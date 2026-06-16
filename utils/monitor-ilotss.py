#!/usr/bin/env python3
# Presently working in 3.9 due to cx_oracle issue with 3.10+
# Version 0.1 TWS April 2026
import argparse
import os
import shutil
import subprocess
import sys
import re
import pickle
import stager_access
import time
import glob
import numpy as np
import pyrap.tables as pt
from pyrap.tables import table
from pathlib import Path
import math
import logging
from assess_calibrators import assess_and_compare

# Require key environment variables to be set before anything else.
_missing_vars = []
RCLONE_CONFIG_DIR = os.environ.get("RCLONE_CONFIG_DIR")
if not RCLONE_CONFIG_DIR:
    _missing_vars.append("RCLONE_CONFIG_DIR  – directory containing rclone macaroon .conf files")

DDF_PIPELINE_CATALOGS = os.environ.get("DDF_PIPELINE_CATALOGS")
if not DDF_PIPELINE_CATALOGS:
    _missing_vars.append("DDF_PIPELINE_CATALOGS  – path to DDF pipeline bootstrap catalogues")

FLOCS_LTA_ROOT = os.environ.get("FLOCS_LTA_ROOT")
if not FLOCS_LTA_ROOT:
    _missing_vars.append("FLOCS_LTA_ROOT  – path to the flocs_lta package directory")

SLURM_ACCOUNT = os.environ.get("SLURM_ACCOUNT")
if not SLURM_ACCOUNT:
    _missing_vars.append("SLURM_ACCOUNT  – SLURM account name for job submission")

SLURM_QUEUE = os.environ.get("SLURM_QUEUE")
if not SLURM_QUEUE:
    _missing_vars.append("SLURM_QUEUE  – SLURM partition/queue name for job submission")

if _missing_vars:
    print(
        "ERROR: The following required environment variables are not set:\n"
        + "\n".join(f"  {v}" for v in _missing_vars)
        + "\n\nExample:\n"
        "  export RCLONE_CONFIG_DIR=/project/lotss/Software/ILoTSS/monitoring/macaroons/\n"
        "  export DDF_PIPELINE_CATALOGS=/project/lotss/Software/ddf-operations/DDF-RUNS/bootstrap-cats/bootstrap-cats\n"
        "  export FLOCS_LTA_ROOT=/project/lotss/Software/ILoTSS/monitoring/flocs-lta/flocs_lta/\n"
        "  export SLURM_ACCOUNT=lotss\n"
        "  export SLURM_QUEUE=normal",
        file=sys.stderr,
    )
    sys.exit(1)

# For approach see https://github.com/LOFAR-VLBI/lofar_vlbi_helpers/tree/main/edfn
#run_linc_calibrator # In flocs run 
#run_linc_target # In flocs run
#run_ddf # In flocs run
#run_delaycal # In flocs run.
#run_ddcal_chunked (faster than default) # In flocs-run§
#1.2" imaging
#facet-subtraction (splitting out MS for each facet)
#imaging (imaging each facet individually)


def slurm_job_status(jobid):
    """Return SLURM job state string or None if unknown."""
    # Try squeue (live)
    try:
        res = subprocess.run(
            ["squeue", "-j", str(jobid), "-h", "-o", "%T"],
            capture_output=True, text=True, check=False
        )
    except FileNotFoundError:
        res = None

    if res and res.returncode == 0 and res.stdout.strip():
        return res.stdout.strip()

    # Fallback to sacct (historical; may require accounting enabled)
    try:
        res2 = subprocess.run(
            ["sacct", "-j", str(jobid), "-n", "-o", "State"],
            capture_output=True, text=True, check=False
        )
    except FileNotFoundError:
        res2 = None

    if res2 and res2.returncode == 0 and res2.stdout.strip():
        # sacct may return multiple lines; take first non-empty and normalize
        for line in res2.stdout.splitlines():
            s = line.strip()
            if s:
                return s.split()[0]  # e.g. COMPLETED, FAILED, RUNNING, PENDING
    return None

def is_job_running(jobid):
    st = slurm_job_status(jobid)
    return st is not None and st.upper() in ("RUNNING","PENDING","CONFIGURING","COMPLETING")


def extract_sbatch_jobid(res):
    text = (getattr(res, "stdout", "") or "") + "\n" + (getattr(res, "stderr", "") or "")
    matches = re.findall(r"Submitted batch job (\d+)", text)
    return int(matches[-1]) if matches else None


def run_with_retries(cmd, attempts=3, delay=60, cwd=None):
    """
    Run `cmd` up to `attempts` times, sleeping `delay` seconds between failures.
    Returns CompletedProcess on success; raises CalledProcessError after final failure.
    """
    last = None
    for attempt in range(1, attempts + 1):
        logging.info("Running (attempt %d/%d): %s", attempt, attempts, cmd)
        last = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
        logging.debug("STDOUT:\n%s", last.stdout.strip())
        logging.debug("STDERR:\n%s", last.stderr.strip())
        if last.returncode == 0:
            logging.info("Command succeeded on attempt %d.", attempt)
            return last
        logging.warning("Command failed (exit %d).", last.returncode)
        if attempt < attempts:
            logging.info("Retrying in %d seconds...", delay)
            time.sleep(delay)
    # all attempts failed -> raise (behaviour like check=True)
    raise subprocess.CalledProcessError(last.returncode, cmd, output=last.stdout, stderr=last.stderr)


def is_corrupt_ms(ms_path):
  if not os.path.isdir(ms_path):
    return True, "path is missing or not a directory"

  required_subtables = (
    "ANTENNA",
    "FIELD",
    "OBSERVATION",
    "SPECTRAL_WINDOW",
    "POLARIZATION",
    "DATA_DESCRIPTION",
  )

  t = None
  try:
    t = pt.table(ms_path, readonly=True, ack=False)
    nrows = t.nrows()
    cols = set(t.colnames())

    # Touch one data cell where possible to catch truncated payloads
    # that still allow table metadata to open.
    if nrows > 0:
      for dcol in ("DATA", "CORRECTED_DATA", "FLOAT_DATA"):
        if dcol in cols:
          t.getcell(dcol, 0)
          break

    for sub in required_subtables:
      sub_path = os.path.join(ms_path, sub)
      st = None
      try:
        st = pt.table(sub_path, readonly=True, ack=False)
        if st.nrows() <= 0:
          return True, f"subtable {sub} has no rows"
        # Read one representative value for lightweight content validation.
        sub_cols = set(st.colnames())
        if sub == "FIELD" and "PHASE_DIR" in sub_cols:
          st.getcell("PHASE_DIR", 0)
        elif sub == "ANTENNA" and "NAME" in sub_cols:
          st.getcell("NAME", 0)
      finally:
        if st is not None:
          st.close()
  except Exception as e:
    return True, str(e)
  finally:
    if t is not None:
      t.close()

  return False, ""


def check_and_remove_corrupt_ms(target):
  """Validate MS directories and remove corrupt ones.

  `target` may be one of:
  - SAS ID (e.g. "123456"): checks L123456/*.MS and L123456/*.ms
  - direct MS path (e.g. /path/to/file.ms): checks that one MS
  - directory path: checks *.MS/*.ms inside that directory

  Returns the number of corrupt MS directories removed.
  """

  target_str = str(target)

  if target_str.lower().endswith(".ms"):
    ms_paths = [target_str]
    checked_scope = target_str
  elif os.path.isdir(target_str):
    ms_paths = sorted(
      glob.glob(os.path.join(target_str, "*.MS"))
      + glob.glob(os.path.join(target_str, "*.ms"))
    )
    checked_scope = target_str
  else:
    ms_dir = f"L{target_str}"
    ms_paths = sorted(
      glob.glob(os.path.join(ms_dir, "*.MS"))
      + glob.glob(os.path.join(ms_dir, "*.ms"))
    )
    checked_scope = ms_dir

  if not ms_paths:
    logging.info("No MS files found to validate in %s.", checked_scope)
    return 0

  removed = 0
  for ms_path in ms_paths:
    is_bad, reason = is_corrupt_ms(ms_path)
    if not is_bad:
      continue
    logging.warning("Corrupt MS detected (%s): %s - removing.", ms_path, reason)
    try:
      shutil.rmtree(ms_path)
      logging.info("Removed corrupt MS: %s", ms_path)
      removed += 1
    except Exception as rm_err:
      logging.warning("Failed to remove corrupt MS %s: %s", ms_path, rm_err)

  if removed:
    logging.info("Removed %d corrupt MS file(s) from %s.", removed, checked_scope)
  else:
    logging.info("No corrupt MS files found in %s.", checked_scope)
  return removed


def download_with_retries(cmd, targid, attempts=3, delay=60, cwd=None):
  """
  Run a download command with retries. Between failed attempts remove any
  intermediate tar files under directory `L{targid}` (pattern `L{targid}/*tar`).

  Returns the CompletedProcess on success, otherwise raises CalledProcessError.
  """
  last = None
  for attempt in range(1, attempts + 1):
    logging.info("Running download (attempt %d/%d): %s", attempt, attempts, cmd)
    last = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    logging.debug("STDOUT:\n%s", last.stdout.strip())
    logging.debug("STDERR:\n%s", last.stderr.strip())
    if last.returncode == 0:
      logging.info("Download succeeded on attempt %d.", attempt)
      targ_tar_pattern = os.path.join(f"L{targid}", "*tar")
      for p in glob.glob(targ_tar_pattern):
        try:
          os.remove(p)
          logging.info("Removed residual tar file %s", p)
        except Exception as e:
          logging.warning("Failed to remove tar file %s: %s", p, e)
      return last
    logging.warning("Download failed (exit %d).", last.returncode)
    if attempt < attempts:
      check_and_remove_corrupt_ms(targid)
      targ_tar_pattern = os.path.join(f"L{targid}", "*tar")
      for p in glob.glob(targ_tar_pattern):
        try:
          os.remove(p)
          logging.info("Removed intermediate file %s", p)
        except Exception as e:
          logging.warning("Failed to remove %s: %s", p, e)
      logging.info("Retrying in %d seconds...", delay)
      time.sleep(delay)
  # all attempts failed -> raise (behaviour like check=True)
  raise subprocess.CalledProcessError(last.returncode, cmd, output=last.stdout, stderr=last.stderr)


def download_with_retries_slurm(cmd, targid, attempts=3, delay=60, cwd=None, cores=1, walltime="120:00:00", script_dir="slurm_jobs", poll_interval=60):
  """
  Submit `cmd` as a SLURM job and wait for completion, with retries.

  Behaviour mirrors `download_with_retries` but uses SLURM submission via
  `submit_slurm_command` and `submit_and_get_jobid`. Between failed attempts
  remove intermediate tar files under directory `L{targid}` (pattern
  `L{targid}/*tar`). Returns a `subprocess.CompletedProcess`-like object on
  success, otherwise raises `subprocess.CalledProcessError`.
  """
  last_state = None
  for attempt in range(1, attempts + 1):
    logging.info("Submitting download as SLURM job (attempt %d/%d): %s", attempt, attempts, cmd)
    # prepare script path
    script_name = os.path.join(script_dir, f"download_L{targid}_attempt{attempt}.sh")
    submit_slurm_command(cmd, out_path=script_name, cores=cores, walltime=walltime)
    try:
      jobid = submit_and_get_jobid(script_name)
    except Exception as e:
      logging.warning("Failed to submit SLURM job: %s", e)
      jobid = None

    if jobid:
      logging.info("Submitted SLURM job %s, waiting for completion...", jobid)
      # wait for job to finish
      while is_job_running(jobid):
        time.sleep(poll_interval)
      last_state = slurm_job_status(jobid)
      logging.info("SLURM job %s finished with state: %s", jobid, last_state)
      if last_state and last_state.upper().startswith("COMPLETED"):
        targ_tar_pattern = os.path.join(f"L{targid}", "*tar")
        for p in glob.glob(targ_tar_pattern):
          try:
            os.remove(p)
            logging.info("Removed residual tar file %s", p)
          except Exception as e:
            logging.warning("Failed to remove tar file %s: %s", p, e)
        return subprocess.CompletedProcess(args=cmd, returncode=0, stdout=f"Submitted via SLURM job {jobid}", stderr="")
    else:
      logging.warning("No job id returned for SLURM submission.")

    logging.warning("SLURM download failed (state=%s).", last_state)
    if attempt < attempts:
      check_and_remove_corrupt_ms(targid)
      targ_tar_pattern = os.path.join(f"L{targid}", "*tar")
      for p in glob.glob(targ_tar_pattern):
        try:
          os.remove(p)
          logging.info("Removed intermediate file %s", p)
        except Exception as e:
          logging.warning("Failed to remove %s: %s", p, e)
      logging.info("Retrying in %d seconds...", delay)
      time.sleep(delay)

  # all attempts failed -> raise
  raise subprocess.CalledProcessError(1, cmd, output="", stderr=f"SLURM job final state: {last_state}")

def get_ms_phasecentre(msfile, field_id=0):
    """
    Get RA and DEC (in degrees) from a Measurement Set.

    Parameters
    ----------
    msfile : str
        Path to the Measurement Set
    field_id : int, optional
        Field index (default = 0)

    Returns
    -------
    ra_deg : float
        Right Ascension in degrees
    dec_deg : float
        Declination in degrees
    """
    t = table(f"{msfile}/FIELD")
    phase_dir = t.getcol("PHASE_DIR")
    t.close()

    ra_rad = phase_dir[field_id, 0, 0]
    dec_rad = phase_dir[field_id, 0, 1]

    return math.degrees(ra_rad), math.degrees(dec_rad)


def make_ddf_parset(output_path, clusterfile=None, exitafter=None, include_inputmodel=True):
    """
    Create a parset-style config file with predefined content.

    Parameters
    ----------
    output_path : str or Path
        Path to the output file.
    clusterfile : str or None
        Value for clusterfile in the [image] section. Omitted when None.
    exitafter : str or None
        When set, adds ``exitafter=<value>`` to the [control] section
        (e.g. 'dirin' for an initial catalogue-only run).
    include_inputmodel : bool
        Whether to include the [inputmodel] section (requires DR3 model files).
    """
    output_path = Path(output_path)

    clusterfile_line = f"clusterfile = {clusterfile}\n" if clusterfile is not None else ""
    exitafter_line = f"exitafter={exitafter}\n" if exitafter is not None else ""
    inputmodel_section = """
[inputmodel]
basedicomodel=DR3-run/image_full_ampphase_di_m.NS
basemaskname=DR3-run/image_full_ampphase_di_m.NS.mask01.fits
baseimagename=DR3-run/image_full_ampphase_di_m.NS
""" if include_inputmodel else ""

    content = f"""[data]
mslist=mslist.txt
full_mslist=big-mslist.txt
colname=DATA

[image]
imsize=20000
robust=-0.15
psf_arcsec=12.0
final_robust=-0.5
final_psf_arcsec=6.0
do_decorr=True
final_rmsfactor=1.0
{clusterfile_line}
[control]
restart=True
bootstrap=False
polcubes=False
stokesv=False
spectral_restored=False
clearcache=False
{exitafter_line}
[masking]
tgss=$$/TGSSADR1_7sigma_catalog.fits
extended_size=2000
thresholds=[15,5,5,5]
rmsfacet = True

[bootstrap]
catalogues=['$$/VLSS.fits','$$/wenss.fits']
names=['VLSS','WENSS']
radii=[40,10]
frequencies=[74e6,327e6]

[solutions]
ndir=45
uvmin=.1
normalize=['None','None','None']
NIterKF=[6, 6, 6, 6, 6, 6, 6]
dt_slow=1.
dt_fast=0.5
dt_di=0.2

[spectra]
do_dynspec=False

[offsets]
method=None

[compression]
compress_polcubes=True
delete_compressed=True
compress_ms=False
{inputmodel_section}"""

    output_path.write_text(content)
    return output_path

def fetch_dr3_model(ms_path, outdir):
  """Fetch LoTSS-DR3 model tarballs for the field and unpack them.

  Returns a set of tarball names that were successfully downloaded and
  extracted (e.g. {'images.tar', 'uv_misc.tar', 'misc.tar'}).  Each tar
  is attempted independently so a partial result is possible.

  'images.tar'   -> input model (image_full_ampphase_di_m.NS.*)
  'uv_misc.tar'  -> cluster catalogue (image_dirin_SSD_m.npy.ClusterCat.npy)
  """
  os.makedirs(outdir, exist_ok=True)
  t_ms = pt.table(ms_path + '/OBSERVATION', readonly=True, ack=False)
  lotssname = str(t_ms[0]['LOFAR_TARGET'][0])
  t_ms.close()
  ddf_conf = os.path.join(RCLONE_CONFIG_DIR, 'maca_sksp_tape_DDF_readonly.conf')
  fetched = set()
  for des in ['images.tar', 'uv_misc.tar', 'misc.tar']:
    logging.info('Fetching LoTSS-DR3 model file %s for %s with rclone...' % (des, lotssname))
    if not os.path.exists(f'{outdir}/{des}'):
      cmd = f'rclone --config={ddf_conf} copy maca_sksp_tape_DDF_readonly:/archive/{lotssname}/{des} {outdir}/'
      res = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
      if res.returncode != 0:
        logging.warning(f"Failed to fetch {des} for {lotssname}\nstderr: {res.stderr.strip()}")
        continue
    logging.info(f"Successfully fetched {des} for {lotssname}")
    res = subprocess.run(f'tar -xf {outdir}/{des} -C {outdir}/', shell=True, check=True, capture_output=True, text=True)
    fetched.add(des)

  if 'images.tar' in fetched:
    with open(f'{outdir}/image_full_ampphase_di_m.NS.DicoModel', 'rb') as f:
        d = pickle.load(f, encoding='latin-1')
    logging.info('Fixing keywords in DicoModel for compatibility with current DDFacet version...')
    d['GD']['Weight']['EnableSigmoidTaper'] = 0
    d['GD']["Weight"]["SigmoidTaperInnerCutoff"] = 0
    d['GD']["Weight"]["SigmoidTaperOuterCutoff"] = 0
    d['GD']["Weight"]["SigmoidTaperInnerRolloffStrength"] = 0.5
    d['GD']["Weight"]["SigmoidTaperOuterRolloffStrength"] = 0.5
    d['GD']["Data"]["Dask"] = False
    with open(f'{outdir}/new.dicomodel', 'wb') as f:
        pickle.dump(d, f)
    os.rename(f'{outdir}/new.dicomodel', f'{outdir}/image_full_ampphase_di_m.NS.DicoModel')
    logging.info('DicoModel keywords updated and saved.')

  return fetched

def _fmt_seconds(s):
  h = s // 3600
  m = (s % 3600) // 60
  return f"{h}h{m}m"
  
def count_srm_lines(path, calib_sas):
  if not path or not os.path.isfile(path):
      return 0
  substr = f"/L{calib_sas}_"
  with open(path, 'r') as fh:
      return sum(1 for line in fh if substr in line)

def extract_flocs_staging_info(text: str):
  """Extract all calibrator SAS IDs and staging IDs for calibrator/target.

  Calibrators are listed in the output as numbered blocks:
    == Closest calibrator observation #N ==
    ...
    SAS ID: <id>

  Returns (calib_sas_list, targ_sas, calib_stage_id, targ_stage_id) where
  calib_sas_list is a list of SAS ID strings ordered by proximity (may be
  empty).  The other three values are strings or None.
  """
  calib_stage = None
  targ_stage = None
  targ_sas = None

  # All calibrator blocks: numbered "== Closest calibrator observation #N =="
  # Capture the SAS ID that immediately follows each block header.
  calib_sas_list = re.findall(
    r"==\s*Closest calibrator observation\s*#?\d*\s*==.*?SAS ID:\s*(\d+)",
    text, re.S
  )

  # Target SAS: usually appears in the initial block before 'Obtaining SURLs'
  pos = text.find("Obtaining SURLs")
  if pos != -1:
    m_t = re.search(r"SAS ID:\s*(\d+)", text[:pos])
    if m_t:
      targ_sas = m_t.group(1)
  # fallback: any SAS ID occurrence
  if not targ_sas:
    m_any = re.search(r"SAS ID:\s*(\d+)", text)
    if m_any:
      targ_sas = m_any.group(1)

  # Staging IDs near the "Staging calibrator data" / "Staging target data" lines
  m_cal = re.search(r"Staging calibrator data.*?Staging request submitted with staging ID\s*(\d+)", text, re.S)
  if m_cal:
    calib_stage = m_cal.group(1)
  m_targ = re.search(r"Staging target data.*?Staging request submitted with staging ID\s*(\d+)", text, re.S)
  if m_targ:
    targ_stage = m_targ.group(1)

  # Fallback: use first/second occurrence if context-based fail
  if not (calib_stage and targ_stage):
    all_ids = re.findall(r"Staging request submitted with staging ID\s*(\d+)", text)
    if len(all_ids) >= 2:
      calib_stage = calib_stage or all_ids[0]
      targ_stage = targ_stage or all_ids[1]
    elif len(all_ids) == 1:
      calib_stage = calib_stage or all_ids[0]

  return calib_sas_list, targ_sas, calib_stage, targ_stage

def extract_sasid_from_srm_list(path: str):
  """Read an SRM list file and return a list of all unique SAS IDs found.

  A SAS ID is the number in patterns like /L123456_SB or /L123456/.
  Returns a list of unique IDs (in order of first appearance),
  or an empty list if the file is missing or nothing matches.
  """
  if not os.path.isfile(path):
    logging.info(f"  [extract_sasid] file not found: {path}")
    return []

  seen = {}  # dict preserves insertion order (Python 3.7+)
  with open(path, 'r') as fh:
    for lineno, line in enumerate(fh, 1):
      m = re.search(r"/L(\d+)", line)
      if m:
        sas = m.group(1)
        if sas not in seen:
          seen[sas] = None
  return list(seen.keys())


def find_surls(srm_file: str, calib_sas: str):
  """Return a list of SURLs from `srm_file` that belong to `calib_sas`.

  Matches any line containing `/L{calib_sas}_` (the standard LOFAR naming
  convention).  Returns an empty list if the file is missing or no lines
  match.
  """
  if not srm_file or not os.path.isfile(srm_file):
    logging.warning(f"find_surls: file not found: {srm_file}")
    return []
  substr = f"/L{calib_sas}_"
  surls = []
  with open(srm_file, 'r') as fh:
    for line in fh:
      line = line.strip()
      if line and substr in line:
        surls.append(line)
  return surls


def check_ms_vs_srm(sas: str, srmlist_path: str, filter_substr: str = None):
  """Compare number of `.MS` entries under L{sas} with SRM list lines.

  - `sas`: SAS id (digits as string)
  - `srmlist_path`: path to srm list file (may be None)
  - `filter_substr`: if provided, only count SRM lines containing this substring

  Returns (ms_count, srm_count).
  """
  import glob
  ms_dir = f"L{sas}"
  ms_count = 0
  if os.path.isdir(ms_dir):
    ms_count = len(glob.glob(os.path.join(ms_dir, "*.MS")))
  srm_count = 0
  if srmlist_path and os.path.isfile(srmlist_path):
    with open(srmlist_path, 'r') as fh:
      for line in fh:
        line = line.strip()
        if not line:
          continue
        if filter_substr:
          if filter_substr in line:
            srm_count += 1
        else:
          srm_count += 1
  return ms_count, srm_count


def submit_slurm_command(cmd: str, out_path: str, cores: int = 1, walltime: str = "50:00:00", job_name: str = None):
  """Write a SLURM jobscript to `out_path` containing `cmd`.

  Parameters:
  - `cmd`: command to run inside the jobscript (string)
  - `out_path`: destination path for the jobscript
  - `cores`: value for `#SBATCH -c`
  - `walltime`: value for `#SBATCH -t` (HH:MM:SS)
  - `job_name`: optional value for `#SBATCH --job-name`

  The script matches the requested template and the resulting file is made
  executable. Returns the path to the written script.
  """
  job_name_line = f"#SBATCH --job-name={job_name}\n" if job_name else ""
  script = f"""#!/bin/bash
#SBATCH -N 1
#SBATCH -c {cores}
#SBATCH -t {walltime}
{job_name_line}##SBATCH --constraint=rome
unset DDF_PIPELINE_DATABASE
export SINGULARITYENV_RCLONE_CONFIG_DIR={RCLONE_CONFIG_DIR}
export SINGULARITYENV_HOSTNAME=spider
export SINGULARITYENV_DDF_PIPELINE_CATALOGS={DDF_PIPELINE_CATALOGS}

{cmd}
"""
  d = os.path.dirname(out_path)
  if d and not os.path.isdir(d):
    os.makedirs(d, exist_ok=True)
  with open(out_path, 'w') as fh:
    fh.write(script)
  try:
    os.chmod(out_path, 0o755)
  except Exception:
    pass
  return out_path


def find_running_job_by_name(job_name: str):
  """Return the job id (int) of a RUNNING or PENDING SLURM job with the given
  name, or None if no such job exists."""
  try:
    res = subprocess.run(
        ["squeue", "--name", job_name, "-h", "-o", "%i"],
        capture_output=True, text=True, check=False
    )
  except FileNotFoundError:
    return None
  if res.returncode == 0:
    for line in res.stdout.splitlines():
      jid = line.strip()
      if jid.isdigit():
        return int(jid)
  return None


def submit_and_get_jobid(script_path: str):
  """Submit `script_path` with `sbatch` and return the job id (int)"""

  res = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
  m = re.search(r"Submitted batch job (\d+)", res.stdout)
  return int(m.group(1))


def run_ddf_pipeline(label, parset, check_glob, workdir, singularity_img, sing_home_dir, sing_bind_dir, max_attempts=3):
  """Submit a DDF pipeline.py job via SLURM and wait for `check_glob` to appear.

  Retries up to `max_attempts` times.  Exits the process on total failure.
  """
  job_name = f"ddf_{label}"
  script = f"slurm_jobs/run_ddf_{label}.sh"
  full_cmd = f"cd {workdir} && singularity exec -H {sing_home_dir} -B {sing_bind_dir} {singularity_img} pipeline.py {parset}"
  for attempt in range(1, max_attempts + 1):
    if glob.glob(check_glob):
      break
    job_id = find_running_job_by_name(job_name)
    if job_id is not None:
      logging.info("DDF %s job '%s' already running/pending as job %d; reusing (attempt %d/%d).", label, job_name, job_id, attempt, max_attempts)
    else:
      logging.info("Submitting DDF %s command (attempt %d/%d): %s", label, attempt, max_attempts, full_cmd)
      submit_slurm_command(full_cmd, out_path=script, cores=61, walltime="120:00:00", job_name=job_name)
      job_id = submit_and_get_jobid(script)
    while is_job_running(job_id):
      logging.info("Waiting for DDF %s job %d to finish (attempt %d/%d)...", label, job_id, attempt, max_attempts)
      time.sleep(60)
    if glob.glob(check_glob):
      break
    logging.warning("DDF %s job %d finished but check_glob '%s' not satisfied (attempt %d/%d).", label, job_id, check_glob, attempt, max_attempts)
    if attempt == max_attempts:
      logging.error("DDF %s failed after %d attempts. Cannot continue.", label, max_attempts)
      sys.exit(1)
    logging.info("Resubmitting DDF %s...", label)


def filter_running_jobs(job_ids):
  """Return subset of `job_ids` that are still running according to `squeue`.

  Simple behaviour: skip falsy entries, call `squeue -j <id> -h`, and include
  the id in the result if `squeue` produces any stdout. If `squeue` is not
  available the function returns an empty list.
  """
  running = []
  for jid in job_ids:
    if not jid:
      continue
    try:
      res = subprocess.run(["squeue", "-j", str(jid), "-h"], capture_output=True, text=True)
    except FileNotFoundError:
      return []
    if res.returncode == 0 and res.stdout.strip():
      running.append(jid)
  return running


def rclone_restore_calibrator(calib_sas, rclone_config=None,
                              rclone_dest="maca_ilotss:LINC_calibrators"):
  """Download and untar a LINC calibrator run from the rclone archive if needed.

  For the given calibrator SAS ID, lists all matching tars on the remote
  (LINC_calibrator_L{calib_sas}_*.tar).  For each remote tar found:
  - skips if the corresponding local directory already exists
  - otherwise downloads the tar and extracts it, then removes the local tar
  """
  if rclone_config is None:
    rclone_config = os.path.join(RCLONE_CONFIG_DIR, "maca_ilotss.conf")
  # List remote tars matching this calibrator SAS ID.
  ls_res = subprocess.run(
      f"rclone --config={rclone_config} lsf {rclone_dest}",
      shell=True, capture_output=True, text=True)
  if ls_res.returncode != 0:
    logging.warning("rclone_restore_calibrator: could not list %s: %s", rclone_dest, ls_res.stderr.strip())
    return

  prefix = f"LINC_calibrator_L{calib_sas}_"
  remote_tars = [f.strip() for f in ls_res.stdout.splitlines()
                 if f.strip().startswith(prefix) and f.strip().endswith('.tar')]

  if not remote_tars:
    logging.info("rclone_restore_calibrator: no remote tar found for L%s.", calib_sas)
    return

  for tar_name in remote_tars:
    local_dir = tar_name[:-4]  # strip .tar
    if os.path.isdir(local_dir):
      logging.info("rclone_restore_calibrator: %s already exists locally, skipping.", local_dir)
      continue

    if os.path.isfile(tar_name):
      logging.info("Local tar %s already exists, skipping download.", tar_name)
    else:
      logging.info("Downloading %s from %s ...", tar_name, rclone_dest)
      res = subprocess.run(
          f"rclone --config={rclone_config} copy {rclone_dest}/{tar_name} .",
          shell=True, capture_output=True, text=True)
      if res.returncode != 0:
        logging.warning("Download failed for %s (exit %d): %s", tar_name, res.returncode, res.stderr.strip())
        continue

    logging.info("Extracting %s ...", tar_name)
    res = subprocess.run(f"tar -xf {tar_name}", shell=True, capture_output=True, text=True)
    if res.returncode != 0:
      logging.warning("Extraction failed for %s (exit %d): %s", tar_name, res.returncode, res.stderr.strip())
      continue
    logging.info("Restored %s.", local_dir)

    try:
      os.remove(tar_name)
      logging.info("Removed local tar %s.", tar_name)
    except Exception as e:
      logging.warning("Could not remove local tar %s: %s", tar_name, e)


def rclone_upload_calibrator(calib_sas, rclone_config=None,
                             rclone_dest="maca_ilotss:LINC_calibrators"):
  """Tar the LINC calibrator run directory and upload it with rclone.

  Skipped silently if:
  - no LINC_calibrator_L{calib_sas}_* directory exists, or
  - cal_solutions.h5 is absent (run did not complete successfully), or
  - the tar already exists on the remote.

  The local tar is removed after a successful upload.
  """
  if rclone_config is None:
    rclone_config = os.path.join(RCLONE_CONFIG_DIR, "maca_ilotss.conf")
  run_dirs = glob.glob(f'LINC_calibrator_L{calib_sas}_*/')
  if not run_dirs:
    logging.warning("rclone_upload_calibrator: no run directory found for L%s, skipping.", calib_sas)
    return
  run_dir = run_dirs[0].rstrip('/')
  sol_h5 = os.path.join(run_dir, 'results_LINC_calibrator', 'cal_solutions.h5')
  if not os.path.exists(sol_h5):
    logging.warning("rclone_upload_calibrator: cal_solutions.h5 missing in %s, skipping.", run_dir)
    return

  tar_name = f'{run_dir}.tar'

  # Skip if already on the remote.
  ls_res = subprocess.run(
      f"rclone --config={rclone_config} ls {rclone_dest}/{tar_name}",
      shell=True, capture_output=True, text=True)
  if ls_res.returncode == 0 and tar_name in ls_res.stdout:
    logging.info("Remote copy of %s already exists, skipping upload.", tar_name)
    return

  # Create tar archive.
  logging.info("Tarring %s -> %s ...", run_dir, tar_name)
  res = subprocess.run(f"tar -cf {tar_name} {run_dir}", shell=True, capture_output=True, text=True)
  if res.returncode != 0:
    logging.warning("Tar failed for %s (exit %d): %s", run_dir, res.returncode, res.stderr.strip())
    return

  # Upload.
  logging.info("Uploading %s to %s ...", tar_name, rclone_dest)
  res = subprocess.run(
      f"rclone --config={rclone_config} copy {tar_name} {rclone_dest}",
      shell=True, capture_output=True, text=True)
  if res.returncode != 0:
    logging.warning("rclone upload failed for %s (exit %d): %s", tar_name, res.returncode, res.stderr.strip())
    return
  logging.info("Uploaded %s successfully.", tar_name)

  # Clean up local tar.
  try:
    os.remove(tar_name)
    logging.info("Removed local tar %s.", tar_name)
  except Exception as e:
    logging.warning("Could not remove local tar %s: %s", tar_name, e)


def rclone_upload_target(targ_sas, rclone_config=None,
                         rclone_dest="maca_ilotss:LINC_target_NL"):
  """Tar the LINC target run directory (excluding large dp3concat MS files) and
  upload it with rclone.

  Skipped silently if:
  - no LINC_target_L{targ_sas}_* directory exists, or
  - the tar already exists on the remote.

  The local tar is removed after a successful upload.
  Files matching ``results_LINC_target/results/*dp3concat`` are excluded from
  the archive to keep its size manageable.
  """
  if rclone_config is None:
    rclone_config = os.path.join(RCLONE_CONFIG_DIR, "maca_ilotss.conf")
  run_dirs = glob.glob(f'LINC_target_L{targ_sas}_*/')
  if not run_dirs:
    logging.warning("rclone_upload_target: no run directory found for L%s, skipping.", targ_sas)
    return
  run_dir = run_dirs[0].rstrip('/')

  tar_name = f'{run_dir}.tar'

  # Skip if already on the remote.
  ls_res = subprocess.run(
      f"rclone --config={rclone_config} ls {rclone_dest}/{tar_name}",
      shell=True, capture_output=True, text=True)
  if ls_res.returncode == 0 and tar_name in ls_res.stdout:
    logging.info("Remote copy of %s already exists, skipping upload.", tar_name)
    return

  # Create tar archive, excluding dp3concat MS files.
  logging.info("Tarring %s -> %s (excluding dp3concat files) ...", run_dir, tar_name)
  exclude = "--exclude='*/results_LINC_target/results/*dp3concat'"
  res = subprocess.run(
      f"tar -cf {tar_name} {exclude} {run_dir}",
      shell=True, capture_output=True, text=True)
  if res.returncode != 0:
    logging.warning("Tar failed for %s (exit %d): %s", run_dir, res.returncode, res.stderr.strip())
    return

  # Upload.
  logging.info("Uploading %s to %s ...", tar_name, rclone_dest)
  res = subprocess.run(
      f"rclone --config={rclone_config} copy {tar_name} {rclone_dest}",
      shell=True, capture_output=True, text=True)
  if res.returncode != 0:
    logging.warning("rclone upload failed for %s (exit %d): %s", tar_name, res.returncode, res.stderr.strip())
    return
  logging.info("Uploaded %s successfully.", tar_name)

  # Clean up local tar.
  try:
    os.remove(tar_name)
    logging.info("Removed local tar %s.", tar_name)
  except Exception as e:
    logging.warning("Could not remove local tar %s: %s", tar_name, e)


def poll_until_output(jid, check_glob, resubmit_cmd, max_resubmits=2, poll_interval=None):
  """
  Wait for SLURM job `jid` to finish, then check whether `check_glob` matches
  any files.  If not, resubmit `resubmit_cmd` and repeat, up to `max_resubmits`
  times.  Returns the path of the first matched file on success, or raises
  RuntimeError if all attempts are exhausted.

  Parameters
  ----------
  jid           : int  – initial SLURM job id
  check_glob    : str  – glob pattern that should match after the job finishes
  resubmit_cmd  : str  – shell command to resubmit the job
  max_resubmits : int  – maximum number of resubmissions (default 2)
  poll_interval : int  – seconds between squeue polls (default: sleep_interval)
  """
  if poll_interval is None:
    poll_interval = sleep_interval

  for attempt in range(max_resubmits + 1):  # attempt 0 = original job
    # Wait for current job to leave the queue (skip if no job was submitted)
    if jid is not None:
      while is_job_running(jid):
        logging.info(f"Job {jid} is still running. Checking again in {poll_interval//60} minutes...")
        time.sleep(poll_interval)
      logging.info(f"Job {jid} has finished (attempt {attempt}).")
    else:
      logging.info(f"No job id for this calibrator (pre-existing run); checking output directly (attempt {attempt}).") 
    matches = glob.glob(check_glob)
    if matches:
      logging.info(f"check_glob '{check_glob}' satisfied: {matches[0]}")
      return matches[0]

    if attempt == max_resubmits:
      raise RuntimeError(
        f"check_glob '{check_glob}' still not satisfied after {max_resubmits} "
        f"resubmit(s). Giving up."
      )

    logging.warning(
      f"check_glob '{check_glob}' not satisfied after job {jid}. "
      f"Resubmitting (attempt {attempt + 1}/{max_resubmits}): {resubmit_cmd}"
    )
    res = subprocess.run(resubmit_cmd, shell=True, check=False, capture_output=True, text=True)
    new_jid = extract_sbatch_jobid(res)
    if new_jid is None:
      raise RuntimeError(f"Resubmit did not return a SLURM job id. stdout: {res.stdout} stderr: {res.stderr}")
    logging.info(f"Resubmitted as job {new_jid}.")
    jid = new_jid


def retry_until_output(cmd, check_glob, restart_rundir_glob, max_attempts=3):
  """
  Run `cmd` synchronously (blocking, e.g. toil) up to `max_attempts` times,
  returning the sorted list of matched paths as soon as `check_glob` is
  satisfied.  Raises RuntimeError after all attempts are exhausted.

  Before every attempt `restart_rundir_glob` is checked: if a tmp run
  directory already exists (including before the very first attempt),
  ``--restart --rundir <dir>`` is appended to `cmd` so a partial run is
  resumed rather than started from scratch.

  Parameters
  ----------
  cmd                  : str  – base shell command to run
  check_glob           : str  – glob that should match after a successful run
  restart_rundir_glob  : str  – glob that locates the tmp run directory used
                                for ``--restart --rundir``
  max_attempts         : int  – total number of attempts (default 3)
  """
  for attempt in range(1, max_attempts + 1):
    rundirs = glob.glob(restart_rundir_glob)
    if rundirs:
      rundir = os.path.abspath(rundirs[0].rstrip('/'))
      _cmd = f"{cmd} --restart --rundir {rundir}"
      logging.info(
        "retry_until_output (attempt %d/%d): tmp dir found, restarting from %s: %s",
        attempt, max_attempts, rundir, _cmd,
      )
    else:
      _cmd = cmd
      logging.info("retry_until_output (attempt %d/%d): %s", attempt, max_attempts, _cmd)
    subprocess.run(_cmd, shell=True, check=False)
    matches = sorted(glob.glob(check_glob))
    if matches:
      logging.info("check_glob '%s' satisfied with %d file(s) on attempt %d.", check_glob, len(matches), attempt)
      return matches
    if attempt < max_attempts:
      logging.warning("check_glob '%s' not satisfied after attempt %d. Retrying...", check_glob, attempt)
  raise RuntimeError(
    f"check_glob '{check_glob}' not satisfied after {max_attempts} attempt(s). Giving up."
  )


def _dp3_uncompress_cmd(ms, outms):
  return (
    f"singularity exec -H {sing_home_dir} -B {sing_bind_dir} {singularity_img} "
    f"DP3 msin={os.path.abspath(ms)} msout={outms} steps=[] "
    f"msout.storagemanager='dysco' msout.uvwcompression=False "
    f"msout.antennacompression=False msout.scalarflags=False"
  )

parser = argparse.ArgumentParser(description="Monitor ILoTSS")
parser.add_argument("--targid", type=int, required=True, help="Target ID (integer)")
parser.add_argument(
    "--exitafter",
    default=None,
    choices=["linccal", "linctarget", "ddfpipeline", "delaycalibration", "ddcalibration"],
    help="Exit cleanly after completing the named stage (default: run all stages).",
)
parser.add_argument(
    "--skipilt",
    action="store_true",
    default=False,
    help="Skip ILT (international baselines): omit --output-fullres-data from LINC target "
         "and automatically exit after the DDF pipeline stage.",
)
parser.add_argument(
    "--skipcleanup",
    action="store_true",
    default=False,
    help="Skip removal of downloaded raw MS data (L{sas}/ directories). "
         "Temporary LINC working directories are still removed regardless.",
)
args, _ = parser.parse_known_args()
targid = args.targid
exitafter = args.exitafter
skipilt = args.skipilt
skipcleanup = args.skipcleanup

sleep_interval = 30 * 60
download_min_fraction = 0.8
sing_pull_dir = os.environ.get("APPTAINER_PULLDIR")
sing_bind_dir = os.getcwd()
sing_home_dir = os.getcwd()
flocs_lta_path = FLOCS_LTA_ROOT
slurm_account = SLURM_ACCOUNT
slurm_queue = SLURM_QUEUE
calib_jobids = {}  # calib_sas -> SLURM job id, or None for pre-existing runs

execdir = f'{os.path.dirname(os.path.abspath(__file__))}/ILoTSS-run-L{targid}'
basedir = f'{os.path.dirname(os.path.abspath(__file__))}'
os.makedirs(execdir, exist_ok=True)
logging.info('Moving to execution directory: %s', execdir)
os.chdir(execdir)

logfile = f"{execdir}/monitor-ilotss.log"
fmt = "%(asctime)s [%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=fmt, handlers=[
    logging.StreamHandler(sys.stdout),
    logging.FileHandler(logfile, mode="a"),
], force=True)

logging.info("Logging to %s", logfile)

text = f"ILoTSS monitor for target {targid}"
width = len(text) + 8
logging.info("#" * width)
logging.info("## " + text + " ##")
logging.info("#" * width)



# Run flocs search and capture output
logging.info("#" * width)
logging.info('Searching the LTA for target and calibrator SURLS')
cmd = f"python {flocs_lta_path}/flocs_search_lta.py --sasid {targid} --freq_end 168.0 --get-surls --num_calibrators 5"
logging.info(f"Running command: {cmd}")
res = run_with_retries(cmd, attempts=3, delay=60)
calib_sas, targ_sas, _, _ = extract_flocs_staging_info(res.stdout)
logging.info(f"Target SAS ID from LTA search: {targ_sas}")
logging.info(f"Calibrator SAS IDs from LTA search: {calib_sas}")

# SRM lists are named srms_<targ_sas>.txt / srms_<targ_sas>_calibrators.txt.
# targ_sas may differ from targid for early-cycle data where LDV renames the SAS ID.
calib_srm_file = f"srms_{targ_sas}_calibrators.txt"
targ_srm_file = f"srms_{targ_sas}.txt"

all_calib_sas_from_srm = extract_sasid_from_srm_list(calib_srm_file) if calib_srm_file else []
targ_sas_from_srm = extract_sasid_from_srm_list(targ_srm_file)[0] if targ_srm_file else None

logging.info(f'Checking SAS IDs from SRM lists: {calib_srm_file} and {targ_srm_file}')
logging.info(f"Target SAS from file: {targ_sas_from_srm}")
logging.info(f"Calibrator SAS IDs from file ({len(all_calib_sas_from_srm)} found): {all_calib_sas_from_srm}")
logging.info(f'Adopting SAS IDs: target={targ_sas_from_srm} calibrators={all_calib_sas_from_srm}')

# Restore any previously archived LINC calibrator runs from the remote.
for _csas in all_calib_sas_from_srm:
    rclone_restore_calibrator(_csas)

targ_num_expected = count_srm_lines(targ_srm_file, targ_sas_from_srm)
logging.info(f'Expected number of files for target ({targ_sas_from_srm}): {targ_num_expected}')
logging.info("#" * width)


logging.info('Checking existing .MS files for target...')
check_and_remove_corrupt_ms(targ_sas_from_srm)
existing_targ_files = glob.glob(f"L{targ_sas_from_srm}/*.MS")
logging.info(f"Existing .MS files for target L{targ_sas_from_srm}: {len(existing_targ_files)} ({len(existing_targ_files)/targ_num_expected:.1%} of expected {targ_num_expected})")
logging.info("#" * width)


for i in range(0,len(all_calib_sas_from_srm)):
  calib_sas_from_srm = all_calib_sas_from_srm[i]
  calib_surls = find_surls(calib_srm_file, calib_sas_from_srm)
  calib_sas_orig = calib_sas[i]
  logging.info(f'Processing calibrator with SAS {calib_sas_from_srm} (original {calib_sas_orig})')

  # If a local tar exists but the directory hasn't been extracted yet, extract it now.
  for _tar in glob.glob(f'LINC_calibrator_L{calib_sas_from_srm}_*.tar'):
    _tar_dir = _tar[:-4]
    if not os.path.isdir(_tar_dir):
      logging.info("Found local tar %s without extracted directory — extracting.", _tar)
      res = subprocess.run(f"tar -xf {_tar}", shell=True, capture_output=True, text=True)
      if res.returncode == 0:
        logging.info("Extracted %s.", _tar_dir)
        try:
          os.remove(_tar)
          logging.info("Removed local tar %s.", _tar)
        except Exception as e:
          logging.warning("Could not remove local tar %s: %s", _tar, e)
      else:
        logging.warning("Extraction of %s failed (exit %d): %s", _tar, res.returncode, res.stderr.strip())

  # If LINC calibrator already completed successfully (cal_solutions.h5 exists), skip
  # downloading and verifying the raw MS data even if L{calib}/ was cleaned up previously.
  linc_calib_already_done = len(glob.glob(
      f'LINC_calibrator_L{calib_sas_from_srm}_*/results_LINC_calibrator/cal_solutions.h5'
  )) > 0
  if linc_calib_already_done:
    logging.info(
        f"Skipping calibrator L{calib_sas_from_srm} download/check: "
        f"cal_solutions.h5 already present from a completed LINC calibrator run."
    )
  else:
    calib_num_expected = count_srm_lines(calib_srm_file, calib_sas_from_srm)
    existing_calib_files = glob.glob(f"L{calib_sas_from_srm}/*.MS")
    num_calib = len(existing_calib_files)
    if calib_num_expected:
      pct_calib = f"{num_calib / calib_num_expected:.1%}"
    else:
      pct_calib = "N/A"
    logging.info(
        f"Existing .MS files for calibrator L{calib_sas_from_srm}: {num_calib} "
        f"({pct_calib} of expected {calib_num_expected})"
    )

    if num_calib >= download_min_fraction * calib_num_expected:
      logging.info(f"Calibrator L{calib_sas_from_srm} already mostly downloaded, skipping staging and download steps.")
    else:
      # Run flocs search and capture output
      logging.info(f'Staging the calibrator {calib_sas_from_srm} (original {calib_sas_orig}) with stager_access')
      # Not using flocs-lta as with products in the compression pipeline it doesnt work because in the LOFAR LTA libraries the compression pipeline is not queriable
      #cmd = f"python {flocs_lta_path}/flocs_search_lta.py --sasid {calib_sas_from_srm} --get-surls --stage-products target --freq_end 168.0"
      #logging.info(f"Running command: {cmd}")
      #res = run_with_retries(cmd, attempts=3, delay=60)
      # Extract IDs from the command output
      #_, targ_sas_dummy, calib_stage, targ_stage = extract_flocs_staging_info(res.stdout)
      calib_stage = stager_access.stage(calib_surls)

      logging.info(
        f"calibrator_stage_id={calib_stage}"
      )

      calib_num_staged = len(stager_access.get_surls_online(calib_stage))
      total_slept = 0
      max_stage_resubmits = 10
      stage_resubmit_count = 0
      skip_calibrator = False
      while calib_num_staged < download_min_fraction * calib_num_expected: # Presently flocs run doesnt do the freq cuts for hte calibrators. Bug to fix.
            logging.info(f"Calibrator stage {calib_stage} status: {calib_num_staged}/{calib_num_expected}; sleeping {sleep_interval//60} minutes")
            time.sleep(sleep_interval)
            total_slept += sleep_interval
            logging.info(f"Total sleep time so far: {_fmt_seconds(total_slept)}")
            calib_num_staged = len(stager_access.get_surls_online(calib_stage))
            stage_status = stager_access.get_status(calib_stage)
            logging.info(f"Staged {calib_num_staged}/{calib_num_expected} files for calibrator stage {calib_stage} (status {stage_status})")
            if stage_status == 'partial success' or stage_status == 'failed':
              if stage_resubmit_count >= max_stage_resubmits:
                logging.warning(
                  f"Staging {calib_stage} is still {stage_status} after {max_stage_resubmits} "
                  f"resubmit(s). Skipping calibrator L{calib_sas_from_srm}."
                )
                skip_calibrator = True
                break
              stage_resubmit_count += 1
              logging.warning(
                f"Staging {calib_stage} is {stage_status}. Resubmitting "
                f"(attempt {stage_resubmit_count}/{max_stage_resubmits})..."
              )
              stager_access.reschedule(calib_stage)
      if skip_calibrator:
        continue

      # Download calibrator data
      logging.info(f'Downloading calibrator data SAS ID {calib_sas_from_srm} (original {calib_sas_orig}) ...')
      cmd = f"python {flocs_lta_path}/flocs_download.py --extract --verification --stage_id {calib_stage}"
      logging.info(f"Running command with SLURM: {cmd}")
      res = download_with_retries_slurm(cmd, targid=calib_sas_from_srm, attempts=3, delay=60)

    # Check download is complete
    # For calibrator SRM list we only count lines that reference the calibrator L{calib_sas}
    check_and_remove_corrupt_ms(calib_sas_from_srm)
    ms_c, srm_c = check_ms_vs_srm(calib_sas_from_srm, calib_srm_file, filter_substr=f"/L{calib_sas_from_srm}_")
    logging.info(f"Calibrator L{calib_sas_from_srm} (original {calib_sas_orig}): .MS files={ms_c}, SRM file entries={srm_c}")
    if ms_c < download_min_fraction * srm_c:
      logging.info(f"FAILING: downloaded number of calibrator files differs significantly from expectations (ms={ms_c} vs srm={srm_c})")
      sys.exit(1)


  # Launch LINC calibrator 
  existingruns = glob.glob(f'tmp.LINC_calibrator_{calib_sas_from_srm}.*/') + glob.glob(f'LINC_calibrator_L{calib_sas_from_srm}_*/') 
  if len(existingruns) == 0:
    logging.info(f'Launching LINC calibrator for SAS ID {calib_sas_from_srm} with flocs-run...')
    cmd = f"flocs-run linc calibrator --slurm-time '4:00:00' --runner cwltool --slurm-cores 32 --scheduler slurm --slurm-account {slurm_account} --slurm-queue {slurm_queue} L{calib_sas_from_srm}"
    logging.info(f"Running command: {cmd}")
    res = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
    calib_jobids[calib_sas_from_srm] = extract_sbatch_jobid(res)
  else:
    logging.info(f"Not running LINC calibrator as found existing LINC calibrator run for SAS ID {calib_sas_from_srm}: {existingruns[0]}")
    calib_jobids[calib_sas_from_srm] = None  # already done; no job to wait for


# Check whether LINC target has already completed (its output MS files exist).
# If so, the raw downloaded data may have been cleaned up on a previous run, but
# there is no need to re-download it.
linc_targ_already_done = (
    len(sorted(glob.glob(
        f'LINC_target_L{targ_sas_from_srm}_*/results_LINC_target/results/*pre-cal.ms'
    ))) > 0.9 * targ_num_expected * 0.1
)
if linc_targ_already_done:
    logging.info(
        f"LINC target output already present for L{targ_sas_from_srm}; "
        "skipping target staging and download steps."
    )

if len(existing_targ_files) < download_min_fraction * targ_num_expected and exitafter != "linccal" and not linc_targ_already_done: # Dont download target if exitafter is linccal or LINC target already ran

  # Run flocs search and capture output
  logging.info('Running flocs_search_lta.py...')
  cmd = f"python {flocs_lta_path}/flocs_search_lta.py --sasid {targid} --get-surls --stage-products target --freq_end 168.0"
  logging.info(f"Running command: {cmd}")
  res = run_with_retries(cmd, attempts=3, delay=60)

  # Extract IDs from the command output
  _, targ_sas_dummy, calib_stage, targ_stage = extract_flocs_staging_info(res.stdout)
  logging.info(
    f"calibrator_stage_id={calib_stage} target_stage_id={targ_stage}"
  )

  targ_num_staged = len(stager_access.get_surls_online(targ_stage))
  total_slept = 0
  max_targ_stage_resubmits = 10
  targ_stage_resubmit_count = 0
  targ_surls = find_surls(targ_srm_file, targ_sas_from_srm)
  while targ_num_staged < download_min_fraction * targ_num_expected:
    logging.info(f"Target stage {targ_stage} status: {targ_num_staged}/{targ_num_expected}; sleeping {sleep_interval//60} minutes")
    time.sleep(sleep_interval)
    total_slept += sleep_interval
    logging.info(f"Total sleep time so far: {_fmt_seconds(total_slept)}")
    targ_num_staged = len(stager_access.get_surls_online(targ_stage))
    logging.info(f"Staged {targ_num_staged}/{targ_num_expected} files for target stage {targ_stage}")
    targ_stage_status = stager_access.get_status(targ_stage)
    if targ_stage_status == 'partial success' or targ_stage_status == 'failed':
      if targ_stage_resubmit_count >= max_targ_stage_resubmits:
        logging.error(
          f"Target staging {targ_stage} is still {targ_stage_status} after "
          f"{max_targ_stage_resubmits} resubmit(s). Cannot continue."
        )
        sys.exit(1)
      targ_stage_resubmit_count += 1
      logging.warning(
        f"Target staging {targ_stage} is {targ_stage_status}. Resubmitting "
        f"(attempt {targ_stage_resubmit_count}/{max_targ_stage_resubmits})..."
      )
      stager_access.reschedule(targ_stage)

  # Download target data
  logging.info(f'Downloading target data SAS ID {targid}...')
  cmd = f"python {flocs_lta_path}/flocs_download.py --extract --verification --stage_id {targ_stage}"
  logging.info(f"Running command with SLURM: {cmd}")
  res = download_with_retries_slurm(cmd, targid=targ_sas_from_srm, attempts=3, delay=60)
else:
  if exitafter == "linccal":
    logging.info('--exitafter linccal: skipping target staging and download steps')
  elif linc_targ_already_done:
    pass  # already logged above
  else:
    logging.info('Target data already mostly downloaded, skipping staging and download steps.')

# Check target download is complete (skip if LINC target already ran and raw data was cleaned up).
if not linc_targ_already_done and exitafter != 'linccal':
  check_and_remove_corrupt_ms(targ_sas_from_srm)
  ms_t, srm_t = check_ms_vs_srm(targ_sas_from_srm, targ_srm_file)
  logging.info(f"Target L{targ_sas_from_srm}: .MS files={ms_t}, SRM list entries={srm_t}")
  if ms_t < download_min_fraction * srm_t:
    logging.error(f"FAILING: target counts differ (ms={ms_t} vs srm={srm_t})")
    sys.exit(1)

### Check the linc calibrator job ids have finished and that some solutions exist.
logging.info('Checking LINC calibrator job ids: %s', calib_jobids)
lincsols = []
for calib, jid in calib_jobids.items():
  if calib is None:
    continue
  resubmit_cmd = (
    f"flocs-run linc calibrator --slurm-time '4:00:00' --runner cwltool "
    f"--slurm-cores 32 --scheduler slurm --slurm-account {slurm_account} "
    f"--slurm-queue {slurm_queue} L{calib}"
  )
  sol_glob = f'LINC_calibrator_L{calib}_*/results_LINC_calibrator/cal_solutions.h5'
  logging.info('Searching for calibrator solutions with glob pattern: %s', sol_glob)
  try:
    sol_path = poll_until_output(jid, sol_glob, resubmit_cmd, max_resubmits=2)
    lincsols.append(sol_path)
    rclone_upload_calibrator(calib)
    tmpdirs = glob.glob(f'tmp.LINC_calibrator_{calib}*/')
    for td in tmpdirs:
      logging.info(f"Removing temporary directory {td} for calibrator L{calib}.")
      shutil.rmtree(td)
  except RuntimeError as e:
    logging.warning(f"Calibrator {calib} failed to produce solutions after all attempts: {e}")
  calib_ms_dir = f"L{calib}"
  if not skipcleanup and os.path.isdir(calib_ms_dir):
    logging.info(f"Removing downloaded calibrator MS files in {calib_ms_dir}/")
    shutil.rmtree(calib_ms_dir)
  elif skipcleanup:
    logging.info(f"--skipcleanup: retaining downloaded MS files in {calib_ms_dir}/.")

    
if not lincsols:
  logging.error("No calibrator solutions found for any calibrator. Cannot continue.")
  sys.exit(1)
logging.info(f"Proceeding with {len(lincsols)} calibrator solution(s): {lincsols}")

_assess_results = assess_and_compare(
    template_path=f'{basedir}/calibrator_medians_cache.h5',
    h5_paths=lincsols,
    n_mad=5.0,
    plot_dir=None,
)
bestcal = f'{execdir}/{_assess_results[0]["h5_path"]}'

logging.info(f"Best calibrator solution selected: {bestcal}")

if exitafter == "linccal":
  logging.info("--exitafter linccal: exiting after LINC calibrator stage.")
  sys.exit(0)

# Launch LINC target
linc_targ_cmd = (
  f"flocs-run linc target --slurm-time '120:00:00' --runner toil --scheduler slurm "
  + ("" if skipilt else "--output-fullres-data ")
  + f"--min-unflagged-fraction 0.05 "
  f"--slurm-account {slurm_account} --slurm-queue {slurm_queue} "
  f"--cal-solutions {bestcal} L{targ_sas_from_srm}"
)
linc_targ_filedir = f'LINC_target_L{targ_sas_from_srm}_*/results_LINC_target/results/'
linc_targ_glob = f'{linc_targ_filedir}*pre-cal.ms'
linc_targ_files = sorted(glob.glob(linc_targ_glob))

if len(linc_targ_files) > 0.9*targ_num_expected*0.1: # Data combined to 10SB blocks
  logging.info(f"Found existing LINC target run with {len(linc_targ_files)} output files for SAS ID {targ_sas_from_srm}.")
else:
  logging.info(f"Launching LINC target for SAS ID {targ_sas_from_srm} with flocs-run...")

  try:
    linc_targ_files = retry_until_output(
      cmd=linc_targ_cmd,
      check_glob=linc_targ_glob,
      restart_rundir_glob=f'tmp.LINC_target_{targ_sas_from_srm}.*/',
      max_attempts=10,
    )
  except RuntimeError as e:
    logging.error(f"LINC target failed after all attempts: {e}")
    sys.exit(1)

if len(linc_targ_files) < 0.9*targ_num_expected*0.1: # Data combined to 10SB blocks
  logging.error(f"Only found {len(linc_targ_files)} LINC target output files, which is less than 90% of the expected {targ_num_expected*0.1}. Cannot continue.")
  sys.exit(1)

logging.info(f"Found {len(linc_targ_files)} LINC target output files for SAS ID {targ_sas_from_srm}.")
# Cleanup
logging.info("Cleaning up temporary files for LINC target stage.")
linctar_tmpdirs = glob.glob(f'tmp.LINC_target_{targ_sas_from_srm}*/')
for ltd in linctar_tmpdirs:
  logging.info(f"Removing temporary directory {ltd} for LINC target.")
  shutil.rmtree(ltd)
targ_ms_dir = f"L{targ_sas_from_srm}"
if not skipcleanup:
  if os.path.isdir(targ_ms_dir):
    logging.info(f"Removing downloaded target MS files in {targ_ms_dir}/")
    shutil.rmtree(targ_ms_dir)
else:
  logging.info(f"--skipcleanup: retaining downloaded MS files in {targ_ms_dir}/.")

rclone_upload_target(targ_sas_from_srm)

if exitafter == "linctarget":
  logging.info("--exitafter linctarget: exiting after LINC target stage.")
  sys.exit(0)


### Doing the DDF-pipeline step.
os.makedirs(f'DDF_{targ_sas_from_srm}', exist_ok=True)
singularity_img = f'{sing_pull_dir}/astronrd_linc_latest.sif'


# Submit DP3 jobs to decompress all LINC target MS files for DDF
active_jobs_ddf = []
for ms in linc_targ_files:
  outms = os.path.abspath(f'DDF_{targ_sas_from_srm}/{os.path.basename(ms)}.uncomp.ms')
  if os.path.exists(outms):
    logging.info(f"Output MS {outms} already exists, skipping.")
    continue
  cmd = _dp3_uncompress_cmd(ms, outms)
  logging.info(f"Submitting DP3 command: {cmd}")
  script = f"DDF_{targ_sas_from_srm}/fix_ms_{os.path.basename(ms)}.sh"
  submit_slurm_command(cmd, out_path=script, cores=1, walltime="4:00:00")
  active_jobs_ddf.append(submit_and_get_jobid(script))
  while len(active_jobs_ddf) > 5:
    active_jobs_ddf = filter_running_jobs(active_jobs_ddf)
    logging.info(f"Throttling submissions; active jobs: {active_jobs_ddf}")
    time.sleep(60)

while active_jobs_ddf:
  active_jobs_ddf = filter_running_jobs(active_jobs_ddf)
  logging.info(f"Waiting for DP3 jobs to finish; active jobs: {active_jobs_ddf}")
  time.sleep(60)

# Check all output files exist and are not corrupt; retry once if missing
for ms in linc_targ_files:
  outms = os.path.abspath(f'DDF_{targ_sas_from_srm}/{os.path.basename(ms)}.uncomp.ms')
  check_and_remove_corrupt_ms(outms)
  if not os.path.exists(outms):
    logging.warning(f"Output MS {outms} not found; submitting retry job...")
    script = f"DDF_{targ_sas_from_srm}/fix_ms_retry_{os.path.basename(ms)}.sh"
    submit_slurm_command(_dp3_uncompress_cmd(ms, outms), out_path=script, cores=1, walltime="4:00:00")
    retry_job_id = submit_and_get_jobid(script)
    logging.info(f"Waiting for retry job {retry_job_id} to finish...")
    while is_job_running(retry_job_id):
      time.sleep(60)
    if not os.path.exists(outms):
      logging.error(f"Output MS {outms} still not found after retry. Cannot continue.")
      sys.exit(1)


# Download existing LoTSS-DR3 models for fast DDF-pipeline.
singularity_img = f'{sing_pull_dir}/ddf.sif'

# Make mslists
workdir = os.path.abspath(f'DDF_{targ_sas_from_srm}')
cmd = f"singularity exec -H {sing_home_dir} -B {sing_bind_dir} {singularity_img} make_mslists.py"
logging.info(f"Running command: {cmd} from working directory {workdir}")
res = subprocess.run(cmd, shell=True, cwd=workdir, check=False, capture_output=True, text=True)

dr3_outdir = f'DDF_{targ_sas_from_srm}/DR3-run'
fetch_dr3_model(outms, dr3_outdir)
has_inputmodel = os.path.exists(f'{dr3_outdir}/image_full_ampphase_di_m.NS.DicoModel')
has_clustercat = os.path.exists(f'{dr3_outdir}/image_dirin_SSD_m.npy.ClusterCat.npy')
logging.info("DR3 file check: has_inputmodel=%s, has_clustercat=%s", has_inputmodel, has_clustercat)

if not has_inputmodel:
  logging.warning("No DR3 input model available. Running full DDF pipeline from scratch.")
  make_ddf_parset(f'DDF_{targ_sas_from_srm}/ddf_full.parset', include_inputmodel=False)
  run_ddf_pipeline(f'full_{targ_sas_from_srm}', 'ddf_full.parset', f'DDF_{targ_sas_from_srm}/summary.txt', workdir, singularity_img, sing_home_dir, sing_bind_dir)
else:
  if not has_clustercat:
    clusterfile_for_reproc = 'image_dirin_SSD_m.npy.ClusterCat.npy'
    if os.path.exists(f'DDF_{targ_sas_from_srm}/{clusterfile_for_reproc}'):
      logging.info("Found cluster catalogue from previous ddf_intial run, reusing for DDF reprocessing.")
    else:
      logging.warning("DR3 cluster catalogue not available but input model is. Running initial DDF pipeline to generate one.")
      make_ddf_parset(f'DDF_{targ_sas_from_srm}/ddf_initial.parset', exitafter='dirin', include_inputmodel=False)
      run_ddf_pipeline(f'initial_{targ_sas_from_srm}', 'ddf_initial.parset', f'DDF_{targ_sas_from_srm}/image_dirin_SSD_m.npy.ClusterCat.npy', workdir, singularity_img, sing_home_dir, sing_bind_dir)
  else:
    clusterfile_for_reproc = 'DR3-run/image_dirin_SSD_m.npy.ClusterCat.npy'

  make_ddf_parset(f'DDF_{targ_sas_from_srm}/VLBI_ddf_reproc.parset', clusterfile=clusterfile_for_reproc, include_inputmodel=True)
  run_ddf_pipeline(f'reproc_{targ_sas_from_srm}', 'VLBI_ddf_reproc.parset', f'DDF_{targ_sas_from_srm}/summary.txt', workdir, singularity_img, sing_home_dir, sing_bind_dir)

if exitafter == "ddfpipeline" or skipilt:
  logging.info("Exiting after DDF pipeline stage%s.",
               " (--skipilt)" if skipilt else " (--exitafter ddfpipeline)")
  sys.exit(0)

logging.info(f'Renaming DDF SOLSDIR')
os.system(f'cd DDF_{targ_sas_from_srm}/SOLSDIR && rename _pre-cal.ms.uncomp.ms .dp3concat *_pre-cal.ms.uncomp.ms')
os.system(f'cp DDF_{targ_sas_from_srm}/DR3-run/image_dirin_SSD_m.npy.ClusterCat.npy DDF_{targ_sas_from_srm}/')


### Run lofar-vlbi-plot
ra_deg, dec_deg = get_ms_phasecentre(linc_targ_files[0])
singularity_img = f'{sing_pull_dir}/astronrd_linc_latest.sif'
cmd = f'singularity exec -H {sing_home_dir} -B {sing_bind_dir} {singularity_img} lofar-vlbi-plot --RA {ra_deg} --DEC {dec_deg} --vlass'
cmd_fallback = f"{cmd} --continue_no_lotss"
used_no_lotss = False

logging.info(f"Running command: {cmd}")
res = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
res_text = ((res.stdout or "") + "\n" + (res.stderr or "")).lower()
lotss_missing = (
  "lotss coverage does not exist" in res_text
  or "contine_no_lotss is set to false" in res_text
  or "continue_no_lotss is set to false" in res_text
)
if res.returncode != 0 or lotss_missing:
  if lotss_missing and res.returncode == 0:
    logging.warning("lofar-vlbi-plot reported missing LoTSS coverage despite zero exit code; retrying with --continue_no_lotss.")
  else:
    logging.warning("lofar-vlbi-plot failed (exit %d); retrying with --continue_no_lotss.", res.returncode)
  if res.stdout.strip():
    logging.warning("lofar-vlbi-plot STDOUT:\n%s", res.stdout.strip())
  if res.stderr.strip():
    logging.warning("lofar-vlbi-plot STDERR:\n%s", res.stderr.strip())
  logging.info(f"Running fallback command: {cmd_fallback}")
  res_fallback = subprocess.run(cmd_fallback, shell=True, check=False, capture_output=True, text=True)
  res_fallback_text = ((res_fallback.stdout or "") + "\n" + (res_fallback.stderr or "")).lower()
  lotss_missing_fallback = (
    "lotss coverage does not exist" in res_fallback_text
    or "contine_no_lotss is set to false" in res_fallback_text
    or "continue_no_lotss is set to false" in res_fallback_text
  )
  if res_fallback.returncode != 0 or lotss_missing_fallback:
    logging.error("lofar-vlbi-plot fallback failed (exit %d).", res_fallback.returncode)
    if res_fallback.stdout.strip():
      logging.error("lofar-vlbi-plot fallback STDOUT:\n%s", res_fallback.stdout.strip())
    if res_fallback.stderr.strip():
      logging.error("lofar-vlbi-plot fallback STDERR:\n%s", res_fallback.stderr.strip())
    sys.exit(1)
  used_no_lotss = True


# Run delay calibration
image_catalogue_arg = "" if used_no_lotss else "--image-catalogue lotss_catalogue.csv "
if used_no_lotss:
  logging.info("lofar-vlbi-plot used --continue_no_lotss; omitting --image-catalogue from delay calibration command.")
delay_cal_cmd = f"flocs-run vlbi delay-calibration --slurm-time '120:00:00' --runner toil --scheduler slurm --ddf-solsdir DDF_{targ_sas_from_srm}/SOLSDIR --ddf-rundir DDF_{targ_sas_from_srm} --do-subtraction --do-validation --ms-suffix 'dp3concat' --delay-calibrator delay_calibrators.csv --mspath {linc_targ_filedir} --slurm-account lotss --slurm-queue normal {image_catalogue_arg} --do-auto-delay-selection --select-best-n-delay-calibrators 1 --starting-skymodel *_vlass.fits --no-use-vlass"
delay_cal_glob = f"VLBI_delay-calibration_{targ_sas_from_srm}_*/results_VLBI_delay-calibration/*ms"
delay_cal_files = sorted(glob.glob(delay_cal_glob))

if len(delay_cal_files) > 0.9*targ_num_expected*0.1: # Data combined to 10SB blocks
  logging.info(f"Found existing VLBI delay calibration run with {len(delay_cal_files)} output files for SAS ID {targ_sas_from_srm}.")
else:
  logging.info(f"Launching VLBI delay calibration for SAS ID {targ_sas_from_srm} with flocs-run...")
  try:
    delay_cal_files = retry_until_output(
      cmd=delay_cal_cmd,
      check_glob=delay_cal_glob,
      restart_rundir_glob=f'tmp.VLBI_delay-calibration.*/',
      max_attempts=10,
    )
  except RuntimeError as e:
    logging.error(f"VLBI delay calibration failed after all attempts: {e}")
    sys.exit(1)

if len(delay_cal_files) < 0.9*targ_num_expected*0.1: # Data combined to 10SB blocks
  logging.error(f"Only found {len(delay_cal_files)} VLBI delay calibration output files, which is less than 90% of the expected {targ_num_expected*0.1}. Cannot continue.")
  sys.exit(1)

logging.info(f"Found {len(delay_cal_files)} VLBI delay calibration output files for SAS ID {targ_sas_from_srm}.")
# Cleanup
logging.info("Cleaning up temporary files for VLBI delay calibration stage.")
delaycal_tmpdirs = glob.glob(f'tmp.VLBI_delay-calibration_{targ_sas_from_srm}*/')
for ltd in delaycal_tmpdirs:
  logging.info(f"Removing temporary directory {ltd} for VLBI delay calibration.")
  shutil.rmtree(ltd)

if exitafter == "delaycalibration":
  logging.info("--exitafter delaycalibration: exiting after VLBI delay calibration stage.")
  sys.exit(0)


# Launch VLBI delay
#flocs-run vlbi delay-calibration \
#--slurm-time "72:00:00" \
#--slurm-queue "normal" \
#--slurm-account lofarvwf \
#--runner toil \
#--scheduler slurm \
#--ddf-solsdir $(realpath ../ddf/SOLSDIR) \
#--ddf-rundir $(realpath ../ddf) \
#--do-subtraction \
#--do-validation \
#--ms-suffix "dp3concat" \
#--apply-delay-solutions \
#--do-auto-delay-selection \
#--delay-calibrator $(realpath delay_calibrators.csv) \
#$(realpath ../target/LINC_target_*/results_LINC_target/results)


# Can run it without slurm. Flocs has restart options.
#cmd = f"flocs-run vlbi delay-calibration --record-toil-stats --scheduler slurm --rundir {rundirs/'rundir'} --outdir {rundirs} --slurm-queue {SLURM_QUEUES} --slurm-time 48:00:00 --slurm-account lofarvlbi --runner toil --delay-calibrator {delay_csv} --ms-suffix dp3concat {linc_target_dir/'results_LINC_target'/'results'}"


logging.info(f'Launching VLBI dd calibration for target SAS ID {targ_sas_from_srm} with flocs-run...')
cmd = f"flocs-run vlbi dd-calibration --slurm-time '120:00:00' --runner cwltool --ddf-solsdir DDF_{targ_sas_from_srm}/SOLSDIR --ddf-rundir DDF_{targ_sas_from_srm} --do-subtraction --do-validation --ms-suffix 'dp3concat' --apply-delay-solutions --do-auto-delay-selection --delay-calibrator delay_calibrators.csv {linc_targ_filedir}"
logging.info(f"Running command: {cmd}")
res = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)

if exitafter == "ddcalibration":
  logging.info("--exitafter ddcalibration: exiting after VLBI DD calibration stage.")
  sys.exit(0)




### THIS DDCAL CAN NOW BE RUN IN FLOCS

# run_ddcal_chunked

#/project/lofarvwf/Share/jdejong/output/EUCLID/edfn_centre/ddcal_dirs/EDFN_directions.csv
#CSV=$1

#SCRIPT_DIR=/home/lofarvwf-jdejong/scripts/lofar_vlbi_helpers/elais_200h
#source /home/lofarvwf-jdejong/scripts/lofar_vlbi_helpers/elais_200h/other/chunk_csv.sh $CSV

# Loop through each CSV chunk file with an index
#ENUMERATE=1
#for CSV_CHUNK in chunk*.csv; do

  # Create run folder
#  RUNFOLDER=chunk_${ENUMERATE}
#  mkdir -p ${RUNFOLDER}
#  mv ${CSV_CHUNK} ${RUNFOLDER}
#  cd ${RUNFOLDER}

  # Submit the job with sbatch
#  sbatch /home/lofarvwf-jdejong/scripts/lofar_vlbi_helpers/edfn/run_ddcal.sh $(realpath "${CSV_CHUNK}")
#  cd ../

  # Increment the enumeration counter
#  ((ENUMERATE++))
#done
