"""
Assess one or more LOFAR calibrator observations against a pre-built template
and report which is best.

Can be used as a standalone script or imported into other code.

Standalone usage:
    python assess_calibrators.py --template calibrator_medians_cache.h5 \\
        --h5 .../L123/cal_solutions.h5 .../L456/cal_solutions.h5

Importable usage:
    from assess_calibrators import assess_and_compare
    results = assess_and_compare(
        template_path='calibrator_medians_cache.h5',
        h5_paths=['/path/to/L123/cal_solutions.h5',
                  '/path/to/L456/cal_solutions.h5'],
        n_mad=5.0,
    )
    # results is a list of dicts, one per h5 file, sorted best-first
"""
import os
import re
import glob
import json
import argparse

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def yy_to_delay_ns(freq_hz, yy_phase):
    """
    Convert a YY polalign phase spectrum to a differential X-Y delay in nanoseconds.
    Masks non-finite samples, unwraps, fits a line φ(ν) = c₀ + c₁ν, and applies
    the radio convention φ(ν) = −2π τ ν  →  τ = −c₁ / (2π).
    Returns np.nan if fewer than 2 finite channels are available.
    """
    mask = np.isfinite(freq_hz) & np.isfinite(yy_phase)
    f, p = freq_hz[mask], yy_phase[mask]
    if len(f) < 2:
        return np.nan
    slope, _intercept = np.polyfit(f, np.unwrap(p), 1)
    return -slope / (2.0 * np.pi) * 1e9   # nanoseconds


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_station_of_interest(ant_name):
    return ant_name.startswith('RS') or ant_name.startswith('CS')


def _get_obs_id(h5_path):
    """Extract the L-number obs ID from the h5 file path."""
    for part in reversed(os.path.abspath(h5_path).split(os.sep)):
        m = re.search(r'(L\d{6,})', part)
        if m:
            return m.group(1)
    return os.path.basename(os.path.dirname(h5_path))


def _get_calibrator_from_h5(h5_path):
    """Read calibrator name directly from cal_solutions.h5."""
    try:
        with h5py.File(h5_path, 'r') as h:
            return h['calibrator/source']['name'][0].decode().upper()
    except Exception:
        return 'Unknown'


# Minimum absolute MAD floor per metric.  Prevents over-flagging when the
# template MAD is near-zero (i.e. a very stable metric across the training
# observations).
MIN_MAD = {
    'clock (ns)':      1.0,    # 1 nanosecond
    'polalign (ns)':   0.5,    # 0.5 ns differential delay
}

# For flagging-fraction metrics (0-100%), flag any station where the observed
# flagging fraction exceeds this absolute threshold (in percent).
FLAG_FRACTION_THRESHOLD = 30.0


def _percentile_flag_check(ant_name, metric, obs_val, ref_vals, issues, all_nsigma):
    """
    Flagging-fraction variant of _flag_check.
    Flags if obs_val (in %) exceeds FLAG_FRACTION_THRESHOLD.
    A continuous pseudo-nsigma = obs_val / FLAG_FRACTION_THRESHOLD is stored
    for scoring (value of 1.0 = exactly at threshold).
    """
    if not np.isfinite(obs_val):
        return
    pseudo_nsigma = obs_val / FLAG_FRACTION_THRESHOLD
    all_nsigma.setdefault(metric, {})[ant_name] = pseudo_nsigma
    if obs_val > FLAG_FRACTION_THRESHOLD:
        # Store in same (metric, obs, ref, scale, nsigma) tuple format as _flag_check
        issues.setdefault(ant_name, []).append(
            (metric, obs_val, FLAG_FRACTION_THRESHOLD, FLAG_FRACTION_THRESHOLD, pseudo_nsigma)
        )


def _mad_flag_check(ant_name, metric, obs_val, med, mad, n_mad, issues, all_nsigma):
    """Flag if deviation from template median exceeds n_mad * MAD (used for bandpass, polalign, clock)."""
    if mad == 0 or not np.isfinite(mad) or not np.isfinite(obs_val):
        return
    scale = max(abs(med), abs(obs_val), 1e-30)
    if mad / scale < 1e-3:
        return
    # Apply a per-metric absolute floor so rock-steady template stations
    # don't get flagged for tiny absolute deviations.
    effective_mad = max(mad, MIN_MAD.get(metric, 0.0))
    nsigma = abs(obs_val - med) / effective_mad
    # Always record for scoring, flag only if above threshold.
    all_nsigma.setdefault(metric, {})[ant_name] = nsigma
    if nsigma > n_mad:
        issues.setdefault(ant_name, []).append((metric, obs_val, med, mad, nsigma))


# ---------------------------------------------------------------------------
# Core assessment function
# ---------------------------------------------------------------------------

def assess_observation(h5_path, template, n_mad=5.0, verbose=True):
    """
    Assess a single cal_solutions.h5 file against the template medians.

    Parameters
    ----------
    h5_path : str
        Path to the cal_solutions.h5 file.
    template : dict
        Dict loaded from the calibrator_medians_cache.pkl produced by build_template.py.
        Keys: obs_to_calibrator, median_bandpass, median_polalign,
              median_clock, median_flagging.
    n_mad : float
        Outlier threshold in units of MAD.
    verbose : bool
        Print the per-antenna report.

    Returns
    -------
    dict with keys:
        obs_id      : str
        calibrator  : str
        issues      : dict  ant_name -> [(metric, obs_val, med, mad, nsigma), ...]
        all_nsigma  : dict  metric -> {ant_name: nsigma}
        h5_path     : str
    """
    median_bandpass = template['median_bandpass']
    median_polalign = template['median_polalign']
    median_clock    = template['median_clock']
    median_flagging = template['median_flagging']

    obs_id     = _get_obs_id(h5_path)
    calibrator = _get_calibrator_from_h5(h5_path)

    if verbose:
        print(f"\nAssessing {obs_id} (calibrator: {calibrator}) against template")
        print(f"Outlier threshold: {n_mad} x MAD")
        print("=" * 80)

    issues = {}
    all_nsigma = {}  # metric -> {ant_name: nsigma}
    observed   = {}  # metric -> {ant_name: observed_value}  (ALL stations, not just flagged)

    def flag(ant, metric, obs_val, med, mad):
        observed.setdefault(metric, {})[ant] = obs_val
        _mad_flag_check(ant, metric, obs_val, med, mad, n_mad, issues, all_nsigma)

    try:
        with h5py.File(h5_path, 'r') as h:
            # ---- Bandpass ----
            if 'calibrator/bandpass' in h:
                bp_ants = [a.decode() if isinstance(a, bytes) else a
                           for a in h['calibrator/bandpass/ant'][...]]
                bp_val  = h['calibrator/bandpass/val'][...]   # (time, freq, ant, pol)
            else:
                bp_ants, bp_val = [], None

            # ---- Polalign ----
            if 'calibrator/polalign' in h:
                pa_ants = [a.decode() if isinstance(a, bytes) else a
                           for a in h['calibrator/polalign/ant'][...]]
                pa_val  = h['calibrator/polalign/val'][...]   # (time, ant, freq, pol)
                pa_freq = h['calibrator/polalign/freq'][...]
            else:
                pa_ants, pa_val, pa_freq = [], None, None

            # ---- Clock ----
            if 'calibrator/clock' in h:
                clk_ants = [a.decode() if isinstance(a, bytes) else a
                            for a in h['calibrator/clock/ant'][...]]
                clk_val  = h['calibrator/clock/val'][...]     # (time, ant), constant in time
            else:
                clk_ants, clk_val = [], None
    except Exception as e:
        print(f"  Error reading {h5_path}: {e}")
        return {'obs_id': obs_id, 'calibrator': calibrator, 'issues': {},
                'all_nsigma': {}, 'h5_path': h5_path}

    # ---- Bandpass ----
    if bp_val is not None and calibrator in {k for d in median_bandpass.values() for k in d}:
        for ai, an in enumerate(bp_ants):
            if not _is_station_of_interest(an): continue
            if an not in median_bandpass or calibrator not in median_bandpass[an]: continue
            entry = median_bandpass[an][calibrator]
            ref_med, ref_mad = entry[1], entry[2]
            sp = np.nanmedian(bp_val[:, :, ai, :], axis=(0, 2))  # (time, freq, ant, pol)
            if len(sp) != len(ref_med): continue
            flag(an, 'bandpass_residual', np.nanmean(np.abs(sp - ref_med)),
                 0.0, np.nanmean(ref_mad) or np.nan)

    # ---- Polalign — YY differential delay (ns) ----
    if pa_val is not None:
        for ai, an in enumerate(pa_ants):
            if not _is_station_of_interest(an): continue
            if an not in median_polalign or calibrator not in median_polalign[an]: continue
            _, ref_vals, ref_mads = median_polalign[an][calibrator]
            yy = pa_val[0, ai, :, 1].astype(float)   # YY at first time slice
            delay_ns = yy_to_delay_ns(pa_freq, yy)
            if np.isnan(delay_ns): continue
            ref_med  = float(np.nanmedian(ref_vals))
            ref_mad  = float(np.nanmedian(ref_mads))
            flag(an, 'polalign (ns)', delay_ns, ref_med, ref_mad)

    # ---- Clock ----
    if clk_val is not None:
        global_clk = np.nanmedian([clk_val[0, i] * 1e9
                                   for i, a in enumerate(clk_ants)
                                   if _is_station_of_interest(a)] or [0.0])
        for ai, an in enumerate(clk_ants):
            if not _is_station_of_interest(an): continue
            if an not in median_clock or calibrator not in median_clock[an]: continue
            _, ref_vals, ref_mads = median_clock[an][calibrator]
            obs_clk_ns = clk_val[0, ai] * 1e9 - global_clk
            flag(an, 'clock (ns)', obs_clk_ns,
                 float(np.nanmedian(ref_vals)), float(np.nanmedian(ref_mads)))

    # ---- JSON flagging fractions (percentage_flagged.final) ----
    parent = os.path.dirname(h5_path)
    json_files = glob.glob(os.path.join(parent, '*_LINC_calibrator_summary.json'))
    if json_files:
        try:
            with open(json_files[0]) as f:
                d = json.load(f)
            final_data = median_flagging.get('final', {})
            for entry in d['metrics']['LINC']['stations']:
                an = entry.get('station', '')
                if not _is_station_of_interest(an): continue
                obs_val = entry.get('percentage_flagged', {}).get('final', np.nan)
                if np.isnan(obs_val): continue
                if an not in final_data or calibrator not in final_data[an]: continue
                _, ref_vals, _ = final_data[an][calibrator]
                observed.setdefault('flag_final', {})[an] = obs_val
                _percentile_flag_check(an, 'flag_final', obs_val,
                                       ref_vals, issues, all_nsigma)
        except Exception as e:
            print(f"  Warning reading JSON: {e}")

    # ---- Print report ----
    if verbose:
        if issues:
            print(f"\nAntennas with outlier solutions (> {n_mad} x MAD):")
            print(f"  {'Station':<14} {'Metric':<22} {'Observed':>12} "
                  f"{'Median':>12} {'MAD':>10} {'N_sigma':>8}")
            print("  " + "-" * 82)
            for an in sorted(issues.keys()):
                for metric, obs_val, med, mad, nsigma in sorted(issues[an], key=lambda x: -x[4]):
                    print(f"  {an:<14} {metric:<22} {obs_val:>12.4f} "
                          f"{med:>12.4f} {mad:>10.4f} {nsigma:>8.1f}")
        else:
            print("  No outlier antennas found.")

    return {'obs_id': obs_id, 'calibrator': calibrator, 'issues': issues,
            'all_nsigma': all_nsigma, 'observed': observed, 'h5_path': h5_path}


# ---------------------------------------------------------------------------
# Scoring and comparison
# ---------------------------------------------------------------------------
#
# How scoring works
# -----------------
# For every station × metric combination, _mad_flag_check / _percentile_flag_check
# compute a CONTINUOUS "n-sigma" value measuring how far the observed value
# deviates from the template.  ALL metrics, including flag_final, contribute
# to the score via this path.
#
#   MAD-based metrics (clock, polalign, bandpass):
#       nsigma = |obs − template_median| / effective_MAD
#     where effective_MAD = max(MAD, MIN_MAD[metric]) to prevent over-sensitivity
#     when the template is very stable.
#
#   Flagging-fraction metric (flag_final):
#       pseudo_nsigma = obs_percent / FLAG_FRACTION_THRESHOLD
#     A value of 1.0 = station is exactly at the warning threshold.
#     A value of 2.0 = station has twice as much flagged data as the threshold.
#     This is stored in all_nsigma exactly like the MAD-based metrics, so
#     flag_final feeds directly into the per-metric score and total score.
#
# Per-metric score
#   = mean(nsigma) across all CS/RS stations for that metric.
# This normalises for station count and rewards observations where the
# *typical* station is close to the template, not just the median one.
#
# Total score
#   = weighted average of per-metric scores, using _WEIGHTS (default all 1.0).
# Lower is better.  Observations are ranked best-first by total score.
#
# NOTE: "Flagged antennas" and "Flagged station count" in the comparison table
# are SEPARATE hard-threshold counts (nsigma > n_mad, default 5) kept for
# human readability.  They do NOT feed into the continuous total score — the
# score already captures the same information in a smoother way via all_nsigma.
#
# Best / margin
#   = (score_2nd − score_best) / score_2nd × 100 %
# ---------------------------------------------------------------------------

# Relative weights for each metric in the total score.  All equal by default.
# Increase a value to make that metric more decisive, e.g. {'clock (ns)': 2.0}.
_WEIGHTS = {}  # all metrics weighted equally


def _score(issues, all_nsigma):
    """
    Compute summary scores from the per-station nsigma dict.

    Returns
    -------
    total : float
        Weighted mean of per-metric mean n-sigmas.  Lower = better.
    n_flagged : int
        Number of stations with at least one flagged metric.
    by_metric_flagged : dict  metric -> int
        Number of flagged stations per metric.
    metric_scores : dict  metric -> float
        Mean n-sigma per metric (used in the comparison table).
    """
    metric_scores = {}
    for metric, ants in all_nsigma.items():
        if ants:
            metric_scores[metric] = float(np.mean(list(ants.values())))
    if metric_scores:
        weights = [_WEIGHTS.get(m, 1.0) for m in metric_scores]
        total = float(np.average(list(metric_scores.values()), weights=weights))
    else:
        total = 0.0
    flagged = set(issues.keys())
    by_metric_flagged = {}
    for an, flags in issues.items():
        for metric, *_ in flags:
            by_metric_flagged[metric] = by_metric_flagged.get(metric, 0) + 1
    return total, len(flagged), by_metric_flagged, metric_scores


def _add_scores(results):
    """Attach score, n_flagged, by_metric, and metric_scores to each result dict.
    Returns a new list sorted best-first (lowest weighted score first).
    """
    out = []
    for r in results:
        total, nants, by_metric_flagged, metric_scores = _score(r['issues'], r['all_nsigma'])
        out.append({**r, 'score': total, 'n_flagged': nants,
                    'by_metric': by_metric_flagged, 'metric_scores': metric_scores})
    return sorted(out, key=lambda x: x['score'])


def compare_results(results):
    """
    Given a list of result dicts (from assess_observation), print a comparison
    table and return the list sorted best-first (lowest weighted score first).

    Parameters
    ----------
    results : list of dict
        Each dict has keys: obs_id, calibrator, issues.

    Returns
    -------
    list of dict, sorted best-first, each with an added 'score' key.
    """
    scored = _add_scores(results)

    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    labels = [f"{r['obs_id']} ({r['calibrator']})" for r in scored]
    col_w = max(20, max(len(l) for l in labels))
    header = f"  {'Metric':<35}" + "".join(f" {l:>{col_w}}" for l in labels)
    print(header)
    print("  " + "-" * (35 + col_w * len(scored) + len(scored)))

    def row(name, vals):
        return f"  {name:<35}" + "".join(f" {v:>{col_w}}" for v in vals)

    print(row("Flagged antennas",              [r['n_flagged'] for r in scored]))
    print(row("Total score (wtd mean n-sigma)", [f"{r['score']:.3f}" for r in scored]))
    print()
    all_metrics = sorted({m for r in scored for m in r['metric_scores']})
    print(row("  Per-metric mean n-sigma", [""] * len(scored)))
    for m in all_metrics:
        w = _WEIGHTS.get(m, 1.0)
        label = f"    {m} (x{w:.0f})"
        print(row(label, [f"{r['metric_scores'].get(m, 0.0):.1f}" for r in scored]))
    print()
    all_flagged_metrics = sorted({m for r in scored for m in r['by_metric']})
    if all_flagged_metrics:
        print(row("  Flagged station count", [""] * len(scored)))
        for m in all_flagged_metrics:
            print(row(f"    {m}", [r['by_metric'].get(m, 0) for r in scored]))

    print()
    best = scored[0]
    if len(scored) > 1 and scored[0]['score'] < scored[1]['score']:
        margin = (scored[1]['score'] - scored[0]['score']) / max(scored[1]['score'], 1) * 100
        print(f"  ✓  Best: {best['obs_id']} ({best['calibrator']})  "
              f"— {margin:.0f}% lower weighted flag score than next best")
    elif len(scored) > 1:
        print("  All observations are equally good.")

    return scored


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_assessment(result, template, plot_dir):
    """
    For a single assessed observation, produce plots showing each station's
    observed value (dot) against the template median ± MAD (error bar).
    Flagged stations are highlighted in red.

    Parameters
    ----------
    result   : dict returned by assess_observation
    template : dict loaded from calibrator_medians_cache.pkl
    plot_dir : str  directory in which to write PNG files
    """
    import math
    os.makedirs(plot_dir, exist_ok=True)

    obs_id    = result['obs_id']
    calibrator = result['calibrator']
    issues    = result['issues']        # ant -> [(metric, obs, med, mad, nsigma), ...]

    flagged_per_metric = {}  # metric -> set of ant_names
    for an, flags in issues.items():
        for metric, *_ in flags:
            flagged_per_metric.setdefault(metric, set()).add(an)

    # ---- clock: two-panel plot (CS / RS + outlier CS) ----
    median_clock = template['median_clock']
    if calibrator in {k for d in median_clock.values() for k in d}:
        all_ants = sorted(an for an in median_clock if calibrator in median_clock[an])
        flagged_clk = flagged_per_metric.get('clock (ns)', set())

        # All observed clock values from the observed dict
        obs_clk = result.get('observed', {}).get('clock (ns)', {})

        # Identify outlier CS stations using template medians (same logic as build_template)
        cs_ants = [an for an in all_ants if an.startswith('CS')]
        cs_meds = {an: float(np.nanmedian(median_clock[an][calibrator][1])) for an in cs_ants}
        if cs_meds:
            bulk_cs_med = float(np.nanmedian(list(cs_meds.values())))
            bulk_cs_mad = float(np.nanmedian([abs(v - bulk_cs_med) for v in cs_meds.values()])) or 1.0
            cs_threshold = max(5.0 * bulk_cs_mad, 10.0)
            outlier_cs = {an for an, v in cs_meds.items() if abs(v - bulk_cs_med) > cs_threshold}
        else:
            outlier_cs = set()

        panels = [
            ('CS (typical)',     lambda an: an.startswith('CS') and an not in outlier_cs),
            ('RS + outlier CS',  lambda an: an.startswith('RS') or an in outlier_cs),
        ]

        fig, axes = plt.subplots(1, 2, figsize=(18, 5))
        for ax, (panel_label, keep) in zip(axes, panels):
            names = [an for an in all_ants if keep(an)]
            if not names:
                ax.set_visible(False)
                continue
            ref_meds, ref_mads = [], []
            for an in names:
                _, vals, mads = median_clock[an][calibrator]
                ref_meds.append(float(np.nanmedian(vals)))
                ref_mads.append(float(np.nanmedian(mads)))

            x = np.arange(len(names))
            ax.errorbar(x, ref_meds, yerr=ref_mads, fmt='none',
                        ecolor='steelblue', elinewidth=1.5, capsize=3, zorder=2,
                        label='Template median ± MAD')
            ax.scatter(x, ref_meds, s=20, color='steelblue', zorder=3)

            for xi, an in enumerate(names):
                if an in obs_clk:
                    color = 'red' if an in flagged_clk else 'orange'
                    ax.scatter(xi, obs_clk[an], s=60, color=color, zorder=4,
                               marker='D',
                               label='Observed (flagged)' if color == 'red' else 'Observed')

            ax.set_xticks(x)
            ax.set_xticklabels(names, rotation=45, ha='right', fontsize=7)
            for tick, an in zip(ax.get_xticklabels(), names):
                is_bad = an in flagged_clk or an in outlier_cs
                if is_bad:
                    tick.set_color('red')
                    tick.set_fontweight('bold')
            ax.set_ylabel('Clock (ns)')
            ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)
            ax.grid(True, axis='y', alpha=0.3)
            # Independent y-limits based on inter-station spread
            centre    = float(np.nanmedian(ref_meds))
            inter_mad = float(np.nanmedian([abs(m - centre) for m in ref_meds])) or 1.0
            obs_floor = max(float(np.nanmedian(ref_mads)), 1.0)
            half      = max(5.0 * inter_mad, 3.0 * obs_floor, 2.0)
            ax.set_ylim(centre - half, centre + half)
            handles, lbls = ax.get_legend_handles_labels()
            ax.legend(dict(zip(lbls, handles)).values(),
                      dict(zip(lbls, handles)).keys(), fontsize=8)
            ax.set_title(f'{panel_label}')

        fig.suptitle(f'{obs_id} ({calibrator}) — Clock (ns)  '
                     f'[CS and RS scaled independently; outlier CS in red]', fontsize=11)
        fig.tight_layout()
        out = os.path.join(plot_dir, f'{obs_id}_clock_ns.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  wrote {out}')

    # ---- polalign and flagging scalar plots ----
    SCALAR_METRICS = [
        ('polalign (ns)', template['median_polalign'], 'Polalign YY delay (ns)'),
    ]
    final_data = template['median_flagging'].get('final', {})
    if final_data:
        SCALAR_METRICS.append(('flag_final', final_data, 'Final flagging (%)'))

    obs_all = result.get('observed', {})

    for metric, median_data, ylabel in SCALAR_METRICS:
        if calibrator not in {k for d in median_data.values() for k in d}:
            continue
        ant_names = sorted(an for an in median_data if calibrator in median_data[an])
        if not ant_names:
            continue
        ref_meds, ref_mads = [], []
        for an in ant_names:
            _, vals, mads = median_data[an][calibrator]
            ref_meds.append(np.nanmedian(vals))
            ref_mads.append(np.nanmedian(mads))

        metric_obs = obs_all.get(metric, {})
        flagged_set = flagged_per_metric.get(metric, set())

        fig, ax = plt.subplots(figsize=(max(10, len(ant_names) * 0.35), 5))
        x = np.arange(len(ant_names))
        # Template
        ax.errorbar(x, ref_meds, yerr=ref_mads, fmt='none',
                    ecolor='steelblue', elinewidth=1.5, capsize=3, zorder=2,
                    label='Template median ± MAD')
        ax.scatter(x, ref_meds, s=20, color='steelblue', zorder=3)
        # All observed values
        obs_x, obs_y, obs_colors = [], [], []
        for xi, an in enumerate(ant_names):
            if an not in metric_obs: continue
            obs_x.append(xi)
            obs_y.append(metric_obs[an])
            obs_colors.append('red' if an in flagged_set else 'orange')
        if obs_x:
            # Normal points first, then outliers on top
            norm_x  = [xi for xi, c in zip(obs_x, obs_colors) if c != 'red']
            norm_y  = [y  for y,  c in zip(obs_y, obs_colors) if c != 'red']
            flag_x  = [xi for xi, c in zip(obs_x, obs_colors) if c == 'red']
            flag_y  = [y  for y,  c in zip(obs_y, obs_colors) if c == 'red']
            if norm_x:
                ax.scatter(norm_x, norm_y, s=40, color='orange', zorder=4,
                           marker='D', label='Observed')
            if flag_x:
                ax.scatter(flag_x, flag_y, s=60, color='red', zorder=5,
                           marker='D', label='Observed (outlier)')

        ax.set_xticks(x)
        ax.set_xticklabels(ant_names, rotation=45, ha='right', fontsize=7)
        for tick, an in zip(ax.get_xticklabels(), ant_names):
            if an in flagged_set:
                tick.set_color('red')
                tick.set_fontweight('bold')
        ax.set_ylabel(ylabel)
        ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)
        ax.grid(True, axis='y', alpha=0.3)
        handles, lbls = ax.get_legend_handles_labels()
        ax.legend(dict(zip(lbls, handles)).values(),
                  dict(zip(lbls, handles)).keys(), fontsize=8)
        ax.set_title(f'{obs_id} ({calibrator}) — {ylabel}')
        fig.tight_layout()
        safe_metric = metric.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
        out = os.path.join(plot_dir, f'{obs_id}_{safe_metric}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  wrote {out}')

    # ---- bandpass: spectrum plot per station ----
    median_bandpass = template['median_bandpass']
    if calibrator in {k for d in median_bandpass.values() for k in d}:
        ant_names = sorted(an for an in median_bandpass if calibrator in median_bandpass[an])
        flagged_bp = flagged_per_metric.get('bandpass_residual', set())
        # load observed spectra from cal_solutions.h5
        h5_path = result.get('h5_path')
        obs_spectra = {}
        if h5_path and os.path.exists(h5_path):
            try:
                with h5py.File(h5_path, 'r') as h:
                    bp_ants = [a.decode() if isinstance(a, bytes) else a
                               for a in h['calibrator/bandpass/ant'][...]]
                    val = h['calibrator/bandpass/val'][...]   # (time, freq, ant, pol)
                for ai, an in enumerate(bp_ants):
                    if not _is_station_of_interest(an): continue
                    obs_spectra[an] = np.nanmedian(val[:, :, ai, :], axis=(0, 2))
            except Exception as e:
                print(f'  Warning reading bandpass for plot: {e}')

        ncols = 4
        nrows = math.ceil(len(ant_names) / ncols)
        fig, axes = plt.subplots(nrows, ncols,
                                  figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
        for idx, an in enumerate(ant_names):
            ax = axes[idx // ncols][idx % ncols]
            entry = median_bandpass[an][calibrator]
            freq, ref_med, ref_mad = entry[0], entry[1], entry[2]
            fmhz = freq / 1e6
            ax.fill_between(fmhz, ref_med - ref_mad, ref_med + ref_mad,
                            color='steelblue', alpha=0.3, label='Template ± MAD')
            ax.plot(fmhz, ref_med, color='steelblue', linewidth=1)
            if an in obs_spectra and len(obs_spectra[an]) == len(ref_med):
                color = 'red' if an in flagged_bp else 'orange'
                ax.plot(fmhz, obs_spectra[an], color=color, linewidth=1,
                        label='Observed' + (' (flagged)' if an in flagged_bp else ''))
            title_color = 'red' if an in flagged_bp else 'black'
            ax.set_title(an, fontsize=8, color=title_color,
                         fontweight='bold' if an in flagged_bp else 'normal')
            ax.tick_params(labelsize=6)
            ax.grid(True, alpha=0.2)
        for idx in range(len(ant_names), nrows * ncols):
            axes[idx // ncols][idx % ncols].set_visible(False)
        fig.suptitle(f'{obs_id} ({calibrator}) — Bandpass vs template', fontsize=12)
        fig.tight_layout()
        out = os.path.join(plot_dir, f'{obs_id}_bandpass.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  wrote {out}')


# ---------------------------------------------------------------------------
# Top-level importable function
# ---------------------------------------------------------------------------

def _load_template(h5_path):
    """Load the HDF5 template cache written by build_template.py.

    Reconstructs the same nested-dict structure previously produced by pickle,
    so that all downstream code is unchanged.
    """
    with h5py.File(h5_path, 'r') as f:
        obs_to_calibrator = json.loads(f.attrs['obs_to_calibrator'])

        # median_bandpass: ant -> cal -> (freq, med, mad, count)
        median_bandpass = {}
        for ant in f['median_bandpass']:
            median_bandpass[ant] = {}
            for cal in f[f'median_bandpass/{ant}']:
                cg = f[f'median_bandpass/{ant}/{cal}']
                median_bandpass[ant][cal] = (
                    cg['freq'][...],
                    cg['med'][...],
                    cg['mad'][...],
                    int(cg.attrs['count']),
                )

        # median_polalign / median_clock: ant -> cal -> (x, vals, mads)
        median_polalign = {}
        for ant in f['median_polalign']:
            median_polalign[ant] = {}
            for cal in f[f'median_polalign/{ant}']:
                cg = f[f'median_polalign/{ant}/{cal}']
                median_polalign[ant][cal] = (cg['x'][...], cg['vals'][...], cg['mads'][...])

        median_clock = {}
        for ant in f['median_clock']:
            median_clock[ant] = {}
            for cal in f[f'median_clock/{ant}']:
                cg = f[f'median_clock/{ant}/{cal}']
                median_clock[ant][cal] = (cg['x'][...], cg['vals'][...], cg['mads'][...])

        # median_flagging: {'final': {ant -> cal -> (x, vals, mads)}}
        median_flagging = {}
        for stage in f['median_flagging']:
            median_flagging[stage] = {}
            for ant in f[f'median_flagging/{stage}']:
                median_flagging[stage][ant] = {}
                for cal in f[f'median_flagging/{stage}/{ant}']:
                    cg = f[f'median_flagging/{stage}/{ant}/{cal}']
                    median_flagging[stage][ant][cal] = (
                        cg['x'][...], cg['vals'][...], cg['mads'][...])

    return dict(
        obs_to_calibrator=obs_to_calibrator,
        median_bandpass=median_bandpass,
        median_polalign=median_polalign,
        median_clock=median_clock,
        median_flagging=median_flagging,
    )


def assess_and_compare(template_path, h5_paths, n_mad=5.0, verbose=True, plot_dir=None):
    """
    Assess one or more cal_solutions.h5 files against a pre-built template
    and return a ranked list.

    Parameters
    ----------
    template_path : str
        Path to the calibrator_medians_cache.h5 produced by build_template.py.
    h5_paths : list of str
        One or more paths to cal_solutions.h5 files.
    n_mad : float
        Outlier threshold in units of MAD (default 5.0).
    verbose : bool
        Print per-antenna reports and comparison table.
    plot_dir : str or None
        If given, write assessment plots (observed vs template) into this directory.

    Returns
    -------
    list of dict, sorted best-first. Each dict has:
        obs_id, calibrator, issues, score, n_flagged, by_metric.
    """
    template = _load_template(h5_path=template_path)

    results = [assess_observation(h5, template, n_mad=n_mad, verbose=verbose)
               for h5 in h5_paths]

    if plot_dir:
        for r in results:
            obs_out = os.path.join(plot_dir, r['obs_id'])
            print(f"\nWriting plots for {r['obs_id']} to {obs_out}/")
            plot_assessment(r, template, obs_out)

    return compare_results(results)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Assess LOFAR calibrator solutions against a pre-built template.')
    parser.add_argument('--template', required=True,
                        help='Path to calibrator_medians_cache.h5 produced by build_template.py.')
    parser.add_argument('--h5', nargs='+', required=True,
                        help='One or more cal_solutions.h5 files to assess and compare.')
    parser.add_argument('--n-mad', type=float, default=5.0,
                        help='Outlier threshold in units of MAD (default: 5.0).')
    parser.add_argument('--plot-dir', default=None,
                        help='If given, write assessment plots into this directory '
                             '(one sub-folder per observation).')
    args = parser.parse_args()

    assess_and_compare(args.template, args.h5, n_mad=args.n_mad,
                       plot_dir=args.plot_dir)


if __name__ == "__main__":
    main()
