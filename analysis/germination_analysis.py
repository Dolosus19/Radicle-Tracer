#!/usr/bin/env python3
"""
Germination analysis for Radicle Tracer per-seed CSVs.

Reads one or more per-seed CSV files exported from Radicle Tracer (index.html),
combines every time point, and writes:

  per_dish_by_time.csv    germination %, emergence %, mean radicle length per dish and time
  per_group_by_time.csv   group means and SDs per time point, with Welch t-tests (2 groups)
  per_dish_t50.csv        time to 50% germination per dish (linear interpolation)
  icc_by_time.csv         intraclass correlation (seeds within dishes), per group and time
  survival_model.txt      discrete-time survival model (time to germination, all images)
  germination_curve.png   mean germination % over time by group
  radicle_length.png      mean radicle length over time by group

Usage
-----
  python germination_analysis.py data/*.csv --out results
  python germination_analysis.py data/*.csv --design design.csv --cutoff 2 --out results

Design file (optional, recommended)
-----------------------------------
A CSV with columns `dish,group` and optionally `run`, one row per dish, e.g.

  dish,group,run
  A1,covered,1
  B1,open,1

Without it, the group is the leading letters of the dish name (A1 -> A) and all
dishes belong to one run. For repeated experiments, give each run its own dish
names (e.g. R1-A1) or list the run in the design file.

Time is read from the `time` column (e.g. "T72h" -> 72 hours).
"""
import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load(paths):
    files = []
    for p in paths:
        files.extend(sorted(glob.glob(p)) or [p])
    if not files:
        sys.exit("No input files found.")
    frames = []
    for f in files:
        d = pd.read_csv(f)
        missing = {"dish", "time", "seed", "status", "length_mm"} - set(d.columns)
        if missing:
            sys.exit(f"{f}: missing columns {sorted(missing)}. Use the per-seed CSV from Radicle Tracer.")
        d["source_file"] = Path(f).name
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    d["hours"] = d["time"].astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0].astype(float)
    if d["hours"].isna().any():
        bad = d.loc[d["hours"].isna(), "time"].unique()[:5]
        sys.exit(f"Could not read a time from: {list(bad)}. Name images like A1_T72h.jpg.")
    d["length_mm"] = pd.to_numeric(d["length_mm"], errors="coerce")
    d.loc[d["status"].astype(str).str.lower() != "traced", "length_mm"] = np.nan
    dup = d.duplicated(["dish", "hours", "seed"], keep="last")
    if dup.any():
        print(f"Note: {dup.sum()} duplicate seed rows found; keeping the last one of each.")
        d = d[~dup]
    unscored = (d["status"].astype(str).str.lower() == "unscored").sum()
    if unscored:
        print(f"Warning: {unscored} seeds are unscored and count as not germinated.")
    return d


def attach_design(d, design_path):
    if design_path:
        des = pd.read_csv(design_path, dtype=str)
        if not {"dish", "group"} <= set(des.columns):
            sys.exit("Design file needs columns dish,group (and optionally run).")
        d = d.merge(des, on="dish", how="left", validate="many_to_one")
        if d["group"].isna().any():
            sys.exit(f"Dishes missing from the design file: {sorted(d.loc[d.group.isna(), 'dish'].unique())}")
        if "run" not in d.columns:
            d["run"] = "1"
    else:
        d["group"] = d["dish"].astype(str).str.extract(r"^([A-Za-z]+)")[0]
        d["run"] = "1"
    d["run"] = d["run"].astype(str)
    d["unit"] = d["run"] + "|" + d["dish"].astype(str)
    return d


def per_dish(d, cutoff, seeds_per_dish):
    def f(x):
        L = x["length_mm"]
        g = L[L > cutoff]
        n = seeds_per_dish or len(x)
        return pd.Series({
            "seeds": n,
            "germinated": int((L > cutoff).sum()),
            "germination_pct": 100 * (L > cutoff).sum() / n,
            "emergence_pct": 100 * L.notna().sum() / n,
            "mean_length_mm": g.mean() if len(g) else np.nan,
            "sd_length_mm": g.std() if len(g) > 1 else np.nan,
            "total_length_mm": L.sum(),
        })
    return (d.groupby(["run", "group", "dish", "unit", "hours"]).apply(f, include_groups=False)
              .reset_index().sort_values(["hours", "run", "group", "dish"]))


def welch(a, b):
    from scipy.stats import ttest_ind
    a, b = a.dropna(), b.dropna()
    if len(a) < 2 or len(b) < 2 or (a.std() == 0 and b.std() == 0):
        return np.nan
    return ttest_ind(a, b, equal_var=False).pvalue


def per_group(pd_dish):
    rows = []
    groups = sorted(pd_dish["group"].unique())
    for h, s in pd_dish.groupby("hours"):
        for g in groups:
            x = s[s.group == g]
            rows.append({"hours": h, "group": g, "dishes": len(x),
                         "germination_pct": x.germination_pct.mean(), "germination_sd": x.germination_pct.std(),
                         "emergence_pct": x.emergence_pct.mean(),
                         "mean_length_mm": x.mean_length_mm.mean(), "length_sd": x.mean_length_mm.std()})
        if len(groups) == 2:
            a, b = (s[s.group == g] for g in groups)
            rows[-1]["p_germination"] = welch(a.germination_pct, b.germination_pct)
            rows[-1]["p_length"] = welch(a.mean_length_mm, b.mean_length_mm)
    return pd.DataFrame(rows)


def monotone_counts(pd_dish):
    """Germinated count per dish and time, forced to never decrease (a seed can't un-germinate)."""
    c = pd_dish.pivot_table(index=["run", "group", "dish", "unit", "seeds"], columns="hours",
                            values="germinated", aggfunc="first").sort_index(axis=1)
    raw = c.copy()
    c = c.cummax(axis=1)
    fixed = (raw.fillna(-1).values < c.fillna(-1).values).any(axis=1).sum()
    if fixed:
        print(f"Note: in {fixed} dish(es) the germinated count went down between sessions; "
              "later counts were raised to the earlier maximum for the survival analysis.")
    return c


def t50(counts):
    hrs = list(counts.columns)
    out = []
    for idx, r in counts.iterrows():
        n = idx[-1]
        half = n / 2
        val = np.nan
        for i in range(1, len(hrs)):
            if r.iloc[i] >= half and r.iloc[i - 1] < half:
                y0, y1 = r.iloc[i - 1], r.iloc[i]
                val = hrs[i - 1] + (half - y0) / (y1 - y0) * (hrs[i] - hrs[i - 1])
                break
        out.append(dict(zip(["run", "group", "dish", "unit", "seeds"], idx), t50_hours=val,
                        final_germination_pct=100 * r.iloc[-1] / n))
    return pd.DataFrame(out)


def survival_model(counts):
    """Discrete-time proportional-hazards model (complementary log-log), all imaging intervals."""
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
    except ImportError:
        return "statsmodels is not installed; skipped. pip install statsmodels"
    hrs = list(counts.columns)
    rows = []
    for idx, r in counts.iterrows():
        run, group, dish, unit, n = idx
        prev = 0
        for i in range(1, len(hrs)):
            risk = n - prev
            if risk <= 0:
                break
            new = int(r.iloc[i] - prev)
            rows.append(dict(run=run, group=group, unit=unit, interval=f"{hrs[i-1]:g}-{hrs[i]:g}",
                             new=new, stay=risk - new))
            prev = r.iloc[i]
    dt = pd.DataFrame(rows)
    # drop intervals in which no seed germinated in any dish (they carry no information)
    active = dt.groupby("interval")["new"].transform("sum") > 0
    dt = dt[active]
    groups = sorted(dt.group.unique())
    formula = "new + stay ~ C(interval) + C(group, Treatment(reference='%s'))" % groups[0]
    if dt.run.nunique() > 1:
        formula += " + C(run)"
    model = smf.glm(formula, data=dt, family=sm.families.Binomial(link=sm.families.links.CLogLog()))
    fit = model.fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dt.unit)[0]})
    lines = ["Discrete-time survival model: time to germination (radicle > cutoff)",
             "Link: complementary log-log (proportional hazards). SEs clustered by dish.",
             f"Reference group: {groups[0]}. Intervals with no germination in any dish were dropped.",
             f"Run included as a block: {'yes' if dt.run.nunique() > 1 else 'no (one run)'}", ""]
    for name in fit.params.index:
        if name.startswith("C(group"):
            g = re.search(r"\[T\.(.+)\]", name).group(1)
            hr = np.exp(fit.params[name])
            lo, hi = np.exp(fit.conf_int().loc[name])
            lines.append(f"Hazard ratio, {g} vs {groups[0]}: {hr:.2f} (95% CI {lo:.2f}-{hi:.2f}), p = {fit.pvalues[name]:.3g}")
    lines += ["", "HR > 1 means faster germination than the reference group.", "", str(fit.summary())]
    return "\n".join(lines)


def icc1(df, col):
    d = df.dropna(subset=[col])
    k = d.groupby("unit")[col].size()
    if len(k) < 2:
        return np.nan
    N, a = k.sum(), len(k)
    gm = d[col].mean()
    means = d.groupby("unit")[col].mean()
    msb = (k * (means - gm) ** 2).sum() / (a - 1)
    msw = ((d[col] - d.groupby("unit")[col].transform("mean")) ** 2).sum() / (N - a)
    k0 = (N - (k ** 2).sum() / N) / (a - 1)
    denom = msb + (k0 - 1) * msw
    return (msb - msw) / denom if denom > 0 else np.nan


def icc_table(d, cutoff):
    d = d.copy()
    d["germ"] = (d.length_mm > cutoff).astype(float)
    d["len_germ"] = d.length_mm.where(d.length_mm > cutoff)
    rows = []
    for (h, g), s in d.groupby(["hours", "group"]):
        rows.append({"hours": h, "group": g, "icc_germination": icc1(s, "germ"), "icc_length": icc1(s, "len_germ")})
    return pd.DataFrame(rows)


def plots(grp, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plots.")
        return
    for col, sd, ylabel, fname, title in [
        ("germination_pct", "germination_sd", "Germinated (%)", "germination_curve.png", "Germination over time"),
        ("mean_length_mm", "length_sd", "Mean radicle length (mm)", "radicle_length.png", "Radicle length over time"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for g, s in grp.groupby("group"):
            s = s.dropna(subset=[col])
            ax.errorbar(s.hours, s[col], yerr=s[sd], marker="o", ms=4, capsize=2, label=str(g))
        ax.set_xlabel("Hours after sowing")
        ax.set_ylabel(ylabel)
        ax.set_title(title + " (mean ± SD of dishes)")
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(out / fname, dpi=150)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", help="per-seed CSV file(s) from Radicle Tracer; wildcards allowed")
    ap.add_argument("--design", help="CSV with dish,group[,run]")
    ap.add_argument("--cutoff", type=float, default=2.0, help="germination cutoff in mm (default 2)")
    ap.add_argument("--seeds-per-dish", type=int, default=None,
                    help="denominator for %% (default: number of seed rows per dish and time)")
    ap.add_argument("--out", default="results", help="output folder (default: results)")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    d = attach_design(load(args.csv), args.design)
    print(f"{len(d)} seed rows · {d.unit.nunique()} dishes · {d.hours.nunique()} time points · "
          f"groups: {', '.join(sorted(d.group.unique()))} · runs: {d.run.nunique()}")

    dish = per_dish(d, args.cutoff, args.seeds_per_dish)
    dish.drop(columns="unit").to_csv(out / "per_dish_by_time.csv", index=False)
    grp = per_group(dish)
    grp.to_csv(out / "per_group_by_time.csv", index=False)
    counts = monotone_counts(dish)
    t = t50(counts)
    t.drop(columns="unit").to_csv(out / "per_dish_t50.csv", index=False)
    icc_table(d, args.cutoff).to_csv(out / "icc_by_time.csv", index=False)
    report = survival_model(counts)
    (out / "survival_model.txt").write_text(report)
    plots(grp, out)

    print("\nT50 (hours to 50% germination), mean ± SD of dishes:")
    for g, s in t.groupby("group"):
        print(f"  {g}: {s.t50_hours.mean():.1f} ± {s.t50_hours.std():.1f} (n = {s.t50_hours.notna().sum()} dishes)")
    print("\n" + report.split("\n\n")[1].strip() if "Hazard ratio" in report else report.splitlines()[0])
    print(f"\nResults written to {out.resolve()}")


if __name__ == "__main__":
    main()
