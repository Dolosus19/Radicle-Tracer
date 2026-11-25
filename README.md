# Radicle Tracer

A browser tool for measuring seed germination from dish photos, plus a Python script that analyses the results over time.

- **Trace** each radicle by clicking along it, from where it leaves the seed coat to the tip. Lengths are calibrated to millimetres from a ruler in the photo.
- **Score** each seed as germinated (radicle over a cutoff, 2 mm by default), below the cutoff, or no radicle, and flag uncertain calls.
- **Export** per-seed and per-dish CSVs, then run `analysis/germination_analysis.py` to get germination curves, time to germination, radicle length and between-group statistics across all time points.

The tool is a single HTML file with no dependencies. Images never leave your computer.

## Quick start

### 1. Open the tool

Use either of these:

- **Locally:** download `index.html` and open it in Chrome, Edge, Firefox or Safari.
- **Online:** enable GitHub Pages for this repository (Settings → Pages → Deploy from branch → `main`, folder `/`). The tool is then at `https://<your-username>.github.io/radicle-tracer/`.

### 2. Name your images

Use `<dish>_T<hours>h.jpg`, for example `A1_T72h.jpg` or `R2-B4_T96h.jpg`. The tool reads the dish and time from the file name, and the analysis script needs both.

### 3. Trace

1. **Open images:** drop a session's photos onto the left panel.
2. **Calibrate:** click two ruler marks a known distance apart (default 80 mm). To reuse one calibration for every image in a session, click *Use this calibration for uncalibrated images*. Calibrating each image on its own ruler is more accurate.
3. **Trace each seed:** click where the radicle leaves the seed coat, click along the radicle, then press **Enter** at the tip. Press **0** if the seed has no radicle, and **U** to flag an uncertain call.
4. **Export:** click *Download per-seed CSV*. Save one per imaging session, or one for everything.

| Key | Action |
| --- | --- |
| Click | Add a point to the current trace (or a calibration mark) |
| Drag / Scroll | Pan / zoom |
| Enter | Finish the seed and go to the next |
| Backspace | Undo the last point |
| 0 | No radicle, next seed |
| U | Toggle uncertain |
| C / T | Calibrate mode / trace mode |
| F | Fit the image to the window |
| 1–9 | Jump to a seed |

### 4. Save your work

Traces autosave in the browser (`localStorage`) and are matched to images by file name and size. Clearing browser data deletes them, so also **Save session file** after each session. To resume, open the same images, then **Open session file**.

*Clear image list* empties the queue so you can start the next time point. It keeps saved traces. *Delete saved traces…* removes them all from the browser.

### 5. Analyse

```bash
pip install -r analysis/requirements.txt
python analysis/germination_analysis.py "data/*.csv" --design design.csv --out results
```

Try it on the included trial data:

```bash
python analysis/germination_analysis.py examples/trial-2026-11/per_seed_all_timepoints.csv \
  --design examples/trial-2026-11/design.csv --out results
```

## Design file

The design file maps each dish to its group and, for repeated experiments, its run:

```csv
dish,group,run
R1-A1,treated,1
R1-B1,control,1
R2-A1,treated,2
```

Without a design file, the group is the leading letters of the dish name (`A1` → `A`) and everything counts as one run. Give dishes unique names across runs, for example with an `R1-` prefix.

## What the analysis produces

| File | Contents |
| --- | --- |
| `per_dish_by_time.csv` | Germination %, emergence % (any radicle), mean and total radicle length for each dish at each time |
| `per_group_by_time.csv` | Group means and SDs at each time, with Welch t-tests on dish means (two groups) |
| `per_dish_t50.csv` | Time to 50% germination for each dish (interpolated) and final germination % |
| `icc_by_time.csv` | Intraclass correlation of seeds within dishes, for germination and length |
| `survival_model.txt` | Time-to-germination model using every image (below) |
| `germination_curve.png`, `radicle_length.png` | Group curves over time, mean ± SD of dishes |

### Statistical approach

- **The dish is the unit.** Seeds in one dish share conditions, so per-group statistics are means of dish values.
- **Time to germination uses all images.** Germinated counts per dish are converted to "new germinations per interval among seeds not yet germinated". These are fitted with a discrete-time proportional-hazards model (binomial GLM, complementary log-log link), with imaging interval and run as factors and standard errors clustered by dish. A hazard ratio above 1 means faster germination than the reference group. Seeds that never germinate count as censored, so differences in final % are included.
- **A seed can't un-germinate.** If a dish's germinated count drops between sessions (a borderline call flipping back), later counts are raised to the earlier maximum, and the script reports how often this happened. Seed numbers don't need to match across sessions.
- **Welch t-tests at each time point** are descriptive and not corrected for multiple comparisons. Pre-register one primary analysis (for example, the hazard ratio) before running a confirmatory experiment.

## CSV formats (exported by the tool)

**Per-seed** (input to the analysis): `image, dish, time, seed, status, length_mm, germinated, uncertain, n_points, note`. `status` is `traced`, `none` or `unscored`. `length_mm` is the traced path length.

**Summary:** one row per image with germination %, mean, SD, min and max length of germinated radicles, number flagged uncertain, and calibration (px/mm and its source).

## Tips from the trial run

- A ruler in every photo, at a fixed camera height, keeps calibration within a few percent.
- Short radicles near the cutoff are the least reliable calls. Flag them with **U** and recheck them, or re-trace a few dishes on another day to measure tracing error.
- Roots longer than the dish (about 90 mm after 6–7 days for *Cucurbita pepo*) are hard to trace accurately. If root length matters, stop measuring length by about 96 h.
- For germination timing, image often around the steep part of the curve (every 6 h from about 42 to 96 h for squash).
- Weigh dishes at each session if moisture could differ between treatments.

See [`docs/protocol.md`](docs/protocol.md) for the imaging and scoring protocol used in the example data.

## Repository layout

```
index.html                         the tracing tool
analysis/germination_analysis.py   analysis script
analysis/requirements.txt
docs/protocol.md                   imaging and scoring protocol
examples/trial-2026-11/            trial data (12 dishes × 15 time points) and design file
```

## License

MIT. See [LICENSE](LICENSE).
