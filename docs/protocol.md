# Imaging and scoring protocol

This is the protocol used for the November 2026 trial run (`examples/trial-2026-11`). Adapt it to your species and setup.

## Setup

- 90 mm Petri dishes, 10 *Cucurbita pepo* seeds per dish in a loose 3-4-3 layout.
- Each dish labelled with marker (e.g. `A1`). For repeated experiments, include the run in the label or file name (`R2-A1`).
- Dish positions moved at random after each imaging session.

Trial-run setups (the two groups in the example data):

| | Group A | Group B |
| --- | --- | --- |
| Paper towels | 2, seeds sandwiched between them | 1, seeds on top |
| Filtered water | 6 mL, once at the start | 3.5 mL, once at the start |

## Imaging

- Photograph each dish from directly above at a fixed height, with a ruler in frame in the same plane as the seeds.
- Every 12 h from 0 to 168 h in the trial. For germination timing, every 6 h from about 42 to 96 h is more informative.
- Name files `<dish>_T<hours>h.jpg`.

## Scoring

- **Germinated:** radicle longer than 2 mm.
- **Length:** measured along the radicle's curve, from where it leaves the seed coat to the tip.
- **Emerged:** any visible radicle, including under 2 mm.
- **No radicle:** press `0`.
- **Uncertain:** press `U` (for example, a tip hidden under the towel, or a radicle you can't tell from the seed coat). Add a note.
- **Calibration:** on the ruler in each image, or once per session applied to all images from that session.

## Known issues from the trial

- Seed numbers were not kept consistent between sessions. The dish-level analysis doesn't need them. If per-seed growth matters, number positions on a template under each dish.
- Sub-2 mm radicles sometimes flipped between "radicle" and "none" between sessions.
- Long, curved roots measured up to several mm shorter on later days (tracing error), and roots beyond about 90 mm outgrew the dish.
- Sandwiched (two-towel) roots can grow under the top towel, which may lead to undercounting.
