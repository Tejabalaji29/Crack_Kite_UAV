# Cranked-kite tailless UAV — OpenVSP model

Parametric conceptual design model of a 900 mm span, 2.0 kg cranked-kite
flying wing for short-range reconnaissance. Built for OpenVSP's Python API so
the geometry, the balance and the trade studies are all driven from one set of
numbers instead of GUI clicks.

## Files

| File | What it does |
|---|---|
| `geometry.py` | Single source of truth. All dimensions plus derived S, AR, MAC, CG window, control-surface area. Run it on its own for the full geometry report. |
| `build_model.py` | Constructs the wing in OpenVSP and writes `cranked_kite_v1.vsp3`. |
| `run_analysis.py` | VSPAERO degen-geom + alpha sweep + stability. Writes a polar CSV and a readable summary. |
| `trade_study.py` | Automated parameter sweeps: crank position, outboard sweep, washout. |

## Setup

OpenVSP's Python bindings must be importable:

```bash
python3 -c "import openvsp"      # must not error
```

If it errors, add your OpenVSP installation's `python/` directory to
`PYTHONPATH`, or pip-install the wheel that ships inside it.

Then put the reflex airfoil next to these scripts:

```bash
# MH 60 from Airfoil Tools, saved as a Selig/Lednicer .dat
curl -o mh60.dat <airfoil source>
```

**The reflex airfoil is not optional.** `build_model.py` will fall back to a
NACA four-series so the model still builds, but a four-series has a strongly
negative pitching moment and a tailless wing cannot trim with it. Any Cm or
neutral-point result from the fallback is meaningless.

## Run order

```bash
python3 geometry.py                    # check the numbers first
python3 build_model.py                 # -> cranked_kite_v1.vsp3
python3 run_analysis.py                # loiter condition
python3 run_analysis.py --condition dash
python3 trade_study.py crank           # headline trade study
python3 trade_study.py sweep
python3 trade_study.py washout
```

Add `--dry-run` to any trade study to check the geometry without spending
solver time. `--dry-run` needs no OpenVSP install.

## Folder contents

```
cranked_kite_uav/
  README.md                     this file
  geometry.py                   configuration + derived properties (run first)
  build_model.py                builds the OpenVSP model, writes the .vsp3
  run_analysis.py               VSPAERO sweep, polars, neutral point
  trade_study.py                crank / sweep / washout studies
  docs/
    UAV_Project_Report.pdf      full project report, 19 pages
    UAV_Project_Report.docx     same report, editable
  figures/
    three_view_wireframe.svg    dimensioned three-view (.png alongside)
    structural_layout.svg       internal structure, plan view (.png alongside)
  outputs/
    study_crank.csv             crank trade study, geometry columns only
```

Files not included, because they must be generated or downloaded locally:

- `mh60.dat` — download from airfoiltools.com, place beside the scripts
- `cranked_kite_v1.vsp3` — produced by `build_model.py`
- `*_polar.csv`, `*_summary.txt` — produced by `run_analysis.py`

## Baseline configuration

| Parameter | Value |
|---|---|
| MTOW | 2.00 kg |
| Span | 900 mm |
| Reference area | 0.2890 m² |
| Aspect ratio | 2.80 |
| MAC | 362.6 mm |
| MAC quarter-chord | 251.2 mm aft of nose |
| Root / crank / tip chord | 520 / 300 / 120 mm |
| Inboard / outboard LE sweep | 42° / 28° |
| Mean quarter-chord sweep | 29.2° |
| Washout at tip | −4° |
| Wing loading | 6.92 kg/m² |
| CG target (6 % SM) | 240 mm aft of nose |
| Acceptable CG window | 233–244 mm |
| Elevon area | 6.94 % of S_ref |
| Loiter / dash speed | 17 / 28 m/s |
| Reynolds at loiter | 4.2 × 10⁵ |

## Trade studies

**Crank position** is the headline study. Span, root chord, tip chord and
reference area are all held fixed; the chord at the crank is solved from the
area constraint, so each variant is a genuinely different break severity
rather than a relabelled straight-taper wing. Crank chord runs 342 mm down to
242 mm as the break moves from 45 % to 70 % of half-span, with mean
quarter-chord sweep rising 27.5° → 31.3°. Plot L/D max and dCm/dCL against
crank position.

**Outboard sweep** 20°–40°: shows the directional-stability versus drag trade.
Plot Cnβ and L/D.

**Washout** 0° to −7°: plot spanwise lift distribution and L/D. This is the
tip-stall protection study.

## Known limitations to state in any report

- Vortex lattice assumes attached inviscid flow. No stall prediction. Treat
  CL beyond roughly 0.9 as fiction.
- At Re ≈ 4 × 10⁵ the real sections suffer laminar separation bubbles that VLM
  cannot see. Cross-check the airfoil polars in XFLR5 at this Reynolds number.
- The thick blended root is poorly represented by VLM. Run the Parasite Drag
  tool for the CD0 buildup and add it to VSPAERO's induced drag.
- The pusher propeller runs in the wing wake. VSPAERO's actuator disc does not
  capture that interaction well; treat the propulsive efficiency as an
  assumption, not a result.

## Verification performed

The closed-form panel formulae in `geometry.py` for reference area, MAC and
MAC spanwise station were checked against brute-force numerical integration
over two million spanwise steps. All three agree to better than 1 × 10⁻⁵.
