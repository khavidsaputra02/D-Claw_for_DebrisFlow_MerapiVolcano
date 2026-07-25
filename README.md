# D-Claw Lahar Merapi Volcano — Cold Lahar (Lahar Dingin) Simulation Pipeline

Simulation and post-processing pipeline for **cold lahar (lahar dingin) flow modeling** on the southern flank of **Mount Merapi**, Indonesia, built on **D-Claw** (the Clawpack/GeoClaw debris-flow extension). The pipeline takes a DEM of the study area, defines a lahar source region from a shapefile, runs a two-phase (solid/fluid) shallow granular-fluid simulation, and produces runout maps, flood-potential maps, sensitivity statistics, and animated visualizations.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Modeling Approach](#modeling-approach)
- [Requirements](#requirements)
- [Installation](#installation)
  - [1. System dependencies](#1-system-dependencies)
  - [2. Clawpack / D-Claw](#2-clawpack--d-claw)
  - [3. Python environment](#3-python-environment)
  - [4. Place this pipeline in a D-Claw example directory](#4-place-this-pipeline-in-a-d-claw-example-directory)
- [Input Data](#input-data)
- [Usage](#usage)
  - [Step 1 — Prepare topography and source inputs (`setinput.py`)](#step-1--prepare-topography-and-source-inputs-setinputpy)
  - [Step 2 — Configure and run the simulation (`setrun.py`)](#step-2--configure-and-run-the-simulation-setrunpy)
  - [Step 3 — Visualize results (`setplot_Selatan.py`)](#step-3--visualize-results-setplotpy)
  - [Step 4 — Extract runout and sensitivity data (`Extract Runout dan Shp file.py`)](#step-4--extract-runout-and-sensitivity-data-extract-runout-dan-shp-filepy)
  - [Step 5 — Generate the flood-potential map (`Peta potensi banjir lahar.py`)](#step-5--generate-the-flood-potential-map-peta-potensi-banjir-laharpy)
  - [Optional — DEM erosion/analysis utility (`DEM Analisis.py`)](#optional--dem-erosionanalysis-utility-dem-analisispy)
- [Key Parameters](#key-parameters)
- [Output Files](#output-files)
- [Sample Result](#sample-result)
- [Troubleshooting](#troubleshooting)
- [Notes on Path Configuration](#notes-on-path-configuration)

---

## Overview

This repository implements an end-to-end workflow to **reproduce the cold lahar (lahar dingin) flood event on the southern flank of Mount Merapi, originating from the "Source Selatan" region**, using the **D-Claw** solver (a two-phase debris-flow model built on top of GeoClaw/AMRClaw within Clawpack).

> **Reproducibility goal:** the entire processing workflow below (Steps 1–5) is built specifically to reproduce this southern-flank lahar event as defined by the included shapefiles — `Source Selatan_1.shp` (the actual lahar source polygon used to initialize the flow) and `Source Selatan.shp` (a matching overlay used for visualization). Running the steps in order, with these two shapefiles and `DEMPredict.tt3` as provided, should regenerate the same southern-source lahar-dingin simulation, runout extent, and hazard map documented in this README. If you intend to model a **different** source region or event, you must supply your own source shapefile(s) and update the corresponding paths/constants noted in each step.

The pipeline:

1. Reads a DEM (`DEMPredict.tt3`) and a source-region shapefile to build the initial basal topography, initial mass fraction, and initial flow depth grids.
2. Runs the D-Claw AMR (adaptive mesh refinement) simulation over the domain.
3. Renders pcolor maps, transects, scatter plots, a timing report, and an MP4 animation of the flow.
4. Extracts the maximum runout distance/extent and writes a summary CSV plus shapefiles of the inundated area.
5. Overlays the simulated lahar-prone extent on a base map with rivers, roads, and settlements to produce a flood-potential (hazard) map.

## Repository Structure

```
D-Claw_Lahar_Merapi_Volcano/
├── setinput.py                       # Build basal_topo / mass_frac / surface_topo (.tt3) + verification plots
├── setrun.py                         # D-Claw run-time parameters (domain, AMR, timestepping, material properties)
├── setplot.py                # visclaw plotting: pcolor maps, transects, scatter plots, MP4 export
├── Extract Runout dan Shp file.py    # Post-run runout/sensitivity extraction -> CSV + shapefiles
├── Peta potensi banjir lahar.py      # Composite lahar flood-potential hazard map (rivers, roads, settlements)
├── DEM Analisis.py                   # Standalone DEM erosion/deposit analysis & colormap utilities
├── DEMPredict.tt3                    # Digital Elevation Model (ESRI ASCII / tt3 format) of the study area
├── Source Selatan/                   # Overlay shapefile (visualization only, not the simulation source)
│   └── Source Selatan.{shp,shx,dbf,prj,cpg,sbn,sbx}
└── Source Selatan_1/                 # Lahar SOURCE shapefile — defines the initial flow region
    └── Source Selatan_1.{shp,shx,dbf,prj,cpg,sbn,sbx}
```

> **Note:** the shapefiles you will need to prepare yourself for your own scenario (rivers, roads, settlements, alternate source polygons) are referenced by filename in `Peta potensi banjir lahar.py` (`sungai.geojson`, `jalan.geojson`, `pemukiman_merapi_scale25ribu.gpkg`, `Laharselatan.shp`, `Laharbarat-baratdaya.shp`) but are **not included** in this archive — see [Input Data](#input-data).

## Modeling Approach

The simulation uses **D-Claw**, a two-phase (solid grain + pore fluid) depth-averaged debris/lahar flow model, on an AMR grid over the DEM domain. Key modeling choices in this pipeline:

- **Domain**: derived directly from the DEM header (`DEMPredict.tt3`), so the D-Claw grid is always pixel-aligned with the DEM.
- **Source region**: a bounding box built from `Source Selatan_1.shp` (EPSG:32749 / UTM Zone 49S), buffered by 50 m, forced to the finest AMR level so the source area is always fully resolved.
- **Initial conditions**: initial lahar depth `h0`, mass fraction `m0`, and surface elevation (`eta`) are painted only inside the source polygon and written out as GeoClaw `qinit` topotype-3 files.
- **Three AMR levels** with refinement ratios `[2, 2, 2]` in space and `[2, 4, 4]` in time — finer levels take more substeps because the lahar front moves fast.
- **Conservative timestepping** (small initial `dt`, capped `dt_max`, low CFL targets) tuned for numerical stability of a fast-moving granular-fluid front.
- **Material parameters** (grain density, fluid density, critical/reference solid fraction, permeability, friction angle, etc.) are set via `rundata.dclaw_data` in `setrun.py`.
- **Gauges** are placed at fixed offsets downstream of the source to record time series of depth/velocity.
- **fgmax** grids record the maximum flow depth/speed reached at every point; **fgout** grids export full time-series snapshots of the flow for animation.

## Requirements

### Core simulation
- **Python** 3.9–3.11 (matches typical Clawpack/D-Claw supported versions)
- **Clawpack** with the **D-Claw** extension (`clawpack.amrclaw`, `clawpack.geoclaw`, `clawpack.clawutil`, `clawpack.visclaw`) — built from source, including the Fortran D-Claw solver
- A **Fortran compiler** (`gfortran`) to build the D-Claw executable
- `make`

### Python packages
| Package | Used for |
|---|---|
| `numpy` | Array/grid handling |
| `scipy` | `RegularGridInterpolator`, `KDTree` (runout extraction) |
| `shapely` | Polygon geometry (source region, overlays) |
| `matplotlib` | All plotting (pcolor maps, hillshade, scatter, colorbars) |
| `geopandas` | Reading/writing shapefiles/GeoPackages for the hazard map |
| `pandas` | Tabular handling in the hazard-map script |
| `pyproj` | Coordinate transforms (`Transformer`) |
| `ffmpeg` (system binary, not pip) | MP4 animation export from `setplot.py` |

## Installation

### 1. System dependencies

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y gfortran make git ffmpeg python3-dev python3-pip
```

### 2. Clawpack / D-Claw

D-Claw is distributed as an extension on top of Clawpack, developed and maintained by Dave George (USGS) at [geoflows/dclaw](https://github.com/geoflows/dclaw) (mirrored from [code.usgs.gov/claw/dclaw](https://code.usgs.gov/claw/dclaw)). Follow the steps below, which mirror the [official D-Claw installation instructions](https://claw.code-pages.usgs.gov/dclaw/src/installation.html).

> D-Claw is supported on **Linux/Unix (including macOS)** only. Windows users should use WSL.

**a. Clone the repositories**

```bash
# Clone Clawpack (provides amrclaw, geoclaw, clawutil, visclaw, pyclaw, ...)
git clone https://github.com/clawpack/clawpack.git
cd clawpack
git submodule init
git submodule update
source pull_all.sh

# Clone D-Claw into the clawpack tree
git clone https://code.usgs.gov/claw/dclaw.git
# (or the GitHub mirror: git clone https://github.com/geoflows/dclaw.git)
```

This places the D-Claw source in a `dclaw/` subfolder of the Clawpack directory, alongside the other Clawpack submodules (GeoClaw, AMRClaw, etc.).

**b. Set the required environment variables**

```bash
export CLAW=/path/to/top/level/clawpack/directory
export PYTHONPATH=/path/to/top/level/clawpack/directory
```

> Both variables **must** be set for D-Claw (and this pipeline) to run — `setrun.py` explicitly checks for `CLAW` at import time and raises an exception if it isn't set. Add both `export` lines to your shell profile (`~/.bashrc` or similar) so they persist across sessions.

**c. Compile D-Claw**

Navigate to one of the D-Claw example directories (`$CLAW/dclaw/examples/<example-name>`) and build a fresh executable:

```bash
cd $CLAW/dclaw/examples/<example-name>
make new
```

Useful Makefile targets (`make help` lists them all):

```text
make .objs        compile object files
make .exe         create the executable
make .data        create .data files from setrun.py
make .output      run the simulation
make output       run the simulation, no dependency checking
make .plots       produce plots (via setplot.py)
make plots        produce plots, no dependency checking
make new          remove all objects, then make .exe
make clean        clean up compilation/html files
make clobber      also clean up output and plot files
```

D-Claw also adds two extra targets on top of the standard Clawpack ones: `make input` (runs `setinput.py`, if present, to preprocess initial conditions) and `make postprocess` (runs `setpostprocess.py`, if present).

### 3. Python environment

D-Claw's own python wrapper only strictly requires `numpy` and `matplotlib`. This pipeline's pre/post-processing scripts add a few more packages on top:

```bash
python3 -m venv .venv
source .venv/bin/activate

pip install numpy matplotlib scipy shapely geopandas pandas pyproj
```

### 4. Place this pipeline in a D-Claw example directory

D-Claw's `make` workflow expects `setrun.py`/`setinput.py`/`setplot.py` (and the executable `xdclaw`) to live together inside an example folder under `$CLAW/dclaw/examples/`. Copy (or symlink) this repository's contents into a new example folder there, e.g.:

```bash
mkdir -p $CLAW/dclaw/examples/merapi_lahar_selatan
cp -r /path/to/D-Claw_Lahar_Merapi_Volcano/* $CLAW/dclaw/examples/merapi_lahar_selatan/
cd $CLAW/dclaw/examples/merapi_lahar_selatan
make new
```

## Input Data

Place the following files at the paths expected by each script (all scripts in this repo resolve most paths **relative to the script's own directory**, except where noted):

> **Key reproducibility inputs:** `Source Selatan_1.shp` and `Source Selatan.shp` are the two files that make this pipeline reproduce *this specific* southern-flank (Source Selatan) lahar-dingin event rather than a generic one. Both are already included in this repository — do not replace them unless you deliberately want to model a different source area.

| File | Required by | Notes |
|---|---|---|
| `DEMPredict.tt3` | `setinput.py`, `setplot.py`, `Peta potensi banjir lahar.py` | ESRI ASCII grid (6-line header), included in this repo |
| `Source Selatan_1/Source Selatan_1.shp` (+ sidecar files) | `setinput.py` | **Lahar source polygon** — defines where flow is initialized; the shapefile this whole pipeline reproduces from |
| `Source Selatan/Source Selatan.shp` (+ sidecar files) | `setinput.py` | Visualization overlay only, does **not** affect the simulation |
| `Laharselatan.shp`, `Laharbarat-baratdaya.shp` | `Peta potensi banjir lahar.py` | Lahar-prone zone polygons for the hazard map — **not included**, supply your own |
| `sungai.geojson` (rivers), `jalan.geojson` (roads) | `Peta potensi banjir lahar.py` | Base-map context layers — **not included** |
| `pemukiman_merapi_scale25ribu.gpkg` (layer `permukiman_area`) | `Peta potensi banjir lahar.py` | Settlement polygons for exposure visualization — **not included** |

All shapefiles/coordinates in this pipeline are in **UTM Zone 49S (EPSG:32749)**.

## Usage

Run the pipeline in the order below from inside the D-Claw example folder (see [installation step 4](#3-python-environment)), with your Clawpack environment activated (`CLAW` and `PYTHONPATH` set, virtualenv active). Steps 1–3 correspond directly to the standard D-Claw workflow (`make input` → `make .data`/`make .output` → `make .plots`); Steps 4–5 are extra post-processing scripts specific to this pipeline.

> **This sequence reproduces the southern-source (Source Selatan) cold lahar event.** Step 1 is where reproducibility is anchored: it consumes `Source Selatan_1.shp` directly to define *where* the lahar initiates and *what shape* the source region has. As long as `Source Selatan_1.shp` (and its sidecar files) and `DEMPredict.tt3` remain unchanged, re-running Steps 1–5 in order will regenerate the same Source Selatan lahar-dingin scenario end-to-end — from initial conditions, through the AMR simulation, to the final runout and flood-potential maps.

### Step 1 — Prepare topography and source inputs (`setinput.py`)

Reads `DEMPredict.tt3` and the **Source Selatan shapefiles** (`Source Selatan_1.shp` as the actual source polygon, `Source Selatan.shp` as the overlay), then writes the three GeoClaw topotype-3 input files D-Claw needs, plus verification plots. **This is the step that ties the whole simulation to the southern-source event** — everything downstream (Steps 2–5) inherits its initial conditions from here.

```bash
python "setinput.py"
# equivalent to D-Claw's standard target:
# make input
```

**What it does:**
- Loads the DEM, flips it so row 0 = south, converts pixel-edge coordinates to pixel-center coordinates.
- Cross-checks the DEM's derived domain bounds against known reference values and prints `OK`/`WARNING` diagnostics.
- Reads **`Source Selatan_1.shp`** as the lahar source polygon — the exact footprint of the southern-flank source area to reproduce — and rasterizes a mass-fraction (`m0`) and initial-depth (`h0`) field inside it.
- Writes:
  - `basal_topo.tt3` — the DEM converted directly to a D-Claw basal topography file
  - `mass_frac.tt3` — initial solid mass fraction grid, non-zero only inside the Source Selatan polygon
  - `surface_topo.tt3` — initial surface elevation (eta) grid
- Generates 4 PNG verification plots: `source_area.png`, `basal_topo_hillshade.png`, `surface_topo.png`, `landslide_depth.png`, and prints summary statistics (elevation range, mean slope, source polygon area, initial depth stats) — use `source_area.png` to visually confirm the source polygon matches the intended Source Selatan extent before proceeding.

To faithfully reproduce the documented Source Selatan event, keep these two physical parameters at their given values (they were calibrated for this event); only change them if you are deliberately testing a different scenario:

```python
m0 = 0.59   # initial solid mass fraction inside the Source Selatan area
h0 = 4.11   # initial lahar depth (meters) at Source Selatan
```

### Step 2 — Configure and run the simulation (`setrun.py`)

`setrun.py` defines all D-Claw run-time parameters and, when executed directly, writes them to the `.data` files the Fortran solver reads.

```bash
# Generate the .data files from setrun.py
make .data

# Build the D-Claw executable (first time only, or after code changes)
make .exe

# Run the simulation (writes output to _output/), regenerating .data if needed
make .output

# Or, to just re-run without dependency checking:
make output
```

**What it configures:**
- Domain bounds and cell counts (derived automatically from the DEM header — do not need manual editing unless you change the DEM).
- Simulation length: `tfinal = 1620.0` s (27 minutes) with output every 10 s (162 frames).
- Stability-tuned timestepping: small initial `dt`, `dt_max = 30 s`, `cfl_desired = 0.25`, `cfl_max = 0.40`.
- 3-level AMR with a forced level-3 flag region over the buffered source polygon, and `regrid_interval = 1` (regrid every step, since the lahar front moves fast).
- Boundary conditions: `extrap` (open/outflow) on all four domain edges.
- GeoClaw physics: gravity, Manning friction (`n = 0.025`), dry tolerance `1e-3`.
- D-Claw material parameters (`rho_s`, `rho_f`, `m0`, `m_crit`, permeability, friction angle `phi = 35°`, viscosity, etc.) via `rundata.dclaw_data`.
- 4 gauges placed downstream of the source at 100/1000/3000/5000 m intervals.
- `fgmax` grid (25 m resolution) for maximum flow depth/speed.
- `fgout` grid at full DEM resolution for animation export (`q_out_vars = [depth, momentum, eta]`).

### Step 3 — Visualize results (`setplot.py`)

Generates the standard Clawpack/visclaw HTML+PNG plot gallery, plus a custom MP4 animation.

```bash
make .plots
# or, without dependency checking:
make plots
```

This calls `setplot.py`'s `setplot()` function internally. To (re)generate just the MP4 from already-produced PNG frames without a full re-plot:

```bash
python setplot.py _plots
```

**Figures produced:**
| Figure # | Name | Content |
|---|---|---|
| 1 | Banjir Lahar Dingin Merapi | Main pcolor map (plan view) + downstream transect panel |
| 2 | Mass Fraction | Solid mass-fraction field over time |
| 6 | Depth | Flow depth field over time |
| 3 | scatter_depth | Depth vs. downstream distance scatter |
| 9 | scatter_speed | Flow speed vs. downstream distance scatter |
| — | timing | Solver performance/timing statistics (`plot_timing_stats`) |
| — | mp4_export | `lahar_utama.mp4` animation assembled from the pcolor + mass-fraction + depth frames |

Settings of note in the script:
- `TFINAL` is automatically imported from `setrun.py` (falls back to `1620.0` s if not found — keep the two scripts in sync).
- `SHOW_KALI_ADEM = True` toggles display of the Kali Adem reference point/marker on the map and transect.
- `MP4_FIGNOS = [1, 2, 6]` controls which figures are combined into the animation; `MOVIE_FPS`, `MOVIE_SLOWDOWN`, and `MP4_CRF` control playback speed and encoding quality.
- HTML output links back to `../README.html` (`plotdata.html_homelink`) — if you rename or move the output folder, update this path.

### Step 4 — Extract runout and sensitivity data (`Extract Runout dan Shp file.py`)

Post-processes the `_output` directory to compute the runout distance from the source and export summary tables/shapefiles.

```bash
python "Extract Runout dan Shp file.py"
```

**What it does:**
- Reads `_TFINAL` from `setrun.py` to synchronize the end-time (minutes) used in the sensitivity summary.
- Uses a KD-tree search over the fgout/fgmax grid to find cells above `H_THRESHOLD = 0.01` m depth.
- Computes runout distance/extent relative to `SOURCE_POINT` (kept in sync with `setrun.py`'s source centroid) and the Kali Adem reference coordinate.
- Writes `hasil_sensitivitas.csv` (sensitivity/runout summary) into the script directory, plus shapefiles of the inundated extent in `EPSG:32749`.

Before running for a new scenario, verify these constants at the top of the script match your current `setrun.py`:

```python
SOURCE_POINT     = (439164.988, 9165667.463)   # must equal (_src_cx, _src_cy) in setrun.py
KALI_ADEM_COORD  = (439540.56, 9161173.13)
KALI_ADEM_RADIUS = 100.0
H_THRESHOLD      = 0.01
```

### Step 5 — Generate the flood-potential map (`flood_potential_map.py`)

Builds a composite hazard/exposure map: DEM hillshade base, lahar-prone zone shapefiles, rivers, roads, and settlement polygons.

```bash
python flood_potential_map.py
```

**Before running**, update `WORK_DIR` and the filenames at the top of the script to point at your local data, or place all inputs alongside the script (the script searches `_SEARCH_DIRS` for each file):

```python
FILE_DEM        = "DEMPredict.tt3"
FILE_SHP_SEL    = "Laharselatan.shp"
FILE_SHP_BAR    = "Laharbarat-baratdaya.shp"
FILE_SUNGAI     = "sungai.geojson"
FILE_JALAN      = "jalan.geojson"
FILE_PERMUKIMAN = "pemukiman_merapi_scale25ribu.gpkg"
LAYER_PERMUKIMAN = "permukiman_area"
WORK_DIR        = r"D:\Tugas Akhir\Peta Potensi Banjir Lahar Dingin"   # <- change this
```

**What it does:**
- Builds a DEM hillshade basemap using `matplotlib.colors.LightSource`.
- Overlays lahar-prone zone polygons (dissolved with `shapely.ops.unary_union` where needed), rivers, and roads.
- Colors settlement polygons (`WARNA_PERMUKIMAN`) to highlight exposure near the hazard zones.
- Reprojects/transforms coordinates as needed via `pyproj.Transformer`.
- Exports the final composite hazard map as a PNG.

### Optional — DEM erosion/analysis utility (`DEM Analisis.py`)

A standalone script (from a separate `radial_slide` example scenario) for analyzing a DEM before/after a simulation — computing an erosion/deposition difference DEM and rendering it with a custom truncated colormap.

**Purpose:** this script cuts (excavates) the source area where the lahar material collapsed/detached, and fills the deposition area produced by the cold lahar flow, directly on the DEM. The goal is to produce a **new, updated DEM** that reflects the topographic changes caused by the lahar event, so the DEM becomes more representative of actual post-event conditions.

This is especially useful for simulating **a sequence of multiple cold lahar events**: instead of reusing the original pre-event DEM for every simulation, each subsequent simulation can use the DEM already updated with the topographic changes (erosion at the source, deposition downstream) left behind by the previous event — making multi-event simulations more physically realistic.

```bash
python "DEM Analisis.py"
```

> This script currently hardcodes absolute paths from the original development machine (see [Notes on Path Configuration](#notes-on-path-configuration)) — edit the `CONFIGURATION` block at the top before running:

```python
OUTPUT_DIR  = "/path/to/your/_output"
DEM_IN      = "/path/to/your/DEMSkenario.tt3"
SHP_IN      = "/path/to/your/source_polygon.shp"
DEM_OUT     = "/path/to/your/DemFinal.tt3"
PLOT_DIR    = "/path/to/your/plot_output"

EROSI_DEPTH = 8.44     # erosion depth assumption (m)
FGOUT_FGNO  = 1
T_TARGET    = 1620.0
H_THRESHOLD = 0.01
```

## Key Parameters

Quick reference for the parameters you're most likely to need to change between scenarios:

| Parameter | Location | Default | Meaning |
|---|---|---|---|
| `m0` | `setinput.py` | `0.59` | Initial solid mass fraction in the source region |
| `h0` | `setinput.py` | `4.11` m | Initial lahar depth in the source region |
| `_TFINAL` | `setrun.py` | `1620.0` s | Total simulation time (also drives frame count) |
| `_FLAG_BUF` | `setrun.py` | `50.0` m | Buffer added around the source polygon for the forced-refinement region |
| `amrdata.amr_levels_max` | `setrun.py` | `3` | Max AMR levels |
| `clawdata.cfl_desired` / `cfl_max` | `setrun.py` | `0.25` / `0.40` | Stability targets for adaptive timestepping |
| `dc.phi` | `setrun.py` | `35` (deg) | Internal friction angle of the granular material |
| `dc.rho_s` / `dc.rho_f` | `setrun.py` | `2500.0` / `1150.0` kg/m³ | Solid grain / pore fluid density |
| `H_THRESHOLD` | `Extract Runout dan Shp file.py` | `0.01` m | Minimum depth counted as "inundated" |
| `MP4_FIGNOS` | `setplot.py` | `[1, 2, 6]` | Figures combined into the MP4 animation |

## Output Files

After a full run (`setinput.py` → `make .output` → `make plots` → post-processing scripts), you should have:

```
basal_topo.tt3, mass_frac.tt3, surface_topo.tt3   # inputs written by setinput.py
source_area.png, basal_topo_hillshade.png,
surface_topo.png, landslide_depth.png             # verification plots
_output/                                          # raw D-Claw solution + fgmax/fgout data
_plots/                                           # visclaw HTML + PNG gallery
_plots/lahar_utama.mp4                            # flow animation
hasil_sensitivitas.csv                            # runout/sensitivity summary
<runout shapefiles>                               # inundated-extent polygons (EPSG:32749)
<flood potential map>.png                          # composite hazard map
```

## Sample Result

The figure below (`frame0162`, the final frame — Figure 1 from `setplot.py`) shows what a successful run of the Source Selatan reproduction pipeline produces at the end of the simulation (`t = 0:27:00`, i.e. `t = 1620 s = tfinal`):

![Sample result: Source Selatan cold lahar simulation, final frame](sample_result_source_selatan.png)

**Top panel — plan-view flow map:**
- Hillshaded DEM of the Merapi summit and southern flank, in UTM Zone 49S (Easting/Northing, meters).
- The **dashed blue outline ("Source Selatan")** marks the lahar source polygon from `Source Selatan_1.shp` — this is where the flow is initialized in Step 1.
- The colored trace running downstream from the source, following the **Kali Gendol** river channel, shows simulated lahar depth (colorbar, "Kedalaman lahar", 0–~15 m) at the final output time.
- The **blue circle "Kali Adem"** marks the fixed downstream reference point (also used by `Extract Runout dan Shp file.py` to report runout relative to this location).
- The **gray transect line "Transek B–B' (Kali Gendol)"** marks the downstream profile line B→B' plotted in the bottom panel.

**Bottom panel — downstream transect (B–B', Sungai Gendol):**
- Elevation profile along the flow path from the source (**B**) to the transect endpoint (**B'**), ~9+ km downstream.
- Colors along the profile classify the flow regime by solid mass fraction `m`: **mud flow** (`m < 0.3`, blue), **hyperconcentrated flow** (`0.3 ≤ m < 0.6`, orange), and **lahar** (`m ≥ 0.6`, dark red).
- The title reports the transect summary statistics at the final time: **`hmax = 6.65 m`** (maximum flow depth along the transect) and **`mmean = 0.62`** (mean solid mass fraction — i.e., the flow is, on average, in the lahar/debris-flow regime by the time it reaches the transect).
- The dashed vertical line marks the **Kali Adem** reference point (~5523 m from the source along the transect), matching the marker in the top panel.

This is the kind of result — a plan-view depth map plus a calibrated downstream transect with flow-regime classification and runout statistics — that Steps 1–3 (`setinput.py` → `setrun.py` → `setplot.py`) reproduce for the Source Selatan event, one PNG frame per output time step (162 frames total for `tfinal = 1620 s`), which are then assembled into `lahar_utama.mp4` (see [Step 3](#step-3--visualize-results-setplotpy)).
