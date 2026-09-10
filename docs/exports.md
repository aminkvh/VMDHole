# Export and import

Each plot tab provides an **Export** menu for the displayed figure and the data
used to draw it. Export the CSV together with a publication figure so processing
and units remain inspectable.

Figures are written directly as EPS. JPEG export uses Ghostscript or
ImageMagick when available and otherwise falls back to EPS; check the status
message and resulting extension.

## Plot exports

| Tab | Principal CSV content |
|---|---|
| Pore Profile | channel coordinate, radius, and selected fill property where applicable |
| Over Time | frame-by-position radius or property matrix |
| Mean Profile | coordinate, mean radius, standard deviation, and contributing-frame count |
| Trends | frame and selected metric |
| Histogram (radius summary) | axial-bin coordinate and uncapped mean/minimum/maximum radius aggregate |
| Hydration | coordinate, relative density, free energy, waters per frame, and available standard deviations |
| Ion Flow | plotted occupancy or passage data and species metadata |

Additional exports include summary metrics, bottleneck residues, unrolled
pore-wall layers, and tunnel lining data. Inspect the CSV header: it is the
authoritative statement of columns and units for that export.

## Save Package

**File > Save Package** writes one folder holding the CSV and figure from each
of the seven plot tabs (Pore Profile, Over Time, Trends, Mean Profile,
Histogram, Hydration, Ion & Water) that has data, so a whole analysis can be
handed on in one piece. It runs the same exporters the per-tab Export menus
run, so the files are identical to the ones those buttons write. It does not
reach the plugin's other exports - bottleneck residues, unrolled pore-wall
layers, tunnel lining, per-opening tables, cavity CSVs - which stay on their
own dialogs.

The folder also contains:

| File | Contents |
|---|---|
| `run_<id>.txt` | every parameter the run used - selection, frames, CPOINT/CVECT, method, engines, sampling, seed, radius file |
| `README.txt` | which tabs were included, which were skipped, and the file list |

A tab is listed as skipped when it had nothing computed, when pressing Abort
stopped the package before reaching it, or when both its exporters ran but
neither actually wrote a file - so a tab in "Included" always has at least one
real file behind it.

## Run identity

Every run is stamped with an id: the date, the time, and a short hash of every
input that decides what HOLE computes. Two runs with different parameters
cannot share an id, and neither can two runs a second apart. It appears in the
`run_<id>.txt` file written beside the results and in a one-line summary in the
Log. Export filenames carry a separate, shorter tag (`_run2`, `_run3`, ...)
that only distinguishes runs sharing one output folder - the first run gets no
tag at all - so a figure and its `run_<id>.txt` are matched by which folder
they sit in, not by a shared string in their names.

## Connolly openings

With **Color** set to `pore_lobes`, each region gear exports the pore or one
lateral opening. The header gear exports all regions. The all-frame table
includes occurrence, dot count, neck, extension, axial position, and azimuth.
An opening absent from a frame has `present=0` and blank measurement cells.

## Filenames

VMDPathFinder suggests a descriptive filename for each export. Confirm the destination
before saving.

## Hydration export

Hydration CSV export always writes the mean density/free-energy profile,
independent of the currently displayed hydration view. Per-frame matrix and
hydrophobicity plot names must not be interpreted as extra columns in that
profile file.

## Tunnel lining

In tunnel mode, **Lining → Export** writes lining data for the selected route.
Standard plot exports remain CSV.

## Saved runs and import

When **Save results** is enabled, each frame has a result directory and the run
root contains provenance/manifest data. Use **File → Import** to restore a saved
HOLE or tunnel calculation without executing the engine again. Imported data
can be plotted and exported, but analyses that need the original trajectory
coordinates, such as ion tracking, also require the matching molecule and frames
to be loaded in VMD.
