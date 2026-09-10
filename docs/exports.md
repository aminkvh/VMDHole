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
| Ion & Water | plotted occupancy or passage data and species metadata |

Additional exports include summary metrics, bottleneck residues, unrolled
pore-wall layers, and tunnel lining data. Inspect the CSV header: it is the
authoritative statement of columns and units for that export.

## Save Package

**File > Save Package** first asks which of the seven plot tabs (Pore Profile,
Over Time, Trends, Mean Profile, Histogram, Hydration, Ion & Water) to include -
a tab without data cannot be ticked - then where to write the folder, and
writes the CSV and figure of each ticked tab, so a whole analysis can be
handed on in one piece. It runs the same exporters the per-tab Export menus
run, so the files are identical to the ones those buttons write. It does not
reach the plugin's other exports - bottleneck residues, unrolled pore-wall
layers, tunnel lining, per-opening tables, cavity CSVs - which stay on their
own dialogs.

The folder also contains:

| File | Contents |
|---|---|
| `run_<id>.txt` | selection, frame summary, pore axis, and pore calculation settings; not a complete record of tunnel or downstream-analysis settings |
| `README.txt` | which tabs were included, which were skipped, and the file list |

Check `README.txt` and the files before sharing the package, especially after
an interrupted export. Save Package is an export bundle, not a reloadable saved
run. Record additional settings from the [minimum reporting set](parameters.md#minimum-reporting-set).

## Run identity

Pore runs receive a timestamp and a short settings hash, recorded in
`run_<id>.txt` and the Log. This label is not a checksum of the input coordinates.
Export filenames use separate suffixes (`_run2`, `_run3`, ...) when runs share
an output folder. Use a separate folder for each analysis to keep figures,
data, and settings together.

## Connolly openings

With **Color** set to `pore_lobes`, each region gear exports the pore or one
lateral opening. The header gear exports all regions. The all-frame table
includes occurrence, dot count, neck, extension, axial position, and azimuth.
An opening absent from a frame has `present=0` and blank measurement cells.
A neck equal to the run's end radius means the route into the opening is wider than the search looks (`open`).

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
root contains provenance/manifest data. Use **File → Load Saved Analysis…** to restore a saved
HOLE or tunnel calculation without executing the engine again. Imported data
can be plotted and exported, but analyses that need the original trajectory
coordinates, such as ion tracking, also require the matching molecule and frames
to be loaded in VMD.
