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

## Connolly openings

With **Color** set to `pore_lobes`, each region gear exports the pore or one
lateral opening. The header gear exports all regions. The all-frame table
includes occurrence, dot count, neck, extension, axial position, and azimuth.
An opening absent from a frame has `present=0` and blank measurement cells.

## Filenames

VMDHole suggests a descriptive filename for each export. Confirm the destination
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
