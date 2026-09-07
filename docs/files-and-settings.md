# Files, settings, and execution

## Result directories

With **Save results** enabled, VMDPathFinder creates a saved run with per-frame data
and a manifest. Use **File → Import** to restore it; do not rearrange files
inside the saved run.

| Setting | Effect |
|---|---|
| Save results | Persist outputs; otherwise use temporary storage |
| Overwrite | Recalculate requested frames after confirmation |
| Keep input PDB | Preserve coordinates used for the calculation |
| Keep visualization | Retain generated VMD objects when results are reset or replaced |

Use a new output directory for a scientifically different parameter set. An
overwrite confirmation protects files; it does not determine whether two runs
are comparable.

## Persistent configuration

**Set default** saves settings for later sessions. Saved runs and exported CSV
files record individual analyses.

## Executable settings

**File → Settings** configures the surface mesher and its grid, `hole`,
`sph_process`, `sos_triangle`, the tunnel engine, and the radius file. The
mesher, the Nelder-Mead search and the Connolly classifier are part of
`sos_triangle` and need no paths of their own. Options that only apply to one
choice appear after it: the grid and neck entries follow the marching-cubes
mesher, the dot density entry sits beside `sos_triangle` and its playback
rows follow it. In HOLE Parameters the Monte Carlo rows, the random seed, SHORTO
and the extra cards follow the Search picker and show for Monte Carlo only. The
Connolly surface engine lives with the HOLE parameters and appears only for
Connolly runs. Use the detected acceleration status to confirm the
selected binaries. See [Installation](installation.md) for recommended builds.


## Performance controls

| Control | Trade-off |
|---|---|
| Parallel jobs | More simultaneous frames; higher CPU and temporary-storage use |
| Prebuild surfaces | Longer initial run; smoother later playback |
| Surface cache | More memory; fewer mesh rebuilds |
| Smoothing | Average the surface over neighbouring analysed frames: *Follow VMD* takes the trajectory-smoothing window of the shown representations, *Off*, or a fixed half-width. The marching mesher averages the frames' distance fields, sos_triangle averages the dot clouds dot by dot; numbers stay per frame |
| Playback stride (`sos_triangle` mesher only) | Draw every Nth triangle while the trajectory plays; the frame at rest is full detail. The marching-cubes mesher draws one full-detail mesh throughout |
| Mean frame cap | Bounds expensive mean/property work on long trajectories |
| Accurate 3D | Better property projection; higher calculation cost |

The automatic job count normally leaves one processor available. On shared
systems, set an explicit value consistent with the scheduler allocation.

## Temporary storage

Frame jobs may use system temporary storage or `/dev/shm` when available.
Large systems and many parallel jobs can require substantial space. If a run is
interrupted abnormally, check the VMD console for the scratch location and
remove only directories known to belong to that stopped run.

## Trajectory preparation

VMDPathFinder uses the coordinates loaded in VMD. For periodic trajectories, make
molecules whole before analysis. Align frames when comparing pore positions or
matching tunnel routes; Tunnel mode enables alignment by default for multi-frame
searches and applies it to the loaded VMD frames. Hydration and ion analyses
also require correct periodic imaging and saved-frame spacing.
