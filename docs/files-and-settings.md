# Files, settings, and execution

## Result directories

Every pore run gets its own folder, `<structure>_<date>-<time>-<hash>`, under
`hole_output_<structure>` beside the structure file (or under the chosen work
directory; under the temp directory when **Save results** is off). It holds the
per-frame data, `run_<id>.txt` with the settings, and a manifest that records
the structure's atom count and radius of gyration. Use **File → Load Saved
Analysis…** on a run folder to restore that run, or on `hole_output_<structure>`
to restore every run in it as a memory each. A run loaded onto a different
structure is flagged before anything is drawn. Do not rearrange files inside a
run folder.

| Setting | Effect |
|---|---|
| Save results | Persist outputs; otherwise use temporary storage |
| Overwrite | Tunnel runs: recalculate requested frames after confirmation (pore runs never overwrite - each gets a new folder) |
| Keep input PDB | Preserve coordinates used for the calculation |
| Keep visualization | Retain generated VMD objects when results are reset or replaced |

Two runs are comparable only if their settings say so: read `run_<id>.txt`.

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
