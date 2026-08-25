# Files, settings, and execution

## Result directories

With **Save results** enabled, VMDHole creates a saved run with per-frame data
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

**File → Settings** configures `hole`, `sph_process`, `sos_triangle`, the tunnel
engine, and the radius file. Use the detected acceleration status to confirm the
selected binaries. See [Installation](installation.md) for recommended builds.


## Performance controls

| Control | Trade-off |
|---|---|
| Parallel jobs | More simultaneous frames; higher CPU and temporary-storage use |
| Prebuild surfaces | Longer initial run; smoother later playback |
| Surface cache | More memory; fewer mesh rebuilds |
| Playback detail | Coarser interactive display; final selected frame remains full detail |
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

VMDHole uses the coordinates loaded in VMD. For periodic trajectories, make
molecules whole before analysis. Align frames when comparing pore positions or
matching tunnel routes; Tunnel mode enables alignment by default for multi-frame
searches and applies it to the loaded VMD frames. Hydration and ion analyses
also require correct periodic imaging and saved-frame spacing.
