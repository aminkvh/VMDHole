# Troubleshooting

Read the VMD console first. VMDHole reports the selected engine, frame, failed
stage, and output path there.

## Plugin does not appear

| Symptom | Check |
|---|---|
| `package require vmdhole` fails | `auto_path` must contain the parent of `vmdhole`; keep `pkgIndex.tcl` beside `vmdhole.tcl` |
| Extension entry is absent | Restart VMD after changing `.vmdrc`, then evaluate `package require vmdhole 1.0` in the Tk Console |

## Pore run fails or finds the wrong path

| Symptom | Likely cause and action |
|---|---|
| Missing atomic radius | Select a radius file covering every selected atom, or revise the selection |
| No readable profile | Put `CPOINT` inside the pore; correct `CVECT`; inspect HOLE output in the frame directory |
| Centreline enters an external groove | Use an explicit two-point vector and a more specific atom selection |
| Profile shifts across frames | Make the molecule whole, align the trajectory, and consider endpoint stabilization |
| Unexpected differences between reruns | Confirm the same seed, executable, atom order, radius file, method, and parameters; a blank seed resolves to `1` |
| Old result would be replaced | Use a new output directory or cancel the overwrite confirmation |

## Surface problems

| Symptom | Action |
|---|---|
| Profile exists but no mesh | Check `sph_process` and `sos_triangle` paths and their console output |
| Connolly produces zero triangles at high density | Use the accelerated surface converter or lower dot density |
| Surface contains holes or spikes | Inspect Dots/Wireframe, adjust dot density, and verify that the centreline and terminal trimming are valid |
| Playback stalls | Prebuild surfaces, increase cache within available memory, or reduce playback detail |

## Tunnel problems

| Symptom | Action |
|---|---|
| No routes | Place the origin inside a buried void; reduce probe/bottleneck thresholds; review Strict Interior and exit constraints |
| Too many similar routes | Enable within-frame clustering or adjust its cutoff |
| Route identity jumps between frames | Make the molecule whole and confirm that Tunnel alignment remains enabled before cross-frame clustering |
| Expected route is absent | Increase ranks per frame, lower the Seen floor, or relax cross-frame maximum deviation |
| Lining selects repeated residues | Check chain and segment identifiers; inspect the generated VMD selection |

## Hydration and ion analysis

| Symptom | Action |
|---|---|
| Hydration Compute finds no water | Use a selection containing one oxygen per explicit water and verify frame/box preparation |
| Water-derived properties are absent | Complete Hydration Compute first; **Per-frame ρ** requires per-frame hydration data |
| Ion Flow shows no species | Verify VMD residue/element naming and that ions enter the analysed radial region |
| Passage and permeation counts differ | Expected: passage includes near-pore entry; permeation requires a complete bulk-to-bulk crossing |
| Permeation warning about sampling or PBC | Re-image the trajectory, use sufficiently frequent saved frames, and treat unsupported triclinic handling as unvalidated |

## Performance and storage

| Symptom | Action |
|---|---|
| A run is unexpectedly slow | Confirm native executable paths; the console identifies Tcl fallbacks |
| Host is overloaded | Reduce Parallel jobs to the allocated CPU count |
| Temporary storage fills | Reduce parallel jobs, choose disk-backed output, and remove only abandoned VMDHole scratch directories |
| Memory grows during playback | Reduce Surface cache and mean frame cap |

When reporting a problem, include the VMDHole and VMD versions, operating
system, complete console message, settings, and a minimal input.
