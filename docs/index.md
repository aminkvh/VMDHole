# VMDHole documentation

VMDHole turns a structure or trajectory loaded in VMD into pore and tunnel
profiles, linked 3D geometry, property annotations, trajectory-wide views, and
exportable figures and tables. Independent frames run in parallel, and native
engines accelerate production calculations.

## Choose the analysis

| If you need to… | Use | Define first | Main results |
|---|---|---|---|
| Characterize one known channel through a structure or trajectory | **Pore (HOLE)** | A point inside the pore (`CPOINT`) and its direction (`CVECT`) | Radius and cross-section profiles, Connolly lateral openings, Capsule geometry, conductance estimates, hydration, ion movement, and permeation |
| Find routes from a buried site to the molecular surface | **Tunnel (MOLE 2)** | A point inside the buried site, or automatic origins | Ranked and clustered routes, bottlenecks, lining residues, property annotations, and route occurrence across frames |

Both workflows keep plots, CSV data, and 3D representations connected to the
current VMD frame. For method derivations and scientific validation, read the
paper and cited methods; these pages focus on operating the plugin and
reporting an analysis.

## Start here

1. [Install VMDHole and native engines](installation.md).
2. Complete the [first pore and tunnel analyses](quickstart.md) (**highly recommended**).
3. Choose a [tutorial](tutorials.md) for the analysis you need.
4. Continue with the matching [workflow](workflows.md).
5. Consult the [parameter reference](parameters.md) for every control and
   default.

## Analyse and interpret

| Page | Scope |
|---|---|
| [Pore workflow](pore-mode.md) | Define the channel, run HOLE, and use the pore analysis views |
| [Tunnel workflow](tunnel-mode.md) | Find, rank, cluster, track, and display routes |
| [Parameter reference](parameters.md) | GUI controls, defaults, radius files, and reporting values |
| [Visualization](visualization.md) | 3D representations, color, material, lining, and playback |
| [Properties](properties.md) | Property definitions, averaging, and picker relationships |
| [Exports](exports.md) | Figures, CSV tables, saved runs, and tunnel lining data |
| [Files and settings](files-and-settings.md) | Output layout, persistence, caching, and performance controls |
| [Scripting](scripting.md) | Tcl and headless operation |

## Verify and troubleshoot

| Page | Scope |
|---|---|
| [Troubleshooting](troubleshooting.md) | Symptom-based diagnosis |
| [Testing](testing.md) | The test tiers, every group and unit test, and the conventions they enforce |
| [References](references.md) | Methods, property scales, and required citations |

## Before interpreting a result

Confirm all of the following:

- The atom selection contains the intended protein, cofactors, and blockers,
  but excludes irrelevant solvent or membrane atoms.
- The radius file covers every selected atom.
- The centreline or tunnel lies in the expected cavity in the 3D view.
- Periodic trajectories are made whole and consistently imaged.
- Frames are aligned when a fixed pore axis or cross-frame tunnel identity is
  used.
- The same method, radii, sampling, seed policy, and trajectory preparation are
  reported for comparisons.
