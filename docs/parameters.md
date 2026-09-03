# Parameter reference

This page lists VMDHole controls and their starting values. Saved settings can
change a default.

## Shared input

| Parameter | Default | Definition |
|---|---:|---|
| Molecule | `top` | VMD molecule to analyse |
| Selection | `protein` | Atoms that define the pore or molecular interior |
| Frames | `now` | `now`, `all`, one frame, `start:end`, or `start:stride:end` |
| Align trajectory | off in pore mode; on for tunnel tracking | RMSD-fit frames to a reference selection before analysis |
| Alignment selection | `protein and name CA` | Atoms used for the fit |
| Reference frame | `0` | Frame used as alignment reference |

## Pore geometry

| Parameter | Default | Definition |
|---|---:|---|
| `CPOINT` | VMD centre of rotation | Initial point inside the pore, in Å |
| `CVECT` | `{0 0 1}` | Direction of the channel axis |
| Stabilize `CPOINT` | off | Carry the point with a local rigid-body fit in every frame |
| Local fit radius | 10 Å | Initial neighbourhood for `CPOINT` or endpoint stabilization |
| Auto-expand fit radius | 15 Å | Larger neighbourhood used if the initial fit has too few atoms |
| Fit RMSD warning | 3 Å | Console-warning threshold for an unstable local fit |
| Track `CPOINT` | off | Re-centre the point on selected atoms near its previous position |
| Track radius | 8 Å | Neighbourhood used by `CPOINT` tracking |
| Stabilize `CVECT` endpoints | off | Carry each two-point endpoint with an independent local fit |
| Exact `CVECT` selection | off | Re-evaluate both point selections in every frame; mutually exclusive with endpoint stabilization |
| Show cues | off | Mark `CPOINT`, `CVECT`, or the tunnel start in the VMD viewer |
| Radius file | configured HOLE file | HOLE `.rad` file containing van der Waals radii for every selected atom |
| `SAMPLE` | 0.25 Å | Distance between successive search planes |
| `ENDRAD` | 15 Å | Radius at which HOLE considers an end open to bulk |
| `SHORTO` | 1 | HOLE output verbosity; `1` suppresses detailed per-slice diagnostics while retaining the profile output |
| `IGNORE` | `HOH WAT TIP SOL` | Residue names excluded from the HOLE calculation |
| Random seed | blank (uses `1`) | Integer seed; blank uses `1` |

### Radius files

The radius file assigns a van der Waals radius to every selected atom. Choose a
HOLE-format `.rad` file in **File → Settings** and confirm that the console
reports no missing atom types. Use a custom file for cofactors, ligands,
modified residues, or nonstandard naming.

Atomic radii are separate from the bare, hydrated, and probe radii used by
**Passability**. Their sources are listed in
[References](references.md#pore-geometry-and-conductance).

### HOLE method and search

| Parameter | Default | Definition |
|---|---:|---|
| Method | Spherical | Spherical probe, Connolly accessible cross-section, or Capsule anisotropic probe |
| Dot density | 15 | Surface sampling density used by Connolly/surface stages |
| Monte Carlo steps | blank (HOLE default 1000) | Optimization steps per search plane |
| Monte Carlo step size | blank (HOLE default 0.1 Å) | Trial displacement scale |
| Monte Carlo temperature | blank (HOLE default 0.1) | Simulated-annealing acceptance parameter |
| Trim Connolly | off/default | Trim terminal Connolly geometry to the analysed pore |
| Hide sideways spill | off | Remove Connolly surface regions classified as lateral spill |
| Margin | 2 Å | Distance beyond the traced pore wall that still belongs to the central pore |
| Opening axial match tolerance | 6 Å | Maximum axial displacement used to match a lateral opening across frames |
| Opening angular match tolerance | 35° | Maximum azimuthal displacement used to match a lateral opening across frames |
| Opening Seen floor | 25% | Minimum percentage of analysed frames required to list a tracked opening |
| Extra cards | blank | Additional HOLE control cards, separated by semicolons |

## Tunnel search

| Parameter | Default | Definition |
|---|---:|---|
| Start point | blank | Buried origin coordinate; may be derived from a selection's centre of geometry (COG) or VMD's centre of rotation (COR) |
| Auto-detect origins | off | Detect candidate internal cavities rather than use only the entered point |
| Probe | 3.0 Å | Probe used in accessible-space construction |
| Interior threshold | 1.25 Å | Minimum clearance used to classify interior voids |
| Origin radius | 5.0 Å | Radius around the requested start used to accept seeds |
| Minimum length | 0 Å | Minimum accepted route length |
| Bottleneck | 1.25 Å | Minimum accepted bottleneck radius |
| Minimum depth | 8 Å | Required graph depth of an interior seed |
| Minimum depth length | 5 Å | Minimum length associated with the depth criterion |
| Surface cover radius | 10 Å | Surface neighbourhood used in terminal classification |
| Auto-origin cover radius | 10 Å | Cover radius used during automatic-origin selection |
| Maximum origins | 5 | Maximum automatically selected origins |
| Bottleneck tolerance | 0 Å | Tolerance applied to the bottleneck filter |
| Maximum similarity | 0.9 | MOLE similarity threshold for redundant routes |
| Weight function | `VoronoiScale` | Route cost: VoronoiScale, LengthAndRadius, Length, or Constant |
| FBL | off | Enable the MOLE FBL option |
| Strict Interior | off | Apply a stricter interior classification; may yield no routes for borderline inputs |

### Tunnel exits and clustering

| Parameter | Default | Definition |
|---|---:|---|
| Custom exit point | blank | User-defined surface target |
| Custom path start/end | blank | Directed path constraint |
| Use custom exits only | off | Exclude automatically found exits |
| Cluster within frame | on | Merge similar routes in one frame |
| Within-frame cutoff | 3 Å | Distance threshold for within-frame route clustering |
| Cross-frame maximum deviation | 12 Å | Largest geometric deviation accepted as one tracked route |
| Ranks per frame | 10 | Highest-ranked routes admitted to cross-frame matching; `0` means all |
| Seen floor | 40% | Minimum frame occupancy shown for a tracked route |
| Align trajectory | on | Fit frames before route matching |
| Pre-mesh budget | 400 | Limit for eagerly prepared route meshes |
| Draft detail | 1 | Tunnel surface sampling stride |
| Tunnel dot density | 15 | Surface sampling density |
| Accurate 3D coloring | off | Project a tunnel property by true 3D surface distance |

## Visualization and property controls

| Parameter | Default | Definition |
|---|---:|---|
| Pore representation | Isosurface | None, Centerline, Dots, Wireframe, or Isosurface |
| Tunnel representation | Isosurface | Isosurface, Wireframe, or Centerline; global or per route |
| Pore surface color | `hole_def` | HOLE radius banding, property, `pore_lat` pore/spill classification, `pore_lobes` individual Connolly openings, or a flat VMD color |
| Tunnel surface color | automatic rank | Route/rank color, selected property, or a flat VMD color |
| Material | Opaque | VMD material applied to the generated representation |
| Playback detail | 4 | Frame stride/detail used while the trajectory slider is moving |
| Synchronize playback | on | Update VMDHole geometry with the VMD frame |
| Pore lining threshold | 3 Å | Maximum atom-to-local-surface distance used to classify lining residues for display and residue-property averaging |
| Property smoothing | 3 Å | Axial smoothing bandwidth |
| Pore-facing only | on | Retain side chains directed toward the lumen where applicable |
| Accurate 3D property projection | on | Project property to the true surface rather than a faster approximation |
| Bottleneck shell | 3 Å | Independent surface-distance cutoff for the bottleneck-residue report |
| Scale bar | on/default | Show property legend in the VMD scene |
| Scale-bar font/corner | saved/default | Legend color and screen position |
| Metrics readout | on/default | Show pore summary and selected ion passability values |
| Metrics species | Water, K, Na, Ca | Species in the on-figure passability summary; Mg, Cl, Li, and Cs are also available |

The 3D surface, Pore Profile Fill, and Mean Profile synchronize a property where
it is available. **Over Time** has an independent property and an explicit
**Compute** gate. Kapcha–Rossky is atom-level in both modes.

## Plot controls

| Tab | Controls |
|---|---|
| Pore Profile | None, Fill, Ellipse fit, or Unrolled; property/layer; ellipse solid/point rendering; swap axes; flip direction |
| Over Time | Radius or Property; HOLE or Ellipse radius source; color scheme; independent property; Compute; flip Y |
| Mean Profile | 2D fill; property; 3D isosurface; color; material; accurate 3D; frame cap; Render smoothly (off); swap/flip |
| Trends | metric; mean overlay; conductivity preset/custom value for conductance; constriction shell |
| Histogram (radius summary) | mean (default), minimum, or maximum radius over 50 fixed spatial bins; swap/flip |

See [Properties](properties.md#unrolled-map-layers) for Unrolled layers.
**Connolly reach** is available for Connolly pores. Plot swap and flip options
change presentation only. The unrolled map requires a native HOLE binary.

## Hydration

| Parameter | Default | Definition |
|---|---:|---|
| Water selection | `water and oxygen` | One representative oxygen atom per water molecule |
| Temperature | 310 K | Temperature in `-RT ln(rho/rho_bulk)` |
| Bin width | 1 Å | Axial density-bin width |
| Density probe cap | configured/default | Maximum radial region used for density; `0` uses the full pore radius |
| Density probe floor | 0 Å | Minimum probe radius; `0` disables it |
| Bulk density | measured; fallback 0.0334 Å⁻³ | Reference density, displayed read-only |
| Gaussian KDE | on | Smooth axial water density |
| KDE bandwidth | 1.4 Å | Smoothing bandwidth; `auto` fits each frame and matches CHAP's software default |
| Poisson energy floor | on | Apply a finite sampling-limited barrier to empty bins |
| Zero energy at mouths | on | Shift the mean mouth free energy to zero |
| CHAP mode | off | Apply CHAP-compatible defaults |
| Fix leucine hydrophobicity | off in CHAP mode | Use corrected leucine value rather than exact CHAP-table compatibility |

Available views are Density, Energy, Hydrophobicity, and Per-frame ρ.
Hydration is pore-only and requires explicit water.

### Reporting a free-energy barrier

Bandwidth controls how strongly the water profile is smoothed and can change a
narrow barrier's height. Use `auto` when matching CHAP's software defaults. Use
a fixed bandwidth when consistent spatial resolution across frames or systems
is required. Compare more than one bandwidth when the conclusion depends on a
narrow gate, and always report the selected value. A floor-limited dry-bin
barrier is a sampling-dependent lower bound.

## Ion Flow, permeation, and passability

| Parameter | Default | Definition |
|---|---:|---|
| Ion Flow view | Occupancy % | Aggregate occupancy map, per-molecule Passage tracks (constriction crossings coloured by direction, drawn on top), or Count vs frame (molecules inside the pore per frame) |
| Passage Show (Water only) | Crossings | Crossings draws only molecules that crossed the constriction; Entered draws every molecule that entered the pore |
| Species | All detected | Filter cached observations by detected ion species; All = every ion type, never water. Water (one oxygen per molecule, from the Hydration water selection) is its own entry |
| Shell | 3 Å | Radial region beyond the mean pore wall included in the map |
| Stride | 1 | Trajectory sampling stride for ion analysis |
| Flip Z | off | Reverse plotted pore direction |
| Permeation bulk planes | automatic | Two bulk boundaries along the per-frame pore axis |
| Time between saved frames | blank | Required to convert counts to rates |
| Applied voltage | blank | Required, with time, to estimate conductance from net transferred charge |
| Passability salt | 150 mM NaCl at 37 °C | Bulk-conductivity preset used by geometric conductance readout |
| Passability custom conductivity | 1.9 S/m | Used when Custom is selected |
| Water probe radius | 1.15 Å | Probe radius used in accessible-volume/passability metrics |
| Readout species | Water, K, Na, Ca | Species shown in the on-figure summary |

Ion Passage means an ion entered the near-pore region. Permeation requires a
complete crossing from one bulk region to the other. Bare/hydrated passability
is a steric comparison, not a free-energy calculation.

## Files and performance

| Parameter | Default | Definition |
|---|---:|---|
| Output directory | automatic | Root for persistent frame directories and manifest |
| Save results | on | Keep results instead of using temporary storage |
| Overwrite | on | Replace prior stored frames after GUI confirmation |
| Keep visualization | on | Retain generated VMD molecules/representations |
| Keep input PDB | off | Preserve exact per-frame coordinates submitted to an engine |
| Prebuild surfaces | off | Generate surfaces before they are displayed |
| Parallel jobs | automatic | Normally available cores minus one |
| Surface cache | 20 frames | Number of prepared surfaces retained in memory |
| Mean-profile frame cap | 1000 | Maximum frames used for expensive mean 3D/property work; `0` uses all frames |

Configure executable paths and the radius file in **File → Settings**. See
[Installation](installation.md) for recommended binaries.

## Minimum reporting set

For a reproducible pore analysis, report the structure/trajectory, frame range
and time stride, periodic imaging and alignment, atom selection, radius file,
`CPOINT`, `CVECT`, method, `SAMPLE`, `ENDRAD`, seed policy, and software/binary
versions. Add dot density and Monte Carlo parameters when changed.

For tunnel analysis, also report origin selection, probe/interior/bottleneck
thresholds, weight function, exit constraints, within-frame clustering,
cross-frame matching parameters, and the **Seen** threshold. For hydration or
ion analyses, report selections, temperature, density reference, binning and
smoothing, trajectory time spacing, bulk planes, voltage, and PBC treatment as
applicable.
