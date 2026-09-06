# `native/nm/` — Nelder-Mead pore search and marching-cubes mesher

Two tools, each a single translation unit:

| tool | what it does |
|---|---|
| `nm_search` | HOLE-style pore profile by per-slice Nelder-Mead instead of Monte Carlo annealing |
| `mesh_csg`  | signed-distance marching cubes for the pore surface: spherical, capsule (`--axis`) and tunnel routes |

`nm_threads.h` gives both, and `conn_lobes`, one thread policy: physical cores, not
hyperthreads, unless `OMP_NUM_THREADS` says otherwise. That is not a
micro-optimisation — on the pore search this machine measured 8 threads at
21 ms and 16 at 368 ms, a 17x cliff, and OpenMP's own default is the 16.

`NM_CELL` overrides the atom-grid cell size in the search (default 4 Å).
Swept on an 18,677-atom channel it changed nothing measurable (350-437 ms
across 3-20 Å at the time, before the thread fix), so it is exposed for
future systems rather than because it pays here.

`NM_MAXSTEP` (default 6 Å) caps how far one call to the per-slice optimizer
may move from its own seed. Without it, a vestibule opening into bulk solvent
is a nearly-flat landscape that plain Nelder-Mead runs away across - every
reflect/expand finds a genuine, if marginal, improvement, so the value-based
`nm_cap` stop never fires while the position drifts arbitrarily far before it
does. Measured on the validated 10-frame corpus: 5 of 10 frames jumped 20-39
Å in a single slice-to-slice step at the mouth, visible in the plugin's Pore
Profile plot as a flat run then a vertical cliff. The default keeps every
bottleneck (min radius) value from that corpus unchanged and reproduces the
deep-pore track exactly except for a handful of sub-0.001 Å rows next to a
correction; see VALIDATION.md.

The other `NM_*` variables (`NM_PROBE_SCALE`, `NM_PROBE_ROUNDS`, `NM_JUMPTOL`,
`NM_CSAMP_MULT`, `NM_FINE_GATE`, `NM_FINE_SCALE`, `NM_SEQ`, `NM_SEQ_PROBE`,
`NM_DEBUG_COARSE`) are investigation instruments; every default reproduces
the validated path byte for byte. `NM_MAXSTEP` is the one exception - its
default is load-bearing, not a no-op.

`nm_search` and `mesh_csg` are wired into the plugin (Search picker under the
HOLE parameters; Surface mesher under Settings > Engines).
Both, and `conn_lobes`, are also compiled INTO `sos_triangle_fast`
(`-DVMDHOLE_MULTICALL`, see `native/build.sh`) and reached as
`sos_triangle --nm-search|--mesh|--conn-lobes ARGS`, so one shipped file
carries them all; the plugin looks there when no standalone build sits
beside it.
`VALIDATION.md` has the numbers and the investigation: the search is 5x
faster than HOLE and lies within HOLE's own seed-to-seed spread on all ten
frames tested - the one frame that looked wrong is one where HOLE itself is
bimodal - and the Connolly pass is byte-identical to HOLE's on the same
centres.

## `nm_search` as the plugin runs it

```
nm_search PDB RAD cx cy cz vx vy vz SAMPLE ENDRAD [options]
  --sph FILE        HOLE-layout hole_out.sph (holcal.f order, ADDEND spheres)
  --tsv FILE        the 8-column hole_profile.tsv the plugin parses
  --conn PROBE GRID HOLE's Connolly pass per slice (GRID 0 = 0.7*PROBE)
  --ignore R1,R2    HOLE's IGNORE card
  --centres FILE    skip the search: "t x y z r" centres (t from CPOINT)
  --quiet           no slice listing on stdout
```

Without options it prints the slices as before. `nm_holeout.h` holds the
emitters and the Connolly port (`concal.f`/`coarea.f`/`addend.f`, taken line
for line from the plugin's bit-exact Tcl port, including HOLE's `dz+dz`
duplicate-test bug). Two exact accelerations sit under it: the duplicate test
is a world-xy hash bounded by the slice's highest z, and the area raster is a
block fill that only runs the per-cell loop on cells a circle can reach, each
grey cell keeping that candidate list for the refinement cycles. Both were
checked byte for byte against the plain loops (108k dots, 315 slices) and
against the Tcl port on six hand-picked slices.

Plugin routes: Spherical + Nelder-Mead runs `nm_search` alone; Connolly with
the fast surface runs it with `--conn`, on its own centres under Nelder-Mead
or on HOLE's (`--centres`, extracted from HOLE's spherical `.sph`) under Monte
Carlo. Capsule always runs HOLE. A control file with cards the engine cannot
take falls back to HOLE with a note.

## `mesh_csg` as the plugin runs it

```
mesh_csg SPH OUT.vmd_plot VOXEL [--tri OUT.tri] [--draw|--dots]
mesh_csg --recolor IN.vmd_plot OUT.vmd_plot COLOUR-OPTIONS
mesh_csg --serve      # one request per stdin line:
                      #   mesh<TAB>SPH<TAB>PLOT<TAB>VOXEL      (addressed form)
                      #   meshdraw<TAB>...                     (--draw form)
                      #   meshdots<TAB>...                     (--dots form)
                      #   recolor<TAB>IN<TAB>OUT<TAB>opt<TAB>val...
                      #   extent<TAB>SPH                       -> "EXT xlo xhi ylo yhi zlo zhi rmax"
                      #   tunnelcluster<TAB>IN<TAB>OUT<TAB>THRESHOLD MAXDEV   (one shipped binary only:
                      #   tunneldist<TAB>IN<TAB>OUT<TAB>WANTMAX     sos_triangle's tunnel kernels, no fork of VMD)
                      # replies "OK ntri" or "ERR ...", "quit" ends it
VOXEL = H             # uniform grid
      | H/HF[/R]      # rectilinear: HF cells in every H cell that touches a
                      # pore sphere of radius < R (default 3H), H elsewhere
--draw                # write HOLE's own "draw ..." records instead of
                      # addressing a molecule - what sos_triangle's
                      # per-triangle recolour modes read
```

Reads `hole_out.sph` directly (occupancy = radius; beta > 0 = pore sphere,
beta = 0 = ADDEND clip sphere), samples `f(p) = max(sdf_pore, -sdf_clip)` on
the grid over the pore spheres, runs marching cubes at iso 0, and drops
triangles the clip term owns, so the mouths stay open exactly as
`sph_process` leaves them. The `.vmd_plot` is in `sos_triangle`'s layout,
with the HOLE red/green/blue groups keyed on the owning sphere's radius.
The refined grid is one grid, so it has no seams, and its origin is snapped
to the HF lattice. Edge crossings are not interpolated from the sampled
field: the surface is a union of spheres, so the crossing on a cell edge
lies exactly on the sphere governing the inside corner, and the quadratic
for it is solved directly. That halves the vertex error at every spacing
and is what lets one mesh (the plugin sends `1.0/0.5`) serve both playback
and the settled view instead of a coarse mesh being replaced by a fine one.
Records are written as `graphics $::VMDHole::_gmol ...` rather than `draw`,
so the plugin can replay a mesh by sourcing the file: `draw` is a Tcl proc
that re-resolves the top molecule on every call, which costs more per
triangle than the drawing. The plugin's parser accepts either form.
The plugin keeps one `--serve` child per session and sends a line per frame,
for the spherical surface and for tunnel routes alike (a route is a sphere
union too);
`--with SPH...` names the other frames of a smoothing window: their fields
are averaged with the centre frame's on one grid and the mean is marched, so
the result is a local average of the surfaces (caps and colour bands come
from the centre frame). The legacy pair gets the same window through
`sos_triangle --sos-smooth OUT RHO CENTRE.sos WITH.sos...`, which moves each
dot of the centre cloud to the mean of itself and its nearest same-facing
dot within RHO in every window cloud.
A capsule run passes `--axis CPOINT CVECT ENDRAD`: HOLE's capsule search
leaves escaped slices in the `.sph` with nothing marking them, and the mesher
keeps a slice only within ENDRAD of the axis, the rule HOLE's profile applies.
sos_triangle serves the lateral and lobe colouring of Connolly clouds and any
frame the mesher fails.
Numbers in VALIDATION.md.

## Colouring, dots and the centreline extent

`--recolor` paints a draw-form mesh by property with the rule sos_triangle's
`--recolor --hydro3d-atoms` applies: pore-lining residues (any atom within
`--thresh` of the nearest centreline sphere's surface; with `--facing 1` the
residue's centre of geometry must be nearer that sphere centre than its CA)
contribute at their centre of geometry, a Gaussian kernel of `--bandwidth`
averages them at each triangle centroid, and the seven-step ramp is applied
to value/range (`--signed 1`) or (value-lo)/(hi-lo). `--values FILE` takes
the contributors as given (`x y z value` rows). Contributor order and the
centroid text match sos_triangle's, and on a 7,624-triangle mesh with the
plugin's own sidecar every triangle came out the same colour. Inside the
plugin this runs in the persistent process: 8-11 ms against 44-141 ms for a
spawned sos_triangle, and the mesh is read once.

`--atoms` also accepts a PDB written by VMD with the residue value in the
B-factor column and 1.0 in occupancy on the CA atoms - one C-speed writepdb
instead of a Tcl loop over every atom (17 ms against 76 on a 11k-atom box).

`--dots` emits one `draw point` per distinct mesh vertex, in HOLE's colour
bands, for the dots display. `extent` (serve only) returns the centreline's
bounding box and largest radius; the plugin's pore-lining sidecar used to
find those by reading the whole .sph in Tcl, 100+ ms on a Connolly frame.
