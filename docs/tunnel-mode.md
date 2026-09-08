# Tunnel workflow

<p align="center"><img src="images/tunnel_3d_side.png" alt="Four tunnels found by the MOLE 2 engine, branching from a buried start point" width="480"></p>

Tunnel mode searches for routes from a buried site to the molecular surface.
Use it for internal cavities, enzyme access paths, and branched egress routes.
Use pore mode instead when one known channel axis is the object of study.

For a publication using Tunnel mode, cite MOLE 2. If route clustering is used,
also cite CAVER 3.0; both references are listed in
[References](references.md#tunnel-analysis).

## 1. Prepare and align

Load the structure or trajectory and choose the atoms that define the molecular
interior. Make periodic structures whole before analysis. For trajectory-wide
route tracking, **Align trajectory** is on by default. Keep it enabled and
choose a stable reference selection unless the trajectory is already aligned.
It fits the loaded frames before searching. Cross-frame clustering compares
route geometry; translation or rotation of an unaligned protein is otherwise
interpreted as route motion.

## 2. Define the origin

The start point should lie in the buried cavity of interest. Enter coordinates,
use a selection's centre of geometry (**COG**) or VMD's centre of rotation
(**COR**), or enable automatic origin detection. A poor origin can return no
routes or routes from the wrong cavity. The **⌖** button opens the on-screen
stick described under pore mode, step 3, to nudge the start point relative to
the current view.

The start point also accepts a **VMD selection** instead of coordinates — its
centre is re-evaluated in every frame, so a residue-defined origin (for example
`resname HEM`) follows the trajectory rather than staying where it was in frame
0. Several origins are given as a `;` list mixing both forms
(`73.8 26.5 26.6; resid 74 and chain A`); each is pinned separately and the
routes found from all of them are merged and de-duplicated, which is MOLE's
pinned multi-origin mode.

Custom exits restrict the search toward known surface regions, and take the same
three forms — coordinates, a selection, or a `;` list of them. A custom path is
defined by start and end points. **Use custom exits only** excludes other exit
candidates; use it only when the biological exit is independently known.

## 3. Set the search criteria

The six controls shown in the main panel are sufficient for most analyses:

| Control | Default | Meaning |
|---|---:|---|
| Probe | 3.0 Å | Probe used to define accessible void space |
| Interior | 1.25 Å | Minimum interior clearance |
| Origin radius | 5.0 Å | Region around the requested start used to seed origins |
| Minimum length | 0 Å | Reject routes shorter than this value |
| Bottleneck | 1.25 Å | Minimum accepted route radius |
| Cluster within frame | on | Merge geometrically similar routes in each frame |

The remaining search, exit, clustering, and rendering controls are available
from the **⚙** next to **MOLE parameters**. They are defined in the
[parameter reference](parameters.md#tunnel-search).
Change one class of parameters at a time and retain the settings with exported
results.

## 4. Run and inspect routes

<p align="center"><img src="images/tunnel_panel.png" alt="Tunnel panel: per-route bottleneck, length, hydrophobicity and charge, with the selected route profile" width="860"></p>

Select **Run Tunnel**. The route table reports:

| Column | Meaning |
|---|---|
| Show | Visibility |
| route/color | Tracked route identifier and display color |
| Rts | Number of route instances represented by the cluster |
| Bneck | Mean bottleneck radius |
| Len | Mean route length |
| Phob | MOLE length-weighted hydrophobicity mean |
| Chg | Mean net formal charge |
| Seen | Percentage of analysed frames containing the tracked route |

Sort by a column to inspect a different property; sorting does not alter route
identity. Expand a row for details. The row gear controls that route's
representation, color, material, and property. The global gear applies display
choices to routes without a per-route override.

Route surfaces are meshed by `mesh_csg` (Settings → Engines → Spherical
mesher), the same marching-cubes mesher the spherical pore uses, since a
route is a union of spheres along its centre line. A route coloured by a
property is meshed by `sos_triangle` instead, because the per-triangle
recolouring reads that program's own mesh records. On a four-route frame of
KcsA, meshing and drawing took 59 ms against 217 ms for `sos_triangle`.

Tunnel properties are Kyte–Doolittle, Wimley–White, Kapcha–Rossky,
Fauchère–Pliska, and the MOLE hydropathy, hydrophobicity, polarity, charge,
ionizable, logP, logD, logS, and mutability fields. See
[Properties](properties.md) for their definitions and citations.

First validate every candidate in 3D. A high-ranked route can still be an
irrelevant solvent-accessible groove, and a route that terminates incorrectly
usually indicates an origin, selection, exit, or interior-classification issue.

## 5. Cluster and track

**Cluster within frame** merges similar candidates produced in one structure.
The bottleneck-row clustering control sets its geometric cutoff.

Cross-frame clustering assigns a persistent route identity to matching routes
from aligned frames. It is controlled by maximum geometric deviation, maximum
ranks considered per frame, and the minimum **Seen** percentage. Restricting
ranks reduces cost but can hide a route that is poorly ranked in some frames.

Treat a low-Seen cluster cautiously in Mean Profile or trend plots: the average
may describe only a small subset of frames. Cross-frame identity is a geometric
classification, not proof that individual solvent molecules use the route.

## 6. Display lining residues

Select a route and open **Lining** to inspect protein residues and HET groups in
contact with it. **Show lining** creates a VMD representation; previous/next
controls step through routes. **Show all** displays the routes that remain after
the current filters.

The lining window exports the selected route's lining data. The standard plot
tabs and CSV exports operate on the selected tracked route. The exported
`tunnel_N.csv` carries per-point `FreeRadius` and `BRadius` alongside the
radius, and the lining window reports the positive and negative residue counts
next to the net charge.

## 7. Available downstream analyses

Tunnel mode supports the radius/profile, Over Time, Mean Profile, Trends,
Histogram, property, lining, and Ion & Water views. Ion & Water measures the
selected route along itself, as distance along the route and distance from
it, so a bent tunnel plots as it is. It does not provide
tunnel hydration, tunnel ellipse fitting, or pore-mode bulk-to-bulk permeation.
Water free-energy and density properties require a pore-mode hydration result.

## 8. Cavities

A **cavity** is the pocket itself - a volume - as distinct from a route, which
is a path out of one. Routes start inside cavities, so the two are views of the
same search: the room and the corridor leaving it.

**Cavities** lists what MOLE found in the displayed frame: type (*Cavity*, or
*Void* when nothing lines it), this frame's volume, the mean and spread over the
frames the cavity was tracked through, how often it was seen, its **max probe**
(the largest sphere that fits inside - whether your ligand fits at all), depth,
and the boundary/inner residue counts. Click a column header to sort. **Residues**
opens the two residue sets with their MOLE properties and a ready-made VMD
selection string.

### Using a cavity to start a search

The hard part of a tunnel run is choosing the origin: a poor one returns no
routes, or routes from the wrong cavity. **Use as start** puts a cavity's start
point into the **Start point** field, so the next run searches from that pocket.

The natural loop is therefore: run once with **Auto-detect origins** (which needs
no start point and returns every cavity), inspect the list, pick the pocket you
care about, then re-run pinned to it.

Two rules are offered, and the status line always says which one was used:

| Rule | Point | Source |
|---|---|---|
| **deepest (MOLE)** | the cavity's deepest point by `DepthLength` | MOLE's own automatic origin, read from the engine rather than recomputed - the point it would have searched from itself |
| **largest sphere (CAVER)** | centre of the largest sphere that fits inside | the rule CAVER Analyst's *Create Starting Point* uses |

Under MD the more robust origin is often not a coordinate at all but a **VMD
selection** in the Start point field, which is re-evaluated every frame and so
follows the protein; the cavity's residue list is a good place to find one.

### Display

Ticking **Draw** meshes the cavity as a sphere-union surface on a **track of its
own**, separate from the routes, so route display settings and cavity ticks do
not disturb each other. **Solid** switches from transparent to opaque, and
**Spheres** draws the clearance spheres themselves instead of a surface over
them. **Show all** / **Hide all** apply to every cavity in the frame.

### What the numbers mean, and what they are not

*Id* is a **tracked** id: the same pocket keeps the same number, the same
colour and the same tick in every frame, matched across frames by centroid
proximity. This matters because MOLE recomputes cavities independently in each
frame and ranks them by volume, so two pockets that swap volume order swap
ranks - keying anything on the rank would mean a ticked cavity silently became
a different pocket on the next frame. That per-frame rank is still shown, in
**Rank here**, and a row reading *absent* is a pocket the displayed frame does
not have (check its *Seen %*). The tracking is this plugin's own: neither MOLE
nor CAVER reports cavity behaviour over a trajectory at all.

Two caveats worth carrying into a figure caption:

- The drawn surface is a **marching-cubes sphere union**. MOLE triangulates the
  boundary facets of its tetrahedra with the corners at *atom centres*, and
  CAVER Analyst renders an analytic solvent-excluded surface. These are three
  different surfaces of the same pocket, so areas and volumes are not
  interchangeable between the programs.
- MOLE's own *Volume* column is not the volume of anything it draws either: it
  is the tetrahedra minus van der Waals caps.

Cavities require the compiled engine; a run made with the pure-Tcl fallback has
none, and the window says so.
