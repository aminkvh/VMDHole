# Visualization

VMDHole creates separate VMD molecules or representations for calculated
geometry. These objects are visual aids tied to numerical results; they do not
modify the analysed structure.

## Pore representations

| Mode | Content | Use |
|---|---|---|
| None | No generated geometry | Batch analysis or fastest plotting |
| Centerline | HOLE sphere centres | Verify the path and endpoints |
| Dots | Surface sample points | Inspect sampling and difficult geometry |
| Wireframe | Mesh edges | Inspect topology without an opaque surface |
| Isosurface | Triangulated pore surface | Presentation and property mapping |
| Ellipse surface or point cloud | Fitted non-circular cross-sections | Inspect slit-like pore shape |

To display the ellipse, choose **Ellipse fit** in **Pore Profile**, then choose
**Solid surface** or **Point cloud** from **Render**. It is available for
Spherical pore analyses, not Connolly or Capsule analyses.

Start with **Centerline**. Interpret profiles or physicochemical maps only after the
path is shown to occupy the intended pore in representative frames.

Connolly and high-density surfaces can exceed the stock HOLE surface
converter's limits. Use the accelerated `sos_triangle`, reduce dot density, or
use Centerline when the mesh is not needed.

## Tunnel representations

Tunnel routes support Isosurface, Wireframe, and Centerline views. Use the
global gear to set defaults and a row gear to override one route. Per-route
controls include visibility, color, property color, material, and
representation. **Show all** applies to routes that pass the current tracking
and Seen filters.

## Color and properties

A pore surface can retain HOLE's radius bands, use a flat VMD color, or use a
supported physicochemical or hydration property. In Connolly mode, `pore_lat`
distinguishes the central pore from lateral extensions. Select `pore_lobes` to
open the inline region table, where each tracked opening can be shown, colored,
annotated, and exported independently. See
[Connolly lateral openings](pore-mode.md#inspect-connolly-lateral-openings) for
the measurements in that table. A tunnel can use automatic route/rank color, a
flat color, or a property. Fixed property-scale limits keep the same meaning
across frames and figures. See [Properties](properties.md) for data sources and
restrictions.

The main surface, Pore Profile Fill, and Mean Profile synchronize a property
where it is available. Over Time is independent and recalculates only after
**Compute**.

## Lining and pore-facing residues

**Lining** creates VMD selections for residues close to the pore or tunnel
surface. The distance is measured from molecular surface to calculated surface,
not from an atom centre to the centreline.

**Pore-facing only** further restricts eligible side chains to those directed
toward the lumen. This filter does not apply to atom-level Kapcha–Rossky values
or water-derived fields.

## Playback and caching

With synchronization enabled, moving the VMD frame slider selects the
corresponding calculated frame. Surface meshes are built lazily and held in an
in-memory cache. Relevant controls are:

- **Prebuild surfaces**: prepare all selected frames before playback.
- **Surface cache**: retain a limited number of meshes.
- **Playback detail**: use a coarser representation during rapid frame changes.
- **Keep visualization**: retain generated VMD objects after a run or reset.

Prebuilding improves interactive playback but increases initial runtime and
disk use. A larger cache reduces rebuilding but increases memory use.

## Mean 3D surface

Mean Profile can revolve the mean radius profile into a trajectory-average 3D
surface. It is not an average Cartesian molecular surface. Property color can
be projected with the fast approximation or **Accurate 3D** mode. **Render
smoothly** applies mesh smoothing to this mean surface and is off by default;
record both choices when comparing figures.

## Scale bar and metrics

The 3D scale bar reports the selected property and its fixed range. Its corner
and font color are presentation settings. The optional metrics readout shows
radius, volume, conductance estimates, and selected steric passability values.
Those annotations do not alter exported numerical data.
