# Tunnel tutorial: cholesterol oxidase

This tutorial uses the bundled cholesterol oxidase structure, PDB 1MXT, to find
routes from a known buried point and then compare them with a whole-structure
origin search.

## Objective

Run one explicit-origin tunnel search, inspect the candidate routes in VMD, and
export the selected route with the inputs needed to repeat the analysis.

## 1. Load and select the structure

Load `vmdhole/1MXT.pdb` in VMD. Open VMDHole, select **Tunnel**, and set:

| Control | Value |
|---|---|
| Molecule | the loaded 1MXT molecule |
| Selection | `protein` |
| Frames | `now` |

The selection defines the molecular interior used by the search. Add a cofactor
or bound ligand only when it should form part of that boundary.

## 2. Set and verify the origin

Leave **Auto-detect origins** off. Click the **Start point** field and paste
`20.4632 0.4374 17.4692`.

Open the **MOLE parameters** gear and enable **Show cues**. Confirm that the
origin marker lies inside the protein. For another system, enter coordinates,
use **COG** with a selection around the buried site, or use **COR** after
centering the VMD view on the site.

## 3. Run the tunnel search

Keep the initial search parameters and select **Run Tunnel**. When the run
finishes, the route table lists the accepted candidates. Select a row to display
that route in VMD.

Inspect each candidate before using its measurements:

- the route begins near the requested origin;
- its centreline remains in accessible internal space;
- it reaches the molecular surface;
- its bottleneck lies on the displayed route rather than an unrelated groove.

If the routes start from the wrong cavity, correct **Start point** or the atom
selection. If no route survives, inspect the origin first and then review the
probe, interior, bottleneck, and minimum-length controls in the
[Tunnel workflow](tunnel-mode.md).

## 4. Inspect geometry and lining

Use the route table to compare rank, bottleneck radius, length, property
summaries, and occurrence. Sorting the table changes its display order but not
route identity.

Select **Lining…** for the current route and enable **Show lining** to display
contacting residues and HET groups. Use the route gear to change its
representation, material, or property coloring without rerunning the search.

## 5. Compare automatic origins

Enable **Auto-detect origins (scan whole structure)** and rerun. The explicit
start point is disabled because it is not used in this mode. VMDHole now scans
qualifying internal cavities and can return routes unrelated to the original
site.

Treat this as a separate analysis question. Compare the starting regions and
3D paths, not only the route ranks. Disable automatic origins before returning
to a site-specific search.

## 6. Export and record the analysis

Export the selected profile figure and CSV from the active plot tab. Use
**Lining… → Export** to save lining data. Retain the saved run or manifest and
record:

- structure, frame, and atom selection;
- explicit start point or automatic-origin setting;
- probe, interior, origin-radius, bottleneck, and length thresholds;
- weight function, exit constraints, and clustering settings;
- VMDHole and tunnel-engine versions.

For a trajectory, make the structure whole and image it consistently. **Align
trajectory** is on by default in Tunnel mode; leave it on for cross-frame
matching unless the trajectory is already aligned. It fits the loaded frames in
VMD. **Seen** is the fraction of analysed frames assigned to the tracked route.

## Citations

Before publishing, open **Help → Guide & Citations… → Citations**. Cite VMD,
VMDHole, HOLE, and MOLE 2 for this workflow. **Cluster within frame** is on by
default, so also cite CAVER 3.0 unless you disable clustering. See
[References](references.md) for the full entries.

Continue with the [Tunnel workflow](tunnel-mode.md).
