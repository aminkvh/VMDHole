# Pore tutorial: gramicidin A

This case study develops a reproducible pore analysis of PDB 1GRM. Its purpose
is to establish the workflow: define the system, verify the pore geometry,
analyse an ensemble, and record enough information to repeat the calculation.

## Objective

Measure the gramicidin channel radius in each NMR model and compare the
profiles. The structure contains five models and no explicit water, ions, or
membrane. Hydration, ion-flow, and permeation analyses therefore do not apply.

## 1. Load and inspect

Load `vmdhole/1GRM.pdb`. Display all atoms and confirm that the two
peptide chains form one continuous channel. Select `all` in VMDHole because the
terminal groups are not included by VMD's `protein` macro.

## 2. Define an axis independent of file orientation

Use the terminal groups at opposite channel mouths:

1. Select **Vector** beside `CVECT`.
2. Set point 1 to `resname ETA and chain A`.
3. Set point 2 to `resname ETA and chain B`.
4. Select **Compute**, inspect the displayed vector, and select **Apply**.
5. Set `CPOINT` to the centre of geometry of `all`, or use the midpoint shown by
   the vector tool.

This definition remains meaningful if the coordinates are rotated. If the
profile is displayed in the opposite direction from the desired convention,
reverse the two points or use the plot flip control.

## 3. Run one model and validate the path

Set **Frames** to `now`, keep the spherical method, and select **Run HOLE**.
Switch the 3D representation to **Centerline** and check that:

- the path lies inside the channel;
- both ends reach bulk-facing openings;
- the path does not jump into an external groove;
- the narrow region is consistent with the molecular structure.

If any check fails, do not tune the plot. Correct the selection, `CPOINT`, or
`CVECT` and rerun.

## 4. Make the run reproducible

Open **HOLE parameters** and confirm the random seed. A blank value resolves to
the fixed seed `1`; enter another integer only when the protocol requires it.
Record:

- structure identifier and model or frame range;
- VMD atom selection;
- radius file;
- `CPOINT`, `CVECT`, `SAMPLE`, `ENDRAD`, and `SHORTO`;
- pore method and dot density;
- Monte Carlo seed, steps, step size, and temperature;
- VMDHole, HOLE, and binary versions.

HOLE uses stochastic optimization. Exact reruns and method comparisons require
the same effective seed and otherwise identical inputs.

## 5. Analyse the five-model ensemble

Set **Frames** to `all` and rerun. For comparisons across a genuine trajectory,
first make the protein whole, image it consistently, and align it to a reference
selection. The five 1GRM models are already suitable for this small example.

Use the analysis tabs in this order:

1. **Pore Profile**: inspect individual frames and locate the constriction.
2. **Over Time**: confirm that the same axial region is compared across frames.
3. **Mean Profile**: report mean radius and variation along the pore.
4. **Trends**: compare minimum radius or pore volume per model.
5. **Histogram**: examine the distribution without implying a time sequence.

Ellipse-based quantities require an additional cross-section calculation. The
**Over Time** tab displays a **Compute** gate before this expensive calculation;
select it deliberately and retain the source label in exported results.

## 6. Add physicochemical context

Select a residue property such as Kyte–Doolittle or Wimley–White. The 3D
surface, Pore Profile fill, and Mean Profile property share one selector. Over
Time has an independent property selector and requires **Compute** after a
property or cutoff change.

Use **Pore-facing only** when the question concerns the physicochemical environment
presented to the lumen. Report the property scale, surface-distance cutoff, and
smoothing bandwidth.

## 7. Export

Export the Pore Profile figure and CSV, Mean Profile CSV, and Trends CSV. Keep
the run directory or manifest with the analysis. A figure without its parameter
set, radius file, frame definition, and trajectory preparation is not fully
reproducible.

## Citations

Before publishing, open **Help → Guide & Citations… → Citations**. Cite
VMDHole, VMD, and HOLE. Add citations for any additional method used, such as
Connolly, ellipse fitting, hydration, or a property scale. See
[References](references.md) for the full entries.

Continue with the [Pore workflow](pore-mode.md) for trajectory operation and
additional pore analysis views.
