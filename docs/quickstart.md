# First pore and tunnel analyses

Use the bundled five-model gramicidin A structure to verify the installation
and produce a pore profile, then use the bundled cholesterol oxidase structure
to find tunnels from internal cavities.

## Pore mode

## 1. Load the structure

In VMD, select **File → New Molecule…**, browse to
`vmdhole/1GRM.pdb`, and select **Load**. Then open
**Extensions → Analysis → VMDHole** and select **Pore**.

## 2. Set the input

Set the following controls:

| Control | Value |
|---|---|
| Molecule | the loaded 1GRM molecule |
| Selection | `all` |
| Frames | `now` |

The `all` selection includes gramicidin's nonstandard terminal residues. For
another system, start with `protein` and add any cofactor, ligand, or blocker
that forms part of the pore wall. Exclude bulk water, membrane atoms, and freely
moving ions unless they are intentionally part of the obstruction.

Open **File → Settings** and confirm that **Radius file** points to a readable
HOLE `.rad` file. Every selected atom, including nonstandard residues and
cofactors, must match a radius rule. Do not use a result if the VMD console
reports missing radii; choose or edit an appropriate radius file first. See
[Radius files](parameters.md#radius-files).

## 3. Check the pore axis

`CPOINT` should lie inside the channel, and `CVECT` should point along it. When
a molecule is activated, VMDHole proposes suitable values from the current
selection and structure. Keep them for this first run if the point and direction
look correct.

Open the **HOLE parameters** gear and enable **Show cues**. Confirm that the
`CPOINT` marker is inside the channel and the `CVECT` arrow follows its long
axis.

To define the gramicidin axis explicitly:

1. Select **Vector** beside `CVECT`.
2. Enter `resname ETA and chain A` for **Point 1**.
3. Enter `resname ETA and chain B` for **Point 2**.
4. Select **Compute**, inspect the displayed direction, and select **Apply**.

The vector sign changes the profile orientation, not the pore being analysed.
For another structure, use **Guess**, **Use Z**, two selections marking the pore
ends, or two coordinates.

## 4. Run the current conformation

Keep the initial spherical-method defaults and select **Run HOLE**. A completed
run should show:

- a radius profile in **Pore Profile**;
- the pore centreline or surface in the VMD display;
- the minimum-radius readout for the current frame.

Inspect the 3D result before using the numbers. The centreline must pass through
the intended channel rather than an external groove. If it does not, correct
`CPOINT`, `CVECT`, or **Selection**, then rerun.

## 5. Compare all five conformers

Set **Frames** to `all` and select **Run HOLE** again. Use:

- **Over Time** to compare position across conformers;
- **Mean Profile** to inspect the mean and spread;
- **Trends** to compare minimum radius, volume, or conductance estimates;
- **Histogram** to inspect the distribution of per-frame radii.

These five experimentally determined NMR conformers show structural variation
and the corresponding variation in pore geometry. Unlike consecutive MD
frames, however, their order does not encode elapsed time.

## Tunnel mode

### 1. Load the structure

Load `vmdhole/1MXT.pdb` in VMD and select **Tunnel** in VMDHole. Set:

| Control | Value |
|---|---|
| Molecule | the loaded 1MXT molecule |
| Selection | `protein` |
| Frames | `now` |

### 2. Find and inspect tunnels

Enter `20.4632 0.4374 17.4692` under **Start point**. Open the **MOLE
parameters** gear, enable **Show cues**, confirm that the marker lies inside the
protein, and select **Run Tunnel**.

The route table lists the tunnels found from qualifying internal cavities.
Select a route and confirm in the VMD display that it begins inside the
structure and reaches the molecular surface. Use **Lining…** to inspect its
contacting residues. Route rank is a search score; visual inspection remains
part of validating a candidate pathway.

## Next steps

Choose the [Pore or Tunnel tutorial](tutorials.md) for the mode you just ran.
Afterward, use the matching [workflow guide](workflows.md) when you need more
than the tutorial settings.
