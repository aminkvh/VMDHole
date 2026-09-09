# Citations and acknowledgements

Use this page when preparing a manuscript, figure caption, data release, or
software acknowledgement. Cite VMDPathFinder and VMD for every analysis. Add
HOLE for Pore mode - every engine and mesher in it implements HOLE's method -
and then only the entries for the other methods or quantities you report. A
Tunnel-mode run that never touches Pore mode does not use HOLE and should not
cite it; the in-plugin citation guide (Help > About > Citations) states the same
rule and is the one to follow.

## Cite for every VMDPathFinder analysis

1. Ahangar, A. A. *VMDPathFinder* (version used) [computer software].
   <https://github.com/aminkvh/VMDPathFinder>
2. Humphrey, W., Dalke, A. & Schulten, K. "VMD: Visual Molecular Dynamics."
   *J. Mol. Graph.* **14**, 33-38 (1996).
   doi:[10.1016/0263-7855(96)00018-5](https://doi.org/10.1016/0263-7855(96)00018-5)
3. Smart, O.S., Neduvelil, J.G., Wang, X., Wallace, B.A. & Sansom, M.S.P.
   "HOLE: A Program for the Analysis of the Pore Dimensions of Ion Channel
   Structural Models." *J. Mol. Graph.* **14**, 354-360 (1996).
   doi:[10.1016/S0263-7855(97)00009-X](https://doi.org/10.1016/S0263-7855(97)00009-X)

Use the VMDPathFinder version shown in **Help → Guide & Citations…**. If a versioned
release DOI is available, use it in place of the repository URL.

## Add citations for the analyses used

### Pore geometry and conductance

- **Connolly surface:** Connolly, M.L. "Analytical Molecular Surface
  Calculation." *J. Appl. Crystallogr.* **16**, 548-558 (1983).
  doi:[10.1107/S0021889883010985](https://doi.org/10.1107/S0021889883010985)
- **Ellipse fit, ellipse volume, or ellipse conductance:** Seiferth, D. &
  Biggin, P.C. "Exploring the Influence of Pore Shape on Conductance and
  Permeation." *Biophys. J.* **123**, 3107-3119 (2024).
  doi:[10.1016/j.bpj.2024.07.010](https://doi.org/10.1016/j.bpj.2024.07.010)
- **Geometric conductance estimate:** Smart, O.S., Breed, J., Smith, G.R. &
  Sansom, M.S.P. "A Novel Method for Structure-Based Prediction of Ion Channel
  Conductance Properties." *Biophys. J.* **72**, 1109-1126 (1997).
  doi:[10.1016/S0006-3495(97)78760-5](https://doi.org/10.1016/S0006-3495(97)78760-5)
- **Access resistance correction:** Hall, J.E. "Access Resistance of a Small
  Circular Pore." *J. Gen. Physiol.* **66**, 531-532 (1975).
  doi:[10.1085/jgp.66.4.531](https://doi.org/10.1085/jgp.66.4.531)
- **Passability using hydrated or bare ion radii:** Nightingale, E.R.
  "Phenomenological Theory of Ion Solvation. Effective Radii of Hydrated Ions."
  *J. Phys. Chem.* **63**, 1381-1387 (1959).
  doi:[10.1021/j150579a011](https://doi.org/10.1021/j150579a011); Pauling, L.
  *The Nature of the Chemical Bond*, 3rd ed. (1960).

### Hydration and pore-wall annotations

- **Hydration, water-density free energy, CHAP mode, or hydrophobic gating:**
  Klesse, G., Rao, S., Sansom, M.S.P. & Tucker, S.J. "CHAP: A Versatile Tool
  for the Structural and Functional Annotation of Ion Channel Pores."
  *J. Mol. Biol.* **431**, 3353-3365 (2019).
  doi:[10.1016/j.jmb.2019.06.003](https://doi.org/10.1016/j.jmb.2019.06.003)
- **VMDPathFinder pore-lining or pore-facing assignments:** cite VMDPathFinder. These
  surface-distance and orientation assignments are not PoreWalker
  classifications. Cite PoreWalker only when making a direct comparison:
  Pellegrini-Calace, M., Maiwald, T. & Thornton, J.M. "PoreWalker."
  *PLoS Comput. Biol.* **5**, e1000440 (2009).
  doi:[10.1371/journal.pcbi.1000440](https://doi.org/10.1371/journal.pcbi.1000440)
- **Automatic hydration KDE bandwidth:** Silverman, B.W. *Density Estimation
  for Statistics and Data Analysis* (1986).
- **Dry-bin lower bounds based on the rule of three:** Hanley, J.A. &
  Lippman-Hand, A. "If Nothing Goes Wrong, Is Everything All Right?"
  *JAMA* **249**, 1743-1745 (1983).

### Tunnel analysis

- **Tunnel mode:** Sehnal, D., Svobodová Vařeková, R., Berka, K., Pravda, L.,
  Navrátilová, V., Banáš, P., Ionescu, C.-M., Otyepka, M. & Koča, J.
  "MOLE 2.0: Advanced Approach for Analysis of Biomacromolecular Channels."
  *J. Cheminform.* **5**, 39 (2013).
  doi:[10.1186/1758-2946-5-39](https://doi.org/10.1186/1758-2946-5-39)
- **Tunnel clustering:** Chovancova, E., Pavelka, A., Benes, P., Strnad, O.,
  Brezovsky, J., Kozlikova, B., Gora, A., Sustr, V., Klvana, M., Medek, P.,
  Biedermannova, L., Sochor, J. & Damborsky, J. "CAVER 3.0: A Tool for the
  Analysis of Transport Pathways in Dynamic Protein Structures." *PLoS Comput.
  Biol.* **8**, e1002708 (2012).
  doi:[10.1371/journal.pcbi.1002708](https://doi.org/10.1371/journal.pcbi.1002708)

### Physicochemical maps

Cite a scale when its values are reported or used to support an interpretation.

| Scale | Reference |
|---|---|
| Kyte-Doolittle | Kyte, J. & Doolittle, R.F. "A Simple Method for Displaying the Hydropathic Character of a Protein." *J. Mol. Biol.* **157**, 105-132 (1982). doi:[10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0) |
| Wimley-White | Wimley, W.C. & White, S.H. "Experimentally Determined Hydrophobicity Scale for Proteins at Membrane Interfaces." *Nat. Struct. Biol.* **3**, 842-848 (1996). doi:[10.1038/nsb1096-842](https://doi.org/10.1038/nsb1096-842) |
| Kapcha-Rossky | Kapcha, L.H. & Rossky, P.J. "A Simple Atomic-Level Hydrophobicity Scale Reveals Protein Interfacial Structure." *J. Mol. Biol.* **426**, 484-498 (2014). doi:[10.1016/j.jmb.2013.09.039](https://doi.org/10.1016/j.jmb.2013.09.039) |
| Grantham polarity | Grantham, R. "Amino Acid Difference Formula to Help Explain Protein Evolution." *Science* **185**, 862-864 (1974). doi:[10.1126/science.185.4154.862](https://doi.org/10.1126/science.185.4154.862) |
| Fauchere-Pliska lipophilicity | Fauchère, J.L. & Pliska, V. "Hydrophobic Parameters Pi of Amino-Acid Side Chains." *Eur. J. Med. Chem.* **18**, 369-375 (1983). |
| Formal charge | Fersht, A.R. *Structure and Mechanism in Protein Science* (1999). |

MOLE-derived tunnel properties, including hydropathy, hydrophobicity, polarity,
charge, ionizable residues, LogP, LogD, LogS, and mutability, are covered by
the MOLE 2.0 citation above.

### Ion movement and trajectory alignment

- **Voltage-driven computational electrophysiology trajectories:** Kutzner, C.,
  Grubmüller, H., de Groot, B.L. & Zachariae, U. "Computational
  Electrophysiology." *Biophys. J.* **101**, 809-817 (2011).
  doi:[10.1016/j.bpj.2011.06.010](https://doi.org/10.1016/j.bpj.2011.06.010)
- **Align trajectory:** acknowledge the VMD RMSD Trajectory Tool and Giorgino,
  T. "Computing 1-D Atomic Densities in Macromolecular Simulations."
  *Comput. Phys. Commun.* **185**, 1109-1114 (2014).
  doi:[10.1016/j.cpc.2013.11.019](https://doi.org/10.1016/j.cpc.2013.11.019)

## Acknowledgements

The screen-space hydrophobicity scale bar adapts techniques from VMD's
`colorscalebar.tcl`, by Wuwei Liang, Dan Wright, John Stone, and Axel
Kohlmeyer.

## Where each idea comes from

Tunnel mode ports MOLE 2's algorithm, and its cavity view borrows interface
ideas from both MOLE and CAVER. This table records which is which, so a method
section can be written without guesswork. Nothing here changes what must be
cited: MOLE for Tunnel mode, CAVER additionally when route clustering is used.

| Feature | Whose idea | Whose algorithm |
|---|---|---|
| Cavity / Void as distinct objects | MOLE | MOLE (Delaunay tetrahedra, depth) |
| Volume and depth columns | MOLE (the only two its own GUI showed) | MOLE |
| Boundary/Inner residues with physicochemical properties | MOLE (present only in its XML, never on screen) | MOLE |
| The drawn cavity surface | both draw one | **this plugin** - marching-cubes sphere union; MOLE uses atom-centre facets, CAVER Analyst an analytic SES |
| "Use as start point" | CAVER Analyst (*Create Starting Point*) | selectable: **MOLE**'s own automatic origin, or **CAVER**'s largest inscribed sphere |
| Colour a cavity by volume / max probe | CAVER Analyst (per-cavity flat colour) | not implemented here |
| Max probe column | CAVER Analyst | max radius over the cavity's own spheres |
| Sortable table, Show/Hide all | CAVER Analyst | interface only |
| Solid / transparent toggle | MOLE (*Solid cavities*) | interface only |
| Spheres view | CAVER Analyst (*Locked Probes*) | the engine's own clearance spheres |
| Snapping a user origin into the cavity | both (MOLE `OriginRadius`, CAVER `rmin`/`dmax`) | MOLE's |
| Reporting the origin actually used | CAVER 3.0 requires it for reproducibility | recorded per frame in the run manifest |
| Cavity tracked across a trajectory | **neither program does this** | this plugin (centroid proximity) |
| Cavity coloured by a property | **neither program does this** | this plugin (per-residue sidecar + the shared recolour kernel) |

Two numerical caveats follow from the third row: our cavity VOLUME is MOLE's own
quantity - the sum of the cavity's tetrahedra minus their van der Waals corner
caps, computed in mole_complex.c and written by the engine - and NOT the volume
of the mesh drawn on screen. The mesh is a marching-cubes sphere union built for
display; its enclosed volume is a different number and is not what any column or
export reports. CAVER Analyst's is a Monte-Carlo estimate over filling balls,
and CAVER Analyst's is a Monte-Carlo estimate over filling balls. They are three
different quantities and should not be compared directly.

### Cavity colouring, specifically

CAVER Analyst has no cavity-owned colouring strategy at all: its colouring
window has tabs for structures, selections and tunnels, and none for cavities.
A cavity there can be given a flat colour per cavity (random, or mapped from its
volume or max probe), or - under *Advanced coloring* - inherit the **protein's
per-atom** colours interpolated over its surface. Colouring a void surface by a
physicochemical property exists in CAVER only for **tunnels**, computed from
atoms within 5 Å of each surface point. MOLE cannot colour a cavity by a
property either; its per-vertex spectrum exporters are unimplemented stubs.

Colouring a cavity surface by a property of its own lining residues is therefore
this plugin's own, and it reuses the same per-residue sidecar and recolour
kernel that colour a route and a pore wall, so the three are one mechanism
rather than three. Neither reference program offers a legend for cavity
colouring; ours follows the scale bar the other property views already use.
