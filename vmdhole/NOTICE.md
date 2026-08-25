# Third-party notices — VMDHole

VMDHole's own code is MIT-licensed (`LICENSE`, repository root). Parts of
`vmdhole.tcl` are derived from HOLE 2 and remain under the Apache License 2.0;
the MIT statement covers VMDHole's own code only.

This file and `LICENSE-Apache-2.0.txt` ship inside `vmdhole/` so that the
licence text travels with the code it covers.

## HOLE 2 — Apache License 2.0

- Upstream: <https://github.com/osmart/hole2>
- Licence text: `LICENSE-Apache-2.0.txt`
- Cite: Smart, O.S., Neduvelil, J.G., Wang, X., Wallace, B.A. & Sansom, M.S.P.
  *HOLE: A program for the analysis of the pore dimensions of ion channel
  structural models.* J. Mol. Graph. **14**, 354–360 (1996).

### Modification notice (Apache 2.0 §4(b))

The inlined pure-Tcl HOLE engine in `vmdhole.tcl` (between the
`BEGIN INLINED HOLE PURE-TCL ENGINE` and `END` markers) is a derived work of
HOLE 2, ported from `holcal.f`, `holeen.f`, `honewp.f`, `hograp.f`, `hcapgr.f`,
`tsatr.f` and the CONNOLLY/CAPSULE routines. It is modified relative to that
source: translated from Fortran 77 to Tcl, restructured to run as a subprocess,
and it refuses control cards it cannot translate rather than ignoring them.

`native/` (not part of the installed plugin) distributes modified HOLE 2
Fortran sources — `hcapen_fast.f`, `coarea_fast.f`, `holcal_par.f`,
`holeen_par.f` — under the Apache licence in `native/LICENSE`.

The plugin also calls unmodified HOLE 2 binaries when installed. No HOLE binary
is included here.

Some ported routines carry legacy Birkbeck College confidentiality headers,
which are also present in the current Apache-licensed upstream. HOLE's author,
Oliver S. Smart, confirmed directly that "the Apache 2.0 licence covers
everything." The headers are retained where they appear in ported code because
they are part of the upstream text, not because their terms apply.

## MOLE 2 — reimplementation of a published method

The tunnel engine reimplements the published method of MOLE 2 (Sehnal et al.
2013; Pravda et al. 2018), validated against MOLE 2's own exported results. No
MOLE source code is included.

## CAVER — cross-snapshot clustering method

Cross-frame tunnel clustering follows the published average-link agglomerative
method of CAVER 3.0 (Chovancova et al., *PLoS Comput Biol* **8**: e1002708,
2012). No CAVER code is included. This attribution must not be removed.
