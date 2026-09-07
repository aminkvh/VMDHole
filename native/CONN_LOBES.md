# `conn_lobes` — fast Connolly lateral-opening classifier

Replaces two pure-Tcl hot loops in the "lobe" coloring path for Connolly
mode (`proc ::VMDPathFinder::_conn_classify_sph`, `_conn_frame_lobes`,
`_split_conn_mesh_by_region`): classifying every `.sph` dot as pore or
lateral and clustering the lateral dots into lobes, and assigning each
triangle of the unified surface mesh to the region its nearest dot belongs
to. Both are read-through-a-file, per-token Tcl loops over tens of
thousands of records; this is the same arithmetic, compiled.

## Measured

18,677-atom Nav trajectory frame, 108,650 dots, 8 lateral lobes (default
margin, no trim/gate filtering — i.e. the untouched default configuration):

| stage | Tcl | native | speedup |
|---|---:|---:|---:|
| classify + cluster | 940-1188 ms | 40 ms | ~25-30x |
| triangle-to-region split (115,169 triangles) | 2,092-2,120 ms | 112 ms | ~19x |
| **full pipeline** (classify, cluster, mesh, split, 9 regions) | 4,690-4,753 ms | 1,719-1,728 ms | **2.7x** |

The remaining ~1.7 s in the native run is `sph_process`/`sos_triangle`
themselves (unchanged, native, building the actual surface geometry) — the
Tcl overhead this tool removes was 3.0-3.5 s of every first-time view of a
Connolly frame with lateral-opening coloring on.

1GRM (575 pore dots, 267 lateral dots, 3 lobes): full pipeline 82-94 ms
either way — small systems were never the problem; the win scales with dot
count.

## Correctness

`test_conn_lobes_engine` (vmdpathfinder/tests/) runs a real HOLE Connolly pass on
1GRM and checks, native vs the Tcl path it replaces:
- pore/lateral/keep dot SETS byte-identical
- escaped-range boundaries agree to 1e-6
- every lobe's z/azimuth/count/escaped-fraction/neck values and dot
  membership agree exactly
- the full region-mesh pipeline's output surfaces are triangle-set
  byte-identical, region by region

Also checked directly on the 108,650-dot Nav fixture: pore/lateral dot sets,
escaped ranges, and all 8 lobes' stats and membership match the Tcl
reference exactly; the union-mesh split assigns the same 77,509 pore /
37,660 lateral triangles as the Tcl split, byte for byte.

## The end-of-channel balls (fixed alongside)

Reported as "blobs at the end of the lateral openings are back". Not the
classifier: HOLE marks its clip spheres (ADDEND records past endrad, escape
dots) with a `LAST-REC-END` line after the record, and `sph_process` draws
a marked sphere only as a cutter that carves the mouth funnel. The
point-cloud reducer the display path runs on clouds over 14,000 points
(`_reduce_conn_sph`) kept only `ATOM` lines, so every marker was lost and
every clip sphere became a drawn one - smooth balls at both ends. That is
why lowering dot density or endrad "fixed" it: below the threshold the
reducer never ran. Measured on a Nav frame (47,349-point cloud): triangles
lying on clip-sphere surfaces 761 with the old reducer, 10 after, 11 in the
unreduced full-resolution mesh.

Fixed at every point the marker could be dropped: the reducer now thins
only the `-999` dots, keeps every centre/escape/clip record, and re-emits
each kept record's marker; the trimmer does the same; the classifier (Tcl
and native, `MARKED` block) records which lines carried one and
`_sph_puts` restores it in every region file. Recipe tags bumped (`_u2`,
`|k1`, `_render2.sph`) so cached surfaces from the old reducer rebuild.
`test_conn_lobes_engine` checks marker parity, region-file restoration and
the reducer. Smooth patches that remain on escaped centre spheres are HOLE's
own geometry (3,161 vs 3,166 triangles at full resolution) - the escape
balloons the Openings list flags at 100% escaped.

## Design notes

- Two subcommands (`classify`, `split`), one binary, matching how the two
  callers use them (`_conn_classify_sph` needs classify+cluster together;
  `_split_conn_mesh_by_region` needs split alone, called after the union
  mesh exists).
- `classify`'s output text reconstructs the exact same dict shape the Tcl
  procs return (`pore`/`lateral`/`keep`/`n_pore`/`n_lat`/`escaped_ranges`,
  plus a `lobes` key `_conn_frame_lobes` reads straight off), so every
  downstream caller — including two-tone coloring, which only calls
  classify, never cluster — is unchanged.
- `split` classifies against a region's `lines` list (written to a small
  temp label `.sph`), never the mesh-building `rsph` file, which for the
  "pore" region also carries centreline/escape points that must not enter
  the classification grid — mirrors the Tcl split function exactly, not
  just approximately.
- A region absent from `out_of` (the caller already has a valid cached
  surface for it) is passed as `--region NAME LABELS -` — registered in
  the classification grid so nearby triangles still resolve correctly, but
  written nowhere, matching Tcl's cached-region skip exactly.
- Every Tcl call site falls back to the original pure-Tcl implementation
  whenever the binary is missing or a run fails (`_conn_classify_native`,
  `_split_conn_mesh_native` return an empty dict / -1) — a stale or absent
  binary degrades speed, never correctness.
- Built and installed the same way as `nm_search` (`sh native/build.sh`,
  discovered next to `sos_triangle`, `state(conn_lobes_exec)` overrides).
