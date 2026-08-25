#!/bin/sh
# Build the self-contained C tools from source. This is the entry point CI
# (.github/workflows/tests.yml) runs on a fresh checkout, and the one the
# plugin's own "No MOLE tunnel engine found" message points users at - so it
# must need nothing beyond a C compiler. The HOLE-side Fortran acceleration
# (patched hole/sph_process) is a separate, optional step that needs a HOLE
# source tree: see build-vmdhole-optimized.sh.
#
# Produces, next to this script:
#   sos_triangle_fast   - surface triangulation, properties, tunnel clustering
#   mole_tunnel_engine  - the MOLE 2 tunnel-search engine
#   hydro_project       - hydration water-projection accelerator
set -e
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
CC="${CC:-cc}"
OPT="${OPT:--O2}"

echo ">> Building sos_triangle_fast ($OPT) ..."
# OpenMP is optional: without it the pragmas are ignored and the binary runs
# serially (Apple clang has no -fopenmp, and CI builds on macOS too).
if $CC $OPT -fopenmp -o "$HERE/sos_triangle_fast" "$HERE/sos_triangle_fast.c" \
        "$HERE/voronoi/vor_predicates.c" "$HERE/voronoi/vor_delaunay.c" -lm -lpthread 2>/dev/null; then
  echo "   (OpenMP enabled)"
else
  $CC $OPT -o "$HERE/sos_triangle_fast" "$HERE/sos_triangle_fast.c" \
        "$HERE/voronoi/vor_predicates.c" "$HERE/voronoi/vor_delaunay.c" -lm -lpthread
  echo "   (serial, no OpenMP)"
fi

echo ">> Building mole_tunnel_engine ($OPT) ..."
# Its own executable, not a mode of sos_triangle: it needs the exact predicates
# at VP_SCALE 1e5 so MOLE's +-0.00005 A general-position jitter survives, while
# sos_triangle_fast needs the default 1e3 and has byte-identity tests pinned to
# it. One binary cannot hold both. Same translation units and flags as
# build-vmdhole-optimized.sh, deliberately - two scripts producing DIFFERENT
# binaries under one name is the trap that script's own comments warn about.
( cd "$HERE" && $CC $OPT -DVP_SCALE=100000.0 -DVP_MAX_COORD=20000000L \
    -o "$HERE/mole_tunnel_engine" \
    mole/mole_main.c mole/mole_tunnel.c mole/mole_lining.c mole/mole_complex.c \
    mole/mole_dh.c mole/mole_rng.c voronoi/vor_delaunay.c voronoi/vor_predicates.c -lm )

echo ">> Building hydro_project ..."
# -ffp-contract=off is NOT optional: without it the compiler may fuse
# multiply-add pairs into single FMA instructions, which round once instead of
# twice like Tcl's expr-by-expr evaluation and silently break bit-identity with
# the Tcl reference at the last ULP. No -Ofast/-ffast-math/-march=native for
# the same reason (see hydration/hydro_project.c's own header).
$CC $OPT -ffp-contract=off -o "$HERE/hydro_project" "$HERE/hydration/hydro_project.c" -lm
echo "   built: $HERE/hydro_project"

echo ">> Done."
