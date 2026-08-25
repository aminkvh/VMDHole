#!/bin/sh
# Build every binary VMDHole uses, in one step. This is the install script:
#   * the VMDHole source patches applied (fast CAPSULE ~40x, fast CONNOLLY area,
#     and PARALLEL CONNOLLY ~6.5x via OpenMP) - byte-identical .sph + all numbers
#     to stock HOLE. See connolly_patches/ (originals backed up to *.vmdhole_orig);
#   * all HOLE binaries compiled with -O2. HOLE's own Makefile sets NO Fortran
#     optimisation flag, so as distributed they are -O0. Measured on 1BL8, same
#     unpatched source, only FFLAGS differing, output byte-identical:
#     2.04x (circular) and 2.91x (CONNOLLY) - see
#     ../paper/benchmarks/results/o2_rebuild.csv;
#   * the stock `sos_triangle` replaced by the optimised `sos_triangle_fast`
#     (~9-10x faster, byte-identical output).
#
# The result is collected into one output folder (default: ./build) so you can
# point VMDHole's executable paths (File -> Settings...) at a single directory.
# Your existing HOLE install is NOT modified.
#
# Usage:
#   ./build-vmdhole-optimized.sh [HOLE2_SRC_DIR] [OUTPUT_DIR] [INSTALL_DIR]
#     HOLE2_SRC_DIR  path to a hole2/src checkout; if omitted, HOLE2 is cloned
#     OUTPUT_DIR     where to place the finished binaries (default: ./build)
#     INSTALL_DIR    also install into this exe dir (backs up originals); or set
#                    HOLE_INSTALL_DIR. e.g. HOLE_INSTALL_DIR="$HOME/hole2/exe"
#
# Env: CC / FC / OPT override the C compiler / Fortran compiler / opt level.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$1"
OUT="${2:-$HERE/build}"
CC="${CC:-cc}"
OPT="${OPT:--O2}"
HOLE_FFLAGS="$OPT -fd-lines-as-comments -fbackslash -std=legacy"
# Default -O2: portable binaries AND byte-identical sos_triangle output (the fork's
# correctness is verified at -O2). For ~7% more speed on THIS machine only you may
# set OPT="-O2 -march=native", but native binaries are NOT portable to other CPUs,
# and -march=native can enable FMA, which changes floating-point rounding and may
# break byte-identity versus the stock binary. Do not ship native-built binaries.
CLEAN_CLONE=""
# Pinned upstream revision. A moving HEAD meant "rebuild the same release six
# months later" could silently produce a DIFFERENT scientific binary from the
# same patch set. Override only deliberately, e.g. HOLE2_REF=master to test
# against current upstream.
HOLE2_REF="${HOLE2_REF:-a8eaf6121ba66625446933f4acd7d6aa336dbb47}"

# 0. Prerequisites ---------------------------------------------------------
# Checked up front and reported together, so a missing toolchain fails in one
# readable message instead of part-way through a clone or a Fortran build.
echo ">> Checking prerequisites ..."
_missing=""
for _t in make python3 git; do
  command -v "$_t" >/dev/null 2>&1 || _missing="$_missing $_t"
done
CC="${CC:-cc}"
FC="${FC:-gfortran}"
command -v "$CC" >/dev/null 2>&1 || _missing="$_missing $CC(C compiler)"
command -v "$FC" >/dev/null 2>&1 || _missing="$_missing $FC(Fortran compiler)"
if [ -n "$_missing" ]; then
  echo "ERROR: missing required tool(s):$_missing" >&2
  echo "       Needed: a POSIX shell, make, Python 3, Git, and C and Fortran" >&2
  echo "       compilers with OpenMP and legacy-Fortran support." >&2
  exit 1
fi
_probe="$(mktemp -d)"
# OpenMP, both compilers: the CONNOLLY/CAPSULE parallel regions need it in
# Fortran, sos_triangle_fast's slice loop in C. Missing OpenMP is not fatal -
# both fall back to serial - so this warns rather than stops.
printf 'program p\n!$ integer :: t\nend program p\n' > "$_probe/omp.f90"
if ! "$FC" -fopenmp -o "$_probe/ompf" "$_probe/omp.f90" >/dev/null 2>&1; then
  echo ">> WARNING: $FC has no working -fopenmp; CONNOLLY/CAPSULE will build serial." >&2
fi
printf 'int main(void){return 0;}\n' > "$_probe/omp.c"
if ! "$CC" -fopenmp -o "$_probe/ompc" "$_probe/omp.c" >/dev/null 2>&1; then
  echo ">> WARNING: $CC has no working -fopenmp; sos_triangle will build serial." >&2
fi
# Legacy Fortran: HOLE 2.2 is fixed-form FORTRAN 77 with lines past column 72.
# Without -std=legacy (or equivalent) gfortran rejects it outright, so this one
# IS fatal - the build cannot proceed.
printf '      PROGRAM P\n      END\n' > "$_probe/f77.f"
if ! "$FC" -std=legacy -c -o "$_probe/f77.o" "$_probe/f77.f" >/dev/null 2>&1; then
  echo "ERROR: $FC cannot compile fixed-form legacy FORTRAN (-std=legacy)." >&2
  echo "       HOLE 2.2 is FORTRAN 77; install gfortran or set FC." >&2
  rm -rf "$_probe"; exit 1
fi
rm -rf "$_probe"
echo ">> Prerequisites OK (make, python3, git, $CC, $FC)."

# 1. HOLE2 source ----------------------------------------------------------
if [ -z "$SRC" ]; then
  CLEAN_CLONE="$(mktemp -d)"
  echo ">> No HOLE2 source given; cloning into $CLEAN_CLONE ..."
  git clone https://github.com/osmart/hole2.git "$CLEAN_CLONE/hole2"
  ( cd "$CLEAN_CLONE/hole2" && git checkout --quiet "$HOLE2_REF" ) || {
    echo "ERROR: could not check out pinned HOLE2 revision $HOLE2_REF." >&2
    echo "       Set HOLE2_REF to override the pin." >&2
    exit 1; }
  echo ">> HOLE2 pinned at $HOLE2_REF"
  SRC="$CLEAN_CLONE/hole2/src"
fi
[ -f "$SRC/Makefile" ] || { echo "ERROR: '$SRC' is not a hole2/src checkout (no Makefile)."; exit 1; }

# 1b. Apply the VMDHole source patches (fast CAPSULE + parallel CONNOLLY) ---
# Idempotent; backs up every file it touches to <file>.vmdhole_orig. Adds
# hcapen_fast/coarea_fast/holcal_par/holeen_par + an RNG counter, and per-file
# -fopenmp Makefile rules (only the parallel-region files - global -fopenmp
# would seg-fault). Requires python3.
command -v python3 >/dev/null 2>&1 || {
  echo "ERROR: python3 is required to apply the CONNOLLY/CAPSULE source patches." >&2
  echo "       Install python3, or apply connolly_patches/ to '$SRC' by hand." >&2
  exit 1; }
echo ">> Applying VMDHole CONNOLLY/CAPSULE patches ..."
python3 "$HERE/connolly_patches/apply_patches.py" "$SRC"

# 2. Build all HOLE binaries with -O2 (per-file -fopenmp added by the patch) -
echo ">> Building HOLE binaries with FFLAGS='$HOLE_FFLAGS' ..."
( cd "$SRC" && make clean >/dev/null 2>&1 || true
  make all FFLAGS="$HOLE_FFLAGS" CFLAGS="$OPT" )
HOLE_EXE="$(cd "$SRC/../exe" && pwd)"

# 3. Build the fast sos_triangle ------------------------------------------
# -fopenmp threads the single-frame ellipse-probe slice loop (--asym-threads N); optional -
# if the compiler lacks OpenMP the pragmas are ignored and it builds/runs serially.
echo ">> Building sos_triangle_fast ($OPT) ..."
# Same translation units and libraries as build.sh, deliberately: this script and
# that one both produce a file called sos_triangle, and without distinct names they produce
# DIFFERENT binaries (this one omitted vor_predicates.c/vor_delaunay.c and
# -lpthread). Only one symbol of those units survives the link, so nothing was
# functionally missing - but "which sos_triangle is running?" is not a question
# a user following the docs should have to ask, and only build.sh's output is the
# one the test suite checks.
if $CC $OPT -fopenmp -o "$HERE/sos_triangle_fast" "$HERE/sos_triangle_fast.c" \
        "$HERE/voronoi/vor_predicates.c" "$HERE/voronoi/vor_delaunay.c" -lm -lpthread 2>/dev/null; then
  echo "   (OpenMP enabled)"
else
  $CC $OPT -o "$HERE/sos_triangle_fast" "$HERE/sos_triangle_fast.c" \
        "$HERE/voronoi/vor_predicates.c" "$HERE/voronoi/vor_delaunay.c" -lm -lpthread
  echo "   (serial, no OpenMP)"
fi

# 4. Collect everything into the output folder ----------------------------
echo ">> Collecting optimised binaries into $OUT ..."
mkdir -p "$OUT"
cp -f "$HOLE_EXE"/* "$OUT"/ 2>/dev/null || true   # all -O2 HOLE binaries
cp -f "$HOLE_EXE/sos_triangle" "$OUT/sos_triangle.stock" 2>/dev/null || true  # keep a backup
cp -f "$HERE/sos_triangle_fast" "$OUT/sos_triangle"       # fast one becomes THE sos_triangle

# 4b. Detection manifest --------------------------------------------------
# VMDHole reads this file (next to the hole binary) to show the per-feature
# green acceleration checks in File -> Settings and to pick Way B (OpenMP)
# thread counts. One "patch=" line per patched source file.
OMP_OK=no
if ldd "$OUT/hole" 2>/dev/null | grep -qi "libgomp"; then OMP_OK=yes; fi
{
  echo "vmdhole_accel_manifest 1"
  echo "built $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  echo "opt $OPT"
  echo "openmp $OMP_OK"
  echo "patch hcapen_fast.f fast-capsule"
  echo "patch coarea_fast.f fast-connolly-area"
  echo "patch holcal_par.f parallel-connolly"
  echo "patch holeen_par.f threadprivate-cache"
  # sph_process's parallel dot cull. apply_patches.py has ALWAYS applied
  # sphqpu_par.f (it is in its patch list and in OMP_FILES), but this manifest
  # never said so - and update_sph_accel_status looks for exactly this line. So
  # a genuinely accelerated sph_process was reported "not accelerated", with
  # its toggle DISABLED, i.e. the user could not even force the serial cull.
  echo "patch sphqpu_par.f parallel-dot-cull"
  echo "patch machine_dep.g77 rng-counter"
  echo "patch Makefile openmp-build"
  # Lets the plugin skip the per-frame writepdb/re-parse round trip: 12.1 ms/frame
  # of ASCII formatting at 18.7k atoms, and it is SERIAL while HOLE runs across
  # every worker, so it is most of the parallel ceiling for the fast pore methods.
  # Lets the plugin hand the tunnel engine a packed coordinate record instead of
  # the 12-column text: the text build costs ~41 ms/frame of Tcl formatting at
  # 18.7k atoms, and that is the rate the whole tunnel search is fed at.
  echo "patch mole_atoms_bin fast-atoms-read"
  echo "patch tsatr_fast.f fast-coord-read"
  echo "patch sos_triangle fast-surface"
  # Provenance: which upstream tree these patches were applied to, and what
  # came out. Without these a manifest recorded only patch NAMES, so two
  # binaries built from different upstream revisions looked identical here.
  echo "hole2_ref $HOLE2_REF"
  if [ -d "$SRC/../.git" ]; then
    echo "hole2_commit $(cd "$SRC/.." && git rev-parse HEAD 2>/dev/null || echo unknown)"
  else
    echo "hole2_commit unknown-not-a-git-checkout"
  fi
  for _b in hole sos_triangle sph_process; do
    [ -f "$OUT/$_b" ] && echo "sha256 $_b $(sha256sum "$OUT/$_b" | cut -d" " -f1)"
  done
} > "$OUT/vmdhole_accel.manifest"
echo ">> Wrote detection manifest: $OUT/vmdhole_accel.manifest (openmp=$OMP_OK)"

# 5. Smoke test the fast sos on a committed fixture ------------------------
if [ -f "$HERE/benchmarks/fixtures/d15.sos" ]; then
  if "$OUT/sos_triangle" -s < "$HERE/benchmarks/fixtures/d15.sos" 2>/dev/null | grep -q "draw "; then
    echo ">> Smoke test: OK."
  else
    echo ">> Smoke test: WARNING - sos_triangle produced no surface." >&2
  fi
fi

# 5b. Build the MOLE 2 tunnel engine ---------------------------------------
# Its own executable, not a mode of sos_triangle: it needs the exact predicates
# at VP_SCALE 1e5 so MOLE's +-0.00005 A general-position jitter survives, while
# sos_triangle_fast and --tunnel-voronoi need the default 1e3 and have
# byte-identity tests pinned to it. One binary cannot hold both.
echo ">> Building mole_tunnel_engine ($OPT) ..."
MOLE_SRC="mole/mole_main.c mole/mole_tunnel.c mole/mole_lining.c mole/mole_complex.c mole/mole_dh.c mole/mole_rng.c voronoi/vor_delaunay.c voronoi/vor_predicates.c"
if ( cd "$HERE" && $CC $OPT -DVP_SCALE=100000.0 -DVP_MAX_COORD=20000000L \
       -o "$OUT/mole_tunnel_engine" $MOLE_SRC -lm ); then
  echo ">> Built mole_tunnel_engine (tunnel mode)."
else
  echo ">> WARNING: mole_tunnel_engine did not build - tunnel mode will have no engine." >&2
fi

# 6. Optionally INSTALL straight into the HOLE exe dir you actually use ----
# So you never have to copy by hand again: set HOLE_INSTALL_DIR (or pass a 3rd
# arg) to the exe dir VMDHole points at, e.g.
#   HOLE_INSTALL_DIR="$HOME/hole2/exe" ./build-vmdhole-optimized.sh
# The current hole/sph_process/sos_triangle there are backed up (timestamped)
# first, then overwritten, and the detection manifest is copied alongside so the
# green acceleration checks light up (VMDHole reads it next to the hole binary).
INSTALL_DIR="${3:-${HOLE_INSTALL_DIR:-}}"
if [ -n "$INSTALL_DIR" ]; then
  if [ ! -d "$INSTALL_DIR" ]; then
    echo ">> HOLE_INSTALL_DIR '$INSTALL_DIR' does not exist - skipping install." >&2
  else
    STAMP="$(date +%Y%m%d-%H%M%S)"
    echo ">> Installing into $INSTALL_DIR (originals backed up as *.bak-$STAMP) ..."
    for b in hole sph_process sos_triangle mole_tunnel_engine; do
      [ -f "$OUT/$b" ] || continue
      [ -f "$INSTALL_DIR/$b" ] && cp -p "$INSTALL_DIR/$b" "$INSTALL_DIR/$b.bak-$STAMP"
      cp -f "$OUT/$b" "$INSTALL_DIR/$b"
    done
    cp -f "$OUT/vmdhole_accel.manifest" "$INSTALL_DIR/vmdhole_accel.manifest"
    echo ">> Installed. Restart VMDHole (or press Detect in Settings) to pick it up."
  fi
fi

[ -n "$CLEAN_CLONE" ] && rm -rf "$CLEAN_CLONE"

cat <<EOF

Done. Every binary VMDHole uses is in:
    $OUT
        hole, sph_process, qpt_conv, ...   (rebuilt with $OPT)
        sos_triangle                        (= the fast build)
        sos_triangle.stock                  (unmodified sos_triangle, kept for comparison)
        mole_tunnel_engine                  (tunnel mode)

Point VMDHole at it in File -> Settings... (set hole / sph_process / sos_triangle
to the files above), or copy them over your existing hole2/exe binaries. VMDHole
picks the thread count and stack size automatically (Way A/B), and shows a green
acceleration check per patched file in Settings.

To install directly into the exe dir you use (no manual copy), re-run with:
    HOLE_INSTALL_DIR="\$HOME/hole2/exe" ./build-vmdhole-optimized.sh
(backs up the old binaries, copies in the new ones + the detection manifest).

Running the accelerated 'hole' by hand for a parallel CONNOLLY frame:
    OMP_NUM_THREADS=8 OMP_STACKSIZE=256M hole < hole.inp
OMP_STACKSIZE must NOT be set below ~4M (worker threads hold a ~2.8MB array);
unset is fine. Verify identical output any time with: ACCEL=$OUT/hole ./verify.sh
EOF
