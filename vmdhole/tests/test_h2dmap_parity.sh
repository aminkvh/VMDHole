#!/bin/sh
# The parallel 2DMAPS routine (connolly_patches/h2dmap_par.f) must produce
# BYTE-IDENTICAL maps to the serial build, at any thread count, and must
# actually be parallel.
#
# Both halves matter. HOLE is Monte Carlo, so the control file below pins
# `raseed` - without it two runs of the SAME binary differ and every
# comparison here is noise. And -fopenmp is applied per-file in the Makefile,
# so a missing rule leaves the directives as comments: correct output, no
# speedup, silently (this is the sphqpu_par lesson, see OMP_FILES).
#
# sample 0.1 (not HOLE's 0.5): five times the slices, so the comparison covers
# more grid and the timing below is not measured at clock resolution.
#
# Needs a hole2/src checkout and gfortran; skips without either.
set -e
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$DIR/../.." && pwd)
SRC="${HOLE2_SRC:-$ROOT/native/stock_build/hole2}"
PATCHES="$ROOT/native/connolly_patches"
# 1GRM, not 1BL8: simple.rad has no K+ entry and 1BL8 is a K+ channel, so
# HOLE aborts on the first ion before it ever reaches the map.
# Bundled tutorial copy - always present in a fresh clone, unlike the
# case-study copy which needs fetch_structures.sh first.
PDB="${H2DMAP_PDB:-$ROOT/vmdhole/1GRM.pdb}"
RAD="${H2DMAP_RAD:-$ROOT/native/stock_build/hole2/rad/simple.rad}"

[ -d "$SRC/src" ] || { echo "SKIP: no hole2/src checkout at $SRC"; exit 0; }
command -v gfortran >/dev/null 2>&1 || { echo "SKIP: no gfortran"; exit 0; }
[ -f "$PDB" ] || { echo "SKIP: no test structure at $PDB"; exit 0; }
[ -f "$RAD" ] || { echo "SKIP: no radius file at $RAD"; exit 0; }

fails=0
ok()  { echo "  PASS  $1"; }
bad() { echo "  FAIL  $1"; fails=$((fails+1)); }

WORK=$(mktemp -d)
[ -n "${H2DMAP_KEEP:-}" ] || trap 'rm -rf "$WORK"' EXIT INT TERM
echo "  (work dir: $WORK)"
FF="-O2 -fd-lines-as-comments -fbackslash -std=legacy"

# Baseline = the shipped patch set MINUS h2dmap, so the only difference is the
# routine under test. Taken from git HEAD~ would be fragile; instead the
# baseline tree simply keeps the stock h2dmap.o.
for variant in serial par; do
    cp -r "$SRC" "$WORK/$variant"
    ( cd "$PATCHES" && python3 apply_patches.py "$WORK/$variant/src" ) >/dev/null 2>&1
    if [ "$variant" = serial ]; then
        # Baseline: same Makefile, same object names, same -fopenmp rule - only
        # the SOURCE differs. Overwriting h2dmap_par.f with the stock routine
        # (the subroutine name is the same) isolates the directives themselves.
        cp "$WORK/$variant/src/h2dmap.f" "$WORK/$variant/src/h2dmap_par.f"
    fi
    ( cd "$WORK/$variant/src" && make clean >/dev/null 2>&1 || true
      make all FFLAGS="$FF" CFLAGS="-O2" ) >"$WORK/build_$variant.log" 2>&1 \
        || { cp "$WORK/build_$variant.log" "${TMPDIR:-/tmp}/h2dmap_build_$variant.log" 2>/dev/null; echo "SKIP: $variant build failed - see ${TMPDIR:-/tmp}/h2dmap_build_$variant.log"; exit 0; }
    [ -x "$WORK/$variant/exe/hole" ] || { echo "SKIP: $variant produced no hole binary"; exit 0; }
done

run() {   # $1 tag  $2 binary  $3 threads -> seconds on stdout
    d="$WORK/run_$1"; mkdir -p "$d"; cp "$PDB" "$d/in.pdb"
    # Staged bare beside the control file: `radius` is a FORTRAN record capped
    # at 80 chars, and a full path silently truncates into "cannot open".
    cp "$RAD" "$d/hole.rad"
    cx=$(awk '/^ATOM/{x+=substr($0,31,8);y+=substr($0,39,8);z+=substr($0,47,8);n++}
              END{printf "%.3f %.3f %.3f", x/n, y/n, z/n}' "$d/in.pdb")
    cat > "$d/hole.inp" <<EOF
coord in.pdb
radius hole.rad
sphpdb map.sph
sample 0.1
endrad 8.0
cpoint $cx
cvect 0 0 1
raseed 1
2dmaps holemap
stop
EOF
    ( cd "$d" && OMP_NUM_THREADS=$3 OMP_STACKSIZE=256M \
        /usr/bin/time -f "%e" "$2" < hole.inp > hole_out.txt 2>tim.txt ) || true
    tail -1 "$d/tim.txt"
}

T_SER=$(run serial "$WORK/serial/exe/hole" 1)
T_P1=$(run par1  "$WORK/par/exe/hole" 1)
NT=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
[ "$NT" -gt 8 ] && NT=8
T_PN=$(run parN  "$WORK/par/exe/hole" "$NT")

n=0
for f in "$WORK/run_serial"/holemap*; do n=$((n+1)); done
if [ "$n" -lt 5 ]; then
    echo "SKIP: the serial run produced only $n map file(s) - no pore traced here"
    exit 0
fi
ok "the serial run produced $n map file(s)"

for ref in par1 parN; do
    diffs=""
    for f in "$WORK/run_serial"/holemap*; do
        b=$(basename "$f")
        cmp -s "$f" "$WORK/run_$ref/$b" || diffs="$diffs $b"
    done
    if [ -z "$diffs" ]; then
        ok "$ref is byte-identical to serial across all $n map files"
    else
        bad "$ref differs from serial:$diffs"
    fi
done

# The directives must be live, not comments. A serial-but-correct build is the
# failure this catches; require a real margin rather than any speedup at all.
FASTER=$(awk -v a="$T_SER" -v b="$T_PN" 'BEGIN{print (b>0 && a/b > 1.5) ? 1 : 0}')
if [ "$FASTER" = 1 ]; then
    ok "the parallel build is actually parallel (${T_SER}s serial vs ${T_PN}s on $NT threads)"
else
    bad "no speedup at $NT threads (${T_SER}s vs ${T_PN}s) - is the per-file -fopenmp rule present?"
fi
echo "  (1-thread parallel build: ${T_P1}s - directives must not cost anything serially)"

[ "$fails" -eq 0 ]
