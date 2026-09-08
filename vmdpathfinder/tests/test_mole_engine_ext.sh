#!/bin/sh
# The three additions to the C MOLE engine (native/mole/mole_main.c) that the
# plugin reads or drives:
#
#   VP records   - a cavity's geometry, one "VP id x y z r" per member
#                  tetrahedron, right after its V/VB/VI block. Centre and
#                  radius are a tunnel P record's own quantities, so a cavity
#                  drawn as a sphere union sits on the tunnel drawn the same way.
#   --origin=    - repeatable pinned origins; each is searched like the
#                  positional ox oy oz and the tunnels are merged as the
#                  auto-origin loop already merges them.
#   --exit=      - repeatable; every exit joins CustomExits.
#   --vdw=El:r   - per-element vdW radius overrides.
#
# test_mole_tcl_port pins the engine to MOLE's own numbers; nothing there
# exercises these flags, and the byte-parity checks there only see that the
# OLD lines are unchanged. This is the positive side: that the new lines and
# flags do what they say. Everything runs on the 1tqn fixture, whose two
# pinned origins (MOLE's test.xml origin and the heme iron) and whose tunnel
# mouths are already characterised by the port's own reference runs.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SRC="$DIR/../../native"
FIX="$DIR/fixtures/mole_atoms_1tqn.txt"
echo "=============================================================="
echo "mole-engine-ext: $SRC/mole/mole_main.c on $(basename "$FIX")"
pass=0; fail=0
ok()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad() { echo "  FAIL  $1"; fail=$((fail+1)); }

CC=${CC:-cc}
command -v "$CC" >/dev/null 2>&1 || { echo "SKIP: no C compiler"; exit 0; }
[ -f "$FIX" ] || { echo "  FAIL  fixture missing: $FIX"; exit 1; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Built from current source, with the engine's own predicate grid - the same
# translation units and flags as native/build.sh.
MOLE_FLAGS="-DVP_SCALE=100000.0 -DVP_MAX_COORD=20000000L"
E="$TMP/engine"
if ! "$CC" -O2 $MOLE_FLAGS -o "$E" "$SRC/mole/mole_main.c" \
        "$SRC/mole/mole_tunnel.c" "$SRC/mole/mole_lining.c" "$SRC/mole/mole_complex.c" \
        "$SRC/mole/mole_dh.c" "$SRC/mole/mole_rng.c" "$SRC/voronoi/vor_delaunay.c" \
        "$SRC/voronoi/vor_predicates.c" -lm 2>"$TMP/err"; then
    bad "mole_tunnel_engine did not build from current source"
    sed 's/^/       /' "$TMP/err"
    echo "  -> $pass passed, $fail failed"
    exit 1
fi

# MOLE's test.xml origin (centroid of A308/A309) and the heme iron - the two
# pinned starts the reference runs use. The first snaps into cavity 1 only,
# the second into cavities 1 and 7.
ORI=$(awk '$6=="A" && ($7==308 || $7==309) {x+=$1;y+=$2;z+=$3;n++} \
           END{printf "%.6f %.6f %.6f", x/n, y/n, z/n}' "$FIX")
ORIC=$(echo "$ORI" | tr ' ' ',')
FE=$(awk '$8=="HEM" && $10=="FE" {print $1, $2, $3; exit}' "$FIX")
FEC=$(echo "$FE" | tr ' ' ',')
DEF="3.0 1.25 8 5.0 0.0 0"
# "cav length" per tunnel, sorted: the identity of a tunnel across runs whose
# ids are renumbered by the (Cavity.Id, Length) sort.
tunset() { awk '$1=="T"{printf "%d %.4f\n", $7, $4}' "$1" | sort; }

# ---- 1. VP: cavity geometry ---------------------------------------------------
MOLE_CAVITY_DEBUG=1 "$E" "$FIX" "$TMP/auto.txt" $DEF > "$TMP/auto.err" 2>&1 \
    || bad "engine failed on the auto-origin run: $(tail -1 "$TMP/auto.err")"
if grep -qxF "# VP id x y z r  (cavity tetrahedra as centre+clearance spheres)" "$TMP/auto.txt"; then
    ok "header documents the VP record"
else
    bad "no '# VP ...' header line"
fi
# Appended after the existing '# V ...' line, not spliced into it.
if awk 'prev ~ /^# V id/ && $0 ~ /^# VP id/ {f=1} {prev=$0} END{exit !f}' "$TMP/auto.txt"; then
    ok "the VP header line follows the V header line"
else
    bad "the VP header line is not directly after the V header line"
fi
nv=$(grep -c '^V ' "$TMP/auto.txt")
nvp=$(grep -c '^VP ' "$TMP/auto.txt")
[ "$nv" -gt 0 ] && [ "$nvp" -gt 0 ] && ok "$nv cavities carry $nvp VP lines" \
    || bad "expected V and VP records (V=$nv VP=$nvp)"
# Per cavity: exactly as many VP lines as the cavity has tetrahedra (the
# engine's own count, printed under MOLE_CAVITY_DEBUG), and every one has a
# positive radius.
awk '$1=="VP"{n[$2]++; if (!($6>0)) badr++} END{for (i in n) print i, n[i]; if (badr) print "BADR", badr}' \
    "$TMP/auto.txt" | sort -n > "$TMP/vp_counts.txt"
awk '$1=="CAVITY"{print $2, $4}' "$TMP/auto.err" | sort -n > "$TMP/cav_counts.txt"
if [ -s "$TMP/cav_counts.txt" ] && cmp -s "$TMP/vp_counts.txt" "$TMP/cav_counts.txt"; then
    ok "every cavity's VP count equals its tetrahedron count ($(tr '\n' ' ' < "$TMP/cav_counts.txt"| sed 's/ $//'))"
else
    bad "VP counts differ from the cavities' tetrahedron counts"
    paste "$TMP/vp_counts.txt" "$TMP/cav_counts.txt" | sed 's/^/       /'
fi
if ! grep -q BADR "$TMP/vp_counts.txt"; then ok "every VP radius is positive"
else bad "$(grep BADR "$TMP/vp_counts.txt" | awk '{print $2}') VP lines with r <= 0"; fi
# Placement: V, VB, VI, then that cavity's VP lines and nothing else until the
# next cavity or the first tunnel. A VP anywhere else, or under another id,
# is a parser trap.
if awk 'BEGIN{state=0; bad=0}
    $1=="V"  { if (state!=0 && state!=3) bad=1; state=1; id=$2; next }
    $1=="VB" { if (state!=1) bad=1; state=2; next }
    $1=="VI" { if (state!=2) bad=1; state=3; next }
    $1=="VP" { if (state!=3 || $2!=id) bad=1; next }
    $1=="T"  { if (state!=3 && state!=0 && state!=9) bad=1; state=9 }
    $1=="P"||$1=="L"||$1=="F"||$1=="H"||$1=="Y"||$1=="W" { if (state!=9) bad=1 }
    END { exit bad }' "$TMP/auto.txt"; then
    ok "VP lines sit immediately after each cavity's V/VB/VI block"
else
    bad "VP lines are out of place"
fi
# The quantity is a tunnel's: the profile spline passes through its first
# control tetrahedron's circumcentre at t=0 and samples the clamped clearance
# there, so a tunnel's first P point must be one of its cavity's VP spheres,
# digit for digit. Checked on the pinned-origin run, whose tunnels come from
# a single, known cavity.
"$E" "$FIX" "$TMP/ori.txt" $DEF $ORI 5.0 > /dev/null 2>&1 \
    || bad "engine failed on the pinned-origin run"
nt=$(grep -c '^T ' "$TMP/ori.txt")
miss=$(awk '$1=="VP"{vp[$2" "$3" "$4" "$5" "$6]=1}
            $1=="T"{cav[$2]=$7; first[$2]=1}
            $1=="P" && first[$2]{first[$2]=0; if (cav[$2]>0 && !(cav[$2]" "$3" "$4" "$5" "$6 in vp)) m++}
            END{print m+0}' "$TMP/ori.txt")
if [ "$nt" -gt 0 ] && [ "$miss" -eq 0 ]; then
    ok "each tunnel's first P point is one of its cavity's VP spheres ($nt tunnels)"
else
    bad "$miss of $nt tunnels start at a point that is not a VP sphere of their cavity"
fi

# ---- 2. --vdw ----------------------------------------------------------------
# A no-op override (the table's own value, in either case) must leave the file
# untouched; a real one must move it. Malformed entries must refuse to run
# rather than search on the built-in radius.
for v in "C:1.61" "c:1.61" "fe:1.7,O:1.45,n:1.55"; do
    "$E" "$FIX" "$TMP/vdw_noop.txt" $DEF --vdw=$v > /dev/null 2>&1
    if cmp -s "$TMP/auto.txt" "$TMP/vdw_noop.txt"; then ok "--vdw=$v (table values) is byte-identical"
    else bad "--vdw=$v changed the output although it equals the table"; fi
done
"$E" "$FIX" "$TMP/vdw_real.txt" $DEF --vdw=C:1.9 > /dev/null 2>&1
if [ -s "$TMP/vdw_real.txt" ] && ! cmp -s "$TMP/auto.txt" "$TMP/vdw_real.txt"; then
    ok "--vdw=C:1.9 changes the result ($(grep -c '^T ' "$TMP/vdw_real.txt") tunnels vs $(grep -c '^T ' "$TMP/auto.txt"))"
else
    bad "--vdw=C:1.9 left the output unchanged"
fi
for v in "XX:1.5" "C:abc" "C" "" "C:1.7," "C:-1" "ZN:1.39,N"; do
    rm -f "$TMP/bad.txt"
    "$E" "$FIX" "$TMP/bad.txt" $DEF "--vdw=$v" > "$TMP/bad.err" 2>&1; rc=$?
    if [ "$rc" -eq 2 ] && [ ! -f "$TMP/bad.txt" ] && grep -q -- "--vdw" "$TMP/bad.err"; then
        ok "--vdw=$v refused (exit 2, no output written)"
    else
        bad "--vdw=$v: exit $rc, expected 2 with a --vdw diagnostic"
    fi
done

# ---- 3. --origin ---------------------------------------------------------------
"$E" "$FIX" "$TMP/o_flag.txt" $DEF --origin=$ORIC > /dev/null 2>&1
if cmp -s "$TMP/ori.txt" "$TMP/o_flag.txt"; then ok "--origin=x,y,z equals the positional origin byte for byte"
else bad "--origin=x,y,z differs from the positional origin"; fi
"$E" "$FIX" "$TMP/o_fe.txt" $DEF $FE 5.0 > /dev/null 2>&1
"$E" "$FIX" "$TMP/o_two.txt" $DEF $ORI 5.0 --origin=$FEC > /dev/null 2>&1
tunset "$TMP/ori.txt" > "$TMP/s1.txt"; tunset "$TMP/o_fe.txt" > "$TMP/s2.txt"
tunset "$TMP/o_two.txt" > "$TMP/s12.txt"
n1=$(wc -l < "$TMP/s1.txt"); n2=$(wc -l < "$TMP/s2.txt"); n12=$(wc -l < "$TMP/s12.txt")
if [ "$n12" -ge "$n1" ] && [ -z "$(comm -23 "$TMP/s1.txt" "$TMP/s12.txt")" ]; then
    ok "a second origin yields a superset ($n1 -> $n12 tunnels)"
else
    bad "two origins lost a tunnel of the first ($n1 -> $n12)"
fi
# Both origins snap into cavity 1, so it gets two origins and each is searched
# on its own, as GetTunnels does per origin: the result is exactly the union.
if sort "$TMP/s1.txt" "$TMP/s2.txt" | cmp -s - "$TMP/s12.txt"; then
    ok "two origins give exactly the union of the two single-origin runs ($n1 + $n2 = $n12)"
else
    bad "two-origin run is not the union of the single-origin runs"
    sort "$TMP/s1.txt" "$TMP/s2.txt" > "$TMP/s1s2.txt"
    diff "$TMP/s12.txt" "$TMP/s1s2.txt" | sed 's/^/       /'
fi
"$E" "$FIX" "$TMP/o_ff.txt" $DEF --origin=$ORIC --origin=$FEC > /dev/null 2>&1
if cmp -s "$TMP/o_two.txt" "$TMP/o_ff.txt"; then ok "positional + --origin equals --origin + --origin"
else bad "positional + --origin differs from two --origin flags"; fi
"$E" "$FIX" "$TMP/o_dup.txt" $DEF $ORI 5.0 --origin=$ORIC > /dev/null 2>&1
if cmp -s "$TMP/ori.txt" "$TMP/o_dup.txt"; then ok "the same point pinned twice is one origin"
else bad "a duplicated origin changed the output"; fi
args=""; i=0
while [ "$i" -lt 16 ]; do args="$args --origin=$i,0,0"; i=$((i+1)); done
rm -f "$TMP/bad.txt"
"$E" "$FIX" "$TMP/bad.txt" $DEF $ORI 5.0 $args > "$TMP/bad.err" 2>&1; rc=$?
if [ "$rc" -eq 2 ] && [ ! -f "$TMP/bad.txt" ]; then ok "a 17th pinned origin is refused (exit 2)"
else bad "17 pinned origins: exit $rc, expected 2"; fi
rm -f "$TMP/bad.txt"
"$E" "$FIX" "$TMP/bad.txt" $DEF --origin=1,2 > /dev/null 2>&1; rc=$?
[ "$rc" -eq 2 ] && ok "--origin=1,2 is refused (exit 2)" || bad "--origin=1,2: exit $rc, expected 2"

# ---- 4. --exit, repeated -------------------------------------------------------
# Two mouths of cavity 1's own tunnels (from the auto-origin run), each of
# which alone yields one tunnel under UseCustomExitsOnly. Together they must
# yield both, each ending at its own exit.
E1="-22.5173,-9.0206,-2.4857"; E2="-29.2520,-31.3578,-13.3185"
"$E" "$FIX" "$TMP/x1.txt" $DEF $ORI 5.0 --exit=$E1 --exitsonly=1 > /dev/null 2>&1
"$E" "$FIX" "$TMP/x2.txt" $DEF $ORI 5.0 --exit=$E2 --exitsonly=1 > /dev/null 2>&1
"$E" "$FIX" "$TMP/x12.txt" $DEF $ORI 5.0 --exit=$E1 --exit=$E2 --exitsonly=1 > /dev/null 2>&1
c1=$(grep -c '^T ' "$TMP/x1.txt"); c2=$(grep -c '^T ' "$TMP/x2.txt"); c12=$(grep -c '^T ' "$TMP/x12.txt")
if [ "$c1" -eq 1 ] && [ "$c2" -eq 1 ] && [ "$c12" -eq 2 ]; then
    ok "two --exit flags yield both tunnels ($c1 + $c2 = $c12)"
else
    bad "two --exit flags: $c1 + $c2 tunnels alone, $c12 together"
fi
tunset "$TMP/x1.txt" > "$TMP/sx1.txt"; tunset "$TMP/x2.txt" > "$TMP/sx2.txt"
tunset "$TMP/x12.txt" > "$TMP/sx12.txt"
if sort "$TMP/sx1.txt" "$TMP/sx2.txt" | cmp -s - "$TMP/sx12.txt"; then
    ok "the two-exit run is the union of the single-exit runs"
else
    bad "the two-exit run is not the union of the single-exit runs"
fi
# Each tunnel ends within 2 A of one of the exits, and each exit is reached.
hit=$(awk -v e1="$E1" -v e2="$E2" 'BEGIN{split(e1,a,","); split(e2,b,",")}
    $1=="P"{lx[$2]=$3; ly[$2]=$4; lz[$2]=$5}
    END{n=0; for (t in lx) {
          d1=sqrt((lx[t]-a[1])^2+(ly[t]-a[2])^2+(lz[t]-a[3])^2)
          d2=sqrt((lx[t]-b[1])^2+(ly[t]-b[2])^2+(lz[t]-b[3])^2)
          if (d1<2) h1=1; else if (d2<2) h2=1; else far++ }
        print (far+0), (h1+0), (h2+0)}' "$TMP/x12.txt")
set -- $hit
if [ "$1" -eq 0 ] && [ "$2" -eq 1 ] && [ "$3" -eq 1 ]; then
    ok "every tunnel ends at one of the two exits, and both exits are reached"
else
    bad "tunnel ends vs exits: $1 far from both, exit1 reached=$2, exit2 reached=$3"
fi
# A surface-only exit alongside a regular one: the SurfaceCavity (C0) route
# and the cavity-1 route both survive, from the same origin.
"$E" "$FIX" "$TMP/xs.txt" $DEF $ORI 5.0 --exit=-22.841,-55.065,-20.445 --exit=$E1 --exitsonly=1 > /dev/null 2>&1
cavs=$(awk '$1=="T"{print $7}' "$TMP/xs.txt" | sort -n | uniq | tr '\n' ' ')
if [ "$cavs" = "0 1 " ]; then ok "surface exit + cavity exit: tunnels from C0 and C1 (cavities: $cavs)"
else bad "surface exit + cavity exit: cavities [$cavs], expected 0 and 1"; fi

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
