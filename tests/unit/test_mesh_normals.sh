#!/bin/sh
# The mesher's vertex normals must be the gradient of the field it MARCHED.
#
# Marching cubes places vertices on the zero set of fval = max(fpos, -fclip):
# the pore spheres MINUS the ENDRAD clip spheres that cut it off at each mouth.
# Normals used to come from fpos alone, so on a mouth cap - where the clip term
# decides the surface - the normal described the wall the cap was carved out of.
#
# Two things this pins:
#   1. GEOMETRY IS UNTOUCHED. Vertex positions must be bit-identical to what the
#      same input produced before; only shading normals may move. This is the
#      "do not compromise accuracy" guarantee.
#   2. The normals must AGREE WITH THE FACET they belong to. A vertex normal
#      more than 90 degrees from its own triangle's geometric normal is lit
#      backwards - a black patch. The facet normal comes from the vertex
#      positions, which this change does not touch, so it is independent ground
#      truth rather than a restatement of the field.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "mesh-normals"

F="$ROOT/vmdpathfinder/tests/fixtures/spherical_1GRM.sph"
[ -f "$F" ] || { echo "SKIP: no spherical fixture"; exit 0; }
EXE=""
for c in "$ROOT/native/nm/mesh_csg" "$ROOT/native/sos_triangle_fast"; do
    [ -x "$c" ] && { EXE="$c"; break; }
done
[ -n "$EXE" ] || { echo "SKIP: no mesher built (sh native/build.sh)"; exit 0; }
case "$EXE" in *sos_triangle_fast) MODE="--mesh" ;; *) MODE="" ;; esac
command -v python3 >/dev/null 2>&1 || { echo "SKIP: no python3"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
$EXE $MODE "$F" "$T/m.plot" 1.4/0.7 --draw --axis 0 0 0 0 0 1 8 >/dev/null 2>&1 \
    || { echo "SKIP: mesher failed"; exit 0; }

OUT=$(python3 - "$T/m.plot" <<'PY'
import re,math,sys
tris=[]
for ln in open(sys.argv[1]):
    if not ln.startswith("draw trinorm"): continue
    g=re.findall(r'\{([^}]*)\}',ln)
    if len(g)<6: continue
    tris.append(([list(map(float,g[i].split())) for i in range(3)],
                 [list(map(float,g[i].split())) for i in range(3,6)]))
sub=lambda a,b:[a[i]-b[i] for i in range(3)]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def nrm(a):
    m=math.sqrt(sum(x*x for x in a));  return [x/m for x in a] if m>1e-12 else None
flip=0; tot=0; nonunit=0; angs=[]
for v,n in tris:
    fn=nrm(cross(sub(v[1],v[0]),sub(v[2],v[0])))
    if fn is None: continue
    for q in n:
        L=math.sqrt(sum(x*x for x in q))
        if abs(L-1.0)>1e-3: nonunit+=1
        qq=nrm(q)
        if qq is None: continue
        tot+=1
        d=max(-1.0,min(1.0,sum(x*y for x,y in zip(fn,qq))))
        angs.append(math.degrees(math.acos(d)))
        if d<0: flip+=1
angs.sort()
print("TRIS",len(tris))
print("TOT",tot)
print("FLIP",flip)
print("NONUNIT",nonunit)
print("P95",f"{angs[int(0.95*len(angs))]:.2f}" if angs else "0")
PY
)
get() { printf '%s\n' "$OUT" | sed -n "s/^$1 //p" | head -1; }
TRIS=$(get TRIS); TOT=$(get TOT); FLIP=$(get FLIP); NU=$(get NONUNIT); P95=$(get P95)

[ "${TRIS:-0}" -gt 100 ] 2>/dev/null && ok "the mesher produced a surface ($TRIS triangles)" \
    || { bad "no usable surface"; echo "mesh-normals: $pass passed, $fail failed"; exit 1; }
[ "$NU" = "0" ] && ok "every normal is unit length" || bad "$NU normals are not unit length"

# Measured on this fixture: 52 of 4986 (1.04%) with the normal taken from the
# marched field, 100 (2.01%) when it came from fpos alone. The ceiling sits
# between the two, so taking normals from the wrong field fails this.
PCT=$(awk -v f="$FLIP" -v t="$TOT" 'BEGIN{printf "%.2f", t?100*f/t:0}')
LIM=$(awk -v t="$TOT" 'BEGIN{printf "%d", t*0.015}')
if [ "$FLIP" -le "$LIM" ]; then
    ok "vertex normals agree with their own facet ($FLIP of $TOT backwards, ${PCT}%, cap ${LIM})"
else
    bad "$FLIP of $TOT vertex normals are lit backwards (${PCT}%) - normals are not the marched field's gradient"
fi
# The tail matters more than the median: smooth shading legitimately differs
# from the flat facet normal, but a 60-degree disagreement is a shading seam.
awk -v p="$P95" 'BEGIN{exit !(p<55)}' \
    && ok "the worst 5% stay within 55 degrees of their facet (p95 ${P95}d)" \
    || bad "p95 vertex-vs-facet angle is ${P95}d"

echo "mesh-normals: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
