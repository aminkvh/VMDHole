#!/bin/sh
# Tunnel Ion Flow must not reference the pore branch's variables.
#
# THE DEFECT THIS GUARDS:
#   The azimuth added for opening attribution was computed from qx/qy/qz, which
#   are assigned ONLY in the pore branch. In tunnel mode `_ion_flow_path_coord`
#   returns z and R and nothing else, so every accepted ion hit an undefined
#   variable and the scan died with a Tcl error. The tests shipped alongside it
#   checked GUI column placement and never executed the calculation.
#
# Source-level because the failing line is inside a 700-line proc that needs a
# loaded trajectory to reach; what matters is that the tunnel branch defines
# whatever the shared code path consumes.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "ion-flow-azimuth"
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

BODY=$(awk '/^proc ::VMDPathFinder::_ion_flow_scan/,/^}/' "$TCL")

# 1. the azimuth must not be taken from q* at the shared append site
if printf '%s' "$BODY" | grep -q 'lset tr_az \$idx \[linsert \[lindex \$tr_az \$idx\] end \$_azv\]'; then
    ok "the ion azimuth is appended from a variable both branches set"
else
    bad "the shared append does not use a per-branch azimuth variable"
fi

# 2. BOTH branches must assign it
# The proc has SEVEN "if {$_flow_tunnel}" branches; the one that matters is the
# per-ion coordinate branch, anchored on its own path_coord call.
TUN=$(printf '%s\n' "$BODY" | awk '/_ion_flow_path_coord .*\] z R/{f=1} f{print; n++} n>8{exit}')
if printf '%s' "$TUN" | grep -q 'set _azv'; then
    ok "the tunnel branch assigns the azimuth it will append"
else
    bad "the tunnel branch does not assign _azv - the pore branch's q* would leak in"
fi
if printf '%s' "$BODY" | grep -c 'set _azv' | grep -qE '^[2-9]'; then
    ok "both branches assign it"
else
    bad "only one branch assigns _azv"
fi

# 3. a tunnel sample must carry NO azimuth, not a fabricated 0.0 - a zero angle
#    would be attributed to whichever opening happens to sit at zero
if printf '%s' "$TUN" | grep -q 'set _azv ""'; then
    ok "a tunnel sample carries an EMPTY azimuth, not a fabricated angle"
else
    bad "the tunnel branch fabricates an azimuth"
fi

# 4. and the consumer must skip empties
OCC=$(awk '/^proc ::VMDPathFinder::_conn_opening_occupancy/,/^}/' "$TCL")
if printf '%s' "$OCC" | grep -q 'if {\$th ne "" && \$t ne ""'; then
    ok "the opening attribution skips samples with no azimuth"
else
    bad "the attribution does not guard against an empty azimuth"
fi

# 5. attribution must be PER FRAME, not one frame's map applied to all
if printf '%s' "$OCC" | grep -q 'dict exists \$permap \$f' && \
   printf '%s' "$OCC" | grep -q 'foreach fr \$result_frames'; then
    ok "each sample is attributed with its OWN frame's opening map"
else
    bad "attribution is not per-frame - one frame's geometry would be applied to all"
fi

# 6. episodes must require consecutive frames
if printf '%s' "$OCC" | grep -q 'f == \$prev_f + 1'; then
    ok "a visit is consecutive frames, so a gap starts a new one"
else
    bad "episodes do not require adjacent frames - gaps would count as dwell"
fi

# 7. water must not be counted as ion traffic
if printf '%s' "$OCC" | grep -q 'species\] eq "Water"} { continue }'; then
    ok "water is excluded from the ion-frame counts"
else
    bad "water contributes to a count the UI calls ion-frames"
fi

echo "ion-flow-azimuth: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
