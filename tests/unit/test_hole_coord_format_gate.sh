#!/bin/sh
# The packed coordinate record may only be handed to a binary that can read it.
#
# THE DEFECT THIS GUARDS:
#   _hole_fast_coord_available used to open the packed path whenever the
#   acceleration manifest merely NAMED the fast-coord-read patch. The manifest
#   says which patches a binary carries, not which record FORMAT it speaks, so a
#   binary built before a format change satisfied the gate and was then handed a
#   record it could not parse. Its reader compares the 8-byte magic and, on a
#   mismatch, falls through to PDB parsing instead of erroring - so the run
#   produced a plausible profile from coordinates nobody supplied.
#
#   This is not hypothetical: the VMDHole -> VMDPathFinder rename changed the
#   magic from VMDHOLEC to VMDPFC01, and only an unrelated manifest RENAME
#   (which happened to disable the path outright) kept the two apart.
#
# The build script now stamps "coordfmt" into the manifest, read out of the
# reader's own source, and the gate requires it to equal what the writer emits.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
READER="$ROOT/native/connolly_patches/tsatr_fast.f"
BUILD="$ROOT/native/build-vmdpathfinder-optimized.sh"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "hole-coord-format-gate"
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

# --- 1. ONE definition of the magic in Tcl, and the writer uses it ----------
nlit=$(grep -c '"VMDPFC01"' "$TCL" || true)
if [ "$nlit" -eq 1 ]; then
    ok "the magic is a literal in exactly one place (_hole_coord_magic)"
else
    bad "the magic appears as a literal $nlit times in the plugin; expected 1"
fi
if grep -q 'puts -nonewline \$fh \[_hole_coord_magic\]' "$TCL"; then
    ok "the writer stamps the shared constant, not its own literal"
else
    bad "_write_hole_coord_bin does not use _hole_coord_magic"
fi

# --- 2. Tcl's magic equals the Fortran reader's ------------------------------
if [ -f "$READER" ]; then
    fmagic=$(sed -n "s/.*VHMAG\.NE\.'\([A-Z0-9]*\)'.*/\1/p" "$READER" | head -1)
    tmagic=$(sed -n 's/.*return "\(VMD[A-Z0-9]*\)".*/\1/p' "$TCL" | head -1)
    if [ -n "$fmagic" ] && [ "$fmagic" = "$tmagic" ]; then
        ok "writer and reader agree on the magic ($tmagic)"
    else
        bad "writer says '$tmagic', reader says '$fmagic'"
    fi
else
    echo "  note: no reader source, skipping the cross-language check"
fi

# --- 3. the build script stamps coordfmt, derived from the reader -----------
if grep -q 'echo "coordfmt' "$BUILD" && grep -q 'tsatr_fast.f' "$BUILD"; then
    ok "the build script stamps coordfmt from the reader's own source"
else
    bad "the build script does not stamp coordfmt"
fi

# --- 4. THE GATE: it must REFUSE a manifest with no/!= coordfmt -------------
# Drives the real proc against three manifests, so this is the behaviour and
# not a grep of the source.
command -v tclsh >/dev/null 2>&1 || { echo "  note: no tclsh, skipping gate exercise"
    echo "hole-coord-format-gate: $pass passed, $fail failed"; [ "$fail" -eq 0 ] || exit 1; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM
mk() { # file, coordfmt-line
    mkdir -p "$T/$1"; : > "$T/$1/hole"
    { echo "vmdpathfinder_accel_manifest 1"; echo "openmp yes"
      echo "patch tsatr_fast.f fast-coord-read"; [ -n "$2" ] && echo "$2"
    } > "$T/$1/vmdpathfinder_accel.manifest"
}
mk good  "coordfmt VMDPFC01"
mk stale ""
mk wrong "coordfmt VMDHOLEC"

out=$(tclsh <<TCLEOF 2>&1
# Only the two procs under test plus the reader they use - sourcing the whole
# plugin needs VMD, which a unit test does not have.
namespace eval ::VMDPathFinder {variable state}
proc ::VMDPathFinder::_accel_manifest_path {} {
    variable state
    return [file join [file dirname \$state(hole_exec)] vmdpathfinder_accel.manifest]
}
$(awk '/^proc ::VMDPathFinder::_read_accel_manifest/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_hole_fast_coord_available/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_hole_coord_magic/,/^}/' "$TCL")
foreach d {good stale wrong} {
    set ::VMDPathFinder::state(hole_exec) "$T/\$d/hole"
    puts "\$d [::VMDPathFinder::_hole_fast_coord_available]"
}
TCLEOF
)
chk() { got=$(echo "$out" | awk -v k="$1" '$1==k{print $2}')
    if [ "$got" = "$2" ]; then ok "$3"; else bad "$3 (got '$got', wanted '$2')"; fi; }
chk good  1 "a matching coordfmt OPENS the packed path"
chk stale 0 "a manifest predating the stamp is refused"
chk wrong 0 "a manifest naming the OLD magic is refused"

echo "hole-coord-format-gate: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
