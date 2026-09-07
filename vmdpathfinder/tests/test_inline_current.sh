#!/bin/sh
# Is the inlined HOLE engine inside vmdpathfinder.tcl still what vmdpathfinder/hole_tcl/ generates?
#
# The plugin ships as ONE script, so vmdpathfinder/hole_tcl/ is the source and the region
# between the sentinels in vmdpathfinder.tcl is a GENERATED copy. Nothing enforced
# that: editing vmdpathfinder/hole_tcl/ without re-running sync_into_plugin.py left the
# shipped engine stale, and editing the inlined region directly left vmdpathfinder/hole_tcl/
# stale - in both cases silently, and in both cases the "source of truth"
# becomes whichever copy the reader happens to open.
#
# This re-inlines into a COPY and requires the result to be byte-identical to
# what is committed. It never writes to the real file.
#
# (the enforceable half. The reviewer's wider
# point, that MOLE should have module sources too, is a refactor and is not
# addressed here: MOLE has no equivalent of vmdpathfinder/hole_tcl/, so its only source is
# the inlined text itself.)
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$DIR/../.." && pwd)
PLUGIN="$REPO/vmdpathfinder/vmdpathfinder.tcl"
SYNC="$REPO/vmdpathfinder/hole_tcl/sync_into_plugin.py"

command -v python3 >/dev/null 2>&1 || { echo "SKIP: no python3"; exit 0; }
[ -f "$SYNC" ]   || { echo "SKIP: vmdpathfinder/hole_tcl/sync_into_plugin.py not found"; exit 0; }
[ -f "$PLUGIN" ] || { echo "  FAIL  vmdpathfinder.tcl not found at $PLUGIN"; exit 1; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp "$PLUGIN" "$TMP/vmdpathfinder.tcl"

if ! python3 "$SYNC" "$TMP/vmdpathfinder.tcl" > "$TMP/log" 2>&1; then
    echo "  FAIL  sync_into_plugin.py could not re-inline"
    sed 's/^/       /' "$TMP/log"
    exit 1
fi

if cmp -s "$PLUGIN" "$TMP/vmdpathfinder.tcl"; then
    echo "  PASS  the inlined HOLE engine matches vmdpathfinder/hole_tcl/ ($(sed -n 's/re-inlined \([0-9]*\) lines.*/\1/p' "$TMP/log") lines)"
    echo "inline_current: 1 passed, 0 failed"
    exit 0
fi

echo "  FAIL  vmdpathfinder.tcl's inlined HOLE engine is NOT what vmdpathfinder/hole_tcl/ generates"
echo "        One of the two was edited without the other. Re-inline with:"
echo "            python3 vmdpathfinder/hole_tcl/sync_into_plugin.py vmdpathfinder/vmdpathfinder.tcl"
echo "        ...but check FIRST which copy holds the change you want - this"
echo "        overwrites the inlined region from vmdpathfinder/hole_tcl/."
diff "$PLUGIN" "$TMP/vmdpathfinder.tcl" | head -25 | sed 's/^/       /'
echo "inline_current: 0 passed, 1 failed"
exit 1
