#!/bin/sh
# Is the inlined HOLE engine inside vmdhole.tcl still what vmdhole/hole_tcl/ generates?
#
# The plugin ships as ONE script, so vmdhole/hole_tcl/ is the source and the region
# between the sentinels in vmdhole.tcl is a GENERATED copy. Nothing enforced
# that: editing vmdhole/hole_tcl/ without re-running sync_into_plugin.py left the
# shipped engine stale, and editing the inlined region directly left vmdhole/hole_tcl/
# stale - in both cases silently, and in both cases the "source of truth"
# becomes whichever copy the reader happens to open.
#
# This re-inlines into a COPY and requires the result to be byte-identical to
# what is committed. It never writes to the real file.
#
# (the enforceable half. The reviewer's wider
# point, that MOLE should have module sources too, is a refactor and is not
# addressed here: MOLE has no equivalent of vmdhole/hole_tcl/, so its only source is
# the inlined text itself.)
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$DIR/../.." && pwd)
PLUGIN="$REPO/vmdhole/vmdhole.tcl"
SYNC="$REPO/vmdhole/hole_tcl/sync_into_plugin.py"

command -v python3 >/dev/null 2>&1 || { echo "SKIP: no python3"; exit 0; }
[ -f "$SYNC" ]   || { echo "SKIP: vmdhole/hole_tcl/sync_into_plugin.py not found"; exit 0; }
[ -f "$PLUGIN" ] || { echo "  FAIL  vmdhole.tcl not found at $PLUGIN"; exit 1; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp "$PLUGIN" "$TMP/vmdhole.tcl"

if ! python3 "$SYNC" "$TMP/vmdhole.tcl" > "$TMP/log" 2>&1; then
    echo "  FAIL  sync_into_plugin.py could not re-inline"
    sed 's/^/       /' "$TMP/log"
    exit 1
fi

if cmp -s "$PLUGIN" "$TMP/vmdhole.tcl"; then
    echo "  PASS  the inlined HOLE engine matches vmdhole/hole_tcl/ ($(sed -n 's/re-inlined \([0-9]*\) lines.*/\1/p' "$TMP/log") lines)"
    echo "inline_current: 1 passed, 0 failed"
    exit 0
fi

echo "  FAIL  vmdhole.tcl's inlined HOLE engine is NOT what vmdhole/hole_tcl/ generates"
echo "        One of the two was edited without the other. Re-inline with:"
echo "            python3 vmdhole/hole_tcl/sync_into_plugin.py vmdhole/vmdhole.tcl"
echo "        ...but check FIRST which copy holds the change you want - this"
echo "        overwrites the inlined region from vmdhole/hole_tcl/."
diff "$PLUGIN" "$TMP/vmdhole.tcl" | head -25 | sed 's/^/       /'
echo "inline_current: 0 passed, 1 failed"
exit 1
