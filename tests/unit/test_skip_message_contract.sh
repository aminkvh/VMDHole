#!/bin/sh
# Invariant: a test that skips the WHOLE GROUP must say so in the form its
# runner can detect.
#
# THE DEFECT THIS GUARDS (verified to go red on the pre-fix tree):
#   vmdhole/tests/run_tests.sh decides a group "checked nothing" by grepping its
#   output for an anchored `^SKIP:`. Under VMDHOLE_RELEASE=1 that is what turns a
#   skipped group into a release-blocking failure.
#
#   hole_tcl_fallback.tcl announced its skip as "  SKIP hole_tcl_fallback: ..."
#   - indented, with the colon after the label rather than after SKIP. That
#   matches neither form, so whenever that group skipped it was silently counted
#   as a PASS and "ALL 20 TEST GROUPS PASSED" overstated what had run.
#
#   run_tests.sh's own header records this same defect being fixed once before,
#   for test_gui_reachable and test_headless_smoke. It recurred because nothing
#   enforced the contract; this enforces it for every group at once.
#
# WHAT COUNTS AS GROUP-LEVEL, precisely:
#   A per-assertion "  SKIP  this one sub-check" is legitimate and deliberately
#   NOT matched by run_tests.sh - run_tests.sh's own comments say so, and
#   flagging those would be wrong. The two are told apart by what follows: a
#   GROUP-level skip abandons the run, so it is followed by an `exit` within a
#   few lines. Only those must use the anchored form.
#
#   Comment lines are ignored entirely. An earlier draft of this test grepped
#   for any line containing SKIP and flagged the files' own documentation -
#   including this header. A test that fails on prose is not reliable.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TDIR="$ROOT/vmdhole/tests"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "skip-contract: $TDIR"
[ -d "$TDIR" ] || { echo "SKIP: no vmdhole/tests directory"; exit 0; }

# 1. The runner still gates the way this test assumes. If run_tests.sh changes
#    its detection, fail loudly rather than keep enforcing a dead contract.
if grep -q '\^SKIP:' "$TDIR/run_tests.sh" 2>/dev/null; then
    ok "run_tests.sh still gates on an anchored ^SKIP:"
else
    bad "run_tests.sh no longer greps for ^SKIP: - this test's premise is stale"
fi

# 2. Every group-level skip must be emitted in the detectable form.
REPORT=$(mktemp)
for f in "$TDIR"/*.sh "$TDIR"/*.tcl; do
    [ -f "$f" ] || continue
    [ "$(basename -- "$f")" = "run_tests.sh" ] && continue
    # Test the EMITTED STRING, not the source line. Most skips are written as
    #   [ -n "$BIN" ] || { echo "SKIP: no vmd on PATH"; exit 0; }
    # so the line does not begin with the command, and a grep pattern like
    #   grep -E '^  (PASS|FAIL|SKIP|...)'
    # mentions SKIP without emitting one. Pull the double-quoted literals out
    # (odd segments of a split on ") and consider only those that BEGIN with a
    # skip announcement - that admits both "SKIP: ..." and the buggy
    # "  SKIP <label>: ..." while ignoring patterns and prose.
    awk -v FNAME="$(basename -- "$f")" '
      { line = $0; buf[NR] = line
        stripped = line; sub(/^[ \t]+/, "", stripped)
        if (stripped ~ /^#/) next                     # comments are not emissions
        n = split(line, seg, "\"")
        for (i = 2; i <= n; i += 2) {                 # even indices are inside quotes
          lit = seg[i]
          probe = lit; sub(/^[ \t]+/, "", probe)
          if (probe ~ /^SKIP/) { skiplit[NR] = lit; skipsrc[NR] = stripped }
        }
      }
      END {
        for (k in skiplit) {
          # Group-level iff the skip is followed by a bare `exit 0` within two
          # lines - the "abandon the run" idiom. A skip that instead exits with
          # a COMPUTED status (`exit $?`, `exit [expr ...]`) has printed its own
          # pass/fail summary first: the group DID run and merely skipped one
          # sub-check, which run_tests.sh deliberately does not count as a group
          # skip. Distinguishing on `exit 0` separates the two exactly; a looser
          # "any exit nearby" rule flagged four legitimate per-assertion skips.
          grouplevel = 0
          for (j = k; j <= k + 2; j++)
            if (j in buf && buf[j] ~ /(^|[ \t;{&|])exit[ \t]+0[ \t]*(;|\}|$)/) grouplevel = 1
          if (!grouplevel) continue                   # per-assertion skip: legitimate
          if (skiplit[k] !~ /^SKIP:/)
            printf "%s|%d|%s\n", FNAME, k, substr(skipsrc[k], 1, 74)
        }
      }
    ' "$f" >> "$REPORT"
done

n_bad=$(wc -l < "$REPORT" | tr -d ' ')
if [ "${n_bad:-0}" -eq 0 ]; then
    ok "every group-level skip message starts with SKIP: at column 0"
else
    bad "$n_bad group-level skip message(s) the release gate cannot see"
    sed 's/^/        /' "$REPORT" | head -8
fi
rm -f "$REPORT"

# 3. The specific regression, pinned by name, so its fix cannot be reverted
#    unnoticed even if the scan above is later loosened.
FB="$TDIR/hole_tcl_fallback.tcl"
if [ -f "$FB" ]; then
    if grep -q 'puts "SKIP: hole_tcl_fallback' "$FB"; then
        ok "hole_tcl_fallback announces its skip in the detectable form"
    else
        bad "hole_tcl_fallback no longer emits 'SKIP: hole_tcl_fallback ...'"
    fi
fi

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
