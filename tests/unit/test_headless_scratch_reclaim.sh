#!/bin/sh
# Scratch reclaim must be reachable without opening the GUI.
#
# THE DEFECT THIS GUARDS:
#   _sweep_stale_tmpdirs had exactly one call site, inside show_gui, and
#   cleanup_temporary_outputs exactly one, inside close_gui. init_executables -
#   the entry point a batch job actually uses - called neither, and
#   docs/scripting.md tells batch authors not to call show_gui. A run killed by
#   the scheduler therefore left its /tmp and /dev/shm scratch behind with
#   nothing in the product that would ever collect it. /dev/shm is RAM.
#
# Checks both halves: that the sweeper is wired into the headless entry point,
# and that it still deletes only what it should once it runs.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/vmdhole/vmdhole.tcl"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "headless-scratch-reclaim: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no vmdhole.tcl"; exit 0; }
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
[ -d /proc ] || { echo "SKIP: sweeper is /proc-gated and this host has none"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

cat > "$T/drv.tcl" <<'TCLEOF'
namespace eval ::VMDHole {}
set SRC [lindex $argv 0]
set SB  [lindex $argv 1]
set fh [open $SRC r]; set src [read $fh]; close $fh

proc lift {src name} {
    set pat "proc ::VMDHole::$name "
    set i [string first $pat $src]
    if {$i < 0} { return "" }
    set buf ""
    foreach ln [split [string range $src $i end] "\n"] {
        append buf $ln "\n"
        if {[info complete $buf]} { return $buf }
    }
    return ""
}

# --- reachability: the headless entry point must call the sweeper -----------
set ie [lift $src init_executables]
if {$ie eq ""} { puts "FATAL: could not lift init_executables"; exit 3 }
if {[string match {*_sweep_stale_tmpdirs*} $ie]} {
    puts "OK reachable: init_executables calls _sweep_stale_tmpdirs"
} else {
    puts "BAD reachable: init_executables never calls _sweep_stale_tmpdirs"
}

# --- one-shot at the CALL SITE, not inside the sweeper ----------------------
# The sweeper itself must stay callable repeatedly: headless_smoke.tcl drives it
# three times to check the age guard from both sides, so a guard inside it would
# make that test's third call a no-op.
set sw [lift $src _sweep_stale_tmpdirs]
if {$sw eq ""} { puts "FATAL: could not lift _sweep_stale_tmpdirs"; exit 3 }
if {[string match {*_swept*} $ie] && ![string match {*_swept*} $sw]} {
    puts "OK oneshot: guarded at the entry point, sweeper itself still repeatable"
} elseif {[string match {*_swept*} $sw]} {
    puts "BAD oneshot: the guard is INSIDE the sweeper - repeat calls become no-ops"
} else {
    puts "BAD oneshot: no guard at the entry point"
}

# --- behaviour: run it against a sandbox, not the real /tmp -----------------
# Redirect the hardcoded base list. If this substitution ever stops matching,
# the test would silently run against the real /tmp - so it is verified.
set sb [string map [list "foreach base {/tmp /dev/shm}" "foreach base [list [list $SB]]"] $sw]
if {$sb eq $sw} { puts "FATAL: could not redirect the sweeper's base list"; exit 3 }
variable ::VMDHole::_swept 0
namespace eval ::VMDHole $sb

set dead 999999
while {[file isdirectory /proc/$dead]} { incr dead }
set live [pid]

set old [expr {[clock seconds] - 7200}]
foreach nm [list vmdhole_scratch_$dead vmdhole_$dead vmdhole_${dead}_f7 \
                 vmdhole_clip_${dead}_1699999999999 vmdhole_scratch_$live] {
    file mkdir [file join $SB $nm]
    file mtime [file join $SB $nm] $old
}
# A fresh dir belonging to a dead pid must survive the one-hour age guard.
file mkdir [file join $SB vmdhole_scratch_${dead}_fresh]

::VMDHole::_sweep_stale_tmpdirs

foreach {nm want} [list \
        vmdhole_scratch_$dead              gone \
        vmdhole_$dead                      gone \
        vmdhole_${dead}_f7                 gone \
        vmdhole_clip_${dead}_1699999999999 gone \
        vmdhole_scratch_$live              kept \
        vmdhole_scratch_${dead}_fresh      kept] {
    set there [file isdirectory [file join $SB $nm]]
    set got [expr {$there ? "kept" : "gone"}]
    if {$got eq $want} { puts "OK sweep: $nm -> $got" } \
                  else { puts "BAD sweep: $nm -> $got, wanted $want" }
}
TCLEOF

mkdir -p "$T/sandbox"
out=$(timeout 60 tclsh "$T/drv.tcl" "$SRC" "$T/sandbox" 2>&1); rc=$?
if [ "$rc" -eq 3 ] || printf '%s\n' "$out" | grep -q '^FATAL'; then
    printf '%s\n' "$out" | sed 's/^/    /'
    echo "  FAIL  the check could not run (see FATAL above)"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
printf '%s\n' "$out" | while IFS= read -r l; do echo "    $l"; done
nok=$(printf '%s\n' "$out" | grep -c '^OK ')
nbad=$(printf '%s\n' "$out" | grep -c '^BAD ')
[ "$nbad" -eq 0 ] && ok "$nok checks: reachable headless, one-shot, and sweeps correctly" \
                  || bad "$nbad of $((nok+nbad)) checks failed"
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
