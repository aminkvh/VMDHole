#!/bin/sh
# Three guards on the headless run path, exercised without VMD.
#
# THE DEFECTS THESE GUARD (each verified to go red on its pre-fix tree):
#
# 1. run_tunnel_analysis interpolated state(tunnel_start) into run.sh and into
#    `exec sh -c` guarded on ELEMENT COUNT alone. "12.3 4.5 nan_typo" is a
#    valid 3-element list that the engine's atof reads as z=0, and
#    {0 0 {1; touch /tmp/x}} is a valid 3-element list that lands as a shell
#    command separator. The start point must be rejected up front unless all
#    three elements are numbers, and shell_quote must neutralise metacharacters
#    at the interpolation sites.
#
# 2. run_tunnel_analysis had nine hand-written `set busy 0` restore points and
#    no catch, so a throw from resolve_output_root (an unwritable output root)
#    left `busy` stuck at 1 and _calc_depth incremented: every later run in
#    that session returned "already in progress" with no output and no error.
#    An error escaping the body must restore busy and balance _end_calc.
#
# 3. run_analysis's launch loop reuses ONE atomselect handle across frames
#    ("$sel frame N; $sel update"). The empty-selection branch counts the frame
#    as failed and continues - it must NOT delete the shared handle, or every
#    frame after the first empty one dies on "invalid command name", which is
#    the mass-discard failure that branch exists to prevent.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/vmdhole/vmdhole.tcl"

echo "headless-run-guards: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no vmdhole.tcl"; exit 0; }
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

cat > "$T/drv.tcl" <<'TCLEOF'
namespace eval ::VMDHole {}
set SRC [lindex $argv 0]; set T [lindex $argv 1]
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
foreach p {run_tunnel_analysis shell_quote} {
    set b [lift $src $p]
    if {$b eq ""} { puts "FATAL: could not lift ::VMDHole::$p"; exit 3 }
    namespace eval ::VMDHole $b
}

# Stubs: enough of the environment to reach the entry checks and the first
# throwing call, nothing more.
proc vmdcon {args} {}
proc molinfo {args} { return 10 }
set ::ncalls_begin 0; set ::ncalls_end 0
proc ::VMDHole::_begin_calc {} { incr ::ncalls_begin }
proc ::VMDHole::_end_calc {}   { incr ::ncalls_end }
proc ::VMDHole::parse_frame_spec {molid spec} { return {0} }
proc ::VMDHole::_tunnel_cfg {} { return {} }
proc ::VMDHole::resolve_output_root {molid kind} { error "simulated: unwritable output root" }

set ::VMDHole::busy 0
array set ::VMDHole::state {
    molid 0  tunnel_auto_origin 0  frame_spec all  status ""
}

proc run_with_seed {seed} {
    set ::VMDHole::state(tunnel_start) $seed
    set ::VMDHole::state(status) ""
    return [catch {::VMDHole::run_tunnel_analysis} ::_err]
}

# --- 1a. non-numeric third element is refused before anything runs ----------
set rc [run_with_seed "12.3 4.5 nan_typo"]
if {$rc == 0 && [string match "*must be three numbers*" $::VMDHole::state(status)]} {
    puts "OK typo seed is refused with a clear status"
} else {
    puts "BAD typo seed: rc=$rc status='$::VMDHole::state(status)'"
}
if {$::VMDHole::busy == 0} { puts "OK busy released after refusal" } else { puts "BAD busy stuck after refusal" }

# --- 1b. a shell-metacharacter element is refused, and nothing executes -----
set marker [file join $T pwned]
set rc [run_with_seed "0 0 {1; touch $marker}"]
if {$rc == 0 && [string match "*must be three numbers*" $::VMDHole::state(status)] && ![file exists $marker]} {
    puts "OK metacharacter seed is refused and no side effect ran"
} else {
    puts "BAD metacharacter seed: rc=$rc marker=[file exists $marker] status='$::VMDHole::state(status)'"
}

# --- 1c. shell_quote really neutralises what the sites interpolate ----------
set q [::VMDHole::shell_quote "1; touch $T/sq_marker"]
catch {exec sh -c "printf %s $q" } out
if {![file exists $T/sq_marker] && $out eq "1; touch $T/sq_marker"} {
    puts "OK shell_quote passes the hostile string through inert"
} else {
    puts "BAD shell_quote: marker=[file exists $T/sq_marker] out='$out'"
}

# --- 2. an error escaping the body restores busy and balances _end_calc -----
set before_diff [expr {$::ncalls_begin - $::ncalls_end}]
set rc [run_with_seed "1.0 2.0 3.0"]
if {$rc == 1 && [string match "*unwritable output root*" $::_err]} {
    puts "OK a valid seed reaches the throwing call and the error propagates"
} else {
    puts "BAD valid seed: rc=$rc err='$::_err'"
}
if {$::VMDHole::busy == 0} {
    puts "OK busy restored after a mid-body throw"
} else {
    puts "BAD busy stuck at $::VMDHole::busy after a mid-body throw"
}
if {($::ncalls_begin - $::ncalls_end) == $before_diff} {
    puts "OK _begin_calc/_end_calc balanced across the throw"
} else {
    puts "BAD calc depth leaked: begin=$::ncalls_begin end=$::ncalls_end"
}

# --- 3. the empty-selection branch must not delete the shared handle --------
set i [string first {lappend _empty_sel_frames $frame} $src]
if {$i < 0} { puts "BAD empty-selection branch not found in run_analysis" } else {
    set j [string first "continue" $src $i]
    set seg [string range $src $i $j]
    if {[string first {$sel delete} $seg] < 0} {
        puts "OK empty-selection branch keeps the shared atomselect alive"
    } else {
        puts "BAD empty-selection branch deletes the shared atomselect"
    }
}
TCLEOF

out=$(timeout 60 tclsh "$T/drv.tcl" "$SRC" "$T" 2>&1); rc=$?
printf '%s\n' "$out" | sed 's/^/    /'
if [ "$rc" -ne 0 ] || printf '%s\n' "$out" | grep -q '^FATAL'; then
    echo "  FAIL  the check could not run (exit $rc)"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
nok=$(printf '%s\n' "$out" | grep -c '^OK ')
nbad=$(printf '%s\n' "$out" | grep -c '^BAD ')
if [ "$nbad" -eq 0 ] && [ "$nok" -ge 8 ]; then
    echo "  PASS  headless run guards hold ($nok checks)"
    echo "  -> $nok passed, 0 failed"
else
    echo "  FAIL  $nbad check(s) failed"
    echo "  -> $nok passed, $nbad failed"; exit 1
fi
