#!/bin/sh
# The two TSV readers must agree. They are the same question asked twice.
#
# THE DEFECT THIS GUARDS (verified to go red on the pre-fix tree):
#   vmdhole.tcl has two readers for hole_profile.tsv. The main-thread one,
#   parse_profile_from_tsv, hands its rows to _resolve_conn_radii, which applies
#   the CONNOLLY convention (Requiv where real, interpolated across HOLE's
#   un-evaluated mid-point rows, spherical-probe radius where the pore has
#   opened to bulk). The worker-thread one, _thr_read_tsv, took column 2
#   verbatim and never called it - even though _thread_parse_initscript copies
#   that helper into the worker specifically "so the CONNOLLY radius convention
#   cannot drift".
#
#   Which reader runs is decided by frame count (>=8 frames, >1 worker, and a
#   Tcl with the Thread package). So on a CONNOLLY trajectory the reported
#   bottleneck radius could change with nothing but frame_spec, silently.
#
# Needs no VMD and no Thread package: both procs are plain Tcl, lifted out of
# the shipped file so the test cannot drift from the code it checks.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/vmdhole/vmdhole.tcl"

echo "tsv-reader-parity: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no vmdhole.tcl"; exit 0; }
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

# A CONNOLLY TSV: col2 = spherical probe radius, col5 = Requiv (the real
# answer), blank col5 on the interleaved mid-point rows HOLE does not evaluate.
# Serial reader -> min_radius 1.20 (Requiv). Buggy threaded reader -> 3.0 (col2).
printf 'coord\tradius\tcen_line_d\tsum\trequiv\tconn_s\trequiv_est\tcap_rad\n'  > "$T/p.tsv"
printf -- '-4.0\t5.0\t-4.0\t0.10\t4.00\t0.10\t4.00\n'                          >> "$T/p.tsv"
printf -- '-3.0\t4.0\t-3.0\t0.20\n'                                            >> "$T/p.tsv"
printf -- '-2.0\t3.0\t-2.0\t0.30\t1.20\t0.30\t1.20\n'                          >> "$T/p.tsv"
printf -- '-1.0\t4.0\t-1.0\t0.40\n'                                            >> "$T/p.tsv"
printf -- '0.0\t5.0\t0.0\t0.50\t4.00\t0.50\t4.00\n'                            >> "$T/p.tsv"

# Quoted heredoc: no shell expansion at all. Paths arrive via argv.
cat > "$T/drv.tcl" <<'TCLEOF'
namespace eval ::VMDHole {}
set SRC [lindex $argv 0]
set TSV [lindex $argv 1]
set fh [open $SRC r]; set src [read $fh]; close $fh

# Lift a proc out of the shipped file by brace balance, so this test always
# checks the code that actually ships. qualified=1 for ::VMDHole::name.
proc lift {src name qualified} {
    set pat [expr {$qualified ? "proc ::VMDHole::$name " : "proc $name "}]
    set i [string first $pat $src]
    if {$i < 0} { return "" }
    set buf ""
    foreach ln [split [string range $src $i end] "\n"] {
        append buf $ln "\n"
        if {[info complete $buf]} { return $buf }
    }
    return ""
}

# Every helper the two readers reach. A missing one must be fatal, not silent:
# an earlier draft of this test let parse_profile_from_tsv throw and still
# reported PASS, because it only counted MISMATCH lines.
foreach p {_resolve_conn_radii _conn_F_from_rows parse_profile_from_tsv} {
    set body [lift $src $p 1]
    if {$body eq ""} { puts "FATAL: could not lift ::VMDHole::$p"; exit 3 }
    namespace eval ::VMDHole $body
}
set body [lift $src _thr_read_tsv 0]
if {$body eq ""} { puts "FATAL: could not lift _thr_read_tsv"; exit 3 }
namespace eval ::VMDHole $body

if {[catch {::VMDHole::parse_profile_from_tsv $TSV 1} a]} {
    puts "FATAL: serial reader threw: $a"; exit 3
}
if {[catch {::VMDHole::_thr_read_tsv $TSV 1} b]} {
    puts "FATAL: threaded reader threw: $b"; exit 3
}
set bd [dict create {*}$b]

set bad 0; set checked 0
foreach key {points min_radius min_coord yvalues rsources} {
    set va [expr {[dict exists $a $key]  ? [dict get $a $key]  : "<missing>"}]
    set vb [expr {[dict exists $bd $key] ? [dict get $bd $key] : "<missing>"}]
    incr checked
    if {$va eq $vb} {
        puts "OK       $key = $va"
    } else {
        puts "MISMATCH $key  serial=$va  threaded=$vb"
        incr bad
    }
}
# Guard the guard: the fixture must actually exercise the convention, i.e. the
# resolved radius must differ from raw column 2. Otherwise this test could pass
# on a fixture where both readers trivially agree.
if {[dict get $a min_radius] != 1.2} {
    puts "FATAL: fixture no longer exercises the CONNOLLY convention (serial min_radius = [dict get $a min_radius], expected 1.2)"
    exit 3
}
puts "CHECKED $checked BAD $bad"
exit [expr {$bad ? 1 : 0}]
TCLEOF

out=$(timeout 120 tclsh "$T/drv.tcl" "$SRC" "$T/p.tsv" 2>&1); rc=$?
echo "$out" | sed 's/^/    /'

if [ "$rc" -eq 3 ] || printf '%s\n' "$out" | grep -q '^FATAL'; then
    echo "  FAIL  the comparison could not run (see FATAL above)"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
if [ "$rc" -ne 0 ]; then
    echo "  FAIL  the serial and threaded TSV readers disagree"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
echo "  PASS  both TSV readers agree on a CONNOLLY profile"
echo "  -> 1 passed, 0 failed"
