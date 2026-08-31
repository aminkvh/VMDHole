#!/bin/sh
# A failed profile parse must not destroy the good hole_profile.tsv beside it.
#
# THE DEFECT THIS GUARDS (verified to go red on the pre-fix tree):
#   Both TSV writers (parse_profile and the worker parser built by
#   _thread_parse_initscript) write to a sibling .part and publish by rename so
#   that a parse which THROWS cannot truncate an existing profile. But the
#   rename ran unconditionally, BEFORE the points==0 test - so a stale
#   hole_out.txt that parses to zero rows still renamed a header-only .part
#   over the good TSV. That is the very outcome the .part scheme was added to
#   prevent: the zero-point return is a failed parse, not a success with an
#   empty table.
#
# Both parsers are exercised: parse_profile is lifted straight out of the
# shipped file, and the worker parser is obtained by calling the REAL
# _thread_parse_initscript and evaluating its script (minus thread::wait), so
# the test cannot drift from the code it checks. No VMD needed.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/vmdhole/vmdhole.tcl"

echo "tsv-publish-on-failure: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no vmdhole.tcl"; exit 0; }
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

# A valid plain-table HOLE output: three rows, min radius 0.8 at coord 0.
cat > "$T/good_out.txt" <<'TXTEOF'
 HOLE run header
 cenxyz.cvec     radius  cen_line_D sum{s/area}
    -2.00000     1.50000    0.10000    0.00100
    -1.00000     1.20000    0.20000    0.00200
     0.00000     0.80000    0.30000    0.00300
 trailer line ends the table
TXTEOF
# A stale/broken output: an error report, no profile table at all.
cat > "$T/bad_out.txt" <<'TXTEOF'
 ***ERROR****
 radius file has no VDWR entry for element QQ
TXTEOF
# One data row: the publish gate is points == 0 exactly, so a single-row
# profile must still publish (guards a gate mutated to any higher threshold).
cat > "$T/one_row_out.txt" <<'TXTEOF'
 cenxyz.cvec     radius  cen_line_D sum{s/area}
     1.00000     2.50000    0.40000    0.00400
TXTEOF

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
foreach p {parse_profile _resolve_conn_radii _hole_error_hint _thread_parse_initscript} {
    set b [lift $src $p]
    if {$b eq ""} { puts "FATAL: could not lift ::VMDHole::$p"; exit 3 }
    namespace eval ::VMDHole $b
}
proc vmdcon {args} {}

# The worker parser, from the real initscript: everything but the thread::wait
# it ends with (there is no worker thread here).
set ws [::VMDHole::_thread_parse_initscript]
regsub {\nthread::wait\n?$} $ws "" ws
if {[string match "*thread::wait*" $ws]} { puts "FATAL: thread::wait not stripped"; exit 3 }
eval $ws

proc readfile {f} { set h [open $f r]; set d [read $h]; close $h; return $d }

# name: label for messages; parsecmd: takes {out tsv} and returns the profile.
proc exercise {name parsecmd} {
    global T
    set dir [file join $T $name]
    file mkdir $dir
    set dest [file join $dir hole_profile.tsv]

    # 1. a successful parse publishes the TSV
    set r [{*}$parsecmd [file join $T good_out.txt] $dest]
    if {[dict get $r valid] == 1 && [dict get $r points] == 3
            && abs([dict get $r min_radius] - 0.8) < 1e-9} {
        puts "OK $name: valid output parses (3 points, min_radius 0.8)"
    } else {
        puts "BAD $name: valid output misparsed: $r"
    }
    if {[dict get $r tsv_file] eq $dest} {
        puts "OK $name: success return names the published destination"
    } else {
        puts "BAD $name: success return names '[dict get $r tsv_file]', not the destination"
    }
    if {[file exists $dest] && ![file exists $dest.part]} {
        puts "OK $name: successful parse published the TSV and removed the .part"
    } else {
        puts "BAD $name: TSV missing or .part left behind after success"
    }
    set good [readfile $dest]

    # 2. a failed parse of a stale output must leave that TSV untouched
    set r [{*}$parsecmd [file join $T bad_out.txt] $dest]
    if {[dict get $r valid] == 0 && [string match "HOLE stopped with an error*" [dict get $r message]]} {
        puts "OK $name: stale output is reported as a failed parse"
    } else {
        puts "BAD $name: stale output not reported as failure: $r"
    }
    if {[dict get $r tsv_file] eq $dest} {
        puts "OK $name: failure return names the destination, not the .part"
    } else {
        puts "BAD $name: failure return names '[dict get $r tsv_file]'"
    }
    if {[readfile $dest] eq $good} {
        puts "OK $name: failed parse left the good TSV byte-identical"
    } else {
        puts "BAD $name: failed parse overwrote the good TSV"
    }
    if {![file exists $dest.part]} {
        puts "OK $name: no .part litter after the failed parse"
    } else {
        puts "BAD $name: .part left behind after the failed parse"
    }

    # 3. a failed parse with NO existing TSV creates none
    set dest2 [file join $dir fresh_profile.tsv]
    {*}$parsecmd [file join $T bad_out.txt] $dest2
    if {![file exists $dest2] && ![file exists $dest2.part]} {
        puts "OK $name: failed parse of a fresh dir creates no TSV at all"
    } else {
        puts "BAD $name: failed parse created a file in a fresh dir"
    }

    # 4. the gate is points == 0 exactly: one row is a success and publishes
    set dest3 [file join $dir one_row_profile.tsv]
    set r [{*}$parsecmd [file join $T one_row_out.txt] $dest3]
    if {[dict get $r valid] == 1 && [dict get $r points] == 1 && [file exists $dest3]} {
        puts "OK $name: a single-row profile publishes (gate is zero, not a threshold)"
    } else {
        puts "BAD $name: single-row profile mishandled: valid=[dict get $r valid] exists=[file exists $dest3]"
    }
}

exercise serial {::VMDHole::parse_profile}
exercise worker {thread_parse_frame}
TCLEOF

out=$(timeout 60 tclsh "$T/drv.tcl" "$SRC" "$T" 2>&1); rc=$?
printf '%s\n' "$out" | sed 's/^/    /'
if [ "$rc" -ne 0 ] || printf '%s\n' "$out" | grep -q '^FATAL'; then
    echo "  FAIL  the check could not run (exit $rc)"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
nok=$(printf '%s\n' "$out" | grep -c '^OK ')
nbad=$(printf '%s\n' "$out" | grep -c '^BAD ')
if [ "$nbad" -eq 0 ] && [ "$nok" -ge 18 ]; then
    echo "  PASS  a failed parse can no longer truncate a published profile ($nok checks)"
    echo "  -> $nok passed, 0 failed"
else
    echo "  FAIL  $nbad check(s) failed"
    echo "  -> $nok passed, $nbad failed"; exit 1
fi
