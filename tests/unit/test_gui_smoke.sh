#!/bin/sh
# GUI smoke + regression suite, driven WITHOUT VMD: the real vmdhole.tcl is
# sourced under plain tclsh+Tk with VMD stubbed, the full widget tree built by
# the real show_gui, and the defects below exercised as scripted user actions.
#
# THE DEFECTS THIS GUARDS (each verified red on its pre-fix tree):
#   close path   - close_gui self-cancelled its own abort (a close mid-run
#                  stopped nothing); stripped four display/property traces
#                  nothing re-adds (controls dead after one reopen); left a
#                  pending overwrite-confirm grab unanswered (all of VMD's Tk
#                  mouse-dead); ignored tunnel results in the keep-vis gate.
#   guards       - Run/Reset/Align/import were deliverable mid-pass (busy
#                  covered only the two Run pipelines, not the bracketed read
#                  passes); align had no bracket at all; run_current_mode's
#                  backstop double-released _calc_depth.
#   molid class  - resolve_molid throws and never returns -1, so a deleted
#                  molecule crashed Passability (then left it permanently,
#                  silently dead), Ion Flow, Permeation and the metrics
#                  readout instead of their own friendly branches; the Mean
#                  tab's Show 3D died reading a never-set _volmode.
#   tunnel/opts  - a stale gear popup wrote edits onto whichever route held
#                  its captured rank; absent-route gear clicks were silent;
#                  "nan" in numeric option fields threw from inside expr
#                  ternaries; Over Time Compute swallowed every error.
#
# Each section runs in a FRESH tclsh so probe stubs cannot leak between them.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
export VMDHOLE_ROOT="$ROOT"

echo "gui-smoke: $ROOT/vmdhole/vmdhole.tcl"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
if ! echo 'if {[catch {package require Tk}]} {exit 1}; exit 0' | tclsh >/dev/null 2>&1; then
    echo "SKIP: no usable Tk/DISPLAY (try xvfb-run)"; exit 0
fi

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM
cat > "$T/harness.tcl" <<'HARNESS_EOF'
# Shared GUI-review harness: source vmdhole.tcl under plain tclsh+Tk with VMD
# stubs and build the real GUI. Usage:  tclsh harness.tcl <script.tcl>
# The child script runs AFTER show_gui with the .vmdhole toplevel built.
package require Tk
wm withdraw .

# ---- VMD stubs --------------------------------------------------------------
set ::STUB_MOLS {}          ;# list of {molid name numframes}
set ::TKDIALOG_LOG {}       ;# every dialog invocation is recorded, not shown
proc vmdcon {args} { lappend ::VMDCON_LOG $args }
set ::VMDCON_LOG {}
proc vmd_install_extension {args} {}
proc molinfo {args} {
    set a0 [lindex $args 0]
    if {$a0 eq "list"} { set r {}; foreach m $::STUB_MOLS {lappend r [lindex $m 0]}; return $r }
    if {$a0 eq "top"}  { return [expr {[llength $::STUB_MOLS] ? [lindex $::STUB_MOLS 0 0] : -1}] }
    foreach m $::STUB_MOLS { if {[lindex $m 0] == $a0} {
        set what [lindex $args 2]
        switch -- $what {
            name      { return [lindex $m 1] }
            numframes { return [lindex $m 2] }
            frame     { return 0 }
            default   { return 0 }
        }
    }}
    error "molinfo: molecule $a0 does not exist"
}
set ::aselc 0
proc atomselect {molid seltext args} {
    if {![llength $::STUB_MOLS]} { error "atomselect: no molecules loaded" }
    set name ::asel[incr ::aselc]
    proc $name {args} [format {
        set a0 [lindex $args 0]
        switch -- $a0 {
            num    { return 0 }
            get    { return {} }
            frame  { return }
            update { return }
            delete { rename %s "" }
            default { return {} }
        }
    } $name]
    return $name
}
foreach c {mol animate graphics display material render axes light imd measure} {
    proc $c {args} { return 0 }
}
proc color {args} { return 0 }
# Tk dialogs: record and return a benign answer instead of blocking
proc tk_messageBox {args} { lappend ::TKDIALOG_LOG [list messageBox $args]; return ok }
proc tk_chooseColor {args} { lappend ::TKDIALOG_LOG [list chooseColor $args]; return "" }
proc tk_getSaveFile {args} { lappend ::TKDIALOG_LOG [list getSaveFile $args]; return "" }
proc tk_getOpenFile {args} { lappend ::TKDIALOG_LOG [list getOpenFile $args]; return "" }
proc tk_chooseDirectory {args} { lappend ::TKDIALOG_LOG [list chooseDirectory $args]; return "" }

# ---- plugin -----------------------------------------------------------------
namespace eval ::VMDHole {}
source [file join $::env(VMDHOLE_ROOT) vmdhole vmdhole.tcl]
# keep config I/O away from the user's real ~/.vmdhole_config
set ::VMDHole::config_file [file join [pwd] .harness_vmdhole_config]
catch {file delete $::VMDHole::config_file}

::VMDHole::show_gui
update

# ---- helpers for review scripts --------------------------------------------
proc walk {w} { set out [list $w]; foreach c [winfo children $w] { lappend out {*}[walk $c] }; return $out }
proc all_widgets {} { walk .vmdhole }
proc invokables {} {
    set out {}
    foreach w [all_widgets] {
        set cls [winfo class $w]
        if {$cls in {Button TButton Checkbutton TCheckbutton Radiobutton TRadiobutton Menubutton TMenubutton}} {
            lappend out $w $cls
        }
    }
    return $out
}


source [lindex $argv 0]
exit 0
HARNESS_EOF

pass=0; fail=0
run_section () {
    name=$1; script=$2
    out=$(cd "$ROOT" && timeout 120 tclsh "$T/harness.tcl" "$script" 2>&1); rc=$?
    printf '%s\n' "$out" | grep -E ' OK$| BAD' | sed 's/^/    /'
    p=$(printf '%s\n' "$out" | grep -c ' OK$'); f=$(printf '%s\n' "$out" | grep -cE ' BAD')
    [ "$rc" -ne 0 ] && { echo "    FAIL $name: harness exit $rc"; f=$((f+1)); }
    pass=$((pass+p)); fail=$((fail+f))
}

cat > "$T/s1.tcl" <<'S1_EOF'
# 1. traces survive close/reopen
set t0 [llength [trace info variable ::VMDHole::state(display_mode)]]
::VMDHole::close_gui
::VMDHole::show_gui
update
set t1 [llength [trace info variable ::VMDHole::state(display_mode)]]
set t2 [llength [trace info variable ::VMDHole::state(hydro_scheme)]]
puts "TRACES before=$t0 after_reopen=$t1 hydro=$t2 [expr {$t1>0 && $t2>0 ? {OK} : {BAD}}]"

# 2. abort survives close only while work is in flight
set ::VMDHole::busy 1; set ::VMDHole::_calc_depth 1
::VMDHole::close_gui
puts "ABORT inflight: raised=[set ::VMDHole::state(abort_requested)] [expr {$::VMDHole::state(abort_requested)==1 ? {OK} : {BAD}}]"
set ::VMDHole::busy 0
::VMDHole::_end_calc
puts "ABORT after unwind: [set ::VMDHole::state(abort_requested)] depth=$::VMDHole::_calc_depth [expr {$::VMDHole::state(abort_requested)==0 && $::VMDHole::_calc_depth==0 ? {OK} : {BAD}}]"
::VMDHole::show_gui; update
::VMDHole::close_gui
puts "ABORT idle close: [set ::VMDHole::state(abort_requested)] [expr {$::VMDHole::state(abort_requested)==0 ? {OK} : {BAD}}]"

# 3. a pending modal answer is forced to Cancel by close
::VMDHole::show_gui; update
unset -nocomplain ::VMDHole::_overwrite_ans
after 200 ::VMDHole::close_gui
set ::VMDHole::_overwrite_ans_pending 1
# simulate the dialog's wait: tkwait returns only if close_gui answers
after 2000 {set ::VMDHole::_overwrite_ans timeout}
tkwait variable ::VMDHole::_overwrite_ans
puts "OVERWRITE answered=[set ::VMDHole::_overwrite_ans] [expr {$::VMDHole::_overwrite_ans eq 0 ? {OK} : {BAD}}]"
S1_EOF
cat > "$T/s2.tcl" <<'S2_EOF'
# depth>0 (bracketed pass) must refuse Run / Reset / Align / molid swap
set ::VMDHole::busy 0
::VMDHole::_begin_calc
set r [::VMDHole::run_analysis]
puts "RUN mid-pass: rc=$r status='$::VMDHole::state(status)' [expr {$r==0 && [string match {*in progress*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"
set ::VMDHole::state(status) ""
::VMDHole::reset_session
puts "RESET mid-pass: [expr {[string match {*in progress*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"
set ::VMDHole::state(status) ""
::VMDHole::do_align_trajectory
puts "ALIGN mid-pass: [expr {[string match {*in progress*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"
::VMDHole::_end_calc

# align brackets busy/_begin_calc and restores on throw
set ::seen_busy -1; set ::seen_depth -1
proc ::VMDHole::align_trajectory {{molid ""}} {
    set ::seen_busy $::VMDHole::busy
    set ::seen_depth $::VMDHole::_calc_depth
    error "boom"
}
set ::TKDIALOG_LOG {}
::VMDHole::do_align_trajectory
puts "ALIGN bracket: during busy=$::seen_busy depth=$::seen_depth after busy=$::VMDHole::busy depth=$::VMDHole::_calc_depth dialogs=[llength $::TKDIALOG_LOG] [expr {$::seen_busy==1 && $::seen_depth==1 && $::VMDHole::busy==0 && $::VMDHole::_calc_depth==0 && [llength $::TKDIALOG_LOG]==1 ? {OK} : {BAD}}]"

# run_current_mode must not double-release an enclosing bracket
proc ::VMDHole::analysis_mode {} { return tunnel }
proc ::VMDHole::run_tunnel_analysis {} {
    # model the real wrapper: restore, then rethrow
    set ::VMDHole::busy 0
    ::VMDHole::_end_calc
    error "engine died"
}
::VMDHole::_begin_calc            ;# outer bracket (depth 1)
set ::VMDHole::busy 1
::VMDHole::_begin_calc            ;# the run's own bracket (depth 2)
::VMDHole::run_current_mode
puts "BACKSTOP: depth=$::VMDHole::_calc_depth (outer bracket must survive) [expr {$::VMDHole::_calc_depth==1 ? {OK} : {BAD}}]"
::VMDHole::_end_calc
S2_EOF
cat > "$T/s3.tcl" <<'S3_EOF'
# state: results empty, no molecules (STUB_MOLS empty) - the delete-and-move-on world
set ::VMDHole::state(molid) 0

# 1. _volmode: ticking Show 3D with no results must not raise, must untick + report
set ::VMDHole::state(show_mean_surface) 1
set rc [catch {::VMDHole::on_show_mean_surface_toggled} err]
puts "VOLMODE rc=$rc show=$::VMDHole::state(show_mean_surface) status='[string range $::VMDHole::state(status) 0 40]' [expr {$rc==0 && $::VMDHole::state(show_mean_surface)==0 && [string match {Mean surface:*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"

# 2. Volume-trap recovery is reachable: persisted Volume mode, no results
set ::VMDHole::state(mean_3d_mode) "Volume"
set ::VMDHole::state(show_mean_surface) 1
set rc [catch {::VMDHole::on_show_mean_surface_toggled} err]
puts "VOLTRAP rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"

# 3. passability: click with dead molid - no raw error, dialog usable after
set rc1 [catch {::VMDHole::show_passability_dialog} e1]
set d .vmdhole.passability
set st1 [expr {[winfo exists $d] ? [wm state $d] : "none"}]
set rc2 [catch {::VMDHole::show_passability_dialog} e2]
set st2 [expr {[winfo exists $d] ? [wm state $d] : "none"}]
puts "PASSABILITY rc1=$rc1 st1=$st1 rc2=$rc2 st2=$st2 [expr {$rc1==0 && $rc2==0 && $st2 eq "normal" ? {OK} : {BAD}}]"

# 4. ion flow with no molecule: friendly branch, not a raw throw
set ::VMDHole::state(status) ""
set rc [catch {::VMDHole::_run_ion_flow} err]
puts "IONFLOW rc=$rc status='$::VMDHole::state(status)' [expr {$rc==0 && $::VMDHole::state(status) eq {Load a molecule first.} ? {OK} : {BAD}}]"

# 5. permeation dialog opens with no molecule
set rc [catch {::VMDHole::show_permeation_dialog} err]
puts "PERMEATION rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"

# 6. metrics readout with dead molid must not throw (canvas + fake args)
canvas .c
set rc [catch {::VMDHole::_draw_metrics_readout .c 0 0 0 200 100} err]
puts "READOUT rc=$rc [expr {$rc==0 || ![string match {*not available*} $err] ? {OK} : "BAD $err"}]"
S3_EOF
cat > "$T/s4.tcl" <<'S4_EOF'
# 1. gear cross-route write is refused
set ::VMDHole::_gear_open_cid 5
array set ::VMDHole::tunnel_xcid {1,2 7}
proc ::VMDHole::_tunnel_display_frame {} { return 1 }
set ::VMDHole::state(status) ""
catch {::VMDHole::_tunnel_gear_set 2 material Glass} err
set wrote [info exists ::VMDHole::tunnel_gear_cid(7,material)]
puts "GEARGUARD wrote=$wrote status='[string range $::VMDHole::state(status) 0 30]' [expr {!$wrote && [string match {*reopen it*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"
# same-route write still works
set ::VMDHole::_gear_open_cid 7
catch {::VMDHole::_tunnel_gear_set 2 material Glass} err
puts "GEARSAME wrote=[info exists ::VMDHole::tunnel_gear_cid(7,material)] [expr {[info exists ::VMDHole::tunnel_gear_cid(7,material)] ? {OK} : "BAD $err"}]"

# 2. absent-route gear click reports and cleans up
unset -nocomplain ::VMDHole::_gear_open_cid
set ::VMDHole::state(status) ""
catch {::VMDHole::_tunnel_gear_click_cid 99} err
puts "GEARABSENT status='[string range $::VMDHole::state(status) 0 40]' [expr {[string match {*absent from the displayed frame*} $::VMDHole::state(status)] ? {OK} : "BAD $err"}]"

# 3. nan option fields fall back instead of throwing
set ::VMDHole::state(conn_pore_margin) nan
set rc [catch {::VMDHole::_conn_pore_margin} v]
puts "MARGIN rc=$rc v=$v [expr {$rc==0 && $v==2.0 ? {OK} : {BAD}}]"
set ::VMDHole::state(conn_lobe_tolz) nan; set ::VMDHole::state(conn_lobe_tola) inf
set rc [catch {::VMDHole::_conn_lobe_tol} v]
puts "LOBETOL rc=$rc v=$v [expr {$rc==0 && $v eq {6.0 35.0} ? {OK} : {BAD}}]"

# 4. Over Time compute failure is reported
proc ::VMDHole::draw_heatmap {} { error "synthetic draw failure" }
set ::VMDHole::state(status) ""
catch {::VMDHole::on_heatmap_property_compute} err
puts "OTCOMPUTE busy=$::VMDHole::busy status='[string range $::VMDHole::state(status) 0 45]' [expr {$::VMDHole::busy==0 && [string match {Over Time compute failed:*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"

# 5. keep_visualization honours tunnel sessions
set ::VMDHole::state(keep_visualization) 1
set ::VMDHole::tunnel_result_frames {0 1}
set ::VMDHole::state(status) ""
::VMDHole::close_gui
puts "KEEPVIS status='[string range $::VMDHole::state(status) 0 24]' [expr {[string match {Visualization kept*} $::VMDHole::state(status)] ? {OK} : {BAD}}]"

# 6. selection-changed path (incl. lining body follow) runs clean
::VMDHole::show_gui; update
set rc [catch {::VMDHole::_on_tunnel_selection_changed} err]
puts "SELCHANGE rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"
S4_EOF

run_section close-path "$T/s1.tcl"
run_section guards     "$T/s2.tcl"
run_section molid      "$T/s3.tcl"
run_section tunnel-opts "$T/s4.tcl"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
