#!/bin/sh
# GUI smoke + regression suite, driven WITHOUT VMD: the real vmdpathfinder.tcl is
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
export VMDPATHFINDER_ROOT="$ROOT"

echo "gui-smoke: $ROOT/vmdpathfinder/vmdpathfinder.tcl"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
if ! echo 'if {[catch {package require Tk}]} {exit 1}; exit 0' | tclsh >/dev/null 2>&1; then
    echo "SKIP: no usable Tk/DISPLAY (try xvfb-run)"; exit 0
fi

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM
export VMDPATHFINDER_HARNESS_TMP="$T"
cat > "$T/harness.tcl" <<'HARNESS_EOF'
# Shared GUI-review harness: source vmdpathfinder.tcl under plain tclsh+Tk with VMD
# stubs and build the real GUI. Usage:  tclsh harness.tcl <script.tcl>
# The child script runs AFTER show_gui with the .vmdpathfinder toplevel built.
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
namespace eval ::VMDPathFinder {}
source [file join $::env(VMDPATHFINDER_ROOT) vmdpathfinder vmdpathfinder.tcl]
# keep config I/O away from the user's real ~/.vmdpathfinder_config
set ::HARNESS_TMP [expr {[info exists ::env(VMDPATHFINDER_HARNESS_TMP)] ? $::env(VMDPATHFINDER_HARNESS_TMP) : [pwd]}]
set ::VMDPathFinder::config_file [file join $::HARNESS_TMP .harness_vmdpathfinder_config]
catch {file delete $::VMDPathFinder::config_file}

::VMDPathFinder::show_gui
update

# ---- helpers for review scripts --------------------------------------------
proc walk {w} { set out [list $w]; foreach c [winfo children $w] { lappend out {*}[walk $c] }; return $out }
proc all_widgets {} { walk .vmdpathfinder }
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
set t0 [llength [trace info variable ::VMDPathFinder::state(display_mode)]]
::VMDPathFinder::close_gui
::VMDPathFinder::show_gui
update
set t1 [llength [trace info variable ::VMDPathFinder::state(display_mode)]]
set t2 [llength [trace info variable ::VMDPathFinder::state(hydro_scheme)]]
puts "TRACES before=$t0 after_reopen=$t1 hydro=$t2 [expr {$t1>0 && $t2>0 ? {OK} : {BAD}}]"

# 2. abort survives close only while work is in flight
set ::VMDPathFinder::busy 1; set ::VMDPathFinder::_calc_depth 1
::VMDPathFinder::close_gui
puts "ABORT inflight: raised=[set ::VMDPathFinder::state(abort_requested)] [expr {$::VMDPathFinder::state(abort_requested)==1 ? {OK} : {BAD}}]"
set ::VMDPathFinder::busy 0
::VMDPathFinder::_end_calc
puts "ABORT after unwind: [set ::VMDPathFinder::state(abort_requested)] depth=$::VMDPathFinder::_calc_depth [expr {$::VMDPathFinder::state(abort_requested)==0 && $::VMDPathFinder::_calc_depth==0 ? {OK} : {BAD}}]"
::VMDPathFinder::show_gui; update
::VMDPathFinder::close_gui
puts "ABORT idle close: [set ::VMDPathFinder::state(abort_requested)] [expr {$::VMDPathFinder::state(abort_requested)==0 ? {OK} : {BAD}}]"

# 3. a pending modal answer is forced to Cancel by close
::VMDPathFinder::show_gui; update
unset -nocomplain ::VMDPathFinder::_overwrite_ans
after 200 ::VMDPathFinder::close_gui
set ::VMDPathFinder::_overwrite_ans_pending 1
# simulate the dialog's wait: tkwait returns only if close_gui answers
after 2000 {set ::VMDPathFinder::_overwrite_ans timeout}
tkwait variable ::VMDPathFinder::_overwrite_ans
puts "OVERWRITE answered=[set ::VMDPathFinder::_overwrite_ans] [expr {$::VMDPathFinder::_overwrite_ans eq 0 ? {OK} : {BAD}}]"
S1_EOF
cat > "$T/s2.tcl" <<'S2_EOF'
# depth>0 (bracketed pass) must refuse Run / Reset / Align / molid swap
set ::VMDPathFinder::busy 0
::VMDPathFinder::_begin_calc
set r [::VMDPathFinder::run_analysis]
puts "RUN mid-pass: rc=$r status='$::VMDPathFinder::state(status)' [expr {$r==0 && [string match {*in progress*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"
set ::VMDPathFinder::state(status) ""
::VMDPathFinder::reset_session
puts "RESET mid-pass: [expr {[string match {*in progress*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"
set ::VMDPathFinder::state(status) ""
::VMDPathFinder::do_align_trajectory
puts "ALIGN mid-pass: [expr {[string match {*in progress*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"
::VMDPathFinder::_end_calc

# align brackets busy/_begin_calc and restores on throw
set ::seen_busy -1; set ::seen_depth -1
proc ::VMDPathFinder::align_trajectory {{molid ""}} {
    set ::seen_busy $::VMDPathFinder::busy
    set ::seen_depth $::VMDPathFinder::_calc_depth
    error "boom"
}
set ::TKDIALOG_LOG {}
::VMDPathFinder::do_align_trajectory
puts "ALIGN bracket: during busy=$::seen_busy depth=$::seen_depth after busy=$::VMDPathFinder::busy depth=$::VMDPathFinder::_calc_depth dialogs=[llength $::TKDIALOG_LOG] [expr {$::seen_busy==1 && $::seen_depth==1 && $::VMDPathFinder::busy==0 && $::VMDPathFinder::_calc_depth==0 && [llength $::TKDIALOG_LOG]==1 ? {OK} : {BAD}}]"

# run_current_mode must not double-release an enclosing bracket
proc ::VMDPathFinder::analysis_mode {} { return tunnel }
proc ::VMDPathFinder::run_tunnel_analysis {} {
    # model the real wrapper: restore, then rethrow
    set ::VMDPathFinder::busy 0
    ::VMDPathFinder::_end_calc
    error "engine died"
}
::VMDPathFinder::_begin_calc            ;# outer bracket (depth 1)
set ::VMDPathFinder::busy 1
::VMDPathFinder::_begin_calc            ;# the run's own bracket (depth 2)
::VMDPathFinder::run_current_mode
puts "BACKSTOP: depth=$::VMDPathFinder::_calc_depth (outer bracket must survive) [expr {$::VMDPathFinder::_calc_depth==1 ? {OK} : {BAD}}]"
::VMDPathFinder::_end_calc
S2_EOF
cat > "$T/s3.tcl" <<'S3_EOF'
# state: results empty, no molecules (STUB_MOLS empty) - the delete-and-move-on world
set ::VMDPathFinder::state(molid) 0

# 1. _volmode: ticking Show 3D with no results must not raise, must untick + report
set ::VMDPathFinder::state(show_mean_surface) 1
set rc [catch {::VMDPathFinder::on_show_mean_surface_toggled} err]
puts "VOLMODE rc=$rc show=$::VMDPathFinder::state(show_mean_surface) status='[string range $::VMDPathFinder::state(status) 0 40]' [expr {$rc==0 && $::VMDPathFinder::state(show_mean_surface)==0 && [string match {Mean surface:*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"

# 2. Volume-trap recovery is reachable: persisted Volume mode, no results
set ::VMDPathFinder::state(mean_3d_mode) "Volume"
set ::VMDPathFinder::state(show_mean_surface) 1
set rc [catch {::VMDPathFinder::on_show_mean_surface_toggled} err]
puts "VOLTRAP rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"

# 3. passability: click with dead molid - no raw error, dialog usable after
set rc1 [catch {::VMDPathFinder::show_passability_dialog} e1]
set d .vmdpathfinder.passability
set st1 [expr {[winfo exists $d] ? [wm state $d] : "none"}]
set rc2 [catch {::VMDPathFinder::show_passability_dialog} e2]
set st2 [expr {[winfo exists $d] ? [wm state $d] : "none"}]
puts "PASSABILITY rc1=$rc1 st1=$st1 rc2=$rc2 st2=$st2 [expr {$rc1==0 && $rc2==0 && $st2 eq "normal" ? {OK} : {BAD}}]"

# 4. ion flow with no molecule: friendly branch, not a raw throw
set ::VMDPathFinder::state(status) ""
set rc [catch {::VMDPathFinder::_run_ion_flow} err]
puts "IONFLOW rc=$rc status='$::VMDPathFinder::state(status)' [expr {$rc==0 && $::VMDPathFinder::state(status) eq {Load a molecule first.} ? {OK} : {BAD}}]"

# 5. permeation dialog opens with no molecule
set rc [catch {::VMDPathFinder::show_permeation_dialog} err]
puts "PERMEATION rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"

# 6. metrics readout with dead molid must not throw (canvas + fake args)
canvas .c
set rc [catch {::VMDPathFinder::_draw_metrics_readout .c 0 0 0 200 100} err]
puts "READOUT rc=$rc [expr {$rc==0 || ![string match {*not available*} $err] ? {OK} : "BAD $err"}]"
S3_EOF
cat > "$T/s4.tcl" <<'S4_EOF'
# 1. gear cross-route write is refused
set ::VMDPathFinder::_gear_open_cid 5
array set ::VMDPathFinder::tunnel_xcid {1,2 7}
proc ::VMDPathFinder::_tunnel_display_frame {} { return 1 }
set ::VMDPathFinder::state(status) ""
catch {::VMDPathFinder::_tunnel_gear_set_from_popup 2 material Glass} err
set wrote [info exists ::VMDPathFinder::tunnel_gear_cid(7,material)]
puts "GEARGUARD wrote=$wrote status='[string range $::VMDPathFinder::state(status) 0 30]' [expr {!$wrote && [string match {*reopen it*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"
# same-route write still works
set ::VMDPathFinder::_gear_open_cid 7
catch {::VMDPathFinder::_tunnel_gear_set_from_popup 2 material Glass} err
puts "GEARSAME wrote=[info exists ::VMDPathFinder::tunnel_gear_cid(7,material)] [expr {[info exists ::VMDPathFinder::tunnel_gear_cid(7,material)] ? {OK} : "BAD $err"}]"

# 2. absent-route gear click reports and cleans up
unset -nocomplain ::VMDPathFinder::_gear_open_cid
set ::VMDPathFinder::state(status) ""
catch {::VMDPathFinder::_tunnel_gear_click_cid 99} err
puts "GEARABSENT status='[string range $::VMDPathFinder::state(status) 0 40]' [expr {[string match {*absent from the displayed frame*} $::VMDPathFinder::state(status)] ? {OK} : "BAD $err"}]"

# 3. nan option fields fall back instead of throwing
set ::VMDPathFinder::state(conn_pore_margin) nan
set rc [catch {::VMDPathFinder::_conn_pore_margin} v]
puts "MARGIN rc=$rc v=$v [expr {$rc==0 && $v==2.0 ? {OK} : {BAD}}]"
set ::VMDPathFinder::state(conn_lobe_tolz) nan; set ::VMDPathFinder::state(conn_lobe_tola) inf
set rc [catch {::VMDPathFinder::_conn_lobe_tol} v]
puts "LOBETOL rc=$rc v=$v [expr {$rc==0 && $v eq {6.0 35.0} ? {OK} : {BAD}}]"

# 4. Over Time compute failure is reported
proc ::VMDPathFinder::draw_heatmap {} { error "synthetic draw failure" }
set ::VMDPathFinder::state(status) ""
catch {::VMDPathFinder::on_heatmap_property_compute} err
puts "OTCOMPUTE busy=$::VMDPathFinder::busy status='[string range $::VMDPathFinder::state(status) 0 45]' [expr {$::VMDPathFinder::busy==0 && [string match {Over Time compute failed:*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"

# 5. keep_visualization honours tunnel sessions
set ::VMDPathFinder::state(keep_visualization) 1
set ::VMDPathFinder::tunnel_result_frames {0 1}
set ::VMDPathFinder::state(status) ""
::VMDPathFinder::close_gui
puts "KEEPVIS status='[string range $::VMDPathFinder::state(status) 0 24]' [expr {[string match {Visualization kept*} $::VMDPathFinder::state(status)] ? {OK} : {BAD}}]"

# 6. selection-changed path (incl. lining body follow) runs clean
::VMDPathFinder::show_gui; update
set rc [catch {::VMDPathFinder::_on_tunnel_selection_changed} err]
puts "SELCHANGE rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"
S4_EOF

cat > "$T/s5.tcl" <<'S5_EOF'
# 1. inf/nan in _num_or fields fall back
set ::VMDPathFinder::state(mole_bottleneck) inf
puts "NUMOR [::VMDPathFinder::_num_or mole_bottleneck 3.0 1] [expr {[::VMDPathFinder::_num_or mole_bottleneck 3.0 1]==3.0 ? {OK} : {BAD}}]"
# 2. water probe survives an emptied/mistyped field
set ::VMDPathFinder::state(metrics_water_probe) ""
set a [::VMDPathFinder::_metrics_water_probe]
set ::VMDPathFinder::state(metrics_water_probe) "1,15"
set b [::VMDPathFinder::_metrics_water_probe]
puts "WPROBE '$a' '$b' [expr {$a==1.15 && $b==1.15 ? {OK} : {BAD}}]"
# 3. show_gui must not clear a live run's abort
::VMDPathFinder::_begin_calc
set ::VMDPathFinder::state(abort_requested) 1
::VMDPathFinder::show_gui
puts "SHOWGUI-LIVE abort=$::VMDPathFinder::state(abort_requested) [expr {$::VMDPathFinder::state(abort_requested)==1 ? {OK} : {BAD}}]"
::VMDPathFinder::_end_calc
# 4. save publishes by rename; no .tmp litter
set cfg $::VMDPathFinder::config_file
::VMDPathFinder::save_config
set tmps [glob -nocomplain "$cfg.tmp*"]
puts "SAVECFG exists=[file exists $cfg] tmps=[llength $tmps] [expr {[file exists $cfg] && ![llength $tmps] ? {OK} : {BAD}}]"
# 5. reopen must not revert a saved rename-atoms setting
set fh [open $cfg w]; puts $fh "hole_fix_atom_names = 1"; puts $fh "ion_radius_fallback = 0"; close $fh
::VMDPathFinder::load_config
puts "MIRROR fix=$::VMDPathFinder::state(hole_fix_atom_names) [expr {$::VMDPathFinder::state(hole_fix_atom_names)==1 ? {OK} : {BAD}}]"
# 6. torn boolean must not blow up close_gui
set fh [open $cfg w]; puts $fh "keep_visualization ="; close $fh
::VMDPathFinder::load_config
set rc [catch {::VMDPathFinder::close_gui} err]
puts "TORNBOOL rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"
# 7. unreadable config: GUI still opens
set ::VMDPathFinder::config_file [file join $::HARNESS_TMP .harness_cfg_dir]
file mkdir $::VMDPathFinder::config_file
set rc [catch {::VMDPathFinder::load_config} err]
puts "UNREADABLE rc=$rc [expr {$rc==0 ? {OK} : "BAD $err"}]"
set ::VMDPathFinder::config_file $cfg
# 8. the weight the combobox writes is what save persists
set ::VMDPathFinder::state(mole_weight_disp) "Length"
::VMDPathFinder::save_config
set fh [open $cfg r]; set t [read $fh]; close $fh
puts "WEIGHTKEY [expr {[string match "*mole_weight_disp = Length*" $t] ? {OK} : {BAD}}]"
# 9. a caller-preset engine path is not persisted by init_executables
set eng [file join $::HARNESS_TMP .harness_fake_hole]
set fh [open $eng w]; puts $fh "#!/bin/sh"; close $fh
file attributes $eng -permissions 0755
file delete -force $cfg
set ::VMDPathFinder::state(hole_exec) $eng
::VMDPathFinder::init_executables
set saved [expr {[file exists $cfg] ? [string match "*$eng*" [read [open $cfg r]]] : 0}]
puts "PRESETSAVE saved=$saved [expr {!$saved ? {OK} : {BAD}}]"
S5_EOF

run_section close-path "$T/s1.tcl"
run_section guards     "$T/s2.tcl"
run_section molid      "$T/s3.tcl"
run_section tunnel-opts "$T/s4.tcl"
run_section config-fields "$T/s5.tcl"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
