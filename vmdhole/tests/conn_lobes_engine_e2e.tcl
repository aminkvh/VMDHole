# End-to-end check of the native conn_lobes engine (classify + cluster +
# triangle-to-region split) against the pure-Tcl reference it replaces, on a
# real Connolly run's own output.
set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc note {msg} { puts "  ....  $msg" }
proc done {} { global pass fail; puts "CL-RESULT pass=$pass fail=$fail"; quit }
set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $root vmdhole 1GRM.pdb]
set RAD  [file join $root native stock_build hole2 rad simple.rad]
cd [file join $here ..]
source vmdhole.tcl
::VMDHole::init_executables
if {[info exists ::env(VMDHOLE_TEST_HOLE)] && $::env(VMDHOLE_TEST_HOLE) ne ""} {
    set ::VMDHole::state(hole_exec) $::env(VMDHOLE_TEST_HOLE)
}
set exe [::VMDHole::tool_path conn_lobes]
if {$exe eq ""} {
    set cand [file join $root native conn_lobes]
    if {[file executable $cand]} { set ::VMDHole::state(conn_lobes_exec) $cand; set exe $cand }
}
if {$exe eq ""} { puts "SKIP: conn_lobes_engine - no conn_lobes binary (sh native/build.sh)"; done }
note "conn_lobes = $exe"
if {![file readable $PDB] || ![file readable $RAD]} { puts "SKIP: conn_lobes_engine - no 1GRM fixture"; done }
if {![file executable $::VMDHole::state(hole_exec)]} { puts "SKIP: conn_lobes_engine - no HOLE binary"; done }

# Produce one real Connolly hole_out.sph to classify against.
set mid [mol new $PDB waitfor all]
set work [file join [::VMDHole::get_temp_base] "vmdhole_cl_[pid]"]
file delete -force $work; file mkdir $work
array set ::VMDHole::state [list molid $mid frame_spec 0 selection all radius_file $RAD \
    cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 12.0 random_seed 1 pore_method connolly \
    display_mode triangulated surface_color pore_lobes dot_density 8 work_dir $work keep_input_pdb 0 \
    save_results 1 extra_cards {} ignore {} hole_fix_atom_names 0 hole_tcl_fallback 0 \
    search_engine mc conn_engine hole prebuild_surfaces 0]
::VMDHole::validate_inputs
if {[catch {::VMDHole::run_analysis} err]} { chk "Connolly run completes" "error: $err" "no error"; done }
chk "one result frame" [llength $::VMDHole::result_frames] 1
if {![llength $::VMDHole::result_frames]} { done }
set f [lindex $::VMDHole::result_frames 0]
set rd [dict get $::VMDHole::results $f run_dir]
set sph [file join $rd hole_out.sph]
chk "hole_out.sph produced" [file exists $sph] 1
if {![file exists $sph]} { done }

set cv {0 0 1}; set cp {0 0 0}; set margin 2.0
set cls_native [::VMDHole::_conn_classify_sph $sph $cv $cp $margin]
chk "native classify ran" [dict exists $cls_native lobes] 1
set lobes_native [::VMDHole::_conn_frame_lobes $cls_native]

rename ::VMDHole::tool_path ::VMDHole::_real_tool_path
proc ::VMDHole::tool_path {name} { if {$name eq "conn_lobes"} { return "" }; ::VMDHole::_real_tool_path $name }
set cls_tcl [::VMDHole::_conn_classify_sph $sph $cv $cp $margin]
set lobes_tcl [::VMDHole::_conn_frame_lobes $cls_tcl]
rename ::VMDHole::tool_path {}
rename ::VMDHole::_real_tool_path ::VMDHole::tool_path

note "n_pore=[dict get $cls_native n_pore] n_lat=[dict get $cls_native n_lat] lobes=[llength $lobes_native]"
chk "pore dot count matches Tcl" [dict get $cls_native n_pore] [dict get $cls_tcl n_pore]
chk "lateral dot count matches Tcl" [dict get $cls_native n_lat] [dict get $cls_tcl n_lat]
chk "pore dot SET matches Tcl" [expr {[lsort [dict get $cls_native pore]] eq [lsort [dict get $cls_tcl pore]]}] 1
chk "lateral dot SET matches Tcl" [expr {[lsort [dict get $cls_native lateral]] eq [lsort [dict get $cls_tcl lateral]]}] 1
chk "keep line SET matches Tcl" [expr {[lsort [dict get $cls_native keep]] eq [lsort [dict get $cls_tcl keep]]}] 1
set _esc_ok 1
foreach rn [dict get $cls_native escaped_ranges] rt [dict get $cls_tcl escaped_ranges] {
    lassign $rn lon hin; lassign $rt lot hit
    if {abs($lon-$lot) > 1e-6 || abs($hin-$hit) > 1e-6} { set _esc_ok 0 }
}
if {[llength [dict get $cls_native escaped_ranges]] != [llength [dict get $cls_tcl escaped_ranges]]} { set _esc_ok 0 }
chk "escaped ranges match Tcl (numerically)" $_esc_ok 1
chk "lobe count matches Tcl" [llength $lobes_native] [llength $lobes_tcl]
set lobe_mismatch 0
foreach ln $lobes_native lt $lobes_tcl {
    lassign $ln zn an nn idxn efn kan kbn
    lassign $lt zt at nt idxt eft kat kbt
    if {$nn != $nt} { incr lobe_mismatch; continue }
    if {abs($zn-$zt) > 1e-6 || abs($an-$at) > 1e-6} { incr lobe_mismatch; continue }
    if {[lsort -integer $idxn] ne [lsort -integer $idxt]} { incr lobe_mismatch }
}
chk "every lobe's stats and membership match Tcl" $lobe_mismatch 0

# HOLE's LAST-REC-END clip markers: both classifiers must report the same set,
# every writer must put them back, and the point-cloud reducer must keep them -
# a clip sphere written without its marker is drawn as a smooth ball.
set src_markers 0
set fh [open $sph r]
while {[gets $fh l] >= 0} { if {[string range $l 0 11] eq "LAST-REC-END"} { incr src_markers } }
close $fh
note "source sph carries $src_markers LAST-REC-END markers"
chk "native and Tcl agree on the marked-record set" \
    [expr {[lsort [dict keys [dict get $cls_native marked]]] eq [lsort [dict keys [dict get $cls_tcl marked]]]}] 1
chk "some records are marked when the source has markers" \
    [expr {$src_markers == 0 || [dict size [dict get $cls_native marked]] > 0}] 1
proc _count_markers {f} { set n 0; set fh [open $f r]; while {[gets $fh l] >= 0} { if {[string range $l 0 11] eq "LAST-REC-END"} { incr n } }; close $fh; return $n }
set _rs [file join $work region_check.sph]
::VMDHole::_write_conn_region_sph $cls_native pore $_rs
set _expect 0
foreach l [concat [dict get $cls_native keep] [dict get $cls_native pore]] { if {[dict exists [dict get $cls_native marked] $l]} { incr _expect } }
chk "region sph restores every marker of its records" [_count_markers $_rs] $_expect
set _red [file join $work reduced_check.sph]
set _nr [::VMDHole::_reduce_conn_sph $sph $_red 0 50]
set _clip 0; set _clipkept 0
set fh [open $sph r]; while {[gets $fh l] >= 0} { if {[string match "ATOM*" $l] && [string trim [string range $l 22 26]] ne "-999"} { incr _clip } }; close $fh
set fh [open $_red r]; while {[gets $fh l] >= 0} { if {[string match "ATOM*" $l] && [string trim [string range $l 22 26]] ne "-999"} { incr _clipkept } }; close $fh
chk "reducer keeps every non-dot (centre/escape/clip) record" $_clipkept $_clip
chk "reducer keeps the markers of what it keeps" [expr {$src_markers == 0 || [_count_markers $_red] > 0}] 1

# region-mesh pipeline: native split vs Tcl split, through the real builder.
proc run_pipeline {label force_tcl} {
    global sph work mid
    if {$force_tcl} {
        rename ::VMDHole::tool_path __saved_engine_path
        proc ::VMDHole::tool_path {name} { if {$name eq "conn_lobes"} { return "" }; __saved_engine_path $name }
    }
    set rundir [file join $work $label]; file delete -force $rundir; file mkdir $rundir
    set cls [::VMDHole::_conn_classify_sph $sph {0 0 1} {0 0 0} 2.0]
    set lobes [::VMDHole::_conn_frame_lobes $cls]
    set regions [list [list pore [dict get $cls pore] red 1]]
    set li 0
    foreach lb $lobes {
        lassign $lb z a n idx ef ka kb
        set lines {}
        foreach i $idx { lappend lines [lindex [dict get $cls lateral] $i] }
        lappend regions [list "lobe$li" $lines blue 1]
        incr li
    }
    set t0 [clock milliseconds]
    set parts [::VMDHole::_build_conn_region_meshes $rundir $cls $regions 8 $sph $mid]
    set ms [expr {[clock milliseconds]-$t0}]
    if {$force_tcl} {
        rename ::VMDHole::tool_path {}
        rename __saved_engine_path ::VMDHole::tool_path
    }
    # parts alternates plot-path, {color material} - key by basename (same
    # tag/dotden naming rule in both runs, only the run dir differs) so the
    # two runs' outputs line up.
    set tris {}
    foreach {plot color} $parts {
        set b [file tail $plot]
        set lines {}
        if {[file exists $plot]} {
            set fh [open $plot r]
            foreach l [split [read $fh] "\n"] { if {[string match "draw trinorm *" $l]} { lappend lines $l } }
            close $fh
        }
        dict set tris $b $lines
    }
    note "$label: [llength $parts] part(s), [expr {$ms}] ms, [dict size $tris] surface(s)"
    return $tris
}
set tris_native [run_pipeline native 0]
set tris_tcl [run_pipeline tclfallback 1]
chk "same surfaces built, native vs Tcl" [lsort [dict keys $tris_native]] [lsort [dict keys $tris_tcl]]
set tri_mismatch 0
dict for {b lines} $tris_native {
    set tl [expr {[dict exists $tris_tcl $b] ? [dict get $tris_tcl $b] : "-missing-"}]
    if {[lsort $lines] ne [lsort $tl]} { incr tri_mismatch; note "  triangle set differs for $b: native=[llength $lines] tcl=[llength $tl]" }
}
chk "every surface's triangle set is BYTE-IDENTICAL to Tcl" $tri_mismatch 0

catch {file delete -force $work}
done
