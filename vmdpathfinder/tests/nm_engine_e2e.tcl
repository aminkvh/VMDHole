# End-to-end check of the Nelder-Mead search engine and the fast Connolly
# surface through the plugin's own run path (validate_inputs + run_analysis),
# against HOLE's Monte Carlo search and HOLE's own conn on the same input.
set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc note {msg} { puts "  ....  $msg" }
proc done {} { global pass fail; puts "NM-RESULT pass=$pass fail=$fail"; quit }
set here [expr {[info exists ::env(VMDPATHFINDER_TEST_DIR)] && $::env(VMDPATHFINDER_TEST_DIR) ne ""
                ? $::env(VMDPATHFINDER_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $root vmdpathfinder 1GRM.pdb]
set RAD  [file join $root native stock_build hole2 rad simple.rad]
cd [file join $here ..]
source vmdpathfinder.tcl
::VMDPathFinder::init_executables
if {[info exists ::env(VMDPATHFINDER_TEST_HOLE)] && $::env(VMDPATHFINDER_TEST_HOLE) ne ""} {
    set ::VMDPathFinder::state(hole_exec) $::env(VMDPATHFINDER_TEST_HOLE)
}

# --- engine selection logic, no binaries needed
foreach {pm se ce fb want} {
    circular mc hole 0 {}
    circular nm hole 0 nm
    circular nm fast 0 nm
    connolly mc hole 0 {}
    connolly mc fast 0 mcconn
    connolly nm fast 0 nmconn
    connolly nm hole 0 nmconn
    capsule  nm fast 0 {}
    circular nm fast 1 {}
} {
    set ::VMDPathFinder::state(pore_method) $pm
    set ::VMDPathFinder::state(search_engine) $se
    set ::VMDPathFinder::state(conn_engine) $ce
    chk "mode for $pm/$se/$ce/fallback=$fb" [::VMDPathFinder::_nm_mode_for_run $fb] $want
}

# --- control-file conn stripping for the MC + fast Connolly route
set tmpd [file join [::VMDPathFinder::get_temp_base] "vmdpathfinder_nm_[pid]"]
file delete -force $tmpd; file mkdir $tmpd
array set ::VMDPathFinder::state [list pore_method connolly extra_cards {conn 1.15 0.5; cutsize 0.9} \
    cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 8.0 radius_file $RAD ignore {}]
::VMDPathFinder::write_control_file [file join $tmpd a.inp] input_frame.pdb hole_out.sph {0 0 0} {0 0 1}
::VMDPathFinder::write_control_file [file join $tmpd b.inp] input_frame.pdb hole_out.sph {0 0 0} {0 0 1} 1
proc _cards {f} { set fh [open $f r]; set t [read $fh]; close $fh; set c {}; foreach l [split $t \n] { lappend c [string tolower [lindex $l 0]] }; return $c }
chk "full control file carries the conn card"  [expr {"conn" in [_cards [file join $tmpd a.inp]]}] 1
chk "stripped control file has no conn card"   [expr {"conn" in [_cards [file join $tmpd b.inp]]}] 0
chk "stripped control file keeps other cards"  [expr {"cutsize" in [_cards [file join $tmpd b.inp]]}] 1
set ::VMDPathFinder::state(extra_cards) {}

# --- engine + fixtures
set exe [::VMDPathFinder::tool_path nm_search]
if {$exe eq ""} {
    set cand [file join $root native nm nm_search]
    if {[file executable $cand]} { set ::VMDPathFinder::state(nm_search_exec) $cand; set exe $cand }
}
if {$exe eq ""} { puts "SKIP: nm_engine - no nm_search binary (sh native/build.sh)"; done }
note "nm_search = $exe"
if {![file readable $PDB] || ![file readable $RAD]} { puts "SKIP: nm_engine - no 1GRM fixture"; done }
if {![file executable $::VMDPathFinder::state(hole_exec)]} { puts "SKIP: nm_engine - no HOLE binary to compare against"; done }

set mid [mol new $PDB waitfor all]
set ncase 0
proc run_case {label settings} {
    global mid tmpd ncase
    incr ncase
    set work [file join $tmpd "case$ncase"]
    file mkdir $work
    catch {::VMDPathFinder::clear_results_for_new_settings}
    array set ::VMDPathFinder::state [list \
        molid $mid frame_spec 0 selection all \
        cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 8.0 \
        random_seed 1 display_mode none work_dir $work keep_input_pdb 0 \
        save_results 1 extra_cards {} ignore {} hole_fix_atom_names 0 hole_tcl_fallback 0]
    array set ::VMDPathFinder::state $settings
    ::VMDPathFinder::validate_inputs
    set t0 [clock milliseconds]
    if {[catch {::VMDPathFinder::run_analysis} err]} {
        chk "$label: run_analysis completes" "error: $err" "no error"
        return {}
    }
    note "$label: run took [expr {[clock milliseconds]-$t0}] ms"
    chk "$label: one result frame" [llength $::VMDPathFinder::result_frames] 1
    if {![llength $::VMDPathFinder::result_frames]} { return {} }
    set f [lindex $::VMDPathFinder::result_frames 0]
    set p [dict get $::VMDPathFinder::results $f profile]
    dict set p run_dir [dict get $::VMDPathFinder::results $f run_dir]
    chk "$label: profile valid" [dict get $p valid] 1
    note "$label: points=[dict get $p points] min_radius=[dict get $p min_radius] at [dict get $p min_coord]"
    return $p
}
proc sph_stats {rd} {
    set f [file join $rd hole_out.sph]
    if {![file exists $f]} { return [list 0 0] }
    set fh [open $f r]; set nc 0; set nd 0
    while {[gets $fh l] >= 0} {
        if {![string match "ATOM*" $l]} continue
        set r [string trim [string range $l 22 25]]
        if {$r eq "-999"} { incr nd } elseif {$r ne "-888"} { incr nc }
    }
    close $fh
    return [list $nc $nd]
}
proc tsv_min_requiv {rd} {
    set f [file join $rd hole_profile.tsv]
    if {![file exists $f]} { return "" }
    set fh [open $f r]; gets $fh; set best ""
    while {[gets $fh l] >= 0} {
        set c [lindex [split $l \t] 4]
        if {[string is double -strict $c] && $c < 1e5 && ($best eq "" || $c < $best)} { set best $c }
    }
    close $fh
    return $best
}

# packed record vs PDB: same engine, same frame, must be the same sphere file
set sel [atomselect $mid all]
set d [file join $tmpd io]; file mkdir $d
$sel writepdb [file join $d f.pdb]
lassign [::VMDPathFinder::_hole_coord_identity $sel [file join $d _id.pdb]] _n _idb
::VMDPathFinder::_write_hole_coord_bin $sel [file join $d f.vhb] $_n $_idb
$sel delete
foreach k {pdb vhb} {
    catch {exec $exe {*}[::VMDPathFinder::tool_args nm_search] [file join $d f.$k] $RAD 0 0 0 0 0 1 0.5 8.0 --sph [file join $d $k.sph] --quiet} _
}
proc _slurp {f} { if {![file exists $f]} { return "" }; set fh [open $f rb]; set c [read $fh]; close $fh; return $c }
chk "engine reads the packed coordinate record" [expr {[string length [_slurp [file join $d vhb.sph]]] > 0}] 1
chk "packed record and PDB give a BYTE-IDENTICAL sphere file" [expr {[_slurp [file join $d pdb.sph]] ne "" && [_slurp [file join $d pdb.sph]] eq [_slurp [file join $d vhb.sph]]}] 1

set A [run_case "MC spherical" {pore_method circular search_engine mc conn_engine hole}]
set B [run_case "NM spherical" {pore_method circular search_engine nm conn_engine hole}]
if {[dict size $A] && [dict size $B]} {
    set rd [dict get $B run_dir]
    chk "NM: engine log kept in the run folder" [file exists [file join $rd nm_search.log]] 1
    chk "NM: no HOLE log (profile comes from the TSV)" [file exists [file join $rd hole_out.txt]] 0
    lassign [sph_stats $rd] nc nd
    chk "NM: sphere file has centreline records" [expr {$nc > 10}] 1
    set d [expr {abs([dict get $A min_radius] - [dict get $B min_radius])}]
    note "bottleneck MC [dict get $A min_radius] vs NM [dict get $B min_radius] (|d| = [format %.3f $d])"
    chk "NM bottleneck within 0.25 A of Monte Carlo" [expr {$d < 0.25}] 1
    set sig [::VMDPathFinder::run_signature $mid all]
    chk "run signature records the search engine" [string match "*search|nm|*" $sig] 1
}
set C [run_case "HOLE Connolly" {pore_method connolly search_engine mc conn_engine hole}]
set D [run_case "MC + fast Connolly" {pore_method connolly search_engine mc conn_engine fast}]
if {[dict size $C] && [dict size $D]} {
    lassign [sph_stats [dict get $C run_dir]] ncC ndC
    lassign [sph_stats [dict get $D run_dir]] ncD ndD
    note "Connolly dots: HOLE $ndC, fast $ndD (centres $ncC / $ncD)"
    chk "fast Connolly kept HOLE's spherical log" [file exists [file join [dict get $D run_dir] hole_sph_out.txt]] 1
    chk "fast Connolly produced dots" [expr {$ndD > 0}] 1
    chk "fast Connolly dot count within 2% of HOLE conn" [expr {$ndC > 0 && abs($ndD - $ndC) <= 0.02 * $ndC + 2}] 1
    set rqC [tsv_min_requiv [dict get $C run_dir]]
    set rqD [tsv_min_requiv [dict get $D run_dir]]
    note "min Requiv: HOLE $rqC, fast $rqD"
    chk "fast Connolly Requiv within 0.05 A of HOLE" [expr {$rqC ne "" && $rqD ne "" && abs($rqC - $rqD) < 0.05}] 1
    chk "MC bottleneck unchanged by the fast Connolly" \
        [expr {abs([dict get $C min_radius] - [dict get $D min_radius]) < 0.011}] 1
}
set E [run_case "NM + fast Connolly" {pore_method connolly search_engine nm conn_engine fast}]
if {[dict size $E]} {
    lassign [sph_stats [dict get $E run_dir]] ncE ndE
    chk "NM + fast Connolly produced dots" [expr {$ndE > 0}] 1
    set rqE [tsv_min_requiv [dict get $E run_dir]]
    note "NM + fast Connolly: dots $ndE, min Requiv $rqE"
    chk "NM + fast Connolly has Requiv columns" [expr {$rqE ne ""}] 1
    if {[dict size $C]} {
        set d [expr {abs([dict get $C min_radius] - [dict get $E min_radius])}]
        chk "NM bottleneck under Connolly within 0.25 A of HOLE" [expr {$d < 0.25}] 1
    }
}
# --- --neck: the widest route through a wall with a round hole
# Carbon atoms on a 1 A lattice in the plane y=0, a hole of radius 5 at the
# origin. The route from y=-6 to the dots at y=+6 must pass the hole, so the
# neck is the hole radius less the atom radius (1.85 in simple.rad). A second
# block with no hole has no route at all.
set nmx [::VMDPathFinder::tool_path nm_search]
if {$nmx ne "" && [file readable $RAD]} {
    set wall [file join $tmpd wall.pdb]
    set fh [open $wall w]
    set n 0
    for {set x -12} {$x <= 12} {incr x} {
        for {set z -12} {$z <= 12} {incr z} {
            if {$x*$x + $z*$z < 25} continue
            incr n
            puts $fh [format "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C" $n $n $x 0.0 $z]
        }
    }
    close $fh
    set solid [file join $tmpd solid.pdb]
    set fh [open $solid w]
    set n 0
    for {set x -12} {$x <= 12} {incr x} {
        for {set z -12} {$z <= 12} {incr z} {
            incr n
            puts $fh [format "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C" $n $n $x 0.0 $z]
        }
    }
    close $fh
    set lobe [file join $tmpd lobe.txt]
    set fh [open $lobe w]
    puts $fh "LOBE 0 -6 0 3"
    puts $fh "0 6 0"; puts $fh "1 6 0"; puts $fh "0 6 1"
    close $fh
    set got ""
    catch {set got [exec $nmx {*}[::VMDPathFinder::tool_args nm_search] $wall $RAD 0 0 0 0 0 1 0.25 22 --neck $lobe --quiet]}
    note "--neck through the hole: $got"
    chk "--neck reports the hole's clearance (5 - 1.85, within 0.35 A)" \
        [expr {[string is double -strict [lindex $got 0]] && abs([lindex $got 0] - 3.15) < 0.35}] 1
    set got ""
    catch {set got [exec $nmx {*}[::VMDPathFinder::tool_args nm_search] $solid $RAD 0 0 0 0 0 1 0.25 22 --neck $lobe --quiet]}
    chk "--neck finds no route through a solid wall" $got "-"
} else {
    note "SKIP --neck: no nm_search binary or radius file"
}

catch {file delete -force $tmpd}
done
