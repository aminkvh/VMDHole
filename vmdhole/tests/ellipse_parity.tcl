# The ellipse probe's C accelerator must agree with its Tcl reference.
#
# The last of the four engine pairs to get parity coverage. HOLE, MOLE and the
# hydration projection each had a C-vs-Tcl test; the ellipse had both paths
# wired and live (_asym_ellipse_c returns "" -> _asym_ellipse_tcl) with nothing
# asserting they produce the same numbers. test_accel_parity.sh does NOT cover
# this: it compares the accelerated binary against STOCK C, never against Tcl.
#
# Only the first few slices are compared. The Tcl fit is ~60 ms/slice here and
# far worse on a real channel, and what is under test is agreement per slice,
# not throughput - the C path is checked over the whole profile by its own
# length assertion.

set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc note {msg} { puts "  ....  $msg" }
proc done {} {
    global pass fail
    puts "ELLIPSE-RESULT pass=$pass fail=$fail"
    quit
}

set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $root vmdhole 1GRM.pdb]
set RAD  [file join $root native stock_build hole2 rad simple.rad]

cd [file join $here ..]
source vmdhole.tcl
# Point the plugin at the engines the suite built before letting it discover.
# find_hole_exe only searches four fixed install paths, so on a tree that builds
# into native/build this group skipped for "no sos_triangle with --asym-ellipse"
# while a perfectly capable one sat in $VMDHOLE_HOLE_EXE_DIR.
if {[info exists ::env(VMDHOLE_HOLE_EXE_DIR)] && $::env(VMDHOLE_HOLE_EXE_DIR) ne ""} {
    foreach {_k _n} {hole_exec hole sph_process_exec sph_process
                     sos_triangle_exec sos_triangle mole_engine_exec mole_tunnel_engine} {
        set _p [file join $::env(VMDHOLE_HOLE_EXE_DIR) $_n]
        if {[file executable $_p]} { set ::VMDHole::state($_k) $_p }
    }
}
::VMDHole::init_executables

if {![file readable $PDB] || ![file readable $RAD]} { puts "  SKIP  no 1GRM fixture"; done }
if {![::VMDHole::asymmetry_c_available] || ![::VMDHole::sos_triangle_has_feature asymellipse]} {
    puts "  SKIP  no sos_triangle with --asym-ellipse"
    done
}

set work [file join [::VMDHole::get_temp_base] "vmdhole_ell_[pid]"]
file mkdir $work
set mid [mol new $PDB waitfor all]
array set ::VMDHole::state [list molid $mid frame_spec 0 selection all \
    radius_file $RAD cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 8.0 \
    random_seed 1 pore_method circular display_mode none work_dir $work \
    save_results 1 extra_cards {} ignore {}]

if {[catch {::VMDHole::run_analysis} rerr]} {
    chk "the HOLE run needed for a .sph completes" "error: $rerr" "no error"
    catch {file delete -force $work}
    done
}

set g [::VMDHole::_asym_gather $mid 0]
if {$g eq ""} {
    # Not a skip: _asym_gather returning "" is how the missing-.sph bug looked,
    # and quietly skipping here is what let it sit unnoticed.
    chk "the run produced a usable .sph centerline" 0 1
    catch {file delete -force $work}
    done
}
lassign $g centers radii atoms u sph
chk "gathered a centerline" [expr {[llength $centers] > 5 ? 1 : 0}] 1

set N 32
set c [::VMDHole::_asym_ellipse_c $atoms $sph $N]
chk "C returns one row per slice" [llength $c] [llength $centers]

set k 3
set t0 [clock milliseconds]
set t [::VMDHole::_asym_ellipse_tcl $atoms [lrange $centers 0 [expr {$k-1}]] \
           [lrange $radii 0 [expr {$k-1}]] $N]
note "Tcl fitted $k slices in [expr {[clock milliseconds]-$t0}] ms"
chk "Tcl returns one row per requested slice" [llength $t] $k

# The C writes its sidecar at float precision, so agreement is ~1e-6, not exact.
set worst_b 0.0
set worst_a 0.0
for {set i 0} {$i < $k} {incr i} {
    lassign [lindex $c $i] cb ca
    lassign [lindex $t $i] tb ta
    set db [expr {abs($cb - $tb)}]
    set da [expr {abs($ca - $ta)}]
    if {$db > $worst_b} { set worst_b $db }
    if {$da > $worst_a} { set worst_a $da }
}
note "worst |db| [format %.3g $worst_b], worst |da| [format %.3g $worst_a]"
chk "minor axis agrees (|d| < 1e-5)" [expr {$worst_b < 1e-5 ? 1 : 0}] 1
chk "major axis agrees (|d| < 1e-5)" [expr {$worst_a < 1e-5 ? 1 : 0}] 1

catch {mol delete $mid}
catch {file delete -force $work}
done
