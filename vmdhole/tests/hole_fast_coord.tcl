# The packed coordinate record must be an ACCELERATOR, not a second answer.
#
# HOLE is an external Fortran program: its only input path is a named coordinate
# file, so every frame gets one written for it. Writing that as a PDB costs
# 12.1 ms/frame of ASCII formatting at 18.7k atoms and runs SERIALLY while HOLE
# itself is spread over every worker, which is most of the parallel ceiling for
# the fast pore methods. The patched reader (tsatr_fast.f) takes a packed record
# instead.
#
# What is actually at risk, and therefore what this checks:
#   * The PDB path quantises coordinates to 0.001 A by writing %8.3f and reading
#     3F8.3. The packed path carries full doubles, so tsatr_fast.f rounds on the
#     way in. Without that the sphere search diverges and the .sph differs - and
#     it would look like a mysterious Monte-Carlo wobble, not a bug.
#   * Anything that cannot take a packed record - stock HOLE, the pure-Tcl
#     engine, a user who ticked "Keep input_frame.pdb" - must still get a PDB.
set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc done {} { global pass fail; puts "FASTCOORD-RESULT pass=$pass fail=$fail"; quit }

set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $root vmdhole 1GRM.pdb]
set RAD  [file join $root native stock_build hole2 rad simple.rad]
cd [file join $here ..]
source vmdhole.tcl
::VMDHole::init_executables
# Lets the suite aim at a freshly built HOLE without touching the user's
# configured one - the byte-identity half is meaningless against a binary that
# cannot read the record.
if {[info exists ::env(VMDHOLE_TEST_HOLE)] && $::env(VMDHOLE_TEST_HOLE) ne ""} {
    set ::VMDHole::state(hole_exec) $::env(VMDHOLE_TEST_HOLE)
}

# --- gates that must hold whatever binary is configured ----------------------
set _sv [expr {[info exists ::VMDHole::state(keep_input_pdb)] ? $::VMDHole::state(keep_input_pdb) : 0}]
set ::VMDHole::state(keep_input_pdb) 1
chk "'Keep input_frame.pdb' forces the PDB path" [::VMDHole::_hole_fast_coord_available] 0
set ::VMDHole::state(keep_input_pdb) $_sv
set _sv_exe $::VMDHole::state(hole_exec)
set ::VMDHole::state(hole_exec) ""
chk "no binary at all means no packed record" [::VMDHole::_hole_fast_coord_available] 0
set ::VMDHole::state(hole_exec) $_sv_exe

if {![file readable $PDB] || ![file readable $RAD]} {
    puts "  ....  no 1GRM fixture - gate checks only"
    done
}
if {![::VMDHole::_hole_fast_coord_available]} {
    # Not a SKIP: the gate returning 0 for a HOLE that cannot read the record is
    # the behaviour stock users depend on, and asserting it is a real check.
    puts "  ....  configured HOLE has no fast-coord-read patch - gate checks only"
    chk "a HOLE without the patch is refused the packed record" \
        [::VMDHole::_hole_fast_coord_available] 0
    done
}

# --- the real thing: same frame, both input paths, identical .sph ------------
set work [file join [::VMDHole::get_temp_base] "vmdhole_fc_[pid]"]
file delete -force $work; file mkdir $work
set mid [mol new $PDB waitfor all]
array set ::VMDHole::state [list molid $mid frame_spec 0 selection all \
    radius_file $RAD cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 8.0 \
    random_seed 1 pore_method circular display_mode none work_dir $work \
    save_results 1 extra_cards {} ignore {} hole_fix_atom_names 0]
set sel [atomselect $mid all]
$sel frame 0
set out {}
foreach mode {pdb vhb} {
    set d [file join $work $mode]
    file mkdir $d
    if {$mode eq "pdb"} {
        set cn input_frame.pdb
        $sel writepdb [file join $d input_frame.pdb]
    } else {
        set pr [::VMDHole::_hole_coord_identity $sel [file join $d _id.pdb]]
        if {[llength $pr] != 2} {
            chk "the fixture can use the packed path" 0 1
            $sel delete; catch {file delete -force $work}; done
        }
        catch {file delete [file join $d _id.pdb]}
        lassign $pr n idb
        set cn input_frame.vhb
        ::VMDHole::_write_hole_coord_bin $sel [file join $d input_frame.vhb] $n $idb
    }
    ::VMDHole::write_control_file [file join $d hole.inp] $cn hole_out.sph {0 0 0} {0 0 1}
    catch {exec sh -c "cd [::VMDHole::shell_quote $d] && [::VMDHole::shell_quote $::VMDHole::state(hole_exec)] < hole.inp > hole_out.txt 2>&1"}
    set c ""
    if {[file exists [file join $d hole_out.sph]]} {
        set f [open [file join $d hole_out.sph] rb]; set c [read $f]; close $f
    }
    lappend out $c
}
lassign $out A B
chk "the PDB path produced a .sph" [expr {[string length $A] > 0}] 1
chk "the packed path produced a .sph" [expr {[string length $B] > 0}] 1
chk "both inputs give a BYTE-IDENTICAL .sph" [expr {$A ne "" && $A eq $B}] 1
$sel delete
catch {file delete -force $work}
done
