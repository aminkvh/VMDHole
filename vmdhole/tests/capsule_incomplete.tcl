# HOLE's own "RUN MAY BE INCOMPLETE" warning must reach the user.
#
# HOLE can stop its axial search early, say so in its output, and still exit 0.
# vmdhole_status.dat records the ENGINE's exit status, so it cannot see this:
# the frame becomes a normal result and nothing on screen says the profile
# stops short of the axis the user set. That is what "the capsule lining is far
# away from the pore" turned out to be - the lining covers the profile, and the
# profile stopped early.
#
# The two literals are HOLE's own, taken from the binary's strings, and there
# are exactly two:
#     WARNING RUN MAY BE INCOMPLETE IN +VE DIRECTION
#     WARNING RUN MAY BE INCOMPLETE IN -VE DIRECTION
#
# Why THIS fixture: the warning is structure-specific, not a knob. Measured over
# {structure} x {radius file} x {sample} x {endrad 5..100}, only this monomer
# WITH amberuni.rad reproduces it - the repo's own 1GRM.pdb does not, at any
# endrad, under either method, and neither does its first NMR model despite
# having the same 272 atoms. A high endrad does NOT force it either. So the case
# is frozen here rather than synthesised, and it repeats 5 runs out of 5.
#
# The control is the point of the test: the SAME structure under the SAME cards
# with only the method changed must produce NO warning. Without it this would
# pass just as well against a detector that fires on every run.

set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc note {msg} { puts "  ....  $msg" }

set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $here fixtures capsule_incomplete_1grm.pdb]
set RAD  [file join $root native stock_build hole2 rad amberuni.rad]
# Honour VMDHOLE_HOLE_EXE_DIR like the rest of the suite; ~/hole2/exe stays the
# install default. Hardcoding $HOME skipped this group on any tree whose HOLE
# lives elsewhere, and run_tests.sh then counted it as a group that passed.
set _exedir [expr {[info exists ::env(VMDHOLE_HOLE_EXE_DIR)] && $::env(VMDHOLE_HOLE_EXE_DIR) ne ""
                   ? $::env(VMDHOLE_HOLE_EXE_DIR) : [file join $::env(HOME) hole2 exe]}]
set REF  [file join $_exedir hole]

cd [file join $here ..]
source vmdhole.tcl

if {![file readable $PDB] || ![file readable $RAD]} {
    puts "  SKIP  no capsule-incomplete fixture / amberuni.rad"
    puts "INCOMPLETE-RESULT pass=$pass fail=$fail"
    quit
}
if {![file executable $REF]} {
    puts "  SKIP  no HOLE binary at $REF (the fallback engine does not emit this warning)"
    puts "INCOMPLETE-RESULT pass=$pass fail=$fail"
    quit
}

set work [file join [::VMDHole::get_temp_base] "vmdhole_incomplete_[pid]"]
file mkdir $work
set mid [mol new $PDB waitfor all]

# Same cards both times; ONLY pore_method changes between the two runs.
proc run_method {meth mid work RAD REF} {
    set w [file join $work $meth]
    file mkdir $w
    array set ::VMDHole::state [list \
        molid $mid frame_spec 0 selection all \
        hole_exec $REF radius_file $RAD \
        cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 5.0 \
        random_seed 1 pore_method $meth display_mode none \
        work_dir $w keep_input_pdb 1 save_results 1 extra_cards {} ignore {}]
    set ::VMDHole::results {}
    set ::VMDHole::result_frames {}
    set ::VMDHole::msg_log {}
    if {[catch {::VMDHole::run_analysis} err]} { return [list error $err] }
    if {![llength $::VMDHole::result_frames]} { return [list norows ""] }
    set f [lindex $::VMDHole::result_frames 0]
    return [list ok [dict get $::VMDHole::results $f run_dir]]
}

# 1. CAPSULE - the case that warns.
lassign [run_method capsule $mid $work $RAD $REF] st cap_dir
chk "capsule run produced a result" $st "ok"

if {$st eq "ok"} {
    # The sidecar has to EXIST as a file: an absent one is indistinguishable
    # from "no warning", which is how a run whose hole_out.txt was already
    # deleted would read. Every pre-existing run dir on disk is in that state,
    # so this must be checked on a run made HERE.
    set side [file join $cap_dir vmdhole_incomplete.dat]
    chk "the job wrote an incomplete sidecar" [file exists $side] 1

    set raw [::VMDHole::_hole_frame_incomplete $cap_dir]
    note "sidecar = \"$raw\""
    chk "...and HOLE's warning is in it" \
        [expr {[string first "MAY BE INCOMPLETE" $raw] >= 0}] 1
    # The DIRECTION is the part that matters on screen: "+VE only" is the shape
    # of a lining that stops halfway up the channel.
    chk "...naming the +VE direction" [::VMDHole::_incomplete_dirs_label $raw] "+"

    # The warning is HOLE's capsule conductance code failing a strict
    # comparison on roundoff (hcapgr.f), not the search stopping: on this very
    # fixture the profile spans both sides of CPOINT. So the sidecar is kept
    # verbatim, but the run must NOT tell the user the search stopped early -
    # the gate keeps a direction only when its side of the .sph is empty.
    lassign [::VMDHole::_sph_side_counts $cap_dir] np nn
    note "stored records: +VE $np, -VE $nn"
    chk "both sides of the capsule .sph hold records" [expr {$np >= 2 && $nn >= 2}] 1
    chk "the gate drops a warning both sides contradict" \
        [::VMDHole::_incomplete_real_dirs $cap_dir $raw] ""
    set surfaced 0
    foreach e $::VMDHole::msg_log {
        if {[string first "stopped its own search early" [lindex $e 2]] >= 0} { set surfaced 1 }
    }
    chk "...so nothing is surfaced for a complete profile" $surfaced 0

    # The gate must still pass a REAL early stop through. Same sidecar text,
    # a .sph whose +VE side is empty: one start record, records on -VE only.
    set gdir [file join $work gateprobe]
    file mkdir $gdir
    set gfh [open [file join $gdir hole_out.sph] w]
    puts $gfh "ATOM      1  QC1 SPH S   0       0.000   0.000   0.000  3.00  0.00"
    puts $gfh "ATOM      1  QC2 SPH S   0       0.500   0.000   0.000  3.00  0.00"
    foreach r {-1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21} {
        puts $gfh [format "ATOM      1  QC1 SPH S%4d       0.000   0.000 %7.3f  3.00  0.00" $r [expr {$r*0.5}]]
        puts $gfh [format "ATOM      1  QC2 SPH S%4d       0.500   0.000 %7.3f  3.00  0.00" $r [expr {$r*0.5}]]
    }
    close $gfh
    lassign [::VMDHole::_sph_side_counts $gdir] gp gn
    chk "side counts read the synthetic .sph" [list $gp $gn] {0 21}
    chk "an empty +VE side keeps the +VE warning" \
        [::VMDHole::_incomplete_real_dirs $gdir "MAY BE INCOMPLETE IN +VE DIRECTION"] \
        "MAY BE INCOMPLETE IN +VE DIRECTION"
    chk "...but not a -VE warning the -VE records contradict" \
        [::VMDHole::_incomplete_real_dirs $gdir \
            "MAY BE INCOMPLETE IN +VE DIRECTION MAY BE INCOMPLETE IN -VE DIRECTION"] \
        "MAY BE INCOMPLETE IN +VE DIRECTION"
}

# 2. SPHERICAL - the control. Same structure, same cards, no warning.
lassign [run_method circular $mid $work $RAD $REF] st2 sph_dir
chk "spherical run produced a result" $st2 "ok"

if {$st2 eq "ok"} {
    set raw2 [::VMDHole::_hole_frame_incomplete $sph_dir]
    note "spherical sidecar = \"$raw2\""
    chk "spherical does NOT warn on the same structure" $raw2 ""
    set surfaced2 0
    foreach e $::VMDHole::msg_log {
        if {[string first "stopped its own search early" [lindex $e 2]] >= 0} { set surfaced2 1 }
    }
    chk "...and nothing is surfaced for it" $surfaced2 0
}

# 3. The SILENT twin: a profile covering only one side of CPOINT, which HOLE
#    does not warn about at all. Driven on a synthetic profile + axis file so it
#    does not need the multi-GB trajectory the real case was found on; the real
#    discrimination (capsule seed 1 flags, seed 2 / spherical / Connolly stay
#    silent) was measured on that trajectory and is recorded in the run summary.
set axdir [file join $work axprobe]
file mkdir $axdir
set afh [open [file join $axdir vmdhole_frame_axis.dat] w]
puts -nonewline $afh "cpoint|0 0 4|cvect|0 0 1"
close $afh
# origin projects to +4.0. Starved on the -ve side: 1 A below, 14 A above.
set p_starved [dict create xvalues {3.0 5.0 10.0 18.0}]
chk "a one-sided profile is flagged, and names the starved side" \
    [::VMDHole::_profile_one_sided $axdir $p_starved] "-"
# The mirror image.
chk "...and the other side too" \
    [::VMDHole::_profile_one_sided $axdir [dict create xvalues {-10.0 -2.0 3.0 4.5}]] "+"
# A balanced profile must stay silent - without this the check would pass
# against a detector that flagged everything.
chk "a balanced profile is NOT flagged" \
    [::VMDHole::_profile_one_sided $axdir [dict create xvalues {-20.0 0.0 10.0 18.0}]] ""
# Lopsided but not starved: 3 A on the short side clears the 2 A floor.
chk "merely lopsided is NOT flagged" \
    [::VMDHole::_profile_one_sided $axdir [dict create xvalues {1.0 6.0 12.0 18.0}]] ""
# No axis file (an imported or older run) - no evidence, so no claim.
chk "no axis file yields no claim" \
    [::VMDHole::_profile_one_sided [file join $work nosuch] $p_starved] ""

# 4. The label helper reads BOTH signs out of HOLE's own wording.
chk "both directions collapse to one label" \
    [::VMDHole::_incomplete_dirs_label \
        "MAY BE INCOMPLETE IN +VE DIRECTION MAY BE INCOMPLETE IN -VE DIRECTION"] "+ and -"
chk "-VE alone" [::VMDHole::_incomplete_dirs_label "MAY BE INCOMPLETE IN -VE DIRECTION"] "-"
chk "no warning yields no label" [::VMDHole::_incomplete_dirs_label ""] ""

catch {mol delete $mid}
catch {file delete -force $work}
puts "INCOMPLETE-RESULT pass=$pass fail=$fail"
# MUST quit: `vmd -e` runs the script then enters its interactive text prompt,
# which blocks on the suite's stdin. Same ending as the other run_analysis tests.
quit
