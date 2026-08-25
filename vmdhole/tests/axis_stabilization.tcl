# Independent check of per-frame axis stabilization.
#
# Builds a trajectory whose frames are the SAME structure rotated by a known,
# increasing angle about a known axis. A correct stabilizer must return a cvect
# that rotates with the structure by that same angle; a broken one returns the
# frame-0 vector unchanged. This is the test the plan demanded, because a
# pre-superposed trajectory cannot distinguish the two.
# Run:  vmd -dispdev text -e tests/axis_stabilization.tcl -args <REPO> <PDB>
#
# ⚠️ Two setup traps this test exists to document, both of which make the
#    stabilizer look BROKEN when it is fine:
#      1. An NMR ensemble must be collapsed to one model first, or the "rotated"
#         frames land on different conformers.
#      2. CVECT Stabilize needs cvect_def_p1/p2 to RESOLVE TO ATOMS. A literal
#         direction has nothing to re-fit, and a slab selection outside the
#         structure's real extent silently selects nothing. Either way the
#         feature no-ops with no warning - check the z range first.
lassign $argv REPO PDB
source [file join $REPO vmdhole vmdhole.tcl]
set ::VMDHole::config_file /tmp/.vmdhole_config_stabtest

set molid [mol new $PDB waitfor all]
# 1GRM is an NMR ensemble - drop every model but the first, or the "rotated"
# frames below would be written onto DIFFERENT conformers and the test would be
# measuring conformational change, not the stabilizer.
if {[molinfo $molid get numframes] > 1} { animate delete beg 1 end -1 $molid }
set nat [molinfo $molid get numatoms]

# ---- build the rotated frames -------------------------------------------------
set NFRAMES 6
set STEP    10.0                    ;# degrees per frame, about +x
set all [atomselect $molid all]
set base [$all get {x y z}]
for {set f 1} {$f < $NFRAMES} {incr f} {
    animate dup frame 0 $molid
    $all frame $f
    $all set {x y z} $base
    $all move [transaxis x [expr {$f * $STEP}] deg]
}
# prove the frames really differ by the imposed rotation
for {set f 0} {$f < $NFRAMES} {incr f} {
    $all frame $f
    set c [measure center $all]
    puts [format "STAB_SETUP f%d center=(%.3f %.3f %.3f)" $f {*}$c]
}
$all delete
puts "STAB_FRAMES [molinfo $molid get numframes] (imposed +${STEP} deg/frame about x)"

namespace eval ::VMDHole {
    variable state
    set state(molid)      $::molid
    set state(selection)  "protein"
    # CPOINT must sit INSIDE the protein or the scoped CA search finds nothing and
    # Stabilize silently falls back to static. CVECT Stabilize additionally needs
    # the two DEFINING POINTS (cvect_def_p1/p2) - a literal direction has nothing
    # to re-fit, which is exactly what _stab_init says.
    set state(cpoint)       "protein"
    set state(cvect)        "0 0 1"
    set state(cvect_def_p1) "protein and z < 1"
    set state(cvect_def_p2) "protein and z > 8"
    set state(track_radius) 8.0
    set state(stab_radius_inner) 10.0
    set state(stab_radius_outer) 15.0
}

proc ang {a b} {
    lassign $a ax ay az; lassign $b bx by bz
    set na [expr {sqrt($ax*$ax+$ay*$ay+$az*$az)}]; set nb [expr {sqrt($bx*$bx+$by*$by+$bz*$bz)}]
    if {$na == 0 || $nb == 0} { return -1 }
    set d [expr {($ax*$bx+$ay*$by+$az*$bz)/($na*$nb)}]
    if {$d > 1.0} {set d 1.0}; if {$d < -1.0} {set d -1.0}
    expr {acos($d) * 57.2957795130823}
}

proc run_mode {label cp_stab cv_stab cp_track cv_exact} {
    global molid NFRAMES STEP
    set ::VMDHole::state(stabilize_cpoint) $cp_stab
    set ::VMDHole::state(stabilize_cvect)  $cv_stab
    set ::VMDHole::state(track_cpoint)     $cp_track
    set ::VMDHole::state(cvect_exact)      $cv_exact
    molinfo $molid set frame 0
    if {[catch {::VMDHole::_stab_init $molid} e]} { puts "  $label: _stab_init FAILED: $e"; return }
    set cv0 ""; set cp0 ""
    puts "  --- $label ---"
    for {set f 0} {$f < $NFRAMES} {incr f} {
        molinfo $molid set frame $f
        if {[catch {::VMDHole::frame_axis $molid $f} r]} { puts "   f$f frame_axis ERROR: $r"; continue }
        lassign $r cp cv
        if {$f == 0} { set cv0 $cv; set cp0 $cp }
        set da [ang $cv0 $cv]
        set dcp [expr {[llength $cp0] == 3 ? \
            sqrt(pow([lindex $cp 0]-[lindex $cp0 0],2)+pow([lindex $cp 1]-[lindex $cp0 1],2)+pow([lindex $cp 2]-[lindex $cp0 2],2)) : -1}]
        puts [format "   f%d  cvect=(%s)  angle_vs_f0=%6.2f deg (expected %5.1f)   |cpoint-cp0|=%7.3f A" \
              $f $cv $da [expr {$f*$STEP}] $dcp]
    }
}

run_mode "STATIC (all off)"        0 0 0 0
run_mode "CVECT Stabilize"          0 1 0 0
run_mode "CPOINT Stabilize"         1 0 0 0
run_mode "CPOINT Track"             0 0 1 0
run_mode "CVECT Exact"              0 0 0 1
puts "STAB_DONE"
quit
