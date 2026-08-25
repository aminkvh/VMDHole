# _hydro_qco_c must return exactly what compute_hydration's own Tcl loop returns.
#
# Pure tclsh - no VMD, no Tk, no trajectory fixture - because the block under
# test is a pure function of its arguments. Sizes match a real solvated frame
# (506 envelope samples, 2213 waters over 3-atom residues); at a much smaller
# envelope the C call's process+I/O overhead makes it the SLOWER of the two,
# which is fine but makes the timing line meaningless.
#
# Compare NUMERICALLY, not as strings: C writes %.17g and Tcl prints
# shortest-round-trip, so identical doubles have different text.

set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
source [file join [file dirname $here] vmdhole.tcl]

namespace eval ::VMDHole {}
# HOLE install location: $VMDHOLE_HOLE_EXE_DIR, else ~/hole2/exe (the install
# layout the docs describe). Nothing machine-specific: if the binaries are not
# there, the C-vs-Tcl comparison below simply exercises the Tcl path alone and
# the group-level guard reports the skip.
set _exedir [expr {[info exists ::env(VMDHOLE_HOLE_EXE_DIR)] && $::env(VMDHOLE_HOLE_EXE_DIR) ne ""
                   ? $::env(VMDHOLE_HOLE_EXE_DIR) : [file join $::env(HOME) hole2 exe]}]
set ::VMDHole::state(sos_triangle_exec) [file join $_exedir sos_triangle]
set ::VMDHole::state(hole_exec) [file join $_exedir hole]

# Deterministic pseudo-random inputs, several residues per molecule.
expr {srand(20260811)}
set penv {}
for {set i 0} {$i < 506} {incr i} {
    lappend penv [list [expr {-15.0 + $i*0.06}] [expr {2.0 + 1.5*sin($i/7.0)}]]
}
set cmin -14.0 ; set cmax 14.0 ; set dcap 0.0
set mx 0.3 ; set my -0.2 ; set mz 0.1
set n [expr {sqrt(0.02*0.02 + 0.03*0.03 + 1.0)}]
set ux [expr {0.02/$n}] ; set uy [expr {0.03/$n}] ; set uz [expr {1.0/$n}]

set wpos {} ; set wres {}
for {set r 0} {$r < 2213} {incr r} {
    for {set a 0} {$a < 3} {incr a} {
        lappend wres $r
        lappend wpos [list [expr {(rand()-0.5)*24}] [expr {(rand()-0.5)*24}] [expr {(rand()-0.5)*34}]]
    }
}

# Reference: the plugin's own Tcl, transcribed from compute_hydration.
proc ref_qco {wres wpos mx my mz ux uy uz cmin cmax penv dcap} {
    if {[llength $wres] == [llength $wpos] && [llength $wpos] > 0} {
        set _wcx [dict create]; set _wcy [dict create]; set _wcz [dict create]; set _wcn [dict create]
        foreach p $wpos rid $wres {
            lassign $p wx wy wz
            dict set _wcx $rid [expr {([dict exists $_wcx $rid] ? [dict get $_wcx $rid] : 0.0) + $wx}]
            dict set _wcy $rid [expr {([dict exists $_wcy $rid] ? [dict get $_wcy $rid] : 0.0) + $wy}]
            dict set _wcz $rid [expr {([dict exists $_wcz $rid] ? [dict get $_wcz $rid] : 0.0) + $wz}]
            dict set _wcn $rid [expr {([dict exists $_wcn $rid] ? [dict get $_wcn $rid] : 0)   + 1}]
        }
        set wpos {}
        foreach rid [dict keys $_wcn] {
            set c [dict get $_wcn $rid]
            lappend wpos [list [expr {[dict get $_wcx $rid]/$c}] \
                               [expr {[dict get $_wcy $rid]/$c}] \
                               [expr {[dict get $_wcz $rid]/$c}]]
        }
    }
    set qco {}
    foreach p $wpos {
        lassign $p wx wy wz
        set dxc [expr {$wx-$mx}]; set dyc [expr {$wy-$my}]; set dzc [expr {$wz-$mz}]
        set co [expr {$dxc*$ux + $dyc*$uy + $dzc*$uz}]
        if {$co < $cmin || $co > $cmax} { continue }
        set rr [::VMDHole::envelope_radius $penv $co]
        if {$rr <= 0} { continue }
        set rr_eff [expr {($dcap > 0 && $rr > $dcap) ? $dcap : $rr}]
        set perp2 [expr {$dxc*$dxc + $dyc*$dyc + $dzc*$dzc - $co*$co}]
        if {$perp2 < 0} { set perp2 0 }
        if {$perp2 <= $rr_eff*$rr_eff} { lappend qco $co }
    }
    return $qco
}

set t0 [clock microseconds]
set a [::VMDHole::_hydro_qco_c $wres $wpos $mx $my $mz $ux $uy $uz $cmin $cmax $penv $dcap]
set t1 [clock microseconds]
set b [ref_qco $wres $wpos $mx $my $mz $ux $uy $uz $cmin $cmax $penv $dcap]
set t2 [clock microseconds]

if {$a eq ""} { puts "  SKIP: hydro_project not found beside sos_triangle"; exit 0 }
puts "C accepted   : [llength $a]"
puts "Tcl accepted : [llength $b]"
set bad 0; set worst 0.0
foreach x $a y $b {
    if {abs($x-$y) > 0.0} {
        incr bad
        set d [expr {abs($x-$y)}]
        if {$d > $worst} { set worst $d }
    }
}
puts "count match  : [expr {[llength $a] == [llength $b] ? {YES} : {NO}}]"
puts "differing    : $bad  (max |delta| $worst)"
set good [expr {[llength $a]==[llength $b] && $bad==0}]
puts [format "  %-52s %s" "C water projection is bit-identical to the Tcl loop" [expr {$good ? "PASS" : "FAIL"}]]
puts [format "time C=%.1f ms  Tcl=%.1f ms" [expr {($t1-$t0)/1000.0}] [expr {($t2-$t1)/1000.0}]]

exit [expr {$good ? 0 : 1}]
