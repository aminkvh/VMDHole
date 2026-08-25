# Dump the pure-Tcl DH triangulation in the reference's canonical shape, so the
# Tcl is compared against MOLE ITSELF rather than only against its own C.
#
# Usage: dh_cells.tcl ATOMS.txt
set here [file dirname [info script]]
# The engine is inlined in vmdhole.tcl now (one shipped script + the
# compiled binaries), so this extracts it instead of sourcing modules.
source [file join [file dirname [file normalize [info script]]] mole_tcl_engine.tcl]
::VMDHole::Mole::read_atoms [lindex $argv 0] A
set piv [::VMDHole::Mole::pivot_points A rad]
set np [expr {[llength $piv] / 3}]
if {[::VMDHole::Mole::dh_build M $piv $np] != 0} { puts stderr "dh_build failed"; exit 1 }
set out {}
for {set t 0} {$t < $M(nt)} {incr t} {
    set cx 0.0; set cy 0.0; set cz 0.0; set vs {}; set nn 0
    for {set k 0} {$k < 4} {incr k} {
        set a [lindex $M(tv) [expr {4*$t+$k}]]
        set b [expr {3*$a}]
        set x [lindex $piv $b]; set y [lindex $piv [expr {$b+1}]]
        set z [lindex $piv [expr {$b+2}]]
        set cx [expr {$cx+$x}]; set cy [expr {$cy+$y}]; set cz [expr {$cz+$z}]
        lappend vs [format "%.3f,%.3f,%.3f" $x $y $z]
        if {[lindex $M(tn) [expr {4*$t+$k}]] < 0} { incr nn }
    }
    lappend out [format "%.4f %.4f %.4f %d %s" [expr {$cx/4}] [expr {$cy/4}] \
                     [expr {$cz/4}] $nn [join $vs ";"]]
}
puts [join $out "\n"]
