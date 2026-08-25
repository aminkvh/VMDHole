# StrictInterior, checked against MOLE's OWN numbers rather than against ours.
#
# MOLE's CLI cannot set StrictInterior - the property is declared once and never
# assigned - so the flag was reached by rebuilding MOLE with it forced true
# (ComplexParameters.cs:82) into a SEPARATE tree. Everything expected here is
# what that build printed.
#
# The point of the flag: it makes MaxClearance the maximum over the adjacent
# EDGES' Clearance, and those are all still zero at that moment because the line
# that would compute them (ComplexComputation.cs:345/348, under the author's own
# "if strict, need to update the edges 1st!") is commented out. Zero is below
# 2 * InteriorThreshold, so the interior removal takes every surviving
# tetrahedron and no cavity is left.
#
# "Nothing came out" is also what a crash produces, so the discriminating check
# is not the zero. It is the SurfaceCavity, which is snapshotted BEFORE the
# interior removal and therefore survives - and grows, because MaxClearance = 0
# also weakens the probe peel to its second term alone. Its volume is a
# non-trivial number that has to come out right.
#
# Usage: test_mole_strict_interior.tcl ATOMS_1BL8.txt

set here [file dirname [info script]]
# The engine is inlined in vmdhole.tcl now (one shipped script + the
# compiled binaries), so this extracts it instead of sourcing modules.
source [file join [file dirname [file normalize [info script]]] mole_tcl_engine.tcl]

set fails 0
proc report {label ok {detail ""}} {
    global fails
    if {!$ok} { incr fails }
    puts [format "  %-52s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
}

# MOLE's own output for 1BL8, from the stock build and from the
# StrictInterior-forced rebuild.
#   channels/voids: the Cavity and Void entries in cavities.xml
#   volume:         the MolecularSurface entry's Volume attribute
set WANT(0) {4 9 24837.879}
set WANT(1) {0 0 26633.820}

# Cavity.Volume = graph.Vertices.Sum(f => f.Volume). The graph is built from
# EDGES, so a surface member with no surviving neighbour is not a vertex and
# contributes nothing.
proc surface_volume {cxVar} {
    upvar 1 $cxVar C
    set v 0.0
    for {set t 0} {$t < $C(nt)} {incr t} {
        if {![lindex $C(surface) $t]} continue
        set deg 0
        for {set k 0} {$k < 4} {incr k} {
            set nb [lindex $C(tn) [expr {4 * $t + $k}]]
            if {$nb >= 0 && [lindex $C(surface) $nb]} { incr deg }
        }
        if {$deg == 0} continue
        set v [expr {$v + [lindex $C(volume) $t]}]
    }
    return $v
}

set atoms [lindex $argv 0]
::VMDHole::Mole::read_atoms $atoms A
set piv [::VMDHole::Mole::pivot_points A rad src pbfac pfree pbb pres restab]
set np [expr {[llength $piv] / 3}]
# MOLE's own triangulation, which is what the C engine uses.
::VMDHole::Mole::dh_build M $piv $np

foreach strict {0 1} {
    lassign $WANT($strict) wch wvd wvol
    ::VMDHole::Mole::complex_build_dh C M $piv $rad 3.0 1.25 8 $strict
    lassign [::VMDHole::Mole::cavities C cav 8 5.0] nc nch nvd
    set vol [format "%.3f" [surface_volume C]]
    report "StrictInterior=$strict: $wch channels, $wvd voids (MOLE's)" \
           [expr {$nch == $wch && $nvd == $wvd}] "(got $nch/$nvd)"
    report "StrictInterior=$strict: SurfaceCavity volume $wvol (MOLE's)" \
           [expr {$vol eq $wvol}] "(got $vol)"
}

exit [expr {$fails ? 1 : 0}]
