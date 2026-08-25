# The Tcl path must write the SAME OUTPUT FILE as the C engine, byte for byte.
#
# Nothing else in this suite covers it. The trace comparisons compare values;
# this compares text, and a formatting difference - a %.4f against a %.3f, a
# missing header field, a different sample count - is invisible to the one and
# fatal to the other, because the plugin's parser reads the file. Until the two
# files match the Tcl path is not a drop-in fallback for the engine.
#
# It also closes the coordinate-interface question recorded in
#: the plugin writes coordinates with `format %.3f` and
# the C engine reads them back through MOLE's own double-rounding parser, so the
# atom table here is deliberately written the way _tunnel_atoms_mole writes one,
# and both sides consume it exactly as the plugin's engine call would.
#
# Usage: test_mole_engine_file.tcl ATOMS.txt C_OUTPUT.txt TCL_OUTPUT.txt [max_pivots] [x y z]

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

set atoms [lindex $argv 0]
set cout  [lindex $argv 1]
set tout  [lindex $argv 2]

::VMDHole::Mole::read_atoms $atoms A
set piv [::VMDHole::Mole::pivot_points A rad src pbfac pfree pbb pres restab]
set np [expr {[llength $piv] / 3}]
if {[llength $argv] > 3 && [lindex $argv 3] > 3 && [lindex $argv 3] < $np} {
    set np [lindex $argv 3]
    set piv [lrange $piv 0 [expr {3 * $np - 1}]]
    set rad [lrange $rad 0 [expr {$np - 1}]]
    set pbfac [lrange $pbfac 0 [expr {$np - 1}]]
    set pfree [lrange $pfree 0 [expr {$np - 1}]]
    set pbb [lrange $pbb 0 [expr {$np - 1}]]
    set pres [lrange $pres 0 [expr {$np - 1}]]
}
# MOLE's own triangulation, which is what the C engine uses.
::VMDHole::Mole::dh_build M $piv $np
::VMDHole::Mole::complex_build_dh C M $piv $rad 3.0 1.25 8
lassign [::VMDHole::Mole::cavities C cav 8 5.0] nc nch nvd
set seed [lrange $argv 4 6]
set tunnels [::VMDHole::Mole::find_tunnels C $cav $piv $rad $np \
                 [list origin $seed origin_radius 5.0] "" $pbfac $pfree]
set lining [dict create]
set lid 0
foreach tn $tunnels {
    incr lid
    dict set lining $lid [dict replace \
        [::VMDHole::Mole::lining [lindex $tn 2] $piv $np $pres $pbb $restab] \
        restab $restab]
}
set cavlines {}
if {$A(has_names)} {
    set rank {}
    for {set i 0} {$i < [llength $cav]} {incr i} {
        set vi [lindex [lindex $cav $i] 3]
        set rr 0
        for {set q 0} {$q < [llength $cav]} {incr q} {
            set cq [lindex $cav $q]
            if {![lindex $cq 2]} continue
            set vq [lindex $cq 3]
            if {$vq > $vi || ($vq == $vi && $q < $i)} { incr rr }
        }
        lappend rank $rr
    }
    for {set c 0} {$c < [llength $cav]} {incr c} {
        lassign [lindex $cav $c] cnt dep hasb vol dlen
        if {!($hasb && $dlen > 5.0 && $dep > 8)} continue
        lassign [::VMDHole::Mole::cavity_residues C $c $pres [llength $restab]] bnd inner
        foreach line [::VMDHole::Mole::format_cavity [expr {[lindex $rank $c] + 1}] \
                          $vol $dep $dlen $bnd $inner $restab] { lappend cavlines $line }
    }
}
::VMDHole::Mole::write_tunnels $tout $A(n) $C(nt) 3.0 1.25 $tunnels \
    $A(has_names) [expr {$A(has_names) ? $lining : {}}] $cavlines

set fa [open $cout r]; set a [read $fa]; close $fa
set fb [open $tout r]; set b [read $fb]; close $fb
set la [split [string trimright $a "\n"] "\n"]
set lb [split [string trimright $b "\n"] "\n"]
report "line count [llength $la]" [expr {[llength $la] == [llength $lb]}] \
       "(Tcl wrote [llength $lb])"
set nbad 0
set first ""
foreach x $la y $lb {
    if {$x ne $y} {
        incr nbad
        if {$first eq ""} { set first "\n         c   $x\n         tcl $y" }
    }
}
report "engine output file byte-identical" [expr {$nbad == 0 && $a eq $b}] \
       "($nbad lines differ$first)"
exit [expr {$fails ? 1 : 0}]
