# The Tcl Delaunay must produce the C's mesh EXACTLY - the same tetrahedra in
# the same slots, the same neighbour links, the same free list.
#
# Not "a Delaunay triangulation of the same points": tetrahedron indices are
# load-bearing downstream (component iteration is index-ordered, and a
# MaxClearance tie is settled by lowest index), so a mesh that differs only in
# slot layout produces different tunnels. Matching the live set as a SET would
# pass while being wrong.
#
# Two stages, checked separately, because a combined comparison cannot say which
# half broke:
#   A  quantised integers in  ->  tetrahedra out      (Bowyer-Watson only)
#   B  atom table in          ->  quantised integers  (jitter + quantisation only)
#
# Usage: test_mole_delaunay.tcl POINTS.txt TETRA.txt [ATOMS.txt]

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

# ---------------------------------------------------------------- reference
set fp [open [lindex $argv 0] r]
gets $fp npt
set points {}
while {[gets $fp line] >= 0} {
    if {[string trim $line] eq ""} continue
    foreach v $line { lappend points $v }
}
close $fp

set ft [open [lindex $argv 1] r]
gets $ft hdr
lassign $hdr ref_nt ref_nfree ref_last
set ref_tv {}; set ref_tnb {}; set ref_tdead {}; set ref_free {}
set row 0
while {[gets $ft line] >= 0} {
    if {[string trim $line] eq ""} continue
    if {$row < $ref_nt} {
        lassign $line d v0 v1 v2 v3 n0 n1 n2 n3
        lappend ref_tdead $d
        lappend ref_tv $v0 $v1 $v2 $v3
        lappend ref_tnb $n0 $n1 $n2 $n3
    } else {
        lappend ref_free $line
    }
    incr row
}
close $ft

# ------------------------------------------------- stage A: Bowyer-Watson
set t0 [clock milliseconds]
set rc [::VMDHole::Mole::dt_build_quantised M $points $npt]
set ms [expr {[clock milliseconds] - $t0}]
report "dt_build_quantised returns 0 on $npt points" [expr {$rc == 0}] "(rc $rc)"

report "tetrahedron slots = $ref_nt" [expr {$M(nt) == $ref_nt}] "(got $M(nt))"

set bad 0; set first ""
for {set i 0} {$i < $ref_nt && $i < $M(nt)} {incr i} {
    set ok 1
    if {[lindex $M(tdead) $i] != [lindex $ref_tdead $i]} { set ok 0 }
    for {set k 0} {$k < 4 && $ok} {incr k} {
        set j [expr {4 * $i + $k}]
        if {[lindex $M(tv) $j] != [lindex $ref_tv $j]} { set ok 0 }
        if {[lindex $M(tnb) $j] != [lindex $ref_tnb $j]} { set ok 0 }
    }
    if {!$ok} {
        incr bad
        if {$first eq ""} {
            set first [format "slot %d: tcl {%s|%s|%s} c {%s|%s|%s}" $i \
                [lindex $M(tdead) $i] [lrange $M(tv) [expr {4*$i}] [expr {4*$i+3}]] \
                [lrange $M(tnb) [expr {4*$i}] [expr {4*$i+3}]] \
                [lindex $ref_tdead $i] [lrange $ref_tv [expr {4*$i}] [expr {4*$i+3}]] \
                [lrange $ref_tnb [expr {4*$i}] [expr {4*$i+3}]]]
        }
    }
}
report "every slot identical (vertices, neighbours, dead)" [expr {$bad == 0}] \
       "($bad slots differ; $first)"
report "free list identical ($ref_nfree entries)" \
       [expr {$M(freet) eq $ref_free}] \
       "(got [llength $M(freet)] entries)"
report "walk start M(last) = $ref_last" [expr {$M(last) == $ref_last}] "(got $M(last))"

# ------------------------------------- stage B: jitter and quantisation only
if {[llength $argv] > 2} {
    ::VMDHole::Mole::read_atoms [lindex $argv 2] A
    set piv [::VMDHole::Mole::pivot_points A]
    set np [expr {[llength $piv] / 3}]
    if {$np > $npt} { set np $npt; set piv [lrange $piv 0 [expr {3 * $npt - 1}]] }
    set got [::VMDHole::Mole::quantise_points $piv $np]
    set nq 0
    foreach a $got b $points { if {$a != $b} { incr nq } }
    report "quantised points identical ($np pivots + 4 corners)" \
           [expr {$got eq $points}] "($nq of [llength $points] differ)"

    # The combined entry point is not exercised by the two-stage checks above,
    # and it hands a caller-named array down through a second upvar level into
    # an `array unset` - a Tcl footgun this project has met before.
    set rc2 [::VMDHole::Mole::dt_build B $piv $np]
    report "dt_build (combined) agrees with the two-stage path" \
           [expr {$rc2 == 0 && $B(nt) == $M(nt) && $B(tv) eq $M(tv)
                  && $B(tnb) eq $M(tnb) && $B(tdead) eq $M(tdead)
                  && $B(freet) eq $M(freet) && $B(last) == $M(last)}] \
           "(rc $rc2, nt $B(nt) vs $M(nt))"
}

puts [format "  %-52s %d ms" "stage A build time" $ms]
exit [expr {$fails ? 1 : 0}]
