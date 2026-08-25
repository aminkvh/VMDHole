# How fast is Tcl at the exact predicates a Bowyer-Watson Delaunay is made of?
# Coordinates are quantised integers, so Tcl's native bignums give exactness
# for free - the question is throughput, not correctness.
proc orient3d {ax ay az bx by bz cx cy cz dx dy dz} {
    set a1 [expr {$ax-$dx}]; set a2 [expr {$ay-$dy}]; set a3 [expr {$az-$dz}]
    set b1 [expr {$bx-$dx}]; set b2 [expr {$by-$dy}]; set b3 [expr {$bz-$dz}]
    set c1 [expr {$cx-$dx}]; set c2 [expr {$cy-$dy}]; set c3 [expr {$cz-$dz}]
    expr {$a1*($b2*$c3-$b3*$c2) - $a2*($b1*$c3-$b3*$c1) + $a3*($b1*$c2-$b2*$c1)}
}
proc insphere {ax ay az bx by bz cx cy cz dx dy dz ex ey ez} {
    set ax [expr {$ax-$ex}]; set ay [expr {$ay-$ey}]; set az [expr {$az-$ez}]
    set bx [expr {$bx-$ex}]; set by [expr {$by-$ey}]; set bz [expr {$bz-$ez}]
    set cx [expr {$cx-$ex}]; set cy [expr {$cy-$ey}]; set cz [expr {$cz-$ez}]
    set dx [expr {$dx-$ex}]; set dy [expr {$dy-$ey}]; set dz [expr {$dz-$ez}]
    set al [expr {$ax*$ax+$ay*$ay+$az*$az}]; set bl [expr {$bx*$bx+$by*$by+$bz*$bz}]
    set cl [expr {$cx*$cx+$cy*$cy+$cz*$cz}]; set dl [expr {$dx*$dx+$dy*$dy+$dz*$dz}]
    set ab [expr {$ax*$by-$bx*$ay}]; set bc [expr {$bx*$cy-$cx*$by}]
    set cd [expr {$cx*$dy-$dx*$cy}]; set da [expr {$dx*$ay-$ax*$dy}]
    set ac [expr {$ax*$cy-$cx*$ay}]; set bd [expr {$bx*$dy-$dx*$by}]
    set abc [expr {$az*$bc-$bz*$ac+$cz*$ab}]; set bcd [expr {$bz*$cd-$cz*$bd+$dz*$bc}]
    set cda [expr {$cz*$da+$dz*$ac+$az*$cd}]; set dab [expr {$dz*$ab+$az*$bd+$bz*$da}]
    expr {($dl*$abc-$cl*$dab)+($bl*$cda-$al*$bcd)}
}
# realistic magnitudes: 1e-5 A grid, protein-sized coordinates
set N 20000
set t0 [clock milliseconds]
for {set i 0} {$i < $N} {incr i} {
    orient3d 1234567 2345678 3456789  1234000 2345000 3456000 \
             1235000 2346000 3457000  1236000 2347000 3458000
}
set t1 [clock milliseconds]
for {set i 0} {$i < $N} {incr i} {
    insphere 1234567 2345678 3456789  1234000 2345000 3456000 \
             1235000 2346000 3457000  1236000 2347000 3458000 \
             1234500 2345500 3456500
}
set t2 [clock milliseconds]
set o [expr {($t1-$t0)*1000.0/$N}]; set s [expr {($t2-$t1)*1000.0/$N}]
puts [format "  orient3d  %.1f us/call" $o]
puts [format "  insphere  %.1f us/call" $s]
# a 3809-point build used 407284 orient3d + 482904 insphere in C
set est [expr {(407284*$o + 482904*$s)/1e6}]
puts [format "\n  predicates alone for one 1MXT-sized build: %.1f s" $est]
puts [format "  (C does the whole triangulation in 0.04 s)"]
puts [format "  ratio, predicates only: %.0fx slower" [expr {$est/0.04}]]
