# The Tcl cavity pipeline must reproduce the C's, quantity by quantity.
#
# Comparing only the cavity counts would be a weak check: those survived several
# of the arithmetic-form errors this port has already hit (normalise-by-multiply,
# the snapshot in RemoveShallowVertices), because a handful of tetrahedra moving
# in or out rarely changes how many components there are. So every per-tetra and
# per-edge quantity is compared, and the doubles are compared EXACTLY - both
# sides evaluate the same IEEE operations in the same order, so anything less
# than exact equality would hide the very differences the C's comments warn about.
#
# Usage: test_mole_complex.tcl ATOMS.txt COMPLEX.txt [max_pivots]

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
set f [open [lindex $argv 1] r]
set ref_rad {}; set ref_T {}; set ref_E {}; set ref_C {}
while {[gets $f line] >= 0} {
    switch -- [lindex $line 0] {
        H { lassign $line _ ref_np ref_nt ref_nsurf ref_nc ref_nch ref_nvd }
        R { lappend ref_rad [lindex $line 1] }
        T { lappend ref_T [lrange $line 1 end] }
        E { lappend ref_E [lrange $line 1 end] }
        C { lappend ref_C [lrange $line 1 end] }
    }
}
close $f

# ------------------------------------------------------------------- build
::VMDHole::Mole::read_atoms [lindex $argv 0] A
set piv [::VMDHole::Mole::pivot_points A rad]
set np [expr {[llength $piv] / 3}]
if {[llength $argv] > 2 && [lindex $argv 2] > 3 && [lindex $argv 2] < $np} {
    set np [lindex $argv 2]
    set piv [lrange $piv 0 [expr {3 * $np - 1}]]
    set rad [lrange $rad 0 [expr {$np - 1}]]
}
report "pivots = $ref_np" [expr {$np == $ref_np}] "(got $np)"

set nbad 0
foreach a $rad b $ref_rad { if {$a != $b} { incr nbad } }
report "vdW radii identical ($np atoms)" [expr {$nbad == 0}] "($nbad differ)"

set t0 [clock milliseconds]
::VMDHole::Mole::dt_build M $piv $np
set tms [expr {[clock milliseconds] - $t0}]
set t0 [clock milliseconds]
::VMDHole::Mole::complex_build C M $piv $rad 3.0 1.25 8
lassign [::VMDHole::Mole::cavities C cav 8 5.0] nc nch nvd
set cms [expr {[clock milliseconds] - $t0}]

report "finite tetrahedra = $ref_nt" [expr {$C(nt) == $ref_nt}] "(got $C(nt))"
report "SurfaceCavity snapshot = $ref_nsurf" [expr {$C(n_surface) == $ref_nsurf}] \
       "(got $C(n_surface))"

# ------------------------------------------- per tetrahedron, every quantity
# Named individually: "something differs" would not say whether the fault is in
# the geometry, the peel, or the depth BFS, and those are separate transcriptions.
array set diff {}
set fields {alive boundary depth depthlen comp cx cy cz vx vy vz volume maxclear
            v0 v1 v2 v3}
foreach fld $fields { set diff($fld) 0 }
set firstbad ""
for {set i 0} {$i < $ref_nt && $i < $C(nt)} {incr i} {
    set r [lindex $ref_T $i]
    lassign $r ri alive bnd dep dlen comp rcx rcy rcz rvx rvy rvz rvol rmax rv0 rv1 rv2 rv3
    set b3 [expr {3 * $i}]
    set b4 [expr {4 * $i}]
    foreach {fld got want} [list \
        alive    [lindex $C(alive) $i]     $alive \
        boundary [lindex $C(boundary) $i]  $bnd \
        depth    [lindex $C(depth) $i]     $dep \
        depthlen [lindex $C(depthlen) $i]  $dlen \
        comp     [lindex $C(comp) $i]      $comp \
        cx  [lindex $C(center) $b3]              $rcx \
        cy  [lindex $C(center) [expr {$b3+1}]]   $rcy \
        cz  [lindex $C(center) [expr {$b3+2}]]   $rcz \
        vx  [lindex $C(vcenter) $b3]             $rvx \
        vy  [lindex $C(vcenter) [expr {$b3+1}]]  $rvy \
        vz  [lindex $C(vcenter) [expr {$b3+2}]]  $rvz \
        volume   [lindex $C(volume) $i]    $rvol \
        maxclear [lindex $C(maxclear) $i]  $rmax \
        v0 [lindex $C(tv) $b4]             $rv0 \
        v1 [lindex $C(tv) [expr {$b4+1}]]  $rv1 \
        v2 [lindex $C(tv) [expr {$b4+2}]]  $rv2 \
        v3 [lindex $C(tv) [expr {$b4+3}]]  $rv3] {
        if {$got != $want} {
            incr diff($fld)
            if {$firstbad eq ""} { set firstbad "tet $i $fld: tcl $got c $want" }
        }
    }
}
set nt_bad 0
foreach fld $fields { incr nt_bad $diff($fld) }
set detail ""
foreach fld $fields { if {$diff($fld)} { append detail " $fld=$diff($fld)" } }
report "every per-tetrahedron quantity identical" [expr {$nt_bad == 0}] \
       "($detail; first: $firstbad)"

# --------------------------------------------------- per edge, every quantity
array set ediff {tn 0 eclear 0 elen 0 eweight 0 evweight 0}
set efirst ""
for {set i 0} {$i < $ref_nt && $i < $C(nt)} {incr i} {
    set r [lindex $ref_E $i]
    for {set k 0} {$k < 4} {incr k} {
        set o [expr {1 + 5 * $k}]
        set e [expr {4 * $i + $k}]
        foreach {fld got want} [list \
            tn       [lindex $C(tn) $e]       [lindex $r $o] \
            eclear   [lindex $C(eclear) $e]   [lindex $r [expr {$o+1}]] \
            elen     [lindex $C(elen) $e]     [lindex $r [expr {$o+2}]] \
            eweight  [lindex $C(eweight) $e]  [lindex $r [expr {$o+3}]] \
            evweight [lindex $C(evweight) $e] [lindex $r [expr {$o+4}]]] {
            if {$got != $want} {
                incr ediff($fld)
                if {$efirst eq ""} { set efirst "tet $i face $k $fld: tcl $got c $want" }
            }
        }
    }
}
set e_bad 0
set edetail ""
foreach fld {tn eclear elen eweight evweight} {
    incr e_bad $ediff($fld)
    if {$ediff($fld)} { append edetail " $fld=$ediff($fld)" }
}
report "every per-edge quantity identical" [expr {$e_bad == 0}] \
       "($edetail; first: $efirst)"

# ---------------------------------------------------------------- cavities
report "components = $ref_nc, channels = $ref_nch, voids = $ref_nvd" \
       [expr {$nc == $ref_nc && $nch == $ref_nch && $nvd == $ref_nvd}] \
       "(got $nc/$nch/$nvd)"
set cbad 0; set cfirst ""
for {set i 0} {$i < $ref_nc && $i < $nc} {incr i} {
    lassign [lindex $ref_C $i] _ rcount rdepth rbnd rvol rdlen
    lassign [lindex $cav $i] count depth bnd vol dlen
    if {$count != $rcount || $depth != $rdepth || $bnd != $rbnd
        || $vol != $rvol || $dlen != $rdlen} {
        incr cbad
        if {$cfirst eq ""} {
            set cfirst "cavity $i: tcl {$count $depth $bnd $vol $dlen} c {$rcount $rdepth $rbnd $rvol $rdlen}"
        }
    }
}
report "every cavity identical (count, depth, volume, length)" [expr {$cbad == 0}] \
       "($cbad differ; $cfirst)"

puts [format "  %-52s %d ms triangulation + %d ms complex" "build time" $tms $cms]
exit [expr {$fails ? 1 : 0}]
