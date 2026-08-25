#!/usr/bin/env tclsh
# ==============================================================================
#  connolly.tcl - the CONNOLLY pore method, ported from hole2/src.
# ==============================================================================
#
#  Ported directly from the Fortran, read source-first:
#    calper.f   the two unit vectors normal to CVECT ("east"/"north")
#    concal.f   the flood-fill that grows a set of probe-sized spheres out
#               from the HOLE centre on one plane
#    coarea.f   the adaptive black/grey/white area estimator that turns that
#               sphere set into an equivalent radius (Requiv)
#
#  What Requiv IS (hole2/doc/old/hole_d04.html, "CONN" section): the radius of
#  a circle with the same area as the part of the plane a CONNR(1)-radius
#  probe sphere can reach, starting from the HOLE centre on that plane.
#
#  Loaded by hole.tcl via `source`; depends on hole::holeen (holeen.f, already
#  in hole.tcl) for every point evaluation - CONCAL is HOLEEN called a few
#  hundred times per plane, nothing more exotic.
#
# ------------------------------------------------------------------------------
#  Fidelity notes - read before changing any number in here
# ------------------------------------------------------------------------------
#  * concal.f's three-dimensional distance-squared, used for its dedup/
#    elimination checks, has a real typo in the shipped source:
#        DIST2 = DIFF1*DIFF1 + DIFF2*DIFF2 + DIFF3+DIFF3
#    (`DIFF3+DIFF3` where a correct port would write `DIFF3*DIFF3`), at all
#    four sites it appears (concal.f:292,401,486,513). Reproduced verbatim
#    below - see dist2_3d_buggy. It is INERT for an axis-aligned CVECT (every
#    point this routine builds shares the same z, so DIFF3=0 either way) but
#    is preserved for the general case, since "fix" would silently produce
#    different numbers than the reference binary for a tilted channel.
#  * coarea.f's own 2D distance-squared (dist2_2d below) has NO such bug -
#    only concal.f's 3D one does.
#  * coarea.f's `ABS(RELCCOORD).LT.0.9*SAMPLE` branch (the "sphere cuts this
#    plane at an angle" case) is unreachable from CONCAL: every point CONCAL
#    hands to COAREA is built from CENTRE via PERPVE/PERPVN displacements
#    only, so it lies exactly on the plane (RELCCOORD ~ 0, machine epsilon)
#    and always takes the "centred on the plane" branch instead. Implemented
#    anyway for fidelity (coarea.f is also called by sph_process for PEG, not
#    ported here, where it is NOT dead) but not exercised by anything below.
#  * SIZE(1) (the coarse grid cell) is hardcoded to 0.25 Angs in coarea.f,
#    independent of CONNR(2) (concal's flood-fill grid spacing) - not a typo,
#    just two different grids for two different jobs.
#  * PI in coarea.f is `2.*ASIN(1.)`, i.e. computed in single precision then
#    widened - about 1e-7 relative error against the double-precision value
#    used below. Invisible at the 4-decimal-place precision HOLE itself
#    prints Requiv to; not reproduced.
#  * Sentinels concal.f/coarea.f actually emit for REQUIV (not invented -
#    read directly off coarea.f/concal.f control flow):
#      - initial HOLE radius < probe radius: REQUIV = that HOLE radius
#        (concal.f "So using HOLE point for calcs" - no flood fill is run)
#      - any accepted point's on-plane circle radius > ENDRAD: REQUIV = 1e6
#        (coarea.f "Escaped so infinite area and Requiv!")
#      - no circle intersects the plane at all: REQUIV = 0.0 (coarea.f's
#        REQUIV is initialised to 0 and the "Cannot find any accessible
#        circles" path returns before ever assigning it - dead in practice
#        for a CONCAL caller, since point 1 is always centred exactly on the
#        plane by construction, but implemented faithfully regardless)
#  * Array bounds: concal.f's flood-fill list (SCOMAX=10000) is fatal in the
#    real program (LERR=.TRUE., the whole HOLE run aborts) - reproduced as a
#    Tcl error. coarea.f's grid-square list (SNMAX=90000) is NOT fatal there:
#    real HOLE prints a WARNING and simply stops refining, keeping the last
#    fully-computed Requiv. Reproduced the same way (see hitcap in the
#    returned dict) since that is what the source actually does.

# ==============================================================================
#  CALPER - east/north unit vectors normal to CVECT              (calper.f)
# ==============================================================================
proc hole::calper {cvect} {
    lassign $cvect vx vy vz
    if {abs($vy) < 1e-6 && abs($vz) < 1e-6} {
        set ex 0.0; set ey 1.0; set ez 0.0
    } else {
        set ex 1.0; set ey 0.0; set ez 0.0
    }
    # remove the component along CVECT, then renormalise
    set cdotp [expr {$vx*$ex + $vy*$ey + $vz*$ez}]
    set ex [expr {$ex - $cdotp*$vx}]
    set ey [expr {$ey - $cdotp*$vy}]
    set ez [expr {$ez - $cdotp*$vz}]
    set en [expr {sqrt($ex*$ex + $ey*$ey + $ez*$ez)}]
    set ex [expr {$ex/$en}]; set ey [expr {$ey/$en}]; set ez [expr {$ez/$en}]
    # PERPVN = PERPVE x CVECT (dCROSS(PERPVE, CVECT, PERPVN) in the Fortran)
    set nx [expr {$ey*$vz - $ez*$vy}]
    set ny [expr {$ez*$vx - $ex*$vz}]
    set nz [expr {$ex*$vy - $ey*$vx}]
    return [list [list $ex $ey $ez] [list $nx $ny $nz]]
}

# 3D distance-squared with concal.f's own DIFF3+DIFF3 typo - see fidelity
# notes above. Do not "fix" this; it changes the answer for tilted CVECT.
proc hole::dist2_3d_buggy {ax ay az bx by bz} {
    set d1 [expr {$ax-$bx}]; set d2 [expr {$ay-$by}]; set d3 [expr {$az-$bz}]
    return [expr {$d1*$d1 + $d2*$d2 + $d3+$d3}]
}

# ==============================================================================
#  CONCAL - flood-fill a set of probe-sized spheres on one plane   (concal.f)
# ==============================================================================
#
#  Returns {points active sconum rad0 fallback escaped} where points is a list
#  of {x y z rad} and active a parallel list of 0/1 flags (SCOTDO's final,
#  post-elimination meaning: "include in the area calculation").
#
#  fallback=1 means the annealed HOLE radius was already below the probe
#  radius - concal.f never runs the flood fill in that case, it just hands
#  the HOLE radius back (see the fidelity notes for why that IS the real
#  Requiv convention, not an approximation of one).

proc hole::concal {centre cvect perpve perpvn n xs ys zs rs endrad probe grid} {
    lassign $centre cx cy cz
    lassign $perpve ex ey ez
    lassign $perpvn nx ny nz
    set endrp3 [expr {$endrad + 3.0}]

    lassign [hole::holeen $cx $cy $cz $n $xs $ys $zs $rs] eng0
    set rad0 [expr {-$eng0}]
    if {$rad0 < $probe} {
        return [list {} {} 0 $rad0 1 0]
    }

    set scomax 10000
    set px [list $cx]; set py [list $cy]; set pz [list $cz]
    set prad [list $rad0]
    set pdo [list 1]
    set sconum 1
    set escaped 0
    set sconxt 0

    while {1} {
        set r [lindex $prad $sconxt]
        set do [lindex $pdo $sconxt]
        if {$r < $endrp3 && $do} {
            lset pdo $sconxt 0
            set bx [lindex $px $sconxt]; set by [lindex $py $sconxt]; set bz [lindex $pz $sconxt]
            for {set ncount -1} {$ncount <= 1 && !$escaped} {incr ncount} {
                for {set ecount -1} {$ecount <= 1 && !$escaped} {incr ecount} {
                    set tx [expr {$bx + $ecount*$grid*$ex + $ncount*$grid*$nx}]
                    set ty [expr {$by + $ecount*$grid*$ey + $ncount*$grid*$ny}]
                    set tz [expr {$bz + $ecount*$grid*$ez + $ncount*$grid*$nz}]
                    set dup 0
                    for {set s 0} {$s < $sconum} {incr s} {
                        set d2 [hole::dist2_3d_buggy $tx $ty $tz \
                                    [lindex $px $s] [lindex $py $s] [lindex $pz $s]]
                        if {$d2 < $grid/1000.0} { set dup 1; break }
                    }
                    if {$dup} { continue }
                    lassign [hole::holeen $tx $ty $tz $n $xs $ys $zs $rs] eng
                    set newrad [expr {-$eng}]
                    if {$newrad >= $probe} {
                        if {$sconum+1 > $scomax} {
                            error "hole::concal: SCOMAX (10000) reached - grid too fine or pore too open"
                        }
                        lappend px $tx; lappend py $ty; lappend pz $tz
                        lappend prad $newrad; lappend pdo 1
                        incr sconum
                        if {$newrad > $endrp3} { set escaped 1 }
                    } else {
                        # "spike" refinement: find where along bx,by,bz -> tx,ty,tz
                        # the pore radius equals the probe radius, by a fixed-point
                        # iteration on the excess (hodotu.f-style, per concal.f).
                        set vxu [expr {$tx-$bx}]; set vyu [expr {$ty-$by}]; set vzu [expr {$tz-$bz}]
                        set vn [expr {sqrt($vxu*$vxu + $vyu*$vyu + $vzu*$vzu)}]
                        set vxu [expr {$vxu/$vn}]; set vyu [expr {$vyu/$vn}]; set vzu [expr {$vzu/$vn}]
                        set delta [expr {0.25*$grid}]
                        set excess 0.0
                        set sx $tx; set sy $ty; set sz $tz
                        for {set cc 0} {$cc < 100} {incr cc} {
                            set sx [expr {$bx + $delta*$vxu}]
                            set sy [expr {$by + $delta*$vyu}]
                            set sz [expr {$bz + $delta*$vzu}]
                            lassign [hole::holeen $sx $sy $sz $n $xs $ys $zs $rs] seng
                            set srad [expr {-$seng}]
                            set excess [expr {$srad - ($probe+0.0001)}]
                            if {abs($excess) < 0.0005} { break }
                            set delta [expr {$delta + $excess}]
                        }
                        if {abs($excess) < 0.0005 && $delta > 0.0} {
                            set dup2 0
                            for {set s 0} {$s < $sconum} {incr s} {
                                set d2 [hole::dist2_3d_buggy $sx $sy $sz \
                                            [lindex $px $s] [lindex $py $s] [lindex $pz $s]]
                                if {$d2 < 0.09} { set dup2 1; break }
                            }
                            if {!$dup2} {
                                if {$sconum+1 > $scomax} {
                                    error "hole::concal: SCOMAX (10000) reached - grid too fine or pore too open"
                                }
                                lappend px $sx; lappend py $sy; lappend pz $sz
                                lappend prad [expr {$excess + $probe + 0.0001}]
                                lappend pdo 0
                                incr sconum
                            }
                        }
                    }
                }
            }
        }
        if {$escaped} break
        # highest-radius still-to-check point propagates first (concal.f 701)
        set maxrad -1e10; set nxt -1
        for {set s 0} {$s < $sconum} {incr s} {
            if {[lindex $pdo $s]} {
                set r2 [lindex $prad $s]
                if {$r2 > $maxrad} { set maxrad $r2; set nxt $s }
            }
        }
        if {$nxt < 0} break
        set sconxt $nxt
    }

    # SCOTDO's meaning flips here: was "needs checking", now "include in area calc"
    set active {}
    for {set s 0} {$s < $sconum} {incr s} { lappend active 1 }

    # elimination pass 1: points inside the escaped end sphere
    set lastIdx [expr {$sconum-1}]
    if {[lindex $prad $lastIdx] > $endrp3} {
        set delim2 [expr {pow([lindex $prad $lastIdx],2)}]
        set lx [lindex $px $lastIdx]; set ly [lindex $py $lastIdx]; set lz [lindex $pz $lastIdx]
        for {set s 0} {$s < $lastIdx} {incr s} {
            set d2 [hole::dist2_3d_buggy $lx $ly $lz [lindex $px $s] [lindex $py $s] [lindex $pz $s]]
            if {$d2 < $delim2} { lset active $s 0 }
        }
    }
    # elimination pass 2: points well within the original HOLE sphere
    if {[lindex $prad 0] > 2.5*$probe} {
        set d [expr {[lindex $prad 0] - 1.5*$probe}]
        set delim2 [expr {$d*$d}]
        set fx [lindex $px 0]; set fy [lindex $py 0]; set fz [lindex $pz 0]
        for {set s 1} {$s < $sconum} {incr s} {
            if {![lindex $active $s]} continue
            set d2 [hole::dist2_3d_buggy $fx $fy $fz [lindex $px $s] [lindex $py $s] [lindex $pz $s]]
            if {$d2 < $delim2} { lset active $s 0 }
        }
    }

    set points {}
    foreach x $px y $py z $pz r $prad { lappend points [list $x $y $z $r] }
    return [list $points $active $sconum $rad0 0 $escaped]
}

# ==============================================================================
#  COAREA - black/grey/white adaptive area estimator                (coarea.f)
# ==============================================================================
#
#  Returns a dict: requiv, sconum, ncycle, area, areag, escaped (1e6
#  sentinel), nocircles (0.0 sentinel), hitcap (SNMAX reached, see notes).

proc hole::coarea {centre perpve perpvn cvect sample endrad points active} {
    lassign $centre cx cy cz
    lassign $cvect vx vy vz
    lassign $perpve ex ey ez
    lassign $perpvn nx ny nz
    set n [llength $points]
    set pi [expr {2.0*asin(1.0)}]
    set root2 [expr {sqrt(2.0)}]

    set circrad {}; set ecoord {}; set ncoord {}
    foreach p $points { lappend circrad 0.0; lappend ecoord 0.0; lappend ncoord 0.0 }
    set act [lrange $active 0 end]

    for {set s 0} {$s < $n} {incr s} {
        if {![lindex $act $s]} continue
        lassign [lindex $points $s] px py pz prad
        set tx [expr {$px-$cx}]; set ty [expr {$py-$cy}]; set tz [expr {$pz-$cz}]
        set relc [expr {$tx*$vx + $ty*$vy + $tz*$vz}]
        if {abs($relc) < 1e-9} {
            set cr $prad
        } elseif {abs($relc) < $prad && abs($relc) < 0.9*$sample} {
            set cr [expr {sqrt($prad*$prad - $relc*$relc)}]
        } else {
            set cr -1e10
            lset act $s 0
        }
        lset circrad $s $cr
        if {$cr > $endrad} {
            # $act at this point is PARTIAL (this loop return matches
            # coarea.f's own early RETURN mid-loop) - concal.f's write block
            # (concal.f:535-556) still runs on an escaped slice using
            # whatever SCOTDO state resulted, so $act must be exposed here
            # too, not just on normal completion - see hole::connolly's own
            # header and vmdhole/hole_tcl/connolly.tcl's write_connolly_sph caller.
            return [dict create requiv 1.0e6 sconum $n ncycle 0 area 0.0 areag 0.0 \
                        escaped 1 nocircles 0 hitcap 0 active $act]
        }
    }

    set emax -1e10; set emin 1e10; set nmax -1e10; set nmin 1e10
    for {set s 0} {$s < $n} {incr s} {
        if {![lindex $act $s] || [lindex $circrad $s] >= $endrad} continue
        lassign [lindex $points $s] px py pz prad
        set tx [expr {$px-$cx}]; set ty [expr {$py-$cy}]; set tz [expr {$pz-$cz}]
        set e [expr {$tx*$ex + $ty*$ey + $tz*$ez}]
        set nn [expr {$tx*$nx + $ty*$ny + $tz*$nz}]
        lset ecoord $s $e; lset ncoord $s $nn
        set cr [lindex $circrad $s]
        if {$e+$cr > $emax} { set emax [expr {$e+$cr}] }
        if {$e-$cr < $emin} { set emin [expr {$e-$cr}] }
        if {$nn+$cr > $nmax} { set nmax [expr {$nn+$cr}] }
        if {$nn-$cr < $nmin} { set nmin [expr {$nn-$cr}] }
    }
    if {abs($emax+1e10) < 0.001} {
        return [dict create requiv 0.0 sconum $n ncycle 0 area 0.0 areag 0.0 \
                    escaped 0 nocircles 1 hitcap 0 active $act]
    }

    # coarse pass: fixed 0.25 Angs cells (independent of CONNR(2))
    set size1 0.25
    set area 0.0
    set greylist {}
    set tripE [expr {int(($emax-$emin+$size1)/$size1)}]
    if {$tripE < 0} { set tripE 0 }
    set edcnt $emin
    for {set ie 0} {$ie < $tripE} {incr ie} {
        set tripN [expr {int(($nmax-$nmin+$size1)/$size1)}]
        if {$tripN < 0} { set tripN 0 }
        set ndcnt $nmin
        for {set jn 0} {$jn < $tripN} {incr jn} {
            set cenE [expr {$edcnt+0.5*$size1}]
            set cenN [expr {$ndcnt+0.5*$size1}]
            set lblack 0; set isblack 0
            for {set s 0} {$s < $n} {incr s} {
                if {![lindex $act $s] || [lindex $circrad $s] >= $endrad} continue
                set cr [lindex $circrad $s]
                set de [expr {$cenE-[lindex $ecoord $s]}]
                set dn [expr {$cenN-[lindex $ncoord $s]}]
                set dist2 [expr {$de*$de+$dn*$dn}]
                set w [expr {$cr-0.5*$root2*$size1}]
                if {$dist2 < $w*$w} { set area [expr {$area+$size1*$size1}]; set isblack 1; break }
                set w [expr {$cr+0.5*$root2*$size1}]
                if {$dist2 < $w*$w} { set lblack 1 }
            }
            if {!$isblack && $lblack} { lappend greylist [list $edcnt $ndcnt] }
            set ndcnt [expr {$ndcnt+$size1}]
        }
        set edcnt [expr {$edcnt+$size1}]
    }

    set areag [expr {double([llength $greylist])*$size1*$size1}]
    set requiv [expr {sqrt(($area+0.5*$areag)/$pi)}]

    # refinement: subdivide each grey cell, classify its 8 boundary points
    set snmax 90000
    set curSize $size1
    set curList $greylist
    set ncycle 0
    set hitcap 0
    while {1} {
        incr ncycle
        set newSize [expr {0.5*$curSize}]
        set newCount 0
        set newList {}
        foreach cell $curList {
            lassign $cell bx by
            set edgeE [list $bx [expr {$bx+0.5*$curSize}] [expr {$bx+$curSize}] \
                            $bx [expr {$bx+$curSize}] \
                            $bx [expr {$bx+0.5*$curSize}] [expr {$bx+$curSize}]]
            set edgeN [list $by $by $by \
                            [expr {$by+0.5*$curSize}] [expr {$by+0.5*$curSize}] \
                            [expr {$by+$curSize}] [expr {$by+$curSize}] [expr {$by+$curSize}]]
            set edblk {0 0 0 0 0 0 0 0}
            set cenE [expr {$bx+0.5*$curSize}]
            set cenN [expr {$by+0.5*$curSize}]
            for {set s 0} {$s < $n} {incr s} {
                if {![lindex $act $s] || [lindex $circrad $s] >= $endrad} continue
                set cr [lindex $circrad $s]
                set de [expr {$cenE-[lindex $ecoord $s]}]
                set dn [expr {$cenN-[lindex $ncoord $s]}]
                set dist2 [expr {$de*$de+$dn*$dn}]
                set w [expr {$cr+0.5*$root2*$curSize}]
                if {$dist2 > $w*$w} continue
                for {set ei 0} {$ei < 8} {incr ei} {
                    set de2 [expr {[lindex $edgeE $ei]-[lindex $ecoord $s]}]
                    set dn2 [expr {[lindex $edgeN $ei]-[lindex $ncoord $s]}]
                    set d2 [expr {$de2*$de2+$dn2*$dn2}]
                    if {$d2 < $cr*$cr} { lset edblk $ei 1 }
                }
                set sumall 0
                foreach v $edblk { incr sumall $v }
                if {$sumall == 8} break
            }
            set sumall 0
            foreach v $edblk { incr sumall $v }
            if {$sumall == 8} {
                set area [expr {$area+$curSize*$curSize}]
            } elseif {$sumall == 0} {
                # white - nothing to do
            } else {
                if {$newCount+4 >= $snmax} {
                    incr newCount 4
                } else {
                    incr newCount 4
                    lappend newList [list $bx $by] [list [expr {$bx+$newSize}] $by] \
                        [list $bx [expr {$by+$newSize}]] [list [expr {$bx+$newSize}] [expr {$by+$newSize}]]
                }
            }
        }
        set areag [expr {double($newCount)*$newSize*$newSize}]
        set rlow [expr {sqrt($area/$pi)}]
        set oldreq $requiv
        set requiv [expr {sqrt(($area+0.5*$areag)/$pi)}]
        if {(abs($requiv-$oldreq) > 0.0005) || ($ncycle < 4)} {
            if {$newCount+4 <= $snmax} {
                set curSize $newSize
                set curList $newList
                continue
            } else {
                set hitcap 1
                break
            }
        } else {
            break
        }
    }
    return [dict create requiv $requiv sconum $n ncycle $ncycle area $area areag $areag \
                escaped 0 nocircles 0 hitcap $hitcap active $act]
}

# ==============================================================================
#  hole::connolly - the public entry point
# ==============================================================================
#
#  -centre {x y z}, -cvect/-perpve/-perpvn unit vectors (perpve/perpvn from
#  hole::calper, computed once outside the per-slice loop, exactly as HOLE
#  itself does), -n/-xs/-ys/-zs/-rs the atom arrays hole::holeen already
#  takes, -endrad, -probe (default 1.15, CONNR(1)'s default), -grid (default
#  0.7*probe, CONNR(2)'s default), -sample (signed, default 1.0 - see the
#  fidelity note above on why its sign is inert here).
#
#  Returns a dict: requiv, sconum, ncycle, fallback, escaped, nocircles,
#  hitcap, points, active - REQUIV alone is concal.f's REQUIV out-parameter;
#  the rest is exposed for testing against the reference binary's own
#  per-plane counts, AND for a .sph writer (write_connolly_sph below): points
#  is concal's own SCOXYZ/SCORAD (a list of {x y z rad}, index 0 always the
#  slice centre itself); active is SCOTDO's FINAL, post-coarea state (what
#  concal.f's own write loop actually filters on, concal.f:544 - "SCOTDO
#  (SCOUNT).OR.(SCORAD(SCOUNT).GT.ENDRAD)") - NOT concal's own
#  pre-elimination active, which coarea.f further narrows. On fallback,
#  points/active are concal's (empty - see hole::concal's own header: no
#  flood fill runs) since coarea never gets called.

proc hole::connolly {args} {
    array set o {
        -probe 1.15 -grid {} -endrad 22.0 -sample 1.0
    }
    array set o $args
    foreach req {-centre -cvect -perpve -perpvn -n -xs -ys -zs -rs} {
        if {![info exists o($req)]} { error "hole::connolly: missing $req" }
    }
    if {$o(-grid) eq ""} { set o(-grid) [expr {0.7*$o(-probe)}] }

    lassign [hole::concal $o(-centre) $o(-cvect) $o(-perpve) $o(-perpvn) \
                 $o(-n) $o(-xs) $o(-ys) $o(-zs) $o(-rs) \
                 $o(-endrad) $o(-probe) $o(-grid)] points active sconum rad0 fallback escaped
    if {$fallback} {
        return [dict create requiv $rad0 sconum 1 ncycle 0 fallback 1 escaped 0 nocircles 0 \
                    hitcap 0 points $points active $active]
    }
    set result [hole::coarea $o(-centre) $o(-perpve) $o(-perpvn) $o(-cvect) $o(-sample) \
                    $o(-endrad) $points $active]
    dict set result fallback 0
    dict set result points $points
    return $result
}

# ==============================================================================
#  write_connolly_sph - CONCAL's own .sph record            (concal.f:535-556)
# ==============================================================================
#
#  Unlike write_sph/write_capsule_sph (hole.tcl/capsule.tcl - one record per
#  slice), CONNOLLY writes EVERY flood-fill sphere point concal.f still
#  considers "in" for that slice (concal.f:544's filter:
#  `SCOTDO(SCOUNT) .OR. SCORAD(SCOUNT).GT.ENDRAD` - PEG_WRITEALL is not a card
#  this port implements, default .FALSE., so always omitted here) - the FIRST
#  point actually written for a slice gets that slice's real IREC, every
#  later point in the same slice gets -999 (concal.f: WRESNO starts as IREC,
#  reset to -999 right after the first successful WRITE - concal.f:536,553).
#  A slice whose annealed radius was already below the probe radius (concal.f
#  fallback - hole::connolly's own `fallback` key) writes NOTHING at all -
#  concal.f's write block is inside the ELSE of that check.
#
#  Both radius fields differ from write_sph's per-slice record: the first is
#  SCORAD(SCOUNT) (that POINT's own HOLE radius, not the slice's), the second
#  is WREQUIV - REQUIV clamped to 999.99 when it exceeds 1000 (concal.f:
#  538-539, an F6.2 field-width clamp, not a new sentinel) - the SAME
#  clamped value on every line of a given slice, computed once.
#
#  ADDEND (hole::addend, shared with write_sph/write_capsule_sph) and the
#  duplicate record-0 rewrite at the +ve/-ve transition are UNCHANGED from
#  the plain path even in CONNOLLY mode: holcal.f's own ADDEND call
#  (holcal.f:781) and record-0 duplicate (holcal.f:803-806) are not gated by
#  CONNR at all - both call plain HOLEEN/WPDBSP on the STORED search state
#  ($discovery, the same list hole::write_sph takes), never CONCAL. Confirmed
#  against the reference-build hole (see vmdhole/hole_tcl/README.md)'s own CONNOLLY+sphpdb output: the duplicate
#  record-0 line has STRRAD duplicated in both radius fields (matching
#  write_sph's plain form), NOT SCORAD/WREQUIV.
#
#  COST NOTE: re-runs hole::connolly per slice independently of any earlier
#  -csv/CONNOLLY pass hole::main already did (that pass only kept `requiv`,
#  not `points`/`active` - see hole::connolly's own header) - CONNOLLY is
#  deterministic (no RNG draws), so this is a wall-clock cost, not a
#  correctness risk, but it roughly DOUBLES the CONNOLLY portion of the run
#  when both -csv and -sph are requested together. Not optimised - see
#  vmdhole/hole_tcl/tests/sph_addend_test.sh for the measured cost on the fixture
#  this project's other CONNOLLY numbers use.
#
#  Byte-identical to the reference-build hole (see vmdhole/hole_tcl/README.md) on 1GRM (cpoint 0 0 4, cvect 0 0 1,
#  sample 0.5, endrad 8.0, raseed 1, conn 1.15 0.2):
#    printf 'sphpdb refn.sph\ncoord 1GRM.pdb\nradius simple.rad\n
#      cpoint 0 0 4\ncvect 0 0 1\nsample 0.5\nendrad 8.0\nraseed 1\n
#      conn 1.15 0.2\n' | the reference-build hole (see vmdhole/hole_tcl/README.md)
#    tclsh hole.tcl -pdb 1GRM.pdb -rad simple.rad -method connolly \
#        -cpoint "0 0 4" -cvect "0 0 1" -sample 0.5 -endrad 8.0 -seed 1 \
#        -probe 1.15 -grid 0.2 -sph tcln.sph
#    diff -q refn.sph tcln.sph      # clean

proc hole::_write_connolly_slice {fh row ncvect perpve perpvn n xs ys zs rs endrad probe grid sample} {
    lassign $row irec t x y z rad a1 a2
    set res [hole::connolly -centre [list $x $y $z] -cvect $ncvect \
                -perpve $perpve -perpvn $perpvn \
                -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
                -endrad $endrad -probe $probe -grid $grid -sample $sample]
    if {[dict get $res fallback]} { return }
    set requiv [dict get $res requiv]
    set wrequiv [expr {$requiv > 1000.0 ? 999.99 : $requiv}]
    set points [dict get $res points]
    set active [dict get $res active]
    set sconum [llength $points]
    set wresno $irec
    for {set s 0} {$s < $sconum} {incr s} {
        lassign [lindex $points $s] px py pz prad
        if {[lindex $active $s] || $prad > $endrad} {
            puts $fh [hole::_sph_atom_line $wresno $px $py $pz $prad $wrequiv]
            if {$prad > $endrad} { puts $fh "LAST-REC-END" }
            set wresno -999
        }
    }
}

proc hole::write_connolly_sph {discovery path args} {
    array set o {-probe 1.15 -grid {}}
    array set o $args
    foreach req {-n -xs -ys -zs -rs -cvect -sample -endrad} {
        if {![info exists o($req)]} { error "hole::write_connolly_sph: missing $req" }
    }
    if {$o(-grid) eq ""} { set o(-grid) [expr {0.7*$o(-probe)}] }
    set n $o(-n); set xs $o(-xs); set ys $o(-ys); set zs $o(-zs); set rs $o(-rs)
    set endrad $o(-endrad)
    set sample [expr {abs($o(-sample))}]
    lassign $o(-cvect) cvx cvy cvz
    set cvn [expr {sqrt($cvx*$cvx + $cvy*$cvy + $cvz*$cvz)}]
    if {$cvn < 1e-9} { error "hole::write_connolly_sph: CVECT is a zero vector" }
    set cvx [expr {$cvx/$cvn}]; set cvy [expr {$cvy/$cvn}]; set cvz [expr {$cvz/$cvn}]
    set ncvect [list $cvx $cvy $cvz]
    lassign [hole::calper $ncvect] perpve perpvn

    set pos {}; set neg {}
    foreach row $discovery {
        if {[lindex $row 0] >= 0} { lappend pos $row } else { lappend neg $row }
    }
    if {[llength $pos] == 0} {
        error "hole::write_connolly_sph: no +ve-direction slice was stored - nothing to write"
    }

    set fh [open $path w]

    foreach row $pos {
        hole::_write_connolly_slice $fh $row $ncvect $perpve $perpvn $n $xs $ys $zs $rs \
            $endrad $o(-probe) $o(-grid) $sample
    }
    if {[llength $pos] < 2} {
        close $fh
        error "hole::write_connolly_sph: only [llength $pos] +ve slice(s) stored - ADDEND's STRCEN(,STRNOP-1) reset needs at least 2; not reproduced for this degenerate case"
    }
    lassign [lindex $pos end-1] li lt lx ly lz lrad la1 la2
    foreach pt [hole::addend [list $lx $ly $lz] $sample $ncvect $perpve $perpvn \
                    $endrad $n $xs $ys $zs $rs] {
        lassign $pt ax ay az arad
        puts $fh [hole::_sph_atom_line -888 $ax $ay $az $arad 0.00]
        puts $fh "LAST-REC-END"
    }

    # holcal.f:803-806 - the duplicate record-0 rewrite is NOT gated by CONNR
    # either; it is the plain WPDBSP form (STRRAD duplicated), not CONCAL's -
    # see this proc's own header.
    lassign [lindex $pos 0] irec0 t0 x0 y0 z0 rad0 a10 a20
    puts $fh [hole::_sph_atom_line $irec0 $x0 $y0 $z0 $rad0 $rad0]

    foreach row $neg {
        hole::_write_connolly_slice $fh $row $ncvect $perpve $perpvn $n $xs $ys $zs $rs \
            $endrad $o(-probe) $o(-grid) $sample
    }
    if {[llength $neg] < 2} {
        close $fh
        error "hole::write_connolly_sph: only [llength $neg] -ve slice(s) stored - ADDEND's STRCEN(,-STRNON+1) reset needs at least 2; not reproduced for this degenerate case"
    }
    lassign [lindex $neg end-1] li lt lx ly lz lrad la1 la2
    foreach pt [hole::addend [list $lx $ly $lz] [expr {-$sample}] $ncvect $perpve $perpvn \
                    $endrad $n $xs $ys $zs $rs] {
        lassign $pt ax ay az arad
        puts $fh [hole::_sph_atom_line -888 $ax $ay $az $arad 0.00]
        puts $fh "LAST-REC-END"
    }

    close $fh
}
