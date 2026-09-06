#!/usr/bin/env tclsh
# ==============================================================================
#  sph_process.tcl - dot-surface generation from a .sph file (sph_process.f)
# ==============================================================================
#
#  Sourced by hole.tcl. Ports the "-sos -dotden D" path of sph_process.f: read
#  a .sph (hole::write_sph / real hole's SPHPDB output), place a template of
#  isotropically-distributed dots on every sphere, drop dots buried inside a
#  neighbouring sphere, and write the surviving dots as a .sos file for
#  hole::sos_triangle.
#
#  Only the -sos default path is ported: no -qpt/-spaghetti (line-drawing
#  output, never read by sos_triangle), no -colour, no -rover/-ignore_small,
#  no -peg, no CAPSULE (.sph QC1/QC2 records) and no -tri (direct triangle
#  output via SPHTRI/TRISPHERE - a completely different, disused code path,
#  see the note in sos_triangle.tcl). All of those are simply not reachable
#  through hole::sph_process's 3-argument signature; feeding them the
#  relevant .sph records raises rather than silently mishandling them.
#
#  Read directly against sph_process.f, sph_process_read.f, ptgen.f, sphqpu.f,
#  and for CAPSULE records sphqpc.f, hocapr.f, hocapd.f.

namespace eval hole {}

proc hole::_f32 {x} {
    # sphqpu.f's RVEC3 is declared REAL (Fortran default = 4-byte): the dot
    # position is truncated to single precision BEFORE the occlusion test and
    # before being written to the .sos file. Skipping this doesn't just cost
    # a few ULP: it can flip which dots survive the occlusion test at a
    # sphere-sphere boundary.
    binary scan [binary format f $x] f y
    return $y
}

proc hole::_ptgen {dotden} {
    # ptgen.f: dots isotropically spaced on a unit sphere as a stack of
    # constant-z rings. PI itself is single-precision in the Fortran (the
    # literal '0.' in ACOS(0.)*2. is a default REAL) then widened to DOUBLE
    # PRECISION on assignment - reproduced here so ring angles match exactly.
    set pi [hole::_f32 [expr {acos(0.0) * 2.0}]]
    set pts {}
    set half [expr {int($dotden) / 2}]
    for {set zc 0} {$zc <= $half} {incr zc} {
        # ZANG = REAL(ZCOUNT)/REAL(DOTDEN)*PI: left-to-right evaluation with
        # PI (the only DOUBLE PRECISION operand) appearing LAST means
        # ZCOUNT/DOTDEN is a single-precision division in the Fortran, THEN
        # widened for the multiply - unlike the X/Y angle below, where "2.*PI"
        # forces double arithmetic from the first operation, so ACOUNT/RNUMB
        # never gets rounded to single on its own. Verified bit-for-bit
        # against a standalone gfortran build of the real ptgen.f (all 152
        # dots at dotden=10, max diff 4e-18 - double-precision noise only).
        set zratio [hole::_f32 [expr {double($zc) / double($dotden)}]]
        set zang [expr {$zratio * $pi}]
        set rcirc [expr {sin($zang)}]
        set zcoord [expr {cos($zang)}]
        # +0.95 rounds up: Fortran REAL->INTEGER assignment truncates toward
        # zero, which for this always-positive value is int().
        set rnumb [expr {int(double(2 * $dotden) * $rcirc + 0.95)}]
        if {$rnumb == 0} {
            lappend pts [list 0.0 0.0 $zcoord]
            lappend pts [list 0.0 0.0 [expr {-$zcoord}]]
        } else {
            for {set ac 0} {$ac < $rnumb} {incr ac} {
                set xcoord [expr {$rcirc * cos(2.0 * $pi * double($ac) / double($rnumb))}]
                set ycoord [expr {$rcirc * sin(2.0 * $pi * double($ac) / double($rnumb))}]
                # Swap x/y every other ring so no great semicircle is left
                # dot-free (ptgen.f's own comment).
                if {$zc % 2 == 0} {
                    set t $xcoord; set xcoord $ycoord; set ycoord $t
                }
                lappend pts [list $xcoord $ycoord $zcoord]
                lappend pts [list $xcoord $ycoord [expr {-$zcoord}]]
            }
        }
    }
    return $pts
}

proc hole::_sph_read {sph_file} {
    # Returns a list of {x y z rad effr last x2 y2 z2} records, parsed at the
    # FIXED column positions sph_process_read.f reads (not whitespace-split -
    # the residue-number field can run into the coordinates with no space).
    # A sphere has x2 y2 z2 equal to x y z. A CAPSULE slice is a QC1 record
    # followed by its QC2 partner (wpdbsp.f): the two cap centres, the capsule
    # radius, and the equal-area radius sphqpc.f computes from them.
    set fh [open $sph_file r]
    set spheres {}
    set pend ""
    while {[gets $fh line] >= 0} {
        if {[string range $line 0 11] eq "LAST-REC-END"} {
            if {[llength $spheres]} {
                set last [lindex $spheres end]
                lset last 5 1
                lset spheres end $last
            }
            continue
        }
        if {[string range $line 0 3] ne "ATOM"} { continue }
        set tag [string range $line 10 21]
        if {$tag ni {"1  QSS SPH S" "1  QC1 SPH S" "1  QC2 SPH S"}} {
            error "hole::sph_process: unrecognized .sph ATOM record: $line"
        }
        set x [string trim [string range $line 30 37]]
        set y [string trim [string range $line 38 45]]
        set z [string trim [string range $line 46 53]]
        set rad [string trim [string range $line 54 59]]
        if {![string is double -strict $x] || ![string is double -strict $y] || \
                ![string is double -strict $z] || ![string is double -strict $rad]} {
            error "hole::sph_process: malformed .sph ATOM record: $line"
        }
        if {$tag eq "1  QC1 SPH S"} { set pend [list $x $y $z $rad]; continue }
        if {$tag eq "1  QC2 SPH S"} {
            if {$pend eq ""} { error "hole::sph_process: QC2 record without its QC1: $line" }
            lassign $pend x1 y1 z1 rad
            set len [expr {sqrt(($x1-$x)*($x1-$x) + ($y1-$y)*($y1-$y) + ($z1-$z)*($z1-$z))}]
            # sphqpc.f: EFFRAD = SQRT(SPRAD*(SPRAD + 0.63661977*length))
            set effr [expr {sqrt($rad*($rad + 0.63661977*$len))}]
            lappend spheres [list $x1 $y1 $z1 $rad $effr 0 $x $y $z]
            set pend ""
            continue
        }
        set effr ""
        if {[string length $line] >= 66} {
            set effr [string trim [string range $line 60 65]]
        }
        # Old .sph format has no effective-radius column, or it reads zero;
        # sph_process_read.f falls back to the real radius in both cases.
        if {$effr eq "" || ![string is double -strict $effr] || abs($effr) < 1e-6} {
            set effr $rad
        }
        lappend spheres [list $x $y $z $rad $effr 0 $x $y $z]
    }
    close $fh
    return $spheres
}

proc hole::_write_sos_line {fh vals} {
    set line ""
    foreach v $vals { append line [format "%12.5f" $v] }
    puts $fh $line
}

proc hole::sph_process {sph_file sos_file dotden {color 0}} {
    # color=1 is sph_process.f's -colour: three radius-banded passes instead of
    # one, each with its own quanta color header, which is what makes the
    # surface read as red/green/blue by pore radius. Bands are RCUT(1)=1.15 and
    # RCUT(2)=2.30 (sph_process.f 405-416); the -Rlow/-Rmid overrides are not
    # exposed because the plugin never passes them.
    if {![string is integer -strict $dotden] || $dotden < 1 || $dotden > 100} {
        error "hole::sph_process: -dotden must be an integer between 1 and 100"
    }
    set spheres [hole::_sph_read $sph_file]
    set n [llength $spheres]
    set templ [hole::_ptgen $dotden]
    set ptno [llength $templ]
    if {$ptno == 0} { error "hole::sph_process: dot template is empty for dotden=$dotden" }

    # Header per pass (SPHCHC) and the SPEFFR band it selects: a sphere is
    # emitted on pass p when RCUT(p-1) < SPEFFR <= RCUT(p). Without -colour
    # there is one pass and RCUT = {-1, 9999}, so every SPEFFR>0 sphere lands in
    # it - the same set the single-pass version emitted.
    if {$color} {
        set cuts {-1.0 1.15 2.30 9999.0}
        set hdrs {{1.0 3.0 -55.0 16.0 0.0 0.0 0.0}
                  {1.0 7.0 -55.0 17.0 0.0 0.0 0.0}
                  {1.0 2.0 -55.0 18.0 0.0 0.0 0.0}}
    } else {
        set cuts {-1.0 9999.0}
        set hdrs {{1.0 7.0 -55.0 17.0 0.0 0.0 0.0}}
    }

    set fh [open $sos_file w]
    set np [llength $hdrs]
    if {[hole::_sph_has_capsule $spheres]} {
        # sph_process.f hands a file with capsule records to sphqpc.f, whose
        # pass loop calls SPHCHC for every pass: under -sos the LAST pass takes
        # the colour -1 end-cap header, so with -colour the high-radius band
        # goes out under it and there is no separate end-cap pass. Reproduced
        # as is: this is what HOLE writes.
        for {set p 0} {$p < $np} {incr p} {
            if {$np == 1} { set hdr [lindex $hdrs 0] } \
            elseif {$p == $np - 1} { set hdr {1.0 -1.0 -1.0 -1.0 0.0 0.0 0.0} } \
            else { set hdr [lindex $hdrs $p] }
            hole::_write_sos_line $fh $hdr
            set lo [lindex $cuts $p]
            set hi [lindex $cuts [expr {$p + 1}]]
            for {set i 0} {$i < $n} {incr i} {
                lassign [lindex $spheres $i] cx cy cz rad effr last
                if {$last} { continue }
                if {$rad <= 0} { continue }
                if {$effr <= $lo || $effr > $hi} { continue }
                hole::_emit_capsule_dots $fh $spheres $i $templ $ptno
            }
        }
        close $fh
        return
    }
    for {set p 0} {$p < $np} {incr p} {
        hole::_write_sos_line $fh [lindex $hdrs $p]
        set lo [lindex $cuts $p]
        set hi [lindex $cuts [expr {$p + 1}]]
        for {set i 0} {$i < $n} {incr i} {
            lassign [lindex $spheres $i] cx cy cz rad effr last
            if {$last} { continue }
            if {$effr <= 0} { continue }
            if {$effr <= $lo || $effr > $hi} { continue }
            hole::_emit_sphere_dots $fh $spheres $i $templ $ptno
        }
    }

    # SOSOUT's extra end-cap pass, quanta color -1 - only the LAST-REC-END
    # spheres, so sos_triangle can cull the triangles that touch them for a
    # sharp edge (sos_triangle.c cull_triangles).
    hole::_write_sos_line $fh {1.0 -1.0 -1.0 -1.0 0.0 0.0 0.0}
    for {set i 0} {$i < $n} {incr i} {
        lassign [lindex $spheres $i] cx cy cz rad effr last
        if {!$last} { continue }
        hole::_emit_sphere_dots $fh $spheres $i $templ $ptno
    }

    close $fh
}

proc hole::_emit_sphere_dots {fh spheres i templ ptno} {
    # sphqpu.f's per-sphere axis swap: SCOUNT is the sphere's 1-based
    # position in the FULL sphere list (all n records, both passes share the
    # same indexing), used to decorrelate dot patterns between neighbouring
    # spheres so their template alignments don't line up into visible seams.
    set scount [expr {$i + 1}]
    set is1 [expr {($scount + 2) % 3}]
    set is2 [expr {$scount % 3}]
    set is3 [expr {($scount + 1) % 3}]

    lassign [lindex $spheres $i] cx cy cz rad effr last

    for {set d 0} {$d < $ptno} {incr d} {
        set p [lindex $templ $d]
        set nx [lindex $p $is1]
        set ny [lindex $p $is2]
        set nz [lindex $p $is3]
        set rx [hole::_f32 [expr {$cx + $rad * $nx}]]
        set ry [hole::_f32 [expr {$cy + $rad * $ny}]]
        set rz [hole::_f32 [expr {$cz + $rad * $nz}]]

        # Occlusion test: reject the dot if it lies inside ANY other sphere.
        # sphqpu.f walks outward from SCOUNT in both directions for the same
        # net effect (first-hit rejects); a full scan over the other n-1
        # spheres gives an identical accept/reject decision.
        set buried 0
        for {set j 0} {$j < [llength $spheres]} {incr j} {
            if {$j == $i} { continue }
            lassign [lindex $spheres $j] ox oy oz orad
            set dx [expr {$rx - $ox}]
            set dy [expr {$ry - $oy}]
            set dz [expr {$rz - $oz}]
            if {[expr {$dx*$dx + $dy*$dy + $dz*$dz}] < [expr {$orad*$orad}]} {
                set buried 1
                break
            }
        }
        if {$buried} { continue }
        hole::_write_sos_line $fh [list 4.0 $rx $ry $rz $nx $ny $nz]
    }
}

proc hole::_sph_has_capsule {spheres} {
    foreach s $spheres {
        lassign $s x y z rad effr last x2 y2 z2
        if {$x2 != $x || $y2 != $y || $z2 != $z} { return 1 }
    }
    return 0
}

proc hole::_capsule_d2 {px py pz x1 y1 z1 x2 y2 z2} {
    # hocapd.f: squared distance from a point to the segment between the two
    # cap centres.
    set ux [expr {$x2-$x1}]; set uy [expr {$y2-$y1}]; set uz [expr {$z2-$z1}]
    set len [expr {sqrt($ux*$ux+$uy*$uy+$uz*$uz)}]
    if {$len > 0.0} { set ux [expr {$ux/$len}]; set uy [expr {$uy/$len}]; set uz [expr {$uz/$len}] }
    set vx [expr {$px-$x1}]; set vy [expr {$py-$y1}]; set vz [expr {$pz-$z1}]
    set plen [expr {$ux*$vx+$uy*$vy+$uz*$vz}]
    if {$plen <= 0.0} { return [expr {$vx*$vx+$vy*$vy+$vz*$vz}] }
    if {$plen >= $len} {
        return [expr {($px-$x2)*($px-$x2)+($py-$y2)*($py-$y2)+($pz-$z2)*($pz-$z2)}]
    }
    set vx [expr {$vx-$plen*$ux}]; set vy [expr {$vy-$plen*$uy}]; set vz [expr {$vz-$plen*$uz}]
    return [expr {$vx*$vx+$vy*$vy+$vz*$vz}]
}

proc hole::_emit_capsule_dots {fh spheres i templ ptno} {
    # sphqpc.f: a dot in each template direction from the capsule's midpoint,
    # at the radius of the capsule surface in that direction (hocapr.f), kept
    # unless it lies inside any other record (hocapd.f). No per-record axis
    # swap and the template direction is written as the normal, as there.
    set pi [expr {2.0*acos(0.0)}]
    lassign [lindex $spheres $i] x1 y1 z1 rad effr last x2 y2 z2
    set cx [expr {0.5*($x1+$x2)}]; set cy [expr {0.5*($y1+$y2)}]; set cz [expr {0.5*($z1+$z2)}]
    set vx [expr {$x1-$x2}]; set vy [expr {$y1-$y2}]; set vz [expr {$z1-$z2}]
    set caplen [expr {sqrt($vx*$vx+$vy*$vy+$vz*$vz)}]
    if {$caplen > 0.0} { set vx [expr {$vx/$caplen}]; set vy [expr {$vy/$caplen}]; set vz [expr {$vz/$caplen}] }
    set half [expr {0.5*$caplen}]
    set alphac [expr {atan2($rad, $half)}]
    set n [llength $spheres]
    for {set d 0} {$d < $ptno} {incr d} {
        lassign [lindex $templ $d] nx ny nz
        # hocapr.f: the capsule's radius along this direction (a sphere when
        # the two centres coincide)
        set cosa [expr {$vx*$nx+$vy*$ny+$vz*$nz}]
        if {$cosa > 1.0} { set cosa 1.0 } elseif {$cosa < -1.0} { set cosa -1.0 }
        set alpha [expr {acos($cosa)}]
        if {$half < 1e-9} {
            set erad $rad
        } elseif {$alpha > $alphac && $alpha < $pi-$alphac} {
            set erad [expr {$rad/sin($alpha)}]
        } else {
            if {$alpha > 0.5*$pi} { set alpha [expr {$pi-$alpha}] }
            set hc [expr {$half*cos($alpha)}]
            set erad [expr {$hc + sqrt($hc*$hc - ($half*$half - $rad*$rad))}]
        }
        set rx [hole::_f32 [expr {$cx + $erad*$nx}]]
        set ry [hole::_f32 [expr {$cy + $erad*$ny}]]
        set rz [hole::_f32 [expr {$cz + $erad*$nz}]]
        set buried 0
        for {set j 0} {$j < $n} {incr j} {
            if {$j == $i} { continue }
            lassign [lindex $spheres $j] ox oy oz orad oeff olast ox2 oy2 oz2
            if {[hole::_capsule_d2 $rx $ry $rz $ox $oy $oz $ox2 $oy2 $oz2] < $orad*$orad} {
                set buried 1
                break
            }
        }
        if {$buried} { continue }
        hole::_write_sos_line $fh [list 4.0 $rx $ry $rz $nx $ny $nz]
    }
}
