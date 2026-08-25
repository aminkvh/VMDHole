#!/usr/bin/env tclsh
# ==============================================================================
#  sos_triangle.tcl - dot-surface triangulation (sos_triangle.c)
# ==============================================================================
#
#  Sourced by hole.tcl. Ports sos_triangle.c's default pipeline
#  (read_cord -> cull_coords -> the advancing-front triangulator ->
#  cull_triangles -> reorder_triangle -> vmd_out) at the -s (smooth,
#  "draw trinorm") setting, which is the ONLY mode the plugin ever invokes
#  (vmdhole.tcl's mesh job always runs "sos_triangle -s") and the one the
#  task's own verification command uses.
#
#  CORRECTION TO THE BRIEF: sos_triangle is sos_triangle.c, a ~1900-line
#  step-by-step ("gift wrapping" / advancing-front) triangulator written by
#  Guy Coates - it has NO connection to trisphere.f. That file's TRISPHERE
#  subroutine (10 hardcoded per-density unit-sphere tessellations) is reached
#  only through SPHTRI, which sph_process.f calls solely for its unused -tri
#  direct-triangle-output option (confirmed by `grep -rln 'call SPHTRI'`: only
#  sph_process.f, and only behind TRIOUT, which hole::sph_process never sets).
#  The default -sos/sos_triangle pipeline this file ports never touches those
#  55,449 lines, so there was nothing there to transcribe or reimplement.
#
#  All mutable state lives in namespace variables (not proc locals passed by
#  upvar): the triangulator is naturally recursive, and an upvar chain only
#  reaches the immediate caller's frame, not the one true top-level store -
#  wrong as soon as a proc calls itself. Reset at the start of every
#  hole::sos_triangle call so back-to-back calls don't see stale state.
#
#  Read directly against sos_triangle.c.

namespace eval hole {
    # cos(102.0) - VOODOO_ANGLE is #define'd as 102.0 but fed straight into
    # C's cos(), which takes RADIANS: this is +0.1016, not cos(102 deg) =
    # -0.208. Reproduced literally, quirk and all.
    variable _voodoo_cos [expr {cos(102.0)}]
    # sos_triangle.c cull_triangles(): default max_vertex_length=5.0 Angs (no
    # -X support - the plugin never passes it).
    variable _max_vertex_len2 25.0

    # dot cloud (post cull_coords), 0-based, parallel arrays for speed
    variable _dx; variable _dy; variable _dz
    variable _dcol; variable _dnx; variable _dny; variable _dnz
    variable _ndots 0

    # base_line tree (calc_tri's nodes), 0-based ids
    variable _na; variable _nb; variable _nc; variable _nz
    variable _nb1a; variable _nb2a
    variable _next_id 0

    # edge_list: key "min,max" -> list of {owner_id owner_field}, append-only
    variable _edges

    variable _tri_list {}
    variable _used
}

proc hole::_sos_read {sos_file} {
    # read_cord(): each line is 7 floats; marker 1.0 = color change (2nd
    # value = quanta color), marker 4.0 = a dot (x y z nx ny nz). sos_triangle
    # parses with C's "%f" (32-bit float) - the text was already float-rounded
    # by sph_process's own '(7F12.5)' output, so re-parsing at float adds no
    # further loss for realistic HOLE coordinate ranges; done anyway for
    # faithfulness and so an externally-supplied .sos gets the same treatment.
    set fh [open $sos_file r]
    set color 1.0
    set raw {}
    while {[gets $fh line] >= 0} {
        set t [string trim $line]
        if {$t eq ""} { continue }
        set n [llength $t]
        if {$n < 4} { continue }
        set v0 [hole::_f32 [lindex $t 0]]
        if {$v0 == 1.0} {
            set color [hole::_f32 [lindex $t 1]]
            continue
        }
        if {$v0 == 4.0} {
            set x [hole::_f32 [lindex $t 1]]
            set y [hole::_f32 [lindex $t 2]]
            set z [hole::_f32 [lindex $t 3]]
            set nx [expr {$n > 4 ? [hole::_f32 [lindex $t 4]] : 0.0}]
            set ny [expr {$n > 5 ? [hole::_f32 [lindex $t 5]] : 0.0}]
            set nz [expr {$n > 6 ? [hole::_f32 [lindex $t 6]] : 0.0}]
            lappend raw [list $x $y $z $color $nx $ny $nz]
        }
    }
    close $fh
    return $raw
}

proc hole::_cull_coords {raw} {
    # cull_coords(): sequential, first-accepted-wins, per-axis L-infinity box
    # test at < 1e-3 (NOT a rounding-bucket test - non-transitive, so this
    # must stay a real "does an accepted point already here match" query).
    # Bucketed on a >=1e-3 grid, searching the 3x3x3 neighbourhood, is an
    # exact reformulation: any pair within the 1e-3 box on all 3 axes must
    # fall in the same or an adjacent cell, and cull_coords only cares
    # whether ANY match exists (the action is the same regardless of which
    # one), so this changes nothing about the result, only its cost.
    set cell 2.0e-3
    set tol 1.0e-3
    array set grid {}
    set dots {}
    set n 0
    foreach r $raw {
        lassign $r x y z color nx ny nz
        set cix [expr {int(floor($x / $cell))}]
        set ciy [expr {int(floor($y / $cell))}]
        set ciz [expr {int(floor($z / $cell))}]
        set dup 0
        for {set dxc -1} {$dxc <= 1 && !$dup} {incr dxc} {
            for {set dyc -1} {$dyc <= 1 && !$dup} {incr dyc} {
                for {set dzc -1} {$dzc <= 1 && !$dup} {incr dzc} {
                    set key "[expr {$cix+$dxc}],[expr {$ciy+$dyc}],[expr {$ciz+$dzc}]"
                    if {![info exists grid($key)]} { continue }
                    foreach idx $grid($key) {
                        lassign [lindex $dots $idx] ox oy oz
                        if {abs($x-$ox) < $tol && abs($y-$oy) < $tol && abs($z-$oz) < $tol} {
                            set dup 1
                            break
                        }
                    }
                }
            }
        }
        if {$dup} { continue }
        lappend dots [list $x $y $z $color $nx $ny $nz]
        lappend grid($cix,$ciy,$ciz) $n
        incr n
    }
    return $dots
}

proc hole::_ekey {x1 x2} { expr {$x1 < $x2 ? "$x1,$x2" : "$x2,$x1"} }

proc hole::_add_edge {x1 x2 owner_id owner_field} {
    variable _edges
    lappend _edges([hole::_ekey $x1 $x2]) [list $owner_id $owner_field]
}

proc hole::_destroy {e1 e2 caller_id caller_field} {
    variable _edges; variable _nb1a; variable _nb2a
    set key [hole::_ekey $e1 $e2]
    if {![info exists _edges($key)] || ![llength $_edges($key)]} { return }
    if {$caller_field eq "b1a"} { set _nb1a($caller_id) 0 } else { set _nb2a($caller_id) 0 }
    lassign [lindex $_edges($key) 0] owner_id owner_field
    # Root-baseline edges are registered with owner_id "" (sos_triangle.c
    # points these at a dead stack local, &dummy, that nothing ever reads
    # again - modelled as a deliberate no-op, not a bug).
    if {$owner_id eq ""} { return }
    if {$owner_field eq "b1a"} { set _nb1a($owner_id) 0 } else { set _nb2a($owner_id) 0 }
}

proc hole::_gen_triangle {a b c} {
    variable _dcol; variable _tri_list; variable _used
    if {$_dcol($a) == $_dcol($b) || $_dcol($a) == $_dcol($c)} {
        set col $_dcol($a)
    } elseif {$_dcol($b) == $_dcol($c)} {
        set col $_dcol($b)
    } else {
        set col $_dcol($a)
    }
    lappend _tri_list [list $a $b $c $col]
    set _used($a) 1; set _used($b) 1; set _used($c) 1
    if {[llength $_tri_list] > 30000} {
        error "hole::sos_triangle: MAX_COORD (30000) polygons exceeded"
    }
}

proc hole::_neighbour {a b z} {
    variable _dx; variable _dy; variable _dz; variable _ndots; variable _voodoo_cos
    set ax $_dx($a); set ay $_dy($a); set az $_dz($a)
    set bx $_dx($b); set by $_dy($b); set bz $_dz($b)
    set zx $_dx($z); set zy $_dy($z); set zz $_dz($z)
    set pmx [expr {0.5*($ax+$bx)}]
    set pmy [expr {0.5*($ay+$by)}]
    set pmz [expr {0.5*($az+$bz)}]
    set zmx [expr {$pmx-$zx}]; set zmy [expr {$pmy-$zy}]; set zmz [expr {$pmz-$zz}]
    set magzm [expr {sqrt($zmx*$zmx+$zmy*$zmy+$zmz*$zmz)}]
    set basedist2 [expr {($ax-$bx)*($ax-$bx)+($ay-$by)*($ay-$by)+($az-$bz)*($az-$bz)}]
    set min_ang 1.0
    set c -1
    for {set i 0} {$i < $_ndots} {incr i} {
        if {$i == $a || $i == $b || $i == $z} { continue }
        set ix $_dx($i); set iy $_dy($i); set iz $_dz($i)
        set mcx [expr {$ix-$pmx}]; set mcy [expr {$iy-$pmy}]; set mcz [expr {$iz-$pmz}]
        set magmc [expr {sqrt($mcx*$mcx+$mcy*$mcy+$mcz*$mcz)}]
        if {$magmc == 0 || $magzm == 0} { continue }
        set dotp [expr {$mcx*$zmx+$mcy*$zmy+$mcz*$zmz}]
        set angle2 [expr {$dotp/($magmc*$magzm)}]
        if {$angle2 <= $_voodoo_cos} { continue }
        set vax [expr {$ix-$ax}]; set vay [expr {$iy-$ay}]; set vaz [expr {$iz-$az}]
        set vbx [expr {$ix-$bx}]; set vby [expr {$iy-$by}]; set vbz [expr {$iz-$bz}]
        set maga [expr {sqrt($vax*$vax+$vay*$vay+$vaz*$vaz)}]
        set magb [expr {sqrt($vbx*$vbx+$vby*$vby+$vbz*$vbz)}]
        if {$maga == 0 || $magb == 0} { continue }
        set dp2 [expr {$vax*$vbx+$vay*$vby+$vaz*$vbz}]
        set angle [expr {$dp2/($maga*$magb)}]
        if {$angle < $min_ang} {
            # Kludge in the original: reject a long thin sliver even if it
            # has the best angle - measured from 'a' only, and does NOT
            # update min_ang, so a later, better-angled point can still win.
            set check_dist [expr {$vax*$vax+$vay*$vay+$vaz*$vaz}]
            if {$check_dist > 9.0*$basedist2} { continue }
            set min_ang $angle
            set c $i
        }
    }
    return $c
}

proc hole::_new_node {a b z} {
    variable _na; variable _nb; variable _nc; variable _nz
    variable _nb1a; variable _nb2a; variable _next_id
    set id $_next_id; incr _next_id
    set _na($id) $a; set _nb($id) $b; set _nz($id) $z; set _nc($id) -1
    set _nb1a($id) 1; set _nb2a($id) 1
    return $id
}

proc hole::_calc_tri {id} {
    variable _na; variable _nb; variable _nc; variable _nz
    variable _nb1a; variable _nb2a
    set a $_na($id); set b $_nb($id); set z $_nz($id)
    set c [hole::_neighbour $a $b $z]
    set _nc($id) $c
    if {$c < 0} { return }
    hole::_gen_triangle $a $b $c
    hole::_destroy $a $c $id b1a
    hole::_destroy $b $c $id b2a
    if {$_nb2a($id)} { hole::_add_edge $b $c $id b2a }
    if {$_nb1a($id)} {
        hole::_add_edge $a $c $id b1a
        set child1 [hole::_new_node $a $c $b]
        hole::_calc_tri $child1
    }
    # Re-read _nb2a($id): a nested destroy() during the base1 recursion above
    # may have just zeroed it (sos_triangle.c's own comment: "check is needed
    # here due to recursive nature of function").
    if {$_nb2a($id)} {
        set child2 [hole::_new_node $c $b $a]
        hole::_calc_tri $child2
    }
}

proc hole::_polygonize {start} {
    variable _dx; variable _dy; variable _dz; variable _ndots
    set sx $_dx($start); set sy $_dy($start); set sz $_dz($start)
    set best -1; set bestd2 0.0; set first 1
    for {set i 0} {$i < $_ndots} {incr i} {
        if {$i == $start} { continue }
        set ddx [expr {$_dx($i)-$sx}]; set ddy [expr {$_dy($i)-$sy}]; set ddz [expr {$_dz($i)-$sz}]
        set d2 [expr {$ddx*$ddx+$ddy*$ddy+$ddz*$ddz}]
        if {$first || $d2 < $bestd2} { set first 0; set bestd2 $d2; set best $i }
    }
    hole::_add_edge $start $best "" ""
    set root [hole::_new_node $start $best $best]
    hole::_calc_tri $root
}

proc hole::_color_name {c} {
    # NUMERIC comparison, not switch -exact: $c came through a float32
    # round-trip and can render as "7.0", which a string switch on "7"
    # silently misses, falling through to the "yellow" default.
    if {$c == 2} { return blue }
    if {$c == 3} { return red }
    if {$c == 7} { return green }
    return yellow
}

proc hole::_vmd_out {tris out_plot} {
    variable _dx; variable _dy; variable _dz
    variable _dnx; variable _dny; variable _dnz
    set fh [open $out_plot w]
    puts $fh "draw delete all"
    set cur ""
    foreach tr $tris {
        lassign $tr a b c col
        if {$col ne $cur} {
            set cur $col
            puts $fh "draw color [hole::_color_name $col]"
        }
        set line "draw trinorm "
        foreach v [list $a $b $c] {
            append line " { " [format "%8.3f" $_dx($v)] " " [format "%8.3f" $_dy($v)] \
                " " [format "%8.3f" $_dz($v)] " } "
        }
        foreach v [list $a $b $c] {
            append line " { " [format "%8.5f" $_dnx($v)] " " [format "%8.5f" $_dny($v)] \
                " " [format "%8.5f" $_dnz($v)] " } "
        }
        append line " "
        puts $fh $line
    }
    close $fh
}

proc hole::sos_triangle {sos_file out_plot} {
    variable _dx; variable _dy; variable _dz
    variable _dcol; variable _dnx; variable _dny; variable _dnz; variable _ndots
    variable _na; variable _nb; variable _nc; variable _nz; variable _nb1a; variable _nb2a
    variable _next_id; variable _edges; variable _tri_list; variable _used

    # Reset ALL namespace state - this proc is not re-entrant/threaded, but
    # IS called more than once per process (e.g. per trajectory frame).
    array unset _dx; array unset _dy; array unset _dz
    array unset _dcol; array unset _dnx; array unset _dny; array unset _dnz
    array unset _na; array unset _nb; array unset _nc; array unset _nz
    array unset _nb1a; array unset _nb2a
    array unset _edges; array unset _used
    set _next_id 0
    set _tri_list {}

    set raw [hole::_sos_read $sos_file]
    if {[llength $raw] > 30000} {
        error "hole::sos_triangle: MAX_COORD (30000) exceeded reading $sos_file"
    }
    set dots [hole::_cull_coords $raw]
    set _ndots [llength $dots]

    set i 0
    foreach d $dots {
        lassign $d x y z c nx ny nz
        set _dx($i) $x; set _dy($i) $y; set _dz($i) $z
        set _dcol($i) $c; set _dnx($i) $nx; set _dny($i) $ny; set _dnz($i) $nz
        incr i
    }

    # _calc_tri recurses roughly once per triangle along the advancing
    # front; measured depth 1833 on 1GRM's 2192-dot mesh. Tcl's default limit
    # (1000) is too low for any real surface, but this runs inside a host
    # process (the VMD plugin) that may have its own expectations - raise it
    # only for the duration of this call and restore whatever it was.
    set _saved_reclimit [interp recursionlimit {}]
    interp recursionlimit {} 200000
    for {set p 0} {$p < $_ndots} {incr p} {
        if {![info exists _used($p)]} {
            hole::_polygonize $p
        }
    }
    interp recursionlimit {} $_saved_reclimit

    # cull_triangles(): drop end-marker (color<0) triangles and
    # overlong-edge triangles.
    variable _max_vertex_len2
    set culled {}
    foreach tr $_tri_list {
        lassign $tr a b c col
        if {$_dcol($a) < 0 || $_dcol($b) < 0 || $_dcol($c) < 0} { continue }
        set lenAB [expr {($_dx($b)-$_dx($a))**2+($_dy($b)-$_dy($a))**2+($_dz($b)-$_dz($a))**2}]
        set lenAC [expr {($_dx($c)-$_dx($a))**2+($_dy($c)-$_dy($a))**2+($_dz($c)-$_dz($a))**2}]
        set lenCB [expr {($_dx($b)-$_dx($c))**2+($_dy($b)-$_dy($c))**2+($_dz($b)-$_dz($c))**2}]
        if {$lenAB > $_max_vertex_len2 || $lenAC > $_max_vertex_len2 || $lenCB > $_max_vertex_len2} { continue }
        lappend culled [list $a $b $c $col]
    }

    # reorder_triangle() (smooth mode only): flip winding so the triangle's
    # own geometric normal agrees with the sum of its vertex normals.
    set reordered {}
    foreach tr $culled {
        lassign $tr a b c col
        set tnx [expr {$_dnx($a)+$_dnx($b)+$_dnx($c)}]
        set tny [expr {$_dny($a)+$_dny($b)+$_dny($c)}]
        set tnz [expr {$_dnz($a)+$_dnz($b)+$_dnz($c)}]
        set abx [expr {$_dx($b)-$_dx($a)}]; set aby [expr {$_dy($b)-$_dy($a)}]; set abz [expr {$_dz($b)-$_dz($a)}]
        set acx [expr {$_dx($c)-$_dx($a)}]; set acy [expr {$_dy($c)-$_dy($a)}]; set acz [expr {$_dz($c)-$_dz($a)}]
        set crx [expr {$aby*$acz-$abz*$acy}]
        set cry [expr {$abz*$acx-$abx*$acz}]
        set crz [expr {$abx*$acy-$aby*$acx}]
        set d [expr {$crx*$tnx+$cry*$tny+$crz*$tnz}]
        if {$d < 0.0} {
            lappend reordered [list $b $a $c $col]
        } else {
            lappend reordered [list $a $b $c $col]
        }
    }

    hole::_vmd_out $reordered $out_plot
    return [list [llength $raw] $_ndots [llength $_tri_list] [llength $culled]]
}
