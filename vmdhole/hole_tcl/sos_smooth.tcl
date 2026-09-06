# hole::sos_smooth - the pure-Tcl port of sos_triangle --sos-smooth.
#  A local average of dot clouds: every dot of the centre frame moves to the
#  mean of itself and its nearest same-facing dot (within rho, normals
#  agreeing) in each window frame; its normal is the renormalised mean of
#  theirs. Header and centreline records are copied as they are. The sums run
#  in the same order as the C, so the two are byte-identical.

namespace eval hole {}

proc hole::_sm_read {path} {
    # dots only: {x y z nx ny nz} per record of type 4
    set fh [open $path r]
    set dots {}
    while {[gets $fh line] >= 0} {
        set v [regexp -all -inline {[-+0-9.eE]+} $line]
        if {[llength $v] < 7 || [lindex $v 0] != 4.0} continue
        lappend dots [lrange $v 1 6]
    }
    close $fh
    return $dots
}

proc hole::_sm_grid {dots cell arrname} {
    # cell -> ascending list of dot indices; returns {ox oy oz}
    upvar 1 $arrname g
    array unset g
    set lo {1e30 1e30 1e30}; set hi {-1e30 -1e30 -1e30}
    foreach d $dots {
        foreach j {0 1 2} {
            set q [lindex $d $j]
            if {$q < [lindex $lo $j]} { lset lo $j $q }
            if {$q > [lindex $hi $j]} { lset hi $j $q }
        }
    }
    if {![llength $dots]} { set lo {0 0 0} }
    set ox [expr {[lindex $lo 0]-$cell}]; set oy [expr {[lindex $lo 1]-$cell}]; set oz [expr {[lindex $lo 2]-$cell}]
    set i 0
    foreach d $dots {
        lassign $d x y z
        set key "[expr {int(($x-$ox)/$cell)}],[expr {int(($y-$oy)/$cell)}],[expr {int(($z-$oz)/$cell)}]"
        lappend g($key) $i
        incr i
    }
    return [list $ox $oy $oz]
}

proc hole::_sm_nearest {dots arrname origin cell x y z nx ny nz rho} {
    upvar 1 $arrname g
    lassign $origin ox oy oz
    set ix [expr {int(($x-$ox)/$cell)}]; set iy [expr {int(($y-$oy)/$cell)}]; set iz [expr {int(($z-$oz)/$cell)}]
    set best [expr {$rho*$rho}]; set bi -1
    foreach dx {-1 0 1} { foreach dy {-1 0 1} { foreach dz {-1 0 1} {
        set key "[expr {$ix+$dx}],[expr {$iy+$dy}],[expr {$iz+$dz}]"
        if {![info exists g($key)]} continue
        foreach i $g($key) {
            lassign [lindex $dots $i] qx qy qz qnx qny qnz
            set ex [expr {$qx-$x}]; set ey [expr {$qy-$y}]; set ez [expr {$qz-$z}]
            set d2 [expr {$ex*$ex+$ey*$ey+$ez*$ez}]
            if {$d2 > $best} continue
            if {$qnx*$nx+$qny*$ny+$qnz*$nz <= 0.0} continue
            if {$d2 < $best || $bi < 0 || $i < $bi} { set best $d2; set bi $i }
        }
    }}}
    return $bi
}

proc hole::sos_smooth {out rho centre with_list} {
    set clouds {}
    set k 0
    foreach w $with_list {
        set dots [hole::_sm_read $w]
        set origin [hole::_sm_grid $dots $rho grid$k]
        lappend clouds [list $dots $origin]
        incr k
    }
    set fin [open $centre r]
    set fout [open $out w]
    set ndots 0
    while {[gets $fin line] >= 0} {
        set v [regexp -all -inline {[-+0-9.eE]+} $line]
        if {[llength $v] < 7 || [lindex $v 0] != 4.0} { puts $fout $line; continue }
        lassign $v _t x y z nx ny nz
        set sx [expr {double($x)}]; set sy [expr {double($y)}]; set sz [expr {double($z)}]
        set snx [expr {double($nx)}]; set sny [expr {double($ny)}]; set snz [expr {double($nz)}]
        set cnt 1
        set j 0
        foreach c $clouds {
            lassign $c dots origin
            set i [hole::_sm_nearest $dots grid$j $origin $rho $x $y $z $nx $ny $nz $rho]
            incr j
            if {$i < 0} continue
            lassign [lindex $dots $i] qx qy qz qnx qny qnz
            set sx [expr {$sx+$qx}]; set sy [expr {$sy+$qy}]; set sz [expr {$sz+$qz}]
            set snx [expr {$snx+$qnx}]; set sny [expr {$sny+$qny}]; set snz [expr {$snz+$qnz}]
            incr cnt
        }
        set sx [expr {$sx/$cnt}]; set sy [expr {$sy/$cnt}]; set sz [expr {$sz/$cnt}]
        set nl [expr {sqrt($snx*$snx+$sny*$sny+$snz*$snz)}]
        if {$nl > 1e-12} { set snx [expr {$snx/$nl}]; set sny [expr {$sny/$nl}]; set snz [expr {$snz/$nl}] } \
        else { set snx $nx; set sny $ny; set snz $nz }
        puts $fout [format "%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" 4.0 $sx $sy $sz $snx $sny $snz]
        incr ndots
    }
    close $fin; close $fout
    return $ndots
}
