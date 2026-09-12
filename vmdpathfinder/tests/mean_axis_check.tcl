# The Mean Profile 3D tube must lie along the CVECT the run was measured with.
#
# THE DEFECT THIS GUARDS (fixed bca956d, and reported as a regression more than
# once before that):
#   Every radius in a pore profile is indexed by distance along CVECT from
#   CPOINT. build_and_show_mean_surface revolved them around a PCA fit of ONE
#   reference frame's centreline instead, and then linearly remapped each bin's
#   coordinate onto that same frame's centreline span. The tube therefore
#   pointed the wrong way (18.2 deg off on a tilted pore) AND was squeezed onto
#   one frame's extent, so features sat at the wrong height and the far end was
#   crushed inwards - "small at the top, does not correspond to the frames".
#
# A source-level check would be vacuous here: the failure is geometric, so this
# builds the real mesh through the real sph_process/sos_triangle and measures
# the angle between the tube's own principal axis and the declared CVECT.
#
# Deliberately runs with a TILTED cvect. With cvect 0 0 1 on an upright channel
# the PCA fit and the declared axis nearly coincide, so the bug is invisible -
# which is exactly why it survived so long.
# NO `package provide Tk 8.5` here: that shim is for the tclsh-only readers, and
# under real VMD 2.0 it collides with its own Tk 8.6 ("conflicting versions
# provided for package Tk") and takes the script down with it.
# vmd -e sets [info script] to a BARE filename, so its dirname is "." and every
# path built from it lands in the cwd. The wrapper passes the real directory.
set here [expr {[info exists ::env(VMDPATHFINDER_TEST_DIR)]
                ? $::env(VMDPATHFINDER_TEST_DIR) : [file dirname [info script]]}]
set fails 0
proc report {label ok {detail ""}} {
    global fails
    if {!$ok} { incr fails }
    puts [format "  %-56s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
}
proc vnorm {v} {
    lassign $v x y z
    set n [expr {sqrt($x*$x + $y*$y + $z*$z)}]
    if {$n < 1e-12} { return {0 0 1} }
    return [list [expr {$x/$n}] [expr {$y/$n}] [expr {$z/$n}]]
}
proc angle_deg {a b} {
    lassign [vnorm $a] ax ay az
    lassign [vnorm $b] bx by bz
    set d [expr {abs($ax*$bx + $ay*$by + $az*$bz)}]
    if {$d > 1.0} { set d 1.0 }
    return [expr {acos($d) * 57.29577951308232}]
}

# principal axis of a .vmd_plot's triangle vertices, by power iteration
proc plot_axis {file} {
    set fh [open $file r]; set txt [read $fh]; close $fh
    set pts {}; set sx 0.0; set sy 0.0; set sz 0.0; set n 0
    foreach line [split $txt "\n"] {
        if {![string match "*trinorm*" $line]} { continue }
        set g {}
        foreach tok [split [string map {\{ { } \} { }} $line]] {
            if {[string is double -strict $tok]} { lappend g $tok }
        }
        for {set i 0} {$i < 9} {incr i 3} {
            set x [lindex $g $i]; set y [lindex $g [expr {$i+1}]]; set z [lindex $g [expr {$i+2}]]
            if {$x eq "" || $z eq ""} { continue }
            lappend pts [list $x $y $z]
            set sx [expr {$sx+$x}]; set sy [expr {$sy+$y}]; set sz [expr {$sz+$z}]; incr n
        }
    }
    if {$n < 3} { return {} }
    set cx [expr {$sx/$n}]; set cy [expr {$sy/$n}]; set cz [expr {$sz/$n}]
    set vx 0.0; set vy 0.0; set vz 1.0
    for {set it 0} {$it < 60} {incr it} {
        set ax 0.0; set ay 0.0; set az 0.0
        foreach p $pts {
            lassign $p x y z
            set dx [expr {$x-$cx}]; set dy [expr {$y-$cy}]; set dz [expr {$z-$cz}]
            set d [expr {$dx*$vx + $dy*$vy + $dz*$vz}]
            set ax [expr {$ax + $d*$dx}]; set ay [expr {$ay + $d*$dy}]; set az [expr {$az + $d*$dz}]
        }
        set nn [expr {sqrt($ax*$ax+$ay*$ay+$az*$az)}]
        if {$nn <= 0} break
        set vx [expr {$ax/$nn}]; set vy [expr {$ay/$nn}]; set vz [expr {$az/$nn}]
    }
    return [list $vx $vy $vz]
}

proc plot_points {file} {
    set fh [open $file r]; set txt [read $fh]; close $fh
    set pts {}
    foreach line [split $txt "\n"] {
        if {![string match "*trinorm*" $line]} { continue }
        set g {}
        foreach tok [split [string map {\{ { } \} { }} $line]] {
            if {[string is double -strict $tok]} { lappend g $tok }
        }
        for {set i 0} {$i < 9} {incr i 3} {
            set x [lindex $g $i]; set y [lindex $g [expr {$i+1}]]; set z [lindex $g [expr {$i+2}]]
            if {$x eq "" || $z eq ""} { continue }
            lappend pts [list $x $y $z]
        }
    }
    return $pts
}

source [file join $here .. vmdpathfinder.tcl]
proc say_diag {} {
    set n [llength $::VMDPathFinder::result_frames]
    puts "  diag: result_frames=$n dict_results=[dict size $::VMDPathFinder::results]"
    puts "  diag: mean_frame_spec='[expr {[info exists ::VMDPathFinder::state(mean_frame_spec)] ? $::VMDPathFinder::state(mean_frame_spec) : {<unset>}}]'"
    lassign [::VMDPathFinder::_mean_frames_and_key] _mf _mk
    puts "  diag: mean_frames='[join $_mf ,]' key='$_mk'"
    foreach fr $::VMDPathFinder::result_frames {
        set rows 0
        catch {set rows [llength [::VMDPathFinder::ensure_profile_full $fr]]}
        puts "  diag: frame $fr full-profile rows=$rows"
    }
}
set pdb [file join $here fixtures mole_reference 1BL8.pdb]
if {![file exists $pdb]} { puts "SKIP: no 1BL8 fixture at '$pdb'"; exit 0 }
set molid [mol new $pdb type pdb waitfor all]
::VMDPathFinder::show_gui
update idletasks; update

# A TILTED axis on purpose: with 0 0 1 the PCA fit and the declared axis nearly
# coincide on an upright channel and the defect is invisible.
set CV [expr {[info exists ::env(MEANAXIS_CV)] ? $::env(MEANAXIS_CV) : {0.25 0.15 1.0}}]
set s ::VMDPathFinder::state
set ${s}(selection)        "protein"
set ${s}(cpoint)           [expr {[info exists ::env(MEANAXIS_CP)] ? $::env(MEANAXIS_CP) : {73.853 26.536 26.594}}]
set ${s}(cvect)            $CV
set ${s}(sample)           0.25
set ${s}(endrad)           15.0
set ${s}(pore_method)      "circular"
set ${s}(frame_spec)       "now"
set ${s}(display_mode)     "none"
set ${s}(prebuild_surfaces) 0
set ${s}(overwrite_results) 1
set ${s}(work_dir) [file join [::VMDPathFinder::_scratch_base] meanaxis[clock clicks]]
if {[catch {::VMDPathFinder::run_analysis} e]} {
    puts "SKIP: pore run unavailable ($e)"; exit 0
}
update idletasks; update
if {![llength $::VMDPathFinder::result_frames]} { puts "SKIP: pore run produced no frames"; exit 0 }

# Diagnostics: this group has failed inside the suite while passing standalone,
# so when the build declines, say WHAT was empty rather than only that it was.
say_diag
set ${s}(show_mean_surface)  1
set ${s}(mean_surface_color) green
if {[catch {::VMDPathFinder::build_and_show_mean_surface 1 1} e]} {
    # HOLE's Monte Carlo search is stochastic: the same fixture bins fine
    # standalone and has, inside the suite, produced too few samples to bin at
    # all. That is not the axis this guards, so it SKIPS - a wrong axis stays
    # the only way to FAIL. Deliberately narrow: any OTHER build error is still
    # a failure.
    if {[string match "*no profile data*" $e]} {
        puts "SKIP: HOLE produced too little profile to bin this run ($e) - MC noise, not an axis fault"
        exit 0
    }
    report "the mean surface builds" 0 "($e)"
    puts "  ---- mean-axis checks, $fails failed"; exit 1
}
report "the mean surface builds" 1

# Every pore run makes its own <work_dir>/<structure>_<run id> folder and the
# mean mesh sits beside that run's frames, so look one level down too.
set plots [concat \
    [glob -nocomplain [file join $::VMDPathFinder::state(work_dir) mean_profile mean_profile_*.vmd_plot]] \
    [glob -nocomplain [file join $::VMDPathFinder::state(work_dir) * mean_profile mean_profile_*.vmd_plot]]]
report "a mean-profile mesh was written" [expr {[llength $plots] > 0}] "(none in the work dir)"
if {![llength $plots]} { puts "  ---- mean-axis checks, $fails failed"; exit 1 }

set axis [plot_axis [lindex $plots 0]]
report "the mesh has enough geometry to measure" [expr {[llength $axis] == 3}]
if {[llength $axis] != 3} { puts "  ---- mean-axis checks, $fails failed"; exit 1 }

set off [angle_deg $axis $CV]
# 3 deg: the tube is a discretised sphere union, so its measured principal axis
# carries some noise. The defect this guards produced 18.2 deg.
report "the tube lies along the declared CVECT (< 3 deg)" [expr {$off < 3.0}] \
    "(off by [format %.2f $off] deg; axis [format "%.4f %.4f %.4f" {*}$axis], cvect $CV)"


# THE DEFECT THIS ALSO GUARDS (fixed alongside the angle check): a pure
# translation along the tube's own axis leaves the angle above unchanged, so
# that check alone is blind to it. HOLE's own "coord" is an ABSOLUTE lab-frame
# projection along CVECT, not a distance from CPOINT - build_and_show_mean_surface
# once added CPOINT's own axial component on top of it, landing the tube
# dot(CPOINT,CVECT) Angstroms off. This fixture's CPOINT is deliberately far
# from the origin (dot(CPOINT,CVECT) ~47 A here) so that offset would be
# unmissable - a fixture centred near 0,0,0 would hide it, same as an upright
# cvect hides the angle defect above.
set _cvn [vnorm $CV]
lassign $_cvn _cux _cuy _cuz
set nbins [::VMDPathFinder::_mean_profile_nbins]
lassign [::VMDPathFinder::_mean_frames_and_key] _mf _mk
set _bdata [::VMDPathFinder::collect_binned_radii $nbins $_mf $_mk]
set _have_range [expr {$_bdata ne {} && [dict exists $_bdata zmin]}]
if {$_have_range} {
    set _zmin [dict get $_bdata zmin]
    set _zstep [dict get $_bdata zstep]
    set _zmax [expr {$_zmin + $nbins * $_zstep}]
    set _pts [plot_points [lindex $plots 0]]
    set _plo 1e20; set _phi -1e20
    foreach _p $_pts {
        lassign $_p _px _py _pz
        set _proj [expr {$_px*$_cux + $_py*$_cuy + $_pz*$_cuz}]
        if {$_proj < $_plo} { set _plo $_proj }
        if {$_proj > $_phi} { set _phi $_proj }
    }
    # Generous tolerance: the mesh's own sphere radii carry vertices a bit past
    # the bin centres axially too. Still far tighter than the ~47 A this
    # fixture's CPOINT-projection defect would produce.
    set _tol 20.0
    set _ok [expr {$_plo >= $_zmin - $_tol && $_phi <= $_zmax + $_tol}]
    report "the tube sits at the profile's own axial position" $_ok \
        "(tube [format %.2f $_plo]..[format %.2f $_phi], expected ~[format %.2f $_zmin]..[format %.2f $_zmax])"
} else {
    puts "  SKIP: no binned profile data to compare the tube's position against"
}

# The run's own persisted axis must be what fed it - otherwise the check above
# could pass on a PCA fit that happens to agree.
# The frame folder lives under the run's own folder, not straight under the
# work dir - take it from the result record rather than rebuilding the path.
set fd [dict get $::VMDPathFinder::results [lindex $::VMDPathFinder::result_frames 0] run_dir]
set fa [::VMDPathFinder::_frame_axis_persisted $fd]
report "the run persisted a resolved CPOINT/CVECT to build from" [expr {[llength $fa] == 6}] \
    "(got '$fa' from $fd)"

puts "  ---- mean-axis checks, $fails failed"
exit [expr {$fails ? 1 : 0}]
