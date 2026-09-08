# The plugin's own MOLE call must give the SAME answer with the binary and
# without it. Run under `vmd -dispdev text` by test_mole_tcl_port.sh.
#
# This tests the GLUE, which nothing else does: the dispatch in
# _tunnel_search_mole, the on-demand loading of the Tcl modules from
# plugin_dir, and _tunnel_atoms_mole feeding both paths. The engine itself is
# checked elsewhere against dumps and against the C's output file; here the
# input comes from VMD and the caller is the plugin.
#
# The binary is taken out of reach by overriding tool_path rather than
# by moving files, so the test cannot damage a working install.
package provide Tk 8.5
source [file join $env(VH) vmdpathfinder.tcl]
mol new $env(PDB) type pdb waitfor all
set ::VMDPathFinder::state(selection) $env(SEL)
if {[info exists env(SOS)]} { set ::VMDPathFinder::state(sos_triangle_exec) $env(SOS) }

# NON-DEFAULT parameters on purpose. Running both paths at MOLE's defaults
# would pass even if _tunnel_cfg were decorative and neither engine ever read
# the panel - which is the same defect class as the dispatcher that was never
# called. A bogus weight name is included: it must fall back to VoronoiScale
# rather than to whatever index a lookup returned.
set ::VMDPathFinder::state(mole_probe) 4.0
set ::VMDPathFinder::state(mole_interior) 1.4
set ::VMDPathFinder::state(mole_weight_disp) LengthAndRadius
set cfg [::VMDPathFinder::_tunnel_cfg]

set t0 [clock milliseconds]
set c [::VMDPathFinder::_tunnel_search_mole 0 0 {} $cfg]
set cms [expr {[clock milliseconds] - $t0}]

rename ::VMDPathFinder::tool_path ::VMDPathFinder::_real_tool_path
proc ::VMDPathFinder::tool_path {name} { if {$name eq "mole_engine"} { return "" }; ::VMDPathFinder::_real_tool_path $name }
set t0 [clock milliseconds]
set t [::VMDPathFinder::_tunnel_search_mole 0 0 {} $cfg]
set tms [expr {[clock milliseconds] - $t0}]

# The engine the DISPATCHER picks. This is its own check, not a detail of the
# one above: _tunnel_search_mole existed and was called by nothing for a whole
# release, so the selector said "mole" while every run did the lattice search.
# Anything but a mole-* tag here means that has come back.
rename ::VMDPathFinder::tool_path {}
rename ::VMDPathFinder::_real_tool_path ::VMDPathFinder::tool_path
set dispatch [::VMDPathFinder::tunnel_search 0 0 {} [::VMDPathFinder::_tunnel_cfg]]
# render_tunnels_for_frame reads tunnel_results; run_tunnel_analysis normally
# fills it, and this harness calls tunnel_search directly.
set ::VMDPathFinder::tunnel_results(0) [lindex $dispatch 1]

set out [open $env(OUT) w]
# The lining must be stored by the same call that stores the tunnels - the
# Lining window reads tunnel_lining($frame) and nothing else populates it.
set nlin 0; set nlayers 0
if {[info exists ::VMDPathFinder::tunnel_lining(0)]} {
    set L $::VMDPathFinder::tunnel_lining(0)
    foreach k [dict keys $L] {
        if {[string match "*.layers" $k]} {
            incr nlin
            incr nlayers [llength [dict get $L $k]]
        }
    }
}
# Property coloring: the spheres handed to colorize_by_sphere_values must
# cover the tunnel and must actually VARY, or the tunnel renders one flat
# color and the feature is decorative. Checked here rather than in the
# renderer because the renderer needs a display.
set psph 0; set pdistinct 0; set prange 0
if {$nlin > 0} {
    set sp [::VMDPathFinder::_tunnel_property_spheres 0 1 hydropathy]
    set psph [llength $sp]
    set vals {}
    foreach e $sp { lappend vals [lindex $e 4] }
    set pdistinct [llength [lsort -unique $vals]]
    lassign [::VMDPathFinder::_tunnel_property_range 0 hydropathy] plo phi
    set prange [expr {$phi > $plo}]
}
puts $out "prop_spheres $psph"
puts $out "prop_distinct $pdistinct"
puts $out "prop_range $prange"
puts $out "lining_tunnels $nlin"
puts $out "lining_layers $nlayers"
puts $out "dispatch_engine [lindex $dispatch 0]"
puts $out "dispatch_tunnels [llength [lindex $dispatch 1]]"
set ntun [llength [lsearch -all [split [lindex $c 1] "\n"] "T *"]]
puts $out "c_status [lindex $c 0]"
puts $out "tcl_status [lindex $t 0]"
puts $out "tunnels $ntun"
# The C engine writes the cavity VP records (and their '# VP' header line)
# ahead of the Tcl engine. While the Tcl side writes none they are set aside
# and REPORTED as vp_pending, so the rest of the text stays under identity;
# once the Tcl engine writes any VP line, every VP line is compared too.
set ct [lindex $c 1]; set tt [lindex $t 1]
set vp_re {^(VP |# VP )}
set nvp_c [llength [lsearch -all -regexp [split $ct "\n"] $vp_re]]
set nvp_t [llength [lsearch -all -regexp [split $tt "\n"] $vp_re]]
if {$nvp_c > 0 && $nvp_t == 0} {
    set ct [join [lsearch -all -inline -not -regexp [split $ct "\n"] $vp_re] "\n"]
}
puts $out "vp_pending [expr {$nvp_c > 0 && $nvp_t == 0 ? $nvp_c : 0}]"
puts $out "identical [expr {$ct eq $tt}]"
puts $out "timing $cms $tms"
puts $out "tcl_diag [lindex $t 2]"
puts $out "cfg_probe [dict get $cfg mole_probe]"
puts $out "cfg_weight [dict get $cfg mole_weight]"
set ::VMDPathFinder::state(mole_weight_disp) "not-a-weight-function"
puts $out "cfg_weight_bogus [dict get [::VMDPathFinder::_tunnel_cfg] mole_weight]"
if {$ct ne $tt} {
    foreach x [split $ct "\n"] y [split $tt "\n"] {
        if {$x ne $y} { puts $out "firstdiff_c $x"; puts $out "firstdiff_tcl $y"; break }
    }
}
close $out
quit
