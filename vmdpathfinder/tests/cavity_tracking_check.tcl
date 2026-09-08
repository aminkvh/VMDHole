package provide Tk 8.5
set here [file dirname [info script]]
source [file join $here .. vmdpathfinder.tcl]
# Cavity identity across frames.
#
# MOLE recomputes cavities independently in every frame and ranks them by
# volume, so a RANK IS NOT AN IDENTITY: two pockets that swap volume order swap
# ranks. Anything user-facing keyed on the rank - the draw tick, the colour, a
# mean volume - then silently follows whichever pocket holds that rank next
# frame. This is the same defect the routes had, and it is invisible to a test
# whose frames are duplicates of one another, so this one builds the awkward
# case directly: two pockets 30 A apart whose volumes swap between two frames.
#
# Pure Tcl, no VMD: _cavity_tracks reads tunnel_lining and nothing else.
set fails 0
proc report {label ok {detail ""}} {
    global fails
    if {!$ok} { incr fails }
    puts [format "  %-56s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
}
proc chk {n g w} { report $n [expr {$g eq $w}] "(got '$g' want '$w')" }
# Two pockets, A near the origin and B 30 A away. Their VOLUMES swap between
# frames, so MOLE's rank-by-volume swaps too - the exact case where a rank is
# not an identity.
proc mk {vol cx} {
    set sph {}
    foreach d {0 1 2} { lappend sph [list [expr {$cx+$d}] 0 0 2.0] }
    return [dict create type Cavity volume $vol depth 5 depthlen 5.0 \
        nboundary 3 ninner 1 bprops {} iprops {} bres {} ires {} spheres $sph \
        origins [list [list 1 [list $cx 0 0] 9.9]]]
}
# frame 0: A is bigger (rank 1); frame 1: B is bigger (rank 1)
set ::VMDPathFinder::tunnel_lining(0) [dict create cav.1 [mk 900.0 0] cav.2 [mk 100.0 30]]
set ::VMDPathFinder::tunnel_lining(1) [dict create cav.1 [mk 900.0 30] cav.2 [mk 100.0 0]]
set ::VMDPathFinder::tunnel_result_frames {0 1}
set ::VMDPathFinder::plot_data_version 1
set tracks [::VMDPathFinder::_cavity_tracks]
puts "tracks: [llength $tracks]"
foreach t $tracks {
    puts [format "  tid %s  centroid %s  ranks %s  seen %.0f%%  mean %.0f" \
        [dict get $t tid] [dict get $t centroid] [dict get $t ids] [dict get $t seen] [dict get $t vol_mean]]
}
chk "the two pockets stay two tracks (not four)" [llength $tracks] 2
# pocket A sits at x~1, pocket B at x~31 - each track must hold one of each
set tA {}; set tB {}
foreach t $tracks {
    if {[lindex [dict get $t centroid] 0] < 15} { set tA $t } else { set tB $t }
}
chk "a track exists for each physical pocket" [expr {[llength $tA] && [llength $tB]}] 1
chk "pocket A is rank 1 in frame 0 but rank 2 in frame 1" \
    [list [dict get $tA ids 0] [dict get $tA ids 1]] {1 2}
chk "pocket B is rank 2 in frame 0 but rank 1 in frame 1" \
    [list [dict get $tB ids 0] [dict get $tB ids 1]] {2 1}
chk "both seen in 100% of frames" \
    [list [dict get $tA seen] [dict get $tB seen]] {100.0 100.0}
chk "each track's mean is the mean of ITS pocket, not of a rank" \
    [list [dict get $tA vol_mean] [dict get $tB vol_mean]] {500.0 500.0}
# a pocket missing from one frame shows as < 100% seen
set ::VMDPathFinder::tunnel_lining(1) [dict create cav.1 [mk 900.0 30]]
set ::VMDPathFinder::plot_data_version 2
set tracks2 [::VMDPathFinder::_cavity_tracks]
set seens {}
foreach t $tracks2 { lappend seens [format %.0f [dict get $t seen]] }
puts "seen with pocket A absent from frame 1: [lsort $seens]"
chk "a pocket absent from a frame is seen in only half of them" [lsort -integer $seens] {50 100}
puts "  ---- cavity-tracking checks, $fails failed"
exit [expr {$fails ? 1 : 0}]
