# The Tcl tunnel stage must reproduce the C's, intermediate by intermediate.
#
# Comparing only the surviving tunnels would be the weakest check yet: a wrong
# spline, a wrong 5-nearest set or a wrong control-path trim all still yield
# plausible tunnels of roughly the right length. So the comparison is against
# every intermediate the C dumps - openings, computed origins, Dijkstra's reach
# and total cost, each raw path, each control path, each profile's 100 sample
# radii and length, the bottleneck verdict and the similarity verdict.
#
# The Tcl driver is ::VMDHole::Mole::find_tunnels, which mirrors mole_main.c -
# the shipped C engine - so what is compared is the real call sequence.
#
# The "V" probe rows test the 5-nearest radius query on its own, at hundreds of
# points rather than the ~100 per tunnel the profiles give. Measured coverage:
# dropping to the ONE nearest atom moves 151 of 600 probes; raising it to the
# TEN nearest - what both method papers say, against the code - moves none, and
# moves no profile radius either. That is a property of the data, not of the
# probes: C/N/O/S vdW radii span ~0.3 A while distances spread far more, so the
# minimum of (distance - vdW) is essentially never attained past the fifth atom.
#
# Usage: test_mole_tunnel.tcl ATOMS.txt TUNNEL.txt [max_pivots] [max_similarity]
#                             [bottleneck] [tolerance]

set here [file dirname [info script]]
# The engine is inlined in vmdhole.tcl now (one shipped script + the
# compiled binaries), so this extracts it instead of sourcing modules.
source [file join [file dirname [file normalize [info script]]] mole_tcl_engine.tcl]

set probes {}
set fails 0
proc report {label ok {detail ""}} {
    global fails
    if {!$ok} { incr fails }
    puts [format "  %-52s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
}

# ---------------------------------------------------------------- reference
set f [open [lindex $argv 1] r]
set ref {}
while {[gets $f line] >= 0} {
    if {[string trim $line] eq ""} continue
    if {[lindex $line 0] eq "H"} {
        lassign $line _ ref_np ref_nt ref_nc ref_nch ref_nvd
        continue
    }
    if {[lindex $line 0] eq "V"} { lappend probes [lrange $line 1 end]; continue }
    lappend ref $line
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
::VMDHole::Mole::dt_build M $piv $np
::VMDHole::Mole::complex_build C M $piv $rad 3.0 1.25 8
lassign [::VMDHole::Mole::cavities C cav 8 5.0] nc nch nvd
report "cavities = $ref_nc, channels = $ref_nch" \
       [expr {$nc == $ref_nc && $nch == $ref_nch}] "(got $nc/$nch)"

# The 5-nearest radius query on its own, at every probe point.
set nrbad 0
foreach pr $probes {
    lassign $pr qx qy qz want
    set got [::VMDHole::Mole::radius_at $qx $qy $qz $piv $rad $np 0]
    if {[format %.17g $got] ne $want} { incr nrbad }
}
report "5-nearest radius identical at [llength $probes] probes" \
       [expr {$nrbad == 0}] "($nrbad differ)"

set maxsim [expr {[llength $argv] > 3 ? [lindex $argv 3] : 0.9}]
set bneck  [expr {[llength $argv] > 4 ? [lindex $argv 4] : 1.25}]
set btol   [expr {[llength $argv] > 5 ? [lindex $argv 5] : 0.0}]
set t0 [clock milliseconds]
set tunnels [::VMDHole::Mole::find_tunnels C $cav $piv $rad $np \
                 [list max_similarity $maxsim bottleneck $bneck \
                       bottleneck_tolerance $btol] trace]
set ms [expr {[clock milliseconds] - $t0}]

# ------------------------------------------ compare the trace, line by line
# Every line is one intermediate. Counting the kinds separately says WHICH stage
# broke - openings, origins, Dijkstra, path, control path, profile or filter -
# rather than only that something did.
array set kind {G origins O openings D dijkstra P path K control-path
                R profile S similarity T tunnels}
array set nbad {}
foreach k [array names kind] { set nbad($kind($k)) 0 }
set firstbad ""
set n [expr {max([llength $ref], [llength $trace])}]
for {set i 0} {$i < $n} {incr i} {
    set a [lindex $trace $i]
    set b [lindex $ref $i]
    if {$a eq $b} continue
    set tag [lindex $b 0]
    if {$tag eq "" || ![info exists kind($tag)]} { set tag [lindex $a 0] }
    if {[info exists kind($tag)]} { incr nbad($kind($tag)) } else { incr nbad(tunnels) }
    if {$firstbad eq ""} {
        set firstbad "line $i:\n         tcl [string range $a 0 150]\n         c   [string range $b 0 150]"
    }
}
report "trace has [llength $ref] lines" \
       [expr {[llength $trace] == [llength $ref]}] \
       "(got [llength $trace])"
set tot 0
set detail ""
foreach k {G O D P K R S T} {
    incr tot $nbad($kind($k))
    if {$nbad($kind($k))} { append detail " $kind($k)=$nbad($kind($k))" }
}
report "every intermediate identical" [expr {$tot == 0}] "($detail\n       $firstbad)"

set reft [lindex $ref end]
report "tunnels: [lrange $reft 1 end]" [expr {[lindex $trace end] eq $reft}] \
       "(got [lrange [lindex $trace end] 1 end])"

puts [format "  %-52s %d ms" "tunnel stage" $ms]
exit [expr {$fails ? 1 : 0}]
