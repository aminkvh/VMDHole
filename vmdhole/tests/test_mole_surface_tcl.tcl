# The pure-Tcl engine must reproduce the SurfaceCavity features, not just the
# cavity ones.
#
# Custom exits and paths were added to the C first, and the byte-identity test
# only exercises the DEFAULT path - a Tcl that ignored --exit and --path
# entirely would still pass it. The plugin falls back to Tcl whenever the binary
# is missing, so a user on that path would silently get no exits and no paths.
#
# Expected values are MOLE's own, read from the committed references rather than
# hard-coded here.
#
# Usage: test_mole_surface_tcl.tcl ATOMS.txt REFDIR_SURFACE_EXIT REFDIR_PATHS

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

# MOLE's own answer, from tunnels.csv's second column.
proc mole_length {dir} {
    set fh [open [file join $dir tunnels.csv] r]
    set rows [split [string trim [read $fh]] "\n"]
    close $fh
    return [lindex [split [lindex $rows 1] ","] 1]
}

lassign $argv atoms refexit refpath refmix

::VMDHole::Mole::read_atoms $atoms A
set piv [::VMDHole::Mole::pivot_points A rad src pbfac pfree pbb pres restab]
set np [expr {[llength $piv] / 3}]
# MOLE's own triangulation, which is what the C engine uses.
::VMDHole::Mole::dh_build M $piv $np
::VMDHole::Mole::complex_build_dh C M $piv $rad 3.0 1.25 8
lassign [::VMDHole::Mole::cavities C cav 8 5.0] nc nch nvd

# The A308/A309 centroid, the origin MOLE's own configs use.
set sx 0.0; set sy 0.0; set sz 0.0; set n 0
for {set i 0} {$i < $A(n)} {incr i} {
    if {[lindex $A(chain) $i] ne "A"} continue
    set r [lindex $A(seq) $i]
    if {$r != 308 && $r != 309} continue
    set b [expr {3 * $i}]
    set sx [expr {$sx + [lindex $A(xyz) $b]}]
    set sy [expr {$sy + [lindex $A(xyz) [expr {$b+1}]]}]
    set sz [expr {$sz + [lindex $A(xyz) [expr {$b+2}]]}]
    incr n
}
set origin [list [expr {$sx/$n}] [expr {$sy/$n}] [expr {$sz/$n}]]

# 1. custom exit into the SurfaceCavity
set want [mole_length $refexit]
set t [::VMDHole::Mole::find_tunnels C $cav $piv $rad $np \
    [list exit_point {-22.841 -55.065 -20.445} exits_only 1 origin $origin] \
    "" $pbfac $pfree $pres $restab]
set got [expr {[llength $t] == 1 ? [format "%.2f" [lindex [lindex $t 0] 0]] : "none"}]
report "Tcl reproduces MOLE's surface-exit tunnel ($want A)" \
       [expr {$got eq $want}] "(got $got)"

# 2. path between two given points
set want [mole_length $refpath]
set t [::VMDHole::Mole::find_tunnels C $cav $piv $rad $np \
    [list path_points [concat $origin {-27.435 0.017 -11.724}]] \
    "" $pbfac $pfree $pres $restab]
set got [expr {[llength $t] == 1 ? [format "%.2f" [lindex [lindex $t 0] 0]] : "none"}]
report "Tcl reproduces MOLE's path ($want A)" [expr {$got eq $want}] "(got $got)"

# 3. a custom exit ALONGSIDE the regular openings (exits_only 0). The Tcl used
#    computed openings only and gated its surface search on exits_only, so the
#    exit contributed nothing - a silent wrong answer on the fallback path,
#    invisible to check 2 above because that one sets exits_only.
set t [::VMDHole::Mole::find_tunnels C $cav $piv $rad $np \
    [list exit_point {-22.841 -55.065 -20.445} origin $origin] \
    "" $pbfac $pfree $pres $restab]
# Exact and IN ORDER: MOLE's TunnelComparer is (Cavity.Id, Length) with the
# SurfaceCavity at Id 0, so a sorted or membership comparison would pass while
# the surface tunnel sat in the wrong place.
set want {}
set fh [open [file join $refmix tunnels.csv] r]
foreach row [lrange [split [string trim [read $fh]] "\n"] 1 end] {
    lappend want [lindex [split $row ","] 1]
}
close $fh
set got {}
foreach tn $t { lappend got [format "%.2f" [lindex $tn 0]] }
report "Tcl reproduces MOLE's mixed custom-exit run, in order" \
       [expr {$got eq $want}] "(got $got want $want)"

exit [expr {$fails ? 1 : 0}]
