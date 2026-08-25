# filter_bottleneck's n<1 branch: MOLE keeps a fully-degenerate tunnel
# UNCONDITIONALLY, derived from FilterBottleneck's own NaN comparisons - see
# the comment on the branch in mole_tunnel.tcl. Mirrors
# native/mole/mole_filter_bottleneck_test.c so both engines are checked the
# same way. No real structure reaches this branch (shortest candidate profile
# across the fixtures is 3.47 A; see check 8b in test_mole_tcl_port.sh).
set here [file dirname [info script]]
# The engine is inlined in vmdhole.tcl now (one shipped script + the
# compiled binaries), so this extracts it instead of sourcing modules.
source [file join [file dirname [file normalize [info script]]] mole_tcl_engine.tcl]

proc make_profile {length} {
    ::VMDHole::Mole::spline_init sx {0.0 1.0} {0.0 100.0} 2
    ::VMDHole::Mole::spline_init sy {0.0 1.0} {0.0 100.0} 2
    ::VMDHole::Mole::spline_init sz {0.0 1.0} {0.0 100.0} 2
    ::VMDHole::Mole::spline_init sr {0.0 1.0} {0.01 0.01} 2
    return [dict create length $length \
        sx [list $sx(n) $sx(t) $sx(y) $sx(m)] sy [list $sy(n) $sy(t) $sy(y) $sy(m)] \
        sz [list $sz(n) $sz(t) $sz(y) $sz(m)] sr [list $sr(n) $sr(t) $sr(y) $sr(m)] \
        sfr {} sbr {}]
}

set bad 0
foreach length {0.01 0.05 0.1 0.12} {
    set prof [make_profile $length]
    set r0 [::VMDHole::Mole::filter_bottleneck $prof 5.0 8.0 0.0]
    puts [format "  length=%.2f tolerance=0.0  -> %s" $length [expr {$r0 ? "keep" : "REJECT"}]]
    if {!$r0} { incr bad }
    set r1 [::VMDHole::Mole::filter_bottleneck $prof 5.0 8.0 2.0]
    puts [format "  length=%.2f tolerance=2.0  -> %s" $length [expr {$r1 ? "keep" : "REJECT"}]]
    if {!$r1} { incr bad }
}
if {$bad} {
    puts "  FAIL: $bad degenerate cases were rejected - MOLE keeps all of them"
    exit 1
}
puts "  PASS: every degenerate profile kept, matching MOLE's derived verdict"
exit 0
