# Chain-collation regression for two sites a fixture-only check cannot
# distinguish: ResidueFlow (_build_flow) and HetResidues (_het_cmp via lsort). Mirrors
# native/mole/mole_lining_regression_test.c - see that file for the case
# derivations (1BL8's THR 107 A/C/D relabelled to a/C/D; 1ERI's DA 6/7 B
# relabelled DA 6 to lowercase b).
set here [file dirname [info script]]
# The engine is inlined in vmdhole.tcl now (one shipped script + the
# compiled binaries), so this extracts it instead of sourcing modules.
source [file join [file dirname [file normalize [info script]]] mole_tcl_engine.tcl]

set bad 0

# ---- Case 1: ResidueFlow ----------------------------------------------------
# restab entries are [name seq chain]; layer dict is {res {..} bb {..}}.
set restab {
    {THR 107 D}
    {THR 107 C}
    {THR 107 a}
}
set layers [list [dict create res {0 1 2} bb {0 0 0}]]
set flow [::VMDHole::Mole::_build_flow layers $restab]
set got {}
foreach f $flow {
    lassign [lindex $restab [lindex $f 0]] rn sq ch
    lappend got "$rn$sq$ch"
}
puts "  ResidueFlow: got [llength $flow] entries, order $got"
if {$got ne {THR107a THR107C THR107D}} {
    puts "  FAIL: want order a, C, D (InvariantCulture) - THR 107 A/C/D relabelled a/C/D"
    incr bad
}

# ---- Case 2: HetResidues -----------------------------------------------------
set hrestab {
    {DA 6 b}
    {DA 7 B}
}
set hout [lsort -command [list ::VMDHole::Mole::_het_cmp $hrestab] {0 1}]
lassign [lindex $hrestab [lindex $hout 0]] n0 s0 c0
lassign [lindex $hrestab [lindex $hout 1]] n1 s1 c1
puts "  HetResidues: got order $n0$s0$c0, $n1$s1$c1"
if {$hout ne {0 1}} {
    puts "  FAIL: want DA6 b, DA7 B unchanged (lowercase sorts before uppercase)"
    incr bad
}

if {$bad} { puts "  FAIL"; exit 1 }
puts "  PASS"
exit 0
