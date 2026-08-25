# MOLE's atom/residue classification, C against Tcl, name by name.
#
# This is the one piece of the port with NO reference coverage. Every other
# check runs through 1tqn, 1BL8, 2ACE or 1AKD - all protein, all single-model,
# none with hydrogens and none with nucleic acid - so the H entry and the whole
# nucleic half of PdbEx.backboneNames are code that the ground-truth runs never
# touch. MOLE cannot be re-run here (no mono), so a parity check against the C
# plus the literal table is what is available, and it is worth having: the two
# engines each hold their own copy of the set, and a copy is a thing that drifts.
#
# Usage: test_mole_classify.tcl C_OUTPUT.txt
#   where C_OUTPUT.txt has "name backbone" / "resn amino" lines from the C.

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

# PdbEx.backboneNames, transcribed from the C# a second time on purpose - if
# both engines drifted the same way this would still catch it.
set BB {C N O H CA P O1P O2P OP1 OP2 O5' C5' C4' O4' C1' C2' C3' O3' O2'}
set AA {ALA ARG ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP
        TYR VAL ASN}
# Names that must NOT be backbone. CB/CG/OG1 are the side-chain atoms the
# lining split turns on; "CA " as calcium is the same string as the alpha
# carbon and MUST still classify as backbone here - the vdW table is where the
# element meaning lives, not this test.
set NOTBB {CB CG CG1 CG2 OG OG1 CD CD1 NE NZ OH SD FE ZN MG HOH XX}
set NOTAA {HEM HOH TIP3 WAT SOL DA DC DG DT U ZN MSE UNK}

set nbb 0; set naa 0
foreach n $BB { if {[::VMDHole::Mole::is_backbone_name $n]} { incr nbb } }
report "all [llength $BB] backbone names classify backbone" \
       [expr {$nbb == [llength $BB]}] "(got $nbb)"
foreach n $NOTBB { if {[::VMDHole::Mole::is_backbone_name $n]} { lappend badbb $n } }
report "side chains, metals and waters are not backbone" \
       [expr {![info exists badbb]}] "([expr {[info exists badbb] ? $badbb : {}}])"
foreach r $AA { if {[::VMDHole::Mole::is_amino_name $r]} { incr naa } }
report "all 20 amino names classify amino" [expr {$naa == 20}] "(got $naa)"
foreach r $NOTAA { if {[::VMDHole::Mole::is_amino_name $r]} { lappend badaa $r } }
report "HET, water, nucleotides and MSE are not amino" \
       [expr {![info exists badaa]}] "([expr {[info exists badaa] ? $badaa : {}}])"

# Case-insensitive: MOLE's sets use StringComparer.OrdinalIgnoreCase, and VMD
# hands back whatever the file had.
report "classification is case-insensitive" \
       [expr {[::VMDHole::Mole::is_backbone_name "ca"] &&
              [::VMDHole::Mole::is_backbone_name "o5'"] &&
              [::VMDHole::Mole::is_amino_name "his"] &&
              ![::VMDHole::Mole::is_amino_name "hem"]}]

# H is in MOLE's backbone set. Easy to "fix" away as a typo; it is not one, and
# it changes FreeRadius on every hydrogen-bearing structure - which in VMD means
# every MD trajectory.
report "H is a backbone name (MOLE's set, not a typo)" \
       [::VMDHole::Mole::is_backbone_name H]

# The C's answers for the same names.
if {[llength $argv] > 0 && [file readable [lindex $argv 0]]} {
    set bad {}
    set ncheck 0
    set fh [open [lindex $argv 0] r]
    foreach line [split [read $fh] "\n"] {
        set f [split [string trim $line]]
        if {[llength $f] != 3} continue
        lassign $f kind name val
        set mine [expr {$kind eq "BB" ? [::VMDHole::Mole::is_backbone_name $name]
                                      : [::VMDHole::Mole::is_amino_name $name]}]
        if {$mine != $val} { lappend bad "$kind $name c=$val tcl=$mine" }
        incr ncheck
    }
    close $fh
    report "C and Tcl agree on all $ncheck names" [expr {![llength $bad]}] \
           "([join [lrange $bad 0 3] {; }])"
} else {
    puts "  SKIP  no C output given"
}
exit [expr {$fails ? 1 : 0}]
