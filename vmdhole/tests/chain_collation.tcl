# Chain-identifier order against MOLE's own collation.
#
# MOLE pins CultureInfo.InvariantCulture (Program.cs:35), so every
# OrderBy(ChainIdentifier) collates rather than compares ordinally. The expected
# string below is not derived - it is what .NET printed under MOLE's own mono
# for OrderBy on exactly this list:
#
#     ids.OrderBy(x => x)  ->  ?,_,0,1,9,a,A,AA,aB,Ab,b,B,c,C,z,Z
#
# Three rules, each of which strcmp gets wrong: punctuation before digits before
# letters; case-insensitive primary; lowercase before uppercase as a tiebreak
# applied only after the whole string compares equal ("aB" before "Ab").
#
# The C's mole_chain_cmp is checked against the same string by the shell wrapper,
# so the two engines are pinned to MOLE rather than to each other.
#
# Usage: chain_collation.tcl
set here [file dirname [info script]]
# The engine is inlined in vmdhole.tcl now (one shipped script + the
# compiled binaries), so this extracts it instead of sourcing modules.
source [file join [file dirname [file normalize [info script]]] mole_tcl_engine.tcl]

set want "?,_,0,1,9,a,A,AA,aB,Ab,b,B,c,C,z,Z"
set got [join [lsort -command ::VMDHole::Mole::_chain_cmp \
                   {A a B b C c Z z 0 1 9 ? _ AA Ab aB}] ","]
if {$got eq $want} {
    puts "  Tcl chain collation = MOLE's InvariantCulture order    PASS"
    exit 0
}
puts "  Tcl chain collation = MOLE's InvariantCulture order    FAIL"
puts "     want $want"
puts "     got  $got"
exit 1
