# The plugin's reader must recover the lining from the engine's file.
#
# Everything else checks the file the engine WRITES. This checks the other side
# of that interface: _tunnel_parse_lining is what any panel would call, and a
# field read from the wrong column is invisible to a byte-comparison of the file
# and fatal to whatever displays it. The expected values are MOLE 2's own
# tunnels.csv, so this is a file -> parser -> MOLE round trip, not a self-check.
#
# Usage: test_mole_lining_parse.tcl ENGINE_OUT.txt MOLE_TUNNELS.csv
package provide Tk 8.5
set here [file dirname [info script]]
source [file join $here .. vmdhole.tcl]

set fails 0
proc report {label ok {detail ""}} {
    global fails
    if {!$ok} { incr fails }
    puts [format "  %-52s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
}

set fh [open [lindex $argv 0] r]; set txt [read $fh]; close $fh
set L [::VMDHole::_tunnel_parse_lining $txt]

set fh [open [lindex $argv 1] r]; set csv [split [string trim [read $fh]] "\n"]; close $fh
# Id,Length,Charge,Ionizable,Hydropathy,Hydrophobicity,Polarity,LogP,LogD,LogS,Mutability
set cols {charge 2 ionizable 3 hydropathy 4 hydrophobicity 5 polarity 6
          logp 7 logd 8 logs 9 mutability 10}

set id 0
set bad {}
foreach row [lrange $csv 1 end] {
    incr id
    set f [split $row ","]
    if {![dict exists $L $id.wprops]} { lappend bad "T$id: no wprops"; continue }
    set w [dict get $L $id.wprops]
    foreach {k i} $cols {
        # MOLE prints two decimals; the engine writes four.
        if {abs([dict get $w $k] - [lindex $f $i]) > 0.005000001} {
            lappend bad "T$id $k [lindex $f $i] vs [dict get $w $k]"
        }
    }
}
report "parsed weighted properties match MOLE's tunnels.csv" \
       [expr {![llength $bad]}] "([join [lrange $bad 0 3] {; }])"
report "parser found all $id tunnels" [expr {$id > 0}]

# Layers and flow must survive the round trip with their structure intact, not
# merely their count: a residue token is "resn:seq:chain:backbone:flowindex" and
# an off-by-one in that split is exactly the kind of thing a count check misses.
set nl 0; set structok 1
for {set i 1} {$i <= $id} {incr i} {
    set ls [dict get $L $i.layers]
    incr nl [llength $ls]
    foreach y $ls {
        foreach r [dict get $y residues] {
            if {![string is integer -strict [dict get $r resid]]
                || [dict get $r backbone] ni {0 1}
                || ![string is integer -strict [dict get $r flow]]
                || [dict get $r flow] >= [llength [dict get $L $i.flow]]
                || [string length [dict get $r resname]] < 1} { set structok 0 }
        }
        if {[dict get $y end] < [dict get $y start]} { set structok 0 }
    }
}
report "every residue token parses into a well-formed record" $structok
report "layers recovered ($nl across $id tunnels)" [expr {$nl > 0}]
exit [expr {$fails ? 1 : 0}]
