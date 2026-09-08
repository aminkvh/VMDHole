# The plugin's reader must recover CAVITIES from the engine's file.
#
# The companion of test_mole_lining_parse.tcl: that one checks the tunnel
# records, this one the cavity block the engine added alongside them (V / VB /
# VI / VP). A cavity read from the wrong column is invisible to a byte
# comparison of the file and fatal to the Cavities window, which is the only
# thing that displays it.
#
# Also covers the two per-tunnel fields the same parser gained: the per-POINT
# FreeRadius/BRadius (P records, fields 6-7) that the MOLE-layout CSV export
# fills its columns from, and npos/nneg in the weighted properties.
#
# Usage: test_mole_cavity_parse.tcl ENGINE_OUT.txt
package provide Tk 8.5
set here [file dirname [info script]]
source [file join $here .. vmdpathfinder.tcl]

set fails 0
proc report {label ok {detail ""}} {
    global fails
    if {!$ok} { incr fails }
    puts [format "  %-52s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
}

set fh [open [lindex $argv 0] r]; set txt [read $fh]; close $fh
set L [::VMDPathFinder::_tunnel_parse_lining $txt]

# What the FILE says, counted independently of the parser.
set nV 0; set nVP 0; set vp_per {}
array set vpc {}
foreach line [split $txt "\n"] {
    set f [split [string trim $line]]
    switch -exact -- [lindex $f 0] {
        V  { incr nV }
        VP { incr nVP; set id [lindex $f 1]
             set vpc($id) [expr {[info exists vpc($id)] ? $vpc($id)+1 : 1}] }
    }
}

# NOT a skip. The only automated caller (test_mole_tcl_port.sh) hands this the
# ground-truth run of a fixture with nine cavities, so an engine that emitted
# none is a regression in the writer, and going quiet about it would leave the
# reader untested while still reporting a pass to the group above.
report "the engine emitted cavity records at all" [expr {$nV > 0}] \
    "(no V lines: this build's writer, or a file from one predating them)"
if {$nV == 0} {
    puts "  ---- cavity-reader checks, 1 failed"
    exit 1
}

# 1. every V record becomes a cavity, ids ascending
set ::VMDPathFinder::tunnel_lining(0) $L
set cavs [::VMDPathFinder::_tunnel_cavities 0]
report "every V record is read back as a cavity" \
    [expr {[dict size $cavs] == $nV}] "(file $nV, parsed [dict size $cavs])"
report "cavity ids come back ascending" \
    [expr {[dict keys $cavs] eq [lsort -integer [dict keys $cavs]]}]

# 2. fields land in the right columns, and each cavity keeps its own spheres
set bad {}
set tot_sph 0
dict for {id cv} $cavs {
    if {[dict get $cv type] ni {Cavity Void}} { lappend bad "cav $id type=[dict get $cv type]" }
    if {![string is double -strict [dict get $cv volume]] || [dict get $cv volume] <= 0} {
        lappend bad "cav $id volume=[dict get $cv volume]"
    }
    if {![string is integer -strict [dict get $cv depth]]} { lappend bad "cav $id depth=[dict get $cv depth]" }
    # MOLE calls a cavity with no boundary residues a Void.
    set nb [dict get $cv nboundary]
    if {[dict get $cv type] eq "Void" && $nb != 0} { lappend bad "cav $id Void with $nb boundary" }
    if {[dict get $cv type] eq "Cavity" && $nb == 0} { lappend bad "cav $id Cavity with no boundary" }
    # the residue LISTS must match the counts the V line declared
    if {[llength [dict get $cv bres]] != $nb} {
        lappend bad "cav $id: $nb boundary declared, [llength [dict get $cv bres]] parsed"
    }
    if {[llength [dict get $cv ires]] != [dict get $cv ninner]} {
        lappend bad "cav $id: [dict get $cv ninner] inner declared, [llength [dict get $cv ires]] parsed"
    }
    set ns [llength [dict get $cv spheres]]
    incr tot_sph $ns
    set want [expr {[info exists vpc($id)] ? $vpc($id) : 0}]
    if {$ns != $want} { lappend bad "cav $id: $want VP lines, $ns spheres parsed" }
    foreach s [dict get $cv spheres] {
        if {[llength $s] != 4} { lappend bad "cav $id: sphere '$s' is not x y z r"; break }
        if {[lindex $s 3] <= 0} { lappend bad "cav $id: non-positive radius [lindex $s 3]"; break }
    }
}
report "cavity fields, residue counts and spheres are consistent" \
    [expr {[llength $bad] == 0}] "([join [lrange $bad 0 2] {; }])"
report "every VP line is assigned to its own cavity" \
    [expr {$tot_sph == $nVP}] "(file $nVP, parsed $tot_sph)"

# 3. the boundary/inner property blocks are read, in MOLE's own column order
set p_bad {}
dict for {id cv} $cavs {
    foreach which {bprops iprops} {
        set p [dict get $cv $which]
        if {![dict size $p]} { continue }
        foreach k {charge ionizable npos nneg hydropathy hydrophobicity polarity
                   logp logd logs mutability} {
            if {![dict exists $p $k]} { lappend p_bad "cav $id $which missing $k"; break }
            if {![string is double -strict [dict get $p $k]]} {
                lappend p_bad "cav $id $which $k=[dict get $p $k]"; break
            }
        }
        # npos/nneg are counts, and charge cannot exceed what they allow.
        if {[dict exists $p npos] && [dict exists $p nneg]} {
            set c [dict get $p charge]
            if {$c > [dict get $p npos] || $c < -[dict get $p nneg]} {
                lappend p_bad "cav $id $which charge $c vs +[dict get $p npos]/-[dict get $p nneg]"
            }
        }
    }
}
report "cavity property blocks parse in MOLE's column order" \
    [expr {[llength $p_bad] == 0}] "([join [lrange $p_bad 0 2] {; }])"

# 4. per-POINT FreeRadius/BRadius, one per P record of that tunnel
set pr_bad {}
array set pcount {}
foreach line [split $txt "\n"] {
    set f [split [string trim $line]]
    if {[lindex $f 0] ne "P" || [llength $f] < 8} continue
    set id [lindex $f 1]
    set pcount($id) [expr {[info exists pcount($id)] ? $pcount($id)+1 : 1}]
}
foreach id [array names pcount] {
    if {![dict exists $L $id.pfree]} { lappend pr_bad "T$id: no per-point freeradius"; continue }
    if {[llength [dict get $L $id.pfree]] != $pcount($id)} {
        lappend pr_bad "T$id: $pcount($id) P records, [llength [dict get $L $id.pfree]] freeradii"
    }
    if {[llength [dict get $L $id.pbr]] != $pcount($id)} {
        lappend pr_bad "T$id: $pcount($id) P records, [llength [dict get $L $id.pbr]] bradii"
    }
}
report "per-point FreeRadius/BRadius kept, one per route point" \
    [expr {[llength $pr_bad] == 0}] "([join [lrange $pr_bad 0 2] {; }])"

# 5. npos/nneg reach the weighted properties the lining window prints
set w_bad {}
foreach k [dict keys $L] {
    if {![string match "*.wprops" $k]} continue
    foreach f {npos nneg} {
        if {![dict exists [dict get $L $k] $f]} { lappend w_bad "$k missing $f" }
    }
}
report "npos/nneg reach the weighted properties" \
    [expr {[llength $w_bad] == 0}] "([join [lrange $w_bad 0 2] {; }])"

puts "  ---- cavity-reader checks, $fails failed"
exit [expr {$fails ? 1 : 0}]
