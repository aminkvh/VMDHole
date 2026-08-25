# Load the pure-Tcl MOLE engine for a plain-tclsh test.
#
# The engine used to ship as eight mole_*.tcl modules beside vmdhole.tcl, and
# every test here sourced whichever subset it needed. It is now INLINED into
# vmdhole.tcl (the plugin ships as one script plus the compiled binaries), so
# those files no longer exist. This file replaces all of those `source` lines.
#
# vmdhole.tcl itself cannot simply be sourced here: it is ~1.5 MB of VMD/Tk
# plugin that calls VMD commands at load time, and this test group runs under
# a bare tclsh with no VMD. So extract exactly the sentinel-delimited engine
# region - which is self-contained ::VMDHole::Mole code, no VMD dependency -
# and source that.
#
# EVERY failure here is fatal on purpose. A silent empty extraction would let
# all 64 checks in test_mole_tcl_port.sh pass while testing nothing at all,
# which is the precise failure mode run_tests.sh's own header warns about: a
# green suite that could not catch the bugs it is named after.

namespace eval ::VMDHole::Mole {}

proc ::_load_inlined_mole_engine {} {
    set here [file dirname [file normalize [info script]]]
    set src  [file join $here .. vmdhole.tcl]
    if {![file readable $src]} {
        error "mole_tcl_engine: cannot read $src"
    }
    set fh [open $src r]
    set text [read $fh]
    close $fh

    set begin "# ===== BEGIN INLINED MOLE PURE-TCL ENGINE ====="
    set end   "# ===== END INLINED MOLE PURE-TCL ENGINE ====="

    set lines [split $text "\n"]
    set b {}; set e {}
    set i 0
    foreach ln $lines {
        if {[string trim $ln] eq $begin} { lappend b $i }
        if {[string trim $ln] eq $end}   { lappend e $i }
        incr i
    }
    # Exactly one of each: a duplicated or missing sentinel means the region is
    # not what this extractor thinks it is, and guessing would be worse than
    # stopping.
    if {[llength $b] != 1 || [llength $e] != 1} {
        error "mole_tcl_engine: expected exactly 1 BEGIN and 1 END sentinel in\
               vmdhole.tcl, found [llength $b] BEGIN and [llength $e] END -\
               the inlined MOLE engine block is missing or malformed"
    }
    set b [lindex $b 0]; set e [lindex $e 0]
    if {$e <= $b + 1} {
        error "mole_tcl_engine: inlined MOLE engine region is empty\
               (BEGIN at line [expr {$b+1}], END at line [expr {$e+1}])"
    }
    set body [join [lrange $lines [expr {$b+1}] [expr {$e-1}]] "\n"]

    set tmpdir [expr {[info exists ::env(TMPDIR)] ? $::env(TMPDIR) : "/tmp"}]
    set tmp [file join $tmpdir "vmdhole_mole_engine_[pid].tcl"]
    set oh [open $tmp w]
    puts $oh $body
    close $oh
    # Sourced (not eval'd) so a syntax error reports a real file and line.
    set code [catch {uplevel #0 [list source $tmp]} err opts]
    file delete -force $tmp
    if {$code} { return -options $opts $err }

    # The extraction can only be trusted if it actually produced the engine.
    if {![llength [info procs ::VMDHole::Mole::find_tunnels]]} {
        error "mole_tcl_engine: sourced the inlined region but\
               ::VMDHole::Mole::find_tunnels is still undefined -\
               extraction produced the wrong text"
    }
}

::_load_inlined_mole_engine
