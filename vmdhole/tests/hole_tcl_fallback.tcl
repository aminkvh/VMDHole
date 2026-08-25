# The pure-Tcl HOLE fallback, checked against the real binary's own profile.
#
# Runs under a bare tclsh: the fallback's plugin-side procs are pure Tcl, so
# they are lifted straight out of vmdhole.tcl rather than reimplemented here -
# a copy would drift and this test would stop covering the shipped code.
#
# What it proves, in order of what would actually break:
#   1. the sentinel extraction + appended driver form a script tclsh can run
#      (that is _hole_tcl_script's own self-check, exercised here);
#   2. every column of the profile TSV matches the binary's, including the
#      mid-point rows, cen_line_D and the running Sum ds/area, which are
#      reconstructed from the search rather than parsed from HOLE's output;
#   3. the .sph is byte-identical to the binary's;
#   4. a card the engine cannot honour REFUSES instead of being ignored.
#
# Skips (does not fail) without ~/hole2/exe/hole: there is nothing to compare to.

set HERE [file dirname [file normalize [info script]]]
set SRC  [file join $HERE .. vmdhole.tcl]
set REF  [file join $::env(HOME) hole2 exe hole]

set pass 0
set fail 0
proc ok {name} { global pass; incr pass; puts "  PASS $name" }
proc no {name why} { global fail; incr fail; puts "  FAIL $name: $why" }
proc chk_eq {name got want} {
    if {$got eq $want} { ok "$name = $got" } else { no $name "got $got, expected $want" }
}

if {![file executable $REF]} {
    puts "  SKIP hole_tcl_fallback: no reference binary at $REF"
    exit 0
}

# ---- lift the fallback's plugin-side procs out of vmdhole.tcl ----------------
# Accumulated with info complete rather than "until a line that is just }": the
# driver body is itself a braced block of procs that close at column 0.
namespace eval ::VMDHole {}
set fh [open $SRC r]; set text [read $fh]; close $fh
set lifted ""
set buf ""
foreach ln [split $text "\n"] {
    if {$buf eq "" && ![string match "proc ::VMDHole::_hole_tcl_*" $ln]} continue
    append buf $ln "\n"
    if {[info complete $buf]} { append lifted $buf; set buf "" }
}
if {[llength [regexp -all -inline {proc ::VMDHole::_hole_tcl_\w+} $lifted]] < 5} {
    puts "  FAIL hole_tcl_fallback: lifted only [regexp -all {proc ::VMDHole::_hole_tcl_} $lifted]\
          fallback procs from vmdhole.tcl"
    exit 1
}
eval $lifted

set TMP [file join [expr {[info exists ::env(TMPDIR)] ? $::env(TMPDIR) : "/tmp"}] \
             "vmdhole_fbtest_[pid]"]
file mkdir $TMP
proc ::VMDHole::get_temp_base {} { return $::TMP }
proc ::vmdcon {args} {}
set ::VMDHole::plugin_dir [file normalize [file join $HERE ..]]

# ---- 1. the engine extracts, the driver parses, the script self-checks -------
set script [::VMDHole::_hole_tcl_script]
if {$script eq ""} {
    no "engine script builds" "_hole_tcl_script returned empty"
    puts "hole_tcl_fallback: $pass passed, [incr fail 0] failed"
    exit 1
}
ok "engine script builds and passes its own self-check"

# ---- 2/3. one frame, fallback vs binary -------------------------------------
# The bundled tutorial structure - present in every clone, unlike the old
# paper/ copy this used to point at.
set PDB [file join $HERE .. 1GRM.pdb]
# The radius file ships as a FIXTURE so this suite is self-contained on a
# fresh clone (it is HOLE 2's own simple.rad, Apache-2.0 - the licence text
# already travels in native/LICENSE). The local build-tree copy is only a
# fallback for dev checkouts that predate the fixture.
set RAD [file join $HERE fixtures simple.rad]
if {![file exists $RAD]} {
    set RAD [file join $HERE .. .. native stock_build hole2 rad simple.rad]
}
if {![file readable $PDB] || ![file readable $RAD]} {
    puts "  SKIP profile comparison: no 1GRM fixture / radius file"
} else {
    file copy -force $PDB [file join $TMP in.pdb]
    file copy -force $RAD [file join $TMP s.rad]
    # The control file the plugin itself writes, minus the per-frame axis: same
    # cards, same order, so _hole_tcl_args_from_inp is parsing the real thing.
    set inp [file join $TMP hole.inp]
    set ih [open $inp w]
    puts $ih "coord  in.pdb"
    puts $ih "radius s.rad"
    puts $ih "sphpdb hole_out.sph"
    puts $ih "sample 0.5"
    puts $ih "endrad 8.0"
    puts $ih "shorto 0"
    puts $ih "cpoint 0 0 0"
    puts $ih "cvect  0 0 1"
    puts $ih "raseed 1"
    puts $ih "stop"
    close $ih

    lassign [::VMDHole::_hole_tcl_args_from_inp $inp] st args
    if {$st ne "ok"} {
        no "control file translates" $args
    } else {
        ok "control file translates to engine flags"
        # The engine takes bare filenames from the control file, exactly as the
        # job pool's run.sh does after cd'ing into the frame's directory.
        set tclsh [::VMDHole::_hole_tcl_exe]
        cd $TMP
        set rc [catch {exec {*}$tclsh $script {*}$args} out]
        if {$rc} {
            no "fallback runs" $out
        } else {
            ok "fallback runs one frame"
            # The binary, same cards, its own .sph under a different name.
            set rh [open [file join $TMP ref.inp] w]
            puts $rh "coord  in.pdb\nradius s.rad\nsphpdb ref.sph\nsample 0.5\nendrad 8.0"
            puts $rh "shorto 0\ncpoint 0 0 0\ncvect  0 0 1\nraseed 1\nstop"
            close $rh
            catch {exec $REF < [file join $TMP ref.inp] > [file join $TMP ref_out.txt] 2>@1}

            # Reference TSV via the SAME extraction the plugin's job pool uses,
            # so a difference here is a difference in the numbers, not in how
            # they were read.
            set rows {}
            set rf [open [file join $TMP ref_out.txt] r]
            set intab 0
            while {[gets $rf line] >= 0} {
                if {!$intab} {
                    if {[string match "*cenxyz.cvec*radius*" $line]} { set intab 1 }
                    continue
                }
                if {[regexp {^\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)} \
                         $line -> a b c d]} {
                    lappend rows [list $a $b $c $d]
                } elseif {[llength $rows]} break
            }
            close $rf

            set tf [open [file join $TMP hole_profile.tsv] r]
            gets $tf
            set mine {}
            while {[gets $tf line] >= 0} {
                set f [split $line "\t"]
                if {[llength $f] >= 4} { lappend mine [lrange $f 0 3] }
            }
            close $tf

            if {[llength $rows] == 0} {
                no "reference profile parsed" "no table in ref_out.txt"
            } elseif {[llength $rows] != [llength $mine]} {
                no "row count matches" "reference [llength $rows], fallback [llength $mine]"
            } else {
                # Numeric, never string: both sides print %.5f here, but a
                # string compare of formatted doubles is the classic way to
                # condemn a correct implementation.
                set worst {0 0 0 0}
                foreach r $rows m $mine {
                    for {set i 0} {$i < 4} {incr i} {
                        set d [expr {abs([lindex $r $i] - [lindex $m $i])}]
                        if {$d > [lindex $worst $i]} { lset worst $i $d }
                    }
                }
                lassign $worst dc dr dcl df
                # 1e-5 is the printed precision of both tables.
                foreach {nm d tol} [list coord $dc 1e-5 radius $dr 1e-4 \
                                         cen_line_d $dcl 1e-3 sum_s_over_area $df 1e-4] {
                    if {$d <= $tol} { ok "TSV $nm matches the binary (max |d| $d)" \
                    } else { no "TSV $nm matches the binary" "max |d| $d > $tol" }
                }
            }

            # .sph: byte-identical, which is what makes the surface pipeline
            # downstream see exactly what the binary would have produced.
            set a [file join $TMP hole_out.sph]
            set b [file join $TMP ref.sph]
            if {![file exists $a] || ![file exists $b]} {
                no ".sph written" "missing [expr {[file exists $a] ? {ref.sph} : {hole_out.sph}}]"
            } else {
                set ah [open $a r]; set at [read $ah]; close $ah
                set bh [open $b r]; set bt [read $bh]; close $bh
                if {$at eq $bt} { ok ".sph is byte-identical to the binary's" \
                } else { no ".sph is byte-identical" "files differ" }
            }
        }
    }
}

# ---- 3b. the surface stages, byte for byte ----------------------------------
# Run at dot density 4: the chain is byte-identical at any density but costs
# ~50 s at the plugin's default 15 and 0.6 s here, and what is under test is
# agreement, not throughput.
set refdir [file join $HERE .. .. hole_tcl reference_bin]
set sphbin [file join $refdir sph_process]
set tribin [file join $refdir sos_triangle]
foreach {v alt} [list sphbin sph_process tribin sos_triangle] {
    if {![file executable [set $v]]} { set $v [file join $::env(HOME) hole2 exe $alt] }
}
set sphin [file join $TMP hole_out.sph]
if {![file executable $sphbin] || ![file executable $tribin]} {
    puts "  SKIP surface comparison: no sph_process / sos_triangle binary"
} elseif {![file exists $sphin]} {
    puts "  SKIP surface comparison: the fallback run produced no .sph"
} else {
    set cwd [pwd]
    cd $TMP
    set tclsh [::VMDHole::_hole_tcl_exe]
    # -color is what the plugin always asks for (the radius-banded red/green/
    # blue surface); the uncolored single-pass form is checked too because the
    # two share the emit loop and only one was ported originally.
    foreach {tag cflag} {plain 0 color 1} {
        set bargs [list -sos -dotden 4]
        if {$cflag} { lappend bargs -color }
        catch {exec {*}[list $sphbin {*}$bargs $sphin b_$tag.sos] >/dev/null 2>@1}
        set rc [catch {exec {*}$tclsh $script --sph-process 4 $cflag $sphin t_$tag.sos} e]
        if {$rc} {
            no ".sos ($tag) built" $e
        } elseif {![file exists "b_$tag.sos"]} {
            no ".sos ($tag) reference built" "sph_process wrote nothing"
        } else {
            set bh [open b_$tag.sos r]; set bt [read $bh]; close $bh
            set th [open t_$tag.sos r]; set tt [read $th]; close $th
            if {$bt eq $tt} { ok ".sos ($tag) is byte-identical to sph_process" \
            } else { no ".sos ($tag) is byte-identical" "files differ" }
        }
    }
    if {[file exists [file join $TMP b_color.sos]]} {
        catch {exec sh -c "[list $tribin] -s < b_color.sos > b.plot 2>/dev/null"}
        set rc [catch {exec {*}$tclsh $script --sos-triangle b_color.sos t.plot} e]
        if {$rc} {
            no ".plot built" $e
        } else {
            set bh [open [file join $TMP b.plot] r]; set bt [read $bh]; close $bh
            set th [open [file join $TMP t.plot] r]; set tt [read $th]; close $th
            if {$bt eq $tt && [string length $bt] > 0} {
                ok ".plot is byte-identical to sos_triangle ([string length $bt] bytes)"
            } else {
                no ".plot is byte-identical" "differ, or empty ([string length $bt] bytes)"
            }
        }
    }
    cd $cwd
}

# ---- 4. a FAILED frame must keep its own diagnostic -------------------------
# The pool deletes the frame's tmpfs directory, and only copies hole_out.txt out
# of it when the TSV came back with no data rows. The awk path gets a header for
# free from its BEGIN block; the fallback branch has to pre-write one, or a
# driver that dies mid-search leaves no TSV at all, `wc -l` errors, the copy
# never happens and the frame vanishes together with the reason it failed.
set hdr_plugin ""
if {[regexp {printf '([^']*)' > hole_profile\.tsv} $text -> h]} { set hdr_plugin $h }
set hdr_driver ""
if {[regexp {puts \$fh "(coord\\tradius[^"]*)"} [::VMDHole::_hole_tcl_driver] -> h]} {
    set hdr_driver $h
}
if {$hdr_plugin eq "" || $hdr_driver eq ""} {
    no "failed-frame header is pre-written" "could not find the header in vmdhole.tcl"
} else {
    # Same header both sides, or the pre-written line is not the one the reader
    # would have got from a successful run.
    if {[string trimright $hdr_plugin "\\n"] eq $hdr_driver} {
        ok "pre-written header matches the driver's"
    } else {
        no "pre-written header matches the driver's" "'$hdr_plugin' vs '$hdr_driver'"
    }
    set fdir [file join $TMP faildir]
    file mkdir $fdir
    set eh [open [file join $fdir empty.pdb] w]; puts $eh "REMARK no atoms\nEND"; close $eh
    set cwd [pwd]
    cd $fdir
    exec sh -c "printf '$hdr_plugin' > hole_profile.tsv"
    set logf [file join $fdir hole_out.txt]
    # A REAL radius file, so the failure under test is the empty PDB and not the
    # missing-radius refusal checked below.
    catch {exec {*}[::VMDHole::_hole_tcl_exe] $script -pdb empty.pdb \
               -rad [file join $TMP s.rad] -sph out.sph -tsv hole_profile.tsv > $logf 2>@1}
    cd $cwd
    set n 0
    if {[file exists [file join $fdir hole_profile.tsv]]} {
        set th [open [file join $fdir hole_profile.tsv] r]
        set n [llength [split [string trimright [read $th] "\n"] "\n"]]
        close $th
    }
    chk_eq "a failed frame leaves a 1-line TSV so its log is kept" $n 1
    # A missing radius file must STOP the run. hole::read_rad_file returns an
    # empty ruleset for one, and every atom then silently falls back to
    # element_radius - a plausible profile with wrong radii, which is worse than
    # no profile. Checked here because the plugin's own ensure_path guarantee
    # does not reach inside the frame's subprocess.
    cd $fdir
    set radlog [file join $fdir rad_out.txt]
    catch {exec {*}[::VMDHole::_hole_tcl_exe] $script -pdb [file join $TMP in.pdb] \
               -rad no_such_file.rad -sph r.sph -tsv r.tsv > $radlog 2>@1}
    cd $cwd
    set rwhy ""
    if {[file exists $radlog]} { set rh [open $radlog r]; set rwhy [read $rh]; close $rh }
    if {[string match "*radius file not readable*" $rwhy] && ![file exists [file join $fdir r.tsv]]} {
        ok "a missing radius file stops the run instead of using element radii"
    } else {
        no "a missing radius file stops the run" \
           "log: [string range $rwhy 0 60] / tsv exists: [file exists [file join $fdir r.tsv]]"
    }
    set why ""
    if {[file exists $logf]} { set lh [open $logf r]; set why [read $lh]; close $lh }
    if {[string match "*no atoms*" $why]} {
        ok "the kept log names the real cause"
    } else {
        no "the kept log names the real cause" "log was: [string range $why 0 80]"
    }
}

# ---- 5. cards the engine cannot honour must REFUSE, not be ignored ----------
# CONNOLLY and CAPSULE used to be here; both are wired now, and
# hole_tcl_pore_methods.tcl checks their tables against the binary. What must
# still refuse is anything the engine has no implementation of at all -
# translating one of those silently answers a question the user did not ask.
foreach {name card} {unknown "2dmaps 1" bad-conn "conn wide"
                     both-methods "capsule\nconn"} {
    set p [file join $TMP refuse.inp]
    set h [open $p w]
    puts $h "coord  in.pdb\nradius s.rad\nsphpdb hole_out.sph\nsample 0.5\nendrad 8.0"
    puts $h $card
    puts $h "stop"
    close $h
    lassign [::VMDHole::_hole_tcl_args_from_inp $p] st why
    if {$st eq "refused"} { ok "refuses $name ($why)" \
    } else { no "refuses $name" "translated it instead of refusing" }
}

# The pore-method cards must now TRANSLATE, and carry their parameters through:
# a `conn 1.4 0.8` that quietly became a default-probe run would be a wrong
# answer that looks right.
foreach {card want} {conn {-method connolly}
                     connolly {-method connolly}
                     capsule {-method capsule}
                     capsul {-method capsule}
                     {conn 1.4 0.8} {-method connolly -probe 1.4 -grid 0.8}
                     {conn 1.4} {-method connolly -probe 1.4}
                     {ignore HOH SOL} {-method spherical -ignore {HOH SOL}}} {
    set p [file join $TMP method.inp]
    set h [open $p w]
    puts $h "coord  in.pdb\nradius s.rad\nsphpdb hole_out.sph\nsample 0.5\nendrad 8.0"
    puts $h $card
    puts $h "stop"
    close $h
    lassign [::VMDHole::_hole_tcl_args_from_inp $p] st args
    set got {}
    foreach k {-method -probe -grid -ignore} {
        set i [lsearch -exact $args $k]
        if {$i >= 0} { lappend got $k [lindex $args [expr {$i+1}]] }
    }
    if {$st eq "ok" && $got eq $want} { ok "'$card' translates to $want" \
    } else { no "'$card' translates to $want" "got $st / $got" }
}

cd $HERE
file delete -force $TMP
puts "hole_tcl_fallback: $pass passed, $fail failed"
if {$fail} { exit 1 }
