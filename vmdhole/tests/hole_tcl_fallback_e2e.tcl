# The pure-Tcl HOLE fallback through the PLUGIN's own run path.
#
# hole_tcl_fallback.tcl checks the pieces; this checks the wiring - that with no
# `hole` binary configured, run_analysis routes to the fallback, the job pool
# runs it, and a real profile lands in the plugin's own results. The two halves
# were the failure mode worth guarding: the pieces all worked while the plugin
# still refused to run, because the fallback has to reproduce the per-frame FILE
# layout the pool's parser reads, not just replace an exec.
#
# The reference is computed HERE from the PDB the plugin itself wrote, not from
# a frozen number: HOLE is Monte Carlo and $sel writepdb is not the input file
# byte for byte, so any hard-coded radius would be measuring the fixture.

set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc note {msg} { puts "  ....  $msg" }

set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $root vmdhole 1GRM.pdb]
set RAD  [file join $root native stock_build hole2 rad simple.rad]
# Honour VMDHOLE_HOLE_EXE_DIR like the rest of the suite; ~/hole2/exe stays the
# install default. Hardcoding $HOME skipped this against any tree whose HOLE
# lives elsewhere, while run_tests.sh reported the group as passed.
set _exedir [expr {[info exists ::env(VMDHOLE_HOLE_EXE_DIR)] && $::env(VMDHOLE_HOLE_EXE_DIR) ne ""
                   ? $::env(VMDHOLE_HOLE_EXE_DIR) : [file join $::env(HOME) hole2 exe]}]
set REF  [file join $_exedir hole]

cd [file join $here ..]
source vmdhole.tcl

if {![file readable $PDB] || ![file readable $RAD]} {
    puts "  SKIP  no 1GRM fixture / radius file"
    puts "E2E-RESULT pass=$pass fail=$fail"
    quit
}
if {[::VMDHole::_hole_tcl_exe] eq ""} {
    puts "  SKIP  no tclsh for the fallback"
    puts "E2E-RESULT pass=$pass fail=$fail"
    quit
}

set work [file join [::VMDHole::get_temp_base] "vmdhole_e2e_[pid]"]
file mkdir $work

set mid [mol new $PDB waitfor all]
array set ::VMDHole::state [list \
    molid $mid frame_spec 0 selection all \
    hole_exec {} radius_file $RAD \
    cpoint {0 0 0} cvect {0 0 1} sample 0.5 endrad 8.0 \
    random_seed 1 pore_method circular display_mode none \
    work_dir $work keep_input_pdb 1 save_results 1 extra_cards {} ignore {}]

# 1. validate_inputs must ROUTE, not raise - a missing binary is only fatal when
#    the fallback is unusable.
if {[catch {::VMDHole::validate_inputs} verr]} {
    chk "validate_inputs routes to the fallback" "error: $verr" "no error"
} else {
    chk "validate_inputs routes to the fallback" $::VMDHole::state(hole_tcl_fallback) 1
}

# 2. The run itself, through the normal job pool.
set t0 [clock milliseconds]
if {[catch {::VMDHole::run_analysis} rerr]} {
    chk "run_analysis completes" "error: $rerr" "no error"
} else {
    note "run took [expr {[clock milliseconds]-$t0}] ms"
    chk "one result frame" [llength $::VMDHole::result_frames] 1
    if {[llength $::VMDHole::result_frames]} {
        set f [lindex $::VMDHole::result_frames 0]
        set p [dict get $::VMDHole::results $f profile]
        chk "profile is valid" [dict get $p valid] 1
        set pts [dict get $p points]
        chk "profile has rows" [expr {$pts > 20 ? 1 : 0}] 1
        note "points=$pts min_radius=[dict get $p min_radius] at [dict get $p min_coord]"

        # 3. Same answer as the binary, on the PDB the plugin actually wrote.
        # The per-frame directory map is local to run_analysis; the run's own
        # root is not, and the frame_%05d layout under it is the plugin's.
        set rd [file join $::VMDHole::state(last_root_dir) [format "frame_%05d" $f]]
        if {![file exists [file join $rd input_frame.pdb]]} { set rd "" }
        if {![file executable $REF]} {
            puts "  SKIP  no reference binary at $REF"
        } elseif {$rd eq ""} {
            chk "found the frame directory" 0 1
        } else {
            set cmp [file join $work cmp]
            file mkdir $cmp
            file copy -force [file join $rd input_frame.pdb] [file join $cmp in.pdb]
            file copy -force $RAD [file join $cmp s.rad]
            set ih [open [file join $cmp c.inp] w]
            puts $ih "coord  in.pdb\nradius s.rad\nsample 0.5\nendrad 8.0\nshorto 0"
            puts $ih "cpoint 0 0 0\ncvect  0 0 1\nraseed 1\nstop"
            close $ih
            set cwd [pwd]
            cd $cmp
            catch {exec $REF < c.inp > out.txt 2>@1}
            cd $cwd
            set rmin ""
            set rf [open [file join $cmp out.txt] r]
            set intab 0
            while {[gets $rf line] >= 0} {
                if {!$intab} {
                    if {[string match "*cenxyz.cvec*radius*" $line]} { set intab 1 }
                    continue
                }
                if {[regexp {^\s*[-+0-9.eE]+\s+([-+0-9.eE]+)\s+[-+0-9.eE]+\s+[-+0-9.eE]+} \
                         $line -> r]} {
                    if {$rmin eq "" || $r < $rmin} { set rmin $r }
                } elseif {$rmin ne ""} break
            }
            close $rf
            if {$rmin eq ""} {
                chk "reference profile parsed" 0 1
            } else {
                set d [expr {abs($rmin - [dict get $p min_radius])}]
                note "reference bottleneck $rmin, fallback [dict get $p min_radius]"
                # Same seed and same input, so this is print precision, not
                # HOLE's 0.0053 A rerun noise.
                chk "bottleneck matches the binary (|d| = [format %.3g $d])" \
                    [expr {$d < 1e-3 ? 1 : 0}] 1
            }
        }
    }
}

# 4. The WHOLE chain with no binaries at all: profile AND a real 3D surface.
#    display_mode triangulated is the one that needs both sph_process and
#    sos_triangle, so this is what proves they are wired rather than merely
#    ported. Low dot density on purpose - the chain is byte-identical at any
#    density (hole_tcl_fallback.tcl checks that) and ~50 s at the default 15.
set work2 [file join [::VMDHole::get_temp_base] "vmdhole_e2e_surf_[pid]"]
file mkdir $work2
array set ::VMDHole::state [list sph_process_exec {} sos_triangle_exec {} \
    display_mode triangulated dot_density 4 work_dir $work2 prebuild_surfaces 0]
if {[catch {::VMDHole::validate_inputs} verr]} {
    chk "validate_inputs allows a surface run with no binaries" "error: $verr" "no error"
} else {
    chk "sph_process routes to the fallback"  $::VMDHole::state(sph_process_exec)  {}
    chk "sos_triangle routes to the fallback" $::VMDHole::state(sos_triangle_exec) {}
    set t1 [clock milliseconds]
    if {[catch {::VMDHole::run_analysis} rerr]} {
        chk "surface run completes" "error: $rerr" "no error"
    } else {
        note "surface run took [expr {[clock milliseconds]-$t1}] ms"
        # The plugin's own two surface procs, on this run's .sph. Called
        # directly rather than through load_surface_for_frame, which also
        # renders into VMD - the wiring under test is the sph_process /
        # sos_triangle commands those procs build, not the drawing.
        set f2 [lindex $::VMDHole::result_frames 0]
        set sph [dict get $::VMDHole::results $f2 sph_file]
        set sos [file join $work2 e2e.sos]
        set plot [file join $work2 e2e.plot]
        if {![file exists $sph]} {
            chk "the run produced a .sph" 0 1
        } elseif {[catch {::VMDHole::run_sph_process $sph $sos 1 4} serr]} {
            chk "run_sph_process works with no binary" "error: $serr" "no error"
        } elseif {[catch {::VMDHole::run_sos_triangle $sos $plot} terr]} {
            chk "run_sos_triangle works with no binary" "error: $terr" "no error"
        } else {
            chk "the surface has real geometry" \
                [::VMDHole::surface_has_geometry $plot] 1
            note "surface [file size $plot] bytes from a [file size $sos]-byte .sos"

            # Tunnel mode writes its mesh job as a /bin/sh script instead of
            # exec'ing directly, so the generated TEXT has to be runnable - a
            # quoting slip there is invisible to every other check here.
            set msh [file join $work2 mesh.sh]
            set mh [open $msh w]
            puts $mh [::VMDHole::_surface_mesh_script 4 $sph \
                          [file join $work2 m.sos] [file join $work2 m.plot]]
            close $mh
            if {[catch {exec sh $msh} merr]} {
                chk "the tunnel mesh script runs" "error: $merr" "no error"
            } else {
                chk "the tunnel mesh script produces geometry" \
                    [::VMDHole::surface_has_geometry [file join $work2 m.plot]] 1
            }

            # display_mode dots takes a different branch: --points is an
            # accelerator-only flag, so with no binary it goes through the Tcl
            # dots_from_trinorm post-step instead.
            set ::VMDHole::state(display_mode) dots
            set dplot [file join $work2 dots.plot]
            if {[catch {::VMDHole::run_sos_triangle $sos $dplot} derr]} {
                chk "dots mode builds without a binary" "error: $derr" "no error"
            } else {
                chk "dots mode: fast --points path correctly unavailable" \
                    [::VMDHole::dots_fast_available] 0
                chk "dots mode produces geometry" \
                    [::VMDHole::surface_has_geometry $dplot] 1
            }
            set ::VMDHole::state(display_mode) triangulated
        }
    }
}
catch {file delete -force $work2}

# 5. The two non-spherical pore methods through run_analysis, still with no
#    binaries. hole_tcl_pore_methods.tcl checks their NUMBERS against the real
#    hole; this checks the plugin can be driven into them at all - the method
#    reaches the engine through a control-file card, so a translation that
#    dropped it would leave a perfectly valid SPHERICAL profile behind and
#    nothing here would look wrong.
#    Both build a surface, and by different routes: CAPSULE's .sph never reaches
#    sph_process at all (create_plot_asset stitches stadium cross-sections in
#    Tcl), while CONNOLLY's ~600-sphere cloud goes through the same
#    sph_process -> sos_triangle chain as the spherical case - cheap here only
#    because create_plot_asset trims and reduces the cloud first.
#    endrad 4.0, not the 8.0 used elsewhere: same code, a third of the slices.
foreach {meth dm} {capsule triangulated connolly triangulated} {
    set w3 [file join [::VMDHole::get_temp_base] "vmdhole_e2e_${meth}_[pid]"]
    file mkdir $w3
    array set ::VMDHole::state [list frame_spec 0 hole_exec {} sph_process_exec {} \
        sos_triangle_exec {} pore_method $meth display_mode $dm work_dir $w3 \
        endrad 4.0 dot_density 4 prebuild_surfaces 0]
    if {[catch {::VMDHole::validate_inputs} verr]} {
        chk "$meth: validate_inputs allows it with no binaries" "error: $verr" "no error"
    } elseif {[catch {::VMDHole::run_analysis} rerr]} {
        chk "$meth: run_analysis completes" "error: $rerr" "no error"
    } else {
        set f3 [lindex $::VMDHole::result_frames 0]
        set p3 [dict get $::VMDHole::results $f3 profile]
        chk "$meth: profile is valid" [dict get $p3 valid] 1
        chk "$meth: profile has rows" [expr {[dict get $p3 points] > 20 ? 1 : 0}] 1
        note "$meth: points=[dict get $p3 points] min_radius=[dict get $p3 min_radius]"
        set s3 [dict get $::VMDHole::results $f3 sph_file]
        chk "$meth: wrote a .sph" [expr {[file exists $s3] && [file size $s3] > 0}] 1
        if {$dm ne "none"} {
            if {[catch {::VMDHole::create_plot_asset [file dirname $s3] $s3 $dm \
                            $mid $f3} a3]} {
                chk "$meth: builds a surface" "error: $a3" "no error"
            } else {
                chk "$meth: builds real geometry" \
                    [expr {[dict get $a3 kind] eq "vmd_plot" &&
                           [::VMDHole::surface_has_geometry [dict get $a3 path]]}] 1
            }
        }
    }
    catch {file delete -force $w3}
}

catch {mol delete $mid}
catch {file delete -force $work}
puts "E2E-RESULT pass=$pass fail=$fail"
# MUST quit: `vmd -e` runs the script and then enters its interactive text
# prompt, which blocks on stdin. Standalone that stdin is closed and the test
# looked fine; under run_tests.sh's $(...) capture it inherits the suite's and
# hangs forever. headless_smoke.tcl ends the same way for the same reason.
quit
