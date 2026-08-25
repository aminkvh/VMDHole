# CONNOLLY and CAPSULE through the pure-Tcl fallback.
#
# The .sph writers for both were ported and checked long before the fallback
# existed; what this covers is the PROFILE TABLE, which the fallback rebuilds
# from the search rather than parsing, and which differs per method:
#   - CAPSULE is hcapgr.f, not hograp.f - the row centre is the capsule
#     mid-point, cen.line.dis and integ.s/(area) both RESTART at the +ve block,
#     and integ takes a whole sample step on sampled rows only.
#   - CONNOLLY adds Requiv / Conn_s/Area / Requiv_estim, filled on sampled rows
#     and left BLANK on mid-points. _resolve_conn_radii reads that blankness as
#     "interpolate here", so the gap pattern is as load-bearing as the numbers.
#
# The reference is the TSV the plugin's OWN parser writes from the binary's
# printout - the two artifacts the plugin actually consumes - not a
# hand-rolled parse, which would only test this file against itself.

set pass 0; set fail 0
proc chk {name got want} {
    global pass fail
    if {$got eq $want} { incr pass; puts "  PASS  $name = $got" } \
    else { incr fail; puts "  FAIL  $name = $got (expected $want)" }
}
proc note {msg} { puts "  ....  $msg" }
proc done {} {
    global pass fail
    puts "POREMETHOD-RESULT pass=$pass fail=$fail"
    quit
}

set here [expr {[info exists ::env(VMDHOLE_TEST_DIR)] && $::env(VMDHOLE_TEST_DIR) ne ""
                ? $::env(VMDHOLE_TEST_DIR) : [pwd]}]
set root [file normalize [file join $here .. ..]]
set PDB  [file join $root vmdhole 1GRM.pdb]
set RAD  [file join $root native stock_build hole2 rad simple.rad]
set REF  [file join $::env(HOME) hole2 exe hole]

cd [file join $here ..]
source vmdhole.tcl

if {![file readable $PDB] || ![file readable $RAD]} { puts "  SKIP  no 1GRM fixture"; done }
if {![file executable $REF]} { puts "  SKIP  no reference binary at $REF"; done }
set exe [::VMDHole::_hole_tcl_exe]
set scr [::VMDHole::_hole_tcl_script]
if {$exe eq "" || $scr eq ""} { puts "  SKIP  no tclsh for the fallback"; done }

# One unit in HOLE's own last printed place: the plain columns are F12.5, the
# three CONNOLLY ones F12.3. Anything real is orders of magnitude bigger - the
# capsule mid-point axis rule, sabotaged, moves the radius column by 0.058.
set TOL {2e-5 2e-5 2e-5 2e-5 2e-3 2e-3 2e-3 2e-5}
set base [file join [::VMDHole::get_temp_base] "vmdhole_pm_[pid]"]

# conn-probe is not a fourth method but the one way a CONNOLLY parameter can
# reach the engine: write_control_file emits a bare `conn`, so a non-default
# probe only ever arrives by the user typing "conn 1.4 0.8" into Extra HOLE
# cards, which suppresses the generated card and flows through verbatim. A
# dropped -probe would run the 1.15 default and produce a profile that looks
# entirely reasonable, so the columns are compared AND the run is required to
# differ from the default-probe one below - otherwise a -probe ignored on both
# sides would pass this vacuously.
#
# ignore-trp is the same shape for the IGNORE card, which is a first-class
# sidebar field rather than an extra card: tsatr.f drops the atom before the
# search, the fallback drops the PDB line, and both must land on the same
# profile. TRP is the pore lining in gramicidin, so dropping it moves the
# bottleneck 0.208 -> 0.534 A and the row count 40 -> 43 - a card silently
# dropped on both sides could not fake that.
set requiv_by_tag [dict create]
set radius_by_tag [dict create]

# guess-* omit CPOINT and CVECT, which is the PLUGIN'S OWN DEFAULT - both
# sidebar fields start empty, so a control file with neither card is the
# ordinary case, not an edge case. HOLE then guesses both (hole.f:623 ->
# cguess.f). Pinning the axis in every case hid a real defect: the fallback
# guessed CPOINT but hard-defaulted CVECT to "0 0 1", so on any structure whose
# pore HOLE calls x or y it silently profiled a different pore. On 1GRM HOLE
# guesses Y, and the fallback returned 94 rows over 23.0 A against the binary's
# 154 over 38.0 A - a perpendicular channel, reported without a warning.
foreach {tag card want_method} {capsule capsule capsule
                                guess-capsule capsule capsule
                                guess-spherical {} spherical
                                connolly conn connolly
                                conn-probe {conn 1.4 0.8} connolly
                                ignore-trp {ignore TRP} spherical
                                spherical {} spherical} {
    set W [file join $base $tag]
    file mkdir $W
    file copy -force $PDB [file join $W in.pdb]
    file copy -force $RAD [file join $W s.rad]
    set fh [open [file join $W hole.inp] w]
    puts $fh "coord in.pdb\nradius s.rad\nsphpdb ref.sph\nsample 0.5\nendrad 8.0\nshorto 0"
    if {[string match "guess-*" $tag]} {
        puts $fh "raseed 1"
    } else {
        puts $fh "cpoint 0 0 0\ncvect 0 0 1\nraseed 1"
    }
    if {$card ne ""} { puts $fh $card }
    puts $fh "stop"
    close $fh

    lassign [::VMDHole::_hole_tcl_args_from_inp [file join $W hole.inp]] st args
    chk "$tag: the control file translates" $st ok
    if {$st ne "ok"} { continue }
    chk "$tag: routed to the right engine" [lindex $args [expr {[lsearch $args -method]+1}]] $want_method

    set cwd [pwd]
    cd $W
    catch {exec $REF < hole.inp > hole_out.txt 2>@1}
    cd $cwd
    set p [::VMDHole::parse_profile [file join $W hole_out.txt] [file join $W ref.tsv]]
    if {![dict get $p valid]} {
        chk "$tag: the binary produced a profile" [dict get $p message] "a profile"
        continue
    }

    set i [expr {[lsearch $args -sph] + 1}]
    set args [lreplace $args $i $i fb.sph]   ;# keep the binary's own .sph to diff
    cd $W
    set t0 [clock milliseconds]
    set err [catch {exec $exe $scr {*}$args} out]
    cd $cwd
    if {$err} { chk "$tag: the fallback ran" "error: $out" "no error"; continue }
    note "$tag: fallback took [expr {[clock milliseconds]-$t0}] ms"

    chk "$tag: .sph byte-identical to the binary" \
        [expr {![catch {exec cmp -s [file join $W ref.sph] [file join $W fb.sph]}]}] 1

    set f [open [file join $W ref.tsv]]; set a [split [string trim [read $f]] \n]; close $f
    set f [open [file join $W hole_profile.tsv]]; set b [split [string trim [read $f]] \n]; close $f
    chk "$tag: same number of profile rows" [llength $b] [llength $a]
    if {[llength $a] != [llength $b]} { continue }

    set worst {0 0 0 0 0 0 0 0}
    set blank 0
    for {set r 1} {$r < [llength $a]} {incr r} {
        set ra [split [lindex $a $r] \t]
        set rb [split [lindex $b $r] \t]
        for {set c 0} {$c < 8} {incr c} {
            set va [lindex $ra $c]; set vb [lindex $rb $c]
            if {$va eq "" && $vb eq ""} continue
            if {$va eq "" || $vb eq ""} { incr blank; continue }
            set d [expr {abs($va - $vb)}]
            if {$d > [lindex $worst $c]} { lset worst $c $d }
        }
    }
    note "$tag: worst |diff| per column = $worst"
    chk "$tag: blank cells fall in the same places" $blank 0
    foreach col {coord radius cen_line_d sum_s_over_area requiv conn_s_over_area requiv_estim cap_rad} \
            d $worst t $TOL {
        chk "$tag: $col agrees (|d| < $t)" [expr {$d < $t ? 1 : 0}] 1
    }
    dict set requiv_by_tag $tag [lmap r [lrange $b 1 end] {lindex [split $r \t] 4}]
    dict set radius_by_tag $tag [lmap r [lrange $b 1 end] {lindex [split $r \t] 1}]
}

# The probe actually reached the search: a bigger Connolly probe reaches less of
# the surface, so its Requiv column must not be the default run's.
if {[dict exists $requiv_by_tag connolly] && [dict exists $requiv_by_tag conn-probe]} {
    chk "a non-default conn probe changes the answer" \
        [expr {[dict get $requiv_by_tag connolly] ne
               [dict get $requiv_by_tag conn-probe] ? 1 : 0}] 1
}
# Same guard for IGNORE: 40 fewer atoms must change the radius profile, or the
# agreement above would only mean the card was dropped on BOTH sides.
if {[dict exists $radius_by_tag spherical] && [dict exists $radius_by_tag ignore-trp]} {
    chk "an ignore list changes the answer" \
        [expr {[dict get $radius_by_tag spherical] ne
               [dict get $radius_by_tag ignore-trp] ? 1 : 0}] 1
}

catch {file delete -force $base}
done
