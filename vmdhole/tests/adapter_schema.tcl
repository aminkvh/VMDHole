# The VMD atom-table adapter, RUN rather than grepped.
#
# _tunnel_atoms_mole is the only producer of the engine's input in the shipped
# plugin, so a schema fault there is invisible to every fixture-based check. A
# grep for the format string cannot see altLoc and insertion code swapped, or
# either emitted as a placeholder; this loads a structure carrying both and
# reads the values back.
set LOG [open $::env(ADAPTER_LOG) w]
proc say {s} { global LOG; puts $LOG $s; flush $LOG }
if {[catch {
    uplevel #0 [list source $::env(VMDHOLE_TCL)]
    set mid [mol new $::env(ADAPTER_PDB) waitfor all]
    set ::VMDHole::state(selection) "all"
    set rows [::VMDHole::_tunnel_atoms_mole $mid 0]
    set fails 0
    proc chk {label ok {detail ""}} {
        global fails
        if {!$ok} { incr fails }
        say [format "  %-52s %s" $label [expr {$ok ? "PASS" : "FAIL $detail"}]]
    }
    chk "adapter emitted [llength $rows] rows" [expr {[llength $rows] == 13}]

    # A missing element used to be taken from the RESIDUE name's first letter -
    # "ALA" gave element "A", not carbon, so MOLE applied the wrong vdW radius
    # and the triangulation moved with it. Derived from the ATOM NAME now, with
    # two-letter elements resolved before the single-letter fallback.
    foreach {nm rn want} {CA ALA C  CA CA Ca  N ALA N  OD1 ASP O  FE HEM Fe
                          CL CL Cl  ZN ZN Zn  NA NA Na  NA ARG N  1HB ALA H} {
        chk "element: name '$nm' in residue '$rn' -> $want" \
            [expr {[::VMDHole::_element_from_atom_name $nm $rn] eq $want}] \
            "(got [::VMDHole::_element_from_atom_name $nm $rn])"
    }

    # "all" must mean ALL. The adapter used to rewrite an exact selection of
    # "all" into "protein", silently dropping nucleic acids, glycans, cofactors,
    # ligands, metals and waters - and for an enzyme tunnel a heme or bound
    # metal frequently defines the cavity wall and the bottleneck, so the tunnel
    # returned would be a different tunnel with no warning. This fixture is
    # SEPARATE because the assertions above are frozen on exact row positions.
    set hid [mol new $::env(ADAPTER_HET_PDB) waitfor all]
    set ::VMDHole::state(selection) "all"
    set hrows [::VMDHole::_tunnel_atoms_mole $hid 0]
    set helem {}
    foreach r $hrows { lappend helem [lindex $r 3] }
    # The FE line in this fixture has a BLANK element column, so the adapter MUST
    # fall back to the atom-name rule. Under the old residue-first-letter rule it
    # became "H" (from HEM) - hydrogen's radius on an iron - which is exactly the
    # silent chemical reinterpretation this guards.
    chk "blank element column resolves from the atom NAME, not the residue" \
        [expr {[lsearch -exact [string toupper $helem] FE] >= 0
               && [lsearch -exact $helem H] < 0}] "(elements: $helem)"
    chk "selection 'all' keeps HETATM cofactors/waters (got [llength $hrows] rows)" \
        [expr {[llength $hrows] == 5 \
               && [lsearch -exact [string toupper $helem] FE] >= 0}] "(elements: $helem)"
    set ::VMDHole::state(selection) "protein"
    set prows [::VMDHole::_tunnel_atoms_mole $hid 0]
    # Asserted as "fewer", not an exact count: what matters is that the two
    # selections differ at all - an exact number would also depend on how VMD
    # classifies a minimal residue, which is not what this check is about.
    chk "selection 'protein' still excludes them" \
        [expr {[llength $prows] < [llength $hrows]}] \
        "(protein=[llength $prows] all=[llength $hrows])"
    set ::VMDHole::state(selection) "all"
    set ncol {}
    foreach r $rows { lappend ncol [llength $r] }
    chk "every row has 12 columns" [expr {[lsort -unique $ncol] eq {12}}] "(got [lsort -unique $ncol])"
    # Column 11 is altLoc, 12 the insertion code - checked by VALUE and by
    # POSITION, so a swap fails.
    set alt {}; set ins {}; set seq {}
    foreach r $rows { lappend alt [lindex $r 10]; lappend ins [lindex $r 11]; lappend seq [lindex $r 6] }
    chk "altLoc column, by position" [expr {$alt eq {- - - - - - - - - A - - A}}] "(got $alt)"
    chk "insertion column, by position" [expr {$ins eq {- - - - A A A A - - - - -}}] "(got $ins)"
    chk "resid column unchanged by the insertion code" \
        [expr {$seq eq {10 10 10 10 10 10 10 10 11 11 11 11 11}}] "(got $seq)"
    # The engine must then accept it and key GLY 10A apart from ALA 10.
    say "  ---- adapter checks, $fails failed"
} err]} { say "  FAIL  the adapter check errored: $err"; say "  ---- adapter checks, 1 failed" }
close $LOG
quit
