#!/usr/bin/env tclsh
# ==============================================================================
#  hole.tcl - HOLE 2 in pure Tcl.  STANDALONE.
# ==============================================================================
#
#  A replication of the HOLE 2 pore-analysis program (Smart, Neduvelil, Wang,
#  Wallace & Sansom 1996; Apache License 2.0) with no compiled dependency.
#  Written standalone so it can be validated on its own against the real
#  binary, but with VMDHole in mind: the entry points, the .sph output and the
#  radius profile all use the shapes vmdhole.tcl already parses, so a working
#  engine can be dropped in as a fallback beside the MOLE one.
#
# ------------------------------------------------------------------------------
#  UPSTREAM NOTICE - retained from the HOLE 2 sources this file was ported from
# ------------------------------------------------------------------------------
#  The Fortran routines listed below carry the following header upstream, and it
#  requires that it travel with any extract. A translation is an extract, so it
#  is reproduced here verbatim:
#
#    ******************************************************************
#    *                                                                *
#    * This software is an unpublished work containing confidential   *
#    * and proprietary information of Birkbeck College. Use,          *
#    * disclosure, reproduction and transfer of this work without the *
#    * express written consent of Birkbeck College are prohibited.    *
#    * This notice must be attached to all copies or extracts of the  *
#    * software.                                                      *
#    *                                                                *
#    * (c) 1995 Oliver Smart & Birkbeck College, All rights reserved  *
#    * (c) 1996 Oliver Smart & Birkbeck College, All rights reserved  *
#    * (c) 1997 Oliver Smart                                          *
#    *                                                                *
#    ******************************************************************
#
#  The HOLE 2 repository (github.com/osmart/hole2) ships a top-level LICENSE
#  granting Apache License 2.0, "Copyright 1993-2023 Oliver Smart", and that
#  grant is the basis on which this port is distributed.
#
#  Oliver S. Smart, asked directly about the contradiction between the two,
#  confirmed in writing: "The Apache 2.0 licence covers everything. You are
#  welcome to use the source code in any way you find useful."
#
#  The header above is therefore reproduced because it asks to be, not because
#  its terms apply - see vmdhole/NOTICE.md.
#
# ------------------------------------------------------------------------------
#
#  Ported from the Apache-licensed Fortran in hole2/src, read directly:
#    holcal.f   the Monte Carlo / simulated-annealing driver
#    holeen.f   the objective function (distance to the nearest atom surface)
#    honewp.f   the perturbation step
#    addend.f   the terminal-sphere grid written to a .sph when -sph is given
#               (hole::addend) - an earlier revision of this file investigated
#               this without porting it; now ported, see hole::write_sph's own
#               header for the byte-identical verification.
#
# ------------------------------------------------------------------------------
#  STATUS - read this before trusting any number out of it
# ------------------------------------------------------------------------------
#  IMPLEMENTED and ported from the real source:
#    * the RNG, exactly (see hole::rng - this is the piece that makes
#      bit-comparable validation possible at all)
#    * PDB input and van der Waals radius assignment (.rad files) - reads the
#      real VDWR card, not a "RADIUS" card that does not exist in any real
#      .rad file (such a gate makes every atom fall back to a hardcoded
#      element radius silently; fixed - see read_rad_file)
#    * HOLEEN, the objective function, including the three-nearest-atom
#      bookkeeping the callers depend on
#    * HONEWP, the in-plane random step
#    * the HOLCAL annealing loop and its linear cooling schedule
#    * CGUESS's CPOINT search (cguess.f) - the CA/all-atom centroid plus its
#      5-cycle grid hill-climb, used when no -cpoint is given
#    * sampling along the axis in both directions with the ENDRAD stop -
#      including that the slice which reaches ENDRAD is NOT stored, matching
#      holcal.f's `IF (-LOWENG.LT.ENDRAD)` storage guard
#    * the .sph writer, BYTE-IDENTICAL to the reference-build hole (see vmdhole/hole_tcl/README.md)'s own sphpdb
#      output - discovery order, real IREC numbering, the ADDEND terminal-
#      sphere grid/LAST-REC-END markers and the record-0 duplicate at the
#      direction transition, not just the stock per-slice records
#      sph_process reads - see hole::write_sph's own header for the exact
#      reproducing command and vmdhole/hole_tcl/tests/sph_addend_test.sh for the
#      durable regression. capsule.tcl/connolly.tcl have their own writers
#      (hole::write_capsule_sph, hole::write_connolly_sph) verified the same
#      way.
#    * the 3D surface: sph_process (sph_process.tcl, dot generation per
#      sph_process.f's -sos/-dotden path) and sos_triangle (sos_triangle.tcl,
#      the advancing-front triangulator in sos_triangle.c - NOT trisphere.f,
#      which that file has no connection to; see sos_triangle.tcl's header).
#      Verified BYTE-IDENTICAL end to end against the real binaries on 1GRM
#      at dotden 10 (out.sph -> .sos -> vmd_plot): 2657/2657 dots, 2192/2192
#      after dedup, 4625/4625 raw triangles, 3280/3280 final triangles, every
#      vertex/normal/color line matching `diff -q`. Getting there required
#      replicating two single-precision truncations most ports would miss -
#      RVEC3 in sphqpu.f is REAL, and ptgen.f's ZANG ratio is evaluated in
#      single precision because PI is not yet in play when it's computed
#      (unlike the X/Y angle just below it, where PI comes first and forces
#      double arithmetic throughout) - see sph_process.tcl. Both stages are
#      O(n) per candidate dot / triangle, same complexity as the C/Fortran,
#      but ~150-200x slower wall-clock in the Tcl interpreter - the same band
#      already measured for the profile search below.
#    * CONNOLLY (concal.f + coarea.f) - see connolly.tcl: the flood fill that
#      grows a set of probe-radius spheres out from each slice's HOLE centre,
#      and the adaptive black/grey/white area estimator that turns them into
#      an equivalent radius (Requiv). Verified against the real binary at the
#      point level (not just the final number): given the SAME full-precision
#      slice centre (concal.f's caller only ever hands it the annealed
#      LOWCEN, never a fresh search), the flood-fill point count, the escape
#      sentinel (Requiv=1e6), the small-pore fallback (Requiv=HOLE radius),
#      and all five area-refinement cycles' Requiv/black-area/grey-area/
#      squares-stored match the reference exactly. NOTE: feeding it the
#      TEXT-OUTPUT-ROUNDED centre (3 decimal places) instead does NOT match
#      past the first cycle - the coarse grid is that sensitive to the
#      centre's own precision, not a porting bug (see vmdhole/hole_tcl/connolly.tcl's
#      header for the DIFF3+DIFF3 typo it deliberately preserves, and for
#      what IS and is not exercised by an axis-aligned CVECT).
#
#  DOES NOT EXIST IN THE REAL SOURCE for the plain spherical calculation this
#  file implements, despite looking like it should - see vmdhole/hole_tcl/refine.tcl:
#    * a post-annealing refinement step. hsbxmi.f exists but is called from
#      exactly one place (holcal.f:493) and is gated `IF (LSPHBX)` - the
#      spherebox option (a SPHPO card), off by default. The plain path's
#      answer for a slice IS the annealing result, unrefined.
#
#  So: this computes a 2D spherical (and, see capsule.tcl, capsule) pore
#  profile. Measured against the real binary on 1GRM (gramicidin) with a
#  pinned seed: bottleneck radius/location and a paired-slice RMS now land
#  inside the real binary's OWN seed-to-seed reproducibility band - see
#  vmdhole/hole_tcl/tests/profile_vs_reference.sh - but this has only been checked
#  on that one structure and CVECT along an axis; do not assume it holds
#  everywhere HOLE does.
#
# ------------------------------------------------------------------------------
#  Usage
# ------------------------------------------------------------------------------
#    tclsh hole.tcl -pdb FILE ?-cpoint "x y z"? ?-cvect "x y z"?
#                   ?-sample 0.25? ?-endrad 22.0? ?-rad FILE? ?-seed N?
#                   ?-sph OUT.sph? ?-csv OUT.csv?
#                   ?-method spherical|capsule|connolly?
#                   ?-probe 1.15? ?-grid N?   (connolly only - CONNR(1)/(2))
#
#  With no -cpoint, CGUESS's own initial search is used (cguess_cpoint),
#  matching what HOLE does when the card is absent - NOT a plain centroid.
# ==============================================================================

namespace eval hole {
    variable VERSION "0.2-connolly"
}

# CONNOLLY (concal.f + coarea.f) lives in its own file - sourced relative to
# this script so `tclsh hole.tcl` works regardless of the caller's cwd.
source [file join [file dirname [info script]] connolly.tcl]

# ==============================================================================
#  1. The random number generator - ported EXACTLY
# ==============================================================================
#
#  This is the piece that decides whether the whole project can be validated.
#
#  holcal.f calls DRAND, and DRAND is NOT in the Apache source release - it
#  comes from a library HOLE links against but does not ship. It is recoverable
#  from the binary: `nm ~/hole2/exe/hole` shows drand_ at 0x22180, and
#  disassembling it shows the body is
#
#      first call:  srand(cseed); if (cseed even) cseed++          [see below]
#      every call:  _gfortran_rand(0)  ->  float in [0,1)
#
#  so DRAND is gfortran's own RAND intrinsic. libgfortran implements that as
#  the Park-Miller "minimal standard" Lehmer generator with Schrage's method:
#
#      seed <- 16807 * seed  (mod 2147483647)
#      value = seed / 2147483647
#
#  reproduced below. That means a Tcl run CAN be made to draw the identical
#  stream as the Fortran one, which is the only way to tell a real porting
#  error from HOLE's own Monte Carlo noise (measured at 0.0053 A between
#  identical reruns of the real binary).
#
#  NOT yet verified against the binary: duran3_ (0x253f0), the random unit
#  vector. The rejection sampler below is the textbook method and consumes 3
#  draws per accepted vector; if a stream comparison ever disagrees, this is
#  the first place to look, and the address above is where to look it up.

namespace eval hole::rng {
    variable seed 1
}

proc hole::rng::srand {s} {
    # Reproduces libgfortran's _gfortran_srand exactly - see hole::rng::rand's
    # header for the second trap this has to clear
    # (objdump --disassemble=_gfortran_srand on the real libgfortran.so.5):
    # a raw store of the given seed, with a single fallback (0 -> a fixed
    # constant) and NO oddification and NO special-casing of any value,
    # including 1. The "oddify" step in the earlier version of this proc is
    # real, but it belongs to HOLE's own drand_ wrapper (machine_dep.f), not
    # to gfortran's SRAND - see hole::rng::seed_like_drand below, which is
    # what a real HOLE run must call instead of this proc directly. This proc
    # stays a literal SRAND reproduction because vmdhole/hole_tcl/tests/rng_stream_test.sh
    # calls it directly and compares against `CALL SRAND(seed)` - changing its
    # semantics would fail that test even if it happened to still work for a
    # full run.
    variable seed
    set s [expr {int($s) & 0xffffffff}]
    if {$s == 0} {
        set seed 0x75bd924
    } else {
        set seed $s
    }
}

proc hole::rng::_schrage_step {s} {
    # One Park-Miller/Lehmer step via Schrage's method (a=16807,
    # m=2147483647=2^31-1): s <- 16807*s mod m. This is the SAME arithmetic
    # _gfortran_irand uses on both of its two call shapes (RAND(0) at
    # objdump offset 2d2fcf, and RAND(nonzero-not-1) at 2d2f80) - factored out
    # so hole::rng::rand and hole::rng::rand_seeded can't drift apart.
    set k [expr {$s / 127773}]
    set r [expr {16807 * ($s - $k * 127773) - 2836 * $k}]
    if {$r < 0} { incr r 2147483647 }
    return $r
}

proc hole::rng::_masked_float {s} {
    # _gfortran_rand (the wrapper _gfortran_irand's
    # result goes through on EVERY call, objdump --disassemble=_gfortran_rand)
    # does NOT return raw_state/2147483647.0. It computes
    # ((raw_state - 1) & 0xFFFFFE00) / 2147483648.0 - i.e. it discards the
    # low 9 bits of (state-1) and divides by 2^31, not 2^31-1. This was
    # invisible for "typical" large seed states (masking off <=511 out of a
    # ~2-billion-scale number is a <=2.4e-7 relative error, under most
    # tolerances) but is a huge relative error for a small state: verified
    # against two independent real gfortran programs -
    #   RAND(1) directly:        state=520932930 (see rand_seeded) -> 0.242578268...
    #   SRAND(1) then RAND(0):   state=16807                       -> 0.0000076293945312 (exactly 2^-17)
    # both match this exact transform and NOTHING ELSE tried (plain
    # state/2147483647.0, state/2147483648.0 unmasked, etc. all fail the
    # second case by orders of magnitude).
    set m [expr {($s - 1) & 0xfffffe00}]
    return [expr {$m / 2147483648.0}]
}

proc hole::rng::rand {} {
    # _gfortran_rand() / Fortran RAND(0): ALWAYS an unconditional "continue"
    # step off the persisted state (objdump: RAND(0)'s argument is a pointer
    # to the value 0, which takes the `test ebx,ebx; je <continue-path>`
    # branch in _gfortran_irand STRAIGHT PAST the "cmp $0x1" special case -
    # that special case is only ever reachable from an EXPLICIT NONZERO
    # argument, i.e. hole::rng::rand_seeded below, never from a plain
    # continuation call. A CURRENT seed that happens to equal 1 here does
    # NOT retrigger anything special - confirmed empirically:
    # vmdhole/hole_tcl/tests/rng_stream_test.sh's SRAND(1)+RAND(0)x50 reference never
    # takes that branch and this proc must not either.
    variable seed
    set seed [_schrage_step $seed]
    return [_masked_float $seed]
}

proc hole::rng::rand_seeded {iseed} {
    # _gfortran_rand() / Fortran RAND(iseed) with an EXPLICIT NONZERO
    # argument - the only call shape that can hit the seed==1 special case
    # (objdump: "cmp $0x1,%ebx; jne <normal Schrage path>; mov
    # $0x1f0cce42,%ebx" at _gfortran_irand+0x2a). Sets the persisted state
    # directly from iseed (a real reseed, unlike hole::rng::rand's
    # continuation), then applies the same masking as every other draw. Not
    # used directly by the search - only by hole::rng::seed_like_drand.
    variable seed
    if {$iseed == 1} {
        set seed 520932930
    } else {
        set seed [_schrage_step $iseed]
    }
    return [_masked_float $seed]
}

proc hole::rng::seed_like_drand {fseed} {
    # Reproduces machine_dep.f's DRAND on its FIRST invocation - the actual
    # seeding path a real HOLE run takes, as opposed to the SRAND-then-RAND(0)
    # shape hole::rng::srand/rand model for the unit test. Traced from
    # `objdump --disassemble=drand_` on the real ~/hole2/exe/hole binary
    # (drand_ at 0x22180, not libgfortran - it is HOLE's own wrapper):
    #   IF (FSEED.NE.0) ISEED=FSEED ELSE CALL GSEED(ISEED)   [time-based;
    #     out of scope here, callers must supply a real seed]
    #   IF (MOD(ISEED,2).EQ.0) ISEED=ISEED+1                 (oddify)
    #   RANDOM = RAND(ISEED)         <- explicit-argument call: THIS is where
    #                                    the seed==1 special case can fire
    #   RANDOM = RAND(0)             <- drand_ makes this SECOND call
    #                                    immediately and unconditionally, and
    #                                    IT OVERWRITES the first call's result
    #                                    before ever returning to the caller
    #                                    (both calls write into the same
    #                                    output slot, back to back, no
    #                                    branch in between)
    # So DRAND's first RETURNED value is a "continue" step applied ONCE on
    # top of the seeded state, and every later DRAND call (holcal.f calls it
    # once per use) is exactly one more "continue" step - i.e. IDENTICAL to
    # hole::rng::rand's own definition. Priming the persisted state to the
    # SEEDED value here (via rand_seeded, discarding ITS return - drand_
    # discards it too) is therefore sufficient: no "is this the first call"
    # flag is needed anywhere else in the port, because the very next
    # hole::rng::rand() call anywhere (honewp, the Metropolis accept check,
    # ...) reproduces DRAND's first returned value bit-for-bit, and every
    # call after that reproduces every later DRAND call the same way.
    set s [expr {int($fseed) & 0x7fffffff}]
    if {$s == 0} { set s 1 }
    if {($s % 2) == 0} { incr s }
    rand_seeded $s
    return
}

proc hole::rng::kick_off {} {
    # FOUND BY DIRECT TRACE COMPARISON, was missing entirely: hole.f:548-549,
    # right after CVECT is unitised/CALPER computed and before HOLCAL is ever
    # called -
    #   "New feature 24 October 1993 write out seed integer used
    #    find one random 3d vector to kick things off"
    #    CALL dURAN3( DISP)
    # - a single dURAN3 call whose result (DISP) is never used for anything;
    # it exists purely to consume 3 DRAND draws once per run, UNCONDITIONALLY
    # (not gated by LCAPS/LSPHBX/SHORTO), before the annealing search's own
    # first honewp call. Every real HOLE run's RNG stream is offset by these
    # 3 draws relative to one that skips this step - which this port did,
    # silently, until caught by a step-by-step trace against an instrumented
    # holcal.f/honewp.f (E24.16 WRITEs, revert after use - see hole_tcl
    # investigation): duran3(DRAND calls 4,5,6) reproduced the real binary's
    # FIRST honewp call's unit vector bit-for-bit, meaning calls 1-3 were
    # consumed by something before honewp - traced to this line. Confirmed
    # against an isolated probe program (FSEED=1, 24 raw DRAND calls linked
    # against the real hole.a): hole::rng::seed_like_drand + hole::rng::rand
    # alone reproduces that probe's stream exactly (the RNG model itself was
    # never wrong) - it is this missing call, not the RNG, that desynced
    # every downstream draw. Callers must invoke this exactly once per run,
    # after seed_like_drand and before the first search call.
    unit_vector
    return
}

proc hole::rng::unit_vector {} {
    # NOT textbook rejection sampling in
    # the unit ball (uniform on the sphere, but a VARIABLE number of draws
    # per call - reject about 47.6% of the time, since a unit ball has
    # volume pi/6 of its bounding [-1,1]^3 cube). The real duran3_ (0x253f0)
    # does no such thing - `objdump --disassemble=duran3_` on the real
    # ~/hole2/exe/hole binary shows a straight-line sequence with NO
    # branches at all: three drand_ calls, each shifted by -0.5 (the exact
    # IEEE754 double 0x3FE0000000000000 read out of .rodata at both
    # 0x3b1f0 and 0x3b2a0), then an unconditional tail-call to duvec2_
    # (0x2dde0, a plain normalise: divide each component by sqrt(sum of
    # squares), no zero-length guard). So this draws a point uniform in the
    # CUBE [-0.5,0.5]^3 and normalises it - NOT uniform on the sphere
    # (biased toward the cube's corners/edges) - but that bias is exactly
    # what the real search has, and reproducing it is the point: this
    # consumes EXACTLY 3 draws, always, never more - unlike rejection
    # sampling, whose variable draw count silently desynchronises the RNG
    # stream from the real binary's from the very first rejected draw
    # onward. Confirmed via honewp_'s own disassembly (0x1bb80): ONE call to
    # duran3_, no retry loop around it either.
    set x [expr {[rand] - 0.5}]
    set y [expr {[rand] - 0.5}]
    set z [expr {[rand] - 0.5}]
    set n [expr {sqrt($x*$x + $y*$y + $z*$z)}]
    return [list [expr {$x/$n}] [expr {$y/$n}] [expr {$z/$n}]]
}

# ==============================================================================
#  2. Input: PDB coordinates and van der Waals radii
# ==============================================================================
#
#  HOLE reads a .rad file mapping (residue, atom) to a radius, with wildcards,
#  and falls back to an element-based guess. The shipped simple.rad /
#  amberuni.rad format is:
#
#      RADIUS <atom-pattern> <residue-pattern> <value>
#
#  patterns may contain ###/*** wildcards. Matching is longest-specific-first,
#  which is what the ordering below reproduces: an exact atom+residue match
#  beats an atom-only match beats the element default.

proc hole::read_rad_file {path} {
    # Does NOT gate on a "RADIUS" keyword, which does
    # not exist in any real HOLE .rad file. tsradr.f reads two record types,
    # 'VDWR' (van der Waals radius, atom+residue wildcard, what HOLEEN needs)
    # and 'BOND' (bond radius, only for the unported MOLQPT sausage plot -
    # skipped here). With the old gate this proc always returned an empty
    # ruleset and every atom silently used element_radius below instead of
    # the file's own values - measured effect on 1GRM: O is 1.65 in
    # simple.rad but 1.40 in element_radius, a 0.25 A error on every oxygen.
    set rules {}
    if {![file exists $path]} { return $rules }
    set fh [open $path r]
    while {[gets $fh line] >= 0} {
        set line [string trim $line]
        if {$line eq "" || [string index $line 0] eq "#"} { continue }
        if {![string match -nocase "VDWR*" $line]} { continue }
        set f [lrange [split [regsub -all {\s+} $line " "] " "] 1 end]
        if {[llength $f] < 3} { continue }
        lassign $f atom res val
        if {![string is double -strict $val]} { continue }
        # tsradr.f reads these as FIXED-WIDTH fields (atom 4 chars, residue 3
        # chars: LINE(6:9), LINE(11:13)) so a trailing blank can be a real
        # part of the wildcard pattern (e.g. "E2? " - GLN's 3-char xplor
        # atom name plus a required blank 4th character). Whitespace-split
        # loses that blank; restore the field width here.
        set atom [string range [format "%-4s" $atom] 0 3]
        set res  [string range [format "%-3s" $res]  0 2]
        lappend rules [list $atom $res $val]
    }
    close $fh
    return $rules
}

# Element fallback, from HOLE's own documentation of the simple.rad defaults.
# Fluorine is present on purpose: its absence from an earlier radius table was
# a real crash ("no pore found" that was actually a missing entry).
#
# NOT radius_for's own unconditional
# fallback - called for ANY atom no VDWR rule matched, even when a real
# radius file WAS supplied. tsatr.f (the real VDWR-rule matcher) has NO
# fallback of its own in that case: if no rule matches, it sets LERR=.TRUE.
# and ABORTS ("Cannot find vdW radius for atom: <name>"), it does not
# substitute anything. Confirmed by running the real binary on a
# K+-containing PDB with simple.rad (K has no VDWR rule in that file): fatal
# error, no profile - silently returning 1.85 instead
# (carbon's radius) for K, producing a full, plausible-looking profile for an
# input the real binary refuses to run at all (see
# vmdhole/hole_tcl/tests/fixtures/1BL8_protein_only.pdb's own REMARK for the worked
# example that found this). Fixed: radius_for below now only reaches this
# proc when NO radius file was supplied at all (rules is empty) - a separate,
# still-unaudited case this session did not investigate (does the real
# binary's own FRADIU='NONE' no-card default abort too, or use some built-in
# table? not checked - every test in this project always passes -rad, so
# this path is not exercised by anything that ships). When a radius file WAS
# supplied and simply does not cover the atom, radius_for now errors instead
# of calling this proc - see radius_for's own comment.
proc hole::element_radius {atname} {
    set a [string trim $atname]
    set e [string toupper [string index $a 0]]
    if {[string is digit -strict $e] && [string length $a] > 1} {
        set e [string toupper [string index $a 1]]
    }
    switch -- $e {
        H { return 1.00 }
        C { return 1.85 }
        N { return 1.75 }
        O { return 1.40 }
        S { return 2.00 }
        P { return 2.10 }
        F { return 1.47 }
        default { return 1.85 }
    }
}

proc hole::radius_for {rules atname resname} {
    # Same fixed-width reasoning as read_rad_file: pad the query name/residue
    # to ATBRK/ATRES's widths (4, 3) so a pattern like "C???" - which needs
    # an exact 4-character match under [string match] - matches "CA" the way
    # LMATCH matches "CA  " in the Fortran. tsradr.f UCASEs the whole rad
    # file but tsatr.f only UCASEs ATRES, not ATBRK; -nocase on both is a
    # harmless superset of that (PDB atom names are already upper case in
    # practice) and the closer of the two to get wrong.
    set an4 [string range [format "%-4s" $atname]   0 3]
    set rn3 [string range [format "%-3s" $resname]  0 2]
    foreach r $rules {
        lassign $r pat_a pat_r val
        if {[string match -nocase $pat_a $an4] && [string match -nocase $pat_r $rn3]} {
            return $val
        }
    }
    if {[llength $rules] > 0} {
        # Does NOT fall through to element_radius
        # unconditionally. tsatr.f has no such fallback when a real VDWR rule
        # list IS in play (rules non-empty here) and none of its patterns
        # match - it sets LERR=.TRUE. and ABORTS with exactly this message
        # (tsatr.f:213-218), it does not substitute a plausible radius. See
        # element_radius's own header for the measured divergence this closes
        # (K+ with simple.rad) and vmdhole/hole_tcl/tests/fixtures/
        # 1BL8_protein_only.pdb's REMARK for the real binary's own abort text.
        error "hole: Cannot find vdW radius for atom: $atname $resname"
    }
    # No radius file was supplied at all - rules is empty for every atom, not
    # just this one, so there is no "no rule matched" case to abort on; keep
    # the permissive element fallback (see its own header for the caveat).
    return [element_radius $atname]
}

proc hole::read_pdb {path {radfile ""}} {
    # Returns {n xs ys zs rs ans} as parallel flat lists - the shape the inner
    # loop wants. Columns are FIXED-WIDTH per the PDB spec, not
    # whitespace-split: a coordinate can run into its neighbour on a large
    # structure. ans (trimmed atom names) is carried along only for CGUESS's
    # "is this a CA" test - callers that only want {n xs ys zs rs} can still
    # lassign the old way, the extra element is just ignored.
    set rules [expr {$radfile ne "" ? [read_rad_file $radfile] : {}}]
    set xs {}; set ys {}; set zs {}; set rs {}; set ans {}
    set fh [open $path r]
    while {[gets $fh line] >= 0} {
        if {![string match "ATOM  *" $line] && ![string match "HETATM*" $line]} { continue }
        set an  [string trim [string range $line 12 15]]
        set rn  [string trim [string range $line 17 19]]
        set x   [string trim [string range $line 30 37]]
        set y   [string trim [string range $line 38 45]]
        set z   [string trim [string range $line 46 53]]
        if {![string is double -strict $x] || ![string is double -strict $y] \
                || ![string is double -strict $z]} { continue }
        lappend xs $x; lappend ys $y; lappend zs $z
        lappend rs [radius_for $rules $an $rn]
        lappend ans $an
    }
    close $fh
    return [list [llength $xs] $xs $ys $zs $rs $ans]
}

# ==============================================================================
#  3. HOLEEN - the objective function                     (holeen.f)
# ==============================================================================
#
#  The energy of a centre is MINUS the distance from it to the nearest atom
#  SURFACE, i.e. -(min over atoms of |c - a| - vdw_a). Minimising the energy
#  therefore maximises the radius of the sphere that fits without touching an
#  atom, which is what HOLE reports. holeen.f is explicit that the sign flip
#  happens only on return ("N.B. ONLY CHANGE SIGN OF ENERGY ON RETURN"), and
#  holcal.f prints "eff radius = -LOWENG" - so the convention here matches.
#
#  The three nearest atoms are tracked because the callers use them: the second
#  and third distances drive the cut-off list rebuild, and the first two name
#  the residues that line the constriction.
#
#  Returns {energy iat1 iat2 iat3 dat2 dat3}, all 0-based atom indices.

proc hole::holeen {cx cy cz n xs ys zs rs} {
    set energy 99999.0
    set iat1 -1; set iat2 -1; set iat3 -1
    set dat2 99999.0; set dat3 99999.0
    for {set i 0} {$i < $n} {incr i} {
        set dx [expr {$cx - [lindex $xs $i]}]
        set dy [expr {$cy - [lindex $ys $i]}]
        set dz [expr {$cz - [lindex $zs $i]}]
        set d2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
        set v  [lindex $rs $i]
        # Compare SQUARED distances against (bound + vdw)^2 before taking any
        # square root - the Fortran does the same, and the sqrt is what this
        # loop would otherwise spend all its time in.
        set b1 [expr {$energy + $v}]
        if {$d2 < $b1*$b1} {
            set d [expr {sqrt($d2) - $v}]
            set iat3 $iat2; set dat3 $dat2
            set iat2 $iat1; set dat2 $energy
            set iat1 $i
            set energy $d
        } else {
            set b2 [expr {$dat2 + $v}]
            if {$d2 < $b2*$b2} {
                set d [expr {sqrt($d2) - $v}]
                set iat3 $iat2; set dat3 $dat2
                set iat2 $i;    set dat2 $d
            } else {
                set b3 [expr {$dat3 + $v}]
                if {$d2 < $b3*$b3} {
                    set iat3 $i
                    set dat3 [expr {sqrt($d2) - $v}]
                }
            }
        }
    }
    # The sign flip, on return only.
    return [list [expr {-$energy}] $iat1 $iat2 $iat3 $dat2 $dat3]
}

# ==============================================================================
#  3a. CGUESS - the initial CPOINT search when no CPOINT is given (cguess.f)
# ==============================================================================
#
#  hole.f only calls CGUESS when CPOINT or CVECT is still at its -55555
#  sentinel; each of the two guesses is then independently gated inside
#  CGUESS by the same per-component check. This ports the CPOINT half only
#  (cguess.f:124-197): find a start point, then hill-climb it. The CVECT
#  guess (cguess.f:200-243, tries X/Y/Z and picks whichever averages the
#  largest radius over an 11-point line) is NOT ported - -cvect always has a
#  concrete default here, so the Fortran's CVECT(1)==-55555 branch can never
#  be the reason a real run took a different path from this one.
#
#  Centroid choice: CA atoms if there are >=3 of them, else all atoms -
#  matches ATBRK(1:2).EQ.'CA' in cguess.f, which for a normal (non-4-char)
#  atom name is exactly the trimmed name read from columns 13-16.
#
#  Grid search: 5 cycles of a 3x3x3 neighbourhood at 1 Angstrom spacing
#  (including the centre itself, offset 0,0,0), strict ">" so a tie keeps the
#  first point found - loop order is X outer, Y middle, Z inner, reproduced
#  exactly since a tie-break otherwise cannot match.

proc hole::cguess_cpoint {n xs ys zs rs ans} {
    set sx 0.0; set sy 0.0; set sz 0.0
    set csx 0.0; set csy 0.0; set csz 0.0; set nca 0
    foreach a $xs b $ys c $zs nm $ans {
        set sx [expr {$sx + $a}]; set sy [expr {$sy + $b}]; set sz [expr {$sz + $c}]
        if {$nm eq "CA"} {
            incr nca
            set csx [expr {$csx + $a}]; set csy [expr {$csy + $b}]; set csz [expr {$csz + $c}]
        }
    }
    if {$nca >= 3} {
        set cx [expr {$csx/$nca}]; set cy [expr {$csy/$nca}]; set cz [expr {$csz/$nca}]
    } else {
        set cx [expr {$sx/$n}]; set cy [expr {$sy/$n}]; set cz [expr {$sz/$n}]
    }
    for {set cycle 0} {$cycle < 5} {incr cycle} {
        set bestrad -1e10
        set bx $cx; set by $cy; set bz $cz
        for {set dx -1} {$dx <= 1} {incr dx} {
            set tx [expr {$cx + double($dx)}]
            for {set dy -1} {$dy <= 1} {incr dy} {
                set ty [expr {$cy + double($dy)}]
                for {set dz -1} {$dz <= 1} {incr dz} {
                    set tz [expr {$cz + double($dz)}]
                    # holeen's first return is -radius (see above); CGUESS
                    # negates it back to a radius before comparing, same as
                    # cguess.f's "TSTENG = -TSTENG".
                    set rad [expr {-[lindex [holeen $tx $ty $tz $n $xs $ys $zs $rs] 0]}]
                    if {$rad > $bestrad} {
                        set bestrad $rad; set bx $tx; set by $ty; set bz $tz
                    }
                }
            }
        }
        set cx $bx; set cy $by; set cz $bz
    }
    return [list $cx $cy $cz]
}

proc hole::cguess_cvect {cx cy cz n xs ys zs rs} {
    # cguess.f:199-243, the half of CGUESS that picks a CHANNEL VECTOR. HOLE
    # guesses CPOINT and CVECT independently, so a run that names one and not
    # the other still gets the other guessed.
    #
    # It only ever returns a CARDINAL axis. For each of x, y, z it walks 11
    # points one Angstrom apart from -5 to +5 through CPOINT, sums HOLEEN's pore
    # radius at each, and keeps the direction with the largest sum - "the way out
    # of the protein is the way the sphere stays biggest". Hence cguess.f's own
    # warning that it "relies on having a reasonably symmetric channel oriented
    # along x, y or z".
    #
    # Ties keep the EARLIER axis: cguess.f tests SUM_PRAD.GT.MAX_PRAD, strictly
    # greater, walking XCOUNT 1->3 from MAX_PRAD = -1E10.
    set best -1e10
    set bestaxis 3
    for {set axis 1} {$axis <= 3} {incr axis} {
        set sum 0.0
        # DO RCOUNT = -5.,5.001,1. - the 5.001 bound is what includes +5.
        for {set i -5} {$i <= 5} {incr i} {
            set t [expr {double($i)}]
            set tx $cx; set ty $cy; set tz $cz
            switch -- $axis {
                1 { set tx [expr {$cx + $t}] }
                2 { set ty [expr {$cy + $t}] }
                3 { set tz [expr {$cz + $t}] }
            }
            # holeen returns -radius; cguess.f negates it back the same way.
            set sum [expr {$sum - [lindex [holeen $tx $ty $tz $n $xs $ys $zs $rs] 0]}]
        }
        if {$sum > $best} { set best $sum; set bestaxis $axis }
    }
    switch -- $bestaxis {
        1 { return [list 1.0 0.0 0.0] }
        2 { return [list 0.0 1.0 0.0] }
        default { return [list 0.0 0.0 1.0] }
    }
}

# ==============================================================================
#  4. HONEWP - the perturbation step                      (honewp.f)
# ==============================================================================
#
#  A random unit vector, projected into the plane normal to the channel axis,
#  re-unitised, then scaled by MCLEN * rand(). Confining the step to that plane
#  is what keeps a sample point on its own slice: the axial coordinate is set
#  by the sampling loop, never by the search.
#
#  The re-unitising and the MCLEN*rand() scaling must NOT collapse into one
#  combined scalar (s = mclen*rand()/n, then dx*s), which
#  is mathematically the same operation as honewp.f/dUVEC2 but NOT the same
#  sequence of roundings. honewp.f calls dUVEC2(DISP) - a PER-COMPONENT
#  division by n (DISP(i)=DISP(i)/LENF, i=1,2,3, three separate divisions) -
#  and only THEN, back in honewp.f, scales the already-normalised vector by
#  MCLEN then WORK (DISP(i)=DISP(i)*MCLEN*WORK, left-to-right per Fortran's
#  same-precedence evaluation). (dx/n)*mclen*work and dx*(mclen*work/n) round
#  differently in the last bit. Confirmed by reading honewp.f AND dUVEC2
#  (ut_vector.f) directly - not a disassembly guess - and this order is worth
#  matching exactly because HONEWP's first call (step 2 of every slice, right
#  after CURENG stops being the guaranteed-accept 1D20 sentinel) is the single
#  most leveraged point in the whole search: a Metropolis annealing walk is
#  chaotic under floating-point perturbation, so a 1-ULP difference here can
#  flip an accept/reject decision and desynchronise the RNG stream for the
#  rest of the slice. Still not sufficient on its own to reach bit-exactness
#  (see vmdhole/hole_tcl/README.md's HOLCAL trace section) but it removes a genuine,
#  source-confirmed mismatch rather than leaving it in.
proc hole::honewp {cx cy cz vx vy vz mclen} {
    lassign [hole::rng::unit_vector] dx dy dz
    set dot [expr {$vx*$dx + $vy*$dy + $vz*$dz}]
    set dx [expr {$dx - $vx*$dot}]
    set dy [expr {$dy - $vy*$dot}]
    set dz [expr {$dz - $vz*$dot}]
    set n [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
    if {$n < 1.0e-12} { return [list $cx $cy $cz] }
    # Per-component division by n first (dUVEC2), THEN *mclen*work
    # left-to-right (honewp.f) - NOT a single combined scalar multiply.
    set dx [expr {$dx/$n}]
    set dy [expr {$dy/$n}]
    set dz [expr {$dz/$n}]
    set work [hole::rng::rand]
    set dx [expr {$dx*$mclen*$work}]
    set dy [expr {$dy*$mclen*$work}]
    set dz [expr {$dz*$mclen*$work}]
    return [list [expr {$cx + $dx}] [expr {$cy + $dy}] [expr {$cz + $dz}]]
}

# ==============================================================================
#  5. The annealing search for ONE slice                  (holcal.f, DO 20)
# ==============================================================================
#
#  Per sample point: start from the previous slice's best centre, then MCSTEP
#  perturb/evaluate/accept steps with a Metropolis criterion and a LINEAR
#  cooling schedule - holcal.f cools by MCKTIN/(0.9*MCSTEP) per step and floors
#  at zero, so the last ~10% of steps are pure downhill.
#
#  Seeding each slice from the previous one's answer is deliberate in the
#  original and matters: it is what makes the profile continuous rather than a
#  set of independent searches that can each find a different local pore.

proc hole::anneal_slice {cx cy cz vx vy vz n xs ys zs rs mcstep mclen mcktin} {
    set curx $cx; set cury $cy; set curz $cz
    set lowx $cx; set lowy $cy; set lowz $cz
    set cureng 1e20
    set loweng 1e20
    set lowiat1 -1; set lowiat2 -1
    set mckt $mcktin
    set cool [expr {$mcktin / (0.9 * double($mcstep))}]
    for {set step 1} {$step <= $mcstep} {incr step} {
        if {$step == 1} {
            set nx $lowx; set ny $lowy; set nz $lowz
        } else {
            lassign [honewp $curx $cury $curz $vx $vy $vz $mclen] nx ny nz
        }
        lassign [holeen $nx $ny $nz $n $xs $ys $zs $rs] neweng iat1 iat2
        if {$neweng < $loweng} {
            set loweng $neweng
            set lowx $nx; set lowy $ny; set lowz $nz
            set lowiat1 $iat1; set lowiat2 $iat2
        }
        set hgood [expr {$neweng - $cureng}]
        if {$hgood < 0.0} {
            set cureng $neweng
            set curx $nx; set cury $ny; set curz $nz
        } elseif {$mckt > 0} {
            # holcal.f's PROB = EXP(-HGOOD/MCKT) has no
            # clamp - it relies on IEEE underflow-to-zero, which Tcl's exp()
            # also does silently (verified: exp(-720.0) => 2.03e-313, no
            # domain error). The old "< -700 ? 0.0 : exp(...)" guard forced an
            # exact 0.0 across roughly (-745,-700) where Fortran's EXP still
            # returns a tiny nonzero double - see capsule.tcl's matching fix
            # for the full account (found via a step-by-step trace of the
            # capsule search against an instrumented holcal.f).
            set arg [expr {-$hgood/$mckt}]
            set prob [expr {exp($arg)}]
            if {[hole::rng::rand] < $prob} {
                set cureng $neweng
                set curx $nx; set cury $ny; set curz $nz
            }
        }
        set mckt [expr {$mckt - $cool}]
        if {$mckt < 0} { set mckt 0.0 }
    }
    # holcal.f does NOT call a
    # steepest-descent refinement (hsbxmi.f) here. It does not, for the
    # plain spherical calculation this file implements. hsbxmi.f is called
    # from exactly one place (holcal.f:493) and it is gated `IF (LSPHBX)` -
    # the spherebox option, which is off unless a SPHPO card is given, and
    # which optimises a DIFFERENT objective (HSBXEN's box area, plus an
    # extra long-axis-angle degree of freedom) that this file has no
    # representation of. Confirmed empirically too: running the real binary
    # on a plain (no SPHPO) input never prints "Applying sd min..." - see
    # vmdhole/hole_tcl/refine.tcl for the full account. So LOWCEN/LOWENG from the
    # annealing loop above IS holcal.f's final answer for this slice; no
    # refinement step is missing here.
    return [list $lowx $lowy $lowz [expr {-$loweng}] $lowiat1 $lowiat2]
}

# ==============================================================================
#  6. The profile: walk the axis in both directions       (holcal.f outer loop)
# ==============================================================================
#
#  From CPOINT, step by SAMPLE along +CVECT until the fitted radius reaches
#  ENDRAD (the pore has opened out into bulk), then the same in -CVECT.
#  The slice that reaches ENDRAD is NOT kept.
#  holcal.f's storage is `IF (-LOWENG.LT.ENDRAD) THEN store ...` (holcal.f
#  ~574) - the slice at or past ENDRAD only prints "This is an end!" and is
#  never written to STRCEN/STRRAD, so it never reaches the printed profile
#  either. Appending it adds a row the reference never has
#  and skews any point-by-point comparison.

proc hole::profile {args} {
    array set o {
        -sample 0.25 -endrad 22.0 -mcstep 1000 -mclen 0.1 -mckt 0.1 -maxsteps 4000
    }
    array set o $args
    set n  $o(-n);  set xs $o(-xs); set ys $o(-ys); set zs $o(-zs); set rs $o(-rs)
    lassign $o(-cpoint) px py pz
    lassign $o(-cvect)  vx vy vz
    set vn [expr {sqrt($vx*$vx + $vy*$vy + $vz*$vz)}]
    if {$vn < 1e-9} { error "hole: CVECT is a zero vector" }
    set vx [expr {$vx/$vn}]; set vy [expr {$vy/$vn}]; set vz [expr {$vz/$vn}]

    set out {}
    set disc {}
    set posc 0
    set negc 0
    set slice0 ""
    foreach dir {1 -1} {
        set cx $px; set cy $py; set cz $pz
        set t 0.0
        for {set k 0} {$k < $o(-maxsteps)} {incr k} {
            # The first slice of the reverse pass repeats t=0, which HOLE also
            # computes twice; it is dropped here rather than emitted twice.
            #
            # Does NOT reposition from CPOINT
            # ($px $py $pz), the raw input point. holcal.f's own reverse-pass
            # reset (holcal.f ~824-827, "back to original point") is
            # `LOWCEN(i) = STRCEN(i,0) + SAMPLE*CVECT(i)` with SAMPLE already
            # negated - STRCEN(,0) is the STORED slice-0 record, i.e. the
            # ANNEALED answer for the first slice, not CPOINT itself. Those
            # two differ (the anneal can walk the in-plane x,y well away from
            # CPOINT even though CPOINT was the seed) - traced via a
            # step-by-step comparison against an instrumented holcal.f: on
            # 1GRM the annealed slice-0 centre was ~0.2 A from CPOINT in x,y,
            # and using CPOINT here fed every negative-direction slice a
            # different seed than the real search, growing to a ~0.03 A
            # radius residual on that whole side while the positive side (not
            # affected by this bug) matched to print precision. $slice0 is
            # set from the very first appended row below (dir 1, k 0).
            if {$dir < 0 && $k == 0} {
                set t [expr {$t - $o(-sample)}]
                if {$slice0 ne ""} {
                    lassign $slice0 s0x s0y s0z
                } else {
                    set s0x $px; set s0y $py; set s0z $pz
                }
                set cx [expr {$s0x + $vx*$t}]
                set cy [expr {$s0y + $vy*$t}]
                set cz [expr {$s0z + $vz*$t}]
                continue
            }
            lassign [anneal_slice $cx $cy $cz $vx $vy $vz $n $xs $ys $zs $rs \
                        $o(-mcstep) $o(-mclen) $o(-mckt)] bx by bz rad a1 a2
            # This slice is the end - do NOT store it (see header comment).
            if {$rad >= $o(-endrad)} { break }
            # ABSOLUTE axial coordinate - the centre projected onto CVECT -
            # because that is what HOLE's own "cenxyz.cvec" column reports,
            # and a profile that cannot be diffed against the reference is
            # not much use for validating this.
            set tabs [expr {$bx*$vx + $by*$vy + $bz*$vz}]
            lappend out [list $tabs $bx $by $bz $rad $a1 $a2]
            # DISCOVERY-ORDER twin of $out, each row tagged with the real
            # STRNOP/-STRNON index (holcal.f:579,602: STRNOP/STRNON start at
            # -1/0 and are incremented BEFORE this row is stored, so the k-th
            # stored +ve row is IREC=k and the k-th stored -ve row is
            # IREC=-(k+1)) - what write_sph needs to reproduce the real
            # binary's -sph record order/numbering, which is NOT sorted-by-t
            # (see write_sph's own header).
            if {$dir == 1} {
                set irec $posc
                incr posc
            } else {
                incr negc
                set irec [expr {-$negc}]
            }
            lappend disc [list $irec $tabs $bx $by $bz $rad $a1 $a2]
            # holcal.f's STRCEN(,0) - the annealed slice-0 record - is only
            # ever this first appended row (dir 1, k 0); see the reverse-pass
            # reset above.
            if {$slice0 eq ""} { set slice0 [list $bx $by $bz] }
            set t [expr {$t + $dir * $o(-sample)}]
            # Next slice starts from THIS slice's answer, displaced along the
            # axis - the continuity the original depends on.
            set cx [expr {$bx + $vx * $dir * $o(-sample)}]
            set cy [expr {$by + $vy * $dir * $o(-sample)}]
            set cz [expr {$bz + $vz * $dir * $o(-sample)}]
        }
    }
    # Two views of the same search: $out sorted by axial coordinate (CSV/
    # min-radius/CONNOLLY consumers, unchanged contract) and $disc in the
    # real binary's own discovery order with its own IREC numbering (the
    # -sph writer, which cannot use the sorted view - see write_sph).
    return [list [lsort -real -index 0 $out] $disc]
}

# ==============================================================================
#  7. Output
# ==============================================================================

#  PORTED: the -sph/ADDEND path and the LAST-REC-END marker (wpdbsp.f,
#  addend.f, holcal.f ~758-806). An earlier revision of this file investigated
#  this without porting it; recorded here (now past tense) so the reasoning
#  survives:
#
#  - wpdbsp.f's OWN "IF (STRRAD(IREC).GT.ENDRAD) WRITE(...) 'LAST-REC-END'"
#    is DEAD CODE at both of its real call sites in holcal.f: the per-slice
#    call (holcal.f:637) only fires inside "IF (-LOWENG.LT.ENDRAD)", and the
#    other (holcal.f:804, IREC=0) re-emits the already-stored, already-below-
#    ENDRAD record-0 sphere when the reverse pass starts ("to have
#    continuous centre line"). STRRAD(IREC) is < ENDRAD at both, always - so
#    write_sph never emits LAST-REC-END for a plain per-slice row; only
#    hole::addend's own records get one.
#  - The marker that actually matters comes from a THIRD, separate source:
#    addend.f, called once per direction (holcal.f:781) right after that
#    direction's search ends. It resets LASCEN to STRCEN(,STRNOP-1) - "the
#    one before as the last is sometimes wonky" (holcal.f:768-778) - DISTINCT
#    from the STRCEN(,0)-based reset hole::profile's own reverse-pass restart
#    uses (that one is for the direction SWITCH; this one is for the
#    direction END) - then evaluates HOLEEN on a 5x5x5 grid of points
#    spanning +-ENDRAD around it, and writes an ATOM record + its own
#    "LAST-REC-END" line for every grid point whose radius exceeds ENDRAD
#    (addend.f:112-149). hole::addend below is that grid search.
#  - Verified against the reference-build hole (see vmdhole/hole_tcl/README.md) with an actual `sphpdb` card (1GRM,
#    cpoint 0 0 0, cvect 0 0 1, sample 0.25, endrad 22.0, raseed 1): the real
#    .sph is 421 lines = 252 ATOM + 169 LAST-REC-END, discovery-ordered
#    (first record at z=0.000, the +ve direction 0..38, its 125-point ADDEND
#    grid, a DUPLICATE of record 0, the -ve direction -1..-43, its own ADDEND
#    grid) - reproduced by this port `diff -q` clean, see write_sph's own
#    header for the exact command and vmdhole/hole_tcl/tests/sph_addend_test.sh for
#    the durable regression test.
#  - capsule's own .sph writer (a QC1/QC2 two-atom-per-slice record,
#    wpdbsp.f:73-83) is now in capsule.tcl (hole::write_capsule_sph) - it
#    reuses this same hole::addend (ADDEND calls plain HOLEEN even in capsule
#    mode; wpdbsp.f/addend.f take no LCAPS-specific path beyond the QC1/QC2
#    record format itself, confirmed by addend.f's own argument list having
#    no capsule-specific inputs at all).
#  - A Connolly .sph writer (concal.f's own SPDBSP write, gated
#    `IF (CONNR(1).LE.0)` at the call site so CONNOLLY runs write there
#    instead of wpdbsp.f) is separate, larger (concal.f:535-556 writes every
#    flood-fill sphere point per slice, not one record) - now in
#    connolly.tcl (hole::write_connolly_sph); see that proc's own header for
#    the byte-identical verification. hole.tcl's CONNOLLY branch calls it,
#    not this write_sph.

# ==============================================================================
#  ADDEND - the terminal-sphere grid at each direction's end       (addend.f)
# ==============================================================================
#
#  A 5(z, along CVECT) x 5(north) x 5(east) = 125-point grid centred on LASCEN
#  (the ADDEND-specific LOWCEN reset - see write_sph's header, NOT the slice
#  this direction actually stored last). Every point's HOLEEN radius is
#  checked; points whose radius is STRICTLY greater than ENDRAD (addend.f:140,
#  ".GT." not ".GE.") are returned for the caller to write out.
#
#  `sample` is passed SIGNED (its sign, not its magnitude, is all this uses -
#  addend.f's own SSIGN is +1 for SAMPLE>0, -1 otherwise, addend.f:105-109) -
#  the caller must pass the ACTUAL per-direction SAMPLE value (positive for
#  the +ve-direction call, negative for the -ve one), not always the
#  magnitude, or every -ve-direction grid point lands on the wrong side.
#
#  Displacements are accumulated as three SEPARATE additions per coordinate
#  (LASCEN, then +Z, then +N, then +E - addend.f:119-133), not one combined
#  sum, to match Fortran's left-to-right evaluation bit for bit; floating
#  addition is not associative, so this only matters at the last-ULP level,
#  but that is exactly the level this project's other traces have found real
#  bugs at (see vmdhole/hole_tcl/README.md's HOLCAL trace section).
#
#  ZCOUNT/NCOUNT/ECOUNT run over small integers (0..4, -2..2, -2..2) and
#  0.5*ENDRAD is a plain double multiply - REAL(ZCOUNT) etc. in the Fortran is
#  an exact int-to-float32-to-double widening for every value these loops
#  actually take, so there is no single-precision trap to replicate here
#  (unlike MCLEN/CAPSULE_PI elsewhere in this project - checked, not assumed).
#
#  Returns a list of {x y z rad} for every accepted point, in the SAME nested
#  order addend.f's own DO loops visit them (Z outermost, then N, then E) -
#  the .sph record order depends on this, not just the point set.

proc hole::addend {lascen sample cvect perpve perpvn endrad n xs ys zs rs} {
    lassign $lascen lcx lcy lcz
    lassign $cvect cvx cvy cvz
    lassign $perpve ex ey ez
    lassign $perpvn nx ny nz
    set ssign [expr {$sample > 0.0 ? 1.0 : -1.0}]
    set half [expr {0.5*$endrad}]
    set out {}
    for {set zcount 0} {$zcount <= 4} {incr zcount} {
        set zdisp [expr {$ssign * double($zcount) * $half}]
        for {set ncount -2} {$ncount <= 2} {incr ncount} {
            set ndisp [expr {double($ncount) * $half}]
            for {set ecount -2} {$ecount <= 2} {incr ecount} {
                set edisp [expr {double($ecount) * $half}]
                set tx $lcx; set ty $lcy; set tz $lcz
                set tx [expr {$tx + $zdisp*$cvx}]
                set ty [expr {$ty + $zdisp*$cvy}]
                set tz [expr {$tz + $zdisp*$cvz}]
                set tx [expr {$tx + $ndisp*$nx}]
                set ty [expr {$ty + $ndisp*$ny}]
                set tz [expr {$tz + $ndisp*$nz}]
                set tx [expr {$tx + $edisp*$ex}]
                set ty [expr {$ty + $edisp*$ey}]
                set tz [expr {$tz + $edisp*$ez}]
                lassign [hole::holeen $tx $ty $tz $n $xs $ys $zs $rs] eng
                set neweng [expr {-$eng}]
                if {$neweng > $endrad} {
                    lappend out [list $tx $ty $tz $neweng]
                }
            }
        }
    }
    return $out
}

# One .sph ATOM line, the format both the plain per-slice record and
# hole::addend's own records share (wpdbsp.f/addend.f both use
# '(A,I4,4X,3F8.3,2F6.2)' after the fixed 'ATOM      1  QSS SPH S' literal -
# byte-verified against the reference-build hole (see vmdhole/hole_tcl/README.md)'s own output, see write_sph).
proc hole::_sph_atom_line {irec x y z rad1 rad2} {
    return [format "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
        1 "QSS" "SPH" "S" $irec $x $y $z $rad1 $rad2]
}

# ==============================================================================
#  write_sph - the STOCK .sph format sph_process reads          (wpdbsp.f/addend.f)
# ==============================================================================
#
#  Takes hole::profile's DISCOVERY-ORDER return (its second list element, NOT
#  the sorted-by-t one CSV/min-radius/CONNOLLY consumers use) - the real
#  binary's own .sph record order is discovery order, not axial-coordinate
#  order, and re-deriving it from a sorted list is not possible in general
#  (see the header note above).
#
#  Byte-identical to the reference-build hole (see vmdhole/hole_tcl/README.md) on the fixture this project's other
#  fidelity claims use (1GRM, cpoint 0 0 0, cvect 0 0 1, sample 0.25, endrad
#  22.0, raseed 1):
#    cd vmdhole/hole_tcl/reference_bin && ./hole < run_sph.inp   # sphpdb ref.sph, same cards (local build - see vmdhole/hole_tcl/README.md)
#    tclsh vmdhole/hole_tcl/hole.tcl -pdb 1GRM.pdb -rad simple.rad -cpoint "0 0 0" \
#        -cvect "0 0 1" -sample 0.25 -endrad 22.0 -seed 1 -sph tcl.sph
#    diff -q ref.sph tcl.sph      # clean
#  See vmdhole/hole_tcl/tests/sph_addend_test.sh for the durable version of this.
#
#  Two degenerate cases are NOT reproduced, matching this project's own
#  "error rather than approximate" precedent (see hole::element_radius):
#  fewer than 2 stored slices in either direction, where ADDEND's own
#  STRCEN(,STRNOP-1)/STRCEN(,-STRNON+1) LASCEN reset would read a slot this
#  port has no value for (index -1, or a slot the OTHER direction owns).
#  Every fixture this project measures against has >=2 stored slices in both
#  directions; not a gap in any reported number, only in inputs nobody has
#  fed this yet.

proc hole::write_sph {discovery path args} {
    array set o {}
    array set o $args
    foreach req {-n -xs -ys -zs -rs -cvect -sample -endrad} {
        if {![info exists o($req)]} { error "hole::write_sph: missing $req" }
    }
    set n $o(-n); set xs $o(-xs); set ys $o(-ys); set zs $o(-zs); set rs $o(-rs)
    set endrad $o(-endrad)
    set sample [expr {abs($o(-sample))}]
    lassign $o(-cvect) cvx cvy cvz
    set cvn [expr {sqrt($cvx*$cvx + $cvy*$cvy + $cvz*$cvz)}]
    if {$cvn < 1e-9} { error "hole::write_sph: CVECT is a zero vector" }
    set cvx [expr {$cvx/$cvn}]; set cvy [expr {$cvy/$cvn}]; set cvz [expr {$cvz/$cvn}]
    set ncvect [list $cvx $cvy $cvz]
    lassign [hole::calper $ncvect] perpve perpvn

    # $discovery is already dir=1 rows followed by dir=-1 rows, each already
    # in the order they were found (hole::profile never reorders within a
    # direction) - splitting on the IREC sign recovers the two per-direction
    # sequences write_sph needs without re-sorting anything.
    set pos {}; set neg {}
    foreach row $discovery {
        if {[lindex $row 0] >= 0} { lappend pos $row } else { lappend neg $row }
    }
    if {[llength $pos] == 0} {
        error "hole::write_sph: no +ve-direction slice was stored - nothing to write"
    }

    set fh [open $path w]

    foreach row $pos {
        lassign $row irec t x y z rad a1 a2
        puts $fh [hole::_sph_atom_line $irec $x $y $z $rad $rad]
    }
    if {[llength $pos] < 2} {
        close $fh
        error "hole::write_sph: only [llength $pos] +ve slice(s) stored - ADDEND's STRCEN(,STRNOP-1) reset (holcal.f:771-773) needs at least 2; not reproduced for this degenerate case (see this proc's header)"
    }
    lassign [lindex $pos end-1] li lt lx ly lz lrad la1 la2
    foreach pt [hole::addend [list $lx $ly $lz] $sample $ncvect $perpve $perpvn \
                    $endrad $n $xs $ys $zs $rs] {
        lassign $pt ax ay az arad
        puts $fh [hole::_sph_atom_line -888 $ax $ay $az $arad 0.00]
        puts $fh "LAST-REC-END"
    }

    # holcal.f:788-806 - "have we looked down from the initial point?": right
    # after the +ve direction's own ADDEND above, HOLE flips SAMPLE negative
    # and RE-WRITES the already-stored record 0 a second time ("to have
    # continuous centre line"), unconditionally, before the -ve direction's
    # search even starts - reproduced literally: $pos's own first row,
    # written again, byte-identical to its first appearance above (STRCEN(,0)
    # is unchanged by the direction flip, so this is not a fresh computation).
    lassign [lindex $pos 0] irec0 t0 x0 y0 z0 rad0 a10 a20
    puts $fh [hole::_sph_atom_line $irec0 $x0 $y0 $z0 $rad0 $rad0]

    foreach row $neg {
        lassign $row irec t x y z rad a1 a2
        puts $fh [hole::_sph_atom_line $irec $x $y $z $rad $rad]
    }
    if {[llength $neg] < 2} {
        close $fh
        error "hole::write_sph: only [llength $neg] -ve slice(s) stored - ADDEND's STRCEN(,-STRNON+1) reset (holcal.f:775-778) needs at least 2; not reproduced for this degenerate case (see this proc's header)"
    }
    lassign [lindex $neg end-1] li lt lx ly lz lrad la1 la2
    foreach pt [hole::addend [list $lx $ly $lz] [expr {-$sample}] $ncvect $perpve $perpvn \
                    $endrad $n $xs $ys $zs $rs] {
        lassign $pt ax ay az arad
        puts $fh [hole::_sph_atom_line -888 $ax $ay $az $arad 0.00]
        puts $fh "LAST-REC-END"
    }

    close $fh
}

proc hole::write_csv {rows path} {
    set fh [open $path w]
    puts $fh "axial_coordinate_angstrom,radius_angstrom,x,y,z"
    foreach row $rows {
        lassign $row t x y z rad
        puts $fh [format "%.4f,%.4f,%.4f,%.4f,%.4f" $t $rad $x $y $z]
    }
    close $fh
}

proc hole::write_connolly_csv {rows path} {
    # rows: {t x y z rad requiv restim}. requiv is CONCAL's REQUIV as-is,
    # including its sentinels (1e6 escaped, 0.0 no accessible circles, or the
    # plain HOLE radius when the small-pore fallback applied - see
    # connolly.tcl's header). restim is hograp.f's RESTIM: requiv itself, or
    # CRATIO*rad when requiv escaped - what the running conductance sum uses.
    set fh [open $path w]
    puts $fh "axial_coordinate_angstrom,radius_angstrom,requiv_angstrom,requiv_estim_angstrom,x,y,z"
    foreach row $rows {
        lassign $row t x y z rad requiv restim
        puts $fh [format "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f" $t $rad $requiv $restim $x $y $z]
    }
    close $fh
}

# ==============================================================================
#  8. Not implemented - these ERROR rather than approximate
# ==============================================================================
#  Silently returning a plausible-looking number for an unported method is the
#  one failure mode that would be worse than not having the file at all.

# CONNOLLY: ported - see connolly.tcl (hole::connolly, hole::concal,
# hole::coarea, hole::calper), sourced near the top of this file.
# REFINE (hsbxmi.f): investigated, not ported - see refine.tcl for why. It is
# spherebox-only in the real source and unreachable from the plain path this
# file implements, so hole::refine_slice below exists only to fail loudly if
# something ever tries to call it, not because a slice is missing a step.
source [file join [file dirname [file normalize [info script]]] refine.tcl]
# CAPSULE: ported - see capsule.tcl (hole::capsule, hole::hcapen,
# hole::anneal_slice_capsule). Sourced from the same directory as this file so
# `tclsh hole.tcl` works regardless of the caller's cwd.
source [file join [file dirname [file normalize [info script]]] capsule.tcl]

# SURFACE: ported - see sph_process.tcl (hole::sph_process) and
# sos_triangle.tcl (hole::sos_triangle). sos_triangle.tcl depends on
# hole::_f32 from sph_process.tcl, so it must be sourced second. Sourced from
# the same directory as this file so `tclsh hole.tcl` works regardless of the
# caller's cwd.
source [file join [file dirname [file normalize [info script]]] sph_process.tcl]
source [file join [file dirname [file normalize [info script]]] sos_triangle.tcl]

# ==============================================================================
#  9. CLI
# ==============================================================================

proc hole::main {argv} {
    # -mclen/-mckt must NOT default to the exact double 0.1.
    # hole.f's own defaults are "MCLEN = 0.1" and "MCKTIN = 0.1" - bare REAL
    # literals (no D0 suffix), which Fortran evaluates in SINGLE precision
    # and then widens to double on assignment, i.e. the real binary's actual
    # default is DOUBLE(FLOAT32(0.1)) = 0.10000000149011612, not exact double
    # 0.1 - a ~1.5e-9 absolute difference that matters because it is fed
    # straight into honewp's per-step displacement scale and the Metropolis
    # PROB denominator, both leveraged points in a chaotic annealing walk.
    # Found the same way as hole::rng::kick_off: a step-by-step trace against
    # an instrumented holcal.f/honewp.f showed HONEWP's step-2 NEWCEN
    # disagreeing by ~1e-8 relative even with the RNG stream and CPOINT
    # exactly aligned; fixing MCLEN/MCKTIN to this value closed it to
    # print-precision noise (see vmdhole/hole_tcl/README.md's HOLCAL trace section).
    # Any card-supplied -mclen/-mckt (or a real MCDISP/MCKT control card) is
    # unaffected - this only corrects the no-card DEFAULT path, which is what
    # "MCLEN 0.1, kT 0.1" in this file's own banner-matching claim means.
    array set o {
        -pdb "" -cpoint "" -cvect "" -sample 0.25 -endrad 22.0
        -rad "" -seed 1 -sph "" -csv "" -mcstep 1000 -mclen 0.10000000149011612
        -mckt 0.10000000149011612
        -method spherical -probe 1.15 -grid {}
    }
    if {[llength $argv] % 2 != 0} { error "hole.tcl: options must come in pairs" }
    array set o $argv
    if {$o(-pdb) eq ""} {
        puts stderr "usage: tclsh hole.tcl -pdb FILE ?-cpoint \"x y z\"? ?-cvect \"x y z\"?"
        puts stderr "                      ?-sample N? ?-endrad N? ?-rad FILE? ?-seed N?"
        puts stderr "                      ?-sph OUT.sph? ?-csv OUT.csv?"
        puts stderr "                      ?-method spherical|capsule|connolly?"
        puts stderr "                      ?-probe 1.15? ?-grid N?  (connolly only, CONNR(1)/CONNR(2))"
        exit 2
    }
    if {$o(-method) ni {spherical capsule connolly}} {
        error "hole.tcl: -method must be spherical, capsule or connolly, got $o(-method)"
    }
    hole::rng::seed_like_drand $o(-seed)
    # hole.f:549 - see hole::rng::kick_off's header. Ordered before CGUESS
    # below to match the real source (hole.f calls this at line 549, CGUESS
    # at line 625), though CGUESS itself draws nothing so the order has no
    # numeric effect - only the search calls after this point are offset.
    hole::rng::kick_off
    lassign [read_pdb $o(-pdb) $o(-rad)] n xs ys zs rs ans
    if {$n == 0} { error "hole.tcl: no atoms read from $o(-pdb)" }

    set cp $o(-cpoint)
    if {$cp eq ""} {
        # NOT a plain centre-of-geometry.
        # Real HOLE does NOT do that when CPOINT is absent - hole.f calls
        # CGUESS (cguess.f), which seeds from a CA (or all-atom) centroid and
        # then hill-climbs it on a 3x3x3/1-Angstrom grid for 5 cycles. See
        # hole::cguess_cpoint above.
        set cp [cguess_cpoint $n $xs $ys $zs $rs $ans]
    }
    # CVECT is guessed SEPARATELY, from whatever CPOINT we ended up with -
    # hole.f:623 enters CGUESS if EITHER card is missing, and cguess.f then
    # tests for each one independently.
    set cv $o(-cvect)
    if {$cv eq ""} {
        set cv [cguess_cvect [lindex $cp 0] [lindex $cp 1] [lindex $cp 2] \
                    $n $xs $ys $zs $rs]
    }

    set t0 [clock milliseconds]
    if {$o(-method) eq "capsule"} {
        lassign [capsule -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
            -cpoint $cp -cvect $cv -sample $o(-sample) -endrad $o(-endrad) \
            -mcstep $o(-mcstep) -mclen $o(-mclen) -mckt $o(-mckt)] rows discovery
        set dt [expr {[clock milliseconds] - $t0}]

        set minr 1e30; set minz 0.0
        foreach row $rows {
            lassign $row t cx cy cz lvx lvy lvz rad a1 a2
            if {$rad < $minr} { set minr $rad; set minz $t }
        }
        puts [format "hole.tcl %s (CAPSULE): %d atoms, %d slices, %d ms" \
            $hole::VERSION $n [llength $rows] $dt]
        puts [format "minimum effective radius %.4f A at axial coordinate %.3f" $minr $minz]
        if {$o(-sph) ne ""} {
            write_capsule_sph $discovery $o(-sph) -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
                -cvect $cv -sample $o(-sample) -endrad $o(-endrad)
            puts "wrote $o(-sph)"
        }
        if {$o(-csv) ne ""} { write_capsule_csv $rows $o(-csv); puts "wrote $o(-csv)" }
        return
    }

    # profile now returns TWO views of the same search - $rows sorted by
    # axial coordinate (unchanged contract; every consumer below still uses
    # this) and $discovery in the real binary's own discovery order/IREC
    # numbering, needed only by write_sph - see hole::profile's own header.
    lassign [profile -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
        -cpoint $cp -cvect $cv -sample $o(-sample) -endrad $o(-endrad) \
        -mcstep $o(-mcstep) -mclen $o(-mclen) -mckt $o(-mckt)] rows discovery
    set dt [expr {[clock milliseconds] - $t0}]

    if {$o(-method) eq "connolly"} {
        # CONN's own centre is holcal.f's LOWCEN for that slice - exactly what
        # $rows already holds (the annealed spherical answer), never a fresh
        # search - see concal.f/holcal.f, and connolly.tcl's header.
        lassign $cv cvx cvy cvz
        set cvn [expr {sqrt($cvx*$cvx + $cvy*$cvy + $cvz*$cvz)}]
        set cvx [expr {$cvx/$cvn}]; set cvy [expr {$cvy/$cvn}]; set cvz [expr {$cvz/$cvn}]
        set ncvect [list $cvx $cvy $cvz]
        lassign [calper $ncvect] perpve perpvn

        set t1 [clock milliseconds]
        set requivs {}
        foreach row $rows {
            lassign $row t x y z rad
            set res [connolly -centre [list $x $y $z] -cvect $ncvect \
                        -perpve $perpve -perpvn $perpvn \
                        -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
                        -endrad $o(-endrad) -probe $o(-probe) -grid $o(-grid) \
                        -sample $o(-sample)]
            lappend requivs [dict get $res requiv]
        }
        set dtconn [expr {[clock milliseconds] - $t1}]

        # CRATIO: mean(requiv/rad) over slices whose Connolly search did NOT
        # escape - hograp.f 233-247, needed BEFORE any RESTIM/CTORES below.
        set cratsum 0.0; set cratno 0
        foreach row $rows requiv $requivs {
            lassign $row t x y z rad
            if {$requiv < 1.0e6} {
                set cratsum [expr {$cratsum + $requiv/$rad}]
                incr cratno
            }
        }
        set cratio [expr {$cratno > 0 ? $cratsum/double($cratno) : 0.0}]

        set pi [expr {2.0*asin(1.0)}]
        set toresi 0.0; set ctores 0.0
        set minr 1e30; set minz 0.0
        set minequiv 1e30
        set connrows {}
        set sample [expr {abs($o(-sample))}]
        set prevrow {}
        foreach row $rows requiv $requivs {
            lassign $row t x y z rad
            if {$prevrow ne ""} {
                lassign $prevrow pt px py pz prad
                set mx [expr {0.5*($px+$x)}]; set my [expr {0.5*($py+$y)}]; set mz [expr {0.5*($pz+$z)}]
                lassign [holeen $mx $my $mz $n $xs $ys $zs $rs] meng
                set mrad [expr {-$meng}]
                set toresi [expr {$toresi + 0.5*$sample/($pi*$mrad*$mrad)}]
            }
            set toresi [expr {$toresi + 0.5*$sample/($pi*$rad*$rad)}]
            # RESTIM: hograp.f 336-343 - substitute CRATIO*rad for an escaped
            # (>=1e6) slice's Requiv so the running sum stays finite.
            if {$requiv > 0.999e6} {
                set restim [expr {$cratio*$rad}]
            } else {
                set restim $requiv
            }
            if {$restim > 0.0} {
                set ctores [expr {$ctores + $sample/($pi*$restim*$restim)}]
            }
            if {$rad < $minr} { set minr $rad; set minz $t }
            if {$requiv < $minequiv} { set minequiv $requiv }
            lappend connrows [list $t $x $y $z $rad $requiv $restim]
            set prevrow $row
        }
        set gmacro [expr {$toresi > 0 ? 1200.0/$toresi : 0.0}]
        set conngmacro [expr {$ctores > 0 ? 1200.0/$ctores : 0.0}]

        puts [format "hole.tcl %s (CONNOLLY): %d atoms, %d slices, %d ms search + %d ms connolly" \
            $hole::VERSION $n [llength $rows] $dt $dtconn]
        puts [format "minimum radius %.4f A at axial coordinate %.3f" $minr $minz]
        puts [format "minimum equivR found: %.4f A" $minequiv]
        puts [format "Gmacro (spherical)= %.4g pS/M   Conn_Gmacro= %.4g pS/M   (1M KCl, rho=1/12 ohm.m)" \
            $gmacro $conngmacro]
        if {$o(-sph) ne ""} {
            write_connolly_sph $discovery $o(-sph) -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
                -cvect $cv -sample $o(-sample) -endrad $o(-endrad) \
                -probe $o(-probe) -grid $o(-grid)
            puts "wrote $o(-sph)"
        }
        if {$o(-csv) ne ""} { write_connolly_csv $connrows $o(-csv); puts "wrote $o(-csv)" }
        return
    }

    set minr 1e30; set minz 0.0
    foreach row $rows {
        lassign $row t x y z rad
        if {$rad < $minr} { set minr $rad; set minz $t }
    }
    puts [format "hole.tcl %s: %d atoms, %d slices, %d ms" \
        $hole::VERSION $n [llength $rows] $dt]
    puts [format "minimum radius %.4f A at axial coordinate %.3f" $minr $minz]
    if {$o(-sph) ne ""} {
        write_sph $discovery $o(-sph) -n $n -xs $xs -ys $ys -zs $zs -rs $rs \
            -cvect $cv -sample $o(-sample) -endrad $o(-endrad)
        puts "wrote $o(-sph)"
    }
    if {$o(-csv) ne ""} { write_csv $rows $o(-csv); puts "wrote $o(-csv)" }
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
    hole::main $argv
}
