#!/usr/bin/env tclsh
# ==============================================================================
#  capsule.tcl - the CAPSULE pore method (hcapen.f, hcapgr.f, holcal.f LCAPS)
# ==============================================================================
#
#  Sourced by hole.tcl, which defines hole::holeen, hole::honewp, hole::rng::*
#  and the plain (spherical) hole::profile - read that file first.
#
#  CAPSULE is not a post-process on top of the spherical search: it replaces
#  the per-slice objective with a STADIUM (two half-spheres joined by a
#  cylinder), described by a PAIR of points (CENTRE, SECCEN a.k.a. LVC in the
#  Fortran). Both points are annealed independently each Monte Carlo step, and
#  the transition between slices carries a real decision (holcal.f:690-747):
#  after moving the found capsule forward by SAMPLE, HOLE tries collapsing it
#  to a single degenerate point at its midpoint and keeps whichever seed - the
#  moved capsule or the collapsed point - gives the better (more negative)
#  energy, biased 2x towards keeping the moved capsule. This is search
#  behaviour, ported here exactly, not a reporting nicety.
#
# ------------------------------------------------------------------------------
#  STATUS
# ------------------------------------------------------------------------------
#  IMPLEMENTED, ported from the real source:
#    * HCAPEN         the capsule objective function (hcapen.f)
#    * the two-point Monte Carlo / simulated-annealing search (holcal.f LCAPS
#      branches), including the moved-vs-midpoint-collapse slice-to-slice seed
#    * the effective-radius profile (the "(sampled)" rows of hcapgr.f's Data
#      table: cenxyz.cvect, eff.rad, cap.rad)
#
#  NOT IMPLEMENTED - these are simply absent from the output, not approximated:
#    * hcapgr.f's "(mid-point)" interpolated rows (extra HCAPEN calls between
#      stored slices, with a choose-the-longer-of-two-vectors mixing rule -
#      a reporting enhancement, not part of the search)
#    * ADDEND-style terminal sphere at the ENDRAD boundary - investigated,
#      not ported; see hole.tcl's write_sph header for the full derivation
#      (it is coupled to the LAST-REC-END marker and the -sph output path,
#      not separable into a quick addition)
#    * the conductance/geometric-factor block at the end of hcapgr.f
#    * HOCAPR/HOCAPD (hocapr.f, hocapd.f) - confirmed by grep to be called
#      only from the dot-surface graphics routines (hodotc.f, sphqpc.f), not
#      from the numeric search/report path, so out of scope here
#
#  A record is stored ONLY if its effective radius is strictly below ENDRAD
#  (holcal.f: "IF (-LOWENG.LT.ENDRAD)") - matching the real output table
#  exactly. This differs from hole::profile in hole.tcl, which appends one
#  extra endrad-exceeding row; that is a pre-existing property of that file,
#  not replicated here.
#
#  VALIDATED against a rebuilt-from-source real binary (1GRM, seed-pinned,
#  raseed 1 and raseed 7 independently) at TWO precisions:
#    - printed (3-decimal) table, via the CLI: eff.rad and cap.rad RMS
#      0.0002-0.0003 A, max 0.0005 A over every paired slice, slice counts
#      matching exactly.
#    - full E24.16 precision, via an instrumented-Fortran trace of the
#      per-slice STORED state (not just the printed table): max|d(eff.rad)|
#      4.4e-16 to 8.9e-16 A, max|d(centre/LVC/dcent)| <= 3.6e-15 A, exact
#      index match (31=31 @ seed 1, 36=36 @ seed 7) - i.e. bit-exact at the
#      machine-precision floor once the CAPSULE_PI fix below is applied. See
#      that fix's own header comment for the full measurement and the
#      reproducing recipe.
#  This search's own logic (the code in THIS file) was never the problem - a
#  step-by-step trace confirmed it matches holcal.f's LCAPS branches exactly,
#  including LOWBRD only ever coming from the post-loop HCAPEN re-evaluation
#  and the moved-vs-midpoint-collapse comparison. The bugs were all in the
#  shared driver capsule.tcl calls into (hole.tcl) or in HCAPEN's own PI
#  constant: a missing RNG warm-up draw (hole::rng::kick_off, hole.f:549),
#  MCLEN/MCKT's default not matching the real binary's actual
#  (single-precision-literal-contaminated) 0.10000000149011612 - not exact
#  double 0.1, and CAPSULE_PI's own single-precision contamination (see
#  below). This comment stays here since a previous revision of this
#  project's own README wrongly speculated the residual was in the capsule
#  search's decision logic itself - the full-precision trace now rules that
#  out directly rather than by elimination.

# ==============================================================================
#  HCAPEN - the capsule objective function                (hcapen.f)
# ==============================================================================
#
#  Returns {energy caprad iat1 iat2 dat2}.
#    energy: -effective_radius, i.e. sqrt(area/pi) with area = pi*r^2 + 2*r*d
#            (d = distance between the two capsule centres, r = raw capsule
#            radius), sign-flipped on return exactly like hole::holeen. If the
#            raw radius r is <= 0 the area transform is skipped (hcapen.f's
#            own "IF (ENERGY.GT.0.)" guard) and energy is just -r.
#    caprad: the raw capsule radius r (positive when open), captured BEFORE
#            the area transform (hcapen.f: "CAPRAD = ENERGY" happens first).
#
#  dat2 (second-nearest clearance, area-transformed with NO positivity guard
#  in the original) is computed for fidelity but is not used by the search or
#  the profile - holcal.f never stores it for the capsule branch.

proc hole::hcapen {cx cy cz sx sy sz n xs ys zs rs} {
    variable CAPSULE_PI
    set ujx [expr {$sx - $cx}]
    set ujy [expr {$sy - $cy}]
    set ujz [expr {$sz - $cz}]
    set dcent [expr {sqrt($ujx*$ujx + $ujy*$ujy + $ujz*$ujz)}]

    # hcapen.f:101 "IF (DCENT.LT.1E-09)" - 1E-09 is a bare REAL literal, so the
    # comparison is really against DOUBLE(FLOAT32(1e-9)) = 9.99999971718068537e-10,
    # not exact double 1e-9 (same class of trap as PI below). Verified against
    # an E24.16 dump of a Fortran `D = 1E-09` assignment: 0.9999999717180685E-09,
    # bit-identical to Tcl's `binary format f 1.0e-9` widened back to double.
    # Provably inert either way for this port's own callers - DCENT is either
    # exactly 0.0 (a post-collapse degenerate slice start) or already several
    # orders of magnitude above 1e-9 from a real MC perturbation - but matched
    # exactly anyway since it costs nothing.
    if {$dcent < 9.99999971718068537e-10} {
        # centre == seccen: hcapen.f falls back to plain HOLEEN (hardcoded
        # cutsize of 5 in the original - hole::holeen has no cutoff list so
        # there is nothing to pass it). holeen already returns -radius, the
        # same convention HCAPEN itself uses on return, so no extra sign flip.
        lassign [hole::holeen $cx $cy $cz $n $xs $ys $zs $rs] energy iat1 iat2 iat3 dat2 dat3
        set caprad [expr {-$energy}]
        return [list $energy $caprad $iat1 $iat2 $dat2]
    }

    set ujx [expr {$ujx/$dcent}]
    set ujy [expr {$ujy/$dcent}]
    set ujz [expr {$ujz/$dcent}]

    set energy 99999.0
    set iat1 -1; set iat2 -1
    set dat2 99999.0
    for {set i 0} {$i < $n} {incr i} {
        set ax [lindex $xs $i]; set ay [lindex $ys $i]; set az [lindex $zs $i]
        set rvx [expr {$ax - $cx}]
        set rvy [expr {$ay - $cy}]
        set rvz [expr {$az - $cz}]
        set rdotu [expr {$rvx*$ujx + $rvy*$ujy + $rvz*$ujz}]
        if {$rdotu < 0} {
            # closest to CENTRE
            set dist [expr {$rvx*$rvx + $rvy*$rvy + $rvz*$rvz}]
        } elseif {$rdotu > $dcent} {
            # closest to SECCEN
            set ddx [expr {$sx - $ax}]; set ddy [expr {$sy - $ay}]; set ddz [expr {$sz - $az}]
            set dist [expr {$ddx*$ddx + $ddy*$ddy + $ddz*$ddz}]
        } else {
            # closest to a point on the centre line
            set rvx [expr {$rvx - $ujx*$rdotu}]
            set rvy [expr {$rvy - $ujy*$rdotu}]
            set rvz [expr {$rvz - $ujz*$rdotu}]
            set dist [expr {$rvx*$rvx + $rvy*$rvy + $rvz*$rvz}]
        }
        set v [lindex $rs $i]
        set b1 [expr {$energy + $v}]
        if {$dist < $b1*$b1} {
            set d [expr {sqrt($dist) - $v}]
            set iat2 $iat1; set dat2 $energy
            set iat1 $i
            set energy $d
        } else {
            set b2 [expr {$dat2 + $v}]
            if {$dist < $b2*$b2} {
                set iat2 $i
                set dat2 [expr {sqrt($dist) - $v}]
            }
        }
    }

    # raw capsule radius, captured before the area transform below
    set caprad $energy
    if {$energy > 0.0} {
        # hcapen.f:197 "ENERGY = PI*ENERGY**2 + 2*ENERGY*DCENT" - Fortran's **
        # binds tighter than *, so this associates as PI*(ENERGY*ENERGY), NOT
        # (PI*ENERGY)*ENERGY. Both are two roundings of the same three inputs
        # but a DIFFERENT pair of roundings, so they can differ in the last
        # bit - the same class of association-order bug already fixed in
        # HONEWP (see that proc's header). Grouped explicitly here to match.
        set area [expr {$CAPSULE_PI*($energy*$energy) + 2.0*$energy*$dcent}]
        set energy [expr {sqrt($area/$CAPSULE_PI)}]
    }
    set energy [expr {-$energy}]

    # dat2's transform has no positivity guard in the original; guard it here
    # only to avoid a Tcl domain error on a negative sqrt argument (dat2 is
    # unused downstream, so a NaN-equivalent sentinel is harmless). Same
    # association-order fix as ENERGY's transform above.
    if {[catch {
        set area2 [expr {$CAPSULE_PI*($dat2*$dat2) + 2.0*$dat2*$dcent}]
        set dat2 [expr {sqrt($area2/$CAPSULE_PI)}]
    }]} {
        set dat2 99999.0
    }

    return [list $energy $caprad $iat1 $iat2 $dat2]
}

# hole.f:363 "PI = 2.*ASIN(1.)" - the literals 2. and 1. are default REAL
# (single precision, no D0 suffix), so ASIN's argument/result and the final
# multiply are all done in float32, THEN widened to DOUBLE PRECISION on
# assignment to PI - the same silent-single-precision-contamination trap as
# MCLEN/MCKTIN's no-card default (see hole.tcl's hole::main) and ptgen.f's
# ring-angle ratio (see sph_process.tcl's hole::_f32/hole::_ptgen). This PI is
# passed into HCAPEN as an argument (hole.f -> holcal.f -> hcapen.f), so every
# capsule area/effective-radius calculation uses DOUBLE(FLOAT32(2*asin(1))) =
# 3.1415927410125732, not true double pi (3.1415926535897931) - a ~2.8e-8
# relative error. Confirmed BIT-IDENTICAL to the real binary by adding a
# TRANSFER-to-INTEGER*8 WRITE right after hole.f:363 in a scratch-copied
# source tree and rebuilding: the real binary itself prints hex
# 400921FB60000000 for PI at runtime (not a proxy - the actual instrumented
# binary), matching Tcl's `binary scan [binary format f [expr
# {2.0*asin(1.0)}]] f y` bit for bit.
#
# Proven to be the ENTIRE source of the ~1e-9-relative
# HCAPEN residual by a FULL per-slice E24.16 trace of the whole two-centre
# search (not just an isolated HCAPEN call): holcal.f's storage site
# (~line 587-598, both STRNOP/STRNON branches) was instrumented to dump
# LOWCEN, LOWLVC, -LOWENG (stored eff.rad) and LOWBRD (cap.rad) per stored
# slice, diffed against a Tcl harness that replays hole::anneal_slice_capsule
# in discovery order (hole_tcl's own trace method - see README's "HOLCAL
# trace" section). 1GRM, cpoint 0 0 4, cvect 0 0 1, sample 0.5, endrad 8.0:
#
#   raseed 1 (31 = 31 slices, exact index match -18..12):
#     WITHOUT this fix (CAPSULE_PI = true double pi): max|d(eff.rad)| =
#       9.405e-08 A (relative 1.4e-8) - centres/LVC/dcent already at 1e-15
#       (i.e. the two-centre SEARCH itself was already bit-exact even
#       without this fix - only the PI-dependent area transform was off).
#     WITH this fix alone (association-order and dcent-threshold left at
#       their pre-fix form): max|d(eff.rad)| = 8.882e-16 A - closed to the
#       machine-precision floor by PI ALONE.
#     WITH all three edits (current form): identical to PI-alone, 8.882e-16.
#   raseed 7 (36 = 36 slices, exact index match): max|d(eff.rad)| = 4.4e-16,
#     max|d(centres/LVC/dcent)| <= 3.6e-15, same pattern.
#
# So the association-order and dcent-threshold edits are CORRECT per source
# (hcapen.f:101,197 - see their own comments above) but MEASURABLY INERT on
# both traced seeds - not "possibly load-bearing", provably so: dcent's
# threshold can only matter for a real dcent within 3e-17 of 1e-9 (the two
# candidate thresholds differ by ~2.8e-17), and every dcent==0.0 case in this
# port's own search (the collapse-to-midpoint reseed in hole::capsule, which
# calls hcapen with centre exactly equal to seccen) is already far below
# either threshold identically, so the choice of threshold cannot change
# which branch is taken. Association order matched bit-for-bit with and
# without the fix at every one of the 67 combined stored slices across both
# seeds; kept for source fidelity, not because it changed any measured value
# here. This also settles Task 2 (the capsule "two-centre separation"
# question): the search state (both centres, not just the radius) is
# bit-exact at the current revision. The "diverges in the two-centre
# separation" language in this file's own history (commits 73cb4ab7,
# 63dcea33) was measured BEFORE the RNG kick-off / MCLEN fixes (02d26f00)
# landed, on a desynchronized RNG stream - stale, not a live contradiction.
#
# Reproducing this measurement: copy native/stock_build/hole2 to a
# scratch dir (never the shared tree), add the TRANSFER/E24.16 WRITEs
# described above to hole.f (after "PI = 2.*ASIN(1.)") and holcal.f (after
# each STRNOP/STRNON storage block), `source ../source.apache && make
# ../exe/hole`, run with `capsul` + the cards above + `raseed N`, and diff
# against a Tcl harness that calls hole::anneal_slice_capsule/hole::hcapen
# directly in discovery order (dir==1: index=k; dir==-1: index=-(k+1)) -
# not through hole::capsule's own CSV output, which is 4-decimal-rounded.
namespace eval hole {
    binary scan [binary format f [expr {2.0*asin(1.0)}]] f _capsule_pi_f32
    variable CAPSULE_PI $_capsule_pi_f32
    unset _capsule_pi_f32
}

# ==============================================================================
#  The annealing search for ONE capsule slice          (holcal.f, DO 20, LCAPS)
# ==============================================================================
#
#  Same Metropolis/cooling structure as hole::anneal_slice, but with a SECOND
#  point (LVC) perturbed by its own HONEWP call each step (holcal.f:293,295 -
#  centre first, then LVC; this order matters for the RNG stream). Only the
#  position of the low-energy LVC is tracked live (holcal.f:358-363); the
#  capsule radius (LOWBRD) is never set inside the loop - it is only produced
#  by the explicit HCAPEN re-evaluation after the loop (holcal.f:509-511),
#  reproduced here as the final call below.
#
#  Returns {lowx lowy lowz lowlvx lowlvy lowlvz effrad caprad iat1 iat2}.

proc hole::anneal_slice_capsule {cx cy cz lvx lvy lvz vx vy vz n xs ys zs rs mcstep mclen mcktin} {
    set curx $cx; set cury $cy; set curz $cz
    set curlvx $lvx; set curlvy $lvy; set curlvz $lvz
    set lowx $cx; set lowy $cy; set lowz $cz
    set lowlvx $lvx; set lowlvy $lvy; set lowlvz $lvz
    set cureng 1e20
    set loweng 1e20
    set mckt $mcktin
    set cool [expr {$mcktin / (0.9 * double($mcstep))}]

    for {set step 1} {$step <= $mcstep} {incr step} {
        if {$step == 1} {
            set nx $lowx; set ny $lowy; set nz $lowz
            set nlvx $lowlvx; set nlvy $lowlvy; set nlvz $lowlvz
        } else {
            lassign [hole::honewp $curx $cury $curz $vx $vy $vz $mclen] nx ny nz
            lassign [hole::honewp $curlvx $curlvy $curlvz $vx $vy $vz $mclen] nlvx nlvy nlvz
        }
        lassign [hole::hcapen $nx $ny $nz $nlvx $nlvy $nlvz $n $xs $ys $zs $rs] \
            neweng newcaprad iat1 iat2 dat2

        if {$neweng < $loweng} {
            set loweng $neweng
            set lowx $nx; set lowy $ny; set lowz $nz
            set lowlvx $nlvx; set lowlvy $nlvy; set lowlvz $nlvz
        }
        set hgood [expr {$neweng - $cureng}]
        if {$hgood < 0.0} {
            set cureng $neweng
            set curx $nx; set cury $ny; set curz $nz
            set curlvx $nlvx; set curlvy $nlvy; set curlvz $nlvz
        } elseif {$mckt > 0} {
            # holcal.f's own PROB = EXP(-HGOOD/MCKT) has
            # no clamp - it relies on IEEE underflow-to-zero for a very
            # negative argument, which Tcl's exp() also does silently (no
            # domain error - verified: exp(-720.0) => 2.03e-313, not an
            # error). The old "< -700 ? 0.0 : exp(...)" guard forced an exact
            # 0.0 across roughly (-745,-700), where Fortran's EXP still
            # returns a tiny nonzero double - a real, if practically
            # unreachable (RANDYN's granularity is ~4.66e-10, far above any
            # of those magnitudes), unfaithfulness to the source.
            set arg [expr {-$hgood/$mckt}]
            set prob [expr {exp($arg)}]
            if {[hole::rng::rand] < $prob} {
                set cureng $neweng
                set curx $nx; set cury $ny; set curz $nz
                set curlvx $nlvx; set curlvy $nlvy; set curlvz $nlvz
            }
        }
        set mckt [expr {$mckt - $cool}]
        if {$mckt < 0} { set mckt 0.0 }
    }

    # holcal.f:509-511 - recompute at (LOWCEN,LOWLVC) purely to obtain LOWBRD
    # (the capsule radius); the reported effective radius still comes from
    # the loop-tracked loweng (holcal.f:590 stores -LOWENG, not this call's
    # energy), kept separate here for the same provenance.
    lassign [hole::hcapen $lowx $lowy $lowz $lowlvx $lowlvy $lowlvz $n $xs $ys $zs $rs] \
        recheck_energy caprad riat1 riat2 rdat2

    return [list $lowx $lowy $lowz $lowlvx $lowlvy $lowlvz [expr {-$loweng}] $caprad $riat1 $riat2]
}

# ==============================================================================
#  hole::capsule - walk the axis in both directions    (holcal.f outer loop)
# ==============================================================================
#
#  Mirrors hole::profile's overall shape (walk +CVECT then -CVECT from
#  CPOINT, stop when the radius reaches ENDRAD) but with capsule state (two
#  points) and the moved-vs-midpoint-collapse transition between slices
#  (holcal.f:690-747).
#
#  The reverse pass anchors at the ACTUAL FOUND record-0 centre/LVC
#  (holcal.f:792-799 use STRCEN(*,0)/STRLVC(*,0)), not the raw input CPOINT.
#  This is a deliberate difference from hole::profile, which restarts the
#  reverse pass from the raw CPOINT as a simplification (see that file's
#  STATUS block); here fidelity to the Fortran was chosen instead, since a
#  seed error at the direction switch would propagate through the whole
#  negative half of the profile.
#
#  Each output row is {t cx cy cz lvx lvy lvz effrad caprad iat1 iat2}, where
#  t is the axial coordinate of the capsule mid-point (0.5*(CENTRE+SECCEN))
#  projected onto CVECT - matching hcapgr.f's CAPCOR/"cenxyz.cvect" column.

proc hole::capsule {args} {
    array set o {
        -sample 0.25 -endrad 22.0 -mcstep 1000 -mclen 0.1 -mckt 0.1 -maxsteps 4000
    }
    array set o $args
    set n $o(-n); set xs $o(-xs); set ys $o(-ys); set zs $o(-zs); set rs $o(-rs)
    lassign $o(-cpoint) px py pz
    lassign $o(-cvect)  vx vy vz
    set vn [expr {sqrt($vx*$vx + $vy*$vy + $vz*$vz)}]
    if {$vn < 1e-9} { error "hole::capsule: CVECT is a zero vector" }
    set vx [expr {$vx/$vn}]; set vy [expr {$vy/$vn}]; set vz [expr {$vz/$vn}]
    set sample [expr {abs($o(-sample))}]

    set out {}
    set disc {}
    set posc 0
    set negc 0
    # holcal.f:246-254 - capsule starts as a DEGENERATE point (both ends at
    # CPOINT).
    set lowx $px; set lowy $py; set lowz $pz
    set lowlvx $px; set lowlvy $py; set lowlvz $pz

    # anchor for the reverse pass, set once the first (+ve, t=0) record is
    # found
    set have_rec0 0
    set rec0x 0; set rec0y 0; set rec0z 0
    set rec0lvx 0; set rec0lvy 0; set rec0lvz 0

    foreach dir {1 -1} {
        set dsample [expr {$dir * $sample}]
        if {$dir == -1} {
            if {!$have_rec0} {
                error "hole::capsule: no slice found in the +ve direction before ENDRAD - cannot anchor the reverse pass"
            }
            # holcal.f:792-799
            set lowx   [expr {$rec0x   + $vx*$dsample}]
            set lowy   [expr {$rec0y   + $vy*$dsample}]
            set lowz   [expr {$rec0z   + $vz*$dsample}]
            set lowlvx [expr {$rec0lvx + $vx*$dsample}]
            set lowlvy [expr {$rec0lvy + $vy*$dsample}]
            set lowlvz [expr {$rec0lvz + $vz*$dsample}]
        }

        for {set k 0} {$k < $o(-maxsteps)} {incr k} {
            lassign [hole::anneal_slice_capsule $lowx $lowy $lowz $lowlvx $lowlvy $lowlvz \
                        $vx $vy $vz $n $xs $ys $zs $rs \
                        $o(-mcstep) $o(-mclen) $o(-mckt)] \
                fx fy fz flvx flvy flvz effrad caprad a1 a2

            # holcal.f: "IF (-LOWENG.LT.ENDRAD)" - only store below ENDRAD
            if {$effrad >= $o(-endrad)} { break }

            set cenx [expr {0.5*($fx+$flvx)}]
            set ceny [expr {0.5*($fy+$flvy)}]
            set cenz [expr {0.5*($fz+$flvz)}]
            set tabs [expr {$cenx*$vx + $ceny*$vy + $cenz*$vz}]
            lappend out [list $tabs $fx $fy $fz $flvx $flvy $flvz $effrad $caprad $a1 $a2]
            # DISCOVERY-ORDER twin, IREC-tagged - same STRNOP/-STRNON
            # numbering as hole::profile's own $disc (holcal.f:579,602); see
            # hole::write_capsule_sph, which needs this and cannot use the
            # sorted $out above (same reasoning as hole::write_sph).
            if {$dir == 1} {
                set irec $posc
                incr posc
            } else {
                incr negc
                set irec [expr {-$negc}]
            }
            lappend disc [list $irec $tabs $fx $fy $fz $flvx $flvy $flvz $effrad $caprad $a1 $a2]

            if {$dir == 1 && !$have_rec0} {
                set rec0x $fx; set rec0y $fy; set rec0z $fz
                set rec0lvx $flvx; set rec0lvy $flvy; set rec0lvz $flvz
                set have_rec0 1
            }

            # ---- seed the NEXT slice (holcal.f:690-747) ----
            set movx   [expr {$fx   + $vx*$dsample}]
            set movy   [expr {$fy   + $vy*$dsample}]
            set movz   [expr {$fz   + $vz*$dsample}]
            set movlvx [expr {$flvx + $vx*$dsample}]
            set movlvy [expr {$flvy + $vy*$dsample}]
            set movlvz [expr {$flvz + $vz*$dsample}]
            lassign [hole::hcapen $movx $movy $movz $movlvx $movlvy $movlvz $n $xs $ys $zs $rs] \
                moveng movcaprad mi1 mi2 md2

            set midx [expr {0.5*($movlvx+$movx)}]
            set midy [expr {0.5*($movlvy+$movy)}]
            set midz [expr {0.5*($movlvz+$movz)}]
            lassign [hole::hcapen $midx $midy $midz $midx $midy $midz $n $xs $ys $zs $rs] \
                mideng midcaprad di1 di2 dd2

            if {$mideng < 2.0*$moveng} {
                # collapse to the degenerate midpoint - much better than the
                # moved-on capsule (holcal.f:722)
                set lowx $midx; set lowy $midy; set lowz $midz
                set lowlvx $midx; set lowlvy $midy; set lowlvz $midz
            } else {
                set lowx $movx; set lowy $movy; set lowz $movz
                set lowlvx $movlvx; set lowlvy $movlvy; set lowlvz $movlvz
            }
        }
    }
    return [list [lsort -real -index 0 $out] $disc]
}

# ==============================================================================
#  Output
# ==============================================================================

proc hole::write_capsule_csv {rows path} {
    set fh [open $path w]
    puts $fh "axial_coordinate_angstrom,eff_radius_angstrom,cap_radius_angstrom,cx,cy,cz,lvx,lvy,lvz"
    foreach row $rows {
        lassign $row t cx cy cz lvx lvy lvz effrad caprad a1 a2
        puts $fh [format "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f" \
            $t $effrad $caprad $cx $cy $cz $lvx $lvy $lvz]
    }
    close $fh
}

# ==============================================================================
#  write_capsule_sph - the CAPSULE .sph record (QC1/QC2)      (wpdbsp.f LCAPS)
# ==============================================================================
#
#  Same discovery-order/duplicate-record-0/ADDEND shape as hole::write_sph
#  (hole.tcl) - see that proc's header for the full derivation - but every
#  stored slice writes TWO atom lines instead of one:
#    QC1  STRCEN (the primary capsule centre, "cx cy cz" here) radius STRBRD
#    QC2  STRLVC (the second capsule centre, "lvx lvy lvz" here) radius STRBRD
#  both with a literal 0.00 second radius field (wpdbsp.f:74-83 - NOT the
#  same value duplicated the way the plain QSS record is; verified against
#  the reference-build hole (see vmdhole/hole_tcl/README.md)'s own -capsul -sphpdb output, see this proc's header
#  test command).
#
#  STRBRD is the RAW capsule radius (hcapen.f's CAPRAD, this port's own
#  `caprad` field - NOT effrad, the area-transformed effective radius CSV
#  output reports) - wpdbsp.f's own comment is explicit about this ("STRBRD
#  is used to store the real capsule radius").
#
#  ADDEND itself (hole::addend, shared with write_sph) takes no LCAPS-aware
#  path in the source - addend.f's own argument list has no capsule-specific
#  input at all, so its terminal grid is plain HOLEEN on LASCEN regardless of
#  LCAPS, and its own emitted records are the ordinary QSS/-888 form, not
#  QC1/QC2 - confirmed directly against the reference-build hole (see vmdhole/hole_tcl/README.md)'s own -capsul
#  -sphpdb output (the ADDEND block there prints 'QSS SPH S-888', not
#  'QCn SPH S-888').
#
#  Byte-identical to the reference-build hole (see vmdhole/hole_tcl/README.md) on 1GRM (cpoint 0 0 4, cvect 0 0 1,
#  sample 0.5, endrad 8.0, raseed 1, capsul):
#    printf 'sphpdb refc.sph\ncoord 1GRM.pdb\nradius simple.rad\n
#      cpoint 0 0 4\ncvect 0 0 1\nsample 0.5\nendrad 8.0\nraseed 1\ncapsul\n' \
#      | the reference-build hole (see vmdhole/hole_tcl/README.md)
#    tclsh hole.tcl -pdb 1GRM.pdb -rad simple.rad -method capsule \
#        -cpoint "0 0 4" -cvect "0 0 1" -sample 0.5 -endrad 8.0 -seed 1 \
#        -sph tclc.sph
#    diff -q refc.sph tclc.sph      # clean
#  See vmdhole/hole_tcl/tests/sph_addend_test.sh for the durable version of this.

# One capsule .sph record pair (QC1 = STRCEN, QC2 = STRLVC, both radius
# STRBRD/0.00 - see hole::write_capsule_sph's header).
proc hole::_qc_sph_lines {fh irec cx cy cz lvx lvy lvz caprad} {
    puts $fh [format "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
        1 "QC1" "SPH" "S" $irec $cx $cy $cz $caprad 0.00]
    puts $fh [format "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
        1 "QC2" "SPH" "S" $irec $lvx $lvy $lvz $caprad 0.00]
}

proc hole::write_capsule_sph {discovery path args} {
    array set o {}
    array set o $args
    foreach req {-n -xs -ys -zs -rs -cvect -sample -endrad} {
        if {![info exists o($req)]} { error "hole::write_capsule_sph: missing $req" }
    }
    set n $o(-n); set xs $o(-xs); set ys $o(-ys); set zs $o(-zs); set rs $o(-rs)
    set endrad $o(-endrad)
    set sample [expr {abs($o(-sample))}]
    lassign $o(-cvect) cvx cvy cvz
    set cvn [expr {sqrt($cvx*$cvx + $cvy*$cvy + $cvz*$cvz)}]
    if {$cvn < 1e-9} { error "hole::write_capsule_sph: CVECT is a zero vector" }
    set cvx [expr {$cvx/$cvn}]; set cvy [expr {$cvy/$cvn}]; set cvz [expr {$cvz/$cvn}]
    set ncvect [list $cvx $cvy $cvz]
    lassign [hole::calper $ncvect] perpve perpvn

    set pos {}; set neg {}
    foreach row $discovery {
        if {[lindex $row 0] >= 0} { lappend pos $row } else { lappend neg $row }
    }
    if {[llength $pos] == 0} {
        error "hole::write_capsule_sph: no +ve-direction slice was stored - nothing to write"
    }

    set fh [open $path w]

    foreach row $pos {
        lassign $row irec t cx cy cz lvx lvy lvz effrad caprad a1 a2
        hole::_qc_sph_lines $fh $irec $cx $cy $cz $lvx $lvy $lvz $caprad
    }
    if {[llength $pos] < 2} {
        close $fh
        error "hole::write_capsule_sph: only [llength $pos] +ve slice(s) stored - ADDEND's STRCEN(,STRNOP-1) reset needs at least 2; not reproduced for this degenerate case"
    }
    lassign [lindex $pos end-1] li lt lcx lcy lcz llvx llvy llvz leffrad lcaprad la1 la2
    foreach pt [hole::addend [list $lcx $lcy $lcz] $sample $ncvect $perpve $perpvn \
                    $endrad $n $xs $ys $zs $rs] {
        lassign $pt ax ay az arad
        puts $fh [hole::_sph_atom_line -888 $ax $ay $az $arad 0.00]
        puts $fh "LAST-REC-END"
    }

    lassign [lindex $pos 0] irec0 t0 cx0 cy0 cz0 lvx0 lvy0 lvz0 effrad0 caprad0 a10 a20
    hole::_qc_sph_lines $fh $irec0 $cx0 $cy0 $cz0 $lvx0 $lvy0 $lvz0 $caprad0

    foreach row $neg {
        lassign $row irec t cx cy cz lvx lvy lvz effrad caprad a1 a2
        hole::_qc_sph_lines $fh $irec $cx $cy $cz $lvx $lvy $lvz $caprad
    }
    if {[llength $neg] < 2} {
        close $fh
        error "hole::write_capsule_sph: only [llength $neg] -ve slice(s) stored - ADDEND's STRCEN(,-STRNON+1) reset needs at least 2; not reproduced for this degenerate case"
    }
    lassign [lindex $neg end-1] li lt lcx lcy lcz llvx llvy llvz leffrad lcaprad la1 la2
    foreach pt [hole::addend [list $lcx $lcy $lcz] [expr {-$sample}] $ncvect $perpve $perpvn \
                    $endrad $n $xs $ys $zs $rs] {
        lassign $pt ax ay az arad
        puts $fh [hole::_sph_atom_line -888 $ax $ay $az $arad 0.00]
        puts $fh "LAST-REC-END"
    }

    close $fh
}
