#!/usr/bin/env tclsh
# ==============================================================================
#  refine.tcl - the post-annealing refinement HOLE applies (hsbxmi.f)
# ==============================================================================
#
#  This file was commissioned to port hsbxmi.f (401 lines) and wire it into
#  hole::anneal_slice's plain spherical path in hole.tcl. Read against the
#  actual call site, that premise does not hold - recorded here so the next
#  reader does not have to re-derive it.
#
#  hsbxmi.f is called from exactly ONE place in the whole source tree:
#
#      holcal.f:484   IF (LSPHBX) THEN
#      holcal.f:493     CALL HSBXMI( LOWCEN, LOWENG, ... )
#      holcal.f:503   ENDIF
#
#  LSPHBX is the SPHERE-BOX option (a SPHPO card), off unless the user asks
#  for it (default LSPHBX=.FALSE., hole.f:399; only rcontr.f:435 sets it
#  true, on a SPHPO card). grepping the full source tree for a refinement
#  call outside that one `IF` block - LCAPS gets its own similar-shaped
#  block right below it, calling HCAPEN, not HSBXMI - turns up nothing:
#
#      $ grep -rn "CALL HSBXMI\|CALL HSBXEN" *.f
#      hsbxmi.f:215,247,262,281,294,320,337:  CALL HSBXEN(...)   <- self-calls
#      holcal.f:307:      CALL HSBXEN( NEWCEN, ...               <- also LSPHBX-gated
#      holcal.f:493:      CALL HSBXMI( LOWCEN, ...               <- the only call
#
#  So for the plain calculation (no SPHPO, no CAPSULE card - the case this
#  file implements and the case in the verification command) holcal.f's
#  DO 10 loop runs straight from the annealing result (LOWCEN, LOWENG) to
#  the "highest radius point found" report with NO minimisation step in
#  between. Confirmed empirically, not just by reading: running the real
#  binary on the plain verification input and grepping its own stdout for
#  the message hsbxmi.f's caller prints before invoking it -
#
#      $ grep -n "Applying sd min" out.txt
#      (no output)
#      $ grep -n "spherebox and capsule options are turned off" out.txt
#       The spherebox and capsule options are turned off.
#
#  - shows the refinement branch never fires. hole::anneal_slice's own
#    return value (in hole.tcl) already IS holcal.f's final answer for a
#    plain slice; nothing is missing there.
#
#  What HSBXMI actually does, for the record, in case SPHPO ever gets
#  ported: crude steepest descent (up to 101 steps, adaptive step length
#  1.2x/0.5x, converges at step<1e-5) over THREE variables - two in-plane
#  displacements plus a rotation angle THETA of the box's long axis - using
#  central-difference derivatives, all through HSBXEN (hsbxen.f), a
#  DIFFERENT objective than HOLEEN: the area of the largest box (not
#  sphere) that fits at a point, biased along a long axis. Porting it
#  faithfully means porting hsbxen.f and the whole LSPHBX state (LVECT,
#  LBOXDM, SBOXDM, BOXRAD) first; there is no version of it that operates
#  on HOLEEN's plain-sphere energy, so it cannot be retrofitted onto
#  hole::anneal_slice's spherical result.
#
#  Per this file's own rule (see hole.tcl section 8): rather than silently
#  no-op or approximate, calling this on the plain path is a hard error.

proc hole::refine_slice {args} {
    error "hole.tcl: refine_slice (hsbxmi.f) is spherebox-only (holcal.f:484\
        IF (LSPHBX)) and is not reachable from the plain spherical path -\
        see vmdhole/hole_tcl/refine.tcl for the derivation. There is nothing to\
        call here for a plain HOLE calculation; do not wire this in."
}
