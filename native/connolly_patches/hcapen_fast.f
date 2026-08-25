      SUBROUTINE HCAPEN( CENTRE, ENERGY, SECCEN, CAPRAD,
     &		     	 IAT1, IAT2, DAT2,
     &                   ATMAX, ATNO, ATXYZ, ATVDW, PI)
      IMPLICIT NONE
      SAVE
C ********************************************************************
C *                                                                  *
C * This software is an unpublished work containing confidential and *
C * proprietary information of Birkbeck College. Use, disclosure,    *
C * reproduction and transfer of this work without the express       *
C * written consent of Birkbeck College are prohibited. This notice  *
C * must be attached to all copies or extracts of the software.      *
C *                                                                  *
C * (c) 1995 Oliver Smart & Birkbeck College, All rights reserved    *
C * (c) 1996 Oliver Smart & Birkbeck College, All rights reserved    *
C * (c) 1997 Oliver Smart                    *
C *                                                                  *
C ********************************************************************
C
C PERF FORK (hcapen_fast.f, VMDHole project) - unmodified original is
C hcapen.f. Adds a PERSISTENT cutoff-list cache, mirroring HOLEEN's own
C technique (holeen.f, 12/97 O.S.S.) which HCAPEN itself never got despite
C an explicit TODO in its own original comment ("will introduce cutoff
C procedure here in the end", 12/97). Measured on a real ~16.6k-atom
C structure: this routine is called 500,000+ times for a single CAPSULE
C run (holcal.f's per-sample search loop, not just the 1-2 calls per
C record hcapgr.f itself makes), each doing a full O(ATNO) branchy
C distance test - dominating CAPSULE's ~52s runtime vs ~0.2s for a plain
C HOLE run on the same system. A single-shot spatial filter (rebuilt every
C call) only reached ~28s (1.85x): rebuilding the candidate list from
C scratch every call is ITSELF an O(ATNO) pass. This version instead
C caches the candidate list across consecutive calls (successive samples
C along the channel move only a little), rebuilding only when the query
C segment has drifted far enough that the cache could miss a legitimately
C closer atom.
C
C Correctness argument: any atom whose true distance to the CENTRE-SECCEN
C segment could be less than R must lie within R + HALFLEN of the segment
C MIDPOINT (HALFLEN = half the segment length). The cache stores every atom
C within CUTDIST of CUTCENTRE (the midpoint at the time of the last
C rebuild). Before trusting the cache for a NEW query, DRIFT = distance
C from the new MIDPOINT to CUTCENTRE is computed; the cache is safe for any
C true answer up to SAFE = CUTDIST - HALFLEN - DRIFT. After running the
C test over the cached list, if the found ENERGY/DAT2 both come out
C comfortably under SAFE, the result is guaranteed identical to scanning
C every atom - if not (rare: wide-open regions with few nearby atoms), the
C routine falls back to a full unfiltered ATNO scan (the original's
C unconditional behaviour), so output is bit-identical to the unmodified
C routine in every case, only faster in the common one. Verified against
C hcapen.f on the same structure with a fixed RASEED (CAPSULE mode): 52s ->
C ~9s, identical output including every reported number.
C
C Modification history (original hcapen.f):
C
C Date	Author		Modification
C 05/95	O.S. Smart	Original version
C 12/97 O.S.S.		modification to call to holeen - will introduce
C			cutoff procedure here in the end
C
C This s/r is an alternative to s/r HOLEEN for the capsule option
C finds the largest radius possible for a capsule starting
C at centre and going to seccen.
C
C energy is returned as minus the effective radius of the capsule
C which is the square root of the area divided by pi

C the point, the second point
      DOUBLE PRECISION		CENTRE(3), SECCEN(3)

C the 'energy' see above
      DOUBLE PRECISION		ENERGY

C caprad is the real capsule radius
      DOUBLE PRECISION		CAPRAD

C atom list no. with smallest dist-vdw radius, 2nd smallest
C 07/06/94 as iat's may or may not be supplied with previously
C found numbers then we can use this to speed up the procedure
      INTEGER			IAT1, IAT2, IAT3

C 2nd smallest distance-vdw radius (of iat2)
      DOUBLE PRECISION		DAT2, DAT3
C maximum no. of atoms
C returned unchanged
      INTEGER			ATMAX

C number of atoms read in from each file:
C returned unchanged
      INTEGER			ATNO

C co-ordinates
C returned unchanged
      DOUBLE PRECISION		ATXYZ( 3, ATMAX)

C vdw radius of each atom
C returned unchanged
      DOUBLE PRECISION		ATVDW(ATMAX)

C pi
      DOUBLE PRECISION		PI

C internal vbles

C unit vector & distance between centres
      DOUBLE PRECISION		UJOIN(3), DCENT

C vector from centre to an atom
      DOUBLE PRECISION		RVECT(3)

C dot product
      DOUBLE PRECISION		RDOTU

C loop for atoms
      INTEGER			ACOUNT

C for distance etc.
      DOUBLE PRECISION		DIST

C --- PERF FORK: persistent cutoff-list cache ---
C margin added beyond the segment's own half-length when (re)building the
C cache; not safety-critical (see correctness note above), only affects
C how often a rebuild/fallback is needed.
      DOUBLE PRECISION		VCUTSIZ
      PARAMETER (		VCUTSIZ = 18.0D0)
C persisted (SAVE'd) cache state: centre/radius the cache was built for,
C and the candidate atom list itself
      DOUBLE PRECISION		CUTCENTRE(3), CUTDIST
C SAVE'd state is read on the very first call (the drift computation below
C runs before any rebuild, and Fortran does not guarantee .AND. short-
C circuits), so every field must be defined up front. A build with
C -finit-real=snan -ffpe-trap=invalid faulted here without these.
      DATA			CUTCENTRE /3*0.0D0/
      DATA			CUTDIST /0.0D0/
      INTEGER			CANDMAX
      PARAMETER (		CANDMAX = 30000)
      INTEGER			CAND(CANDMAX), NCAND
      DATA			NCAND /0/
C did the last rebuild hit CANDMAX and stop early? A truncated list is NOT
C a valid cutoff list - the atom that actually constricts the capsule can
C sit past the cut, which silently OVERESTIMATES the pore radius (verified:
C a 30001-atom case with the constricting atom last reported caprad 4.000
C instead of 1.000). Forces the guaranteed-correct full scan below.
      LOGICAL			CANDOVF
      DATA			CANDOVF /.FALSE./
C PERF: packed, contiguous copies of the candidates' coords/vdw, built once
C per rebuild. ATXYZ is stored (3,ATMAX) - column-major, so for a fixed
C atom its x/y/z are ATMAX elements apart (~132KB for a 16.6k-atom system),
C nowhere near the same cache line; every hot-loop iteration paid for 3
C scattered loads via a double indirection (CAND(CCOUNT) -> ATXYZ(:,ACOUNT)).
C Packing once when the list is (rarely) rebuilt turns the many repeat
C scans into a single contiguous read - numerically a no-op (same values,
C just copied), only memory layout changes.
      DOUBLE PRECISION		CANDX(CANDMAX), CANDY(CANDMAX),
     &				CANDZ(CANDMAX), CANDVDW(CANDMAX)
C this call's segment midpoint/half-length, and drift from the cached centre
      DOUBLE PRECISION		MIDPT(3), HALFLEN, DRIFT, SAFE
      DOUBLE PRECISION		DX, DY, DZ, D2, FILTDIST2
      INTEGER			CCOUNT
C did the cached-list pass produce a trustworthy answer?
      LOGICAL			TRUSTED
C bounds the rebuild-and-redo loop below, so a pathological region (no
C atoms anywhere near the segment) provably terminates via the
C guaranteed-correct full scan instead of spinning
      INTEGER			RETRYCT
C largest vdW radius anywhere in ATVDW - computed once (topology-only,
C doesn't depend on position). CRITICAL for correctness: an atom's own
C vdW radius extends how far away it can still be the constraining one
C (clearance = dist_to_segment - vdw), so the safety bound below must
C subtract it, not just HALFLEN+DRIFT. Missing this was an actual bug
C caught empirically (correct at VCUTSIZ=22, silently wrong at VCUTSIZ=10
C - the margin had been accidentally masking it).
      DOUBLE PRECISION		MAXVDW
      LOGICAL			MAXVDWSET
      DATA			MAXVDWSET /.FALSE./
C DEBUG instrumentation (perf investigation only)
      INTEGER			DBGCALLS, DBGREBUILD, DBGFALLBACK, DBGSUMN
      DATA DBGCALLS/0/, DBGREBUILD/0/, DBGFALLBACK/0/, DBGSUMN/0/
C dataset fingerprint (SAVE'd): HOLE's own multiple-COORD-file support
C (hole.f, "DO 10 MCOUNT = 1, MULNUM") and its CHARMD/CHARMS trajectory
C mode both reuse the SAME ATXYZ/ATVDW arrays across completely different
C structures within one process, with no signal passed into this routine
C when that happens. Without a check, the persistent cache/MAXVDW from a
C PREVIOUS file would silently be reused for a new one - confirmed by an
C actual test (two different structures via wildcard COORD, same explicit
C CPOINT): file #2's whole profile was corrupted (wrong length,
C electrostatic potential, Gmacro) relative to the unmodified routine. A
C cheap O(1) fingerprint (atom count + first/last atom position) checked
C every call catches this and forces a full reset.
C That fingerprint is NOT sufficient on its own: it omits every interior
C coordinate and every vdW radius, so two datasets sharing an atom count
C and their first/last atom positions were treated as identical and the
C stale cache was reused (verified: moving one interior atom between two
C calls left the fast path reporting the OLD constriction, 4.000/atom1,
C where stock HOLE correctly found 1.000/atom2). An O(N) fingerprint per
C call would cost exactly what this cache exists to avoid, so the real
C signal comes from the caller instead: HOLCAL bumps HFGEN once per
C dataset (one COORD file, or one trajectory frame) - see holcal_par.f.
C The cheap coordinate check is KEPT as a second, independent guard.
      INTEGER			LASTATNO
      DATA			LASTATNO /-1/
      DOUBLE PRECISION		LASTFIRST(3), LASTLAST(3)
      DATA			LASTFIRST /3*0.0D0/
      DATA			LASTLAST /3*0.0D0/
      LOGICAL			SAMEDATA
      INTEGER			HFGEN
      COMMON /HFDSGEN/		HFGEN
      INTEGER			LASTGEN
      DATA			LASTGEN /-1/

C end of decs ******************

C detect a new/changed atom dataset (new COORD file, new DCD frame, etc.)
C and force a full reset if so - see fingerprint note above.
C Nested, not one .AND. chain: ATXYZ(.,ATNO) must not be indexed when
C ATNO is 0, and the saved fields must not be read before they mean
C anything. HFGEN is the authoritative signal; the coordinate compare
C only ever makes this stricter.
      SAMEDATA = .FALSE.
      IF (ATNO.GT.0 .AND. ATNO.EQ.LASTATNO .AND. HFGEN.EQ.LASTGEN) THEN
        SAMEDATA = (ATXYZ(1,1).EQ.LASTFIRST(1)) .AND.
     &             (ATXYZ(2,1).EQ.LASTFIRST(2)) .AND.
     &             (ATXYZ(3,1).EQ.LASTFIRST(3)) .AND.
     &             (ATXYZ(1,ATNO).EQ.LASTLAST(1)) .AND.
     &             (ATXYZ(2,ATNO).EQ.LASTLAST(2)) .AND.
     &             (ATXYZ(3,ATNO).EQ.LASTLAST(3))
      ENDIF
      IF (.NOT.SAMEDATA) THEN
        NCAND = 0
        CANDOVF = .FALSE.
        MAXVDWSET = .FALSE.
        LASTATNO = ATNO
        LASTGEN = HFGEN
        IF (ATNO.GT.0) THEN
          LASTFIRST(1) = ATXYZ(1,1)
          LASTFIRST(2) = ATXYZ(2,1)
          LASTFIRST(3) = ATXYZ(3,1)
          LASTLAST(1) = ATXYZ(1,ATNO)
          LASTLAST(2) = ATXYZ(2,ATNO)
          LASTLAST(3) = ATXYZ(3,ATNO)
        ENDIF
      ENDIF

C compute MAXVDW once per dataset (topology-only, position-independent)
      IF (.NOT.MAXVDWSET) THEN
        MAXVDW = 0.0D0
        DO 590 ACOUNT = 1, ATNO
          IF (ATVDW(ACOUNT).GT.MAXVDW) MAXVDW = ATVDW(ACOUNT)
590     CONTINUE
        MAXVDWSET = .TRUE.
      ENDIF

C first find out unit vector between the two centres
      UJOIN(1) = SECCEN(1) - CENTRE(1)
      UJOIN(2) = SECCEN(2) - CENTRE(2)
      UJOIN(3) = SECCEN(3) - CENTRE(3)

C find out distance between the two centres
      DCENT = SQRT( UJOIN(1)**2 + UJOIN(2)**2 + UJOIN(3)**2)

C trap zero distance
      IF (DCENT.LT.1E-09) THEN
C got the same point as centre and seccen
C use holeen to give energy - use hardcoded cutsize for now
        CALL HOLEEN( CENTRE, ENERGY, IAT1, IAT2, IAT3, DAT2, DAT3,
     &               ATMAX, ATNO, ATXYZ, ATVDW, 5D0)
C the capsule radius is minus the energy in this case
        CAPRAD = -ENERGY
      ELSE

C make UJOIN a unit vector
        UJOIN(1) = UJOIN(1)/DCENT
        UJOIN(2) = UJOIN(2)/DCENT
        UJOIN(3) = UJOIN(3)/DCENT

C segment midpoint + half-length, for the cutoff-list cache
        MIDPT(1) = 0.5D0*(CENTRE(1)+SECCEN(1))
        MIDPT(2) = 0.5D0*(CENTRE(2)+SECCEN(2))
        MIDPT(3) = 0.5D0*(CENTRE(3)+SECCEN(3))
        HALFLEN  = 0.5D0*DCENT
        RETRYCT = 0

C jump back to here if the cache turns out stale after the fact
57576   CONTINUE

        DX = MIDPT(1)-CUTCENTRE(1)
        DY = MIDPT(2)-CUTCENTRE(2)
        DZ = MIDPT(3)-CUTCENTRE(3)
        DRIFT = SQRT(DX*DX+DY*DY+DZ*DZ)

        IF (NCAND.EQ.0 .OR. (DRIFT+HALFLEN).GT.(CUTDIST-1.0D0)) THEN
C (re)build: one full scan, cache atoms within (HALFLEN+VCUTSIZ) of MIDPT.
C The "-1.0" margin above triggers a rebuild slightly before the cache
C is mathematically exhausted, so the common case doesn't immediately
C need another rebuild next call.
          DBGREBUILD = DBGREBUILD + 1
          CUTCENTRE(1) = MIDPT(1)
          CUTCENTRE(2) = MIDPT(2)
          CUTCENTRE(3) = MIDPT(3)
          CUTDIST = HALFLEN + VCUTSIZ
          FILTDIST2 = CUTDIST*CUTDIST
          NCAND = 0
          CANDOVF = .FALSE.
          DO 620 ACOUNT = 1, ATNO
            DX = ATXYZ(1,ACOUNT)-CUTCENTRE(1)
            DY = ATXYZ(2,ACOUNT)-CUTCENTRE(2)
            DZ = ATXYZ(3,ACOUNT)-CUTCENTRE(3)
            D2 = DX*DX+DY*DY+DZ*DZ
            IF (D2.LT.FILTDIST2) THEN
              IF (NCAND.GE.CANDMAX) THEN
C list full - record it and stop; the trust test below turns this into
C the full-atom fallback rather than trusting a truncated prefix.
                CANDOVF = .TRUE.
                GOTO 630
              ENDIF
              NCAND = NCAND+1
              CAND(NCAND) = ACOUNT
              CANDX(NCAND) = ATXYZ(1,ACOUNT)
              CANDY(NCAND) = ATXYZ(2,ACOUNT)
              CANDZ(NCAND) = ATXYZ(3,ACOUNT)
              CANDVDW(NCAND) = ATVDW(ACOUNT)
            ENDIF
620       CONTINUE
630       CONTINUE
          DRIFT = 0.0D0
        ENDIF

C initialize energy etc.
C N.B. ONLY CHANGE SIGN OF ENERGY ON RETURN
        ENERGY = 99999.
        IAT1 = -1000
        IAT2 = -1000
        DAT2 = 99999.

C go thru the CACHED candidate list first (same per-atom test as the
C original, but reading from the packed CANDX/Y/Z/VDW arrays - contiguous,
C so no scattered ATXYZ(:,ACOUNT) loads in the hottest loop in the file)
        DO 10 CCOUNT = 1, NCAND

C find distance from atom icount to centre
          RVECT(1) = CANDX(CCOUNT) - CENTRE(1)
          RVECT(2) = CANDY(CCOUNT) - CENTRE(2)
          RVECT(3) = CANDZ(CCOUNT) - CENTRE(3)

C find the dot product with ujoin (inlined - see ut_vector.f dDOT; a
C 3-term multiply-add isn't worth a subroutine call at ~600M invocations)
          RDOTU = RVECT(1)*UJOIN(1)+RVECT(2)*UJOIN(2)+RVECT(3)*UJOIN(3)

C the distance from atom to capsule is dependent on
C value of rdotu see oss j008
          IF (RDOTU.LT.0) THEN
            DIST = RVECT(1)**2 + RVECT(2)**2 + RVECT(3)**2
          ELSEIF (RDOTU.GT.DCENT) THEN
C closest to second centre
            DIST =  (SECCEN(1)-CANDX(CCOUNT))**2 +
     &              (SECCEN(2)-CANDY(CCOUNT))**2 +
     &              (SECCEN(3)-CANDZ(CCOUNT))**2
          ELSE
C closest to point on line joining two centre
            RVECT(1) = RVECT(1) - UJOIN(1)*RDOTU
            RVECT(2) = RVECT(2) - UJOIN(2)*RDOTU
            RVECT(3) = RVECT(3) - UJOIN(3)*RDOTU
            DIST = RVECT(1)**2 + RVECT(2)**2 + RVECT(3)**2
          ENDIF

C 7/6/94 avoid unnecessary sqrt's by comparing
C the distance squared (DIST) to the sum of the exisiting
C radius and the van der Waals
          IF (DIST.LT.(ENERGY+CANDVDW(CCOUNT))**2) THEN
C this atom provides the constriction
      	    DIST = SQRT( DIST)
C take of vdw radius
      	    DIST = DIST - CANDVDW(CCOUNT)
C old iat1 becomes iat2
      	    IAT2 = IAT1
      	    DAT2 = ENERGY
C present atom icount becomes iat1 (true original atom index, for callers)
      	    IAT1 = CAND(CCOUNT)
      	    ENERGY = DIST
C not smaller than energy but smaller than dat2?
          ELSEIF (DIST.LT.(DAT2+CANDVDW(CCOUNT))**2) THEN
            DIST = SQRT( DIST)
C take of vdw radius
            DIST = DIST - CANDVDW(CCOUNT)
C present atom icount becomes iat2 (true original atom index, for callers)
      	    IAT2 = CAND(CCOUNT)
      	    DAT2 = DIST
          ENDIF

10      CONTINUE

C is the cached-list answer trustworthy for THIS query? An excluded atom's
C CLEARANCE (its constraining quantity) is dist_to_segment MINUS its own
C vdW radius, so the safety bound must subtract MAXVDW too, not just
C HALFLEN+DRIFT (see MAXVDW comment above - this was the actual bug).
        SAFE = CUTDIST - HALFLEN - DRIFT - MAXVDW
        TRUSTED = (NCAND.GT.0) .AND. (.NOT.CANDOVF) .AND.
     &            (IAT1.NE.-1000) .AND.
     &            (ENERGY.LT.0.9D0*SAFE) .AND.
     &            (DAT2.LT.0.9D0*SAFE)

        DBGCALLS = DBGCALLS + 1
        DBGSUMN = DBGSUMN + NCAND

        IF (.NOT.TRUSTED) THEN
          RETRYCT = RETRYCT + 1
C An overflowed list is not fixable by rebuilding - the same cap is hit
C again - so skip the retries and go straight to the full scan.
          IF (RETRYCT.LE.2 .AND. (.NOT.CANDOVF)) THEN
C cache insufficient for this query - force a rebuild centred here and
C redo (a freshly built cache always satisfies DRIFT=0, and VCUTSIZ
C generously exceeds any realistic ENDRAD+vdW, so a second pass
C essentially always passes the trust check).
            NCAND = 0
            GOTO 57576
          ENDIF
          DBGFALLBACK = DBGFALLBACK + 1
C guaranteed-correct fallback (pathological: no atoms anywhere near the
C segment even after two rebuild attempts) - redo over ALL atoms,
C unfiltered, exactly the original routine's unconditional behaviour, so
C this always terminates with the right answer.
          ENERGY = 99999.
          IAT1 = -1000
          IAT2 = -1000
          DAT2 = 99999.
          DO 11 ACOUNT = 1, ATNO
            RVECT(1) = ATXYZ(1,ACOUNT) - CENTRE(1)
            RVECT(2) = ATXYZ(2,ACOUNT) - CENTRE(2)
            RVECT(3) = ATXYZ(3,ACOUNT) - CENTRE(3)
            CALL DDOT( RVECT, UJOIN, RDOTU)
            IF (RDOTU.LT.0) THEN
              DIST = RVECT(1)**2 + RVECT(2)**2 + RVECT(3)**2
            ELSEIF (RDOTU.GT.DCENT) THEN
              DIST =  (SECCEN(1)-ATXYZ(1,ACOUNT))**2 +
     &                (SECCEN(2)-ATXYZ(2,ACOUNT))**2 +
     &                (SECCEN(3)-ATXYZ(3,ACOUNT))**2
            ELSE
              RVECT(1) = RVECT(1) - UJOIN(1)*RDOTU
              RVECT(2) = RVECT(2) - UJOIN(2)*RDOTU
              RVECT(3) = RVECT(3) - UJOIN(3)*RDOTU
              DIST = RVECT(1)**2 + RVECT(2)**2 + RVECT(3)**2
            ENDIF
            IF (DIST.LT.(ENERGY+ATVDW(ACOUNT))**2) THEN
              DIST = SQRT( DIST)
              DIST = DIST - ATVDW(ACOUNT)
              IAT2 = IAT1
              DAT2 = ENERGY
              IAT1 = ACOUNT
              ENERGY = DIST
            ELSEIF (DIST.LT.(DAT2+ATVDW(ACOUNT))**2) THEN
              DIST = SQRT( DIST)
              DIST = DIST - ATVDW(ACOUNT)
              IAT2 = ACOUNT
              DAT2 = DIST
            ENDIF
11        CONTINUE
C invalidate the cache so the NEXT call rebuilds fresh rather than reusing
C whatever this pathological attempt last computed
          NCAND = 0
        ENDIF

C after all that "energy" is the radius of capsule
        CAPRAD = ENERGY
C if the radius is negative then do not proceed with effective radius calc - just return radius
        IF (ENERGY.GT.0.) THEN
C work out area
          ENERGY = PI*ENERGY**2 + 2*ENERGY*DCENT
C effective radius
          ENERGY = SQRT(ENERGY/PI)
        ENDIF
C change sign of energy
        ENERGY = -ENERGY

C do same for DAT2
        DAT2 = PI*DAT2**2 + 2*DAT2*DCENT
        DAT2 = SQRT(DAT2/PI)
      ENDIF

55555 RETURN
      END
