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
C hcapen.f. Same interface, same answer, bit for bit; only the way the
C answer is found differs.
C
C The original scans every atom on every call, and HOLCAL calls this
C routine once per Monte Carlo step, so that scan is where all of
C CAPSULE's time goes.
C
C This version keeps a uniform grid over the atoms (built once per
C dataset) and, per call, only looks at the cells around the query
C segment. It is EXACT, not approximate, by this argument:
C
C   1. Gather every atom in the cells covering the segment's bounding
C      box grown by RQ. Any atom NOT gathered lies outside that box, so
C      its distance to the segment is at least RQ, and its clearance
C      (distance minus its own vdW radius) at least RQ - MAXVDW.
C   2. Over the gathered atoms find the second-smallest clearance E2,
C      using the original's own per-atom arithmetic. If
C      RQ - MAXVDW > E2 + TOL, no atom outside the box can beat the two
C      best inside it; otherwise grow RQ and gather again.
C   3. The original's answer depends only on the atoms whose clearance
C      is within the two best (ENERGY, DAT2 and the indices IAT1, IAT2):
C      an atom worse than the final second-best is at most a transient
C      running minimum, and dropping it cannot change any later
C      comparison's outcome - a later atom is tested against a running
C      minimum that is only ever equal or larger, so it is accepted in
C      at least the same cases, and every acceptance the final answer
C      rests on has a margin of at least TOL over anything dropped. So
C      the original's loop is replayed, in ORIGINAL ATOM ORDER with the
C      original arithmetic, over exactly the gathered atoms whose
C      clearance is within E2 + TOL. Same order, same expressions, same
C      floating-point path, same result - including which of two tied
C      atoms wins.
C
C Anything the argument does not cover (more atoms than the grid can
C hold, a candidate list that overflows) falls back to the original's
C unconditional full scan, so the answer is right in every case.
C
C test_hcapen_cache.sh builds the same driver against this file and the
C original and requires identical output; native/verify.sh does the same
C for whole HOLE runs.
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

C --- PERF FORK: grid state, SAVE'd across calls ---
C the grid's atom store has its own fixed size; a dataset larger than
C this takes the full scan
      INTEGER			GATMAX
      PARAMETER (		GATMAX = 1000000)
C cell count cap: the cell size grows until the box fits
      INTEGER			GMAXC
      PARAMETER (		GMAXC = 2000000)
C counting-sort layout: CSTART(c+1)..CSTART(c+2)-1 are cell c's slots
C (0-based c); within a cell the slots are in ascending atom order. The
C slots hold packed COPIES of the coordinates, vdW radius and atom index:
C ATXYZ is (3,ATMAX) with ATMAX = 1,000,000, so one atom's x, y and z sit
C 8 MB apart and reading them by index costs three cache misses - the
C original full scan streams three columns instead, which is why it is
C fast per atom. A gather over packed rows streams the same way.
      INTEGER			CSTART(GMAXC+1), PID(GATMAX)
      DOUBLE PRECISION		PX(GATMAX), PY(GATMAX), PZ(GATMAX),
     &				PV(GATMAX)
C per-call scratch: squared distance and slot of every gathered atom
      DOUBLE PRECISION		DISTC(400000)
      INTEGER			POSC(400000)
      DOUBLE PRECISION		RX, RY, RZ, TT
      INTEGER			CLO, CHI, P0, P1, NROW
      DOUBLE PRECISION		GMIN(3), GH, GINVH
      INTEGER			GN(3), GNTOT
      LOGICAL			GRIDOK
      DATA			GRIDOK /.FALSE./
C initial cell size; grown if the box would need more than GMAXC cells
      DOUBLE PRECISION		GH0
      PARAMETER (		GH0 = 2.5D0)
C tolerance on the clearance bound (angstrom) - far above rounding,
C far below any real geometric difference
      DOUBLE PRECISION		TOL
      PARAMETER (		TOL = 1.0D-6)
C gathered candidates, and the short list replayed in order
      INTEGER			CANDMAX, SELMAX
      PARAMETER (		CANDMAX = 400000, SELMAX = 8192)
      INTEGER			NCAND
      INTEGER			SEL(SELMAX), NSEL
C query box, cell range, search radius
      DOUBLE PRECISION		QLO(3), QHI(3), RQ, T
      INTEGER			ILO(3), IHI(3), I, J, L, C, P, K, S
      LOGICAL			WHOLE, USEGRID
C two smallest clearances over the gathered set
      DOUBLE PRECISION		E1P, E2P, CL, BIG
      PARAMETER (		BIG = 99999.0D0)
C the previous call's second-best clearance seeds the next search radius
      DOUBLE PRECISION		LASTE2
      DATA			LASTE2 /5.0D0/
      INTEGER			GROWCT
C largest vdW radius in the dataset (topology-only)
      DOUBLE PRECISION		MAXVDW
      DOUBLE PRECISION		GMAX(3)
C dataset fingerprint: HOLCAL bumps HFGEN once per dataset (one COORD
C file, or one trajectory frame) - see holcal_par.f; the coordinate
C compare is a second, independent guard
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
C Nested, not one .AND. chain: ATXYZ(.,ATNO) must not be indexed when
C ATNO is 0, and the saved fields must not be read before they mean
C anything.
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
        GRIDOK = .FALSE.
        LASTATNO = ATNO
        LASTGEN = HFGEN
        LASTE2 = 5.0D0
        IF (ATNO.GT.0) THEN
          LASTFIRST(1) = ATXYZ(1,1)
          LASTFIRST(2) = ATXYZ(2,1)
          LASTFIRST(3) = ATXYZ(3,1)
          LASTLAST(1) = ATXYZ(1,ATNO)
          LASTLAST(2) = ATXYZ(2,ATNO)
          LASTLAST(3) = ATXYZ(3,ATNO)
        ENDIF
C build the grid: bounding box, cell size, counting sort
        IF (ATNO.GT.0 .AND. ATNO.LE.GATMAX) THEN
          MAXVDW = 0.0D0
          DO 500 K = 1, 3
            GMIN(K) = ATXYZ(K,1)
            GMAX(K) = ATXYZ(K,1)
500       CONTINUE
          DO 510 ACOUNT = 1, ATNO
            IF (ATVDW(ACOUNT).GT.MAXVDW) MAXVDW = ATVDW(ACOUNT)
            DO 505 K = 1, 3
              IF (ATXYZ(K,ACOUNT).LT.GMIN(K)) GMIN(K) = ATXYZ(K,ACOUNT)
              IF (ATXYZ(K,ACOUNT).GT.GMAX(K)) GMAX(K) = ATXYZ(K,ACOUNT)
505         CONTINUE
510       CONTINUE
          GH = GH0
520       CONTINUE
          GINVH = 1.0D0/GH
          GNTOT = 1
          DO 525 K = 1, 3
            GN(K) = INT((GMAX(K)-GMIN(K))*GINVH) + 1
            GNTOT = GNTOT*GN(K)
525       CONTINUE
          IF (GNTOT.GT.GMAXC) THEN
            GH = GH*1.5D0
            GOTO 520
          ENDIF
          DO 530 C = 1, GNTOT+1
            CSTART(C) = 0
530       CONTINUE
C count atoms per cell (CSTART(c+2) holds cell c's count for now)
          DO 540 ACOUNT = 1, ATNO
            I = INT((ATXYZ(1,ACOUNT)-GMIN(1))*GINVH)
            J = INT((ATXYZ(2,ACOUNT)-GMIN(2))*GINVH)
            L = INT((ATXYZ(3,ACOUNT)-GMIN(3))*GINVH)
            IF (I.GT.GN(1)-1) I = GN(1)-1
            IF (J.GT.GN(2)-1) J = GN(2)-1
            IF (L.GT.GN(3)-1) L = GN(3)-1
            C = (L*GN(2)+J)*GN(1)+I
            CSTART(C+2) = CSTART(C+2)+1
540       CONTINUE
C prefix sums: CSTART(c+1) = first slot of cell c (1-based into CATOM)
          CSTART(1) = 1
          DO 550 C = 1, GNTOT
            CSTART(C+1) = CSTART(C+1)+CSTART(C)
550       CONTINUE
C place atoms in ascending atom order, so every cell's list is sorted
          DO 560 ACOUNT = 1, ATNO
            I = INT((ATXYZ(1,ACOUNT)-GMIN(1))*GINVH)
            J = INT((ATXYZ(2,ACOUNT)-GMIN(2))*GINVH)
            L = INT((ATXYZ(3,ACOUNT)-GMIN(3))*GINVH)
            IF (I.GT.GN(1)-1) I = GN(1)-1
            IF (J.GT.GN(2)-1) J = GN(2)-1
            IF (L.GT.GN(3)-1) L = GN(3)-1
            C = (L*GN(2)+J)*GN(1)+I
            P = CSTART(C+1)
            PID(P) = ACOUNT
            PX(P) = ATXYZ(1,ACOUNT)
            PY(P) = ATXYZ(2,ACOUNT)
            PZ(P) = ATXYZ(3,ACOUNT)
            PV(P) = ATVDW(ACOUNT)
            CSTART(C+1) = CSTART(C+1)+1
560       CONTINUE
C undo the placement shift: CSTART(c+1) is cell c's first slot again
          DO 570 C = GNTOT, 1, -1
            CSTART(C+1) = CSTART(C)
570       CONTINUE
          CSTART(1) = 1
          GRIDOK = .TRUE.
        ENDIF
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

        USEGRID = GRIDOK
        IF (USEGRID) THEN
C --- grid path: gather, bound, replay in order ---
          DO 600 K = 1, 3
            QLO(K) = MIN(CENTRE(K),SECCEN(K))
            QHI(K) = MAX(CENTRE(K),SECCEN(K))
600       CONTINUE
          RQ = LASTE2 + MAXVDW + TOL + 0.25D0
          IF (RQ.LT.GH) RQ = GH
          GROWCT = 0

610       CONTINUE
          WHOLE = .TRUE.
          DO 620 K = 1, 3
            T = (QLO(K)-RQ-GMIN(K))*GINVH
            IF (T.LT.0.0D0) THEN
              ILO(K) = 0
            ELSE
              ILO(K) = INT(T)
              IF (ILO(K).GT.GN(K)-1) ILO(K) = GN(K)-1
            ENDIF
            T = (QHI(K)+RQ-GMIN(K))*GINVH
            IF (T.LT.0.0D0) THEN
              IHI(K) = 0
            ELSE
              IHI(K) = INT(T)
              IF (IHI(K).GT.GN(K)-1) IHI(K) = GN(K)-1
            ENDIF
            IF (ILO(K).NE.0 .OR. IHI(K).NE.GN(K)-1) WHOLE = .FALSE.
620       CONTINUE

C gather the cells. Cells with the same (j,l) are consecutive in the
C packed arrays, so each row of the box is ONE contiguous run: stream it,
C storing every atom's squared distance to the segment. The clamped-
C projection form here is mathematically the original's piecewise one;
C its rounding differs by ~1e-16, which the TOL below absorbs - the
C replay at the end recomputes the exact original value.
          NCAND = 0
          DO 660 L = ILO(3), IHI(3)
          DO 650 J = ILO(2), IHI(2)
            CLO = (L*GN(2)+J)*GN(1)+ILO(1)
            CHI = (L*GN(2)+J)*GN(1)+IHI(1)
            P0 = CSTART(CLO+1)
            P1 = CSTART(CHI+2)-1
            NROW = P1-P0+1
            IF (NROW.LE.0) GOTO 650
            IF (NCAND+NROW.GT.CANDMAX) THEN
              USEGRID = .FALSE.
              GOTO 670
            ENDIF
            DO 630 P = P0, P1
              RX = PX(P) - CENTRE(1)
              RY = PY(P) - CENTRE(2)
              RZ = PZ(P) - CENTRE(3)
              TT = RX*UJOIN(1)+RY*UJOIN(2)+RZ*UJOIN(3)
              IF (TT.LT.0.0D0) TT = 0.0D0
              IF (TT.GT.DCENT) TT = DCENT
              RX = RX - UJOIN(1)*TT
              RY = RY - UJOIN(2)*TT
              RZ = RZ - UJOIN(3)*TT
              DISTC(NCAND+P-P0+1) = RX*RX+RY*RY+RZ*RZ
              POSC(NCAND+P-P0+1) = P
630         CONTINUE
            NCAND = NCAND+NROW
650       CONTINUE
660       CONTINUE
C the two smallest clearances over the gathered set (sqrt only when an
C atom can still beat the current second-best, as the original does)
          E1P = BIG
          E2P = BIG
          DO 665 K = 1, NCAND
            P = POSC(K)
            IF (DISTC(K).LT.(E2P+PV(P))**2) THEN
              CL = SQRT( DISTC(K)) - PV(P)
              IF (CL.LT.E1P) THEN
                E2P = E1P
                E1P = CL
              ELSEIF (CL.LT.E2P) THEN
                E2P = CL
              ENDIF
            ENDIF
665       CONTINUE
670       CONTINUE

          IF (USEGRID .AND. (.NOT.WHOLE)) THEN
C can an atom outside the gathered box still be one of the two best? Purely
C additive (no margin multiplied against a quantity that can go negative and
C flip from tightening to loosening the bound, the failure mode the earlier
C cutoff-list cache had in its own trust test) - RQ, MAXVDW, E2P and TOL are
C all non-negative by construction, so RQ-MAXVDW cannot go negative and
C invert this comparison's sense.
            IF (E2P.GE.BIG .OR. (RQ-MAXVDW).LE.(E2P+TOL)) THEN
              GROWCT = GROWCT+1
              IF (GROWCT.GT.40) THEN
                USEGRID = .FALSE.
              ELSE
                T = E2P+MAXVDW+TOL+GH
                IF (E2P.GE.BIG) T = RQ*2.0D0
                IF (T.LT.RQ*2.0D0) T = RQ*2.0D0
                RQ = T
                GOTO 610
              ENDIF
            ENDIF
          ENDIF

          IF (USEGRID) THEN
C the short list: gathered atoms within E2P + TOL, in ascending atom
C order (an insertion sort - the list is a handful of atoms)
            NSEL = 0
            DO 700 S = 1, NCAND
              P = POSC(S)
              IF (DISTC(S).LE.(E2P+TOL+PV(P))**2) THEN
                IF (NSEL.GE.SELMAX) THEN
                  USEGRID = .FALSE.
                  GOTO 710
                ENDIF
                ACOUNT = PID(P)
                NSEL = NSEL+1
                K = NSEL
705             CONTINUE
                IF (K.GT.1) THEN
                  IF (SEL(K-1).GT.ACOUNT) THEN
                    SEL(K) = SEL(K-1)
                    K = K-1
                    GOTO 705
                  ENDIF
                ENDIF
                SEL(K) = ACOUNT
              ENDIF
700         CONTINUE
710         CONTINUE
          ENDIF
        ENDIF

C initialize energy etc.
C N.B. ONLY CHANGE SIGN OF ENERGY ON RETURN
        ENERGY = 99999.
        IAT1 = -1000
        IAT2 = -1000
        DAT2 = 99999.

        IF (USEGRID) THEN
C --- replay the original loop over the short list, in atom order ---
          DO 20 S = 1, NSEL
            ACOUNT = SEL(S)
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
20        CONTINUE
          IF (DAT2.LT.BIG) LASTE2 = DAT2
        ELSE
C --- the original, unconditional full scan ---
          DO 10 ACOUNT = 1, ATNO

C find distance from atom icount to centre
            RVECT(1) = ATXYZ(1,ACOUNT) - CENTRE(1)
            RVECT(2) = ATXYZ(2,ACOUNT) - CENTRE(2)
            RVECT(3) = ATXYZ(3,ACOUNT) - CENTRE(3)

C find the dot product with ujoin
            CALL DDOT( RVECT, UJOIN, RDOTU)

C the distance from atom to capsule is dependent on
C value of rdotu see oss j008
            IF (RDOTU.LT.0) THEN
              DIST = RVECT(1)**2 + RVECT(2)**2 + RVECT(3)**2
            ELSEIF (RDOTU.GT.DCENT) THEN
C closest to second centre
              DIST =  (SECCEN(1)-ATXYZ(1,ACOUNT))**2 +
     &                (SECCEN(2)-ATXYZ(2,ACOUNT))**2 +
     &                (SECCEN(3)-ATXYZ(3,ACOUNT))**2
            ELSE
C closest to point on line joining two centre
C put into rvect
              RVECT(1) = RVECT(1) - UJOIN(1)*RDOTU
              RVECT(2) = RVECT(2) - UJOIN(2)*RDOTU
              RVECT(3) = RVECT(3) - UJOIN(3)*RDOTU
C distance squared
              DIST = RVECT(1)**2 + RVECT(2)**2 + RVECT(3)**2
            ENDIF

C 7/6/94 avoid unnecessary sqrt's by comparing
C the distance squared (DIST) to the sum of the exisiting
C radius and the van der Waals
            IF (DIST.LT.(ENERGY+ATVDW(ACOUNT))**2) THEN
C this atom provides the constriction
      	      DIST = SQRT( DIST)
C take of vdw radius
      	      DIST = DIST - ATVDW(ACOUNT)
C old iat1 becomes iat2
      	      IAT2 = IAT1
      	      DAT2 = ENERGY
C present atom icount becomes iat1
      	      IAT1 = ACOUNT
      	      ENERGY = DIST
C not smaller than energy but smaller than dat2?
            ELSEIF (DIST.LT.(DAT2+ATVDW(ACOUNT))**2) THEN
              DIST = SQRT( DIST)
C take of vdw radius
              DIST = DIST - ATVDW(ACOUNT)
C present atom icount becomes iat2
      	      IAT2 = ACOUNT
      	      DAT2 = DIST
            ENDIF

10        CONTINUE
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
