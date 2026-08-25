C Regression driver for the two HCAPEN cutoff-cache defects found in the
C 2026-08-12 code review. Calls HCAPEN directly (not through HOLE) because
C both defects live in SAVE'd state that only shows up across calls.
C
C Built twice by test_hcapen_cache.sh - once against hcapen_fast.f, once
C against the stock hcapen.f - and the two outputs must agree.
C
C   CASE 1  overflow: more than CANDMAX(=30000) atoms inside the cutoff
C           sphere, with the REAL constricting atom last. A truncated
C           candidate list silently overestimates the capsule radius.
C   CASE 2  stale reuse: two calls whose atom count and FIRST/LAST atom
C           coordinates are identical, but with an interior atom moved.
C           The O(1) fingerprint cannot see the move; only the caller's
C           dataset generation (HFGEN) can.
      PROGRAM TCACHE
      IMPLICIT NONE
      INTEGER           ATMAX
      PARAMETER (       ATMAX = 40000)
      DOUBLE PRECISION  ATXYZ(3,ATMAX), ATVDW(ATMAX)
      DOUBLE PRECISION  CENTRE(3), SECCEN(3), ENERGY, CAPRAD, DAT2
      INTEGER           IAT1, IAT2, ATNO, I
      DOUBLE PRECISION  PI
      INTEGER           HFGEN
      COMMON /HFDSGEN/  HFGEN

      PI = 3.14159265358979D0
      HFGEN = 0

C ---- CASE 1: 30001 candidates, constricting atom LAST ----------------
C All atoms sit on a shell 4.0 A from the origin (clearance 4.0 - 0.0),
C except the last, which sits 1.0 A away. A correct search finds 1.0.
      ATNO = 30001
      DO 10 I = 1, ATNO-1
        ATXYZ(1,I) = 4.0D0*COS(DBLE(I)*0.001D0)
        ATXYZ(2,I) = 4.0D0*SIN(DBLE(I)*0.001D0)
        ATXYZ(3,I) = 0.0D0
        ATVDW(I)   = 0.0D0
10    CONTINUE
      ATXYZ(1,ATNO) = 1.0D0
      ATXYZ(2,ATNO) = 0.0D0
      ATXYZ(3,ATNO) = 0.0D0
      ATVDW(ATNO)   = 0.0D0
      CENTRE(1) = 0.0D0
      CENTRE(2) = 0.0D0
      CENTRE(3) = 0.0D0
      SECCEN(1) = 0.0D0
      SECCEN(2) = 0.0D0
      SECCEN(3) = 0.1D0
      HFGEN = 1
      CALL HCAPEN( CENTRE, ENERGY, SECCEN, CAPRAD, IAT1, IAT2, DAT2,
     &             ATMAX, ATNO, ATXYZ, ATVDW, PI)
      WRITE(*,900) 'CASE1', CAPRAD, IAT1

C ---- CASE 2: same count, same first/last atom, interior atom moved ---
C Atom 1 and atom ATNO are the fingerprint HCAPEN used to sample; atom 2
C is the one that actually moves in and constricts the capsule.
      ATNO = 100
      DO 20 I = 1, ATNO
        ATXYZ(1,I) = 4.0D0
        ATXYZ(2,I) = DBLE(I)*0.001D0
        ATXYZ(3,I) = 0.0D0
        ATVDW(I)   = 0.0D0
20    CONTINUE
      HFGEN = 2
      CALL HCAPEN( CENTRE, ENERGY, SECCEN, CAPRAD, IAT1, IAT2, DAT2,
     &             ATMAX, ATNO, ATXYZ, ATVDW, PI)
      WRITE(*,900) 'CASE2A', CAPRAD, IAT1
C Move ONLY atom 2 - first and last atoms, and the count, are untouched.
      ATXYZ(1,2) = 1.0D0
      ATXYZ(2,2) = 0.0D0
      HFGEN = 3
      CALL HCAPEN( CENTRE, ENERGY, SECCEN, CAPRAD, IAT1, IAT2, DAT2,
     &             ATMAX, ATNO, ATXYZ, ATVDW, PI)
      WRITE(*,900) 'CASE2B', CAPRAD, IAT1

900   FORMAT(A6,' caprad=',F8.3,' iat1=',I7)
      END
