C ********************************************************************
C * tsatr_fast.f - VMDHole patch of tsatr.f.                          *
C *                                                                  *
C * Adds ONE alternate input path: when the coordinate file's name    *
C * ends in .vhb, read a packed binary record VMDHole wrote straight  *
C * from VMD's atom selection, instead of re-parsing 80-column ASCII  *
C * that VMD had just finished formatting. Everything downstream -    *
C * the HTEST/FE name resolution, UCASE, the IGNRES skip, both LMATCH *
C * radius lookups, OATNO - is the stock code, unchanged.             *
C *                                                                  *
C * Coordinates are rounded to 0.001 A on the way in. The PDB path    *
C * quantises them by writing %8.3f and reading 3F8.3, so rounding    *
C * here is what keeps the two paths bit-identical rather than merely *
C * close - HOLE's sphere search is chaotic enough that 1e-12 A of    *
C * extra precision would move the .sph output.                       *
C *                                                                  *
C * Layout (little-endian, native, same host writes and reads it):    *
C *   "VMDHOLEC" 8 bytes | version int4 | natoms int4                 *
C *   natoms x char5   PDB cols 13-17 (HTEST then the 4-char name)    *
C *   natoms x char3   PDB cols 18-20 (residue name)                  *
C *   natoms x char1   PDB col  22    (chain)                         *
C *   natoms x int4    PDB cols 23-26 (residue number, as HOLE's I4)  *
C *   natoms x real8   x, then all y, then all z                      *
C * The identity blocks are copied out of a PDB VMDHole wrote once    *
C * per run, so they are the SAME bytes the ASCII path would have     *
C * produced - only the per-frame coordinates are packed fresh.       *
C ********************************************************************
      SUBROUTINE TSATR(  SIN, NOUT, LERR, LDBUG, LVDW, LBND,
     &	ATMAX, ATNO, ATBRK, ATRES, ATCHN, ATRNO, ATXYZ, ATVDW, ATBND,
     &	MAXLST, BNDNO, BNDBRK, BNDR, VDWNO, VDWBRK, VDWRES, VDWR,
     &  IGNRES, OATNO)
      IMPLICIT NONE
C ********************************************************************
C *                                                                  *
C * This software is an unpublished work containing confidential and *
C * proprietary information of Birkbeck College. Use, disclosure,    *
C * reproduction and transfer of this work without the express       *
C * written consent of Birkbeck College are prohibited. This notice  *
C * must be attached to all copies or extracts of the software.      *
C *                                                                  *
C * (c) 1993 Oliver Smart & Birkbeck College, All rights reserved    *
C * (c) 1996 Oliver Smart & Birkbeck College, All rights reserved    *
C *                                                                  *
C ********************************************************************
C
C Modification history:
C
C Date	Author		Modification
C 12/93	O.S. Smart	Original public release in HOLE suite beta1.0
C 11/95 O.S.S.		Addition of IGNRES string for residue types to be
C			ignored, OATNO to store original atom numbers
C 03/96 O.S.S.		Support for proper hydrogen naming vble HTEST
C
C
C this s/r read pdb atom records from file opened to stream SIN
C and set up vdW and bond radii for each atom (if LVDW and LBND are
C set .true.).
C Adapted to read chain identifier Nov '93

C co-ord file input stream no.
C (return unchanged)
      INTEGER			SIN

C output stream no (to user)
C (return unchanged)
      INTEGER			NOUT

C maximum no. of atoms
C set up as parameter in program (so don't change!)
      INTEGER			ATMAX

C number of atoms read in:
      INTEGER			ATNO

C atom names found in Brookhaven file:
      CHARACTER*4		ATBRK(ATMAX)

C residue name in brook
      CHARACTER*3		ATRES(ATMAX)

C chain identifier in brookhaven pdb file
      CHARACTER*1		ATCHN(ATMAX)

C integer residue no.
      INTEGER			ATRNO(ATMAX)

C co-ordinates
      DOUBLE PRECISION		ATXYZ(3,ATMAX)

C vdw and bond radii of each atoms
      DOUBLE PRECISION		ATVDW(ATMAX), ATBND(ATMAX)

C maximum no of entries in lists
C set up as parameter in program (so don't change!)
      INTEGER			MAXLST

C bond radius list
C (return unchanged)
      INTEGER			BNDNO
      CHARACTER*4		BNDBRK(MAXLST)
      DOUBLE PRECISION		BNDR(MAXLST)

C vdW radius list
C (return unchanged)
      INTEGER			VDWNO
      CHARACTER*4		VDWBRK(MAXLST)
      CHARACTER*3		VDWRES(MAXLST)
      DOUBLE PRECISION		VDWR(MAXLST)


C string listing residue types to be ignored
      CHARACTER*80		IGNRES

C need to record on the initial read of the pdb file the
C original atom numbers of the each atom from the pdb file
C - as we ignore some residues in the read.
C oatno(0) is the total number of original atoms in the pdb file
C oatno(1) is the original atom number of the stored atom#1 etc.
      INTEGER                   OATNO(0:ATMAX)

C logical vairable - produce debug output
C (return unchanged)
      LOGICAL			LDBUG

C if an error found is found in s/r - set lerr true
C and program will stop
C (return unchanged)
      LOGICAL			LERR

C set up bond list if LBND, vdW if LVDW
C (return unchanged)
      LOGICAL			LBND, LVDW

C line to read info
      CHARACTER*132		LINE

C count index
      INTEGER			JCOUNT

C function to match character strings
      LOGICAL			LMATCH

C one character test vble to see whether we have funny
C hydrogen type record
      CHARACTER*1		HTEST


C --- VMDHole binary-coordinate path -------------------------------
C name of the file connected to SIN, tested for the .vhb suffix
      CHARACTER*512             VHFNAM
C true when that suffix is present
      LOGICAL                   VHBIN
C header fields and loop/bookkeeping indices
      INTEGER                   VHVER, VHNAT, VHI, VHJ, VHIOS, VHL
C magic, and the scratch the packed 5-char name field is read into
      CHARACTER*8               VHMAG
      CHARACTER*5               VHNM5

C end of decs ******************

C initialize variables
      ATNO = 0
      OATNO(0) = 0
      LERR = .FALSE.

C === VMDHole: binary coordinate record? ==========================
      VHBIN = .FALSE.
      VHFNAM = ' '
      INQUIRE( SIN, NAME = VHFNAM)
      VHL = LEN(VHFNAM)
10001 CONTINUE
      IF ((VHL.GT.1) .AND. (VHFNAM(VHL:VHL).EQ.' ')) THEN
        VHL = VHL - 1
        GOTO 10001
      ENDIF
      IF (VHL.GE.4) THEN
        IF (VHFNAM(VHL-3:VHL).EQ.'.vhb') VHBIN = .TRUE.
      ENDIF
      IF (.NOT.VHBIN) GOTO 10        

C reopen the same file for raw byte access - the caller opened it
C formatted, and INQUIRE(NAME) still reports it afterwards, which is
C all hole.f does with the unit before it closes it.
      CLOSE( SIN)
      OPEN( UNIT= SIN, FILE= VHFNAM(1:VHL), ACCESS= 'STREAM',
     &      FORM= 'UNFORMATTED', STATUS= 'OLD', IOSTAT= VHIOS)
      IF (VHIOS.NE.0) THEN
        WRITE(NOUT,*) '***ERROR***', CHAR(7)
        WRITE(NOUT,*) 'Cannot reopen binary co-ord file: ',
     &                VHFNAM(1:VHL)
        LERR = .TRUE.
        GOTO 55555
      ENDIF
      READ( SIN, IOSTAT= VHIOS) VHMAG
      IF ((VHIOS.NE.0) .OR. (VHMAG.NE.'VMDHOLEC')) THEN
        WRITE(NOUT,*) '***ERROR***', CHAR(7)
        WRITE(NOUT,*) 'Not a VMDHole binary co-ord file: ',
     &                VHFNAM(1:VHL)
        LERR = .TRUE.
        GOTO 55555
      ENDIF
      READ( SIN) VHVER
      READ( SIN) VHNAT
      IF (VHVER.NE.1) THEN
        WRITE(NOUT,*) '***ERROR***', CHAR(7)
        WRITE(NOUT,*) 'Unsupported binary co-ord version ', VHVER
        LERR = .TRUE.
        GOTO 55555
      ENDIF
      IF (VHNAT.GT.ATMAX) THEN
        WRITE(NOUT,*) '***ERROR***', CHAR(7)
        WRITE(NOUT,*) 'Have exeeded array bound ATMAX ', ATMAX
        WRITE(NOUT,*) '(the maximum no. of atoms allowed)'
        LERR = .TRUE.
        GOTO 55555
      ENDIF

C ---- slurp the blocks positionally, atom i into slot i ----------
C HTEST resolution is done here, inline, exactly as the ASCII branch
C does it right after its READ.
      DO 10010 VHI = 1, VHNAT
        READ( SIN) VHNM5
        ATBRK(VHI) = VHNM5(2:5)
        IF (VHNM5(1:1).EQ.'H') ATBRK(VHI) = VHNM5(1:1)//ATBRK(VHI)
        IF ((ATBRK(VHI).EQ.'E  ').AND.(VHNM5(1:1).EQ.'F'))
     &                                      ATBRK(VHI) = 'FE  '
10010 CONTINUE
      DO 10020 VHI = 1, VHNAT
        READ( SIN) ATRES(VHI)
10020 CONTINUE
      DO 10030 VHI = 1, VHNAT
        READ( SIN) ATCHN(VHI)
10030 CONTINUE
      DO 10040 VHI = 1, VHNAT
        READ( SIN) ATRNO(VHI)
10040 CONTINUE
      DO 10060 VHJ = 1, 3
        DO 10050 VHI = 1, VHNAT
          READ( SIN) ATXYZ(VHJ,VHI)
C match the PDB path's %8.3f / 3F8.3 round trip
          ATXYZ(VHJ,VHI) = DNINT( ATXYZ(VHJ,VHI)*1000.0D0)/1000.0D0
10050   CONTINUE
10060 CONTINUE

C ---- compact in place, applying IGNRES + the radius lookups -----
      OATNO(0) = VHNAT
      ATNO = 0
      DO 10100 VHI = 1, VHNAT
        CALL UCASE(ATRES(VHI))
        IF (INDEX(IGNRES, ATRES(VHI)).NE.0) GOTO 10100
        ATNO = ATNO + 1
        ATBRK(ATNO) = ATBRK(VHI)
        ATRES(ATNO) = ATRES(VHI)
        ATCHN(ATNO) = ATCHN(VHI)
        ATRNO(ATNO) = ATRNO(VHI)
        ATXYZ(1,ATNO) = ATXYZ(1,VHI)
        ATXYZ(2,ATNO) = ATXYZ(2,VHI)
        ATXYZ(3,ATNO) = ATXYZ(3,VHI)
        IF (LBND) THEN
          DO 10070 JCOUNT = 1, BNDNO
            IF (LMATCH( '?', 4, BNDBRK(JCOUNT), ATBRK(ATNO))) THEN
              ATBND(ATNO) = BNDR(JCOUNT)
              GOTO 10071
            ENDIF
10070     CONTINUE
          WRITE(NOUT,*) '***ERROR***', CHAR(7)
          WRITE(NOUT,*) 'Cannot find bond radius for atom:'
          WRITE(NOUT,*) ATBRK(ATNO), ATRES(ATNO), ATRNO(ATNO)
          LERR = .TRUE.
          GOTO 55555
        ENDIF
10071   CONTINUE
        IF (LVDW) THEN
          DO 10080 JCOUNT = 1, VDWNO
            IF ( LMATCH( '?', 4, VDWBRK(JCOUNT), ATBRK(ATNO)) .AND.
     &           LMATCH( '?', 3, VDWRES(JCOUNT), ATRES(ATNO))  ) THEN
              ATVDW(ATNO) = VDWR(JCOUNT)
              GOTO 10081
            ENDIF
10080     CONTINUE
          WRITE(NOUT,*) '***ERROR***', CHAR(7)
          WRITE(NOUT,*) 'Cannot find vdW radius for atom:'
          WRITE(NOUT,*) ATBRK(ATNO), ATRES(ATNO), ATRNO(ATNO)
          LERR = .TRUE.
          GOTO 55555
        ENDIF
10081   CONTINUE
        IF (LDBUG) WRITE(NOUT,*) 'debug Atom: ',
     &    ATBRK(ATNO), ATRES(ATNO), ATRNO(ATNO),
     &    (ATXYZ(JCOUNT,ATNO), JCOUNT = 1, 3)
        IF (LDBUG) WRITE(NOUT,*) '           ',
     &    ' number ', ATNO, ATBND(ATNO), ATVDW(ATNO)
        OATNO(ATNO) = VHI
10100 CONTINUE
      GOTO 55555
C === end VMDHole binary path =====================================



C read from file until 'atom' record is found
C (do until loop to 55555 - the return statement)
10    CONTINUE
        READ( SIN, '(A)', END= 55555) LINE
C October 1993 read HETATM's as well as atoms
	IF ( (LINE(1:4).NE.'ATOM') .AND. 
     &       (LINE(1:6).NE.'HETATM') ) GOTO 10

C 'ATOM' found ___read in line
C one more atom
	ATNO = ATNO + 1
C increment original number of atoms
        OATNO(0) = OATNO(0) + 1

C make sure that we have not exceeded maximum 
C no of atoms which can be read in
	IF (ATNO.GT.ATMAX) THEN
	  WRITE(NOUT,*) '***ERROR***', CHAR(7)
	  WRITE(NOUT,*) 'Have exeeded array bound ATMAX ', ATMAX
	  WRITE(NOUT,*) '(the maximum no. of atoms allowed)'
	  LERR = .TRUE.
	  GOTO 55555
	ENDIF

C read line in pdb format
C reading four character brookhaven atomname but if 
C the standard can put an "H" preceeding this -
C read into HTEST
	READ(LINE(5:80), '(8X,A1,A4,A3,1X,A1,I4,4X,3F8.3)')
     &    HTEST,
     &	  ATBRK(ATNO), ATRES(ATNO), ATCHN(ATNO), ATRNO(ATNO),
     &	  (ATXYZ(JCOUNT,ATNO), JCOUNT = 1, 3)
        IF (HTEST.EQ.'H') ATBRK(ATNO) = HTEST//ATBRK(ATNO)
C 22/5/98 support for FE records of heme
        IF ((ATBRK(ATNO).EQ.'E  ').AND.(HTEST.EQ.'F'))
     &                                      ATBRK(ATNO) = 'FE  '

C should we ignore this residue?
C Check to see whether residue type appears in ignres
C comparison should be between upper case strings
        CALL UCASE(ATRES(ATNO))
        IF (INDEX(IGNRES, ATRES(ATNO)).NE.0) THEN
          ATNO = ATNO - 1
C (leave oatno(0) alone)
C read next atom
          GOTO 10
        ENDIF

C setup bond radius
	IF (LBND) THEN
C go thru' bond radius lists
	  DO 20 JCOUNT = 1, BNDNO
C use function lmatch to match character strings
C 4 is the no. of characters to be matched and ? is the wildcard
C character
	    IF (LMATCH( '?', 4, BNDBRK(JCOUNT), ATBRK(ATNO))) THEN
	       ATBND(ATNO) = BNDR(JCOUNT)
	       GOTO 201
	    ENDIF
20	  CONTINUE
C error cannot find bond radius for atom
	  WRITE(NOUT,*) '***ERROR***', CHAR(7)
	  WRITE(NOUT,*) 'Cannot find bond radius for atom:'
	  WRITE(NOUT,*) ATBRK(ATNO), ATRES(ATNO), ATRNO(ATNO)
	  LERR = .TRUE.
	  GOTO 55555

	ENDIF

C have found bond radius for atom
201	CONTINUE

C setup van der Waals radius
	IF (LVDW) THEN
C go thru' bond radius lists
	  DO 30 JCOUNT = 1, VDWNO
C use function lmatch to match character strings
C 4 is the no. of characters to be matched and ? is the wildcard
C character HERE we must match 3*C residue name as well.
	    IF ( LMATCH( '?', 4, VDWBRK(JCOUNT), ATBRK(ATNO)) .AND.
     &		 LMATCH( '?', 3, VDWRES(JCOUNT), ATRES(ATNO))  ) THEN
	       ATVDW(ATNO) = VDWR(JCOUNT)
	       GOTO 301
	    ENDIF
30	  CONTINUE
C error cannot find vdW radius for atom
	  WRITE(NOUT,*) '***ERROR***', CHAR(7)
	  WRITE(NOUT,*) 'Cannot find vdW radius for atom:'
	  WRITE(NOUT,*) ATBRK(ATNO), ATRES(ATNO), ATRNO(ATNO)
	  LERR = .TRUE.
	  GOTO 55555

	ENDIF

C have found bond radius for atom
301	CONTINUE

C debug output
	IF (LDBUG) WRITE(NOUT,*) 'debug Atom: ',
     &	  ATBRK(ATNO), ATRES(ATNO), ATRNO(ATNO),
     &	  (ATXYZ(JCOUNT,ATNO), JCOUNT = 1, 3)
	IF (LDBUG) WRITE(NOUT,*) '           ',
     &	  ' number ', ATNO, ATBND(ATNO), ATVDW(ATNO)

C record original atom number
        OATNO(ATNO) = OATNO(0)

C read next entry:
      GOTO 10

55555 RETURN
      END
