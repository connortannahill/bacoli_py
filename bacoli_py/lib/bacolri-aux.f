C bacolri-aux.f (supporting routines for BACOLRI, the source code
C the source code contained in the file bacolri.f).

C This file contains:

c (i) The LAMPACK software package: P. Keast. FORTRAN package for
c solving certain almost block diagonal matrices.

C (ii) From the book, Hairer, E. and Wanner, G., Solving Ordinary
c Equations II. Stiff and Differential-Algebraic Problems. 
c TODO: Finish this citation.
c            radau5, radcor

c (iii) From the paper, de Boor, C., Package for Calculating with
c B-splines, c SIAM J. Numer. Anal., vol. 14, no. 3, June 1977,
c pp. 441-472: 
c            bsplvd, bsplvn, interv.
c However, the versions of these routines appearing in BACOLI 
c have been modified or rewritten to some extent. 
c Similar versions of BSPLVD and INTERV are available
c at www.netlib.org/pppack and a similar version of BSPLVN is
c available at www.netlib.org/slatec. As mentioned above, software
c from SLATEC is in the public domain. It's disclaimer may be found 
c at: www.netlib.org/slatec/src/aaaaaa.f. There does not appear to
c be a license associated with the PPPACK collection.

c (iv) From the EISPACK (www.netlib.org/eispack/) numerical software
c collection: 
c            imtql1, imtql2, pythag. 
c There does not appear to be a license associated with the 
c EISPACK collection.

c (v) From Patrick Keast: 
c            gauleg. 
c Written by Pat Keast, Dalhousie University. No licensing 
c information available.

c (vi) From the LINPACK (www.netlib.org/linpack/) numerical software
c collection: 
c            decomc, decomr, estrad, slvrad.
c We have modified this source code slightly: declarations like
c array(1) were changed to array(*). There does not appear to
c be a license associated with the LINPACK collection.

c (vii) From the BLAS (www.netlib.org/blas/) numerical software
c collection: 
c            daxpy, dcopy, dscal.
c These routines are subject to the following legal restrictions:
c see www.netlib.org/blas/faq.html#2

        SUBROUTINE CCRCMPri(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
        implicit none
C
C***************************************************************
C
C  C C R C M PDECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - COMPLEX(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - COMPLEX(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - COMPLEX(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C               PIVOT  - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C               PIVOT  - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 1, 2003.
c
c-----------------------------------------------------------------------
C***************************************************************
C
c-----------------------------------------------------------------------
        integer n,nrwtop,novrlp,nrwblk,nclblk,nbloks,nrwbot,iflag,
     &          nrwtp1,nrowel,nrwel1,nvrlp0,i,iplus1,ipvt,j,l,incr,k,
     &          kplus1,jplus1,jminn,loop,incrj,iplusn,incrn,irwblk,
     &          ipvblk,jrwblk
c-----------------------------------------------------------------------
        DOUBLE COMPLEX TOPBLK,ARRAY,BOTBLK
        DOUBLE COMPLEX ROWPIV,ROWMLT,COLPIV,SWAP,COLMLT
        DOUBLE PRECISION ROWMAX,COLMAX,TEMPIV,ZERO,PIVMAX
        INTEGER PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*)
        DATA ZERO/0.0/
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        IFLAG = 0
        PIVMAX = ZERO
        NRWTP1 = NRWTOP+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
        IF(N.NE.NBLOKS*NRWBLK+NOVRLP)GO TO 10
        IF(NOVRLP.NE.NRWTOP+NRWBOT)GO TO 10
        IF(NCLBLK.NE.NOVRLP+NRWBLK)GO TO 10
        IF(NOVRLP.GT.NRWBLK)GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        GO TO 100
10      CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IFLAG = 1
        RETURN
100     CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
        DO 190 I = 1,NRWTOP
           IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IPVT = I
           COLMAX = ABS(TOPBLK(I,I))
           DO 110 J = IPLUS1,NOVRLP
              TEMPIV = ABS(TOPBLK(I,J))
              IF(TEMPIV.LE.COLMAX)GO TO 110
                 IPVT = J
                 COLMAX = TEMPIV
110        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
           PIVMAX = MAX(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           PIVOT(I) = IPVT
           IF(IPVT.EQ.I)GO TO 140
              DO 120 L = I,NRWTOP
                 SWAP = TOPBLK(L,IPVT)
                 TOPBLK(L,IPVT) = TOPBLK(L,I)
                 TOPBLK(L,I) = SWAP
120           CONTINUE
              DO 130 L = 1,NRWBLK
                 SWAP = ARRAY(L,IPVT,1)
                 ARRAY(L,IPVT,1) = ARRAY(L,I,1)
                 ARRAY(L,I,1) = SWAP
130           CONTINUE
140        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           COLPIV = TOPBLK(I,I)
           DO 180 J = IPLUS1,NOVRLP
              COLMLT = TOPBLK(I,J)/COLPIV
              TOPBLK(I,J) = COLMLT
              IF(IPLUS1.GT.NRWTOP)GO TO 160
                 DO 150 L = IPLUS1,NRWTOP
                    TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
150              CONTINUE
160           CONTINUE
              DO 170 L = 1,NRWBLK
                 ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
170           CONTINUE
180        CONTINUE
190     CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
        INCR = 0
        DO 395 K = 1,NBLOKS
           KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
           DO 270 J = NRWTP1,NRWBLK
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = ABS(ARRAY(JMINN,J,K))
              LOOP = JMINN+1
              DO 210 I = LOOP,NRWBLK
                 TEMPIV = ABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.ROWMAX)GO TO 210
                 IPVT = I
                 ROWMAX = TEMPIV
210           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO  1000
              PIVMAX = MAX(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 230
                 DO 220 L = J,NCLBLK
                    SWAP = ARRAY(IPVT,L,K)
                    ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
                    ARRAY(JMINN,L,K) = SWAP
220              CONTINUE
230           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = ARRAY(JMINN,J,K)
              DO 240 I = LOOP,NRWBLK
                 ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
240           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 260 L = JPLUS1,NCLBLK
                 ROWMLT = ARRAY(JMINN,L,K)
                 DO 250 I = LOOP,NRWBLK
                    ARRAY(I,L,K) = ARRAY(I,L,K)
     *                                -ROWMLT*ARRAY(I,J,K)
250              CONTINUE
260           CONTINUE
270        CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
           DO 390 I = NRWEL1,NRWBLK
              IPLUSN = I+NRWTOP
              IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = IPLUSN
              COLMAX = ABS(ARRAY(I,IPVT,K))
              LOOP = IPLUSN+1
              DO 310 J = LOOP,NCLBLK
                 TEMPIV = ABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.COLMAX)GO TO 310
                 IPVT = J
                 COLMAX = TEMPIV
310           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = MAX(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRN = INCR+IPLUSN
              PIVOT(INCRN) = INCR+IPVT
              IRWBLK = IPLUSN-NRWBLK
              IF(IPVT.EQ.IPLUSN)GO TO 340
                 DO 315 L = I,NRWBLK
                    SWAP = ARRAY(L,IPVT,K)
                    ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
                    ARRAY(L,IPLUSN,K) = SWAP
315              CONTINUE
                 IPVBLK = IPVT-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 330
                    DO 320 L = 1,NRWBLK
                       SWAP = ARRAY(L,IPVBLK,KPLUS1)
                       ARRAY(L,IPVBLK,KPLUS1)
     *                                 = ARRAY(L,IRWBLK,KPLUS1)
                       ARRAY(L,IRWBLK,KPLUS1) = SWAP
320                 CONTINUE
                    GO TO 340
330              CONTINUE
                 DO 335 L = 1,NRWBOT
                    SWAP = BOTBLK(L,IPVBLK)
                    BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
                    BOTBLK(L,IRWBLK) = SWAP
335              CONTINUE
340           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              COLPIV = ARRAY(I,IPLUSN,K)
              DO 380 J = LOOP,NCLBLK
                 COLMLT = ARRAY(I,J,K)/COLPIV
                 ARRAY(I,J,K) = COLMLT
                 IF(I.EQ.NRWBLK)GO TO 350
                    DO 345 L = IPLUS1,NRWBLK
                       ARRAY(L,J,K) = ARRAY(L,J,K)
     *                                -COLMLT*ARRAY(L,IPLUSN,K)
345                 CONTINUE
350              CONTINUE
                 JRWBLK = J-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 370
                    DO 360 L = 1,NRWBLK
                       ARRAY(L,JRWBLK,KPLUS1) =
     *                                  ARRAY(L,JRWBLK,KPLUS1)
     *                         -COLMLT*ARRAY(L,IRWBLK,KPLUS1)
360                 CONTINUE
                    GO TO 380
370              CONTINUE
                 DO 375 L = 1,NRWBOT
                    BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)
     *                              -COLMLT*BOTBLK(L,IRWBLK)
375              CONTINUE
380           CONTINUE
390        CONTINUE
           INCR = INCR + NRWBLK
395     CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 500
           DO 470 J = NRWTP1,NVRLP0
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = ABS(BOTBLK(JMINN,J))
              LOOP = JMINN+1
              DO 410 I = LOOP,NRWBOT
                 TEMPIV = ABS(BOTBLK(I,J))
                 IF(TEMPIV.LE.ROWMAX) GO TO 410
                 IPVT = I
                 ROWMAX = TEMPIV
410           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = MAX(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 430
                 DO 420 L = J,NOVRLP
                    SWAP = BOTBLK(IPVT,L)
                    BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
                    BOTBLK(JMINN,L) = SWAP
420              CONTINUE
430           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = BOTBLK(JMINN,J)
              DO 440 I = LOOP,NRWBOT
                 BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
440           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 460 L = JPLUS1,NOVRLP
                 ROWMLT = BOTBLK(JMINN,L)
                 DO 450 I = LOOP,NRWBOT
                    BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
450              CONTINUE
460           CONTINUE
470        CONTINUE
500     CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
        IF(PIVMAX+ABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
1000    CONTINUE
        IFLAG = -1
        RETURN
        END

        SUBROUTINE CCRSLVri(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B)
        implicit none
C
C***************************************************************
C
C  C C R S L V  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - COMPLEX(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - COMPLEX(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - COMPLEX(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - COMPLEX(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - COMPLEX(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C                    X - COMPLEX(N)
C                         THE SOLUTION VECTOR
C
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 1, 2003.
c
c-----------------------------------------------------------------------
C***************************************************************
C
c-----------------------------------------------------------------------
        integer nrwtop,novrlp,nrwblk,nclblk,nbloks,nrwbot,nrwtp1,nrwbk1,
     &          nvrlp1,nrwtp0,nrwbt1,nrowel,nrwel1,nvrlp0,nblks1,nbktop,
     &          j,loop,i,incr,k,incrtp,incrj,incri,jpivot,jrwtop,ll,
     &          nrwbtl,l,l1,iplusn,incrn,ipvtn,nrwell,ipvti
c-----------------------------------------------------------------------
        DOUBLE COMPLEX TOPBLK,ARRAY,BOTBLK,B
        DOUBLE COMPLEX DOTPRD,XJ,XINCRJ,BINCRJ,SWAP
        INTEGER PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           B(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              XJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*XJ
110           CONTINUE
120        CONTINUE
130     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
250              CONTINUE
260           CONTINUE
270        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -B(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ)
440           CONTINUE
              B(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = B(INCRN)
                 B(INCRN) = B(IPVTN)
                 B(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -B(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J)
510        CONTINUE
           B(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = B(I)
                 B(I) = B(IPVTI)
                 B(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
        RETURN
        END

c=======================================================================
c     RADAU5 DAE SOLVER
c=======================================================================
      SUBROUTINE RADAU5(N,FCN,X,Y,XEND,H,RTOL,ATOL,ITOL,JAC,MAS,SOLOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,CWORK,IDID,
     &                  f, fvec, derivf, bndxa, difbxa, bndxb, difbxb,
     &                  uinit, vec)
      implicit none
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
C     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
C                     M*Y'=F(X,Y).
C     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA)
C     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT.
C     CF. SECTION IV.8
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND
C              E-MAIL:  Ernst.Hairer@math.unige.ch
C                       Gerhard.Wanner@math.unige.ch
C
C     THIS CODE IS PART OF THE BOOK:
C         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
C         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
C
C     VERSION OF JULY 9, 1996
C        (small correction April 14, 2000)
C
C     INPUT PARAMETERS
C     ----------------
C     N           DIMENSION OF THE SYSTEM
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),F(N)
C                    F(1)=...   ETC.
C                 RPAR, IPAR (SEE BELOW)
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
C                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y.
C                 THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(N,X,Y,DFY,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),DFY(*)
C                    DFY(1,1)= ...
C
C     ----   MAS HAS ANALOG MEANINGS      -----
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
C
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
C                 MATRIX M.
c-----------------------------------------------------------------------
c                 SUBROUTINE MAS(AM,RPAR,IPAR)
c                 the mass-matrix AM is stored in the ABD form.
c-----------------------------------------------------------------------
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION.
C                 IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (N,X,Y,ITOL,ATOL,RTOL,
C                                       RPAR,IPAR,IRTRN)
C                    DOUBLE PRECISION X,Y(N),CONT(LRC)
C                    ....
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM.
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE OF THE CODE
C                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE
C                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
C                 WORK(21),..,WORK(LWORK) SERVE AS WORKING SPACE
C                 FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST
C                            3*lenpd+12*N+20
C
C     LWORK       DECLARED LENGTH OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..,
C                 IWORK(20) TO ZERO BEFORE CALLING.
C                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA.
C                 "LIWORK" MUST BE AT LEAST 3*N+20.
C
C     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.
C
C ----------------------------------------------------------------------
C
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
C              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
C              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
C              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
C              IT DOES NOT WORK FOR BANDED JACOBIAN,
C              AND NOT FOR IMPLICIT SYSTEMS.
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP.
C              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7.
C
C    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION
C              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD.
C              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED.
C              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS
C              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN
C              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.).
C              DEFAULT IS IWORK(4)=0.
C
C       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
C       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
C       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
C       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER.
C       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
C       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
C
C    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR
C              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
C              DEFAULT IWORK(5)=N.
C
C    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0.
C
C    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0.
C
C    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY
C              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
C              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL
C              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1.
C              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS;
C              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES
C              OFTEN SLIGHTLY FASTER RUNS
C
C ----------
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
C              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
C              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER
C              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO
C              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.
C              DEFAULT 0.001D0.
C
C    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
C              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER.
C              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0)
C
C    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE
C              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A
C              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR
C              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE
C              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS
C              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD.
C              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 .
C
C    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WORK(8) <= HNEW/HOLD <= WORK(9)
C              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0
C
C-----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS
C     -----------------
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND).
C
C     Y(N)        NUMERICAL SOLUTION AT X
C
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
C   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
C   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                      OR NUMERICALLY)
C   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
C   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
C   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
C                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
C   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES
C   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH
C                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS,
C                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED
C-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 1, 2003.
c
c-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS
C *** *** *** *** *** *** *** *** *** *** *** *** ***
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,ITOL,LWORK,LIWORK,IDID,IWORK(LIWORK),IPAR(*)

        external                f
        external                fvec
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
        external                uinit

      DOUBLE PRECISION X,XEND,H,Y(N),WORK(LWORK)
      DOUBLE PRECISION ATOL(*),RTOL(*),RPAR(*)
      DOUBLE COMPLEX CWORK(*)
      logical vec
c-----------------------------------------------------------------------
c local variables
      INTEGER NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,NMAX,NIT,
     &        IEZ1,IEZ2,IEZ3,IEY0,IESCAL,IEF1,IEF2,IEF3,IECON,IEJAC,
     &        IEMAS,IEE1,IEE2R,ISTORE,IEIP1,IEIP2,IECZ2
      DOUBLE PRECISION UROUND,EXPM,QUOT,SAFE,THET,TOLST,FNEWT,QUOT1,
     &                 QUOT2,HMAX,FACL,FACR
c-----------------------------------------------------------------------
c loop indices
      INTEGER I
c-----------------------------------------------------------------------

      LOGICAL ARRET,STARTN,PRED
      EXTERNAL FCN,JAC,MAS,SOLOUT
c-----------------------------------------------------------------------
      INTEGER NCONTI,NPDE,NINT,KCOL,NSZJAC
      PARAMETER (NCONTI = 2)
      NPDE = IPAR(1)
      KCOL = IPAR(2)
      NINT = IPAR(3)
      NSZJAC = NPDE * NPDE * (NCONTI + NINT * KCOL * (KCOL + NCONTI)
     &       + NCONTI)
c-----------------------------------------------------------------------
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS
C *** *** *** *** *** *** ***
       NFCN=0
       NJAC=0
       NSTEP=0
       NACCPT=0
       NREJCT=0
       NDEC=0
       NSOL=0
       ARRET=.FALSE.
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0
      IF (WORK(1).EQ.0.0D0) THEN
         UROUND=1.0D-16
      ELSE
         UROUND=WORK(1)
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- CHECK AND CHANGE THE TOLERANCES
      EXPM=2.0D0/3.0D0
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          ELSE
              QUOT=ATOL(1)/RTOL(1)
              RTOL(1)=0.1D0*RTOL(1)**EXPM
              ATOL(1)=RTOL(1)*QUOT
          END IF
      ELSE
          DO I=1,NPDE
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          ELSE
              QUOT=ATOL(I)/RTOL(I)
              RTOL(I)=0.1D0*RTOL(I)**EXPM
              ATOL(I)=RTOL(I)*QUOT
          END IF
          END DO
      END IF
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF (IWORK(2).EQ.0) THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF (NMAX.LE.0) THEN
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS
      IF (IWORK(3).EQ.0) THEN
         NIT=7
      ELSE
         NIT=IWORK(3)
         IF (NIT.LE.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS
      IF(IWORK(4).EQ.0)THEN
         STARTN=.FALSE.
      ELSE
         STARTN=.TRUE.
      END IF
C -------- PRED   STEP SIZE CONTROL
      IF(IWORK(8).LE.1)THEN
         PRED=.TRUE.
      ELSE
         PRED=.FALSE.
      END IF
C --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION
      IF (WORK(2).EQ.0.0D0) THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
      IF (WORK(3).EQ.0.D0) THEN
         THET=0.001D0
      ELSE
         THET=WORK(3)
         IF (THET.GE.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(3)=',WORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
      TOLST=RTOL(1)
      IF (WORK(4).EQ.0.D0) THEN
         FNEWT=MAX(10*UROUND/TOLST,MIN(0.03D0,TOLST**0.5D0))
      ELSE
         FNEWT=WORK(4)
         IF (FNEWT.LE.UROUND/TOLST) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(4)=',WORK(4)
            ARRET=.TRUE.
         END IF
      END IF
C --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST.
      IF (WORK(5).EQ.0.D0) THEN
         QUOT1=1.D0
      ELSE
         QUOT1=WORK(5)
      END IF
      IF (WORK(6).EQ.0.D0) THEN
         QUOT2=1.2D0
      ELSE
         QUOT2=WORK(6)
      END IF
      IF (QUOT1.GT.1.0D0.OR.QUOT2.LT.1.0D0) THEN
         WRITE(6,*)' CURIOUS INPUT FOR WORK(5,6)=',QUOT1,QUOT2
         ARRET=.TRUE.
      END IF
C -------- MAXIMAL STEP SIZE
      IF (WORK(7).EQ.0.D0) THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(7)
      END IF
C -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(8).EQ.0.D0)THEN
         FACL=5.D0
      ELSE
         FACL=1.D0/WORK(8)
      END IF
      IF(WORK(9).EQ.0.D0)THEN
         FACR=1.D0/8.0D0
      ELSE
         FACR=1.D0/WORK(9)
      END IF
      IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT WORK(8,9)=',WORK(8),WORK(9)
            ARRET=.TRUE.
         END IF
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEZ1=21
      IEZ2=IEZ1+N
      IEZ3=IEZ2+N
      IEY0=IEZ3+N
      IESCAL=IEY0+N
      IEF1=IESCAL+N
      IEF2=IEF1+N
      IEF3=IEF2+N
      IECON=IEF3+N
      IEJAC=IECON+4*N
      IEMAS=IEJAC+NSZJAC
      IEE1=IEMAS+NSZJAC
c-----------------------------------------------------------------------
c ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN CWORK -----
      IEE2R=1
      IECZ2=IEE2R+NSZJAC
c-----------------------------------------------------------------------
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEE1+NSZJAC-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP1=21
      IEIP2=IEIP1+N
C --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEIP2+N-1
      IF (ISTORE.GT.LIWORK) THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
c      write(*,*) 'calling radcor'
      CALL RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,JAC,MAS,SOLOUT,
     &   IDID,NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,STARTN,
     &   PRED,FACL,FACR, WORK(IEZ1),WORK(IEZ2),WORK(IEZ3),WORK(IEY0),
     &   WORK(IESCAL),WORK(IEF1),WORK(IEF2),WORK(IEF3),WORK(IEJAC),
     &   WORK(IEE1),CWORK(IEE2R),WORK(IEMAS),IWORK(IEIP1),IWORK(IEIP2),
     &   WORK(IECON),NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR,
     &   CWORK(IECZ2),f, fvec, derivf, bndxa, difbxa, bndxb, difbxb,
     &   uinit, vec)
c      write(*,*) 'leaving radcor'

      IWORK(14)=NFCN
      IWORK(15)=NJAC
      IWORK(16)=NSTEP
      IWORK(17)=NACCPT
      IWORK(18)=NREJCT
      IWORK(19)=NDEC
      IWORK(20)=NSOL
C -------- RESTORE TOLERANCES
      EXPM=1.0D0/EXPM
      IF (ITOL.EQ.0) THEN
              QUOT=ATOL(1)/RTOL(1)
              RTOL(1)=(10.0D0*RTOL(1))**EXPM
              ATOL(1)=RTOL(1)*QUOT
      ELSE
          DO I=1,NPDE
              QUOT=ATOL(I)/RTOL(I)
              RTOL(I)=(10.0D0*RTOL(I))**EXPM
              ATOL(I)=RTOL(I)*QUOT
          END DO
      END IF
C ----------- RETURN -----------
      RETURN
      END
C ***********************************************************
C
      SUBROUTINE RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,MAS,SOLOUT,IDID,NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,
     &   QUOT2,NIT,STARTN,PRED,FACL,FACR,Z1,Z2,Z3,Y0,SCAL,F1,F2,F3,
     &   FJAC,E1,E2R,FMAS,IP1,IP2,CONT,NFCN,NJAC,NSTEP,NACCPT,
     &   NREJCT,NDEC,NSOL,RPAR,IPAR,CZ2,f, fvec, derivf, bndxa, difbxa, 
     &   bndxb, difbxb, uinit, vec)
      implicit none
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR RADAU5
C     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED
C ----------------------------------------------------------
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 21, 2003.
c
c-----------------------------------------------------------------------
C ----------------------------------------------------------
C         DECLARATIONS
C ----------------------------------------------------------
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,ITOL,IDID,NMAX,NIT,NFCN,NJAC,NSTEP,NACCPT,NREJCT,
     &        NDEC,NSOL,IP1(N),IP2(N),IPAR(*)
      DOUBLE PRECISION X,XEND,HMAX,H,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,
     &                 FACL,FACR,Y(N),Z1(N),Z2(N),Z3(N),Y0(N),SCAL(N),
     &                 F1(N),F2(N),F3(N),CONT(4*N)
c-----------------------------------------------------------------------
c the dimensions of the following arrays have been changed.
      DOUBLE PRECISION E1(*), FJAC(*), FMAS(*)
      DOUBLE COMPLEX E2R(*), CZ2(N)
c-----------------------------------------------------------------------
      DOUBLE PRECISION RTOL(*),ATOL(*),RPAR(*)
      LOGICAL STARTN,PRED
      EXTERNAL FCN,JAC,MAS,SOLOUT
      external f, fvec, derivf, bndxa, difbxa, bndxb, difbxb, uinit
c-----------------------------------------------------------------------
c variables in common expression
      INTEGER NN,NN2,NN3,NN4
      DOUBLE PRECISION XSOL,HSOL,C2M1,C1M1
      COMMON /CONRA5/NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1
c-----------------------------------------------------------------------
c local variables
      INTEGER NPDE,NINT,KCOL
c-----------------------------------------------------------------------
      INTEGER NSING,IRTRN,N2,N3,IER,NEWT
      DOUBLE PRECISION SQ6,C1,C2,C1MC2,DD1,DD2,DD3,U1,ALPH,BETA,CNO,T11,
     &                 T12,T13,T21,T22,T23,T31,TI11,TI12,TI13,TI21,TI22,
     &                 TI23,TI31,TI32,TI33,POSNEG,HMAXN,HOLD,FACCON,
     &                 CFAC,HHFAC,FAC1,ALPHN,
     &                 BETAN,XPH,C1Q,C2Q,C3Q,AK1,AK2,AK3,Z1I,Z2I,Z3I,
     &                 THETA,A1,A2,A3,DYNO,DENOM,THQ,DYNOLD,THQOLD,
     &                 DYTH,QNEWT,F1I,F2I,F3I,ERR,FAC,QUOT,HNEW,FACGUS,
     &                 HACC,ERRACC,AK,ACONT3,HOPT,QT
      LOGICAL REJECT,FIRST,CALJAC
      LOGICAL LAST, vec
c-----------------------------------------------------------------------
c loop indices
      INTEGER I,J
c-----------------------------------------------------------------------
      NPDE = IPAR(1)
      KCOL = IPAR(2)
      NINT = IPAR(3)
c-----------------------------------------------------------------------
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
C --------- DUPLIFY N FOR COMMON BLOCK CONT -----
      NN=N
      NN2=2*N
      NN3=3*N
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
c-----------------------------------------------------------------------
c      write(*,*) 'calling mas'
      CALL MAS(FMAS,RPAR,IPAR)
c      write(*,*) 'leaving mas'
c-----------------------------------------------------------------------
C ---------- CONSTANTS ---------
      SQ6=SQRT(6.D0)
      C1=(4.D0-SQ6)/10.D0
      C2=(4.D0+SQ6)/10.D0
      C1M1=C1-1.D0
      C2M1=C2-1.D0
      C1MC2=C1-C2
      DD1=-(13.D0+7.D0*SQ6)/3.D0
      DD2=(-13.D0+7.D0*SQ6)/3.D0
      DD3=-1.D0/3.D0
      U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0
      ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0
      BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*SQRT(3.D0)/60.D0
      CNO=ALPH**2+BETA**2
      U1=1.0D0/U1
      ALPH=ALPH/CNO
      BETA=BETA/CNO
      T11=9.1232394870892942792D-02
      T12=-0.14125529502095420843D0
      T13=-3.0029194105147424492D-02
      T21=0.24171793270710701896D0
      T22=0.20412935229379993199D0
      T23=0.38294211275726193779D0
      T31=0.96604818261509293619D0
      TI11=4.3255798900631553510D0
      TI12=0.33919925181580986954D0
      TI13=0.54177053993587487119D0
      TI21=-4.1787185915519047273D0
      TI22=-0.32768282076106238708D0
      TI23=0.47662355450055045196D0
      TI31=-0.50287263494578687595D0
      TI32=2.5719269498556054292D0
      TI33=-0.59603920482822492497D0
      POSNEG=SIGN(1.D0,XEND-X)
      HMAXN=MIN(ABS(HMAX),ABS(XEND-X))
      IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
      H=MIN(ABS(H),HMAXN)
      H=SIGN(H,POSNEG)
      HOLD=H
      REJECT=.FALSE.
      FIRST=.TRUE.
      LAST=.FALSE.
      IF ((X+H*1.0001D0-XEND)*POSNEG.GE.0.D0) THEN
         H=XEND-X
         LAST=.TRUE.
      END IF
      HOPT=H
      FACCON=1.D0
      CFAC=SAFE*(1+2*NIT)
      NSING=0
c-----------------------------------------------------------------------
      N2=2*N
      N3=3*N
      IF (ITOL.EQ.0) THEN
          DO I=1,N
             SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
          END DO
      ELSE
          DO J=1,N/NPDE
             DO I=1,NPDE
                SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I+(J-1)*NPDE))
             END DO
          END DO
      END IF
      HHFAC=H
c      write(*,*) 'calling fcn'
      CALL FCN(N,X,Y,Y0,RPAR,IPAR,f,fvec,derivf,bndxa,bndxb, vec)
c      write(*,*) 'leaving fcn'
      NFCN=NFCN+1
C --- BASIC INTEGRATION STEP
  10  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
      NJAC=NJAC+1
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
c-----------------------------------------------------------------------
c      write(*,*) 'calling jac'
      CALL JAC(N,X,Y,FJAC,RPAR,IPAR, f, derivf, bndxa, difbxa, bndxb, 
     &         difbxb)
c      write(*,*) 'leaving jac'
c-----------------------------------------------------------------------
      CALJAC=.TRUE.
  20  CONTINUE
C --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS
      FAC1=U1/H
      ALPHN=ALPH/H
      BETAN=BETA/H
c-----------------------------------------------------------------------
c      write(*,*) 'calling decomr'
      CALL DECOMR(N,NPDE,NINT,KCOL,FJAC,FMAS,FAC1,E1,IP1,IER)
c      write(*,*) 'leaving decomr'
c-----------------------------------------------------------------------
      IF (IER.NE.0) GOTO 78
c-----------------------------------------------------------------------
c      write(*,*) 'calling decomc'
      CALL DECOMC(N,NPDE,NINT,KCOL,FJAC,FMAS,ALPHN,BETAN,E2R,IP2,IER)
c      write(*,*) 'leaving decomc'
c-----------------------------------------------------------------------
      IF (IER.NE.0) GOTO 78
      NDEC=NDEC+1
  30  CONTINUE
      NSTEP=NSTEP+1
      IF (NSTEP.GT.NMAX) GOTO 178
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND) GOTO 177
c-----------------------------------------------------------------------
      do 31 i = 1, npde
         scal(i) = scal(i)/hhfac
   31 continue
c      do 32 i = npde*(nint*kcol+1)+1, npde*(nint*kcol+2)
c         scal(i) = scal(i)/hhfac
c   32 continue
c      do 33 i = npde*(nint*kcol+2)+1, npde*(nint*kcol+3)
c         scal(i) = scal(i)/hhfac
c   33 continue
      do 34 i = n-npde+1, n
         scal(i) = scal(i)/hhfac
   34 continue
c-----------------------------------------------------------------------
      XPH=X+H
C *** *** *** *** *** *** ***
C  STARTING VALUES FOR NEWTON ITERATION
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN) THEN
         DO I=1,N
            Z1(I)=0.D0
            Z2(I)=0.D0
            Z3(I)=0.D0
            F1(I)=0.D0
            F2(I)=0.D0
            F3(I)=0.D0
         END DO
      ELSE
         C3Q=H/HOLD
         C1Q=C1*C3Q
         C2Q=C2*C3Q
         DO I=1,N
            AK1=CONT(I+N)
            AK2=CONT(I+N2)
            AK3=CONT(I+N3)
            Z1I=C1Q*(AK1+(C1Q-C2M1)*(AK2+(C1Q-C1M1)*AK3))
            Z2I=C2Q*(AK1+(C2Q-C2M1)*(AK2+(C2Q-C1M1)*AK3))
            Z3I=C3Q*(AK1+(C3Q-C2M1)*(AK2+(C3Q-C1M1)*AK3))
            Z1(I)=Z1I
            Z2(I)=Z2I
            Z3(I)=Z3I
            F1(I)=TI11*Z1I+TI12*Z2I+TI13*Z3I
            F2(I)=TI21*Z1I+TI22*Z2I+TI23*Z3I
            F3(I)=TI31*Z1I+TI32*Z2I+TI33*Z3I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
            NEWT=0
            FACCON=MAX(FACCON,UROUND)**0.8D0
            THETA=ABS(THET)
  40        CONTINUE
            IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
            DO I=1,N
               CONT(I)=Y(I)+Z1(I)
            END DO
c            write(*,*) 'calling fcn1'
            CALL FCN(N,X+C1*H,CONT,Z1,RPAR,IPAR,f,fvec,derivf,bndxa,
     &               bndxb, vec)
            DO I=1,N
               CONT(I)=Y(I)+Z2(I)
            END DO
c            write(*,*) 'calling fcn2'
            CALL FCN(N,X+C2*H,CONT,Z2,RPAR,IPAR,f,fvec,derivf,bndxa,
     &               bndxb, vec)
            DO I=1,N
               CONT(I)=Y(I)+Z3(I)
            END DO
c            write(*,*) 'calling fcn3'
            CALL FCN(N,XPH,CONT,Z3,RPAR,IPAR,f,fvec,derivf,bndxa,bndxb,
     &               vec)
            NFCN=NFCN+3
c            write(*,*) 'all fcn called'
C ---     SOLVE THE LINEAR SYSTEMS
           DO I=1,N
              A1=Z1(I)
              A2=Z2(I)
              A3=Z3(I)
              Z1(I)=TI11*A1+TI12*A2+TI13*A3
              Z2(I)=TI21*A1+TI22*A2+TI23*A3
              Z3(I)=TI31*A1+TI32*A2+TI33*A3
           END DO
c-----------------------------------------------------------------------
c           write(*,*) 'calling slvrad'
        CALL SLVRAD(N,NPDE,KCOL,NINT,FMAS,FAC1,ALPHN,BETAN,E1,E2R,
     &              Z1,Z2,Z3,F1,F2,F3,CONT,IP1,IP2,CZ2)
c           write(*,*) 'leaving slvrad'
c-----------------------------------------------------------------------
            NSOL=NSOL+1
            NEWT=NEWT+1
            DYNO=0.D0
            DO I=1,N
               DENOM=SCAL(I)
               DYNO=DYNO+(Z1(I)/DENOM)**2+(Z2(I)/DENOM)**2
     &          +(Z3(I)/DENOM)**2
            END DO
            DYNO=SQRT(DYNO/N3)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
            IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
                THQ=DYNO/DYNOLD
                IF (NEWT.EQ.2) THEN
                   THETA=THQ
                ELSE
                   THETA=SQRT(THQ*THQOLD)
                END IF
                THQOLD=THQ
                IF (THETA.LT.0.99D0) THEN
                    FACCON=THETA/(1.0D0-THETA)
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
                    IF (DYTH.GE.1.0D0) THEN
                         QNEWT=MAX(1.0D-4,MIN(20.0D0,DYTH))
                         HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
                         H=HHFAC*H
                         REJECT=.TRUE.
                         LAST=.FALSE.
                         IF (CALJAC) GOTO 20
                         GOTO 10
                    END IF
                ELSE
                    GOTO 78
                END IF
            END IF
            DYNOLD=MAX(DYNO,UROUND)
            DO I=1,N
               F1I=F1(I)+Z1(I)
               F2I=F2(I)+Z2(I)
               F3I=F3(I)+Z3(I)
               F1(I)=F1I
               F2(I)=F2I
               F3(I)=F3I
               Z1(I)=T11*F1I+T12*F2I+T13*F3I
               Z2(I)=T21*F1I+T22*F2I+T23*F3I
               Z3(I)=T31*F1I+    F2I
            END DO
            IF (FACCON*DYNO.GT.FNEWT) GOTO 40
C --- ERROR ESTIMATION
c      write(*,*) 'calling estrad'
      CALL ESTRAD (N,NPDE,KCOL,NINT,FMAS,H,DD1,DD2,DD3,FCN,NFCN,
     &          Y0,Y,X,E1,Z1,Z2,Z3,CONT,F1,F2,IP1,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR, f, fvec, derivf, bndxa, 
     &          difbxa, bndxb, difbxb, uinit, vec)
c      write(*,*) 'leaving estrad'
C --- COMPUTATION OF HNEW
C --- WE REQUIRE .2<=HNEW/H<=8.
      FAC=MIN(SAFE,CFAC/(NEWT+2*NIT))
      QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC))
      HNEW=H/QUOT
C *** *** *** *** *** *** ***
C  IS THE ERROR SMALL ENOUGH ?
C *** *** *** *** *** *** ***
      IF (ERR.LT.1.D0) THEN
C --- STEP IS ACCEPTED
         FIRST=.FALSE.
         NACCPT=NACCPT+1
         IF (PRED) THEN
C       --- PREDICTIVE CONTROLLER OF GUSTAFSSON
            IF (NACCPT.GT.1) THEN
               FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE
               FACGUS=MAX(FACR,MIN(FACL,FACGUS))
               QUOT=MAX(QUOT,FACGUS)
               HNEW=H/QUOT
            END IF
            HACC=H
            ERRACC=MAX(1.0D-2,ERR)
         END IF
         HOLD=H
         X=XPH
         DO I=1,N
            Y(I)=Y(I)+Z3(I)
            Z2I=Z2(I)
            Z1I=Z1(I)
            CONT(I+N)=(Z2I-Z3(I))/C2M1
            AK=(Z1I-Z2I)/C1MC2
            ACONT3=Z1I/C1
            ACONT3=(AK-ACONT3)/C2
            CONT(I+N2)=(AK-CONT(I+N))/C1M1
            CONT(I+N3)=CONT(I+N2)-ACONT3
         END DO
         IF (ITOL.EQ.0) THEN
             DO I=1,N
                SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
             END DO
         ELSE
             DO J=1,N/NPDE
                DO I=1,NPDE
                   SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I+(J-1)*NPDE))
                END DO
             END DO
         END IF
         CALL SOLOUT(N,X,Y,ITOL,RPAR,IPAR,IRTRN)
         IF (IRTRN.NE.0) GOTO 179
         CALJAC=.FALSE.
         IF (LAST) THEN
            H=HOPT
            IDID=1
            RETURN
         END IF
         CALL FCN(N,X,Y,Y0,RPAR,IPAR,f,fvec,derivf,bndxa,bndxb, vec)
         NFCN=NFCN+1
         HNEW=POSNEG*MIN(ABS(HNEW),HMAXN)
         HOPT=HNEW
         HOPT=MIN(H,HNEW)
         IF (REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
         REJECT=.FALSE.
         IF ((X+HNEW/QUOT1-XEND)*POSNEG.GE.0.D0) THEN
            H=XEND-X
            LAST=.TRUE.
         ELSE
            QT=HNEW/H
            HHFAC=H
            IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) GOTO 30
            H=HNEW
         END IF
         HHFAC=H
         IF (THETA.LE.THET) GOTO 20
         GOTO 10
      ELSE
C --- STEP IS REJECTED
         REJECT=.TRUE.
         LAST=.FALSE.
         IF (FIRST) THEN
             H=H*0.1D0
             HHFAC=0.1D0
         ELSE
             HHFAC=HNEW/H
             H=HNEW
         END IF
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         IF (CALJAC) GOTO 20
         GOTO 10
      END IF
C --- UNEXPECTED STEP-REJECTION
  78  CONTINUE
      IF (IER.NE.0) THEN
          NSING=NSING+1
          IF (NSING.GE.5) GOTO 176
      END IF
      H=H*0.5D0
      HHFAC=0.5D0
      REJECT=.TRUE.
      LAST=.FALSE.
      IF (CALJAC) GOTO 20
      GOTO 10
C --- FAIL EXIT
 176  CONTINUE
      WRITE(6,979)X
      WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER
      IDID=-4
      RETURN
 177  CONTINUE
      WRITE(6,979)X
      WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H
      IDID=-3
      RETURN
 178  CONTINUE
      WRITE(6,979)X
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'
      IDID=-2
      RETURN
C --- EXIT CAUSED BY SOLOUT
 179  CONTINUE
 979  FORMAT(' EXIT OF RADAU5 AT X=',E18.4)
      IDID=2
      RETURN
      END

C=======================================================================
C     DEBOOR, BSPLINE PACKAGE
C=======================================================================
      SUBROUTINE BSPLVDri ( XT, K, X, ILEFT, VNIKX, NDERIV )
      implicit none
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C CALCULATES THE VALUE AND THE FIRST NDERIV-1 DERIVATIVES OF ALL
C B-SPLINES WHICH DO NOT VANISH AT X.  THE ROUTINE FILLS THE TWO-
C DIMENSIONAL ARRAY VNIKX(J,IDERIV), J=IDERIV, ... ,K WITH NONZERO
C VALUES OF B-SPLINES OF ORDER K+1-IDERIV, IDERIV=NDERIV, ... ,1, BY
C REPEATED CALLS TO BSPLVN.
C
C LAST MODIFIED BY RONG WANG, DEC 13, 2006.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  BSPLVN
C USER ROUTINES CALLED..     NONE
C CALLED BY..                COLPNT,INITAL,VALUES
C FORTRAN FUNCTIONS USED..   DBLE,MAX
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      INTEGER K,NDERIV,ILEFT
      DOUBLE PRECISION X
      DOUBLE PRECISION XT(*),VNIKX(K,NDERIV)
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER KO,IDERIV,IDERVM,KMD,JM1,IPKMD,JLOW
      DOUBLE PRECISION A(20,20)
      DOUBLE PRECISION FKMD,DIFF,V
C-----------------------------------------------------------------------
C CONSTANT
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.D0)
      PARAMETER (ONE  = 1.D0)
C-----------------------------------------------------------------------
C LOOP INDICES
      INTEGER I,J,M,L
C-----------------------------------------------------------------------
      KO = K + 1 - NDERIV
      CALL BSPLVNri(XT,KO,1,X,ILEFT,VNIKX(NDERIV,NDERIV))
      IF (NDERIV .LE. 1) GO TO 130
      IDERIV = NDERIV
      DO 20 I=2,NDERIV
        IDERVM = IDERIV-1
        DO 10 J=IDERIV,K
          VNIKX(J-1,IDERVM) = VNIKX(J,IDERIV)
   10   CONTINUE
        IDERIV = IDERVM
        CALL BSPLVNri(XT,0,2,X,ILEFT,VNIKX(IDERIV,IDERIV))
   20 CONTINUE
      DO 40 I=1,K
        DO 30 J=1,K
          A(I,J) = ZERO
   30   CONTINUE
        A(I,I) = ONE
   40 CONTINUE
      KMD = K
      DO 120 M=2,NDERIV
        KMD = KMD - 1
        FKMD =  DBLE(KMD)
        I = ILEFT
        J = K
   50   CONTINUE
        JM1 = J-1
        IPKMD = I + KMD
        DIFF = XT(IPKMD) -XT(I)
        IF (JM1 .NE. 0) THEN
          IF (DIFF .NE. ZERO) THEN
            DO 60 L=1,J
              A(L,J) = (A(L,J) - A(L,J-1))/DIFF*FKMD
   60       CONTINUE
          ENDIF
          J = JM1
          I = I - 1
          GO TO 50
        ENDIF
        IF (DIFF .NE. ZERO) THEN
          A(1,1) = A(1,1)/DIFF*FKMD
        ENDIF
        DO 110 I=1,K
          V = ZERO
          JLOW = MAX(I,M)
          DO 100 J=JLOW,K
            V = A(I,J)*VNIKX(J,M) + V
  100     CONTINUE
          VNIKX(I,M) = V
  110   CONTINUE
  120 CONTINUE
  130 RETURN
      END

      SUBROUTINE BSPLVNri( XT, JHIGH, INDEX, X, ILEFT, VNIKX )
      implicit none
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT THE
C POINT X OF ORDER MAX(JHIGH,(J+1)(INDEX-1)) FOR THE BREAKPOINT SEQ-
C UENCE XT.  ASSUMING THAT XT(ILEFT) .LE. X .LE. XT(ILEFT+1), THE ROUT-
C INE RETURNS THE B-SPLINE VALUES IN THE ONE DIMENSIONAL ARRAY VNIKX.
C
C LAST MODIFIED BY RONG WANG, JAN 8, 2001.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  NONE
C USER ROUTINES CALLED..     NONE
C CALLED BY..                BSPLVD
C FORTRAN FUNCTIONS USED..   NONE
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      DOUBLE PRECISION XT(*),X,VNIKX(*)
      INTEGER JHIGH,INDEX,ILEFT
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER IPJ,IMJP1,JP1,JP1ML
      DOUBLE PRECISION VMPREV,VM
C-----------------------------------------------------------------------
C CONSTANT
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.D0)
      PARAMETER (ONE  = 1.D0)
      DOUBLE PRECISION DELTAM(20),DELTAP(20)
      INTEGER J
C-----------------------------------------------------------------------
C LOOP INDICE
      INTEGER L
C-----------------------------------------------------------------------
      DATA J/1/,DELTAM/20*0.D+0/,DELTAP/20*0.D+0/

      IF(INDEX.EQ.1) THEN
        J = 1
        VNIKX(1) = ONE
        IF (J .GE. JHIGH) GO TO 40
      ENDIF
   20 CONTINUE
      IPJ = ILEFT+J
      DELTAP(J) = XT(IPJ) - X
      IMJP1 = ILEFT-J+1
      DELTAM(J) = X - XT(IMJP1)
      VMPREV = ZERO
      JP1 = J+1
      DO 30 L=1,J
        JP1ML = JP1-L
        VM = VNIKX(L)/(DELTAP(L) + DELTAM(JP1ML))
        VNIKX(L) = VM*DELTAP(L) + VMPREV
        VMPREV = VM*DELTAM(JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      J = JP1
      IF (J .LT. JHIGH) GO TO 20
   40 RETURN
      END

      SUBROUTINE INTERVri ( XT, LXT, X, ILEFT, MFLAG, ILO)
      implicit none
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C COMPUTES LARGEST ILEFT IN (1,LXT) SUCH THAT XT(ILEFT) .LE. X.  THE
C PROGRAM STARTS THE SEARCH FOR ILEFT WITH THE VALUE OF ILEFT THAT WAS
C RETURNED AT THE PREVIOUS CALL (AND WAS SAVED IN THE LOCAL VARIABLE
C ILO) TO MINIMIZE THE WORK IN THE COMMON CASE THAT THE VALUE OF X ON
C THIS CALL IS CLOSE TO THE VALUE OF X ON THE PREVIOUS CALL.  SHOULD
C THIS ASSUMPTION NOT BE VALID, THEN THE PROGRAM LOCATES ILO AND IHI
C SUCH THAT XT(ILO) .LE. X .LT. XT(IHI) AND, ONCE THEY ARE FOUND USES
C BISECTION TO FIND THE CORRECT VALUE FOR ILEFT.
C
C LAST MODIFIED BY RONG WANG, JAN 9, 2001.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  NONE
C USER ROUTINES CALLED..     NONE
C CALLED BY..                COLPNT,INITAL,VALUES
C FORTRAN FUNCTIONS USED..   NONE
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      INTEGER LXT,ILEFT,MFLAG
      DOUBLE PRECISION XT(LXT),X
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER ILO,IHI,ISTEP,MIDDLE
C-----------------------------------------------------------------------
      IF(MFLAG.EQ.-2) ILO = 1
      IHI = ILO + 1
      IF (IHI .LT. LXT) GO TO 20
      IF (X .GE. XT(LXT)) GO TO 110
      IF (LXT .LE. 1) GO TO 90
      ILO = LXT - 1
      GO TO 21
   20 IF (X .GE. XT(IHI)) GO TO 40
   21 IF (X .GE. XT(ILO)) GO TO 100
C-----------------------------------------------------------------------
C NOW X .LT. XT(IHI).  FIND LOWER BOUND.
C-----------------------------------------------------------------------
      ISTEP = 1
   31 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO .LE. 1) GO TO 35
      IF (X .GE. XT(ILO)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 31
   35 ILO = 1
      IF (X .LT. XT(1)) GO TO 90
      GO TO 50
C-----------------------------------------------------------------------
C NOW X .GE. XT(ILO).  FIND UPPER BOUND.
C-----------------------------------------------------------------------
   40 ISTEP = 1
   41 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .GE. LXT) GO TO 45
      IF (X .LT. XT(IHI)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 41
   45 IF (X .GE. XT(LXT)) GO TO 110
      IHI = LXT
C-----------------------------------------------------------------------
C NOW XT(ILO) .LE. X .LT. XT(IHI).  NARROW THE INTERVAL.
C-----------------------------------------------------------------------
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO) GO TO 100
C-----------------------------------------------------------------------
C NOTE..  IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1.
C-----------------------------------------------------------------------
      IF (X .LT. XT(MIDDLE)) GO TO 53
      ILO = MIDDLE
      GO TO 50
   53 IHI = MIDDLE
      GO TO 50
C-----------------------------------------------------------------------
C SET OUTPUT AND RETURN.
C-----------------------------------------------------------------------
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END



C=======================================================================
C     EISPACK Routines
C=======================================================================
      SUBROUTINE IMTQL1ri(N,D,E,IERR) 
      implicit none
C 
      INTEGER I,J,L,M,N,II,MML,IERR 
      DOUBLE PRECISION D(N),E(N) 
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG 
C 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL1, 
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON, 
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE. 
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971). 
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC 
C     TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD. 
C 
C     ON INPUT 
C 
C        N IS THE ORDER OF THE MATRIX. 
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX. 
C 
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY. 
C 
C      ON OUTPUT 
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN 
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND 
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE 
C          THE SMALLEST EIGENVALUES. 
C 
C        E HAS BEEN DESTROYED. 
C 
C        IERR IS SET TO 
C          ZERO       FOR NORMAL RETURN, 
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN 
C                     DETERMINED AFTER 30 ITERATIONS. 
C 
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) . 
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, 
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
C 
C     THIS VERSION DATED APRIL 1983. 
C 
C     ------------------------------------------------------------------ 
C 
      IERR = 0 
      IF (N .EQ. 1) GO TO 1001 
C 
      DO 100 I = 2, N 
         E(I-1) = E(I) 
  100 CONTINUE
C 
      E(N) = 0.0D0 
C 
      DO 290 L = 1, N 
         J = 0 
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT .......... 
  105    DO 110 M = L, N 
            IF (M .EQ. N) GO TO 120 
            TST1 = ABS(D(M)) + ABS(D(M+1)) 
            TST2 = TST1 + ABS(E(M)) 
            IF (TST2 .EQ. TST1) GO TO 120 
  110    CONTINUE 
C 
  120    P = D(L) 
         IF (M .EQ. L) GO TO 215 
         IF (J .EQ. 30) GO TO 1000 
         J = J + 1 
C     .......... FORM SHIFT .......... 
         G = (D(L+1) - P) / (2.0D0 * E(L)) 
         R = PYTHAG(G,1.0D0) 
         G = D(M) - P + E(L) / (G + SIGN(R,G)) 
         S = 1.0D0 
         C = 1.0D0 
         P = 0.0D0 
         MML = M - L 
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... 
         DO 200 II = 1, MML 
            I = M - II 
            F = S * E(I) 
            B = C * E(I) 
            R = PYTHAG(F,G) 
            E(I+1) = R 
            S = F / R 
            C = G / R 
            G = D(I+1) - P 
            R = (D(I) - G) * S + 2.0D0 * C * B 
            P = S * R 
            D(I+1) = G + P 
            G = C * R - B 
  200    CONTINUE 
C 
         D(L) = D(L) - P 
         E(L) = G 
         E(M) = 0.0D0 
         GO TO 105 
C     .......... ORDER EIGENVALUES .......... 
  215    IF (L .EQ. 1) GO TO 250 
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- .......... 
         DO 230 II = 2, L 
            I = L + 2 - II 
            IF (P .GE. D(I-1)) GO TO 270 
            D(I) = D(I-1) 
  230    CONTINUE 
C 
  250    I = 1 
  270    D(I) = P 
  290 CONTINUE 
C 
      GO TO 1001 
C     .......... SET ERROR -- NO CONVERGENCE TO AN 
C                EIGENVALUE AFTER 30 ITERATIONS .......... 
 1000 IERR = L 
 1001 RETURN 
      END 
      SUBROUTINE IMTQL2ri(NM,N,D,E,Z,IERR) 
      implicit none
C 
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR 
      DOUBLE PRECISION D(N),E(N),Z(NM,N) 
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG 
C 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2, 
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON, 
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE. 
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971). 
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS 
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD. 
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO 
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS 
C     FULL MATRIX TO TRIDIAGONAL FORM. 
C 
C     ON INPUT 
C 
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
C          DIMENSION STATEMENT. 
C 
C        N IS THE ORDER OF THE MATRIX. 
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX. 
C 
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY. 
C 
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE 
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS 
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN 
C          THE IDENTITY MATRIX. 
C 
C      ON OUTPUT 
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN 
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT 
C          UNORDERED FOR INDICES 1,2,...,IERR-1. 
C 
C        E HAS BEEN DESTROYED. 
C 
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC 
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE, 
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED 
C          EIGENVALUES. 
C 
C        IERR IS SET TO 
C          ZERO       FOR NORMAL RETURN, 
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN 
C                     DETERMINED AFTER 30 ITERATIONS. 
C 
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) . 
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, 
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
C 
C     THIS VERSION DATED APRIL 1983. 
C 
C     ------------------------------------------------------------------ 
C 
      IERR = 0 
      IF (N .EQ. 1) GO TO 1001 
C 
      DO 100 I = 2, N 
         E(I-1) = E(I) 
  100 CONTINUE
C 
      E(N) = 0.0D0 
C 
      DO 240 L = 1, N 
         J = 0 
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT .......... 
  105    DO 110 M = L, N 
            IF (M .EQ. N) GO TO 120 
            TST1 = ABS(D(M)) + ABS(D(M+1)) 
            TST2 = TST1 + ABS(E(M)) 
            IF (TST2 .EQ. TST1) GO TO 120 
  110    CONTINUE 
C 
  120    P = D(L) 
         IF (M .EQ. L) GO TO 240 
         IF (J .EQ. 30) GO TO 1000 
         J = J + 1 
C     .......... FORM SHIFT .......... 
         G = (D(L+1) - P) / (2.0D0 * E(L)) 
         R = PYTHAG(G,1.0D0) 
         G = D(M) - P + E(L) / (G + SIGN(R,G)) 
         S = 1.0D0 
         C = 1.0D0 
         P = 0.0D0 
         MML = M - L 
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... 
         DO 200 II = 1, MML 
            I = M - II 
            F = S * E(I) 
            B = C * E(I) 
            R = PYTHAG(F,G) 
            E(I+1) = R 
            S = F / R 
            C = G / R 
            G = D(I+1) - P 
            R = (D(I) - G) * S + 2.0D0 * C * B 
            P = S * R 
            D(I+1) = G + P 
            G = C * R - B 
C     .......... FORM VECTOR .......... 
            DO 180 K = 1, N 
               F = Z(K,I+1) 
               Z(K,I+1) = S * Z(K,I) + C * F 
               Z(K,I) = C * Z(K,I) - S * F 
  180       CONTINUE 
C 
  200    CONTINUE 
C 
         D(L) = D(L) - P 
         E(L) = G 
         E(M) = 0.0D0 
         GO TO 105 
  240 CONTINUE 
C     .......... ORDER EIGENVALUES AND EIGENVECTORS .......... 
      DO 300 II = 2, N 
         I = II - 1 
         K = I 
         P = D(I) 
C 
         DO 260 J = II, N 
            IF (D(J) .GE. P) GO TO 260 
            K = J 
            P = D(J) 
  260    CONTINUE 
C 
         IF (K .EQ. I) GO TO 300 
         D(K) = D(I) 
         D(I) = P 
C 
         DO 280 J = 1, N 
            P = Z(J,I) 
            Z(J,I) = Z(J,K) 
            Z(J,K) = P
  280    CONTINUE
C 
  300 CONTINUE
C 
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN 
C                EIGENVALUE AFTER 30 ITERATIONS .......... 
 1000 IERR = L 
 1001 RETURN
      END

C=======================================================================
C     GAULEG
C=======================================================================
      SUBROUTINE GAULEGri(N, NSQ, PTS, WTS, WORK, FLAG)
      implicit none
C $ID: GAULEG.F,V 1.11 1992/06/25 15:09:31 KEAST EXP $
C LAST MODIFIED BY RONG WANG, 2001/03/08
C
C     GAULEG RETURNS THE POINTS AND WEIGHTS FOR GAUSS-LEGENDRE
C     QUADRATURE OR GAUSS-LOBATTO QUADRATURE OVER THE INTERVAL
C     [-1,1] OR [0,1].
C
C     ON INPUT:
C
C        N     : THE NUMBER OF GAUSS-LEGENDRE POINTS.
C        NSQ   : EQUAL TO N*N, TO HANDLE STORAGE FOR EIGENVECTORS.
C        PTS   : DOUBLE PRECISION (N).
C        WTS   : DOUBLE PRECISION (N).  IF FLAG (SEE BELOW) IS 1 OR 3,
C                WTS IS USED ONLY FOR TEMPORARY WORKSPACE.
C        WORK  : DOUBLE PRECISION (NSQ), WORK SPACE FOR CALL TO IMTQL2,
C                IF WEIGHTS ARE REQUIRED. IF FLAG IS 1 OR 3, WORK IS
C                NOT REFERENCED, AND MAY BE DECLARED AS SCALAR IN THE
C                CALLING PROGRAM.
C        FLAG  : SPECIFIES WHETHER WEIGHTS ARE ALSO REQUIRED, AND
C                WHETHER GAUSS-LEGENDRE OR LOBATTO POINTS ARE WANTED.
C                   FLAG  = 1: GAUSS-LEGENDRE POINTS ONLY OVER [-1,1];
C                         = 2: GAUSS-LEGENDRE POINTS ONLY OVER [0,1];
C                         = 3: GAUSS-LEGENDRE POINTS AND WEIGHTS OVER
C                              [-1,1];
C                         = 4: GAUSS-LEGENDRE POINTS AND WEIGHTS OVER
C                              [0,1];
C                         = 5: LOBATTO POINTS ONLY OVER [-1,1];
C                         = 6: LOBATTO POINTS ONLY OVER [0,1];
C                         = 7: LOBATTO POINTS AND WEIGHTS OVER [-1,1];
C                         = 8: LOBATTO POINTS AND WEIGHTS OVER [0,1];
C                FOR ANY OTHER VALUE, THE DEFAULT IS GAUSS-LEGENDRE
C                POINTS AND WEIGHTS OVER [-1,1].
C
C     ON OUTPUT:
C
C        FOR FLAG <> 1 OR 2:
C          PTS : PTS(I) IS THE ITH GAUSS-LEGENDRE POINT IN [-1,1],
C                PTS(I) < PTS(I+1), I = 1,2,..,N-1.
C          WTS : WTS(I) IS THE ITH GAUSS-LEGENDRE WEIGHT IF FLAG <> 1.
C
C        FOR FLAG = 3 OR 4:
C          PTS : PTS(I) IS THE ITH LOBATTO POINT IN [-1,1],
C                PTS(I) < PTS(I+1), I = 1,2,..,N-1.
C                CLEARLY, PTS(1) = -1.0, PTS(N) = 1.0.
C          WTS : WTS(I) IS THE ITH LOBATTO WEIGHT IF FLAG = 4.
C
C        WORK  : WORK, USED TO STORE EIGENVECTORS, UNREFERENCED IF FLAG
C                IS 1 OR 3.
C
C     SUBROUTINES USED:
C
C        IMQTL1: EISPACK ROUTINE TO COMPUTE THE EIGENVALUES OF A
C                SYMMETRIC TRIDIAGONAL MATRIX.
C
C        IMQTL2: EISPACK ROUTINE TO COMPUTE THE EIGNEVECTORS AND
C                EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX.
C
C     FUNCTIONS USED:
C
C        PYTHAG: EISPACK FUNCTION TO COMPUTE EUCLIDEAN NORM.
C
C     INTRINSIC FUNCTIONS USED:
C
C        SQRT, DBLE.
C
C     VERSION: JUNE 22 1992, PAT KEAST.
C
C     DECLARATIONS:
C
C        PARAMETERS:
C
      INTEGER          N, NSQ, FLAG
      DOUBLE PRECISION PTS(N), WTS(N), WORK(NSQ)
C
C        LOCAL VARIABLES:

      INTEGER          IFAIL, J, NM2
      DOUBLE PRECISION FOUR, THREE, TWO, ONE, ZERO
*     .. EXTERNAL FUNCTIONS ..
      EXTERNAL         IMTQL2ri
      PARAMETER ( FOUR = 4.0D0, THREE = 3.0D0, TWO = 2.0D0,
     *            ONE = 1.0D0, ZERO = 0.0D0 )

      DO 10 J = 1,N
         PTS(J) = ZERO
   10 CONTINUE

      DO 20 J = 1,NSQ
         WORK(J) = ZERO
   20 CONTINUE

      DO 30 J = 1,N
         WORK((J-1)*N+J) = ONE
   30 CONTINUE

      IF ( FLAG .LE.4 .OR. FLAG .GT. 8 ) THEN
C        GAUS-LEGENDRE POINTS AND WEIGHTS.
         DO 40 J = 1,N-1
            WTS(J+1) = DBLE(J)/SQRT(DBLE(4*J*J)-ONE)
   40    CONTINUE

         IF ( FLAG .EQ. 1 .OR. FLAG .EQ. 2 ) THEN
C           COMPUTE ONLY THE GAUSS-LEGENDRE POINTS OVER [-1,1].
            CALL IMTQL1ri(N, PTS, WTS, IFAIL )
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 2 ) THEN
               DO 45 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
   45          CONTINUE
            ENDIF
         ELSE
C           COMPUTE BOTH POINTS AND WEIGHTS OVER [-1,1].
            IFAIL = 1
C
            CALL IMTQL2ri(N, N, PTS, WTS, WORK, IFAIL)

            DO 50 J = 1,N
               WTS(J) = TWO*WORK((J-1)*N+1)*WORK((J-1)*N+1)
   50       CONTINUE
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 4 ) THEN
               DO 55 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
                  WTS(J) = WTS(J)/TWO
   55          CONTINUE
            ENDIF
         ENDIF
C
      ELSE
C        THE LOBATTO POINTS AND WEIGHTS.

C        FIRST, COMPUTE THE ORDER N-2 JACOBI POINTS AND/OR WEIGHTS.
         NM2 = N-2
         DO 60 J = 1,NM2-1
            WTS(J+1) = SQRT(DBLE(J*(J+2))/DBLE((2*J+1)*(2*J+3)))
   60    CONTINUE

         IF ( FLAG .EQ. 5 .OR. FLAG .EQ. 6) THEN
C           COMPUTE ONLY THE GAUSS-LOBATTO POINTS OVER [-1,1].
            CALL IMTQL1ri(NM2, PTS(2), WTS, IFAIL )
            PTS(1) = -ONE
            PTS(N) = ONE
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 6 ) THEN
               DO 65 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
   65          CONTINUE
            ENDIF
         ELSE
C           COMPUTE BOTH POINTS AND WEIGHTS.
            IFAIL = 1
C
            CALL IMTQL2ri(NM2, NM2, PTS(2), WTS, WORK, IFAIL)
            PTS(1) = -ONE
            PTS(N) = ONE

            DO 70 J = 2,N-1
               WTS(J) = (FOUR/THREE)*WORK((J-2)*NM2+1)*WORK((J-2)*NM2+1)
     *                  /(ONE - PTS(J)*PTS(J))
   70       CONTINUE
            WTS(1) = ZERO
            DO 80 J = 2,N-1
               WTS(1) = WTS(1) - WTS(J)
   80       CONTINUE
            WTS(1) = ONE + WTS(1)/TWO
            WTS(N) = WTS(1)
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 8 ) THEN
               DO 85 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
                  WTS(J) = WTS(J)/TWO
   85          CONTINUE
            ENDIF
         ENDIF
C
         RETURN
      ENDIF
*
      RETURN
      END

C=======================================================================
C     LAPACK Routines
C=======================================================================

C ***********************************************************
C
      SUBROUTINE DECOMC(N,NPDE,NINT,KCOL,FJAC,FMAS,ALPHN,BETAN,E2R,
     &                  IP2,IER)
      implicit none
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 23, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,IP2(N),IER
      INTEGER NPDE, NINT, KCOL
      DOUBLE PRECISION FJAC(*),FMAS(*),ALPHN,BETAN
      DOUBLE COMPLEX E2R(*),BB
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J
      INTEGER NEQ,NSIZTB,NSIZBK,KCOLTM,ITEMP,IABDTP,IABDBK,IABDBT
c-----------------------------------------------------------------------
c Subroutines Called:
c                               ccrcmp
c
c-----------------------------------------------------------------------

      BB = CMPLX(ALPHN,BETAN)
      ITEMP = 0
c      DO 100 I = 1, 2
      i = 1
         KCOLTM = KCOL+I-1
         NSIZTB = NPDE*NPDE*NCONTI
         NSIZBK = NPDE*NPDE*KCOLTM*(KCOLTM+NCONTI)*NINT

         IABDTP = ITEMP  + 1
         IABDBK = IABDTP + NSIZTB
         IABDBT = IABDBK + NSIZBK

         DO 10 J = 1,NSIZTB
            E2R(IABDTP-1+J)=BB*FJAC(IABDTP-1+J)
   10    CONTINUE

         DO 20 J = 1, NSIZBK
            E2R(IABDBK-1+J) = -FJAC(IABDBK-1+J)+BB*FMAS(IABDBK-1+J)
   20    CONTINUE

         DO 30 J = 1,NSIZTB
            E2R(IABDBT-1+J)=BB*FJAC(IABDBT-1+J)
   30    CONTINUE

         NEQ = NPDE*(KCOLTM*NINT+NCONTI)

c     BACOLRI --> BACOLRILAM
         CALL CCRCMPri(NEQ,E2R(IABDTP),NPDE,2*NPDE,E2R(IABDBK),
     &               KCOLTM*NPDE,(KCOLTM+NCONTI)*NPDE,NINT,E2R(IABDBT),
     &               NPDE,IP2(1+(I-1)*NPDE*(KCOL*NINT+NCONTI)),IER)
C          CALL CLAMDEC(NEQ,E2R(IABDTP),NPDE,2*NPDE,E2R(IABDBK),
C      &               KCOLTM*NPDE,(KCOLTM+NCONTI)*NPDE,NINT,E2R(IABDBT),
C      &               NPDE,IP2(1+(I-1)*NPDE*(KCOL*NINT+NCONTI)),IER)

         ITEMP  = NSIZTB + NSIZBK + NSIZTB

c  100 CONTINUE
C
      RETURN
      END

C ******************************************
C     VERSION OF SEPTEMBER 18, 1995
C ******************************************
C
      SUBROUTINE DECOMR(N,NPDE,NINT,KCOL,FJAC,FMAS,FAC1,E1,IP1,IER)
      implicit none
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 5, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,IP1(N),IER
      INTEGER NPDE, NINT, KCOL
      DOUBLE PRECISION FJAC(*),FMAS(*),FAC1,E1(*)
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J
      INTEGER NEQ,NSIZTB,NSIZBK,KCOLTM,ITEMP,IABDTP,IABDBK,IABDBT
c-----------------------------------------------------------------------
c Subroutines Called:
c                               crdcmp
c                               daxpy
c
c-----------------------------------------------------------------------

      ITEMP = 0
c      DO 100 I = 1, 1
      i = 1
         KCOLTM = KCOL+I-1
         NSIZTB = NPDE*NPDE*NCONTI
         NSIZBK = NPDE*NPDE*KCOLTM*(KCOLTM+NCONTI)*NINT

         IABDTP = ITEMP  + 1
         IABDBK = IABDTP + NSIZTB
         IABDBT = IABDBK + NSIZBK
         DO 10 J = 1,NSIZTB
            E1(IABDTP-1+J) = ZERO
   10    CONTINUE

         CALL DAXPYri(NSIZTB, FAC1, FJAC(IABDTP), 1, E1(IABDTP), 1)

         DO 20 J = 1, NSIZBK
            E1(IABDBK-1+J) = - FJAC(IABDBK-1+J)
   20    CONTINUE

         CALL DAXPYri(NSIZBK, FAC1, FMAS(IABDBK), 1, E1(IABDBK), 1)

         DO 30 J = 1,NSIZTB
            E1(IABDBT-1+J) = ZERO
   30    CONTINUE

         CALL DAXPYri(NSIZTB, FAC1, FJAC(IABDBT), 1, E1(IABDBT), 1)

         NEQ = NPDE*(KCOLTM*NINT+NCONTI)

c     BACOLRI --> BACOLRILAM
         CALL LAMDECri(NEQ,E1(IABDTP),NPDE,2*NPDE,E1(IABDBK),
     &               KCOLTM*NPDE,
     &               (KCOLTM+NCONTI)*NPDE,NINT,E1(IABDBT),NPDE,
     &               IP1(1+(I-1)*NPDE*(KCOL*NINT+NCONTI)),IER)

         ITEMP  = NSIZTB + NSIZBK + NSIZTB
c  100 CONTINUE
C
      RETURN
      END
C ***********************************************************
C
      SUBROUTINE ESTRAD(N,NPDE,KCOL,NINT,FMAS,H,DD1,DD2,DD3,FCN,
     &          NFCN,Y0,Y,X,E1,Z1,Z2,Z3,CONT,F1,F2,IP1,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR, f, fvec, derivf, bndxa,
     &          difbxa, bndxb, difbxb, uinit, vec)
      implicit none
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 30, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
        double precision        one
        parameter              (one = 1.0D0)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,NPDE,KCOL,NINT,NFCN,IP1(N),IPAR(*)
      DOUBLE PRECISION FMAS(*),H,DD1,DD2,DD3,Y0(N),Y(N),X,E1(*),Z1(N),
     &                 Z2(N),Z3(N),CONT(N),F1(N),F2(N),SCAL(N),ERR,
     &                 FAC1,RPAR(*),XIN(N)
      LOGICAL FIRST,REJECT
      EXTERNAL FCN
      external f, fvec, derivf, bndxa, difbxa, bndxb, difbxb, uinit
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J,K,M,III,II,KK,MM,NPDTP1,
     &        NPDBK1,NPDBT1
      DOUBLE PRECISION HEE1,HEE2,HEE3,SUM
      logical vec
c-----------------------------------------------------------------------
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
C
      DO 10 I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
   10 CONTINUE

      iii = npde*npde*kcol
      do 50 i = 1, nint
         do 40 k = 1, kcol
            do 30 m = 1, npde
               sum = zero
               ii = npde+(i-1)*npde*kcol+(k-1)*npde+m
               do 20 j = 1, kcol + nconti
                  kk = 1+(i-1)*iii*(kcol+nconti)+(j-1)*iii
     &               +(k-1)*npde+npde*npde*nconti
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  sum = sum + fmas(kk) * f1(mm)
   20          continue
               f2(ii)   = sum
               cont(ii) = sum + y0(ii)
   30       continue
   40    continue
   50 continue
      do 60 i = 1, npde
         f2(i) = zero
         cont(i) = - y0(i) * fac1
   60 continue
      do 70 i = 1, npde
         ii = i+npde*(kcol*nint+1)
         f2(ii) = zero
         cont(ii) = - y0(ii) * fac1
   70 continue

      npdtp1 = 1
      npdbk1 = npdtp1 + npde * npde * nconti
      npdbt1 = npdbk1 + npde * npde * nint * kcol * (kcol + nconti)

c      write(*,*) 'calling crslve'
c     BACOLRI --> BACOLRILAM
      call dcopyri(n,cont,1,xin,1)
      call lamsolri(e1(npdtp1), npde, 2*npde, e1(npdbk1), kcol*npde,
     &           (kcol+nconti)*npde, nint, e1(npdbt1), npde, ip1, xin,
     &           cont)
c      write(*,*) 'leaving crslve'

      err = zero

      do 90 i = 1, n
         err = err + (cont(i)/scal(i))**2
   90 continue

      err = max(sqrt(err/n), 1.d-10)

      if (err .lt. one) return

      if (first .or. reject) then

         do 100 i = 1, n
            cont(i) = y(i) + cont(i)
  100    continue
c            write(*,*) 'calling fcn'
         call fcn(n, x, cont, f1, rpar, ipar, f, fvec,  derivf, bndxa,
     &            bndxb, vec)
c            write(*,*) 'leaving fcn'
         nfcn = nfcn + 1
         do 110 i = 1, n
            cont(i) = f1(i) + f2(i)
  110    continue
         do 120 i = 1, npde
            cont(i) = - cont(i) * fac1
  120    continue
         do 150 i = n-npde+1, n
            cont(i) = - cont(i) * fac1
  150    continue

c     BACOLRI --> BACOLRILAM
         call dcopyri(n,cont,1,xin,1)
         call lamsolri(e1(npdtp1), npde, 2*npde, e1(npdbk1), kcol*npde,
     &               (kcol+nconti)*npde, nint, e1(npdbt1), npde, ip1,
     &               xin, cont)

         err = zero

         do 160 i = 1, n
            err = err + (cont(i)/scal(i))**2
  160    continue

         err = max(sqrt(err/n), 1.d-10)

      endif

      RETURN
      END
C ***********************************************************
C
      SUBROUTINE SLVRAD(N,NPDE,KCOL,NINT,FMAS,FAC1,ALPHN,BETAN,
     &                  E1,E2R,Z1,Z2,Z3,F1,F2,F3,CONT,IP1,IP2,CZ2)
      implicit none
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 23, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,NPDE,KCOL,NINT,IP1(N),IP2(N)
      DOUBLE PRECISION FMAS(*),FAC1,ALPHN,BETAN,E1(*),
     &                 Z1(N),Z2(N),Z3(N),F1(N),F2(N),F3(N),CONT(N),
     &                 XIN(N)
      DOUBLE COMPLEX E2R(*),CZ2(N),CXIN(N)
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J,K,M,L,III,II,KK,MM,NPDTP1,
     &        NPDBK1,NPDBT1
      DOUBLE PRECISION S1,S2,S3,BB

c-----------------------------------------------------------------------

      iii = npde*npde*kcol
      do 40 i = 1, nint
         do 30 k = 1, kcol
            do 20 m = 1, npde
               s1 = zero
               s2 = zero
               s3 = zero
               ii = npde+(i-1)*npde*kcol+(k-1)*npde+m
               do 10 j = 1, kcol + nconti
                  kk = 1+(i-1)*iii*(kcol+nconti)+(j-1)*iii
     &               +(k-1)*npde+npde*npde*nconti
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  bb = fmas(kk)
                  s1 = s1 - bb * f1(mm)
                  s2 = s2 - bb * f2(mm)
                  s3 = s3 - bb * f3(mm)
   10          continue
               z1(ii)   = z1(ii) + s1 * fac1
               z2(ii)   = z2(ii) + s2 * alphn - s3 * betan
               cont(ii) = z3(ii) + s3 * alphn + s2 * betan
   20       continue
   30    continue
   40 continue
      do 50 i = 1, npde
         z1(i) = - z1(i) * fac1
         cont(i) = z3(i)
   50 continue
      do 60 i = 1, npde
         ii = i+npde*(kcol*nint+1)
         z1(ii) = - z1(ii) * fac1
         cont(ii) = z3(ii)
   60 continue

      do 80 i = 1, n
         cz2(i) = cmplx(z2(i),cont(i))
   80 continue

      do 90 i = 1, npde
         cz2(i) = - cmplx(alphn,betan) * cz2(i)
   90 continue

      do 120 i = n-npde+1, n
         cz2(i) = - cmplx(alphn,betan) * cz2(i)
  120 continue

      npdtp1 = 1
      npdbk1 = npdtp1 + npde * npde * nconti
      npdbt1 = npdbk1 + npde * npde * nint * kcol * (kcol + nconti)

c     BACOLRI --> BACOLRILAM
      call dcopyri(N, z1, 1, xin, 1)
      call lamsolri(e1(npdtp1), npde, 2*npde, e1(npdbk1), kcol*npde,
     &            (kcol+nconti)*npde, nint, e1(npdbt1), npde, ip1, xin,
     &            z1)

c     BACOLRI --> BACOLRILAM
      call ccrslvri(e2r(npdtp1), npde, 2*npde, e2r(npdbk1), kcol*npde,
     &            (kcol+nconti)*npde, nint, e2r(npdbt1), npde, ip2, cz2)
      call ccopy(N, cz2, 1, cxin, 1)
c      write(*,*) cxin(5)
c      cxin = cz2
C       call clamsol(e2r(npdtp1), npde, 2*npde, e2r(npdbk1), kcol*npde,
C      &             (kcol+nconti)*npde, nint, e2r(npdbt1), npde, ip2,
C      &             cxin, cz2)

      do 130 i = 1, n
         z2(i) = dble(cz2(i))
         z3(i) = dimag(cz2(i))
  130 continue

      return
      end

C=======================================================================
C     BLAS Routines
C=======================================================================
      subroutine daxpyri(n,da,dx,incx,dy,incy)
      implicit none
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine dcopyri(n,dx,incx,dy,incy)
      implicit none
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      subroutine dscalri(n,da,dx,incx)
      implicit none
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

c=======================================================================
c     LAMPACK Routines
c=======================================================================
      SUBROUTINE LAMSOLri(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,X)
      implicit none
C
C***************************************************************
C
C  L A M S O L  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  L A M D E C.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  L A M D E C
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  L A M D E C
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  L A M D E C
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  L A M D E C
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - DOUBLE PRECISION(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C                    X - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,X,B
        DOUBLE PRECISION DOTPRD,BJ,XINCRJ,BINCRJ,SWAP
        INTEGER PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*),X(*)
        INTEGER NRWTOP, NOVRLP, NRWBLK, NCLBLK, NBLOKS, NRWBOT,
     *          NRWTP1, NRWBK1, NRWTP0, NRWBT1, NROWEL, NRWEL1,
     *          NVRLP0, NBLKS1, NBKTOP, NBKTP0, J, LOOP, INCR,
     *          K, INCRTP, INCRI, JPIVOT, JRWTOP, NRWBTL, L, 
     *          L1, IPLUSN, INCRN, IPVTN, NRWELL, IPVTI, NVRLP1,
     *          I, INCRJ, LL
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
        NBKTP0 = NBKTOP-1
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
       IF(NRWTOP.EQ.1)GO TO 135
        DO 130 J = 1,NRWTP0
              BJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*BJ
110           CONTINUE
130     CONTINUE
135    CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              BJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
       IF(NRWBK1.GT.NBKTP0)GO TO 275
           DO 270 J = NRWBK1,NBKTP0
                 INCRJ = INCR+J
                 JRWTOP = J -NRWTOP
                 BJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BJ
250              CONTINUE
270        CONTINUE
275        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           BJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*BJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           X(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -X(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*X(INCRJ)
440           CONTINUE
              X(INCRN) = DOTPRD/ARRAY(I,IPLUSN,K)
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = X(INCRN)
                 X(INCRN) = X(IPVTN)
                 X(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -X(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              X(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -X(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*X(J)
510        CONTINUE
           X(I) = DOTPRD/TOPBLK(I,I)
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = X(I)
                 X(I) = X(IPVTI)
                 X(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
        RETURN
        END

      SUBROUTINE LAMDECri(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *            NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
      implicit none
C
C***************************************************************
C
C  L A M D E C DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK
        DOUBLE PRECISION ROWMAX,ROWPIV,ROWMLT,COLMAX,COLPIV
        DOUBLE PRECISION SWAP,COLMLT,PIVMAX,ZERO,TEMPIV
        INTEGER N, NRWTOP, NOVRLP, NRWBLK, NCLBLK, NBLOKS, NRWBOT,
     *          IFLAG
        INTEGER PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*)
        DATA ZERO/0.0D0/
C
        INTEGER NRWTP1, NROWEL, NRWEL1, NVRLP0, I, IPLUS1, IPVT, J,
     *          L, K, KPLUS1, JPLUS1, JMINN, LOOP, INCRJ, IPLUSN, 
     *          INCRN, IRWBLK, IPVBLK, JRWBLK, INCR
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        IFLAG = 0
        PIVMAX = ZERO
        NRWTP1 = NRWTOP+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
        IF(N.NE.NBLOKS*NRWBLK+NOVRLP)GO TO 10
        IF(NOVRLP.NE.NRWTOP+NRWBOT)GO TO 10
        IF(NCLBLK.NE.NOVRLP+NRWBLK)GO TO 10
        IF(NOVRLP.GT.NRWBLK)GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        GO TO 100
10      CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IFLAG = 1
        RETURN
100     CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP ROW ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
        DO 190 I = 1,NRWTOP
           IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IPVT = I
           COLMAX = ABS(TOPBLK(I,I))
           DO 110 J = IPLUS1,NOVRLP
              TEMPIV = ABS(TOPBLK(I,J))
              IF(TEMPIV.LE.COLMAX)GO TO 110
                 IPVT = J
                 COLMAX = TEMPIV
110        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
           PIVMAX = MAX(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           PIVOT(I) = IPVT
           IF(IPVT.EQ.I)GO TO 140
              DO 120 L = I,NRWTOP
                 SWAP = TOPBLK(L,IPVT)
                 TOPBLK(L,IPVT) = TOPBLK(L,I)
                 TOPBLK(L,I) = SWAP
120           CONTINUE
              DO 130 L = 1,NRWBLK
                 SWAP = ARRAY(L,IPVT,1)
                 ARRAY(L,IPVT,1) = ARRAY(L,I,1)
                 ARRAY(L,I,1) = SWAP
130           CONTINUE
140        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM ROW
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           COLPIV = TOPBLK(I,I)
          IF(IPLUS1.GT.NRWTOP)GO TO 143
             DO 142 L = IPLUS1,NRWTOP
                TOPBLK(L,I) = TOPBLK(L,I)/COLPIV
142          CONTINUE
143       CONTINUE
          DO 144 L = 1,NRWBLK
             ARRAY(L,I,1) = ARRAY(L,I,1)/COLPIV
144       CONTINUE
           DO 180 J = IPLUS1,NOVRLP
              COLMLT = TOPBLK(I,J)
              IF(IPLUS1.GT.NRWTOP)GO TO 160
                 DO 150 L = IPLUS1,NRWTOP
                    TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
150              CONTINUE
160           CONTINUE
              DO 170 L = 1,NRWBLK
                 ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
170           CONTINUE
180        CONTINUE
190     CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
        INCR = 0
        DO 395 K = 1,NBLOKS
           KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
           DO 270 J = NRWTP1,NRWBLK
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = ABS(ARRAY(JMINN,J,K))
              LOOP = JMINN+1
              DO 210 I = LOOP,NRWBLK
                 TEMPIV = ABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.ROWMAX)GO TO 210
                 IPVT = I
                 ROWMAX = TEMPIV
210           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO  1000
              PIVMAX = MAX(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 230
                 DO 220 L = J,NCLBLK
                    SWAP = ARRAY(IPVT,L,K)
                    ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
                    ARRAY(JMINN,L,K) = SWAP
220              CONTINUE
230           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = ARRAY(JMINN,J,K)
              DO 240 I = LOOP,NRWBLK
                 ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
240           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 260 L = JPLUS1,NCLBLK
                 ROWMLT = ARRAY(JMINN,L,K)
                 DO 250 I = LOOP,NRWBLK
                    ARRAY(I,L,K) = ARRAY(I,L,K)
     *                                -ROWMLT*ARRAY(I,J,K)
250              CONTINUE
260           CONTINUE
270        CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP ROW ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
           DO 390 I = NRWEL1,NRWBLK
              IPLUSN = I+NRWTOP
              IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = IPLUSN
              COLMAX = ABS(ARRAY(I,IPVT,K))
              LOOP = IPLUSN+1
              DO 310 J = LOOP,NCLBLK
                 TEMPIV = ABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.COLMAX)GO TO 310
                 IPVT = J
                 COLMAX = TEMPIV
310           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = MAX(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRN = INCR+IPLUSN
              PIVOT(INCRN) = INCR+IPVT
              IRWBLK = IPLUSN-NRWBLK
              IF(IPVT.EQ.IPLUSN)GO TO 340
                 DO 315 L = I,NRWBLK
                    SWAP = ARRAY(L,IPVT,K)
                    ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
                    ARRAY(L,IPLUSN,K) = SWAP
315              CONTINUE
                 IPVBLK = IPVT-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 330
                    DO 320 L = 1,NRWBLK
                       SWAP = ARRAY(L,IPVBLK,KPLUS1)
                       ARRAY(L,IPVBLK,KPLUS1)
     *                                 = ARRAY(L,IRWBLK,KPLUS1)
                       ARRAY(L,IRWBLK,KPLUS1) = SWAP
320                 CONTINUE
                    GO TO 340
330              CONTINUE
                 DO 335 L = 1,NRWBOT
                    SWAP = BOTBLK(L,IPVBLK)
                    BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
                    BOTBLK(L,IRWBLK) = SWAP
335              CONTINUE
340           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM ROW
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              COLPIV = ARRAY(I,IPLUSN,K)
             IF(I.EQ.NRWBLK)GO TO 342
                DO 341 L = IPLUS1,NRWBLK
                   ARRAY(L,IPLUSN,K) = ARRAY(L,IPLUSN,K)/COLPIV
341             CONTINUE
342          CONTINUE
             IF(K.EQ.NBLOKS)GO TO 344
                DO 343 L = 1,NRWBLK
                   ARRAY(L,IRWBLK,KPLUS1) = ARRAY(L,IRWBLK,KPLUS1)/
     *                                         COLPIV
343             CONTINUE
                GO TO 346
344          CONTINUE
                DO 345 L = 1,NRWBOT
                   BOTBLK(L,IRWBLK) = BOTBLK(L,IRWBLK)/COLPIV
345             CONTINUE
346          CONTINUE
              DO 380 J = LOOP,NCLBLK
                 COLMLT = ARRAY(I,J,K)
                 IF(I.EQ.NRWBLK)GO TO 350
                    DO 347 L = IPLUS1,NRWBLK
                       ARRAY(L,J,K) = ARRAY(L,J,K)
     *                                -COLMLT*ARRAY(L,IPLUSN,K)
347                 CONTINUE
350              CONTINUE
                 JRWBLK = J-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 370
                    DO 360 L = 1,NRWBLK
                       ARRAY(L,JRWBLK,KPLUS1) =
     *                                  ARRAY(L,JRWBLK,KPLUS1)
     *                         -COLMLT*ARRAY(L,IRWBLK,KPLUS1)
360                 CONTINUE
                    GO TO 380
370              CONTINUE
                 DO 375 L = 1,NRWBOT
                    BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)
     *                              -COLMLT*BOTBLK(L,IRWBLK)
375              CONTINUE
380           CONTINUE
390        CONTINUE
           INCR = INCR + NRWBLK
395     CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 500
           DO 470 J = NRWTP1,NVRLP0
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = ABS(BOTBLK(JMINN,J))
              LOOP = JMINN+1
              DO 410 I = LOOP,NRWBOT
                 TEMPIV = ABS(BOTBLK(I,J))
                 IF(TEMPIV.LE.ROWMAX) GO TO 410
                 IPVT = I
                 ROWMAX = TEMPIV
410           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = MAX(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 430
                 DO 420 L = J,NOVRLP
                    SWAP = BOTBLK(IPVT,L)
                    BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
                    BOTBLK(JMINN,L) = SWAP
420              CONTINUE
430           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = BOTBLK(JMINN,J)
              DO 440 I = LOOP,NRWBOT
                 BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
440           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 460 L = JPLUS1,NOVRLP
                 ROWMLT = BOTBLK(JMINN,L)
                 DO 450 I = LOOP,NRWBOT
                    BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
450              CONTINUE
460           CONTINUE
470        CONTINUE
500     CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
        IF(PIVMAX+ABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
1000    CONTINUE
        IFLAG = -1
        RETURN
        END

*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX CX(*),CY(*)
*       ..
*
*
*> Purpose:
*  =============
*>
*>
*>
*>    CCOPY copies a vector x to a vector y.
*>
*
*  Arguments:
*  ==========
*
*> N
*> 
*>         N is INTEGER
*>         number of elements in input vector(s)
*>
*> CX
*>          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*>
*> INCX
*>          INCX is INTEGER
*>         storage spacing between elements of CX
*>
*> CY
*>          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*>
*> INCY
*>          INCY is INTEGER
*>         storage spacing between elements of CY
*
*  Authors:
*  ========
*
*> Univ. of Tennessee
*> Univ. of California Berkeley
*> Univ. of Colorado Denver
*> NAG Ltd.
*
*> November 2017
*
*> complex_blas_level1
*
*> Further Details:
*  ================
*>
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*>
*  =====================================================================
      SUBROUTINE ccopy(N,CX,INCX,CY,INCY)
      implicit none
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX CX(*),CY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            cy(i) = cx(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            cy(iy) = cx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END
