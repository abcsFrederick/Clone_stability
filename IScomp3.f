      PROGRAM ISCOMP3
      PARAMETER (MAXSITE=20000,MAXIS=20000,MAXTIME=5)
      CHARACTER*300 LINE,BL300
      CHARACTER*80 INFILE1,OUTFILE,PREFIX
      CHARACTER*30 P3L(MAXSITE),USEL(MAXSITE),BL30,CH30
      CHARACTER*20 HEADER(MAXTIME),BL20,CH20
      CHARACTER*10 CH10A,CH10B,BL10
      CHARACTER*3 CH3A,CH3B,CH3C
      CHARACTER*1 BLNK,HTAB,CH1
      INTEGER HAVE(MAXSITE),ICLONE1(MAXIS),ICLONE2(MAXIS),IHT(10),
     *ICIN(10),LHEAD(MAXTIME)
      REAL COUNT(MAXSITE,2),USE(MAXSITE,2),WGHT(MAXSITE),WGHT0(MAXSITE),
     *AX(5),BX(5),YCA(5),VAR(MAXIS),DRELE(5)
      REAL*8 SEED,COUNTIN(100,MAXTIME),TOTAL(MAXTIME),
     *ROULET(100,MAXTIME)
      DATA BL300/' '/
      DATA BL10/'          '/
      DATA BL20/'                    '/
      DATA BL30/'                              '/
      DATA BLNK/' '/
C..This program uses a roulette wheel selection procedure to select sets
C....of 2000 insertion sites from the 100 largest clones and calculates
C....their differences. This is repeated 10000 times.
      CALL GETARG(1,INFILE1)
      CALL GETARG(2,PREFIX)
      CALL GETARG(3,CH1)
      CALL GETARG(4,CH20)
      READ(CH1,'(I1)') NCOL
      IF(NCOL.GT.MAXTIME) THEN
        WRITE(*,'(A)') 'ERROR...NEED TO INCREASE MAXTIME'
        STOP
      ENDIF
      READ(CH20,*) SEED
      NPICK=100
      HTAB=CHAR(9)
      LENGTH=2020
C..READ INFILE1
      OPEN(UNIT=1,FILE=INFILE1,FORM='FORMATTED')
      REWIND(1)
      READ(1,'(A)') LINE
      IHT(1)=0
      DO 1 I=1,NCOL
      IHT(I+1)=INDEX(LINE,HTAB)
      IF(IHT(I+1).EQ.0) IHT(I+1)=INDEX(LINE,BL10)
      LINE(IHT(I+1):IHT(I+1))=BLNK
    1 CONTINUE
      DO 101 I=1,NCOL
      HEADER(I)=BL10
      LHEAD(I)=IHT(I+1)-IHT(I)-1
      HEADER(I)(1:LHEAD(I))=LINE(IHT(I)+1:IHT(I+1)-1)
  101 CONTINUE
      DO 103 I=1,100
      READ(1,*)(COUNTIN(I,J),J=1,NCOL)
  103 CONTINUE
      CLOSE(1)
C..BUILD THE NCOL ROULETTE WHEELS
      DO 109 J=1,NCOL
      TOTAL(J)=0.0D0
      DO 105 I=1,100
      TOTAL(J)=TOTAL(J)+COUNTIN(I,J)
  105 CONTINUE
      TOTAL(J)=1.0D0/TOTAL(J)
      ROULET(1,J)=COUNTIN(1,J)*TOTAL(J)
      DO 107 I=2,100
      ROULET(I,J)=ROULET(I-1,J)+COUNTIN(I,J)*TOTAL(J)
  107 CONTINUE
      ROULET(100,J)=1.0D0
  109 CONTINUE
C..EXAMINE THE GROUPS ONE AT A TIME AND PAIRWISE
      LUKE=ICHAR('A')-1
      DO 148 JC=1,NCOL
      DO 146 IC=JC,NCOL
C
      LPRE=INDEX(PREFIX,BLNK)
      OUTFILE=PREFIX
      LUKE=LUKE+1
      OUTFILE(LPRE:LPRE)=CHAR(LUKE)
      OUTFILE(LPRE+1:LPRE+4)='.txt'
      OPEN(UNIT=2,FILE=OUTFILE,FORM='FORMATTED')
      REWIND(2)
      IJ=INDEX(INFILE1,BLNK)-1
      WRITE(2,'(6A)') 'Comparing ',HEADER(JC)(1:LHEAD(JC)),
     *' versus ',HEADER(IC)(1:LHEAD(IC)),
     *' from ',INFILE1(1:IJ)
      WRITE(2,'(A,I3,A)') 'Selecting top ',NPICK,' clones from each.'
      WRITE(2,'(A,F12.0)') 'SEED = ',SEED
C..START RUNNING THE SIMULATIONS
      WRITE(2,'(/A3,4(A1,F3.1))') '0.0',(HTAB,0.5*FLOAT(I),I=1,4)
      DO 44 IPERM=1,10000
      DO 14 J=1,2
      DO 12 I=1,100
      COUNT(I,J)=0.0
   12 CONTINUE
   14 CONTINUE
      DO 16 I=1,100
      HAVE(I)=0
   16 CONTINUE
C..FILL IN COUNT(I,1)
      CALL SURAND(SEED,LENGTH,VAR)
      DO 113 I=1,2000
      DO 111 K=1,100
      IF(ROULET(K,JC).GE.VAR(I)) THEN
        COUNT(K,1)=COUNT(K,1)+1
        GOTO 113
      ENDIF
  111 CONTINUE
      WRITE(*,'(A)') 'Problem with COUNT(I,1)'
      WRITE(*,'(10F8.5)')(ROULET(K,JC),K=1,100)
      WRITE(*,'(/F8.5)') VAR(I)
      COUNT(100,1)=COUNT(100,1)+1
  113 CONTINUE
C..FILL IN COUNT(I,2)
      CALL SURAND(SEED,LENGTH,VAR)
      DO 117 I=1,2000
      DO 115 K=1,100
      IF(ROULET(K,IC).GE.VAR(I)) THEN
        COUNT(K,2)=COUNT(K,2)+1
        GOTO 117
      ENDIF
  115 CONTINUE
      WRITE(*,'(A)') 'Problem with COUNT(I,2)'
      WRITE(*,'(10F8.5)')(ROULET(K,IC),K=1,100)
      WRITE(*,'(/F8.5)') VAR(I)
      COUNT(100,2)=COUNT(100,2)+1
  117 CONTINUE
C..SELECT THE NPICK CLONES FOR COMPARISON
      DO 24 IP=1,NPICK
      CM=-1.0
      IM=0
      DO 22 I=1,100
      IF(COUNT(I,1).GT.CM.AND.HAVE(I).EQ.0) THEN
        CM=COUNT(I,1)
        IM=I
      ENDIF
   22 CONTINUE
      USE(IP,1)=COUNT(IM,1)
      USE(IP,2)=COUNT(IM,2)
      HAVE(IM)=1
      IF(IP.EQ.NPICK) CL=COUNT(IM,1)
   24 CONTINUE
      NUSE=NPICK
C..LOOK FOR TIES NOT INCLUDED SO FAR
      DO 26 I=1,100
      IF(HAVE(I).EQ.0.AND.COUNT(I,1).EQ.CL) THEN
        NUSE=NUSE+1
        USEL(NUSE)=P3L(I)
        USE(NUSE,1)=COUNT(I,1)
        USE(NUSE,2)=COUNT(I,2)
        HAVE(I)=1
      ENDIF
   26 CONTINUE
C..LOOK FOR LARGEST CLONES FROM INFILE2
      DO 30 IP=1,NPICK
      CM=-1.0
      IM=0
      DO 28 I=1,100
      IF(HAVE(I).NE.2.AND.COUNT(I,2).GT.CM) THEN
        CM=COUNT(I,2)
        IM=I
      ENDIF
   28 CONTINUE
      IF(HAVE(IM).EQ.0) THEN
        NUSE=NUSE+1
        USEL(NUSE)=P3L(IM)
        USE(NUSE,1)=COUNT(IM,1)
        USE(NUSE,2)=COUNT(IM,2)
      ENDIF
      HAVE(IM)=2
      IF(IP.EQ.NPICK) CL=COUNT(IM,2)
   30 CONTINUE
C..LOOK FOR TIES WITH CL
      NA2=0
      DO 32 I=1,100
      IF(HAVE(I).NE.2.AND.COUNT(I,2).EQ.CL) THEN
        NA2=NA2+1
        IF(HAVE(I).EQ.0) THEN
          NUSE=NUSE+1
          USEL(NUSE)=P3L(I)
          USE(NUSE,1)=COUNT(I,1)
          USE(NUSE,2)=COUNT(I,2)
        ENDIF
        HAVE(I)=2
      ENDIF
   32 CONTINUE
C      WRITE(2,'(A,I4)') '#Clones from Set 2 is ',NPICK+NA2
C      WRITE(2,'(A,I4,A)') 'There are ',NUSE,' clones overall'
C..BUILD WGHT0(NUSE)
      DO 133 I=1,NUSE
      WGHT0(I)=0.5*(USE(I,1)+USE(I,2))
      IF(WGHT0(I).EQ.0.0) WGHT0(I)=1.0
  133 CONTINUE
C..START THE NONLINEAR FITS
      DO 42 IZ=1,5
      Z=0.5*FLOAT(IZ-1)
      DO 36 I=1,NUSE
      WGHT(I)=WGHT0(I)**Z
   36 CONTINUE
      W=0.0
      WX=0.0
      WY=0.0
      WXX=0.0
      WXY=0.0
      DO 38 I=1,NUSE
      W=W+WGHT(I)
      WX=WX+WGHT(I)*USE(I,1)
      WY=WY+WGHT(I)*USE(I,2)
      WXX=WXX+WGHT(I)*USE(I,1)*USE(I,1)
      WXY=WXY+WGHT(I)*USE(I,1)*USE(I,2)
   38 CONTINUE
      B=(WX*WY-W*WXY)/(WX*WX-W*WXX)
      A=(WY-B*WX)/W
      AX(IZ)=A
      BX(IZ)=B
      RMSE=0.0
      AAE=0.0
      RELE=0.0
      DEN=1.0/(B+(1.0/B))
      DO 40 I=1,NUSE
      XDP=(USE(I,2)+(USE(I,1)/B)-A)*DEN
      YDP=B*XDP+A
      DIS=SQRT((XDP-USE(I,1))**2+(YDP-USE(I,2))**2)
      RMSE=RMSE+DIS**2
      AAE=AAE+DIS
      RELE=RELE+DIS**2/WGHT0(I)
   40 CONTINUE
      RMSE=SQRT(RMSE/FLOAT(NUSE))
      AAE=AAE/FLOAT(NUSE)
      RELE=RELE/FLOAT(NUSE)
      DRELE(IZ)=RELE
   42 CONTINUE
      WRITE(2,'(F10.5,4(A1,F10.5))') DRELE(1),(HTAB,DRELE(I),I=2,5)
C..DONE WITH A PERMUTATION
   44 CONTINUE
      CLOSE(2)
  146 CONTINUE
  148 CONTINUE
      END
      SUBROUTINE SURAND(SEED,N,X)
      REAL*8 SEED,A,XM,XMI,ONE
      DIMENSION X(*)
      DATA A,XM,ONE/16807.0D0,2147483647.0D0,1.0D0/
      IF(N.LT.1) RETURN
      XMI=ONE/XM
      DO 10 I=1,N
      SEED=A*SEED
      IX=INT(SEED*XMI)
      SEED=SEED-FLOAT(IX)*XM
      X(I)=SEED*XMI
   10 CONTINUE
      RETURN
      END
