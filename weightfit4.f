      PROGRAM WEIGHTFIT4
      PARAMETER (MAXSITE=10000)
      CHARACTER*80 LINE,INFILE,OUTFILE
      CHARACTER*30 P3L(MAXSITE),USEL(MAXSITE),BL30,CH30,SET1,SET2
      CHARACTER*3 CH3
      CHARACTER*1 BLNK,HTAB,CH1
      INTEGER HAVE(MAXSITE),USEI(MAXSITE,2),LSET1,LSET2
      REAL COUNT(MAXSITE,2),USE(MAXSITE,2),WGHT(MAXSITE),WGHT0(MAXSITE),
     *AX(5),BX(5),YCA(5)
      DATA BL30/'                              '/
      DATA BLNK/' '/
C..THIS PROGRAM DIFFERS FROM WEIGHTFIT IN THAT THERE IS A SINGLE INPUT FILE
C..THAT CONTAINS pseudo3LTR, Set1, Set2
      CALL GETARG(1,INFILE)
      CALL GETARG(2,OUTFILE)
      CALL GETARG(3,CH3)
      READ(CH3,*) NPICK
      HTAB=CHAR(9)
C
      OPEN(UNIT=2,FILE=OUTFILE,FORM='FORMATTED')
      REWIND(2)
      IJ=INDEX(INFILE,BLNK)-1
C      WRITE(2,'(2A)') 'Comparing data in ',INFILE(1:IJ)
C..READ INFILE1
      OPEN(UNIT=1,FILE=INFILE,FORM='FORMATTED')
      REWIND(1)
      NHAVE=0
      READ(1,'(A)') LINE
C      WRITE(*,'(A80)') LINE
      L1=INDEX(LINE,HTAB)
      LINE(L1:L1)=BLNK
      L2=INDEX(LINE,HTAB)
      LINE(L2:L2)=BLNK
      L3=INDEX(LINE,HTAB)
      IF(L3.EQ.0) L3=INDEX(LINE,BL30(1:10))
      WRITE(2,'(6A)') 'Comparing ',LINE(L1+1:L2-1),' versus ',
     *LINE(L2+1:L3-1),' from ',INFILE(1:IJ)
      SET1=BL30
      LSET1=L2-L1-1
      SET1(1:LSET1)=LINE(L1+1:L2-1)
      SET2=BL30
      LSET2=L3-L2-1
      SET2(1:LSET2)=LINE(L2+1:L3-1)
      WRITE(2,'(A,I3,A)') 'Selecting top ',NPICK,' clones from each.'
    2 READ(1,'(A)',ERR=4,END=4) LINE
      IJ=INDEX(LINE,HTAB)
      IF(IJ.LT.2) GOTO 4
      CH30=BL30
      CH30(1:IJ-1)=LINE(1:IJ-1)
      LINE(1:IJ)=BL30(1:IJ)
      READ(LINE,*) IC1,IC2
      IF(IC1.EQ.0.AND.IC2.EQ.0) GOTO 2
      NHAVE=NHAVE+1
      P3L(NHAVE)=CH30
      COUNT(NHAVE,1)=IC1
      COUNT(NHAVE,2)=IC2
      HAVE(NHAVE)=0
      GOTO 2
    4 CLOSE(1)
C..LOAD DATA FROM SET1 TO USE(NUSE,1) AND USEL(NUSE)
      DO 14 IP=1,NPICK
      CM=0.0
      IM=0
      DO 12 I=1,NHAVE
      IF(COUNT(I,1).GT.CM.AND.HAVE(I).EQ.0) THEN
        CM=COUNT(I,1)
        IM=I
      ENDIF
   12 CONTINUE
      USEL(IP)=P3L(IM)
      USE(IP,1)=COUNT(IM,1)
      USE(IP,2)=COUNT(IM,2)
      HAVE(IM)=1
      IF(IP.EQ.NPICK) CL=COUNT(IM,1)
   14 CONTINUE
      NUSE=NPICK
C..LOOK FOR TIES NOT INCLUDED SO FAR
      DO 16 I=1,NHAVE
      IF(HAVE(I).EQ.0.AND.COUNT(I,1).EQ.CL) THEN
        NUSE=NUSE+1
        USEL(NUSE)=P3L(I)
        USE(NUSE,1)=COUNT(I,1)
        USE(NUSE,2)=COUNT(I,2)
        HAVE(I)=1
      ENDIF
   16 CONTINUE
      WRITE(2,'(3A,I4)') '#Clones from ',SET1(1:LSET1),' is ',NUSE
C..LOOK FOR LARGEST CLONES FROM SET2
      DO 20 IP=1,NPICK
      CM=0.0
      IM=0
      DO 18 I=1,NHAVE
      IF(HAVE(I).NE.2.AND.COUNT(I,2).GT.CM) THEN
        CM=COUNT(I,2)
        IM=I
      ENDIF
   18 CONTINUE
      IF(HAVE(IM).EQ.0) THEN
        NUSE=NUSE+1
        USEL(NUSE)=P3L(IM)
        USE(NUSE,1)=COUNT(IM,1)
        USE(NUSE,2)=COUNT(IM,2)
      ENDIF
      HAVE(IM)=2
      IF(IP.EQ.NPICK) CL=COUNT(IM,2)
   20 CONTINUE
C..LOOK FOR TIES WITH CL
      NA2=0
      DO 22 I=1,NHAVE
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
   22 CONTINUE
      WRITE(2,'(3A,I4)') '#Clones from ',SET2(1:LSET2),' is ',NPICK+NA2
      WRITE(2,'(A,I4,A)') 'There are ',NUSE,' clones overall'
C..BUILD WGHT0(NUSE)
C..FIRST SCALE ALL DATA SO THAT THE AVERAGE OF EACH SET IS 20.0
      TOT1=0.0
      TOT2=0.0
      DO 24 I=1,NUSE
      USEI(I,1)=NINT(USE(I,1))
      USEI(I,2)=NINT(USE(I,2))
      TOT1=TOT1+USE(I,1)
      TOT2=TOT2+USE(I,2)
   24 CONTINUE
      SC1=20.0*FLOAT(NUSE)/TOT1
      SC2=20.0*FLOAT(NUSE)/TOT2
      DO 101 I=1,NUSE
      USE(I,1)=SC1*USE(I,1)
      USE(I,2)=SC2*USE(I,2)
      WGHT0(I)=0.5*(USE(I,1)+USE(I,2))
  101 CONTINUE
C..START THE LINEAR FITS, UNWEIGHTED AND WEIGHTED
      WRITE(2,'(/99A)') 
     *'Weight',HTAB,'Intercept',HTAB,'Slope',HTAB,'RMSE',HTAB,'AAE',
     *HTAB,'RELE'
      DO 32 IZ=1,5
      Z=0.5*FLOAT(IZ-1)
      DO 26 I=1,NUSE
      WGHT(I)=WGHT0(I)**Z
   26 CONTINUE
      W=0.0
      WX=0.0
      WY=0.0
      WXX=0.0
      WXY=0.0
      DO 28 I=1,NUSE
      W=W+WGHT(I)
      WX=WX+WGHT(I)*USE(I,1)
      WY=WY+WGHT(I)*USE(I,2)
      WXX=WXX+WGHT(I)*USE(I,1)*USE(I,1)
      WXY=WXY+WGHT(I)*USE(I,1)*USE(I,2)
   28 CONTINUE
      B=(WX*WY-W*WXY)/(WX*WX-W*WXX)
      A=(WY-B*WX)/W
      AX(IZ)=A
      BX(IZ)=B
      RMSE=0.0
      AAE=0.0
      RELE=0.0
      DEN=1.0/(B+(1.0/B))
      DO 30 I=1,NUSE
      XDP=(USE(I,2)+(USE(I,1)/B)-A)*DEN
      YDP=B*XDP+A
      DIS=SQRT((XDP-USE(I,1))**2+(YDP-USE(I,2))**2)
      RMSE=RMSE+DIS**2
      AAE=AAE+DIS
      RELE=RELE+DIS**2/WGHT0(I)
   30 CONTINUE
      RMSE=SQRT(RMSE/FLOAT(NUSE))
      AAE=AAE/FLOAT(NUSE)
      RELE=RELE/FLOAT(NUSE)
      WRITE(2,'(F3.1,5(A1,F10.5))') Z,HTAB,A,HTAB,B,HTAB,RMSE,HTAB,
     *AAE,HTAB,RELE
   32 CONTINUE
C..WRITE OUT THE DATA
      WRITE(2,'(/99A)') 'pseudo3LRT',HTAB,SET1(1:LSET1),HTAB,
     *SET2(1:LSET2),
     *HTAB,'CALC2(0.0)',HTAB,'CALC2(0.5)',HTAB,'CALC2(1.0)',
     *HTAB,'CALC2(1.5)',HTAB,'CALC2(2.0)'
      DO 36 I=1,NUSE
      DO 34 IZ=1,5
      YCA(IZ)=AX(IZ)+BX(IZ)*USE(I,1)
      IF(YCA(IZ).LT.0.0) YCA(IZ)=0.0
   34 CONTINUE
      IJ=INDEX(USEL(I),BLNK)-1
      WRITE(2,'(A,7(A1,F7.2))') USEL(I)(1:IJ),HTAB,USE(I,1),HTAB,
     *USE(I,2),(HTAB,YCA(IZ),IZ=1,5)
   36 CONTINUE
C
      WRITE(2,'(//5A)') 'pseudo3LRT',HTAB,SET1(1:LSET1),HTAB,
     *SET2(1:LSET2)
      DO 38 I=1,NUSE
      IJ=INDEX(USEL(I),BLNK)-1
      WRITE(2,'(A,2(A1,I8))') USEL(I)(1:IJ),HTAB,USEI(I,1),HTAB,
     *USEI(I,2)
   38 CONTINUE
      CLOSE(2)
      END
