PROGRAM Domcalc
************************************************************************
* Magnetic antidot simulation - XY Model v.6.1 *
* --------------------------------------------- *
* This code is an variation of the XY model, metropolis algorithm. *
* It simulates the time-evolution of an array of variable spin *
* elements at various temperatures with periodic antidots. *
************************************************************************
* This version includes the following features:
* Energy & magnetization effects to the first nearest neighbor
* Susceptibility, heat capacity & 4th order cumulant
* Demagnetization, external field & anisotropy are enabled
* Magnetization statistics of subdomains
* flip frequency monitoring
************************************************************************
************************************************************************
* Copyright (C) 2006 Brian Weir
*
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the
* Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the
* Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
* 02111-1307 USA
************************************************************************
Define initial variables
INTEGER ITRNUM, L, K, DLOOP, NEG, I, J, RNGSED, N, SINC
INTEGER STPNUM, Q, T, TP, IDX, MSNAP, LX, LY
C Define array variables
REAL*8 EX, EY, SX(10000), SY(10000), NNX, NNY
C Define calculation variables
REAL*4 RAND, TMPSTP
REAL*8 TMPKTJ, OP, ENET, M, PROB, ENETT
REAL*8 PI, TWOPI
CHARACTER*12 FNAME, LVAR
LOGICAL MFLIP, TMPFLIP, CTFLIP
DIMENSION ENET(6000000), OP(6000000)
C Define variables for susceptibility, heat capacity and cumulant
56
REAL*8 X, C, SUMOP, SUMMAG, SUME, SQE, SQOP, U, SQQOP
REAL*8 ENSUM, MOPSUM, CEN, COP, SQRNOP, RANOP, SUMX, SQX
REAL*8 MX, MY, SX2, SY2, EI, DE, E
REAL*8 SQENET(6000000)
C Define quantum anti-dot properties
INTEGER DOTX, DOTY, SPACX, SPACY
INTEGER NN, SN, NT
COMMON /DOT/ DOTX, DOTY, SPACX, SPACY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX, SY, SN(100000)
COMMON /EN/ ENET, SQENET, K
C Define complex Hamiltonian variables
REAL*8 HEXTX, HEXTY, D, EISO, ISODIR, HANGL, HEXT
COMMON /VAR/ EISO, HEXTX, HEXTY
C Extra info variables
REAL*4 RND1, RND2
REAL*8 SX1, SY1
INTEGER CHANGE(100000)
REAL*8 MBIN1X(110), MBIN1Y(110), MBIN2X(110)
REAL*8 MBIN2Y(110), MBIN3X(110), MBIN3Y(110)
REAL*8 MBIN1(110), MBIN2(110), MBIN3(110)
REAL*8 X1, X2, X3, BINSIZE(110)
COMMON /DOM1/ BINSIZE,MBIN1,MBIN1X,MBIN1Y,MBIN2,MBIN2X,MBIN2Y,
% MBIN3,MBIN3X,MBIN3Y,X1,X2,X3,TMPKTJ, ITRNUM, NT
C Define demagnetization calculation variables
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMAG, HDMAG
REAL*8 DMDISTH(20,2), DMDISTV(20,2)
INTEGER IDM
COMMON /DMG/ DMA, DMB, DMC, DMD, DME,
% MXA, MXB, MXC, MXD, MXE,
% MYA, MYB, MYC, MYD, MYE
COMMON /DMG2/ DMAG, HDMAG, IDM
COMMON /DIST/ DMDISTH, DMDISTV
C More variables for spin pictures
INTEGER NPASSPICT, ANS5, IP, ANS4,ANS3
CHARACTER*12 FSPINPICT
c CHANGE tracks how often a site flips
57
C Ask for user defined variables
C Define the array
PRINT *, 'Enter array length: '
READ (*,*) LX
PRINT *, 'Enter array width: '
READ (*,*) LY
C set initial values
ANS5 = 0
ANS4=0
ANS3=0
IP = 0
ITRNUM = 1
N = LX*LY
NT = N
ENETT = -2*N
TMPNUM = 20
PI = 3.141592653589793d0
TWOPI = 2.0D0*PI
HEXT = 0
EISO = 0
DMAG = 0
IDM = 1
fspinpict = 'spins.txt'
C Initialize the entire array to initial energy, +1
DO 10 I = 1, N
SX(I) = 1
SY(I) = 0
SN(I) = 1
CHANGE(I) = 0
C NOTE: I is used for columns while J is used for rows
10 CONTINUE
c Initialize the bins for unit cell statistics
DO I = 1, 102
MBIN1X(I) = 0
MBIN1Y(I) = 0
MBIN1(I) = 0
MBIN2X(I) = 0
MBIN2Y(I) = 0
MBIN2(I) = 0
MBIN3X(I) = 0
MBIN3Y(I) = 0
MBIN3(I) = 0
c Store statistics in 'bins'
c bins range from '-1' to '1'
58
c each bin is magnitude '0.02' wide, 100 total
BINSIZE(I) = -1 + 0.02*(I-1)
END DO
C Call NEARBR to initialize the nearest neighbor array
CALL NEARBR
C Ask for file
PRINT *, 'Name the file to record to'
READ (*,*) FNAME
CLOSE (UNIT = 1)
C Ask for multiple temperatures
PRINT *, 'Run a range of temperatures? (Y/N)'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
TMPFLIP = .TRUE.
PRINT *, 'Enter initial temperature in kT/J'
READ (*,*) TMPKTJ
PRINT *, 'Enter temperature step size'
READ (*,*) TMPSTP
PRINT *, 'How many steps?'
READ (*,*) STPNUM
ELSE
TMPFLIP = .FALSE.
PRINT *, 'Enter temperature in kT/J'
READ (*,*) TMPKTJ
TMPSTP = 0
STPNUM = 1
END IF
C Ask for anti-dot generation
PRINT *, 'Use anti-dots? (Y/N)'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
PRINT *, 'Use square dots?'
READ (*,*) LVAR
PRINT *, 'Dot diameter(X)?'
READ (*,*) DOTX
IF (LVAR .EQ. 'y') THEN
DOTY = DOTX
ELSE
PRINT *, 'Dot height(y)?'
READ (*,*) DOTY
END IF
PRINT *, 'Dot X spacing?'
59
READ (*,*) SPACX
PRINT *, 'Dot Y spacing?'
READ (*,*) SPACY
CALL DOTGEN
NT = 0
ENETT = 0
DO I = 1, N
NT = NT + SN(I)
ENETT = ENETT - SN(I)*(SN(NN(I, 1)) + SN(NN(I, 2)))
END DO
END IF
E = ENETT
EI = ENETT
PRINT *, 'Enter number of passes'
READ (*,*) ITRNUM
PRINT *, 'Enter RNG seed'
READ (*,*) RNGSED
C This allows the system to reach some equilibrium before recording data
PRINT *, 'Skip how many passes before recording?'
READ (*,*) DLOOP
C Ask for complex Hamiltonian variables
PRINT *, 'Use anisotropy? (y/n)'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
PRINT *, 'Anisotropy energy?'
READ (*,*) EISO
C PRINT *, 'Direction of the anisotropy? (+1/-1)'
C READ (*,*) ISODIR
ELSE
EISO = 0
END IF
PRINT *, 'Use external field? (y/n)'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
PRINT *, 'Enter the magnitude of the field'
READ (*,*) HEXT
PRINT *, 'Enter the angle of the field'
READ (*,*) HANGL
HANGL = (HANGL * PI) / 180
HEXTX = HEXT*DCOS(HANGL)
HEXTY = HEXT*DSIN(HANGL)
ELSE
60
HEXT = 0
END IF
PRINT *, 'Include demagnitization? (y/n)'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
PRINT *, 'Enter the magnitude of the field'
READ (*,*) DMAG
C Call the remaining 4 nearest neighbor arrays
CALL NEARBR2
CALL NEARBR3
CALL NEARBR4
CALL NEARBR5
C precalculate neighbor distances
CALL DMDIST
ELSE
DMAG = 0
END IF
write(*,*) 'Print out spin pictures?'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
ANS5 = 1
write(*,*) 'Enter lattice pass number for pictures:'
read(*,*) npasspict
IP = NPASSPICT
c write(*,*) 'Enter file name for pictures:'
c read(*,*) fspinpict
fspinpict = 'spins.txt'
else
ANS5 = 0
endif
write(*,*) 'Record flip-rate statistics?'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
ANS4 = 1
else
ANS4=0
end if
write(*,*) 'Record domain data?'
READ (*,*) LVAR
IF (LVAR .EQ. 'y') THEN
ANS3 = 1
else
ANS3=0
end if
61
************************************************************************
* Write additional information to file:
OPEN (UNIT = 1, FILE = 'data.txt')
C Write parameters
WRITE (1, *) 'Lattice size:', LX, 'x', LY
WRITE (1, *) 'Number of passes:', ITRNUM
WRITE (1, *) 'Number of skipped passes:', DLOOP
WRITE (1, *) 'Number generator seed:', RNGSED
C Write dot parameters
WRITE (1, *) 'Antidot geometry:', DOTX, DOTY, SPACX, SPACY
C Write variables
WRITE (1, *) 'Anisotropy Magnitude:', EISO
WRITE (1, *) 'External Field Magnitude:', HEXT
WRITE (1, *) 'External Field Angle:', HANGL
WRITE (1, *) 'Demagnitization Magnitude:', DMAG
C Write output headings
close(unit=1,status='save')
OPEN (UNIT = 1, FILE = FNAME)
WRITE (1, 15) 'TMPKTJ', 'OP', 'X', 'M', 'MX', 'MY'
15 FORMAT(6(5X, A6))
c CLOSE(UNIT = 1)
C Prepare domain data
IF(ans3 .eq.1) then
close(unit=1,status='save')
OPEN (UNIT = 1, FILE = 'domain.txt')
WRITE (1, 16) 'Bin', 'OP1', 'M1X', 'M1Y', 'OP2', 'M2X',
% 'M2Y', 'OP3', 'M3X', 'M3Y', 'X1', 'X2', 'X3'
16 FORMAT(13(5X, A6))
close(unit=1,status='save')
open(unit=1,file=fname, STATUS='OLD')
end if
c Prepare spin statistics data
IF(ans3 .eq.1) then
close(unit=1,status='save')
OPEN (UNIT = 1, FILE = fspinpict)
write(1,17) 'i,','sn(i),','sx(i),','sy(i),'
17 FORMAT(10(5X, A6))
close(unit=1,status='save')
OPEN (UNIT = 1, FILE = FNAME, STATUS = 'OLD')
endif
********************
C RAN2 requires a negative seed, so we'll just correct it instead of asking
RNGSED = -ABS(RNGSED)
********************
PRINT *, 'Randomizing system, please wait...'
62
C Here's a loop to start up with DLOOP skipped passes
DO K=1, DLOOP
DO 30 I=1,N
IDM = I
SX1 = SX(I)
SY1 = SY(I)
C Make sure holes are not counted
IF (SN(I) .EQ. 1) THEN
NNX = 0
NNY = 0
DO J=1,4
NNX = NNX + SX(NN(I,J))
NNY = NNY + SY(NN(I,J))
END DO
CALL RAN2(RNGSED, RAND)
RND1 = RAND
SX2 = DCOS(RND1*TWOPI)
SY2 = DSIN(RND1*TWOPI)
EX= SX2*NNX
EY= SY2*NNY
CALL DMAGC
EI = -(NNX*SX(I) + NNY*SY(I))
% -EISO * (SX(I)**2)
% -HEXTX*SX(I) -HEXTY*SY(I)
% -HDMAG
SX(I) = SX2
SY(I) = SY2
CALL DMAGC
E = -(EX + EY)
% -EISO * (SX(I)**2)
% -HEXTX*SX(I)- HEXTY*SY(I)
% -HDMAG
DE = E - EI
IF (DE .LT. 0) THEN
SY(I) = SY2
SX(I) = SX2
ENETT = ENETT + DE
ELSE
CALL RAN2(RNGSED, RAND)
63
IF (RAND .LE. EXP(-DE / TMPKTJ)) THEN
SY(I) = SY2
SX(I) = SX2
ENETT = ENETT + DE
ELSE
SX(I) = SX1
SY(I) = SY1
END IF
END IF
END IF
30 CONTINUE
ENDDO
********************
C We begin our temperature outer loop now
DO Q = 1, STPNUM
PRINT *, 'Evaluating temperature', TMPKTJ
********************
C Reset sums before starting next temperature
SUMMAG = 0
SUMOP = 0
SUME = 0
SQOP = 0
SQE = 0
SQQOP = 0
********************
C Now we evaluate the system ITRNUM times
DO K = 1, ITRNUM
C K is a placeholder variable
C Here begins the logic loop for evaluating each item
C we will use row-wise processing because its easier for neighbors
C We include nearest-neighbors here, but
DO 20 I=1,N
IDM = I
SX1 = SX(I)
SY1 = SY(I)
C Make sure holes are not counted
IF (SN(I) .EQ. 1) THEN
NNX = 0
NNY = 0
DO J=1,4
64
NNX = NNX + SX(NN(I,J))
NNY = NNY + SY(NN(I,J))
END DO
C Now, we evaluate the energy
CALL RAN2(RNGSED, RAND)
RND1 = RAND
SX2 = DCOS(RND1*TWOPI)
SY2 = DSIN(RND1*TWOPI)
EX= SX2*NNX
EY= SY2*NNY
CALL DMAGC
EI = -(NNX*SX(I) + NNY*SY(I))
% -EISO * (SX(I)**2)
% -HEXTX*SX(I)- HEXTY*SY(I)
% -HDMAG
SX(I) = SX2
SY(I) = SY2
CALL DMAGC
E = -(EX + EY)
% -EISO * (SX(I)**2)
% -HEXTX*SX(I)- HEXTY*SY(I)
% -HDMAG
DE = E - EI
C call random number here. RAND is a number from 0->1
C for high temperatures, S is more likely to flip
C allow for temperature effects here
IF (DE .LT. 0) THEN
SY(I) = SY2
SX(I) = SX2
ENETT = ENETT + DE
CHANGE(I) = CHANGE(I) + 1
ELSE
CALL RAN2(RNGSED, RAND)
IF (RAND .LE. EXP(-DE / TMPKTJ)) THEN
SY(I) = SY2
SX(I) = SX2
ENETT = ENETT + DE
CHANGE(I) = CHANGE(I) + 1
ELSE
SX(I) = SX1
SY(I) = SY1
END IF
END IF
END IF
20 CONTINUE
65
********************
C Resolve Magnetization, M = ( POS - NEG ) / MSIZE**2
M = 0
MX = 0
MY = 0
DO I = 1, N
MX = MX + SX(I)
MY = MY + SY(I)
END DO
M = SQRT(MX**2 + MY**2)
M = M / NT
C The order parameter is the absolute value of the magnetization
OP(K) = ABS(M)
C Set the Energy and Order Parameter averages now
SUMMAG = SUMMAG + M
SUMOP = SUMOP + OP(K)
SQOP = SQOP + OP(K)*OP(K)
SQQOP = SQQOP + OP(K)**4
c Recalculate energies for the full lattice
E = 0
SUME = 0
SQE = 0
ENET(K)=0
SQENET(K)=0
DO I=1,N
IDM = I
NNX = 0
NNY = 0
DO J=1,2
NNX = NNX + SX(NN(I,J))
NNY = NNY + SY(NN(I,J))
END DO
CALL DMAGC
E = -(NNX*SX(I) + NNY*SY(I))
% -EISO * (SX(I)**2)
% -(HEXTX*SX(I) + HEXTY*SY(I))
% -HDMAG
ENET(K) = ENET(K) + E
SQENET(K) = SQENET(K) + E*E
c this ends the energy recalculation here
end do
c take spin picture data if necessary
66
if(ans5.eq.1) then
if((K/NPASSPICT).eq.0) then
close(unit=1,status='save')
open(unit=1,file=fspinpict, STATUS='OLD')
do i=1,N
write(1,42) i,sn(i),sx(i),sy(i)
42 format(2x,i5,',',i2,',',f9.5,',',f9.5)
enddo
close(unit=1,status='save')
open(unit=1,file=fname, STATUS='OLD')
endif
endif
IF(ANS3.EQ.1) THEN
CALL DOMAIN
end if
C this ends the iteration loop here, K
ENDDO
DO K=1, ITRNUM
SUME = SUME+ENET(K)
SQE = SQE+ENET(K)*ENET(K)
END DO
SUME = SUME/ITRNUM
SQE = SQE/ITRNUM
********************
C Calculate the susceptibility, heat capacity and cumulant
X = (NT / TMPKTJ) * ( (SQOP/ITRNUM) - (SUMOP/ITRNUM)**2 )
C = ( (SQE) - (SUME)**2 ) / (NT * TMPKTJ * TMPKTJ )
U = 1 - ( (SQQOP/ITRNUM) / (3 * ((SQOP/ITRNUM)**2)) )
C Store data per temperature
OPEN (UNIT = 1, FILE = FNAME, STATUS='OLD')
WRITE (1, 70) TMPKTJ, (SUMOP/ITRNUM), X, M, MX, MY
70 FORMAT(6(2X, E12.5))
C Print domain data
IF(ANS3.EQ.1) THEN
close(unit=1,status='save')
OPEN (UNIT = 1, FILE = 'domain.txt', STATUS='OLD')
DO I =1,101
WRITE (1, 72) BINSIZE(I), ABS(MBIN1(I)), MBIN1X(I), MBIN1Y(I),
% ABS(MBIN2(I)), MBIN2X(I), MBIN2Y(I), ABS(MBIN3(I)),
% MBIN3X(I), MBIN3Y(I), X1, X2, X3
72 FORMAT(13(2X, E12.5))
67
END DO
close(unit=1,status='save')
open(unit=1,file=fname, STATUS='OLD')
END IF
C End temperature loop
TMPKTJ = TMPKTJ + TMPSTP
END DO
C Add rate-of-change data
if(ANS4.EQ.1)THEN
close(unit=1,status='save')
OPEN (UNIT = 1, FILE = 'roc.txt', STATUS='OLD')
DO I = 1, LY
WRITE(1,80) (CHANGE(LX*(I-1)+J),J=1,LX)
80 FORMAT(2X, 40(I7,1X))
END DO
CLOSE(UNIT = 1, STATUS='SAVE')
end if
C End Program
END
************************************************************************
* NEARBR: a nearest-neighbor indexing array subroutine *
************************************************************************
subroutine NEARBR
C Sitenum is the lattice location, 1-N
C nnidx is the nn location, 1-4 clockwise
INTEGER I, J, WX, WY, NR
COMMON /NN/ WX, WY, NR(10000, 4)
C Here, W corresponds to L, and NR() to NN()
DO I = 1, WX*WY
C Set locations of nearest neighbors
NR(I, 1) = I - WX
NR(I, 2) = I + 1
NR(I, 3) = I + WX
NR(I, 4) = I - 1
C Correct if an edge
C right side
IF (MOD(I,WX) .EQ. 0) THEN
NR(I, 2) = I - WX + 1
C left side
ELSE IF (MOD((I-1),WX) .EQ. 0) THEN
NR(I, 4) = I + WX - 1
END IF
C Top
68
IF (I .LE. WX) THEN
NR(I, 1) = I + WX*(WY-1)
C Bottom
ELSE IF (I .GT. WX*(WY-1)) THEN
NR(I, 3) = I - WX*(WY-1)
END IF
END DO
RETURN
END
************************************************************************
* DOTGEN: an anti-dot generator *
************************************************************************
SUBROUTINE DOTGEN
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
REAL*8 SX, SY
INTEGER SN
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
VAR = 0
NUMX = INT( LX / (X + SPX) )
NUMY = INT( LY / (Y + SPY) )
PRINT *, "Generating holes..."
C All dots begin in the upper left corner
C Take each point in the original dot and copy it at regular intervals
DO A = 1, Y
C Repeat each y point in the first dot NUMY times
DO B = 1, NUMY
C Now, repeat each x point NUMX times
DO C = 1, X
DO D = 1, NUMX
VAR = (A + (B-1)*(Y+SPY) -1)*LX + C + (X+SPX)*(D-
1)
SN(VAR) = 0
SX(VAR) = 0
SY(VAR) = 0
END DO
END DO
69
END DO
END DO
END
************************************************************************
* RAN2: a random number generator that actually works *
************************************************************************
Subroutine RAN2(idum,rn)
c from Press in Comp. Phys. v6 522 (1992)
Parameter (IM1=2147483563, IM2=2147483399, AM=1./IM1,
& IMM1=IM1-1, IA1=40014, IA2=40692, IQ1=53668, IQ2=52774,
& IR1=12211, IR2=3791, NTAB=32, NDIV=1+IMM1/NTAB,
& EPS=1.2e-7, RNMX=1.-EPS)
Integer idum2, j, k, iv(NTAB), iy
Save iv, iy, idum2
Data idum2/123456789/, iv/NTAB*0/, iy/0/
if(idum.le.0) then !initialize
idum=max(-idum,1)
idum2=idum
do j=NTAB+8,1,-1
k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if(idum.lt.0) idum=idum+IM1
if(j.le.NTAB) iv(j)=idum
enddo
iy=iv(1)
endif
k=idum/IQ1 !start here when not initializing
idum=IA1*(idum-k*IQ1)-k*IR1
if(idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if(idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if(iy.lt.1) iy=iy+IMM1
rn=min(AM*iy,RNMX)
return
end
************************************************************************
70
* DOMAIN: A domain magnetization tracking subroutine *
* *
* ___________ *
* | 4 | 1 | *
* |____|____| *
* | 3 | 2 | *
* |____|____| *
* *
************************************************************************
SUBROUTINE DOMAIN
INTEGER X, Y, SPX, SPY
INTEGER I, J, K, A, B
INTEGER LX, LY, NN
REAL*8 SX, SY
REAL*8 MX1, MY1, MX2, MY2, MX3, MY3
REAL*8 M1, M2, M3
REAL*8 MBIN1X(110), MBIN1Y(110), MBIN2X(110)
REAL*8 MBIN2Y(110), MBIN3X(110), MBIN3Y(110)
REAL*8 MBIN1(110), MBIN2(110), MBIN3(110)
REAL*8 BINSIZE(110)
INTEGER SN, NT
REAL*8 SQOP1, SQOP2, SQOP3
REAL*8 X1, X2, X3, TMPKTJ, S1, S2, S3
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DOM1/ BINSIZE,MBIN1,MBIN1X,MBIN1Y,MBIN2,MBIN2X,MBIN2Y,
% MBIN3,MBIN3X,MBIN3Y, X1, X2, X3, TMPKTJ, ITRNUM, NT
A = 0
B = 0
MX1 = 0
MY1 = 0
MX2 = 0
MY2 = 0
MX3 = 0
MY3 = 0
S1 = 0
S2 = 0
S3 = 0
J = 1
X1 = 0
X2 = 0
X3 = 0
K = 0
71
C outer loop checks each row
do I=1, LY
c A is the first location of the current row
A = (I-1)*LX + 1
c inner loop checks each cell in the row
do J=1, LX
c B is the current location in the lattice
B = (I-1)*LX + J
IF (SN(A) .EQ. 0) THEN
c this is a row with a hole
c hole rows contain holes (#4) and cell #1
IF (SN(B) .EQ. 1) THEN
MX1 = SX(B)
MY1 = SY(B)
M1 = SQRT(MX1**2 + MY1**2)
DO K = 1,101
IF(M1 .GE. BINSIZE(K) .AND. M1 .LT. BINSIZE(K+1)) THEN
MBIN1(K) = MBIN1(K)+1
END IF
IF(MX1 .GE. BINSIZE(K) .AND. MX1 .LT. BINSIZE(K+1)) THEN
MBIN1X(K) = MBIN1X(K)+1
END IF
IF(MY1 .GE. BINSIZE(K) .AND. MY1 .LT. BINSIZE(K+1)) THEN
MBIN1Y(K) = MBIN1Y(K)+1
END IF
end do
ELSE
END IF
ELSE
c this is a row without a hole
c non-hole rows contain cell #s 2 & 3
IF (SN(J) .EQ. 0) THEN
MX3 = SX(B)
MY3 = SY(B)
M3 = SQRT(MX3**2 + MY3**2)
DO K = 1,101
IF(M3 .GE. BINSIZE(K).AND. M3 .LT. BINSIZE(K+1)) THEN
MBIN3(K) = MBIN3(K)+1
END IF
IF(MX3 .GE. BINSIZE(K).AND. MX3 .LT. BINSIZE(K+1)) THEN
MBIN3X(K) = MBIN3X(K)+1
END IF
IF(MY3 .GE. BINSIZE(K).AND. MY3 .LT. BINSIZE(K+1)) THEN
MBIN3Y(K) = MBIN3Y(K)+1
END IF
END DO
ELSE
MX2 = SX(B)
MY2 = SY(B)
72
M2 = SQRT(MX2**2 + MY2**2)
DO K = 1,101
IF(M2 .GE. BINSIZE(K).AND. M2 .LT. BINSIZE(K+1)) THEN
MBIN2(K) = MBIN2(K)+1
END IF
IF(MX2 .GE. BINSIZE(K).AND. MX2 .LT. BINSIZE(K+1)) THEN
MBIN2X(K) = MBIN2X(K)+1
END IF
IF(MY2 .GE. BINSIZE(K).AND. MY2 .LT. BINSIZE(K+1)) THEN
MBIN2Y(K) = MBIN2Y(K)+1
END IF
END DO
END IF
END IF
end do
end do
RETURN
END
************************************************************************
* NEARBR2: A 2nd nearest neighbor arranging subroutine *
* *
************************************************************************
SUBROUTINE NEARBR2
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
INTEGER SN, I, J, N
REAL*8 SX, SY
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMAG, HDMAG, U
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DMG/ DMA, DMB, DMC, DMD, DME,
% MXA, MXB, MXC, MXD, MXE,
% MYA, MYB, MYC, MYD, MYE
COMMON /DMG2/ DMAG, HDMAG
73
I = 1
J = 1
N = LX*LY
C Designations for neighbor locations follow this pattern:
c DM[position]([location], [site number])
c where position is A-E in order of neighbor distance
c and location is in order, clockwise from the top
c and distance is edge distance
C Convert 1st nearest neighbors to new format
DO I=1,N
DMA(I,1) = NN(I,1)
DMA(I,2) = NN(I,2)
DMA(I,3) = NN(I,3)
DMA(I,4) = NN(I,4)
DO J = 5,8
DMA(I,J) = 0
END DO
ENDDO
C Record the site location of each neighbor
DO I=1,N
DMB(I,1) = I - LX + 1
DMB(I,2) = I + LX + 1
DMB(I,3) = I + LX - 1
DMB(I,4) = I - LX -1
C Correct the neighbor locations
C top
IF (I .LE. LX) THEN
DMB(I,1) = DMB(I,1) + N
DMB(I,4) = DMB(I,4) + N
ENDIF
C bottom
IF (I .GT. LX*(LY-1)) THEN
DMB(I,2) = DMB(I,2) - N
DMB(I,3) = DMB(I,3) - N
ENDIF
C left side
IF (MOD((I-1),LX) .EQ. 0) THEN
DMB(I,4) = DMB(I,4) + LX
DMB(I,3) = DMB(I,3) + LX
ENDIF
C right side
IF (MOD(I,LX) .EQ. 0) THEN
DMB(I,1) = DMB(I,1) - LX
DMB(I,2) = DMB(I,2) - LX
74
ENDIF
ENDDO
RETURN
END
************************************************************************
* NEARBR3: A 3rd nearest neighbor arranging subroutine *
* *
************************************************************************
SUBROUTINE NEARBR3
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
INTEGER SN, I, J, N
REAL*8 SX, SY
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMAG, HDMAG, U
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DMG/ DMA, DMB, DMC, DMD, DME,
% MXA, MXB, MXC, MXD, MXE,
% MYA, MYB, MYC, MYD, MYE
COMMON /DMG2/ DMAG, HDMAG
I = 1
J = 1
N = LX*LY
C Designations for neighbor locations follow this pattern:
c DM[position]([location], [site number])
c where position is A-E in order of neighbor distance
c and location is in order, clockwise from the top
c and distance is edge distance
C Record the site location of each neighbor
DO I=1,N
DMC(I,1) = I - 2*LX
DMC(I,2) = I + 2
DMC(I,3) = I + 2*LX
75
DMC(I,4) = I - 2
C Correct the neighbor locations
C top
IF (I .LE. 2*LX) THEN
DMC(I,1) = DMC(I,1) + N
ENDIF
C bottom
IF (I .GT. (N - 2*LX)) THEN
DMC(I,3) = DMC(I,3) - N
ENDIF
C left side
IF (MOD((I-1),LX) .EQ. 0) THEN
DMC(I,4) = DMC(I,4) + LX
ELSE IF (MOD((I-2),LX) .EQ. 0) THEN
DMC(I,4) = DMC(I,4) + LX
ENDIF
C right side
IF (MOD(I,LX) .EQ. 0) THEN
DMC(I,2) = DMC(I,2) - LX
ELSE IF (MOD(I+1,LX) .EQ. 0) THEN
DMC(I,2) = DMC(I,2) - LX
ENDIF
ENDDO
RETURN
END
************************************************************************
* NEARBR4: A 4th nearest neighbor arranging subroutine *
* *
************************************************************************
SUBROUTINE NEARBR4
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
INTEGER SN, I, J, N
REAL*8 SX, SY
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMAG, HDMAG, U
COMMON /DOT/ X, Y, SPX, SPY
76
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DMG/ DMA, DMB, DMC, DMD, DME,
% MXA, MXB, MXC, MXD, MXE,
% MYA, MYB, MYC, MYD, MYE
COMMON /DMG2/ DMAG, HDMAG
I = 1
J = 1
N = LX*LY
C Designations for neighbor locations follow this pattern:
c DM[position]([location], [site number])
c where position is A-E in order of neighbor distance
c and location is in order, clockwise from the top
c and distance is edge distance
C Record the site location of each neighbor
DO I=1,N
DMD(I,1) = I - 2*LX + 1
DMD(I,2) = I - LX + 2
DMD(I,3) = I + LX + 2
DMD(I,4) = I + 2*LX +1
DMD(I,5) = I + 2*LX - 1
DMD(I,6) = I + LX - 2
DMD(I,7) = I - LX - 2
DMD(I,8) = I - 2*LX - 1
C Correct the neighbor locations
C top
IF (I .LE. LX) THEN
DMD(I,1) = DMD(I,1) + N
DMD(I,2) = DMD(I,2) + N
DMD(I,7) = DMD(I,7) + N
DMD(I,8) = DMD(I,8) + N
ENDIF
IF ((I .GT. LX) .AND. (I .LE. 2*LX)) THEN
DMD(I,1) = DMD(I,1) + N
DMD(I,8) = DMD(I,8) + N
ENDIF
C bottom
IF (I .GT. LX*(LY-1)) THEN
DMD(I,3) = DMD(I,3) - N
DMD(I,4) = DMD(I,4) - N
DMD(I,5) = DMD(I,5) - N
DMD(I,6) = DMD(I,6) - N
ENDIF
IF ((I .GT. LX*(LY-2)) .AND. (I .LT. N-LX)) THEN
DMD(I,4) = DMD(I,4) - N
77
DMD(I,5) = DMD(I,5) - N
ENDIF
C left side
IF (MOD((I-1),LX) .EQ. 0) THEN
DMD(I,5) = DMD(I,5) + LX
DMD(I,6) = DMD(I,6) + LX
DMD(I,7) = DMD(I,7) + LX
DMD(I,8) = DMD(I,8) + LX
ENDIF
IF (MOD((I-2),LX) .EQ. 0) THEN
DMD(I,6) = DMD(I,6) + LX
DMD(I,7) = DMD(I,7) + LX
ENDIF
C right side
IF (MOD(I,LX) .EQ. 0) THEN
DMD(I,1) = DMD(I,1) - LX
DMD(I,2) = DMD(I,2) - LX
DMD(I,3) = DMD(I,3) - LX
DMD(I,4) = DMD(I,4) - LX
ENDIF
IF (MOD((I+1),LX) .EQ. 0) THEN
DMD(I,2) = DMD(I,2) - LX
DMD(I,3) = DMD(I,3) - LX
ENDIF
ENDDO
RETURN
END
************************************************************************
* NEARBR5: A 5th nearest neighbor arranging subroutine *
* *
************************************************************************
SUBROUTINE NEARBR5
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
INTEGER SN, I, J, N
REAL*8 SX, SY
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMAG, HDMAG, U
78
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DMG/ DMA, DMB, DMC, DMD, DME,
% MXA, MXB, MXC, MXD, MXE,
% MYA, MYB, MYC, MYD, MYE
COMMON /DMG2/ DMAG, HDMAG
I = 1
J = 1
N = LX*LY
C Designations for neighbor locations follow this pattern:
c DM[position]([location], [site number])
c where position is A-E in order of neighbor distance
c and location is in order, clockwise from the top
c and distance is edge distance
C Record the site location of each neighbor
DO I=1,N
DME(I,1) = I - 2*LX + 2
DME(I,2) = I + 2*LX + 2
DME(I,3) = I + 2*LX - 2
DME(I,4) = I - 2*LX - 2
C Correct the neighbor locations
C top
IF (I .LE. LX) THEN
DME(I,1) = DME(I,1) + N
DME(I,4) = DME(I,4) + N
ENDIF
IF ((I .GT. LX) .AND. (I .LE. 2*LX)) THEN
DME(I,1) = DME(I,1) + N
DME(I,4) = DME(I,4) + N
ENDIF
C bottom
IF (I .GT. LX*(LY-1)) THEN
DME(I,2) = DME(I,2) - N
DME(I,3) = DME(I,3) - N
ENDIF
IF ((I .GT. LX*(LY-2)) .AND. (I .LT. N-LX)) THEN
DME(I,2) = DME(I,2) - N
DME(I,3) = DME(I,3) - N
ENDIF
C left side
IF (MOD((I-1),LX) .EQ. 0) THEN
DME(I,3) = DME(I,3) + LX
DME(I,4) = DME(I,4) + LX
79
ENDIF
IF (MOD((I-2),LX) .EQ. 0) THEN
DME(I,3) = DME(I,3) + LX
DME(I,4) = DME(I,4) + LX
ENDIF
C right side
IF (MOD(I,LX) .EQ. 0) THEN
DME(I,1) = DME(I,1) - LX
DME(I,2) = DME(I,2) - LX
ENDIF
IF (MOD((I+1),LX) .EQ. 0) THEN
DME(I,1) = DME(I,1) - LX
DME(I,2) = DME(I,2) - LX
ENDIF
ENDDO
RETURN
END
************************************************************************
* DMAGC: A demagnetization calculator to the 5th nearest neighbor *
* *
************************************************************************
SUBROUTINE DMAGC
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
INTEGER SN, I, J, K, N, IDM
REAL*8 SX, SY
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMDISTH(20,2), DMDISTV(20,2)
REAL*8 DMAG, HDMAG, UHX, UHY, UVX, UVY
REAL*8 X1,Y1,X2,Y2
REAL*8 MNETH(10000,20), MNETV(10000,20)
C Variables are named where the MX series is the magnetization x-components
C and the MY series is the magnetization y-components of the neighbors
C DMDIST records the magnitude of the distances to each neighbor
C positive and negative signs will be added later to make the data more manageable
80
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DMG/ DMA, DMB, DMC, DMD, DME,
% MXA, MXB, MXC, MXD, MXE,
% MYA, MYB, MYC, MYD, MYE
COMMON /DMG2/ DMAG, HDMAG,IDM
COMMON /DIST/ DMDISTH, DMDISTV
I = IDM
J = 1
K = 1
N = LX*LY
UXH = 0
UYH = 0
UXV = 0
UYV = 0
HDMAG = 0
c record the component of the magnetization of each neighbor from its site
IF (DMAG .ne. 0) THEN
DO J = 1, 4
MXA(I, J)=SX(DMA(I, J))
MYA(I, J)=SY(DMA(I, J))
MXB(I, J)=SX(DMB(I, J))
MYB(I, J)=SY(DMB(I, J))
MXC(I, J)=SX(DMC(I, J))
MYC(I, J)=SY(DMC(I, J))
MXD(I, J)=SX(DMD(I, J))
MYD(I, J)=SY(DMD(I, J))
MXE(I, J)=SX(DME(I, J))
MYE(I, J)=SY(DME(I, J))
END DO
DO J = 5,8
MXD(I, J)=SX(DMD(I, J))
MYD(I, J)=SY(DMD(I, J))
END DO
C For the purpose of organization, we'll reorder these
C into the sums we will use in our calculation
c This section will be very long
MNETH(I,1)=-MYE(I,4)+MYD(I,7)
MNETH(I,2)=-MYD(I,8)+MYB(I,4)
MNETH(I,3)=-MYC(I,1)+MYA(I,1)
MNETH(I,4)=-MYD(I,1)+MYB(I,1)
MNETH(I,5)=-MYE(I,1)+MYD(I,2)
MNETH(I,6)=-MYD(I,7)+MYC(I,4)
MNETH(I,7)=-MYB(I,4)+MYA(I,4)
MNETH(I,8)=-MYA(I,1)+SY(I)
81
MNETH(I,9)=-MYB(I,1)+MYA(I,2)
MNETH(I,10)=-MYD(I,2)+MYC(I,2)
MNETH(I,11)=-MYC(I,4)+MYD(I,6)
MNETH(I,12)=-MYA(I,4)+MYB(I,3)
MNETH(I,13)=-SY(I)+MYA(I,3)
MNETH(I,14)=-MYA(I,2)+MYB(I,2)
MNETH(I,15)=-MYC(I,2)+MYD(I,3)
MNETH(I,16)=-MYD(I,6)+MYE(I,3)
MNETH(I,17)=-MYB(I,3)+MYD(I,5)
MNETH(I,18)=-MYA(I,3)+MYC(I,3)
MNETH(I,19)=-MYB(I,2)+MYD(I,4)
MNETH(I,20)=-MYD(I,3)+MYE(I,2)
MNETV(I,1)=MXE(I,4)-MXD(I,8)
MNETV(I,2)=MXD(I,8)-MXC(I,1)
MNETV(I,3)=MXC(I,1)-MXD(I,1)
MNETV(I,4)=MXD(I,1)-MXE(I,1)
MNETV(I,5)=MXD(I,7)-MXB(I,4)
MNETV(I,6)=MXB(I,4)-MXA(I,1)
MNETV(I,7)=MXA(I,1)-MXB(I,1)
MNETV(I,8)=MXB(I,1)-MXD(I,2)
MNETV(I,9)=MXC(I,4)-MXA(I,4)
MNETV(I,10)=MXA(I,4)-SX(I)
MNETV(I,11)=SX(I)-MXA(I,2)
MNETV(I,12)=MXA(I,2)-MXC(I,2)
MNETV(I,13)=MXD(I,6)-MXB(I,3)
MNETV(I,14)=MXB(I,3)-MXA(I,3)
MNETV(I,15)=MXA(I,3)-MXB(I,2)
MNETV(I,16)=MXB(I,2)-MXD(I,3)
MNETV(I,17)=MXE(I,3)-MXD(I,5)
MNETV(I,18)=MXD(I,5)-MXC(I,3)
MNETV(I,19)=MXC(I,3)-MXD(I,4)
MNETV(I,20)=MXD(I,4)-MXE(I,2)
C Calculate the net demagnetization effect on the site.
c note that the magnetization is the sum of the neighbors
c on either side of the edge being considered
c 20 horizontal and 20 vertical edges
DO J = 1,20
c sum the horizontal edges, x components
UXH = UXH + MNETH(I,J)*DMDISTH(J,1)
c sum the horizontal edges, y components
UYH = UYH + MNETH(I,J)*DMDISTH(J,2)
c sum the vertical edges, x components
UXV = UXV + MNETV(I,J)*DMDISTV(J,1)
c sum the vertical edges, y components
UYV = UYV + MNETV(I,J)*DMDISTV(J,2)
82
END DO
c return the scaled demagnetization energy
HDMAG = DMAG * (SX(I)*(UXH+UXV)+SY(I)*(UYH+UYV))
END IF
RETURN
END
************************************************************************
* DMDIST: A neighbor distance calculator *
* *
************************************************************************
SUBROUTINE DMDIST
INTEGER A, B, C, D
INTEGER NUMX, NUMY, VAR
INTEGER X, Y, SPX, SPY
INTEGER L, NN
INTEGER SN, I, J, K, N
REAL*8 SX, SY
REAL*8 DMA(10000,10), DMB(10000,10), DMC(10000,10), DMD(10000,10),
% DME(10000,10)
REAL*8 MXA(10000,10), MXB(10000,10), MXC(10000,10), MXD(10000,10),
% MXE(10000,10)
REAL*8 MYA(10000,10), MYB(10000,10), MYC(10000,10), MYD(10000,10),
% MYE(10000,10)
REAL*8 DMDISTH(20,2), DMDISTV(20,2)
REAL*8 UHX, UHY, UVX, UVY
REAL*8 X1,Y1,X2,Y2
COMMON /DOT/ X, Y, SPX, SPY
COMMON /NN/ LX, LY, NN(10000,4)
COMMON /NN2/ SX(10000), SY(10000), SN(100000)
COMMON /DIST/ DMDISTH, DMDISTV
x1=0
x2=0
y1=0
y2=0
N = LX*LY
C Record the distances to each neighbor edge
C the distance array is in order like the lattice array
c note that the distances are the same for several edges;
c this is due to symmetry
C Also, x1 and y1 are Xmax and Ymax, respectively
C As such, x2 and y2 are thus Xmin and Ymin
C note that to simplify data use, we are calculating the distance
C component of the demagnetization energy here instead of just
83
C recording x1,x2,y1,y2
C this makes the later calculation much simpler to follow
C also the index L in DMDISTH(K,L) is the component, i.e. 1=x, 2=y
C First calculate the distances for the horizontal components
Y1 = 1.5
DO J = 1,4
X1 = -2.5
DO I = 1,5
X2 = X1 + 1
K = 5*(J-1)+I
DMDISTH(K,1) = (ABS(y1)**(-1))*((x2/y1)**2+1)**(-0.5D0)
% -(ABS(y1)**(-1))*((x1/y1)**2+1)**(-0.5D0)
DMDISTH(K,2) = ( -(x2/(y1**2))*((x2/y1)**2+1)**(-0.5D0)
% +(x1/(y1**2))*((x1/y1)**2+1)**(-0.5D0) )*(y1/abs(y1))
c the term (y1/abs(y1)) corrects for vector direction
X1 = X2
END DO
Y1 = Y1 - 1
END DO
C Now calculate the vertical segments
Y1 = 1.5
DO J = 1, 5
Y2 = Y1 + 1
X1 = -1.5
DO I = 1, 4
K = 4*(J-1)+I
DMDISTV(K,1) = ( -(y2/(x1**2))*((y2/x1)**2+1)**(-0.5D0)
% +(y1/(x1**2))*((y1/x1)**2+1)**(-0.5D0) )*(x1/abs(x1))
DMDISTV(K,2) = (ABS(x1)**(-1))*((y2/x1)**2+1)**(-0.5D0)
% -(ABS(x1)**(-1))*((y1/x1)**2+1)**(-0.5D0)
X1 = X1 + 1
END DO
Y1 = Y1 - 1
END DO
RETURN
END
************************************************************************
