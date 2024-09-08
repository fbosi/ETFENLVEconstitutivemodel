!****************************************************************************	
*USER SUBROUTINES
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

	INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME  
      DIMENSION STRESS(NTENS),STATEV(NSTATV), 
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
	 
	DOUBLE PRECISION PROPS
	 
	REAL*16 CTEMD, CTETD, CTE3, CTEMDold, CTETDold, CTE3old    
      REAL*16 Ea, Tr, T0, Temperature, Poisson, V, DPKeqDt, maxEqS
      INTEGER NProny, counter,istat
      REAL*16 tauj(21),Etresh
      REAL*16 S110, S210, S120, S130, S220, S230, S660 
      REAL*16 S11j(21),S12j(21),S21j(21),S13j(21),S22j(21),S23j(21),S66j(21)
      REAL vE(3)
	REAL vEold(3)


      REAL*16 PK1pold(NTENS)
      REAL*16 q11old(21),q12old(21),q13old(21)
	REAL*16 q21old(21),q22old(21),q23old(21)
	REAL*16 q66old(21)
      
      REAL*16 aT, dtr, Et(3) 
      REAL*16 S11, S12, S13, S22, S23, S66
      REAL*16 f11, f12, f13, f21, f22, f23, f66
      REAL*16 DFGRD1T(3,3), I(3,3), E(3,3), PK1(3,3), C(3,3), det, J
      REAL*16 dPK1dE(3,3), F(3,3)
      
	REAL*16  omega(3),Ueigval(3),eigvec(3,3)
      REAL*16  U(3,3)
	
      REAL*16 Ep(3), Eppr(3), dif(3)
      
      REAL*16 lambdap(2), lambdat, theta, detF, dive, oldEqStr
      
      REAL*16 temp1(3,3), temp2(3,3), sigma(3,3)
	
	 
!*************************** PROPS DESIGNATION *************************** 
!	 User-specified array of material constants associated with this user material:
!
!     PROPS  1 - 7  : CTEMD % Coefficient of Thermal Expansion MD
!     PROPS 8 - 14  : CTETD % Coefficient of Thermal Expansion TD
!     PROPS 15 - 21 : CTE3  % Coefficient of Thermal Expansion out of plane
!     PROPS 22 - 23 : Arrhenius   % Arrhenius Constants 
!     PROPS 24 - 44 : tauj  % Retardation Times
!     PROPS 45      : D01    % Prony First Coefficient MD direction
!     PROPS 46 - 66 : Dj1    % Prony Coefficients MD direction
!     PROPS 67      : V    % V activation volume
!     PROPS 68      : D01    % Prony First Coefficient MD direction
!     PROPS 69 - 89 : Dj1    % Prony Coefficients MD direction
!     PROPS 90      : D01    % Prony First Coefficient MD direction
!     PROPS 91 - 111 : Dj1    % Prony Coefficients MD direction

!***************************  STATEV DESIGNATION *************************** 
!***************************    State variable   *************************** 
!
!     STATEV   1 -   3 : PK1old  % PK1 stresses at the beginning of increment - sigma3=0 
!     STATEV   4 - 151 : qold    % hereditary integrals at the beginning of increment
!     STATEV 152- 155 : E       % Mech strain at the beginning of increment - eps3 different from 0 
!     STATEV 156 : E1tot         % Total strain MD
!     STATEV 157 : E2tot         % Total strain TD
!     STATEV 158 : E3tot         % Total strain 3
!     STATEV 159 : T0            % temperature of the step
!     STATEV 160 : DsEq/dt	     % stress increment in time
!     STATEV 161 : masSeq	     % maximum equivalent stress for Eyring
!     STATEV 162 : Seq	         % equivalent stress
!     STATEV 163 : eeq	         % equivalent strain
!     STATEV 164 : eReq	         % equivalent strain rate
!
!***************************   Material parameters  **************************** 
                       
!	Coefficients of thermal expansion - CTE - current step -use T in Celsius!
	  Temperature = TEMP + 273.15
      CTEMD=PROPS(1)*Temperature**6+PROPS(2)*Temperature**5+PROPS(3)*Temperature**4
     & +PROPS(4)*Temperature**3+PROPS(5)*Temperature**2+PROPS(6)*Temperature**1
     & +PROPS(7)
	  CTETD=PROPS(8)*Temperature**6+PROPS(9)*Temperature**5+PROPS(10)*Temperature**4
     & +PROPS(11)*Temperature**3+PROPS(12)*Temperature**2+PROPS(13)*Temperature**1
     & +PROPS(14)
	CTE3=PROPS(15)*Temperature**6+PROPS(16)*Temperature**5+PROPS(17)*Temperature**4
     & +PROPS(18)*Temperature**3+PROPS(19)*Temperature**2+PROPS(20)*Temperature**1
     & +PROPS(21)

	 
!   Arrhenius constant 
      Ea = PROPS(22)
      Tr = PROPS(23)

!   Number of Prony terms
	NProny=21 	
		
!   Obtain compliance Prony coefficients 
	  Poisson = 0.43
      do j=1,NProny
        tauj(j)=PROPS(23+j)
		S11j(j)=PROPS(45+J)
		S12j(j)=-Poisson*(PROPS(45+j)+PROPS(69+j))/2
		S21j(j)=-Poisson*(PROPS(45+j)+PROPS(69+j))/2
        S22j(j)=PROPS(69+j)
		S66j(j)=PROPS(91+j)
	    S13j(j)=S12j(j)
        S23j(j)=S12j(j)
      end do
	  	
!   Obtain instantaneous compliance 
	  S110=PROPS(45)
	  S120=-Poisson*(PROPS(45)+PROPS(68))/2
	  S210=-Poisson*(PROPS(45)+PROPS(68))/2
	  S220=PROPS(68)
	  S660=PROPS(90)*((1+Poisson))
	  S130=-Poisson*(PROPS(45)+PROPS(68))/2
	  S230=-Poisson*(PROPS(45)+PROPS(68))/2
	  
 
!****************** solution-dependent variables **************************

!     Temperature of the previous step      
	  if ((KSTEP == 1 .AND. KINC == 1))  then
	     T0 = TEMP+273.15
		 STATEV(159) = T0
	  else
		 T0 = STATEV(159)
	  end if
	  
!	Coefficients of thermal expansion - CTE - old step

      CTEMDold=PROPS(1)*T0**6+PROPS(2)*T0**5+PROPS(3)*T0**4
     & +PROPS(4)*T0**3+PROPS(5)*T0**2+PROPS(6)*T0**1
     & +PROPS(7)
	  CTETDold=PROPS(8)*T0**6+PROPS(9)*T0**5+PROPS(10)*T0**4
     & +PROPS(11)*T0**3+PROPS(12)*T0**2+PROPS(13)*T0**1
     & +PROPS(14)
	  CTE3old=PROPS(15)*T0**6+PROPS(16)*T0**5+PROPS(17)*T0**4
     & +PROPS(18)*T0**3+PROPS(19)*T0**2+PROPS(20)*T0**1
     & +PROPS(21)	

!     Add thermal strain
      Et(1)=(CTEMD-CTEMDold)*0.000000000000001
      Et(2)=(CTETD-CTETDold)*0.000000000000001
      Et(3)=(CTE3-CTE3old)*0.000000000000001		
	  
!     First Piola-Kirchhoff stresses at the beginning of increment
      do j=1,3
		PK1pold(j) = STATEV(j)
	  end do 

	! Eyring variables	
	  V = PROPS(67) 
	  DPKeqDt = STATEV(160)
	  maxEqS = STATEV(161)
      Etresh = 0.d0

!     Heredity integrals at the beginning of increment
	counter = 4 ! first position of qold in STATEV
	do j = 1,NProny
    	  q11old(j) = STATEV(counter) 
    	  q12old(j) = STATEV(counter+NProny)
    	  q13old(j) = STATEV(counter+2*NProny)
    	  q21old(j) = STATEV(counter+3*NProny)
    	  q22old(j) = STATEV(counter+4*NProny)
    	  q23old(j) = STATEV(counter+5*NProny)
    	  q66old(j) = STATEV(counter+6*NProny)
	      counter = counter+1   
	  
	end do


!****************** Obtain Nominal Strain ***************************************

!     compute inplane nominal strain at the end of increment 
!     3-axis of deformation gradient F is always normal to current membrane surface
!     current deformation gradient F for membrane elements is always symmetric 

      call onem(U)
      F=DFGRD1
	  DFGRD1T=transpose(DFGRD1)             
      C=matmul(DFGRD1T,DFGRD1)         
      call onem(I)   
	  E = 0.d0
	  
! recall subroutine to compute the principal values of U and E
	call spectral(C,omega,eigvec,istat)
	
! Calculate the principal values of U and E
      Ueigval(1) = sqrt(omega(1))
      Ueigval(2) = sqrt(omega(2))
      Ueigval(3) = sqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
	
! Calculate the complete tensors U and E
	U = matmul(matmul(eigvec,U),transpose(eigvec))
	E=U-I   
      
!****************** Update stress ***************************************	

! Strain rate planar
		
      if (KSTEP == 1 .AND. KINC == 1)  then
		DEPLDt = 0.d0
      else 
		  vE(1) = E(1,1)-Et(1)
		  vE(2) = E(2,2)-Et(2)
		  vE(3) = E(1,2)-Et(3)
		  vEold(1) = STATEV(152)
		  vEold(2) = STATEV(153)
		  vEold(3) = STATEV(155)

		  LocIndex = MAXLOC(vE, DIM=1)
		  if (DTIME >0) then
			DEPLDt = (vE(LocIndex)-vEold(LocIndex))/DTIME
	      else
			DEPLDt = 0.d0
		  end if  
		  
	  end if
!	Arrhenius Shift factor

      aT=10**((-Ea/(2.303*8.31446261815324))*(-1/Temperature+1/Tr))
	  
!    Equivalent stress
	  sEq = (0.5*((STATEV(1)-STATEV(2))**2+(STATEV(1))**2+(STATEV(2))**2)+3*STATEV(3)**2)**0.5

!	Eyring law switch for unloading

	  if (DEPLDt < Etresh) then
	      sEyring = maxEqS*1e-6
	  else
	      sEyring = sEq*1e-6
	  end if
	  x = ((V*sEyring)/(Temperature*8.31446262*1e-3))
	  arg = sinh(x)
      if (arg<1e-10) then
		  arg = x + 1/6*x**3 + 1/120*x**5;
      end if
         aSigma = (V*sEyring/(Temperature*arg*8.31446262*1e-3))
 
      if (KSTEP == 1 .AND. KINC == 1)  then
		aSigma = 1
      end if
	  
	  if ((sEq<=0) .AND. (DPKeqDt >= 0)) then 
 		aSigma = 1
      end if     	

!     Reduced time increment
	  dtr=DTIME/(aT*aSigma)
		  
!     compute S and f
      S11 = S110
	  S21 = S210
      S12 = S120
      S13 = S130
      S22 = S220
      S23 = S230
      S66 = S660
      f11=0 
      f12=0 
      f13=0 
      f21=0
      f22=0
      f23=0
      f66=0
	  !counter2 = 201
      do j = 1, NProny
            stq=dtr/tauj(j)
                if (stq<1.0D-10) then
                    stp=1-stq/2+stq**2/6-stq**3/24
                else
                    stp=(1-exp(-stq))/stq
                end if
        S11=S11+S11j(j)-S11j(j)*stp
        S12=S12+S12j(j)-S12j(j)*stp
		S21=S21+S21j(j)-S21j(j)*stp
        S13=S13+S13j(j)-S13j(j)*stp
        S22=S22+S22j(j)-S22j(j)*stp
        S23=S23+S23j(j)-S23j(j)*stp
        S66=S66+S66j(j)-S66j(j)*stp
        f11=f11+S11j(j)*(exp(-stq)*q11old(j)-stp*PK1pold(1))
        f12=f12+S12j(j)*(exp(-stq)*q12old(j)-stp*PK1pold(2))
        f13=f13+S13j(j)*(exp(-stq)*q13old(j)-stp*PK1pold(1))
        f21=f21+S21j(j)*(exp(-stq)*q21old(j)-stp*PK1pold(1))
        f22=f22+S22j(j)*(exp(-stq)*q22old(j)-stp*PK1pold(2))
        f23=f23+S23j(j)*(exp(-stq)*q23old(j)-stp*PK1pold(2))
        f66=f66+S66j(j)*(exp(-stq)*q66old(j)-stp*PK1pold(3)) 
      end do

!   if (DTEMP .NE. 0) then
	  E(1,1)=E(1,1)-Et(1)
      E(2,2)=E(2,2)-Et(2)
	  E(3,3)=E(3,3)-Et(3)
!      end if

!     PK1 stresses in material directions 
      det=S11*S22-S12*S21
      PK1(1,1)=(S22*(E(1,1)+f11+f12) 
     & -S21*(E(2,2)+f21+f22))/det
      PK1(2,2)=(-S12*(E(1,1)+f11+f12)
     & +S11*(E(2,2)+f21+f22))/det
      PK1(3,3)=0.0
      PK1(1,2)=(E(1,2)+f66)/S66
	  PK1(1,3)=0.0
	  PK1(2,1)=(E(1,2)+f66)/S66
	  PK1(2,3)=0.0
	  PK1(3,1)=0.0
      PK1(3,2)=0.0

	  sEqNew = (0.5*((PK1(1,1)-PK1(2,2))**2+(PK1(1,1))**2+(PK1(2,2))**2)+3*PK1(1,2)**2)**0.5
	  STATEV(160) = (sEqNew-sEq)/DTIME 
	  DEPLDt = 100
	  if ((STATEV(160)>0) .AND. (DEPLDt > Etresh)) then
		STATEV(161) = sEqNew 
	  else 
	    STATEV(161) = maxEqS 
      end if

	  
!     Nominal strains in material directions
      E(3,3)=S13*PK1(1,1)-f13+S23*PK1(2,2)-f23

	!eq Stress, eq strain old and current,strainRate,InstElModulus and yield boolean
      STATEV(162) = sEqNew
      oldEqStr = STATEV(163)
      STATEV(163) = (1/(1+Poisson))*(0.5*((E(1,1)-E(2,2))**2+(E(1,1)-E(3,3))**2+(E(2,2)-E(3,3))**2)+3*E(1,2)**2)**0.5
	
	
      call mdet(F,J)
      if (J.le.0.d0) then
	  write(*,*) '** ERROR **'
      write(*,*) 'Problem in kinematics, J.le.0, J='
      write(*,*) J
      write(*,*) '** ERROR **'
        call xit
      end if	  
	  
      temp2=matmul(PK1,DFGRD1T)
      sigma=temp2/J
 

!******************************************************************************

      STRESS(1) = sigma(1,1)
      STRESS(2) = sigma(2,2)
      STRESS(3) = sigma(1,2)
	  
!************************   update tangent stiffness   ************************

!     Tangent stiffness in material directions
      det=S11*S22-S12*S21
      dPK1dE(1,1)=S22/det
      dPK1dE(1,2)=-S12/det
      dPK1dE(1,3)=0.0
      dPK1dE(2,1)=-S21/det
      dPK1dE(2,2)=S11/det
      dPK1dE(2,3)=0.0
      dPK1dE(3,1)=0.0
      dPK1dE(3,2)=0.0
      dPK1dE(3,3)=1/S66    

!     Push forward tangent stiffness
      DDSDDE=dPK1dE    
	  


!********************  Update solution-dependent variables  ******************** 

!     Update PK2 stresses at the beginning of increment
      STATEV(1) = PK1(1,1)
      STATEV(2) = PK1(2,2)
      STATEV(3) = PK1(1,2)

!     Update hereditary integrals at the end of increment
      counter = 4
      do j = 1, NProny	
	   stq=dtr/tauj(j)
                if (stq<1.0D-10) then
                    stp=1-stq/2+stq**2/6-stq**3/24
                else
                    stp=(1-exp(-stq))/stq
                end if
        STATEV(counter)=exp(-stq)*q11old(j)+
     & (PK1(1,1)-PK1pold(1))*stp
    	  STATEV(counter+NProny)=exp(-stq)*q12old(j)+
     & (PK1(2,2)-PK1pold(2))*stp
        STATEV(counter+2*NProny)=exp(-stq)*q13old(j)+ 
     & (PK1(1,1)-PK1pold(1))*stp
        STATEV(counter+3*NProny)=exp(-stq)*q21old(j)+ 
     & (PK1(1,1)-PK1pold(1))*stp
        STATEV(counter+4*NProny)=exp(-stq)*q22old(j)+ 
     & (PK1(2,2)-PK1pold(2))*stp
        STATEV(counter+5*NProny)=exp(-stq)*q23old(j)+ 
     & (PK1(2,2)-PK1pold(2))*stp
        STATEV(counter+6*NProny)=exp(-stq)*q66old(j)+ 
     & (PK1(1,2)-PK1pold(3))*stp
        counter=counter+1
      end do

!     Update Nominal strain at the end of increment
      STATEV(152)= E(1,1) 
      STATEV(153)= E(2,2) 
      STATEV(154)= E(3,3) 
      STATEV(155)= E(1,2) 
	
!     Total strains 	
	STATEV(156)=E(1,1)+Et(1)    
	STATEV(157)=E(2,2)+Et(2)    
	STATEV(158)=E(3,3)+Et(3)   

!********************  Update energy   **************************************	

!     Update elastic strain energy in SSE
      SSE = 0.0
    
!     Update creep dissipation in SCD
      SCD = 0.0

!     Update plastic dissipation in SPD
	SPD = 0.0
	 
      RETURN

      END
!****************************************************************************	
!****************************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!****************************************************************************
      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      INCLUDE 'ABA_PARAM.INC'
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*16 D(3),V(3,3),A(3,3),K(3,3)

      K = A
      !
      call jacobi(K,3,np,D,V,nrot,istat)
	
      return
      end
!****************************************************************************
!****************************************************************************     
      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      INCLUDE 'ABA_PARAM.INC'
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*16 A(np,np),D(np),V(np,np),B(nmax),Z(nmax)
      real*16 sm,tresh,G,T,H,theta,S,C,tau
       
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
	  B(ip) = A(ip,ip)
	  D(ip) = B(ip)
	  Z(ip) = 0.d0
      end do
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,100
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
	      sm = sm + abs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*abs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the 
              !  off-diagonal element is small.
              !                                                                                                                                
	      if ((i.gt.4).and.(abs(D(ip))+G.eq.abs(D(ip))).and.  
     &	     (abs(D(iq))+G.eq.abs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (abs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (abs(H)+G.eq.abs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
	          T =A(ip,iq)/H
	        else
	          theta = 0.5d0*H/A(ip,iq)
	          T =1.d0/(abs(theta)+sqrt(1.d0+theta**2.d0))
	          if (theta.lt.0.d0) T = -T
	        end if
	        C = 1.d0/sqrt(1.d0 + T**2.d0)
	        S = T*C
	        tau = S/(1.d0 + C)
	        H = T*A(ip,iq)
	        Z(ip) = Z(ip) - H
	        Z(iq) = Z(iq) + H
	        D(ip) = D(ip) - H
	        D(iq) = D(iq) + H
	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
		!		
	        do j=1,ip-1
	          G = A(j,ip)
	          H = A(j,iq)
	          A(j,ip) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations P < J < Q
                !
	        do j=ip+1,iq-1
	          G = A(ip,j)
	          H = A(j,iq)
	          A(ip,j) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations Q < J <= N
                !
	        do j=iq+1,n
                  G = A(ip,j)
	          H = A(iq,j)
	          A(ip,j) = G - S*(H + G*tau)
	          A(iq,j) = H + S*(G - H*tau)
	        end do
	        do j = 1,n
	          G = V(j,ip)
	          H = V(j,iq)
	          V(j,ip) = G - S*(H + G*tau)
	          V(j,iq) = H + S*(G - H*tau)
	        end do
	        nrot = nrot + 1
              end if
	    end do
	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
	  do ip=1,n
	    B(ip) = B(ip) + Z(ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  end do
	end do

      !  If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '100 iterations in Jacobi should never happen'
      istat = 0
      
      return
      end subroutine jacobi	
!****************************************************************************
!****************************************************************************
      subroutine mdet(A,dett)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*16  A(3,3),dett


      dett = A(1,1)*A(2,2)*A(3,3)  
     &     + A(1,2)*A(2,3)*A(3,1)  
     &     + A(1,3)*A(2,1)*A(3,2)  
     &     - A(3,1)*A(2,2)*A(1,3)   
     &     - A(3,2)*A(2,3)*A(1,1)   
     &     - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
!****************************************************************************
!****************************************************************************
      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      INCLUDE 'ABA_PARAM.INC'
	!
      !
      integer i,j
      !
      real*16 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.d0
            else
              A(i,j) = 0.d0
            end if
         end do
      end do

      return
      end 
!****************************************************************************
