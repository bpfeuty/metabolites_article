!! --gfortran simu.f90 -lseulexé -lminpack -llapack; ./a.out  
!! !units micromolar and seconds !!!

module GLOB
  implicit none
  integer,parameter::nvar=8,npar=19,nreg=5,nreac=12,nplot=50000
  integer :: nr_max,nmax,numit,nrf
  real(8), dimension(nvar,nplot) :: xplot(nvar,nplot),tplot(nplot)
  real(8) :: pi=acos(-1.d0),l2=log(2D0),l10=log(1D1)
  integer, parameter :: max_val=10**3
  real(8) :: p0,k0,k_nad,k_g6pd,k_pgd,k_tkt,k_gly1,k_gly3,k_nuc,k_gapd,k_pfk,k_pfb,ksat
  real(8) :: keq_tkt,keq_gly1,keq_gly3,rr(nreg),ss(2)
  type indiv  
     sequence
     real(8) :: k(npar),sc
  end type indiv
  type(indiv) :: p,pb
  type(indiv),dimension(10)::pfam

contains

  
  subroutine sys_ox(v,dvdt,t) !D,X,Y 
    real(8),dimension(nvar),intent(in)::v
    real(8),dimension(nvar),intent(out)::dvdt
    real(8),intent(in)::t
    real(8)::pox1,pox2,ptkt,pnuc,pgly1,pgly2,pgap,pnad,pgly3
    real(8)::G6P,R5P,GAP,NADPH,F6P,PG,FBP,H2O2
    G6P=v(1); F6P=v(2); FBP=v(3); GAP=v(4); PG=v(5); R5P=v(6); NADPH=v(7); H2O2=v(8);
    
    pox1=k_g6pd*G6P*(1-NADPH)/(1.+rr(1)*NADPH);
    pox2=k_pgd*PG*(1-NADPH)/(1.+rr(2)*PG);
    ptkt=k_tkt*(R5P**3-F6P**2*GAP/keq_tkt)/(1+(R5P**3+F6P**2*GAP)*ksat);
    pgly1=k_gly1*(G6P-F6P/keq_gly1)/(1.+rr(3)*PG)
    pgly2=k_pfk*F6P/(1+rr(4)*H2O2)-k_pfb*FBP;     
    pgly3=k_gly3*(FBP-GAP/keq_gly3)
    pgap=k_gapd*GAP/(1+rr(5)*H2O2);
    pnuc=k_nuc*R5P;  pnad=k_nad*NADPH*H2O2;              
    dvdt(1)=p0-pox1-pgly1;      !G6P
    dvdt(2)=pgly1-pgly2+2*ptkt; !F6P
    dvdt(3)=pgly2-pgly3;        !FBP
    dvdt(4)=2*pgly3-pgap+ptkt;  !gap
    dvdt(5)=pox1-pox2;          !PG
    dvdt(6)=pox2-pnuc-3*ptkt;   !r5P
    dvdt(7)=pox1+pox2-pnad;     !NADPH
    dvdt(8)=k0-pnad;     !NADPH
  end subroutine sys_ox

  !! routine d'integration numerique qui utilise SEULEX !!
  subroutine evolSEU(tbeg,tend,X,idid,tol)
    real(8), intent(in) :: tol
    real(8)             :: tbeg,tend
    real(8), intent(inout) :: X(nvar)
    integer, intent(out) :: idid
    real(8),save::h
    integer,parameter :: lrcont=20,licont=10
    integer :: itol,ijac,ifcn,iout
    integer :: mljac=Nvar,mujac=Nvar,imas=0
    integer :: MUMAS=Nvar,MLMAS=Nvar
    integer, parameter :: lwork=2*(Nvar+10)*Nvar+61+20
    real(8) :: work(lwork),rpar
    integer, parameter :: liwork=2*Nvar+32+20
    integer :: iwork(liwork),ipar
    real(8) :: Rtol,Atol
    external fcn_HW,solout,mass,jac_HW !! code externe utilisé
    itol=0 
    rtol=tol
    atol=1d-6  	!! précision abolue
    iout=1 	!! appelle soulout
    ijac=0 	!! difference finie pour le jacobien -> =1  utilise JAC
    ifcn=0 	!! F(X,Y) INDEPENDENT OF X (AUTONOMOUS) -> =1 nonautonomous

    work=0d0
    !! work(2)=0.5D0 !! mettre en commentaire si adaptatif
    iwork=0
    iwork(1)=0 		!! pas de transfo hessenberg du jac -> =1 sinon
    iwork(2)=10**8 	!! nb de step max
    iwork(4)=3  	!! step size sequence
    call SEULEX(Nvar,FCN_HW,IFCN,tbeg,X,tend,H,&
         &RTOL,ATOL,ITOL,&
         &JAC_HW ,IJAC,MLJAC,MUJAC,&
         &MASS ,IMAS,MLMAS,MUMAS,&
         &SOLOUT,IOUT,&
         &WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    nr_max=min(iwork(17),nplot) !! nombre de pas
    !idid
  end subroutine evolSEU

end  module GLOB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_algebra   !! EXPM, EIGEN, JACSOL
  use GLOB; implicit none;
contains
  FUNCTION EXPM(A,n) ! Matrix exponential EXPM(A, n) = exp(A)
    INTEGER i,n
    REAL(kind=kind(0.d0)),DIMENSION(n,n)::C,expm,A
    C=A;
    EXPM=A;
    DO i=1,n
       expm(i,i)=expm(i,i)+1.0
    END DO
    DO i=2,16
       C=MATMUL(C, A)/I
       expm=expm+C
    END DO
  END FUNCTION EXPM
  
  function inv(A) result(Ainv)
    real(8), dimension(:,:), intent(in) :: A
    real(8), dimension(size(A,1),size(A,2)) :: Ainv

    real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv

  subroutine eigen(x,eig)  ! calculate eigenvalues using LAPACK; x->jac
    integer::N,LWORK,LDVR,LDVL,LDA,info;
    PARAMETER (N=NVAR, LWORK=108, LDVR=N, LDVL=N, LDA=N)
    character::JOBVL,JOBVR
    real(kind=kind(0.d0))::VR(LDVR,N),VL(LDVL,N),WR(N),WI(N),WORK(LWORK)
    real(kind=kind(0.d0)),intent(out)::eig(N,2)
    real(kind=kind(0.d0)),intent(in)::x(LDA,N)
    JOBVL='N';JOBVR='N';
    call DGEEV(JOBVL,JOBVR,N,x,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
    eig(:,1)=WR(:); eig(:,2)=WI(:);
  end subroutine eigen

  subroutine eigenvec(x,eig)  !vl x->jac
    ! calculate eigenvector using LAPACK
    integer::N,LWORK,LDVR,LDVL,LDA,info,j8(1),j7,jj,j88;;
    PARAMETER (N=NVAR, LWORK=108, LDVR=N, LDVL=N, LDA=N)
    character::JOBVL,JOBVR
    real(kind=kind(0.d0))::WR(N),WI(N),WORK(LWORK),VR(LDVL,N),yor(1),yoi(1),VL(LDVR,N),le(LDVR,N),ri(LDVL,N)
    real(kind=kind(0.d0)),intent(out)::eig(N,2)
    real(kind=kind(0.d0)),intent(in)::x(LDA,N) ! JACOBIAN !
    JOBVL='V';JOBVR='V';
    call DGEEV(JOBVL,JOBVR,N,x,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
    !eig(:,1)=WR(:); eig(:,2)=WI(:);
    do j7=1,n
       j8=maxloc(WR); j88=j8(1); yor=WR(j8); yoi=WI(j8); 
       eig(j7,1)=yor(1); eig(j7,2)=yoi(1); WR(j8)=-1000.;
       le(:,j7)=VL(:,j88);ri(:,j7)=VR(:,j88);
    enddo
  end subroutine eigenvec

  subroutine jacsol(x,fjac) ! COMPUTE JACOBIAN with MINPACK
    integer n,ldfjac,iflag,ml,mu,j1
    parameter (n=nvar,ldfjac=nvar)
    real(kind=kind(0.d0))::epsfcn,t
    REAL(kind=kind(0.d0))::fvec(n),wa1(n),wa2(n) !fvec sort de fcn
    REAL(kind=kind(0.d0)),intent(in)::x(n)
    REAL(kind=kind(0.d0)),intent(out)::fjac(n,n)
    external FCN
    ml=n-1; mu=n-1; epsfcn=0.001;
    !do j1=1,n
    !   call fcn_HW(n,x,fvec,iflag)
    !enddo
    call fdjac1(FCN,n,x,fvec,fjac,n,iflag,ml,mu,epsfcn,wa1,wa2)
    ! MINPACK; FDJAC1 estimates an N by N jacobian matrix using forward differences.
  end subroutine jacsol

  subroutine syssol(x) !compute jacobian analytically
    use glob;  implicit none
    integer::n,info,lwa
    parameter (n=nvar, lwa=n*(3*n+13))
    real(kind=kind(0.d0))::wa(lwa),tol,fvec(n)
    real(kind=kind(0.d0))::x(n)
    external FCN
    tol=1d-8;
    call hybrd1(FCN,n,x,fvec,tol,info,wa,lwa)
    ! HYBRD - HYBRD1 HYBRJ HYBRJ1
    ! subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
  end subroutine syssol

  subroutine flux(v,yy,zz)
    !nv is 4 or 5  (if PG) nr=10 ou 11; nm=6;
    !reaction number 
    real(8),intent(in)::v(nvar);
    real(8),intent(out)::yy(15),zz(nreac)
    real(8)::pox1,pox2,ptkt,pnuc,pgly1,pgly2,pgap,pnad
    real(8)::G6P,PG,R5P,F6P,GAP,NADPH,FBP,H2O2
    integer::i,j
    G6P=v(1); PG=v(5); R5P=v(6); F6P=v(2); FBP=v(3); GAP=v(4); NADPH=v(7); H2O2=v(8)
    yy(1)=p0;
    yy(2)=k_gly1*G6P/(1+rr(3)*PG);
    yy(3)=k_gly1*F6P/keq_gly1/(1+rr(3)*PG);
    yy(4)=k_pfk*F6P/(1+rr(4)*H2O2); yy(5)=k_pfb*FBP;
    yy(6)=k_gly3*FBP; yy(7)=k_gly3*GAP/keq_gly3 ! ref?
    yy(8)=k_gapd*GAP/(1+rr(5)*H2O2);
    yy(11)=k_nuc*R5P;
    yy(9)=k_g6pd*G6P*(1-NADPH)/(1.+rr(1)*NADPH);
    yy(10)=k_pgd*PG*(1-NADPH)/(1.+rr(2)*PG);
    yy(12)=k_tkt*R5P**3/(1+(R5P**3+F6P**2*GAP)*ksat);
    yy(13)=k_tkt*F6P**2*GAP/keq_tkt/(1+(R5P**3+F6P**2*GAP)*ksat);
    yy(14)=k_nad*NADPH*H2O2;
    yy(15)=k0
    zz(1)=yy(1);
    zz(2)=yy(2)-yy(3); !gpi
    zz(3)=yy(4); !pfk
    zz(4)=yy(5); !pfb
    zz(5)=yy(6)-yy(7); !ald
    zz(6)=yy(8); !gapd
    zz(7)=yy(9); !g6pd
    zz(8)=yy(10); !pgd
    zz(9)=yy(11); !prp
    zz(10)=yy(12)-yy(13); !tkt
    zz(11)=yy(14); !knad
    zz(12)=yy(15); !knad
  end subroutine flux
  
end module module_algebra

!! PROGRAMME PRINCIPAL !! 
program principal !!REF
  use glob; use module_algebra; implicit none;
  integer,parameter::nstat=10000 !6000 to 60000
  real(8),dimension(15)::fl,fl0
  real(8),dimension(nvar)::cs11
  real(8),dimension(nreac)::flr,flr0,tmp_r1,tmp_r2,flr0b,test,cj11

  real(8),dimension(npar)::par,par0,ranp2,ranp
  real(8)::tmp0(6,7),tmp1(nvar,4),mul_ccc_k2(6,4),cjeq(6),cgeq(6)
  real(8)::dpar(npar,nstat),dx(nvar,nstat)
  real(8)::dfl(nreac,nstat)
  real(8),dimension(nreac,nreac)::fcc,invfcc,fcc2 !,dj,invdj
  real(8)::ccc(nvar,nreac),dccc(nvar,nreac,nstat),dfcc(nreac,nreac,nstat)
  real(8)::covpg(npar),mv1,mv2
  real(8),dimension(nvar)::xx,x0,xs
  integer::nn(8),j1,j2,idid,j3,j4,j5,nend,j6,nnp(6)
  integer::INFO,IPIV(nreac),jk(4),itmp
  real(8)::ti,tf,x1,eta,gam,tmpdv(nvar),pno,tmpvar2,rreff(2)
  integer,parameter::lwork=50
  real(8)::work(lwork),somme,disk(3,30),e1,cj1,cj0,e0

  open(112,file='fig4B.res'); 
  open(113,file='fig4Ba.res'); 
  open(114,file='fig4Bb.res'); 
  open(115,file='fig4Bc.res');
  open(116,file='fig4Bd.res');
  open(117,file='correl1.res');
  open(118,file='correl2.res'); 

  call random_number(ranp);
  do j6=1,4,1
     !! j6=1 noreg
     !! j6=2 r1=5
     !! j6=3 r1=5; r5=2
     !! j6=4 r1=5; r4=2
     jk=0.; write(*,*) j6
     do j3=1,nstat !0
1       call random_number(ranp); par=10**(2D0*(ranp-0.5D0)); par(16:19)=ranp(16:19); rr=0.;
        if (par(1)<0.) goto 2
        if (j6==1) then
           rr=0.; 
           par(2:15)=(/1.,1.,0.5,1.,1.,2.,2.,0.5,0.5,4./3.,0.5,2.,1.,2./)*par(2:15); !!for figure 2
        endif
        if (j6==2) then
           rr((/1,2,3,5/))=(/5.,0.,0.,0./)
           par(2:15)=(/1.,1.,0.5,1.,1.,8.,2.,0.5,0.5,4./3.,0.5,2.,1.,2./)*par(2:15);
           !par(2:15)=(/1.,0.5,1.5,1.,1.,4.,2.,0.02,0.5,4./3.,2.,1.5,3./)*par(2:15);
           !rr(1)=4./3.; par(2:15)=(/1.,1.,0.5,1.,1.,4.,2.,0.5,0.5,4./3.,2.,1.,2./)*par(2:15);
        endif
        if (j6==3) then !j'ai mus prp tres bas !!
           rr((/1,2,3,5/))=(/5.,0.,0.,2./)
           !par(2:15)=(/1.,1.,0.5,1.,10.,4.,2.,0.2,0.5,4./3.,2.,1.,3./)*par(2:15); 
           par(2:15)=(/0.5,1.,0.5,0.5,2.,8.,2.,0.05,0.2,4./3.,1.,2.,1.,2./)*par(2:15); !! AIDE 
        endif
        if (j6==4) then
           rr((/1,4,3,5/))=(/5.,2.,0.,0./) 
           !par(2:15)=(/4.,1.,0.5,1.,2.,4.,2.,0.2,0.5,4./3.,2.,1.,3./)*par(2:15);
           !par(2:15)=(/1.,0.5,1.,1.,2.,4.,2.,0.05,0.5,4./3.,2.,1.,2./)*par(2:15);  !25
           !par(2:15)=(/1.,1.,0.5,1.,4.,8.,2.,0.05,0.5,4./3.,2.,1.,2./)*par(2:15);  !55
           par(2:15)=(/1.,2.,0.5,1.,1.,8.,2.,0.05,0.5,4./3.,1.,2.,1.,2./)*par(2:15);  !55
        endif
        par(1)=1D0; dpar(1:15,j3)=log(par(1:15))/log(1D1); dpar(16:19,j3)=par(16:19); 
        par0=par; call pconv; ! parameter initialization
        ti=0D0; tf=5D3; xx=0.5D0; !initial conditions
        call evolSEU(ti,tf,xx,idid,1.d-6); !numerical integration
        call sys_ox(xx,tmpdv,0D0); if (sum(tmpdv)>1D-4) goto 1 !derivative calculation, check steady state
        x0=xx; dx(:,j3)=x0; 
        call flux(x0,fl0,flr0); !flux calculation 
        test=flr0;
        dfl(:,j3)=flr0; !flux storage
        !if (maxval(fl0)>2D1) goto 1 
        
        !!!!! CONTROL COEFFICIENT !!!!!!!
        j1=12; !index of perturbed parameter
        par(12)=par0(12)*1.02; !perturbaiton of 2\%
        call pconv;  
        xx=x0; call syssol(xx); call sys_ox(xx,tmpdv,0D0); !use of a solver to find steady state concentration
        call flux(xx,fl,flr); !flux calculation 
        do j2=1,nreac
           dfcc(j2,j1,j3)=(flr(j2)-flr0(j2))/(par(j1)-par0(j1))*par0(j1)/flr0(j2)
        enddo
        do j2=1,nvar
           dccc(j2,j1,j3)=(xx(j2)-x0(j2))/(par(j1)-par0(j1))*par0(j1)/x0(j2)
        enddo
        
        if (((x0(7)-0.48)*(0.52-x0(7))>0).and.(dccc(7,12,j3)<-0.04)) then !!! condition of nadph=nadp+
           write(112,*)  dfl(2,j3),dccc(7,12,j3),dccc(2,12,j3),fl0(3)
           if ((dccc(2,12,j3)>0.).and.(fl0(3)<1.)) write(113,*)  dfl(2,j3),dccc(7,12,j3),fl0(2)+fl0(3)
           if ((dccc(2,12,j3)>0.).and.(fl0(3)>1.)) write(114,*)  dfl(2,j3),dccc(7,12,j3),fl0(2)+fl0(3)
           if ((dccc(2,12,j3)>0.).and.(fl0(3)>1.)) then
              write(117,*) dccc((/2,7/),12,j3)
           else
              write(118,*) dccc((/2,7/),12,j3)
           endif
           if ((dccc(2,12,j3)<0.).and.(fl0(3)<1.)) write(115,*)  dfl(2,j3),dccc(7,12,j3),fl0(2)+fl0(3)
           if ((dccc(2,12,j3)<0.).and.(fl0(3)>1.)) write(116,*)  dfl(2,j3),dccc(7,12,j3),fl0(2)+fl0(3)
        endif
        
     enddo !parameter sampling
2    nend=j3 
  enddo !j6

contains

  subroutine pconv !!11 parametres
    p0=par(1); k_gly1=par(2); k_pfk=par(3); k_pfb=par(4); k_gly3=par(5); k_gapd=par(6);
    k_g6pd=par(7); k_pgd=par(8); k_nuc=par(9); k_tkt=par(10); k_nad=par(11);
    k0=par(12); ksat=0.; keq_gly1=par(13); keq_tkt=par(14); keq_gly3=par(15);
  end subroutine pconv

end program principal

!!!!!!!!!!!EXTERNAL SUBROUTINES!!!!!!!!!!!!!!!
subroutine solout (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
  use glob
  implicit none
  integer :: N,NR,LRC,LIC,IRTRN,ipar                               
  real(8) :: X,Y(N),RC(LRC),IC(LIC),rpar,XOLD
  real(4) :: t,val,diff
  real(8) :: F(N)

  !!! stocke l'évolution temporelle dans un tableau
  if (NR.LE.nplot) then
     tplot(NR)=X;
     xplot(:,NR)=Y(:);
  endif
  nmax=NR
!  if ((X>0).and.(Y(nvar).LT.0.5)) then !! --- CRITERE D'ARRET de evol_SEU---!!YO!!
!     Y=-1; IRTRN=-1; return;
 ! endif
  return
end subroutine solout

!!$ ---------------------------------------
!!$   Routines nécessaires à la communication avec 
!!$   l'intégrateur numérique
!!$----------------------------------------

SUBROUTINE FCN(n,x,fvec,iflag)
  use glob; implicit none;
  integer::iflag,n
  real(8)::t,x(n),fvec(n)
  call sys_ox(x,fvec,t) 
end SUBROUTINE FCN

SUBROUTINE FCN_HW(N,t,xHW,F,rpar,ipar)
  use glob
  implicit none
  integer :: n,ipar
  real(8) :: t,xHW(n),f(n),rpar
  call sys_ox(XHW,f,t) 
end SUBROUTINE FCN_HW

subroutine jac_HW(N,X,Y,DFY,LDFY,RPAR,IPAR)
  use glob
  implicit none
  integer :: N,LDFY,LRC,LIC,IRTRN,ipar                                       
  real(8) :: X,Y(N),DFY(LDFY,N),rpar
  !call jacobien(Y,DFY,X)
  return
end subroutine jac_HW

subroutine mass
  return
end subroutine mass
