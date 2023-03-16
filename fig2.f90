!! --gfortran simu.f90 -lseulex -lminpack -llapack; ./a.out  
!! !units micromolar and seconds 

module GLOB
  implicit none
  integer,parameter::nvar=7,npar=18,nreg=4,nmod=7,nreac=11,nm=4,nplot=50000
  real(8),dimension(13,7)::K1(13,7),N1(6,13) ! Kern sans Nadph et directional
  real(8),dimension(nreac,nm)::K2,K2x ! Kern sans Nadph et reversible
  real(8),dimension(nm,nreac)::K2T ! Kern sans Nadph et reversible
  real(8),dimension(nvar,nreac)::N2 ! stochio sans NADPH et reversible
  integer :: nr_max,nmax,numit,nrf
  real(8), dimension(nvar,nplot) :: xplot(nvar,nplot),tplot(nplot)
  real(8) :: pi=acos(-1.d0),l2=log(2D0),l10=log(1D1)
  integer, parameter :: max_val=10**3
  real(8) :: p0,k_nad,k_g6pd,k_pgd,k_tkt,k_gly1,k_gly3,k_nuc,k_gapd,k_pfk,k_pfb,ksat
  real(8) :: keq_tkt,keq_gly1,keq_gly3,rr(nreg),ss(2)
  type indiv  
     sequence
     real(8) :: k(npar),sc
  end type indiv
  type(indiv) :: p,pb
  type(indiv),dimension(10)::pfam

contains

  subroutine matdef
    ! hK,G1,G2,PFK,PFB,A1,A2,GAPd,G6d,PGd,PRP,Tk1,Tk2 (13)
    N1(1,:)=(/1.,-1.,1.,0.,0.,0.,0.,0.,-1.,0.,0.,0.,0./);  !G6P
    N1(2,:)=(/0.,1.,-1.,-1.,1.,0.,0.,0.,0.,0.,0.,2.,-2./);  !F6P
    N1(3,:)=(/0.,0.,0.,1.,-1.,-1.,1.,0.,0.,0.,0.,0.,0./);  !FBP
    N1(4,:)=(/0.,0.,0.,0.,0.,2.,-2.,-1.,0.,0.,0.,1.,-1./);  !GAP
    N1(5,:)=(/0.,0.,0.,0.,0.,0.,0.,0.,1.,-1.,0.,0.,0./);   !PG
    N1(6,:)=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,-1.,-3.,3./);  !R5P
    K1=0.;  K1((/1,2,4,6/),1)=1.; k1(8,1)=2.; ! GLYC
    K1((/1,2,4,6,11,13/),2)=(/1.,1.,0.2,0.2,1.2,0.4/);  ! PROL
    K1((/1,3,5,7,9,10,12/),3)=(/1.,5.,1.,1.,6.,6.,2./); ! CYCLING
    K1((/2,3/),4)=(/1.,1./);   !GPI
    K1((/4,5/),5)=(/1.,1./);   !PFK
    K1((/6,7/),6)=(/1.,1./);   !ALD
    K1((/12,13/),7)=(/1.,1./); !TKT
    ! hK,G1,PFK,PFB,A1,6!GAPd,G6d,PGd,PRP,Tk1 (10)
    N2(1,:)=(/1.,-1.,0.,0.,0.,0.,-1.,0.,0.,0.,0./);  !G6P
    N2(2,:)=(/0.,1.,-1.,1.,0.,0.,0.,0.,0.,2.,0./);  !F6P
    N2(3,:)=(/0.,0.,1.,-1.,-1.,0.,0.,0.,0.,0.,0./);  !FBP
    N2(4,:)=(/0.,0.,0.,0.,2.,-1.,0.,0.,0.,1.,0./);  !GAP
    N2(5,:)=(/0.,0.,0.,0.,0.,0.,1.,-1.,0.,0.,0./);   !PG
    N2(6,:)=(/0.,0.,0.,0.,0.,0.,0.,1.,-1.,-3.,0./);  !R5P
    N2(7,:)=(/0.,0.,0.,0.,0.,0.,1.,1.,0.,0.,-1./)
    K2=0.; K2((/1,2,3,5/),1)=1.; K2(6,1)=2.; !GLYC
    K2((/1,2,3,5,9,10/),2)=(/1.,1.,0.2,0.2,1.2,-0.4/); ! PROL
    K2((/1,2,4,5,7,8,10,11/),3)=(/1.,-5.,1.,-1.,6.,6.,2.,12./); ! CYCL
    K2((/3,4/),4)=(/1.,1./);   !PFK
    K2T=transpose(K2);
   
  end subroutine matdef
  
  subroutine sys_ox(v,dvdt,t) !D,X,Y 
    real(8),dimension(nvar),intent(in)::v
    real(8),dimension(nvar),intent(out)::dvdt
    real(8),intent(in)::t
    real(8)::pox1,pox2,ptkt,pnuc,pgly1,pgly2,pgap,pnad,pgly3
    real(8)::G6P,R5P,GAP,NADPH,F6P,PG,FBP
    G6P=v(1); F6P=v(2); FBP=v(3); GAP=v(4); PG=v(5); R5P=v(6); NADPH=v(7); 
    
    pox1=k_g6pd*G6P*(1-NADPH)/(1.+rr(1)*NADPH);
    pox2=k_pgd*PG*(1-NADPH)/(1.+rr(2)*NADPH);
    ptkt=k_tkt*(R5P**3-F6P**2*GAP/keq_tkt)/(1+(R5P**3+F6P**2*GAP)*ksat);
    pgly1=k_gly1*(G6P-F6P/keq_gly1)/(1.+rr(3)*PG)
    pgly2=(1.+rr(4)*F6P)*k_pfk*F6P-k_pfb*FBP;     
    pgly3=k_gly3*(FBP-GAP/keq_gly3)
    pgap=k_gapd*GAP!/(1+rr(4)*k_nad);
    pnuc=k_nuc*R5P;  pnad=k_nad*NADPH;              
    dvdt(1)=p0-pox1-pgly1;      !G6P
    dvdt(2)=pgly1-pgly2+2*ptkt; !F6P
    dvdt(3)=pgly2-pgly3;        !FBP
    dvdt(4)=2*pgly3-pgap+ptkt;  !gap
    dvdt(5)=pox1-pox2;          !PG
    dvdt(6)=pox2-pnuc-3*ptkt;   !r5P
    dvdt(7)=pox1+pox2-pnad;     !NADPH
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

  subroutine flux(v,yy,zz,ww)
    !nv is 4 or 5  (if PG) nr=10 ou 11; nm=6;
    !reaction number 
    real(8),intent(in)::v(nvar);
    real(8),intent(out)::yy(14),zz(nmod),ww(nm)
    real(8)::pox1,pox2,ptkt,pnuc,pgly1,pgly2,pgap,pnad,yyred(nmod),yyred2(nm)
    real(8),dimension(nmod,nmod)::efm,efminv,efmm
    real(8),dimension(nm,nm)::efm4,efminv4,efmm4
    real(8)::G6P,PG,R5P,F6P,GAP,NADPH,FBP
    integer::i,j
    G6P=v(1); PG=v(5); R5P=v(6); F6P=v(2); FBP=v(3); GAP=v(4); NADPH=v(7);

    yy(1)=p0; yy(2)=k_gly1*G6P/(1+rr(3)*PG); yy(3)=k_gly1*F6P/keq_gly1/(1+rr(3)*PG);
    yy(4)=(1+F6P*rr(4))*k_pfk*F6P; yy(5)=k_pfb*FBP;
    yy(6)=k_gly3*FBP; yy(7)=k_gly3*GAP/keq_gly3 ! ref?
    yy(8)=k_gapd*GAP !/(1+rr(4)*k_nad);
    yy(11)=k_nuc*R5P;
    yy(9)=k_g6pd*G6P*(1-NADPH)/(1.+rr(1)*NADPH);
    yy(10)=k_pgd*PG*(1-NADPH)/(1.+rr(2)*NADPH);
    yy(14)=k_nad*NADPH;
    yy(12)=k_tkt*R5P**3/(1+(R5P**3+F6P**2*GAP)*ksat);
    yy(13)=k_tkt*F6P**2*GAP/keq_tkt/(1+(R5P**3+F6P**2*GAP)*ksat);;
    !2,3; 4,5; 6,7; 12,13
    !decomposition yy j'enelver phi_ox7, phi_gap6, phi_nuc8 & phi_gly2-5
    yyred=yy((/1,2,3,4,6,12,13/));
    !hk,gpi+,gpi-,pfk,ald,tkt+,tkt-
    efm(:,1)=(/1.,1.,0.,1.,1.,0.,0./);    !glycolytic mode
    efm(:,2)=(/1.,1.,0.,0.2,0.2,0.,0.4/); !prolif mode    
    efm(:,3)=(/1.,0.,5.,0.,0.,2.,0./);    !cycling
    efm(:,4)=(/0.,1.,1.,0.,0.,0.,0./);    ! gly1 futile cycle
    efm(:,5)=(/0.,0.,0.,1.,0.,0.,0./);    ! gly2 futile cycle
    efm(:,6)=(/0.,0.,0.,0.,1.,0.,0./);    ! gly3 futile cycle
    efm(:,7)=(/0.,0.,0.,0.,0.,1.,1./);    ! nonox futile cycle
    efminv=inv(efm); zz=matmul(efminv,yyred);

    yyred2(1)=yy(1);
    yyred2(2)=yy(2)-yy(3);
    yyred2(3)=yy(4);
    yyred2(4)=yy(12)-yy(13);

    !hk,gpi,pfk/ald,tkt
    efm4(:,1)=(/1.,1.,1.,0./);    !glycolytic mode
    efm4(:,2)=(/1.,1.,0.2,-0.4/); !prolif mode    
    efm4(:,3)=(/1.,-5.,0.,2./);   !cycling
    efm4(:,4)=(/0.,0.,1.,0./);    !gly2 futile cycle
    efminv4=inv(efm4);
    ww=matmul(efminv4,yyred2)/yy(1);
    
    zz=matmul(efminv,yyred)/yy(1);
  end subroutine flux
  
  subroutine fluxrev(x,y) 
    real(8),dimension(14),intent(in)::x
    real(8),dimension(11),intent(out)::y
    y(1)=x(1);
    y(2)=x(2)-x(3); !gpi
    y(3)=x(4); !pfk
    y(4)=x(5); !pfb
    y(5)=x(6)-x(7); !ald
    y(6)=x(8); !gapd
    y(7)=x(9); !g6pd
    y(8)=x(10); !pgd
    y(9)=x(11); !prp
    y(10)=x(12)-x(13); !tkt
    y(11)=x(14); !knad
  end subroutine fluxrev
end module module_algebra

!! PROGRAMME PRINCIPAL !! 
program principal !!REF
  use glob; use module_algebra; implicit none;
  integer,parameter::nstat=10000 
  real(8),dimension(14)::fl,fl0
  real(8),dimension(nreac)::flr,flr0,tmp_r1,tmp_r2,flr0b,test
  real(8),dimension(npar)::par,par0,ranp2,ranp
  real(8),dimension(nmod)::yefm0(nmod),yefm(nmod)
  real(8),dimension(nm)::wefm0,wefm
  real(8)::tmp0(6,7),tmp1(nvar,4),mul_ccc_k2(6,4)
  real(8),dimension(nm,nreac+1)::mean,var
  real(8)::wp(nm,nreac+1,5),cum(nm,nreac+1,200),wpfcc(nreac,5),cumfcc(nreac,200)
  real(8)::wpcc(nvar,nreac,5),cumcc(nvar,nreac,200),wpddcc(nvar,5),cumddcc(nvar,200)
  real(8),dimension(nvar,nreac,200)::disccc   !! distributions
  real(8),dimension(nreac,200)::disfcc   !! distributions
  real(8),dimension(3,2,nreac,200)::dd  !!3cas, 2mode,
  real(8),dimension(2,nvar,200)::ddcc  !!3cas, 2mode,
  real(8)::cov(nm,nreac,nreac),cov1(nreac,npar+3+nreac),covx(3,nreac,nvar)
  real(8)::covp(3,nreac,npar),covw(3,nreac,nm),covjr(3,nreac,4),covcc(3,nreac,nvar)
  real(8)::zm(2,npar+3+nreac),zmdo(4,2,npar+3+nreac) !covariance
  real(8)::dpar(npar,nstat),dx(nvar,nstat),dwefm(5,nstat)
  real(8)::dfl(nreac+1,nstat),dflr(4,nstat),dwefmdo(8,5,nstat) !storage
  real(8)::dpardo(8,npar,nstat),dfldo(8,nreac+1,nstat),dflrdo(8,4,nstat) !storage
  real(8),dimension(nreac,nreac)::fcc,invfcc,fcc2 !,dj,invdj
  real(8)::ccc(nvar,nreac),dccc(nvar,nreac,nstat),dfcc(nreac,nreac,nstat)
  real(8)::covpg(npar),mv1,mv2
  !real(8),dimension(nreac,nm)::tmp1_nmnr,tmp2_nmnr
  real(8),dimension(nm,nm)::dw
  real(8),dimension(nvar)::xx,x0,xs
  integer::nn(8),j1,j2,idid,j3,j4,j5,nend,j6,nnp(6)
  integer::INFO,IPIV(nreac),jdis(nm,nreac),itmp,jtr(10),jdis2(nvar,nreac),jk(4)
  real(8)::ti,tf,x1,eta,gam,tmpvar,tmpmean,tmpdv(nvar),pno,tmpmean2,tmpvar2,rreff(2)
  integer,parameter::lwork=50
  real(8)::work(lwork),dndg

  open(10,file='cx_fig2.res');   !sampling
  open(11,file='cw_fig2.res');   !sampling
  open(12,file='cxmv_fig2.res'); !meanvar for cx
  open(13,file='wp_cc1.res');    !whiskecrplot for cs
  open(73,file='wp_c.res');      !whiskecrplot for cs
  open(173,file='wp_cx.res');    !whiskecrplot for cs
  open(79,file='wp_j.res');        !whiskerplot for cj only g6pd 
  open(179,file='wp_jx.res');      !minimum and maximu (grey on figure) 
  !open(21,file='all.res');
  open(30,file='param504208.dat'); !black  subdomain in figure 2 
  open(31,file='param851005.dat'); !red    ""
  open(32,file='param108505.dat'); !orange ""
  open(33,file='param473817.dat'); !blue   "" 
  open(34,file='param601030.dat'); !yellow ""
  open(35,file='param106030.dat'); !green  ""
  !these files have been generated by a random sampling procedure
  !paramXXYYZZ where XX corresponds to J_GAP/J_HK and YY to J_PPP/6J_HK
  nnp=(/4517,6198,5183,1603,1625,771/);

  call matdef; !stochio, kernel
  if (0>1) then
     tmp0=matmul(N1,K1); tmp1=matmul(N2,K2);
     write(*,'(50f4.1)') tmp1  !write(*,'(100f6.3)') matmul(K2T,K2)
  endif

  mean=0D0; var=0D0;
  do j6=1,6,1
     jk=0.; 
     cum=0D0; disccc=0.; nn=0; cumddcc=0D0; dd=0.; ddcc=0.; disfcc=0.; 
     do j3=1,nnp(j6)
1       call random_number(ranp); ranp2=2D0*(ranp-0.5D0); par=10**(ranp2); rr=0.; par(15:18)=ranp(15:18);
        read(29+j6,*) par(1:18);
        if (par(1)<0.) goto 2
        par(1)=1D0; dpar(1:14,j3)=log(par(1:14))/log(1D1); dpar(15:18,j3)=par(15:18); 
        par0=par; call pconv; ti=0D0; tf=5D3; xx=0.5D0; call evolSEU(ti,tf,xx,idid,1.d-6); 
        call sys_ox(xx,tmpdv,0D0); if (sum(tmpdv)>1D-4) goto 1
        x0=xx; dx(:,j3)=x0;
        call flux(x0,fl0,yefm0,wefm0); call fluxrev(fl0,flr0); !flr0b=matmul(K2,wefm0);
        test=flr0; K2x=K2;
        dfl(1:nreac-3,j3)=flr0((/2,3,4,5,6,7,9,10/));
        dfl(nreac-2,j3)=fl0(2)+fl0(3);  dfl(nreac-1,j3)=fl0(4)+fl0(5);
        dfl(nreac,j3)=fl0(6)+fl0(7);  dfl(nreac+1,j3)=fl0(12)+fl0(13);
        dfl(:,j3)=min(dfl(:,j3),5.); dflr(:,j3)=dfl(nreac-2:nreac+1,j3);
        dwefm(1:3,j3)=wefm0(1:3); dwefm(4,j3)=1./6.*wefm0(3)+5./6.*wefm0(1); dwefm(5,j3)=1./6.*wefm0(3)+5./6.*wefm0(2);
        call sgels ('N', nreac, nm, 1, K2x, nreac , test, nreac, work, lwork, INFO)
        !A(m=nr,n=nm) X(n=nm)=B(n=nr); trans=N ou T
        if (wefm0(3)>1./6.) then
           jk(3)=jk(3)+1.
        else
           jk(4)=jk(4)+1.
        endif
        do j1=1,nreac,1 !calculs 
           par=par0; par(j1)=par0(j1)*1.05D0; call pconv;  
           !ti=0.; tf=5D2; xx=x0; call evolSEU(ti,tf,xx,idid,1.d-6);
           xx=x0;call syssol(xx); call sys_ox(xx,tmpdv,0D0); !if (sum(tmpdv)>1D-2) goto 1
           call flux(xx,fl,yefm,wefm); call fluxrev(fl,flr);       
           do j2=1,nreac
              dfcc(j2,j1,j3)=(flr(j2)-flr0(j2))/(par(j1)-par0(j1))*par0(j1)/flr0(j2)
              if (j2==7) then
                 itmp=int(25D0*(dfcc(7,j1,j3)+4)+1); itmp=min(itmp,200); itmp=max(1,itmp)
                 disfcc(j1,itmp)= disfcc(j1,itmp)+1.
              endif
           enddo
           do j2=1,nvar
              dccc(j2,j1,j3)=(xx(j2)-x0(j2))/(par(j1)-par0(j1))*par0(j1)/x0(j2)
              jdis2(j2,j1)=int(25D0*(dccc(j2,j1,j3)+4)+1); !200
              if (jdis2(j2,j1)<1) jdis2(j2,j1)=1
              if (jdis2(j2,j1)>200) jdis2(j2,j1)=200
              itmp=jdis2(j2,j1)
              disccc(j2,j1,itmp)=disccc(j2,j1,itmp)+1D0
              write(10,*) real(j1)+0.1*real(j2-3),dccc(j2,j1,j3) ! cx_fig2.res
              if (j1==11) then
                 if (wefm0(3)>1./6.)  then
                    ddcc(1,j2,itmp)=ddcc(1,j2,itmp)+1D0
                 else
                    ddcc(2,j2,itmp)=ddcc(1,j2,itmp)+1D0
                 endif
              endif
              if (j2==5) then
                 if (wefm0(3)>1./6.)  then
                    dd(1,1,j1,itmp)=dd(1,1,j1,itmp)+1D0
                 else
                    dd(1,2,j1,itmp)=dd(1,2,j1,itmp)+1D0
                 endif
              endif
              if (j2==7) then
                 if (wefm0(3)>1./6.)  then
                    dd(2,1,j1,itmp)=dd(2,1,j1,itmp)+1D0
                 else
                    dd(2,2,j1,itmp)=dd(2,2,j1,itmp)+1D0
                 endif
              endif
           enddo

        enddo
     enddo 

2    nend=j3 
     write(*,*) j6,nend

!!!! %%%%%%%%%%%% WHISKER PLOT FOR FCC %%%%%%%%%%%%%%% OUI
     do j3=1,nreac
        cumfcc(j3,1)=disfcc(j3,1)/real(nend);
        wpfcc=-2D1; 
        do j1=2,200 !trace des wp.dat
           cumfcc(j3,j1)=disfcc(j3,j1)/real(nend)+cumfcc(j3,j1-1);
           if ((cumfcc(j3,j1)>0.1).and.(wpfcc(j3,4)<-1D1)) wpfcc(j3,4)=((real(j1)-0.5)/25D0-4D0)
           if ((cumfcc(j3,j1)>0.25).and.(wpfcc(j3,2)<-1D1)) wpfcc(j3,2)=((real(j1)-0.5)/25D0-4D0)
           if ((cumfcc(j3,j1)>0.5).and.(wpfcc(j3,1)<-1D1)) wpfcc(j3,1)=((real(j1)-0.5)/25D0-4D0)
           if ((cumfcc(j3,j1)>0.75).and.(wpfcc(j3,3)<-1D1)) wpfcc(j3,3)=((real(j1)-0.5)/25D0-4D0)
           if ((cumfcc(j3,j1)>0.9).and.(wpfcc(j3,5)<-1D1)) wpfcc(j3,5)=((real(j1)-0.5)/25D0-4D0)
        enddo
        write(79,'(12f8.3)') real(j3)+0.1*real(j6),wpfcc(j3,:);
        mv1=minval(dfcc(7,j3,1:nend)); mv2=maxval(dfcc(7,j3,1:nend));
        write(179,*) real(j3)+0.1*real(j6),(mv1+mv2)/2.,mv2-(mv1+mv2)/2.
     enddo
     write(79,*)  ' ';

!!!! %%%%%%%%%%%% WHISKER PLOT FOR CS (wpcc) %%%%%%%%%%%%%%%
     wpcc=-2D1;
     do j2=1,nvar
        do j3=1,nreac
           cumcc(j2,j3,1)=disccc(j2,j3,1)/real(nend);
           do j1=2,200 !trace des wp.dat
              cumcc(j2,j3,j1)=disccc(j2,j3,j1)/real(nend)+cumcc(j2,j3,j1-1);
              if ((cumcc(j2,j3,j1)>0.1).and.(wpcc(j2,j3,4)<-1D1)) wpcc(j2,j3,4)=((real(j1)-0.5)/25D0-4D0)
              if ((cumcc(j2,j3,j1)>0.25).and.(wpcc(j2,j3,2)<-1D1)) wpcc(j2,j3,2)=((real(j1)-0.5)/25D0-4D0)
              if ((cumcc(j2,j3,j1)>0.5).and.(wpcc(j2,j3,1)<-1D1)) wpcc(j2,j3,1)=((real(j1)-0.5)/25D0-4D0)
              if ((cumcc(j2,j3,j1)>0.75).and.(wpcc(j2,j3,3)<-1D1)) wpcc(j2,j3,3)=((real(j1)-0.5)/25D0-4D0)
              if ((cumcc(j2,j3,j1)>0.9).and.(wpcc(j2,j3,5)<-1D1)) wpcc(j2,j3,5)=((real(j1)-0.5)/25D0-4D0)
           enddo
           if (j3==11) then
              mv1=minval(dccc(j2,j3,1:nend)); mv2=maxval(dccc(j2,j3,1:nend));
              write(13,*) real(j2)+0.1*real(j6),wpcc(j2,j3,:) !for nadph perturbation
              write(173,*) real(j2)+0.1*real(j6),(mv1+mv2)/2.,mv2-(mv1+mv2)/2.
           endif
           if (j2==5) then ! Whisker plot for PG
              mv1=minval(dccc(j2,j3,1:nend)); mv2=maxval(dccc(j2,j3,1:nend));
              write(73,*) real(j3)+0.1*real(j6),wpcc(j2,j3,:) !for PG changes
              !write(173,*) real(j3)+0.1*real(j6),(mv1+mv2)/2.,mv2-(mv1+mv2)/2.
           endif
        enddo
     enddo
     write(13,*) ' '; write(73,*) ' ';
enddo

contains

  subroutine pconv !!11 parametres
    p0=par(1); k_gly1=par(2); k_pfk=par(3); k_pfb=par(4); k_gly3=par(5); k_gapd=par(6);
    k_g6pd=par(7); k_pgd=par(8); k_nuc=par(9); k_tkt=par(10); k_nad=par(11); ksat=0.;
    !keq_gly1=1D0; keq_tkt=1D0; keq_gly3=1D0;
    keq_gly1=par(12); keq_tkt=par(13); keq_gly3=par(14);
  end subroutine pconv

end program principal

  
!!!!!!!!!!!EXTERNAL SUBROUTINES!!!!!!!!!!!!!!!

!!!!%%%SOLOUT

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
  call sys_ox(XHW,f,t) !! changer
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
