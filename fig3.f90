!! Compilation and execution command :
!! > gfortran simu.f90 -lseulex2 -lminpack -llapack; ./a.out  
!! Units: micromolar and seconds
!! 18 fevrier 2021
module GLOB
  implicit none
  integer,parameter::nvar=7,npar=18,nreg=5,nmod=7,nreac=11,nm=4,nplot=5000
  real(8),dimension(nm,nreac)::K2T ! Kern sans Nadph et reversible
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
  
  subroutine sys_ox(v,dvdt,t) !D,X,Y 
    real(8),dimension(nvar),intent(in)::v
    real(8),dimension(nvar),intent(out)::dvdt
    real(8),intent(in)::t
    real(8)::pox1,pox2,ptkt,pnuc,pgly1,pgly2,pgap,pnad,pgly3
    real(8)::G6P,R5P,GAP,NADPH,F6P,PG,FBP
    G6P=v(1); F6P=v(2); FBP=v(3); GAP=v(4); PG=v(5); R5P=v(6); NADPH=v(7); 
    
    pox1=k_g6pd*G6P*(1-NADPH)/(1.+rr(1)*NADPH);
    pox2=k_pgd*PG*(1-NADPH)/(1.+rr(2)*NADPH+rr(4)*PG);
    ptkt=k_tkt*(R5P**3-F6P**2*GAP/keq_tkt)/(1+(R5P**3+F6P**2*GAP)*ksat);
    pgly1=k_gly1*(G6P-F6P/keq_gly1)/(1.+rr(3)*PG)
    pgly2=(1.)*k_pfk*F6P-k_pfb*FBP;     
    pgly3=k_gly3*(FBP-GAP/keq_gly3)
    pgap=k_gapd*GAP/(1+rr(5)*k_nad);
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

  !! compute flux values from variable values
  subroutine flux(v,yy)
    !reaction number 
    real(8),intent(in)::v(nvar);
    real(8),intent(out)::yy(14)
    real(8)::pox1,pox2,ptkt,pnuc,pgly1,pgly2,pgap,pnad,yyred(nmod),yyred2(nm)
    real(8)::G6P,PG,R5P,F6P,GAP,NADPH,FBP
    integer::i,j
    G6P=v(1); PG=v(5); R5P=v(6); F6P=v(2); FBP=v(3); GAP=v(4); NADPH=v(7);

    yy(1)=p0;
    yy(2)=k_gly1*G6P/(1+rr(3)*PG);
    yy(3)=k_gly1*F6P/keq_gly1/(1+rr(3)*PG);
    yy(4)=(1+F6P*rr(4))*k_pfk*F6P; yy(5)=k_pfb*FBP;
    yy(6)=k_gly3*FBP; yy(7)=k_gly3*GAP/keq_gly3 ! ref?
    yy(8)=k_gapd*GAP/(1+rr(5)*k_nad);
    yy(11)=k_nuc*R5P;
    yy(9)=k_g6pd*G6P*(1-NADPH)/(1.+rr(1)*NADPH);
    yy(10)=k_pgd*PG*(1-NADPH)/(1.+rr(2)*NADPH+rr(4)*PG);
    yy(14)=k_nad*NADPH;
    yy(12)=k_tkt*R5P**3/(1+(R5P**3+F6P**2*GAP)*ksat);
    yy(13)=k_tkt*F6P**2*GAP/keq_tkt/(1+(R5P**3+F6P**2*GAP)*ksat);;
    !2,3; 4,5; 6,7; 12,13
    !decomposition yy j'enelver phi_ox7, phi_gap6, phi_nuc8 & phi_gly2-5
    yyred=yy((/1,2,3,4,6,12,13/));
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
  integer,parameter::nstat=30000 !sampling size
  real(8),dimension(14)::fl,fl0
  real(8),dimension(nvar)::cs11
  real(8),dimension(nreac)::flr,flr0,cj11

  real(8),dimension(npar)::par,par0,ranp
  real(8)::cjeq(6),cgeq(6)
  real(8)::dpar(npar,nstat),dx(nvar,nstat)
  real(8)::dfl(nreac+1,nstat),dflr(4,nstat) 
  real(8),dimension(nreac,nreac)::fcc,fcc2 !,dj,invdj
  real(8)::ccc(nvar,nreac),dccc(nvar,nreac,nstat),dfcc(nreac,nreac,nstat)
  real(8),dimension(nvar)::xx,x0,xs
  integer::j1,j2,idid,j3,nend,j6
  integer::INFO,IPIV(nreac)
  real(8)::ti,tf,x1,tmpdv(nvar)
  integer,parameter::lwork=50
  real(8)::work(lwork),dndg,e1,cj1,cj0,e0
  open(11,file='fig3a.res');
  open(12,file='fig3b.res');
  open(13,file='fig3c.res');
  open(14,file='fig3d.res');
  open(15,file='fig3e.res');
  open(16,file='fig3f.res');

  open(21,file='fig3_6pg_a.res');
  open(22,file='fig3_6pg_b.res');
  open(23,file='fig3_6pg_c.res');
  open(24,file='fig3_6pg_d.res');
  open(25,file='fig3_6pg_e.res');
  open(26,file='fig3_6pg_f.res');

  open(31,file='fig3_g6p_a.res');
  open(32,file='fig3_g6p_b.res');
  open(33,file='fig3_g6p_c.res');
  open(34,file='fig3_g6p_d.res');
  open(35,file='fig3_g6p_e.res');
  open(36,file='fig3_g6p_f.res');
  open(40,file='fig3F_corr.res');
  open(41,file='fig3F_corr_or.res');
  open(42,file='fig3F_corr_re.res');
  open(43,file='fig3F_corr_bl.res');

  !j6=1 : no regulation
  !j6=2 : r1=2 
  !j6=3 : r3=2
  !j6=4 : r3=2

  do j6=1,6,1
     write(*,*) j6
     do j3=1,nstat !XX
1       call random_number(ranp); !ranp2 entre 1 et -1
        par=10**(4.D0*(ranp-0.5D0)); !par entre 10^-2 and 10^2
        par(15:18)=ranp(15:18); !entre 0 et 1 ;
        rr(4)=0.;
        if (par(1)<0.) goto 2
        if (j6==1) then
           rr=0.; 
           par(2:14)=(/1.,1.,0.5,1.,1.,2.,2.,0.5,0.5,4./3.,2.,1.,2./)*par(2:14); !!for figure 2
        endif
        if (j6==2) then
           rr((/1,2,3,5/))=(/2.,0.,0.,0./) !*par(15:18);  !si tirage au sort de ri
           par(2:14)=(/1.,1.,0.5,1.,1.,4.,2.,0.5,0.5,4./3.,2.,1.,2./)*par(2:14);
           !par(2:14)=(/1.,0.5,1.5,1.,1.,4.,2.,0.02,0.5,4./3.,2.,1.5,3./)*par(2:14);
        endif

        if (j6==3) then 
           rr((/1,2,3,5/))=(/0.,0.,2.,0./) 
           par(2:14)=(/2.,1.,0.5,1.,1.,2.,2.,0.5,0.5,4./3.,2.,1.,2./)*par(2:14);
           !par(2:14)=(/2.5,1.,0.5,1.,1.,2.,4.,0.5,0.5,4./3.,2.,1.,2./)*par(2:14); ???
        endif

        if (j6==4) then
           rr((/1,2,3,5/))=(/2.,0.,2.,0./)
           !rr(1)=4./3.; rr(3)=3.D0; 
           par(2:14)=(/2.,1.,0.5,1.,1.,4.,2.,0.2,0.5,4./3.,2.,1.,2./)*par(2:14); ! ??? 
           !par(2:14)=(/2.5,1.,0.5,1.,1.,2.,4.,0.5,0.5,4./3.,2.,1.,3./)*par(2:14);
        endif
        if (j6==5) then !j'ai mus prp tres bas !!
           rr((/1,2,3,5/))=(/2.,2.,2.,0./) 
           !par(2:14)=(/1.,1.,0.5,1.,10.,4.,2.,0.2,0.5,4./3.,2.,1.,3./)*par(2:14); 
           par(2:14)=(/2.,1.,0.5,1.,1.,4.,4.,0.5,0.5,4./3.,2.,1.,2./)*par(2:14); ! ??? 
        endif

        par(1)=1D0; par(2)=1D0; dpar(1:14,j3)=log(par(1:14))/log(1D1); dpar(15:18,j3)=par(15:18);
        par0=par; call pconv; !initialisation des parametres
        ti=0D0; tf=5D3; xx=0.5D0; ! condition initiales
        call evolSEU(ti,tf,xx,idid,1.d-6);  !simulation de l'etat de reference
        call sys_ox(xx,tmpdv,0D0); !calcul de dvdt
        if (sum(tmpdv)>1D-5) goto 1 ! check
        x0=xx; dx(:,j3)=x0; !stock value of concentration
        call flux(x0,fl0); !calcul les flux
        !write(*,'(10f10.5)') fl0
        call fluxrev(fl0,flr0);  !calcul les flux net
        dfl(1:nreac-3,j3)=flr0((/2,3,4,5,6,7,9,10/));
        dfl(nreac-2,j3)=fl0(2)+fl0(3);  dfl(nreac-1,j3)=fl0(4)+fl0(5);
        dfl(nreac,j3)=fl0(6)+fl0(7);  dfl(nreac+1,j3)=fl0(12)+fl0(13);
        if (maxval(dfl(:,j3))>2D1) goto 1  !refuse les valeurs de flux net trop grandes 
        dflr(:,j3)=dfl(nreac-2:nreac+1,j3);

!!!!!!!!!!!!!!!!!%%%%%%%% CONTROL COEFFICIENT
        do j1=1,nreac,1 
           par=par0; par(j1)=par0(j1)*1.01D0; call pconv;  
           !ti=0.; tf=5D2; xx=x0; call evolSEU(ti,tf,xx,idid,1.d-6);
           xx=x0;
           call syssol(xx);
           call sys_ox(xx,tmpdv,0D0); !if (sum(tmpdv)>1D-2) goto 1
           call flux(xx,fl); 
           call fluxrev(fl,flr);
           do j2=1,nreac !calcyl des dfcc
              dfcc(j2,j1,j3)=(flr(j2)-flr0(j2))/(par(j1)-par0(j1))*par0(j1)/flr0(j2)
           enddo
           do j2=1,nvar !calcul des dccc
              dccc(j2,j1,j3)=(xx(j2)-x0(j2))/(par(j1)-par0(j1))*par0(j1)/x0(j2)
           enddo
        enddo

        cj11=dfcc(:,11,j3);  cs11=dccc(:,11,j3);
        call equations(x0,fl0,cs11,cj11,cjeq,cgeq); !whats for ?
        e1=-x0(7)*(1+rr(1))/(1-x0(7))/(1+rr(1)*x0(7))
        e0=-x0(7)/(1-x0(7));
        cj0=(-e0+dccc(5,11,j3))/(1-e0)
        cj1=(-e1+dccc(1,11,j3))/(1-e1)
        
        if ((dccc(5,11,j3)<0).and.(dccc(1,11,j3)<0)) then
           write(10+j6,*) dccc(5,11,j3),dfcc(7,11,j3)/x0(7),dccc(1,11,j3)
           write(43,*) dfl(6,j3),dfcc(7,11,j3)/x0(7),fl0(2) !JPPP
        endif
        if ((dccc(5,11,j3)>0).and.(dccc(1,11,j3)<0)) then
           write(20+j6,*) dccc(5,11,j3),dfcc(7,11,j3)/x0(7),dccc(1,11,j3)
           write(42,*) dfl(6,j3),dfcc(7,11,j3)/x0(7),fl0(2) !JPPP
        endif
        if ((dccc(5,11,j3)>0).and.(dccc(1,11,j3)>0)) then
           write(30+j6,*) dccc(5,11,j3),dfcc(7,11,j3)/x0(7),dccc(1,11,j3)
           write(41,*) dfl(6,j3),dfcc(7,11,j3)/x0(7),fl0(2) !JPPP
        endif
        write(40,*) dfl(6,j3),dfcc(7,11,j3)/x0(7),fl0(2) !JPPP
     enddo !end nstat
2    nend=j3 !IMPORTANT
  enddo

contains

  subroutine pconv !11 parameters of the model
    p0=par(1); k_gly1=par(2); k_pfk=par(3); k_pfb=par(4); k_gly3=par(5); k_gapd=par(6);
    k_g6pd=par(7); k_pgd=par(8); k_nuc=par(9); k_tkt=par(10); k_nad=par(11); ksat=0.;
    keq_gly1=par(12); keq_tkt=par(13); keq_gly3=par(14);
  end subroutine pconv

  subroutine equations(vx,vj,vcs,vcj,vk,vg)
    real(8),intent(in)::vx(nvar),vj(14),vcs(nvar),vcj(nreac)
    real(8),intent(out)::vk(6),vg(6)
    real(8)::e2,e0,e1,e3,kt,jt
    e0=-vx(7)/(1-vx(7))
    e1=-vx(7)*(1+rr(1))/(1-vx(7))/(1+rr(1)*vx(7))
    e2=-rr(3)*vx(5)/(1+rr(3)*vx(5))
    kt=k_gly1/(1+rr(3)*vx(5))
    e3=vx(1)*kt; jt=vj(2)-vj(3)
    vk(1)=-e1*e3; !K1
    vk(2)=-e0*e2*(vj(2)-vj(3))*(1+rr(2)*vx(5)); !K2
    vk(3)=vcs(2)*vj(3); !K3
    vk(4)=vj(9)+fl(2)*(1-e1)+(vj(2)-vj(3))*e2*(1-e0)*(1+rr(2)*vx(5));
    vk(5)=sum(vk(1:3))/vk(4)
    vg=0D0;
    vg(1)=e1*vj(9); !0
    vg(2)=e2*(vj(2)-vj(3))*(e1-e0); !>0
    vg(3)=vcs(2)*vj(3)*(1-e1) !<0
    vg(4)=vk(4) 
    vg(5)=sum(vg(1:3))/vk(4)
    vk(6)=-e1*(jt+fl(3)*(1-vcs(2)/e1))/(1-e1*jt+fl(3)*(1-e1))
  end subroutine equations
end program principal
  
!!!!!!!!!!! EXTERNAL SUBROUTINES !!!!!!!!!!!!!!!

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
