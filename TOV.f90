Program TOV_solver
    ! This code works for solving the Tolman-Oppenheimer-Volkov (TOV) equation
    !      & tidal deformability by 4th order adaptive-step Runge-Kutta method.
    !       
    !                                      Written by Chencan Wang, 04/10/2021  
    !
    ! Running this code needs provided EOS table written in following order:
    ! |-------------------------------------------------------------------!
    ! | num. dens. | ene. dens. [MeV/fm^3] | pre [MeV/fm^3] | cs^2 (dP/de)|
    ! |-------------------------------------------------------------------! 
    !    Be sure that pre is monotonous to ene and the num. dens.
    !    
    !   
    !
    !
    !   Physical constants (SI)
    !   G (gravitational constant) = 6.67259 x 10^-11 m^3/(s^2  kg)
    !   c (speed of light) = 299792.458 km 
    !   1 Mo (solar mass) = 1.98892 x 10^30 kg 
    !   1 MeV/c^2 = 1.7826619069 x 10^-30 kg =  0.8962964357 x 10^-60 Mo
    !
    !   ** Derived quantities (set c = 1 )
    !      1 Mev/fm^3 =  8.962964357 x 10^-7 Mo/km^3
    !      1 G/c^2 = 1.47662444342 km/Mo  
    !
    implicit none 
    !==================================================! 
    !        Input EOS given by <filename>             !
    !==================================================!
    character(len=30) :: filename = 'data.dat' 
    ! Show explicit neutron star profile if providen .TRUE..
    logical :: Details = .False.!.True. 
    !   In general, no further modification is needed 
    !       for the contents below.
    !==================================================!
    !
    !
    !
    integer       ::  nmax,  neos,  n, np, nwrite 
    !-----------------------------------------------!
    !                 Parameters                    !  
    !-----------------------------------------------!
    !  cdRmax --- max radius step [km]  
    !  cR1M --- G/c^2(in unit [km/Mo], Mo the solar mass)
    !  ceos2ns --- transfer the unit of density [MeV/fm^3] into [Mo/km^3]
    real(kind=16) :: chmax, cRmax, cpi, cR1M, ceos2ns, ceps, chmin
    !  NS is a abbreviation of 'neutron star'.
    parameter(      nmax=1000000,                          &
                    cpi=3.1415926535898d0,                 &    
                    cR1M=1.47662444342d0,                  &
                    ceos2ns=4.0*cpi*8.962964357e-7,        &
    !======================================================!        
    !               Adjustable parameters                  ! 
                    chmax=0.01d0, cRmax = 20., chmin =1d-5,&
    !               The ending pressure should be in accordance with the given EOS.
                    np=100,  ceps = 1d-6)  
    !======================================================!
    ! Initial pressure 
    real(kind=16)   :: pres(np), pre_out  
    !             1     2     3     4     5     6     7     8     9    10
    data pres  /640,  620,  600,  580,  560,  540,  520,  500,  480,  460,   &
                440,  420,  400,  380,  360,  340,  320,  300,  280,  260,   &
                240,  220,  200,  180,  160,  140,  120,  100,   90,   80,   &
                 70,   60,   58,   56,   54,   52,   50,   48,   46,   44,   &
                 42,   40,   35,   30,   25,   20,   18,   16,   14,   12,   &
                 10,    9,    8,    7,    6,    5,   4,    3,    2,    1,    &
                  0.5    /
    !-----------------------------------------------!
    !                NS properties                  ! 
    !-----------------------------------------------!
    !  Rns[Radius], Mns[Mass], Pns[Pressure], Ens[Energy], Dns[Number], Cns[SpeedofSound], Yns[function y(r)]. 
    real(kind=16) :: Rns(nmax), Mns(nmax), Pns(nmax), Ens(nmax), Dns(nmax), Cns(nmax), Yns(nmax)
    !  The Providen EOS,                   
    real(kind=16), allocatable :: rho(:), ene(:), pre(:), cs2(:)
    character(len=10) :: cpres
    !-----------------------------------------!
    !      Embedded Runge-Kutta method        !
    !       with Cash-Karp Parameters         ! 
    !-----------------------------------------!
    real(kind=16) :: a2, a3, a4, a5, a6
    real(kind=16) :: b21, b32, b31, b41, b42, b43, b54, b53, b52, b51, &
                     b65, b64, b63, b62, b61, c1, c2, c3, c4, c5, c6,  &
                     d1, d2, d3, d4, d5, d6   
    parameter(a2 = 0.2,  a3 = 0.3,  a4 = 0.6,  a5 = 1, a6 = .875)
    parameter(b21 = .2,  b31 = .075, b32 = .225, b41=.3, b42 = -.9, b43 = 1.2,          &
              b51 =-11./54, b52 = 2.5, b53 = -70./27, b54 = 35./27,  b61 = 1631./55296, &
              b62 = 175./512,  b63 = 575./13824,  b64 = 44275./110592,  b65 = 253./4096 )
    parameter(c1 = 37./378,  c3 = 250./621, c4 = 125./594,  c6 = 512./1771)
    parameter(d1 = c1 - 2825./27648,  d3 = c3 - 18575./48384,  d4 = c4 - 13525./55296,  &
              d5 = -277./14336,   d6 = c6 - .25 )    
    real(kind=16)  ::  h, htry, FDrv1(3), FDrv2(3), FDrv3(3), &
                                FDrv4(3), FDrv5(3), FDrv6(3), difF(3)
    real(kind=16)  ::  xr, xm, xy, xp, xe, xc ! Temporary record
    
    !-------------------!
    ! Read the EoS data !
    !-------------------!
    open(unit = 1 , file=trim(filename), status='old')
    !--------------------------------------------!
    read(1,*) ! Read the title 
    neos = 0
    do while(.not. eof(1))
        neos = neos + 1 
        !   Read      den       ene          pre         cs2 
        Read(1,*)  Rns(neos), Ens(neos),  Pns(neos),  Cns(neos)
    end do
    close(1)
    ! It is natural to set the exiting pressure as the lowest one in the input table.
    pre_out = minval(Pns)  
    if(pre_out<=0.d0) pre_out = 1d-12
    ! Establish the EOS 
    allocate(rho(neos), ene(neos), pre(neos), cs2(neos))
    rho = Rns(1:neos)
    ene = Ens(1:neos)
    pre = Pns(1:neos)
    cs2 = Cns(1:neos)
    
  
    open(unit = 0, file='RvsM.d')
    write(0,200)
      
     
    do nwrite= 1,  np
        ! Clear the values of every NS table 
        do n=1,nmax
            Rns(n)=0.d0 
            Mns(n)=0.d0
            Pns(n)=0.d0
            Ens(n)=0.d0
            Dns(n)=0.d0
            Cns(n)=0.d0
            Yns(n)=0.d0
        end do 
        !
        !==============================!
        ! Initial value:               !  
        !    Central NS Pressure (R=0) !              
        !==============================!
        xp = pres(nwrite)
        if(xp<=pre_out) exit
        Pns(1) = xp 
        Ens(1) = Fpolint(xp,pre,ene)
        Dns(1) = Fpolint(xp,pre,rho)
        Cns(1) = Fpolint(xp,pre,cs2)
        ! NS properties 
        Rns(1) = 0.d0
        Mns(1) = 0.d0
        Yns(1) = 2.0 ! For tidal deformability 
        ! 
        n = 1 
        htry = chmax
        do while(Pns(n)>=pre_out.and. Rns(n)<cRmax)            
1           if(htry>chmax) htry = chmax
            if(htry<chmin) htry = chmin
            h = htry
            xr = Rns(n) 
            xp = Pns(n) 
            xe = Ens(n)
            xc = Cns(n) 
            xm = Mns(n) 
            xy = Yns(n)
            ! First step
            FDrv1 = h*Deriv(xr, xm, xy, xp, xe, xc) 
            ! Second step
            xr = Rns(n) + a2*h 
            xp = Pns(n) + b21*FDrv1(1) 
            xm = Mns(n) + b21*FDrv1(2) 
            xy = Yns(n) + b21*FDrv1(3) 
            xe = Fpolint(xp,pre,ene)
            xc = Fpolint(xp,pre,cs2)
            FDrv2 = h*Deriv(xr, xm, xy, xp, xe, xc) 
            ! Third step
            xr = Rns(n) + a3*h 
            xp = Pns(n) + b31*FDrv1(1) + b32*FDrv2(1)
            xm = Mns(n) + b31*FDrv1(2) + b32*FDrv2(2)
            xy = Yns(n) + b31*FDrv1(3) + b32*FDrv2(3)
            xe = Fpolint(xp,pre,ene)
            xc = Fpolint(xp,pre,cs2)
            FDrv3 = h*Deriv(xr, xm, xy, xp, xe, xc) 
            ! Forth step 
            xr = Rns(n) + a4*h 
            xp = Pns(n) + b41*FDrv1(1) + b42*FDrv2(1) + b43*FDrv3(1)
            xm = Mns(n) + b41*FDrv1(2) + b42*FDrv2(2) + b43*FDrv3(2)
            xy = Yns(n) + b41*FDrv1(3) + b42*FDrv2(3) + b43*FDrv3(3)
            xe = Fpolint(xp,pre,ene)
            xc = Fpolint(xp,pre,cs2)
            FDrv4 = h*Deriv(xr, xm, xy, xp, xe, xc)  
            ! Fifth step
            xr = Rns(n) + a5*h 
            xp = Pns(n) + b51*FDrv1(1) + b52*FDrv2(1) + b53*FDrv3(1) + b54*FDrv4(1)
            xm = Mns(n) + b51*FDrv1(2) + b52*FDrv2(2) + b53*FDrv3(2) + b54*FDrv4(2)
            xy = Yns(n) + b51*FDrv1(3) + b52*FDrv2(3) + b53*FDrv3(3) + b54*FDrv4(3)
            xe = Fpolint(xp,pre,ene)
            xc = Fpolint(xp,pre,cs2)
            FDrv5 = h*Deriv(xr, xm, xy, xp, xe, xc) 
            ! Sixth step
            xr = Rns(n) + a6*h 
            xp = Pns(n) + b61*FDrv1(1) + b62*FDrv2(1) + b63*FDrv3(1) + b64*FDrv4(1) + b65*FDrv5(1)
            xm = Mns(n) + b61*FDrv1(2) + b62*FDrv2(2) + b63*FDrv3(2) + b64*FDrv4(2) + b65*FDrv5(2)
            xy = Yns(n) + b61*FDrv1(3) + b62*FDrv2(3) + b63*FDrv3(3) + b64*FDrv4(3) + b65*FDrv5(3)
            xe = Fpolint(xp,pre,ene)
            xc = Fpolint(xp,pre,cs2)
            FDrv6 = h*Deriv(xr, xm, xy, xp, xe, xc) 
            ! The difference 
            DifF = d1*FDrv1 + d3*FDrv3 + d4*FDrv4 + d5*FDrv5 + d6*FDrv6
            xe = maxval(abs(DifF))
            !
            if(h>chmin .and. xe>ceps) then ! Calculate again if given accuracy is not satisfied.
                htry = .9*(ceps/xe)**.25*h 
                go to 1 
            else 
                htry =  (ceps/xe)**.2*h
            end if 
            !
            ! Now step once 
            Rns(n+1) = Rns(n) + h
            Pns(n+1) = Pns(n) + c1*FDrv1(1) + c3*FDrv3(1) + c4*FDrv4(1) + c6*FDrv6(1)
            Mns(n+1) = Mns(n) + c1*FDrv2(2) + c3*FDrv3(2) + c4*FDrv4(2) + c6*FDrv6(2)
            Yns(n+1) = Yns(n) + c1*FDrv2(3) + c3*FDrv3(3) + c4*FDrv4(3) + c6*FDrv6(3)
            Ens(n+1) = Fpolint(Pns(n+1),pre,ene)
            Dns(n+1) = Fpolint(Pns(n+1),pre,rho)
            Cns(n+1) = Fpolint(Pns(n+1),pre,cs2)  
            n = n + 1 
        end do 
        !
        neos = n-1
        xc = cR1M*Mns(neos)/Rns(neos)
        xy = Yns(neos) 
        xe = 1.6*xc**5*(1-2*xc)**2*(2-xy+2*xc*(xy-1))/     &
             (6*xc*(2-xy+xc*(5*xy-8))+4*xc**3*(13-11*xy+   &
             xc*(3*xy-2) + 2*xc**2*(1+xy))+3*(1-2*xc)**2*  &
             (2-xy+2*xc*(xy-1))*log(1-2*xc))            ! k2
        xr = 2.*xe/(3*xc**5)  
        write(0,20) Rns(neos), Mns(neos), Dns(1), Pns(1), Ens(1), xc,  xy,  xe,  xr
        !
        if(Details) then 
            write(cpres,'(f10.2)')  pres(nwrite)
            open(unit=nwrite, file='NSdetails_P='//trim(adjustl(cpres))//'.d')
            write(nwrite,100)
            !
            do n = 1, neos
                xc = sign(Cns(n), sqrt(abs(Cns(n))))
                if(xc<=0) write(*,'(a,f10.4,a)') &
                    'Non-positive sound speed at R=',Rns(n),' km.'
                write(nwrite,10) Rns(n), Mns(n), Dns(n), Pns(n), Ens(n), &
                      xc, Yns(n)
            end do 
            close(nwrite)
        end if        
    end do 
    
10  format(7e15.4)
20  format(f10.4,7f12.4,f16.4)
    
100 format(5x,'Radius',10x,'Mass',13x, 'n',12x,'pre',12x,'ene',13x,'cs',14x,'y')    
200 format(4x,'Radius',7x,'Mass',9x,'nc',9x,'pre',9x,'ene',11x,'C',10x,'yR',10x,'k2',12x,'Lam')
    contains 
    
    
    function Deriv(r,m,y,p,e,c) result(Yout)
    !   The derivative part of R.H.S. of TOV equation &
    !       the tidal y function 
    !       r [km],  m [Mo],  y & c [1], p & e [MeV/fm^3] 
        real(kind=16)  :: r, m, y, p, e, c 
        real(kind=16)  :: Yout(3) 
        real(kind=16)  :: te, tp, td, tr, rr 
        
        rr = r*r 
        te = ceos2ns*rr*e   ! 4 pi r^2 e  [Mo/km]
        tp = ceos2ns*rr*p   ! 4 pi r^2 P  [Mo/km]
        tr = r/cR1M         ! Change R [km] to M [Mo]
        if(c<ceps) c = p/e  ! A stupid way to remove the singularity 
        if(r==0) then 
            Yout(1) = 0.d0 
            Yout(3) = 0.d0 
        else 
            td = 1./(1.-2*m/tr)
            Yout(1) = -(e+p)*(m+r*tp)/(r*(tr-2*m)) ! [MeV/fm^3  /  km]
            Yout(3) = -(  y**2 + y*td*(1.-cR1M*(te-tp))      & 
                      +(cR1M*(5*te+9*tp + (te+tp)/c) -6)*td &
                      -4*(td*(cR1M*tp+m/tr))**2  )/r        ! [km^(-1)]
        end if 
        Yout(2) = te ![Mo/km]
        return 
    end function 
 

    function Fpolint(x,xx,yy,np) result(y)
        ! One-dimensional interpolation by polynomials
        !    x --- input x to be evaluated 
        !    xx, yy --- the data arrays {(x,y)}
        !    np --- optional, the power for evalution
        !           default np=4   
        implicit none 
        real*16  :: x, y, xx(:), yy(:)
        integer, optional :: np 
        ! if present, this function will perform np order Polynomial evalution
        integer :: nv, i, jl, ju, jm
        real*16, allocatable :: xa(:), ya(:)
        
        nv = size(xx)! size of vector xx 
        ! Locate x by bisection
        jl = 0
        ju = nv+1 
        ! Initialize the lower and upper limit 
        do while(ju-jl>1) 
            jm = (ju+jl)/2
            ! simultaneously establishing
            if((xx(nv)>=xx(1)) == (x>=xx(jm))) then 
                jl = jm
            else 
                ju = jm 
            end if  
        end do 
       ! if(jl>nv.or.ju<1) then
       !    y = 0.d0; return
       ! end if 
            
        ! Now jl is the lower bound to locate x
        i = 4
        if(present(np)) i = np
        jm = 0.5*(i+1) 
        allocate(xa(i),ya(i))
        if(jl<=jm) then 
            xa = xx(1:i); ya = yy(1:i)
        else if (jl>=nv-jm) then
            xa = xx(1+nv-i:nv); ya = yy(1+nv-i:nv) 
        else 
            xa = xx(jl-jm+1:jl+jm); ya = yy(jl-jm+1:jl+jm)
        end if
        ! Begin to do interpolation by Neville's algorithm
        do jm = 1, i-1
            do ju = 1, i-jm
                ya(ju) = (ya(ju)*(x-xa(ju+jm)) &
                         +ya(ju+1)*(xa(ju)-x))/(xa(ju)-xa(ju+jm))
            end do
        end do
        y = ya(1)
        deallocate(xa,ya)
        return 
    end function 
    
    end program 