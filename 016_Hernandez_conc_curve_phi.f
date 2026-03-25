c Replicating Hernandez 1 with curvature..... real space curvature calculated , time increment..... ellipsoid eqn... protein coming in the middle

      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer  NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv, nuc
      double precision bx,by,bz, c_x, c_y,c_z, eps, r, h, a_x, a_y, a_z
      double precision nuu, r_ellip, r2, F, a
      parameter(Nx=64,Ny=64,Nz=64,NT=100000,snap_int=1000, c_max=3.d0)
      parameter(pi=acos(-1.d0), dt=0.001d0, grutol=0.2, lambda=0.1d0)
      parameter(r0=10,c_h=2.d0,Dc=0.1d0,As=1.0d0,Af=20.0d0, tol=0.1d0)
      parameter(Du=0.5d0, sig=4.d0, beta=0.0d0, lambda1=2.5d0)
      parameter(Lx=10.d0*pi, Ly=10.d0*pi, Lz=10.d0*pi, c_ini=1.0d0)
      parameter(alpha=5.0)

      double precision, allocatable :: u(:,:,:), c(:,:,:), fc(:,:,:)
      double precision, allocatable :: kx(:), ky(:), kz(:),ux(:,:,:)
      double precision, allocatable :: uy(:,:,:), uz(:,:,:)
      double precision, allocatable :: iku(:,:,:), fu(:,:,:)
      double precision, allocatable :: Curve(:,:,:), n_x(:,:,:)
      double precision, allocatable :: n_y(:,:,:), n_z(:,:,:)
      double precision, allocatable :: cc(:,:,:), ciku(:,:,:)
      
 
      complex*16, allocatable :: ck(:,:,:), fch(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ikuhy(:,:,:), ikuhz(:,:,:), uh(:,:,:)
      complex*16, allocatable :: fuh(:,:,:), Curveh(:,:,:)
      complex*16, allocatable :: cch(:,:,:), ikbcu(:,:,:), cikuh(:,:,:)

      integer*8 :: plan
      integer :: NxC, it, isnap, l
      complex*16, parameter :: iota = (0.d0, 1.d0)
      complex*16, parameter :: minus_iota = (0.d0, -1.d0)
      character(len=64) :: fname

      call execute_command_line('rm -f step_*.dat')
      call execute_command_line('rm -f step_*.vtu')
      call execute_command_line('rm -f cutc_step_*.dat')
      call execute_command_line('rm -f cutc_step_*.vtu')
      call execute_command_line('rm -f cute_step_*.dat')
      call execute_command_line('rm -f cute_step_*.vtu')          
c     Wavenumbers
      dx=Lx/Nx
      dy=Ly/Ny
      dz=Lz/Nz


      NxC=Nx/2+1

      norm=1.d0/(dble(Nx*Ny*Nz))
c      write(*,*)norm
      
      bx=2.0d0*pi/Lx
      by=2.0d0*pi/Ly
      bz=2.0d0*pi/Lz

      
      allocate(kx(NxC))
      allocate(ky(Ny))
      allocate(kz(Nz))


      do 12 i=1, NxC
      	kx(i)=bx*dble(i-1)
 12   continue
 
      do 21 i=1,Nz
      	if (i.le.(Nz/2+1)) then
      		kz(i)=bz*dble(i-1)
      	else
      		kz(i)=bz*dble(i-1-Nz)
      	endif
 21    continue
 
      do 11 i=1,Ny
      	if (i.le.(Ny/2+1)) then
      		ky(i)=by*dble(i-1)
      	else
      		ky(i)=by*dble(i-1-Ny)
      	endif
 11   continue


 
c     Creating initial condition-arrays
      allocate(u(Nx, Ny, Nz))
      allocate(c(Nx, Ny, Nz))
      allocate(fc(Nx, Ny, Nz))
      allocate(ux(Nx, Ny, Nz))
      allocate(uy(Nx, Ny, Nz))
      allocate(uz(Nx, Ny, Nz))
      allocate(iku(Nx, Ny, Nz))
      allocate(fu(Nx, Ny, Nz))
      allocate(Curve(Nx, Ny, Nz))
      allocate(n_x(Nx, Ny, Nz))
      allocate(n_y(Nx, Ny, Nz))
      allocate(n_z(Nx, Ny, Nz))
      allocate(cc(Nx, Ny, Nz))
      allocate(ciku(Nx, Ny, Nz))
      
              
      allocate(ck(NxC, Ny, Nz))
      allocate(fch(NxC, Ny, Nz))
      allocate(uh(NxC, Ny, Nz))
      allocate(fuh(NxC, Ny, Nz))
      allocate(Curveh(NxC, Ny, Nz))
      allocate(ikbcu(NxC, Ny, Nz))
      allocate(cch(NxC, Ny, Nz))
      allocate(cikuh(NxC, Ny, Nz))

      
C***************Initiation of fields*************************
c Initial phi field
      c_x=(Nx*1.d0+1)/2
      c_y=(Ny*1.d0+1)/2
      c_z=(Nz*1.d0+1)/2
      
      
c For elipsoidal only
      a_x = 25.d0    ! semi-axis length along x
      a_y = 15.d0    ! semi-axis length along y
      a_z = 15.d0    ! semi-axis length along z

c For Nephroidal & Concave Plain

      a=25    
       do 13 i=1,Nx
        do 14 j=1,Ny
         do 15 k=1,Nz

c***************************************************************
c    Circular Initial Field
c         r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
c          u(i,j,k)=-tanh(alpha*(r-r0))
c***************************************************************

c***************************************************************
c    Cubic Initial Field 
! cube half width (distance from center to face)
c           a = r0         ! use r0 as cube half-length or rename as needed

c           dx = abs(dble(i) - c_x) 
c           dy = abs(dble(j) - c_y) 
c           dz = abs(dble(k) - c_z) 
           
c           if ((dx.le.a).and.(dy.le.a).and.(dz.le.a)) then
c             u(i,j,k)=1.d0
c           else
c             u(i,j,k)=-1.d0
c          endif
c******************************************************************

c*************************************************************
c    *********** Ellipsoidal Initial Field ***********

            xr = (dble(i) - c_x) / a_x
            yr = (dble(j) - c_y) / a_y
            zr = (dble(k) - c_z) / a_z
            r_ellip = sqrt(xr**2 + yr**2 + zr**2)
            u(i,j,k) = -tanh(alpha * (r_ellip - 1.d0))

c            if (r_ellip .le. 1.d0) then
c               u(i,j,k) =  1.d0     ! inside ellipsoid
c            else
c               u(i,j,k) = -1.d0     ! outside ellipsoid
c            endif

c**************Nephroidal**************

c       xr = (dble(i) - c_x)/a
c       yr = (dble(j) - c_y)/a
c       zr = (dble(k) - c_z)/a

c       r2 = xr**2 + yr**2 + zr**2

c       F = (r2 -0.8**2)**3 - 0.8 * yr**2

c       u(i,j,k) = -tanh(alpha * F)

c            if (F.le.0.d0) then
c               u(i,j,k) =  1.d0     ! inside nephroid
c            else
c              u(i,j,k) = -1.d0     ! outside nephroid
c            endif           


c*************Concave Plain**************
c        xr = (dble(i) - c_x)/a
c        yr = (dble(j) - c_y)/a
c        zr = (dble(k) - c_z)/a

c        F = zr - (xr**2 + yr**2)

c       u(i,j,k) = -tanh(alpha * F)

c            if (F.le.0.d0) then
c               u(i,j,k) =  1.d0     ! inside concave
c            else
c              u(i,j,k) = -1.d0     ! outside concave
c           endif           
c**********************************************
 15     continue
 14    continue
 13   continue
 
 
 
c*************************Initial Run************************************************

       do 65 it=1, 1
       
       do 59 i=1,Nx
        do 60 j=1,Ny
         do 61 k=1, Nz
          
           fu(i,j,k)=Du*(4*As*u(i,j,k)*(u(i,j,k)**2-1)*c(i,j,k)**2
     1       *(c(i,j,k)-c_max)**2+Af*2*u(i,j,k)*c(i,j,k)**2)
          
 61      continue    
 60     continue    
 59    continue
 
       call fft_forward(Nx, Ny, Nz, u, uh)
       call fft_forward(Nx, Ny, Nz, fu, fuh)

       do 62 i=1, NxC
         do 63 j=1, Ny
          do 64 k=1, Nz      	    
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)


            nuu=Du*2*sig*k2
            uh(i,j,k)=exp(-nuu*k2*dt)*(uh(i,j,k)-k2*fuh(i,j,k)*dt)
 64      continue
 63     continue
 62    continue

 65    continue

       call fft_backward(Nx, Ny, Nz, uh, u)
       u=u*norm
c*******************************************************************






      call compute_grad_u(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
c3      call bending_curvature(Nx,Ny,Nz,uh,kx,ky,kz,iku,Curve,tol,
c3     1 norm)
c5      call compute_curvature(Nx, Ny, Nz, NxC,uh,kx,ky,kz,tol,norm,
c5     1     ux, uy, uz, iku, Curve)                       
c4      call make_paraboloid_curvature(Nx,Ny,Nz,u,iku,
c4     1 c_x,c_y,c_z,a_x,a_y,a_z,Curve,tol)
        call compute_curvature_real(Nx,Ny,Nz,NxC,uh,dx,dy,dz,
     1 tol,norm,ux,uy,uz,iku,Curve)
c    Initial Concentration Field
      do 23 i=1,Nx
       do 24 j=1,Ny
        do 25 k=1,Nz
          h=(dble(k)-c_z)
          r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
c          if ((h.ge.0).and.(abs(iku(i,j,k)).ge.grutol)) then
c          if ((r.ge.(r0-0.5)).and.(r.le.(r0+0.5))) then
           if (abs(u(i,j,k)).le.0.9) then
              c(i,j,k)=c_ini
           else
              c(i,j,k)=0
           endif
c2           c(i,j,k)=c_ini
 25     continue
 24    continue
 23    continue
       

c**********Saving the initial file*****************************
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm

      open(8, file='step_0.dat')

      do 41 i=1,Nx
        do 52 j=1, Ny
         do 70 k=1,Nz
          write(8,*)i, j, k, u(i,j,k), c(i,j,k), Curve(i,j,k),   
     1      sqrt(iku(i,j,k))
 70     continue
 52    continue
 41   continue
      close(8)


      call fft_forward(Nx, Ny, Nz, c, ck)
      call fft_forward(Nx, Ny, Nz, u, uh)

      


c**************************************************************

c*****************Time Marching********************************

      do 3 it=1, NT
    
      call fft_backward(Nx, Ny, Nz, ck, c)
      call fft_backward(Nx, Ny, Nz, uh, u)
      
      c=c*norm
      u=u*norm
             
       do 19 i=1,Nx
        do 20 j=1,Ny
         do 22 k=1, Nz

          fc(i,j,k)=Dc*(2*As*(((u(i,j,k)**2)-1)**2)*((c(i,j,k)-c_max)
     1    *c(i,j,k)*(2*c(i,j,k)-c_max))+Af*2*(u(i,j,k)**2)*c(i,j,k))
             
c1           fu(i,j,k)=Du*(4*As*u(i,j,k)*(u(i,j,k)**2-1)*c(i,j,k)**2
c1     1       *(c(i,j,k)-c_max)**2+Af*2*u(i,j,k)*c(i,j,k)**2)
     
           cc(i,j,k)=Curve(i,j,k)*c(i,j,k) 
             
c           ciku(i,j,k)=Curve(i,j,k)*sqrt(iku(i,j,k))       
     
 22      continue    
 20     continue    
 19    continue
 
       call fft_forward(Nx, Ny, Nz, fc, fch)
       call fft_forward(Nx, Ny, Nz, c, ck)
       
       call fft_forward(Nx, Ny, Nz, u, uh)
c1       call fft_forward(Nx, Ny, Nz, fu, fuh)

       call fft_forward(Nx, Ny, Nz, cc, cch)
c      call fft_forward(Nx, Ny, Nz, ciku, cikuh)
 
 
       do 26 i=1, NxC
         do 27 j=1, Ny
          do 28 k=1, Nz      	    
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)


      	    nuc=Dc*lambda*k2*2*As
      	    
            ck(i,j,k)=exp(-k2*nuc*dt)*(ck(i,j,k)-k2*fch(i,j,k)*dt
     1        +lambda1*(-k2)*cch(i,j,k)*dt)
            

c1            nuu=Du*2*sig*k2
c1            uh(i,j,k)=(exp(-nuu*k2*dt))*(uh(i,j,k)-k2*fuh(i,j,k)*dt)
c     1      -beta*cikuh(i,j,k)*dt)
 28      continue
 27     continue
 26    continue
 
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
 
      call compute_grad_u(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
     
c3      call bending_curvature(Nx,Ny,Nz,uh,kx,ky,kz,iku,Curve,tol,
c3     1 norm)
     
c      call compute_curvature(Nx, Ny, Nz, NxC,uh,kx, ky, kz, tol,norm, 
c     1     ux, uy, uz, iku, Curve)

c4       call make_paraboloid_curvature(Nx,Ny,Nz,u,iku,
c4     1 c_x,c_y,c_z,a_x,a_y,a_z,Curve,tol)
         
c***********************************************
       
     

c********Saving files*******

      if (mod(it,snap_int).eq.0) then      
      
       
       call fft_backward(Nx, Ny, Nz, ck, c)
       call fft_backward(Nx, Ny, Nz, uh, u)
       c=c*norm
       u=u*norm
       
      
       write(fname, '("step_", I0, ".dat")') it
       open(8, file=fname, status='replace')      
       do 1 i=1,Nx
        do 4 j=1,Ny
         do 5 k=1, Nz
          write(8,*)i, j, k, u(i,j,k), c(i,j,k), Curve(i,j,k), 
     1      sqrt(iku(i,j,k))
 5       continue
 4      continue
 1     continue
       close(8)
      
       call fft_forward(Nx, Ny, Nz, c, ck)
       call fft_forward(Nx, Ny, Nz, u, uh)
      endif
      

 3    continue


c *******************************************************
      call dfftw_cleanup()
      stop
      end
      
c*******************************************************************      
      subroutine fft_forward(Nx, Ny, Nz, in, out)
      implicit none
      include 'fftw3.f'
      integer, intent(in) :: Nx, Ny, Nz
      real*8, intent(inout) :: in(Nx,Ny,Nz)
      complex*16, intent(inout) :: out(Nx/2+1,Ny,Nz)
      integer*8 :: plan

      call dfftw_plan_dft_r2c_3d(plan,Nx,Ny,Nz,in,out,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,in,out)
      call dfftw_destroy_plan(plan)
      end subroutine fft_forward
c****************************************************************      
      subroutine fft_backward(Nx, Ny, Nz, in, out)
      implicit none
      include 'fftw3.f'
      integer, intent(in) :: Nx, Ny, Nz
      complex*16, intent(inout) :: in(Nx/2+1,Ny,Nz)
      real*8, intent(inout) :: out(Nx,Ny,Nz)
      integer*8 :: plan

      call dfftw_plan_dft_c2r_3d(plan,Nx,Ny,Nz,in,out,FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(plan, in, out)
      call dfftw_destroy_plan(plan)
      end subroutine fft_backward
      
c*****************************************************************      
      subroutine compute_grad_u(Nx,Ny,Nz,NxC,u,uh,kx,ky,kz
     1       ,ux,uy,uz,iku,norm)
      implicit none
      include 'fftw3.f'

      integer, intent(in) :: Nx, Ny, Nz, NxC
      double precision :: u(Nx,Ny,Nz)
      complex*16, intent(inout) :: uh(NxC,Ny,Nz)
      double precision :: kx(NxC), ky(Ny), kz(Nz)
      double precision :: ux(Nx,Ny,Nz), uy(Nx,Ny,Nz)
      double precision :: uz(Nx,Ny,Nz)
      double precision :: iku(Nx,Ny,Nz)     
      real*8, intent(in) :: norm             

      complex*16 :: uhx(NxC,Ny,Nz), uhy(NxC,Ny,Nz)
      complex*16 :: uhz(NxC,Ny,Nz)
      integer :: i, j, k
      complex*16, parameter :: iota = (0.d0, 1.d0)


      call fft_forward(Nx, Ny, Nz, u, uh)


      do 51 i = 1, NxC
        do 52 j = 1, Ny
            do 53 k = 1, Nz
                uhx(i,j,k) = iota * kx(i) * uh(i,j,k)
                uhy(i,j,k) = iota * ky(j) * uh(i,j,k)
                uhz(i,j,k) = iota * kz(k) * uh(i,j,k)
 53         continue
 52      continue
 51    continue 

      call fft_backward(Nx, Ny, Nz, uhx, ux)
      call fft_backward(Nx, Ny, Nz, uhy, uy)
      call fft_backward(Nx, Ny, Nz, uhz, uz)

      ux = ux * norm
      uy = uy * norm
      uz = uz * norm


      do 54 i = 1, Nx
        do 55 j = 1, Ny
            do 56 k = 1, Nz
                iku(i,j,k) = ux(i,j,k)**2 + uy(i,j,k)**2 
     1           + uz(i,j,k)**2
 56         continue
 55      continue
 54    continue 

      end subroutine compute_grad_u



c*****************************Fourier Space Curvature********************************************
      subroutine compute_curvature(Nx, Ny, Nz, NxC, uh, kx,ky,kz,tol,
     1 norm, ux, uy, uz, iku, Curve)
      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer :: Nx, Ny, Nz
      double precision  :: kx(NxC), ky(Ny), kz(Nz)
      double precision  :: norm, mag
      double precision :: ux(Nx,Ny,Nz), uy(Nx,Ny,Nz)
      double precision :: uz(Nx,Ny,Nz)
      double precision :: iku(Nx,Ny,Nz)
      double precision :: Curve(Nx,Ny,Nz)
      
      double precision :: n_x(Nx,Ny,Nz), n_y(Nx,Ny,Nz)
      double precision :: n_z(Nx,Ny,Nz)
      double precision :: ikaux(Nx,Ny,Nz), ikauy(Nx,Ny,Nz)
      double precision :: ikauz(Nx,Ny,Nz), u(Nx, Ny, Nz)

      complex*16 :: n_xh(NxC,Ny,Nz), n_yh(NxC,Ny,Nz)
      complex*16 :: n_zh(NxC,Ny,Nz), uh(NxC, Ny, Nz)
      complex*16 :: ikauxh, ikauyh, ikauzh
      complex*16 :: ikauh(NxC,Ny,Nz)
      
      

      integer :: i, j, k
      complex*16, parameter :: iota = (0.d0, 1.d0)
c************************************************************


      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
      
      
      do 38  i = 1, Nx
        do 39 j = 1, Ny
            do 40 k = 1, Nz
            
		mag=sqrt(iku(i,j,k))
                if (abs(u(i,j,k)).lt.1) then
c                 if (mag.ge.tol) then
                    n_x(i,j,k) = ux(i,j,k) / mag
                    n_y(i,j,k) = uy(i,j,k) / mag
                    n_z(i,j,k) = uz(i,j,k) / mag
                else

                    n_x(i,j,k) = 0.d0
                    n_y(i,j,k) = 0.d0
                    n_z(i,j,k) = 0.d0
               endif


 40         continue
 39     continue
 38    continue

c      call fft_backward(Nx, Ny, Nz, uh, u)
c      u=u*norm

c      open(300,file="normals.dat")
c       do i=1,Nx
c        do j=1,Ny
c          do k=1,Nz
c           write(300,*) i,j,k,u(i,j,k), n_x(i,j,k), n_y(i,j,k), 
c     1      n_z(i,j,k)
c          end do
c        end do
c      end do
c      close(300)

c      call fft_forward(Nx, Ny, Nz, u, uh)      

      call fft_forward(Nx, Ny, Nz, n_x, n_xh)
      call fft_forward(Nx, Ny, Nz, n_y, n_yh)
      call fft_forward(Nx, Ny, Nz, n_z, n_zh)
      


      do 42 i = 1, NxC
        do 43 j = 1, Ny
            do 44 k = 1, Nz
                ikauxh = iota * kx(i) * n_xh(i,j,k)
                ikauyh = iota * ky(j) * n_yh(i,j,k)
                ikauzh = iota * kz(k) * n_zh(i,j,k)
                ikauh(i,j,k)=ikauxh+ikauyh+ikauzh                
 44         continue
 43     continue
 42   continue
 


      call fft_backward(Nx, Ny, Nz, ikauh, Curve)
      
      do 138  i = 1, Nx
        do 139 j = 1, Ny
            do 140 k = 1, Nz
            
c		if (sqrt(iku(i,j,k)).le.tol) Curve(i,j,k)=0.d0
              if (abs(u(i,j,k)).gt.0.2) Curve(i,j,k)=0.d0
 140     continue                
 139     continue
 138    continue

      Curve=-1.d0*Curve*norm

      call fft_forward(Nx, Ny, Nz, u, uh)
      end subroutine compute_curvature
c**********************************************************************

c************Paraboloid Curvature************************
      subroutine make_paraboloid_curvature(Nx,Ny,Nz,uh,iku,
     1 c_x,c_y,c_z,a_x,a_y,a_z,Curve,tol)

      implicit double precision(a-h,o-z)

      integer Nx,Ny,Nz, i,j,k
      double precision u(Nx,Ny,Nz), iku(Nx,Ny,Nz), Curve(Nx,Ny,Nz)
      double precision c_x,c_y,c_z,tol,x,y,z,kappa0,gamma

      kappa0 = 15.0d0      ! peak curvature
      gamma  =1.d0     ! paraboloid strength


      do 71 i=1,Nx
       do 72 j=1,Ny
        do 73 k=1,Nz

          x = (dble(i) - c_x)/a_x
          y = (dble(j) - c_y)/a_y
          z = (dble(k) - c_z)/a_z

          if (sqrt(iku(i,j,k)) .ge. tol) then
            Curve(i,j,k) = gamma + kappa0*(x**2)
          else
            Curve(i,j,k) = 0.d0

          endif

 73     continue
 72    continue
 71   continue

      end subroutine make_paraboloid_curvature
      
      
      
      
      
      
      
      
      
c**********Make bending curvature*********
      subroutine bending_curvature(Nx,Ny,Nz,uh,kx,ky,kz,iku,Curve,tol,
     1 norm)      
     
      implicit double precision(a-h,o-z)

      integer Nx,Ny,Nz, i,j,k
      double precision eps, kx(Nx/2+1), ky(Ny), kz(Nz), k2u(Nx,Ny,Nz)
      double precision u(Nx,Ny,Nz), iku(Nx,Ny,Nz), Curve(Nx,Ny,Nz)
      double precision norm
      complex*16 k2uh(Nx/2+1, Ny, Nz), uh(Nx/2+1, Ny, Nz)
      
      eps=0.1
      
      do 77 i=1, Nx/2+1
       do 78 j=1,Ny
        do 79 k=1,Nz
         k2=(kx(i)**2+ky(j)**2+kz(k)**2)
         k2uh(i,j,k)=k2*uh(i,j,k)
 79     continue
 78    continue
 77   continue
      
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
      call fft_backward(Nx, Ny, Nz, k2uh, k2u)
      k2u=k2u*norm
      
      do 74 i=1,Nx
       do 75 j=1,Ny
        do 76 k=1,Nz
           Curve(i,j,k)=(-u(i,j,k)+u(i,j,k)**3-eps**2*k2u(i,j,k))
 76     continue
 75    continue
 74   continue
      
      call fft_forward(Nx, Ny, Nz, u, uh)

      
      end subroutine bending_curvature
      
c*****************************************************************










c************Real Space Curvature**********************

      subroutine compute_curvature_real(Nx, Ny, Nz, NxC, uh, dx,dy,dz,
     1 tol,norm, ux, uy, uz, iku, Curve)
      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer :: Nx, Ny, Nz,NxC
      double precision :: dx, dy, dz, tol,kappa
      double precision  :: norm, mag
      double precision :: ux(Nx,Ny,Nz), uy(Nx,Ny,Nz)
      double precision :: uz(Nx,Ny,Nz)
      double precision :: iku(Nx,Ny,Nz)
      double precision :: Curve(Nx,Ny,Nz)
      
      double precision :: n_x(Nx,Ny,Nz), n_y(Nx,Ny,Nz)
      double precision :: n_z(Nx,Ny,Nz)
      double precision :: ikaux(Nx,Ny,Nz), ikauy(Nx,Ny,Nz)
      double precision :: ikauz(Nx,Ny,Nz), u(Nx, Ny, Nz)

      complex*16 :: n_xh(NxC,Ny,Nz), n_yh(NxC,Ny,Nz)
      complex*16 :: n_zh(NxC,Ny,Nz), uh(NxC, Ny, Nz)
      complex*16 :: ikauxh, ikauyh, ikauzh
      complex*16 :: ikauh(NxC,Ny,Nz)
      
      

      integer :: i, j, k
      complex*16, parameter :: iota = (0.d0, 1.d0)
      
      logical valid_p, valid_m
c************************************************************


      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
      
      
      do 381  i = 1, Nx
        do 391 j = 1, Ny
            do 401 k = 1, Nz
            
		mag=sqrt(iku(i,j,k))
c                if (abs(u(i,j,k)).lt.tol) then
                 if (mag.ge.tol) then

                    n_x(i,j,k) = ux(i,j,k) / mag
                    n_y(i,j,k) = uy(i,j,k) / mag
                    n_z(i,j,k) = uz(i,j,k) / mag
                else

                    n_x(i,j,k) = 0.d0
                    n_y(i,j,k) = 0.d0
                    n_z(i,j,k) = 0.d0
               endif

 401         continue
 391     continue
 381    continue



c      open(300,file="normals.dat")
c       do i=1,Nx
c        do j=1,Ny
c          do k=1,Nz
c           write(300,*) i,j,k,u(i,j,k), n_x(i,j,k), n_y(i,j,k), 
c     1      n_z(i,j,k)
c          end do
c        end do
c      end do
c      close(300)

c Interface only curvature

c Central Difference Formula

c       do 80 i = 2, Nx-1
c        do 81 j = 2, Ny-1
c         do 82 k = 2, Nz-1
c          mag=sqrt(iku(i,j,k)) 
c           if (mag.le.0.01) then
c            Curve(i,j,k)=0
c           else
c            dnx_dx = (n_x(i+1,j,k) - n_x(i-1,j,k)) / (2.d0*dx)
c            dny_dy = (n_y(i,j+1,k) - n_y(i,j-1,k)) / (2.d0*dy)
c            dnz_dz = (n_z(i,j,k+1) - n_z(i,j,k-1)) / (2.d0*dz)
c            kappa = -(dnx_dx + dny_dy + dnz_dz)
c            Curve(i,j,k)=kappa/2.d0*(tanh(20*(u(i,j,k)+0.7))
c     1                   +tanh(20*(-u(i,j,k)+0.7)))
c           endif
c 82      continue
c 81     continue
c 80    continue

c-------------------3 point formula with boundary correction-----------------------------
       do 83 i = 3, Nx-2
       	do 84 j = 3, Ny-2
       	 do 85 k = 3, Nz-2
       	   mag=sqrt(iku(i,j,k))
       	   if (mag.le.tol) then
       	    Curve(i,j,k)=0
       	   else
c================X-direction========================
            valid_p=((sqrt(iku(i+1,j,k)).ge.tol))
            valid_m=((sqrt(iku(i-1,j,k)).ge.tol))
                 
            if (valid_p.and.valid_m) then
             dnx_dx = (n_x(i+1,j,k) - n_x(i-1,j,k)) / (2.d0*dx)
            else if (valid_p.and.sqrt(iku(i+2,j,k)).ge.tol) then
             dnx_dx = (-3*n_x(i,j,k) + 4*n_x(i+1,j,k) - n_x(i+2,j,k))
     1        / (2.d0*dx)
            else if (valid_m.and.sqrt(iku(i-2,j,k)).ge.tol) then
             dnx_dx = (3*n_x(i,j,k) - 4*n_x(i-1,j,k) + n_x(i-2,j,k))
     1       / (2.d0*dx)
            else
             dnx_dx=0.d0
            endif
c===============Y-direction=======================
            valid_p=((sqrt(iku(i,j+1,k)).ge.tol))
            valid_m=((sqrt(iku(i,j-1,k)).ge.tol))
                 
            if (valid_p.and.valid_m) then
             dny_dy = (n_y(i,j+1,k) - n_y(i,j-1,k)) / (2.d0*dy)
            else if (valid_p.and.sqrt(iku(i,j+2,k)).ge.tol) then
             dny_dy = (-3*n_y(i,j,k) + 4*n_y(i,j+1,k) - n_y(i,j+2,k))
     1        / (2.d0*dy)
            else if (valid_m.and.sqrt(iku(i,j-2,k)).ge.tol) then
             dny_dy = (3*n_y(i,j,k) - 4*n_y(i,j-1,k) + n_y(i,j-2,k))
     1       / (2.d0*dy)
            else
             dny_dy=0.d0
            endif                  
c===============Z-direction=======================
            valid_p=((sqrt(iku(i,j,k+1)).ge.tol))
            valid_m=((sqrt(iku(i,j,k-1)).ge.tol))
                 
            if (valid_p.and.valid_m) then
             dnz_dz = (n_z(i,j,k+1) - n_z(i,j,k-1)) / (2.d0*dz)
            else if (valid_p.and.sqrt(iku(i,j,k+2)).ge.tol) then
             dnz_dz = (-3*n_z(i,j,k) + 4*n_z(i,j,k+1) - n_z(i,j,k+2))
     1        / (2.d0*dz)
            else if (valid_m.and.sqrt(iku(i,j,k-2)).ge.tol) then
             dnz_dz = (3*n_z(i,j,k) - 4*n_z(i,j,k-1) + n_z(i,j,k-2))
     1       / (2.d0*dz)
            else
            dnz_dz=0.d0
            endif
            kappa=-(dnx_dx+dny_dy+dnz_dz)
            Curve(i,j,k)=kappa
            
           endif
           
 85      continue
 84     continue
 83    continue
c-------------------------------------------------------------------------------------
       call fft_forward(Nx, Ny, Nz, u, uh)
       end subroutine compute_curvature_real

c**************************************************************
