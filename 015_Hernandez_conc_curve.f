c Replicating Hernandez 1 with curvature..... two curvature calculated , time increment... no phi dynamics... nephroid eqn

      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer  NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv, nuc
      double precision bx,by,bz, c_x, c_y,c_z, eps, r, h, a_x, a_y, a_z
      double precision nuu, r_ellip, r2, F, a
      parameter(Nx=32,Ny=32,Nz=32,NT=100000,snap_int=1000, c_max=3.d0)
      parameter(pi=acos(-1.d0), dt=0.001d0, grutol=0.2, lambda=0.01d0)
      parameter(r0=6,c_h=2.d0,Dc=0.1d0,As=1.0d0,Af=20.0d0, tol=0.1d0)
      parameter(Du=0.5d0, sig=4.d0, beta=0.1d0, lambda1=0.05d0)
      parameter(Lx=10.d0*pi, Ly=10.d0*pi, Lz=10.d0*pi, c_ini=0.3d0)
      parameter(alpha=20.0)

      double precision, allocatable :: u(:,:,:), c(:,:,:), fc(:,:,:)
      double precision, allocatable :: kx(:), ky(:), kz(:),ux(:,:,:)
      double precision, allocatable :: uy(:,:,:), uz(:,:,:)
      double precision, allocatable :: iku(:,:,:), fu(:,:,:)
      double precision, allocatable :: Curve(:,:,:), abs_ux(:,:,:)
      double precision, allocatable :: abs_uy(:,:,:), abs_uz(:,:,:)
      double precision, allocatable :: cc(:,:,:)
      
 
      complex*16, allocatable :: ck(:,:,:), fch(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ikuhy(:,:,:), ikuhz(:,:,:), uh(:,:,:)
      complex*16, allocatable :: fuh(:,:,:), Curveh(:,:,:)
      complex*16, allocatable :: cch(:,:,:), ikbcu(:,:,:)

      integer*8 :: plan
      integer :: NxC, it, isnap, l
      complex*16, parameter :: iota = (0.d0, 1.d0)
      complex*16, parameter :: minus_iota = (0.d0, -1.d0)
      character(len=64) :: fname

      call execute_command_line('rm -f step_*.dat')
      call execute_command_line('rm -f step_*.vtu')
          
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
      allocate(abs_ux(Nx, Ny, Nz))
      allocate(abs_uy(Nx, Ny, Nz))
      allocate(abs_uz(Nx, Ny, Nz))
      allocate(cc(Nx, Ny, Nz))
      
              
      allocate(ck(NxC, Ny, Nz))
      allocate(fch(NxC, Ny, Nz))
      allocate(uh(NxC, Ny, Nz))
      allocate(fuh(NxC, Ny, Nz))
      allocate(Curveh(NxC, Ny, Nz))
      allocate(ikbcu(NxC, Ny, Nz))
      allocate(cch(NxC, Ny, Nz))

      
C***************Initiation of fields*************************
c Initial phi field
      c_x=(Nx*1.d0+1)/2
      c_y=(Ny*1.d0+1)/2
      c_z=(Nz*1.d0+1)/2
      
      
c For elipsoidal only
      a_x = 12.d0    ! semi-axis length along x
      a_y = 8.d0    ! semi-axis length along y
      a_z = 8.d0    ! semi-axis length along z

c For Nephroidal

      a=10      
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

c          dx = abs(dble(i) - c_x) 
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

c            xr = (dble(i) - c_x) / a_x
c            yr = (dble(j) - c_y) / a_y
c            zr = (dble(k) - c_z) / a_z
c            r_ellip = sqrt(xr**2 + yr**2 + zr**2)
c            u(i,j,k) = -tanh(alpha * (r_ellip - 1.d0))

c            if (r_ellip .le. 1.d0) then
c               u(i,j,k) =  1.d0     ! inside ellipsoid
c            else
c               u(i,j,k) = -1.d0     ! outside ellipsoid
c            endif

c**************Nephroidal**************

       xr = (dble(i) - c_x)/a
       yr = (dble(j) - c_y)/a
       zr = (dble(k) - c_z)/a

       r2 = xr**2 + yr**2 + zr**2

       F = (r2 -0.8**2)**3 - 0.8 * yr**2

       u(i,j,k) = -tanh(alpha * F)

c            if (F.le.0.d0) then
c               u(i,j,k) =  1.d0     ! inside ellipsoid
c            else
c              u(i,j,k) = -1.d0     ! outside ellipsoid
c            endif           
 15     continue
 14    continue
 13   continue
 
 
 
c*************************Initial Run************************************************

       do 65 it=1, 10000
       
       do 59 i=1,Nx
        do 60 j=1,Ny
         do 61 k=1, Nz
          
           fu(i,j,k)=-Du*(4*As*u(i,j,k)*(u(i,j,k)**2-1)*c(i,j,k)**2
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


            nuu=Du*2*sig
            uh(i,j,k)=exp(-nuu*k2*dt)*(uh(i,j,k)+k2*fuh(i,j,k)*dt)
 64      continue
 63     continue
 62    continue

 65    continue

       call fft_backward(Nx, Ny, Nz, uh, u)
       u=u*norm
c*******************************************************************






      call compute_grad_u(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
      call compute_curvature(Nx, Ny, Nz, NxC,uh,kx,ky,kz,tol,norm,
     1     ux, uy, uz, iku, Curve)                       
c      call compute_curvature_real(Nx,Ny,Nz,uh,dx,dy,dz,
c     1     ux,uy,uz,iku,Curve)
c    Initial Concentration Field
      do 23 i=1,Nx
       do 24 j=1,Ny
        do 25 k=1,Nz
          h=(dble(k)-c_z)
          r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
c          if ((h.ge.0).and.(abs(iku(i,j,k)).ge.grutol)) then
c          if ((h.ge.0).and.(r.ge.(r0-0.5)).and.(r.le.(r0+0.5))) then
c              c(i,j,k)=10
c          else
c              c(i,j,k)=0
c          endif
           c(i,j,k)=c_ini
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
             
c           fu(i,j,k)=-Du*(4*As*u(i,j,k)*(u(i,j,k)**2-1)*c(i,j,k)**2
c     1       *(c(i,j,k)-c_max)**2+Af*2*u(i,j,k)*c(i,j,k)**2)
     
           cc(i,j,k)=Curve(i,j,k)*c(i,j,k)          
     
 22      continue    
 20     continue    
 19    continue
 
       call fft_forward(Nx, Ny, Nz, fc, fch)
       call fft_forward(Nx, Ny, Nz, c, ck)
       
       call fft_forward(Nx, Ny, Nz, u, uh)
c       call fft_forward(Nx, Ny, Nz, fu, fuh)

       call fft_forward(Nx, Ny, Nz, cc, cch)
 
 
       do 26 i=1, NxC
         do 27 j=1, Ny
          do 28 k=1, Nz      	    
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)


      	    nuc=Dc*lambda*k2*2*As
      	    
            ck(i,j,k)=exp(-k2*nuc*dt)*(ck(i,j,k)-k2*fch(i,j,k)*dt
     1        +lambda1*(-k2)*cch(i,j,k)*dt)
            
c            ikbcu(i,j,k)=iota*(kx(i)+ky(j)+kz(k))*beta
c     1        *Curveh(i,j,k)*u(i,j,k)
c            nuu=Du*2*sig
c            uh(i,j,k)=exp(-nuu*k2*dt)*(uh(i,j,k)+k2*fuh(i,j,k)*dt)
 28      continue
 27     continue
 26    continue
 
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
 
      call compute_grad_u(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
     
      call compute_curvature(Nx, Ny, Nz, NxC,uh,kx, ky, kz, tol,norm, 
     1     ux, uy, uz, iku, Curve)
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

    ! differentation in Fourier space (ik * uhat)
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



c*****************************Fourier Space********************************************
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
      
      double precision :: abs_ux(Nx,Ny,Nz), abs_uy(Nx,Ny,Nz)
      double precision :: abs_uz(Nx,Ny,Nz)
      double precision :: ikaux(Nx,Ny,Nz), ikauy(Nx,Ny,Nz)
      double precision :: ikauz(Nx,Ny,Nz), u(Nx, Ny, Nz)

      complex*16 :: abs_uxh(NxC,Ny,Nz), abs_uyh(NxC,Ny,Nz)
      complex*16 :: abs_uzh(NxC,Ny,Nz), uh(NxC, Ny, Nz)
      complex*16 :: ikauxh, ikauyh, ikauzh
      complex*16 :: ikauh(NxC,Ny,Nz)
      
      

      integer :: i, j, k
      complex*16, parameter :: iota = (0.d0, 1.d0)
c************************************************************

      do 38  i = 1, Nx
        do 39 j = 1, Ny
            do 40 k = 1, Nz
            
		mag=sqrt(iku(i,j,k))
                if (mag.gt.tol) then
                    abs_ux(i,j,k) = ux(i,j,k) / mag
                    abs_uy(i,j,k) = uy(i,j,k) / mag
                    abs_uz(i,j,k) = uz(i,j,k) / mag
                else

                    abs_ux(i,j,k) = 0.d0
                    abs_uy(i,j,k) = 0.d0
                    abs_uz(i,j,k) = 0.d0
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
c           write(300,*) i,j,k,u(i,j,k), abs_ux(i,j,k), abs_uy(i,j,k), 
c     1      abs_uz(i,j,k)
c          end do
c        end do
c      end do
c      close(300)

c      call fft_forward(Nx, Ny, Nz, u, uh)      

      call fft_forward(Nx, Ny, Nz, abs_ux, abs_uxh)
      call fft_forward(Nx, Ny, Nz, abs_uy, abs_uyh)
      call fft_forward(Nx, Ny, Nz, abs_uz, abs_uzh)
      


      do 42 i = 1, NxC
        do 43 j = 1, Ny
            do 44 k = 1, Nz
                ikauxh = iota * kx(i) * abs_uxh(i,j,k)
                ikauyh = iota * ky(j) * abs_uyh(i,j,k)
                ikauzh = iota * kz(k) * abs_uzh(i,j,k)
                ikauh(i,j,k)=ikauxh+ikauyh+ikauzh                
 44         continue
 43     continue
 42   continue
 


      call fft_backward(Nx, Ny, Nz, ikauh, Curve)

      Curve=-1.d0*Curve*norm


      end subroutine compute_curvature
c**********************************************************************
