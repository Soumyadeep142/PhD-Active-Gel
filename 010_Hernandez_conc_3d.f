c Replicating Hernandez 1

      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer  NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv, nu
      double precision bx,by,bz, c_x, c_y,c_z, eps, r, N1, N2, N3, h
      parameter(Nx=16,Ny=16,Nz=16,NT=200000,snap_int=1000,c_max=2.d0)
      parameter(pi=acos(-1.d0), dt=0.01d0, grutol=0.2, lambda=0.1d0)
      parameter(alpha=0.1d0, r0=6,c_h=2.d0,Dc=0.1d0,As=0.5d0,Af=0.5d0)
      parameter(Lx=10.d0*pi, Ly=10.d0*pi, Lz=10.d0*pi, c_ini=20)

      double precision, allocatable :: u(:,:,:), c(:,:,:), fk(:,:,:)
      double precision, allocatable :: kx(:), ky(:), kz(:),ikux(:,:,:)
      double precision, allocatable :: ikuy(:,:,:), ikuz(:,:,:)
      double precision, allocatable :: iku(:,:,:)
 
      complex*16, allocatable :: ck(:,:,:), fkh(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ikuhy(:,:,:), ikuhz(:,:,:), uh(:,:,:)

      integer*8 :: plan
      integer :: NxC, Ny, Nz, it, isnap, l
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
      allocate(fk(Nx, Ny, Nz))
      allocate(ikux(Nx, Ny, Nz))
      allocate(ikuy(Nx, Ny, Nz))
      allocate(ikuz(Nx, Ny, Nz))
      allocate(iku(Nx, Ny, Nz))
  
      
              
      allocate(ck(NxC, Ny, Nz))
      allocate(fkh(NxC, Ny, Nz))
      allocate(uh(NxC, Ny, Nz))
      allocate(ikuhx(NxC, Ny, Nz))
      allocate(ikuhy(NxC, Ny, Nz))
      allocate(ikuhz(NxC, Ny, Nz))
      

      
C***************Initiation of fields*************************
c Initial phi field
      c_x=(Nx*1.d0+1)/2
      c_y=(Ny*1.d0+1)/2
      c_z=(Nz*1.d0+1)/2
      
       do 13 i=1,Nx
        do 14 j=1,Ny
         do 15 k=1,Nz

          r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
          u(i,j,k)=-tanh(2*(r-r0))
 15     continue
 14    continue
 13   continue



c    FFTW Plans      
      call fft_forward(Nx, Ny, Nz, u, uh)
  

      do 42 i=1,NxC
       do 117 j=1,Ny
        do 118 k=1,Nz      
         ikuhx(i,j,k)=iota*kx(i)*uh(i,j,k)
         ikuhy(i,j,k)=iota*ky(j)*uh(i,j,k)
         ikuhz(i,j,k)=iota*kz(k)*uh(i,j,k)
 118     continue
 117    continue
 42    continue
      call fft_backward(Nx, Ny, Nz, ikuhx, ikux)     
      call fft_backward(Nx, Ny, Nz, ikuhy, ikuy)
      call fft_backward(Nx, Ny, Nz, ikuhz, ikuz)
      
      ikux=ikux*norm
      ikuy=ikuy*norm
      ikuz=ikuz*norm
     
       do 44 i=1,Nx
       do 45 j=1,Ny
        do 46 k=1,Nz
           iku(i,j,k)=ikux(i,j,k)**2+ikuy(i,j,k)**2+ikuz(i,j,k)**2 
 46     continue
 45    continue
 44    continue
                       
c    Initial Concentration Field
      do 23 i=1,Nx
       do 24 j=1,Ny
        do 25 k=1,Nz
          h=(dble(k)-c_z)
          r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
c          if ((h.ge.0).and.(abs(iku(i,j,k)).ge.grutol)) then
          if ((h.ge.0).and.(r.ge.(r0-0.5)).and.(r.le.(r0+0.5))) then
              c(i,j,k)=c_ini
          else
              c(i,j,k)=0
          endif
c           c(i,j,k)=1
 25     continue
 24    continue
 23    continue
       

c**********Saving the initial file*****************************
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm

      open(8, file='step_0.dat')

      do 41 i=1,Nx
        do 52 j=1, Ny
         do 61 k=1,Nz
          write(8,*)i, j, k, u(i,j,k), c(i,j,k)     
 61     continue
 52    continue
 41   continue
      close(8)


      call fft_forward(Nx, Ny, Nz, c, ck)

      


c**************************************************************


c*****************Time Marching********************************

      do 3 it=1, NT
    
      call fft_backward(Nx, Ny, Nz, ck, c)
      c=c*norm

       
       do 19 i=1,Nx
        do 20 j=1,Ny
         do 22 k=1, Nz

           fk(i,j,k)=Dc*(2*As*((u(i,j,k)**2-1)**2)*((c(i,j,k)-c_max)
     1       *c(i,j,k)*(2*c(i,j,k)-c_max))+Af*2*u(i,j,k)**2*c(i,j,k))
 22      continue    
 20     continue    
 19    continue
 
       call fft_forward(Nx, Ny, Nz, fk, fkh)
       call fft_forward(Nx, Ny, Nz, c, ck)
 
       do 26 i=1, NxC
         do 27 j=1, Ny
          do 28 k=1, Nz      	    
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)
      	    nu=Dc*lambda*k2*2*As
            ck(i,j,k)=exp(-k2*nu*dt)*(ck(i,j,k)-k2*fkh(i,j,k)*dt)
 28      continue
 27     continue
 26    continue

c*****************************


       
      

c********Saving files*******

      if (mod(it,snap_int).eq.0) then      
      
       
       call fft_backward(Nx, Ny, Nz, ck, c)
       c=c*norm
      
       write(fname, '("step_", I0, ".dat")') it
       open(8, file=fname, status='replace')      
       do 1 i=1,Nx
        do 4 j=1,Ny
         do 5 k=1, Nz
          write(8,*)i, j, k, u(i,j,k), c(i,j,k)
 5       continue
 4      continue
 1     continue
       close(8)
       

c       call fft_forward(Nx, u, uh)
       call fft_forward(Nx, Ny, Nz, c, ck)
      endif
      

 3    continue
c *******************************************************
      call dfftw_cleanup()
      stop
      end
      
      
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
