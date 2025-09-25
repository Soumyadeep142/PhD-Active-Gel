      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer Nx, Ny, Nz, NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv
      double precision bx,by,bz, c_x, c_y,c_z, eps, c, r
      parameter(Nx=16, Ny=16, Nz=16, NT=5000, snap_int=100)
      parameter(pi=acos(-1.d0), dt=0.0001d0)
      parameter(alpha=0.1d0, r0=4)
      parameter(Lx=2.d0*pi, Ly=2.d0*pi, Lz=2.d0*pi)
      double precision, allocatable :: kx(:), ky(:), kz(:)
      double precision, allocatable :: u(:,:,:),p(:,:,:)
      double precision, allocatable :: pc(:,:,:), ikuc(:,:,:)
      double precision, allocatable :: uc(:,:,:), v(:), iku(:,:,:)
      double precision, allocatable :: ikux(:,:,:), ikuy(:,:,:)
      double precision, allocatable :: ikuz(:,:,:)
      complex*16, allocatable :: uh(:,:,:), uhc(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ph(:,:,:), ikuhy(:,:,:), ikuhz(:,:,:)
      integer*8 :: plan
      integer :: NxC, it, isnap, l
      complex*16, parameter :: iota = (0.d0, 1.d0)
      complex*16, parameter :: minus_iota = (0.d0, -1.d0)
      character(len=64) :: fname

      call execute_command_line('rm -f step_*.dat')
          
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
      
      allocate(kx(NxC), ky(Ny), kz(Nz))
      
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

      do 12 i=1, NxC
      	kx(i)=bx*dble(i-1)
 12   continue
 
 
c     Creating initial condition-arrays
      allocate(u(Nx, Ny, Nz))
      allocate(uc(Nx,Ny,Nz))
      allocate(p(Nx, Ny, Nz))
      allocate(iku(Nx,Ny,Nz))
      allocate(ikux(Nx,Ny,Nz))
      allocate(ikuy(Nx,Ny,Nz))
      allocate(ikuz(Nx,Ny,Nz))
      allocate(pc(Nx, Ny, Nz))
      allocate(ikuc(Nx,Ny,Nz))
      allocate(ikuhx(NxC,Ny,Nz))
      allocate(ikuhy(NxC,Ny,Nz))
      allocate(ikuhz(NxC,Ny,Nz))

      u=0.d0
      p=0.d0
      pc=0.d0
      iku=0.d0
      ikuc=0.d0
      

c Initial phi field
      c_x=(Nx*1.d0+1)/2
      c_y=(Ny*1.d0+1)/2
      c_z=(Nz*1.d0+1)/2      
      do 13 i=1,Nx
       do 14 j=1,Ny
        do 15 k=1,Nz
         r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
         u(i,j,k)=0.5d0*(1.d0-tanh(r-r0))         
 15     continue
 14    continue
 13   continue
     
            

c    FFTW Plans      
      allocate(uh(NxC, Ny, Nz))
      allocate(uhc(NxC, Ny, Nz))
      allocate(ph(NxC, Ny, Nz))

      call fft_forward(Nx, Ny, Nz, u, uh)

      
      
c    FFTW initial condition      

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
       
       iku=dsqrt(ikux**2+ikuy**2+ikuz**2)
       
c Initial psi field
      open(8, file='step_0.dat')
      c=1
      p=iku**2*c
      do 41 i=1,Nx
        do 52 j=1, Ny
         do 61 k=1,Nz
          write(8,*)i, j, k, c, p(i,j,k), u(i,j,k), iku(i,j,k)      
 61     continue
 52    continue
 41   continue
      close(8)
      
      call fft_forward(Nx, Ny, Nz, p, ph)  
      
c   Velocity 
      allocate(v(3))
      v(1)=10
      v(2)=0
      v(3)=0
      
      eps=0.d0
c   Time Marching
      open(8,file='adv_c', status='replace')   
      do 3 it=1, NT
      

c    FFTW initial condition      

      do 4 i=1,NxC
       do 17 j=1,Ny
        do 18 k=1,Nz
      
         kv=kx(i)*v(1)+ky(j)*v(2)+kz(k)*v(3)
         uh(i,j,k)=uh(i,j,k)*(1.d0-iota*kv*dt)
         ph(i,j,k)=ph(i,j,k)*(1.d0-iota*kv*dt)
         ikuhx(i,j,k)=iota*kx(i)*uh(i,j,k)
         ikuhy(i,j,k)=iota*ky(j)*uh(i,j,k)
         ikuhz(i,j,k)=iota*kz(k)*uh(i,j,k)

 18     continue
 17    continue
 4    continue


c   Saving files

      if (mod(it,snap_int).eq.0) then      
       call fft_backward(Nx, Ny, Nz, ph, p)
       p=p*norm
      
       call fft_backward(Nx, Ny, Nz, ikuhx, ikux)     
       call fft_backward(Nx, Ny, Nz, ikuhy, ikuy)
       call fft_backward(Nx, Ny, Nz, ikuhz, ikuz)
       
       ikux=ikux*norm
       ikuy=ikuy*norm
       ikuz=ikuz*norm
       
       iku=dsqrt(ikux**2+ikuy**2+ikuz**2)
       
       call fft_backward(Nx, Ny, Nz, uh, u)
       u=u*norm
      
       write(fname, '("step_", I0, ".dat")') it
       open(8, file=fname, status='replace')      
       do 1 i=1,Nx
        do 2 j=1, Ny
         do 31 k=1,Nz
          d=abs(iku(i,j,k))
          c=p(i,j,k)/((d*d))        
          write(8,*)i, j, k, c, p(i,j,k), u(i,j,k), iku(i,j,k)
 31      continue
 2      continue
 1     continue
       close(8)
       
       call fft_forward(Nx, Ny, Nz, p, ph)     
       call fft_forward(Nx, Ny, Nz, u, uh)
      endif
      
      
 3    continue
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

