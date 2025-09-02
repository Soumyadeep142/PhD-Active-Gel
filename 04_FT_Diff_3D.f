      implicit double precision(a-h, o-z)
      include 'fftw3.f'
      
      integer Nx, Ny, Nz, NT, snap_int
      real*8 pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm
      double precision :: ksq
      parameter(Nx=64, Ny=64, Nz=64, NT=500, snap_int=10)
      parameter(pi=acos(-1.d0), dt=0.001d0)
      parameter(alpha=0.1d0)
      parameter(Lx=2.d0*pi, Ly=2.d0*pi, Lz=2.d0*pi)
      double precision, allocatable :: kx(:), ky(:), kz(:)
      double precision, allocatable :: u(:,:,:)
      double precision, allocatable :: uc(:,:,:)
      complex*16, allocatable :: uh(:,:,:), uhc(:,:,:)
      integer*8 :: plan_fwd, plan_bwd
      integer :: NxC, it, isnap, l

c     Wavenumbers
      dx=Lx/Nx
      dy=Ly/Ny
      dz=Lz/Nz
c      write(*,*)dx, dy,dz
      Nxc=Nx/2+1
      write(*,*) NxC
      norm=1.d0/(dble(Nx*Ny*Nz))
      
      bx=2.0d0*pi/Lx
      by=2.0d0*pi/Ly
      bz=2.0d0*pi/lz
      
      allocate(kx(Nxc), ky(Ny), kz(Nz))
      
      do 2 i=1,Nz
      	if (i.le.(Nz/2+1)) then
      		kz(i)=bz*dble(i-1)
      	else
      		kz(i)=bz*dble(i-1-Nz)
      	endif
 2    continue
 
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
      
      u(1,1,1)=1.d0/(dx*dy*dz)

c    FFTW Plans      
      allocate(uh(NxC, Ny, Nz))
      allocate(uhc(NxC, Ny, Nz))
      call dfftw_plan_dft_r2c_3d(plan_fwd,Nx,Ny,Nz,u,uh,FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_3d(plan_bwd,Nx,Ny,Nz,uhc,uc,FFTW_ESTIMATE)
      
      open(8, file='u_3d_data.dat', status='replace')
      open(9, file='kx_mode.dat', status='replace')
      open(10, file='ky_mode.dat', status='replace')
      open(11, file='kz_mode.dat', status='replace')
      isnap=0
      
      call dfftw_execute_dft_r2c(plan_fwd, u, uh)

c    Time Marching
      do 3 it=1,NT

c    FFTW initial condition      

         do 4 i=1,NxC
          do 17 j=1,Ny
           do 18 k=1,Nz
            ksq=kx(i)**2+ky(j)**2+kz(k)**2
            uh(i,j,k)=uh(i,j,k)*(1.d0-alpha*ksq*dt)
 18        continue
 17       continue
 4       continue
      uhc=uh
      
                   
      if (mod(it, snap_int).eq.0) then
         isnap=isnap+1
         file_t=it*dt

         write(9,*) file_t, abs(uhc(3,3,4))
c         write(10,*) file_t, abs(uhc(i,j,k))
c         write(11,*) file_t, abs(uhc(i,j,k))

         
         call dfftw_execute_dft_c2r(plan_bwd, uhc, uc)
         uc=uc*norm
         write(8, '(2F20.5)') file_t, uc(1,1,1)

         
         
      endif
          
 3    continue
      close(8)
      close(9)
      close(10)
      close(11)
      call dfftw_destroy_plan(plan_fwd)
      call dfftw_destroy_plan(plan_bwd)
      call dfftw_cleanup()
      
      
      
      deallocate(u,kx,ky,kz,uh,uc,uhc)
      stop
      end
