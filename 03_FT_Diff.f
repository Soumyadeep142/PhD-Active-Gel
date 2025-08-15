      implicit double precision(a-h, o-z)
      include 'fftw3.f'
      
      integer N, NT, snap_int
      real*8 pi, kappa, dt,L,dx
      
      parameter(N=128, NT=500, snap_int=10)
      parameter(pi=acos(-1.d0), kappa=0.1, dt=0.01)
      parameter(L=2*pi,x_ini=0., alpha=0.1d0)
      double precision, allocatable :: x_array(:), kvec(:), u(:)
      double precision, allocatable :: u_copy(:)
      complex*16, allocatable :: uh(:), uh_snap(:,:)
      integer*8 :: plan_fwd, plan_bwd, it
      character(len=50) :: fname
c     Creating x-arrays

      allocate(x_array(N))
      allocate(u(N))
      x=x_ini
      dx=L/N  
      do 1 i=1,N
        x=x_ini+(i-1)*dx
        x_array(i)=x
        u(i)=dexp(-(x-L/2.d0)**2/0.1d0)
 1    continue
      u_copy=u

         
      
      
c     Wavenumbers
      allocate(kvec(N/2))
      do 2 i=1,N/2+1
       kvec(i)=2.d0*pi*(i-1)/(N*dx)
 2    continue
      
      allocate(uh(N/2+1))
c    FFTW Plans
      call dfftw_plan_dft_r2c_1d(plan_fwd, N, u, uh, FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_1d(plan_bwd, N, uh, u, FFTW_ESTIMATE)
      
c    FFTW initial condition      
      call dfftw_execute_dft_r2c(plan_fwd, u, uh)
c    Time Marching
      T=0.d0
      isnap=0
      nsnaps = NT / snap_int
      allocate(uh_snap(N/2+1, nsnaps))
      do 3 it=1,NT
         T=T+dt
         do 4 j=1,N/2+1
           decay=dexp(-alpha*kvec(j)**2*dt)
           uh(j)=decay*uh(j)
 4       continue
      if (mod(it, snap_int).eq.0) then
         isnap=isnap+1
         uh_snap(:, isnap)=uh
      endif
 3    continue
      
      do 7 isnap=1, nsnaps
       uh(:)=uh_snap(:, isnap)
       call dfftw_execute_dft_c2r(plan_bwd, uh, u)

c      Normalise
       do 5 i=1,N
         U(i)=U(i)/dble(N)
 5     continue
 
c     Write the output file
       write(fname,'("u_t",I4.4,".dat")') isnap*snap_int
       open(8, file=fname)
       do 6 i=1,N
         write(8,'(2E20.10)')x_array(i), u(i)
 6     continue      
       close(8)
 7     continue
      call dfftw_destroy_plan(plan_fwd)
      call dfftw_destroy_plan(plan_bwd)
      call dfftw_cleanup()
      
      
      
      deallocate(x_array, u, u_copy, kvec, uh)      
      stop
      end
