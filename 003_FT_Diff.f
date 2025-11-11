      implicit double precision(a-h, o-z)
      include 'fftw3.f'
      
      integer N, NT, snap_int
      real*8 pi, kappa, dt,L,dx
      
      parameter(N=32, NT=500, snap_int=1)
      parameter(pi=acos(-1.d0), dt=0.01d0)
      parameter(L=2*pi,x_ini=0, alpha=0.1d0)
      double precision, allocatable :: x_array(:), kvec(:), u(:)
      double precision, allocatable :: uc(:)
      complex*16, allocatable :: uh(:), uhc(:), uh_snap(:,:)
      integer*8 :: plan_fwd, plan_bwd, it
      character(len=50) :: fname
c     Creating x-arrays

      allocate(x_array(N))
      allocate(u(N))
      allocate(uc(N))
      x=x_ini
      dx=L/N  
      do 1 i=1,N
c        x=x_ini+(i-1)*dx
c        x_array(i)=x
c        u(i)=dexp(-(x-L/2.d0)**2/0.1d0)
        if (i.eq.(N/2+1)) then
         u(i)=1.d0/dx
        else
         u(i)=0.d0
        endif
 1    continue
c      write(*,*) u(N)       
      
c     Wavenumbers
      allocate(kvec(N/2+1))
      do 2 i=1,N/2+1
       kvec(i)=2.d0*pi*(i-1)/(N*dx)
 2    continue
      
      allocate(uh(N/2+1))
      allocate(uhc(N/2+1))
c    FFTW Plans
      call dfftw_plan_dft_r2c_1d(plan_fwd, N, u, uh, FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_1d(plan_bwd, N, uh, u, FFTW_ESTIMATE)      
      
      call dfftw_execute_dft_r2c(plan_fwd, u, uh)
      uhc=uh
c      write(*,*) uhc
      call dfftw_destroy_plan(plan_fwd)      
c    Time Marching
      T=0.d0
      isnap=0
      nsnaps = NT / snap_int
      allocate(uh_snap(N/2+1, nsnaps))
      do 3 it=1,NT
         T=T+dt
         do 4 j=1,N/2+1
           uh(j) = uh(j) * (1.d0-alpha * (kvec(j)**2) * dt)
c           uh(j) = (uhc(j)) * dexp(- alpha * (kvec(j)**2) * T)
 4       continue
      if (mod(it, snap_int).eq.0) then
         isnap=isnap+1
         uh_snap(:, isnap)=uh
      endif
      

 3    continue
      
      do 7 isnap=1, nsnaps
       uh(:)=uh_snap(:, isnap)

       open(9, file='uh_mode.dat')
       open(8, file='u_data.dat')
       do 17 i=1,N/2+1
        file_t=isnap*snap_int*dt
        if (kvec(i).eq.5) then
         
         write(9,*)file_t, real(uh(i)), aimag(uh(i)), abs(uh(i))
        endif
 17    continue  
      
       call dfftw_plan_dft_c2r_1d(plan_bwd, N, uh, uc, FFTW_ESTIMATE)
       call dfftw_execute_dft_c2r(plan_bwd, uh, uc)
       call dfftw_destroy_plan(plan_bwd)
       
      uc=uc/N
 
c     Write the output file
c       write(fname,'("u_t",I4.4,".dat")') isnap*snap_int
c       open(8, file=fname)
c       do 6 i=1,N
c         xc=x_array(i)-L/2
c         if (xc.eq.0) then 
c         write(*,*)xc,i
         write(8,'(2F20.5)')file_t, uc(N/2+1)
c         endif
c 6     continue      
c       close(8)
       

 7    continue


      call dfftw_cleanup()
      
      
      
      deallocate(x_array, u, kvec, uh)      
      stop
      end
