c Replicating Hernandez 1

      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer  NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv, nu
      double precision bx,by,bz, c_x, c_y,c_z, eps, r, N1, N2, N3, h
      parameter(Nx=32, NT=50000, snap_int=100, c_max=2.d0)
      parameter(pi=acos(-1.d0), dt=0.01d0, grutol=1e-1, lambda=5.d0)
      parameter(alpha=0.1d0, r0=8, c_h=0.d0, Dc=0.2d0,As=0.1d0,Af=0.2d0)
      parameter(Lx=10.d0*pi)

      double precision, allocatable :: u(:), c(:), fk(:),  kx(:)
 
      complex*16, allocatable :: ck(:), fkh(:)

      integer*8 :: plan
      integer :: NxC, it, isnap, l
      complex*16, parameter :: iota = (0.d0, 1.d0)
      complex*16, parameter :: minus_iota = (0.d0, -1.d0)
      character(len=64) :: fname

      call execute_command_line('rm -f step_*.dat')
          
c     Wavenumbers
      dx=Lx/Nx


      NxC=Nx/2+1

      norm=1.d0/(dble(Nx))
c      write(*,*)norm
      
      bx=2.0d0*pi/Lx

      
      allocate(kx(NxC))


      do 12 i=1, NxC
      	kx(i)=bx*dble(i-1)
 12   continue
      write(*,*) kx
 
c     Creating initial condition-arrays
      allocate(u(Nx))
      allocate(c(Nx))
      allocate(fk(Nx))
  
      
              
      allocate(ck(NxC))
      allocate(fkh(NxC))

      
C***************Initiation of fields*************************
c Initial phi field
      c_x=(Nx*1.d0+1)/2
    
       do 13 i=1,Nx

          r=sqrt((dble(i)-c_x)**2)
          u(i)=-tanh(2*(r-r0))

 13   continue


c    Initial Concentration Field
      do 23 i=1,Nx
c          if (i.gt.Nx/2) then     
c          	c(i)=1.d0
c          else
c          	c(i)=0.5
c          endif
          c(i)=1.d0
 23    continue

c**********Saving the initial file*****************************
      open(8, file='step_0.dat')

      do 41 i=1,Nx

          write(8,*)i, u(i), c(i)     

 41   continue
      close(8)


      call fft_forward(Nx, c, ck)
      
c initialize previous RHSh to zero


c**************************************************************


c*****************Time Marching********************************

      do 3 it=1, NT
    
      call fft_backward(Nx, ck, c)
      c=c*norm

       
       do 19 i=1,Nx

           fk(i)=Dc*(2*As*((u(i)**2-1)**2)*((c(i)-c_max)*c(i)
     1           *(2*c(i)-c_max))+Af*2*u(i)**2*c(i))
     
     
 19    continue
 
       call fft_forward(Nx, fk, fkh)
       call fft_forward(Nx, c, ck)
 
       do 26 i=1,NxC      	    
      	    k2=(kx(i)**2)
      	    nu=Dc*lambda*k2*2*As
c      	    write(*,*) nu(i)
            ck(i)=exp(-k2*nu*dt)*(ck(i)-k2*fkh(i)*dt)
c            write(*,*)Dc*k2*fkh(i)*dt
 26    continue

c*****************************


       
      

c********Saving files*******

      if (mod(it,snap_int).eq.0) then      
      
       
       call fft_backward(Nx, ck, c)
       c=c*norm
      
       write(fname, '("step_", I0, ".dat")') it
       open(8, file=fname, status='replace')      
       do 1 i=1,Nx

          write(8,*)i, c(i)

 1     continue
       close(8)
       

c       call fft_forward(Nx, u, uh)
       call fft_forward(Nx, c, ck)

      endif
      

 3    continue
c *******************************************************
      call dfftw_cleanup()
      stop
      end
      
      
      subroutine fft_forward(N, in, out)
      implicit none
      include 'fftw3.f'
      integer, intent(in) :: N
      real*8, intent(inout) :: in(N)
      complex*16, intent(inout) :: out(N/2+1)
      integer*8 :: plan

      call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,in,out)
      call dfftw_destroy_plan(plan)
      end subroutine fft_forward
      
      subroutine fft_backward(N, in, out)
      implicit none
      include 'fftw3.f'
      integer, intent(in) :: N
      complex*16, intent(inout) :: in(N/2+1)
      real*8, intent(inout) :: out(N)
      integer*8 :: plan

      call dfftw_plan_dft_c2r_1d(plan,N,in,out,FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(plan, in, out)
      call dfftw_destroy_plan(plan)
      end subroutine fft_backward

