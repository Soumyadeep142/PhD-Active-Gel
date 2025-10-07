c This code does concentration,velocity as a field only

      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer Nx, Ny, Nz, NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv
      double precision bx,by,bz, c_x, c_y,c_z, eps, r, N1, N2, N3
      parameter(Nx=16, Ny=16, Nz=16, NT=5000, snap_int=100)
      parameter(pi=acos(-1.d0), dt=0.0001d0)
      parameter(alpha=0.1d0, r0=4)
      parameter(Lx=2.d0*pi, Ly=2.d0*pi, Lz=2.d0*pi)
      double precision, allocatable :: kx(:), ky(:), kz(:), cvx(:,:,:)
      double precision, allocatable :: u(:,:,:), c(:,:,:), cvy(:,:,:)
      double precision, allocatable :: ikuc(:,:,:), zi(:,:,:)
      double precision, allocatable :: uc(:,:,:), iku(:,:,:), cvz(:,:,:)
      double precision, allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
      double precision, allocatable :: ikcvx(:,:,:), ikcvy(:,:,:)
      double precision, allocatable :: ikcvz(:,:,:), k2c(:,:,:)
      double precision, allocatable :: ikux(:,:,:), ikuy(:,:,:)
      double precision, allocatable :: ikuz(:,:,:), N(:,:,:), RHS(:,:,:)
      double precision, allocatable :: ikcx(:,:,:), ikcy(:,:,:)
      double precision, allocatable :: ikcz(:,:,:), kzx(:,:,:)
      double precision, allocatable :: kzy(:,:,:), kzz(:,:,:)
      double precision, allocatable :: ikcv(:,:,:), cvkz(:,:,:)
      double precision, allocatable :: ikckz(:,:,:), cdtlndp(:,:,:)
      
      complex*16, allocatable :: vxc(:,:,:), vyc(:,:,:), vzc(:,:,:)
      complex*16, allocatable :: uh(:,:,:), vuh(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ph(:,:,:), ikuhy(:,:,:), ikuhz(:,:,:)
      complex*16, allocatable :: ck(:,:,:), ikcvxh(:,:,:), ikcvyh(:,:,:)
      complex*16, allocatable :: ikcvzh(:,:,:), k2ch(:,:,:), zih(:,:,:)
      complex*16, allocatable :: ikzx(:,:,:), ikzy(:,:,:), cvxh(:,:,:)
      complex*16, allocatable :: ikzz(:,:,:) , ikcxh(:,:,:), cvyh(:,:,:)
      complex*16, allocatable :: ikcyh(:,:,:), ikczh(:,:,:), cvzh(:,:,:)
      complex*16, allocatable :: cdtlndph(:,:,:), cvkzh(:,:,:)
      complex*16, allocatable :: ikckzh(:,:,:), RHSh(:,:,:)
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
      allocate(iku(Nx,Ny,Nz))
      allocate(ikux(Nx,Ny,Nz))
      allocate(ikuy(Nx,Ny,Nz))
      allocate(ikuz(Nx,Ny,Nz))
      allocate(vx(Nx, Ny, Nz))
      allocate(vy(Nx, Ny, Nz))
      allocate(vz(Nx, Ny, Nz))
      allocate(N(Nx,Ny,Nz))
      allocate(c(Nx,Ny,Nz))
      allocate(ikcvx(Nx,Ny,Nz))
      allocate(ikcvy(Nx,Ny,Nz))
      allocate(ikcvz(Nx,Ny,Nz))                  
      allocate(zi(Nx, Ny, Nz))      
      allocate(kzx(Nx, Ny, Nz))
      allocate(kzy(Nx, Ny, Nz))
      allocate(kzz(Nx, Ny, Nz))      
      allocate(cvkz(Nx, Ny, Nz))
      allocate(ikcv(Nx, Ny, Nz))
      allocate(ikcz(Nx, Ny, Nz))
      allocate(k2c(Nx, Ny, Nz))
      allocate(ikcx(Nx, Ny, Nz))
      allocate(ikcy(Nx, Ny, Nz))
      allocate(ikckz(Nx, Ny, Nz))   
      allocate(cdtlndp(Nx, Ny, Nz))
      allocate(RHS(Nx, Ny, Nz))
      allocate(cvx(Nx, Ny, Nz))
      allocate(cvy(Nx, Ny, Nz))     
      allocate(cvz(Nx, Ny, Nz))    
      
      allocate(ikuhx(NxC,Ny,Nz))
      allocate(ikuhy(NxC,Ny,Nz))
      allocate(ikuhz(NxC,Ny,Nz))                 
      allocate(vxc(NxC, Ny, Nz))
      allocate(vyc(NxC, Ny, Nz))
      allocate(vzc(NxC, Ny, Nz))
      allocate(uh(NxC, Ny, Nz))
      allocate(vuh(NxC, Ny, Nz))                  
      allocate(ck(NxC,Ny,Nz))
      allocate(ikcvxh(NxC, Ny, Nz))
      allocate(ikcvyh(NxC, Ny, Nz))
      allocate(ikcvzh(NxC, Ny, Nz))            
      allocate(zih(NxC, Ny, Nz))          
      allocate(ikzx(NxC, Ny, Nz))
      allocate(ikzy(NxC, Ny, Nz))
      allocate(ikzz(NxC, Ny, Nz))
      allocate(k2ch(NxC, Ny, Nz))
      allocate(ikcxh(NxC, Ny, Nz))
      allocate(ikcyh(NxC, Ny, Nz))
      allocate(ikczh(NxC, Ny, Nz))
      allocate(cdtlndph(NxC, Ny, Nz))
      allocate(cvxh(NxC, Ny, Nz))
      allocate(cvyh(NxC, Ny, Nz))
      allocate(cvzh(NxC, Ny, Nz))
      allocate(cvkzh(NxC, Ny, Nz))
      allocate(ikckzh(NxC, Ny, Nz))
      allocate(RHSh(NxC, Ny, Nz))
      
C***************Initiation of fields*************************
c Initial phi field
      c_x=(Nx*1.d0+1)/2
      c_y=(Ny*1.d0+1)/2
      c_z=(Nz*1.d0+1)/2      
       do 13 i=1,Nx
        do 14 j=1,Ny
         do 15 k=1,Nz
          r=sqrt((dble(i)-c_x)**2+(dble(j)-c_y)**2+(dble(k)-c_z)**2)
          u(i,j,k)=0.5d0*(1.d0-tanh(r-r0))
 15      continue
 14    continue
 13   continue
 
c    FFTW Plans      
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
      
      iku=ikux**2+ikuy**2+ikuz**2            

      open(8, file='step_0.dat')

      do 41 i=1,Nx
        do 52 j=1, Ny
         do 61 k=1,Nz
          write(8,*)i, j, k, u(i,j,k)     
 61     continue
 52    continue
 41   continue
      close(8)
      

      
c   Initial Velocity Field
      do 9 i=1,Nx
       do 10 j=1,Ny
        do 16 k=1,Nz
           vx(i,j,k)=0
           vy(i,j,k)=0
           vz(i,j,k)=0
 16     continue
 10    continue
 9    continue
 

c    Initial Concentration Field
      do 23 i=1,Nx
       do 24 j=1,Ny
        do 25 k=1,Nz
           c(i,j,k)=1
 25     continue
 24    continue
 23    continue

      call fft_forward(Nx, Ny, Nz, c, ck)
      call fft_forward(Nx, Ny, Nz, vx, vxc)
      call fft_forward(Nx, Ny, Nz, vy, vyc)
      call fft_forward(Nx, Ny, Nz, vz, vzc)
c**************************************************************


c*****************Time Marching********************************

      do 3 it=1, NT

c   Calculating |grad phi|                  
      do 4 i=1,NxC
       do 17 j=1,Ny
        do 18 k=1,Nz
         ikuhx(i,j,k)=iota*kx(i)*uh(i,j,k)
         ikuhy(i,j,k)=iota*ky(j)*uh(i,j,k)
         ikuhz(i,j,k)=iota*kz(k)*uh(i,j,k)     
 18     continue
 17    continue
 4    continue
 
      call fft_backward(Nx, Ny, Nz, ikuhx, ikux)     
      call fft_backward(Nx, Ny, Nz, ikuhy, ikuy)
      call fft_backward(Nx, Ny, Nz, ikuhz, ikuz)
      
      ikux=ikux*norm
      ikuy=ikuy*norm
      ikuz=ikuz*norm
      
      iku=ikux**2+ikuy**2+ikuz**2
      
      do 19 i=1,Nx
      	do 20 j=1,Ny
      	  do 22 k=1,Nz
      	    N1=vx(i,j,k)*ikux(i,j,k)
      	    N2=vy(i,j,k)*ikuy(i,j,k)
      	    N3=vz(i,j,k)*ikuz(i,j,k)
      	    N(i,j,k)=N1+N2+N3
 22       continue
 20     continue
 19   continue
 
      call fft_forward(Nx, Ny, Nz, N, vuh)
      
      do 5 i=1,NxC
      	do 6 j=1,Ny
      	 do 7 k=1,Nz
      	   uh(i,j,k)=uh(i,j,k)-vuh(i,j,k)*dt
 7       continue
 6      continue
 5    continue

c*********Terms***********

       do 29 i=1,Nx
        do 30 j=1,Ny
         do 32 k=1,Nz
           zi(i,j,k)=log(iku(i,j,k))
 32      continue
 30     continue
 29    continue
        
      
      
      call fft_forward(Nx, Ny,  Nz, zi, zih)
      
    
      do 26 i=1,NxC
      	do 27 j=1,Ny
      	  do 28 k=1,Nz
      	  
      	    k2=-(kx(i)**2+ky(j)**2+kz(k)**2)
      	    cdtlndph(i,j,k)=k2*ck(i,j,k)*2*uh(i,j,k)*vuh(i,j,k)
      	 
            ikcvxh(i,j,k)=iota*ck(i,j,k)*kx(i)*vxc(i,j,k)
            ikcvyh(i,j,k)=iota*ck(i,j,k)*ky(j)*vyc(i,j,k)
            ikcvzh(i,j,k)=iota*ck(i,j,k)*kz(k)*vzc(i,j,k)
            
            k2ch(i,j,k)=-(kx(i)**2+ky(j)**2+kz(k)**2)*ck(i,j,k)
            
            ikzx(i,j,k)=iota*kx(i)*zih(i,j,k)
            ikzy(i,j,k)=iota*ky(j)*zih(i,j,k)
            ikzz(i,j,k)=iota*kz(k)*zih(i,j,k)
            
            ikcxh(i,j,k)=iota*kx(i)*ck(i,j,k)
            ikcyh(i,j,k)=iota*ky(j)*ck(i,j,k)
            ikczh(i,j,k)=iota*kz(k)*ck(i,j,k)
            
            cvxh(i,j,k)=ck(i,j,k)*vxc(i,j,k)
            cvyh(i,j,k)=ck(i,j,k)*vyc(i,j,k)
            cvzh(i,j,k)=ck(i,j,k)*vzc(i,j,k)
 28       continue
 27     continue
 26    continue

       call fft_backward(Nx, Ny, Nz, cdtlndph, cdtlndp)
       
       call fft_backward(Nx, Ny, Nz, ikcvxh, ikcvx)
       call fft_backward(Nx, Ny, Nz, ikcvyh, ikcvy)
       call fft_backward(Nx, Ny, Nz, ikcvzh, ikcvz)
       
       call fft_backward(Nx, Ny, Nz, k2ch, k2c)
       
       call fft_backward(Nx, Ny, Nz, ikzx, kzx)
       call fft_backward(Nx, Ny, Nz, ikzy, kzy)
       call fft_backward(Nx, Ny, Nz, ikzz, kzz)
       
       call fft_backward(Nx, Ny, Nz, ikcxh, ikcx)
       call fft_backward(Nx, Ny, Nz, ikcyh, ikcy)
       call fft_backward(Nx, Ny, Nz, ikczh, ikcz)
       
       call fft_backward(Nx, Ny, Nz, cvxh, cvx)
       call fft_backward(Nx, Ny, Nz, cvyh, cvy)
       call fft_backward(Nx, Ny, Nz, cvzh, cvz)

       cdtlndp=cdtlndp*norm
             
       ikcvx=ikcvx*norm
       ikcvy=ikcvy*norm
       ikcvz=ikcvz*norm
       
       k2c=k2c*norm
       
       kzx=kzx*norm
       kzy=kzy*norm
       kzz=kzz*norm    
          
       ikcx=ikcx*norm
       ikcy=ikcy*norm
       ikcz=ikcz*norm
       
       cdtlndp=cdtlndp/iku
       
       cvx=cvx*norm
       cvy=cvy*norm
       cvz=cvz*norm

       do 33 i=1, Nx
       do 34 j=1, Ny
       do 35 k=1, Nz
       	
       ikcv(i,j,k)=sqrt(ikcvx(i,j,k)**2+ikcvy(i,j,k)**2+ikcvz(i,j,k)**2)
       cvkz(i,j,k)=cvx(i,j,k)*kzx(i,j,k)+cvy(i,j,k)*kzy(i,j,k)+
     1               cvz(i,j,k)*kzz(i,j,k)     
       ikckz(i,j,k)=ikcx(i,j,k)*kzx(i,j,k)+ikcy(i,j,k)*kzy(i,j,k)
     1                +ikcz(i,j,k)*kzz(i,j,k)       
 35    continue
 34    continue
 33    continue

       do 36 i=1,Nx
       	do 37 j=1,Ny
       	  do 38 k=1,Nz
            RHS(i,j,k)=cdtlndp(i,j,k)+ikcv(i,j,k)+k2c(i,j,k)+cvkz(i,j,k)
     1                 +ikckz(i,j,k)  
 38       continue
 37     continue
 36    continue
 
c   Time Marching in Concentration field
 
       call fft_forward(Nx, Ny, Nz, RHS, RHSh)
       
       do 39 i=1, NxC
        do 40 j=1, Ny
          do 42 k=1, Nz
           ck(i,j,k)=ck(i,j,k)+RHSk(i,j,k)
 42       continue
 40     continue
 39    continue
c*****************************


       
      

c********Saving files*******

      if (mod(it,snap_int).eq.0) then      
       
       call fft_backward(Nx, Ny, Nz, uh, u)
       u=u*norm
       
       call fft_backward(Nx, Ny, Nz, ck, c)
       c=c*norm
      
       write(fname, '("step_", I0, ".dat")') it
       open(8, file=fname, status='replace')      
       do 1 i=1,Nx
        do 2 j=1, Ny
         do 31 k=1,Nz
          write(8,*)i, j, k, u(i,j,k), c(i,j,k)
 31      continue
 2      continue
 1     continue
       close(8)
       

       call fft_forward(Nx, Ny, Nz, u, uh)
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

