c Replicating Hernandez 1 with curvature..... real space curvature calculated , time increment..... ellipsoid eqn... protein given constant in the middle by hand, now doing couple

      implicit double precision(a-h, o-z)
      include 'fftw3.f'

      integer  NT, snap_int
      double precision pi, dt,Lx, Ly, Lz, dx, dy, dz, k2, norm, kv, nuc
      double precision bx,by,bz, c_x, c_y,c_z, eps, r, h, a_x, a_y, a_z
      double precision nuu, r_ellip, r2, F, a,lambda, lambda1, lambda2
      parameter(Nx=64,Ny=64,Nz=64,NT_i=30000,snap_int=1000, c_max=3.d0)
      parameter(pi=acos(-1.d0), dt=0.001d0, grutol=0.2, lambda=0.05d0)
      parameter(r0=10,c_h=2.d0,Dc=0.1d0,As=1.0d0,Af=5.0d0, tol=0.1d0)
      parameter(Du=1.0d0, sig=0.75d0, beta=0.0d0, lambda1=0.4d0)
      parameter(Lx=10.d0*pi, Ly=10.d0*pi, Lz=10.d0*pi, c_ini=3.0d0)
      parameter(alpha=5.0, lambda2=0.4d0, Dc1=0.01d0, NT_f=100000)

      double precision, allocatable :: u(:,:,:), c(:,:,:), fc(:,:,:)
      double precision, allocatable :: kx(:), ky(:), kz(:),ux(:,:,:)
      double precision, allocatable :: uy(:,:,:), uz(:,:,:)
      double precision, allocatable :: iku(:,:,:), fu(:,:,:)
      double precision, allocatable :: Curve(:,:,:), n_x(:,:,:)
      double precision, allocatable :: n_y(:,:,:), n_z(:,:,:)
      double precision, allocatable :: cc(:,:,:), ciku(:,:,:)
      
 
      complex*16, allocatable :: ck(:,:,:), fch(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ikuhy(:,:,:), ikuhz(:,:,:), uh(:,:,:)
      complex*16, allocatable :: fuh(:,:,:), Curveh(:,:,:), uh_l(:,:,:)
      complex*16, allocatable :: cch(:,:,:), ikbcu(:,:,:), cikuh(:,:,:)

      integer*8 :: plan
      integer :: NxC, it, isnap, l
      complex*16, parameter :: iota = (0.d0, 1.d0)
      complex*16, parameter :: minus_iota = (0.d0, -1.d0)
      character(len=1024) :: fname
      character(len=256) :: folder
      character(len=512) :: cmd
      character*32 sNx, sdt, slambda, sAs, sAf, sDc
      character*32 sDu, ssig, sl1, salpha, sl2, sDc1

      call execute_command_line('rm -f step_*.dat')
      call execute_command_line('rm -f step_*.vtu')
      call execute_command_line('rm -f cutc_step_*.dat')
      call execute_command_line('rm -f cutc_step_*.vtu')
      call execute_command_line('rm -f cute_step_*.dat')
      call execute_command_line('rm -f cute_step_*.vtu')  

c************Folder Creation*******************************
c --- convert numbers to strings ---
      cmd = 'mkdir -p 018_params' 
      call system(cmd)
      
      write(sNx,'(I10)') Nx
      write(sdt,'(F10.5)') dt
      write(slambda,'(F10.5)') lambda
      write(sAs,'(F10.5)') As
      write(sAf,'(F10.5)') Af
      write(sDc,'(F10.5)') Dc
      write(sDu,'(F10.5)') Du
      write(ssig,'(F10.5)') sig
      write(sl1,'(F10.5)') lambda1
      write(salpha,'(F10.5)') alpha
      write(sl2,'(F10.5)') lambda2
      write(sDc1,'(F10.5)') Dc1

c --- remove spaces ---
      call strip_spaces(sNx)
      call strip_spaces(sdt)
      call strip_spaces(slambda)
      call strip_spaces(sAs)
      call strip_spaces(sAf)
      call strip_spaces(sDc)
      call strip_spaces(sDu)
      call strip_spaces(ssig)
      call strip_spaces(sl1)
      call strip_spaces(salpha)
      call strip_spaces(sl2)
      call strip_spaces(sDc1)

c --- build folder name (SAFE) ---
       folder = '018_params'

       folder = trim(folder)//'/N_'//trim(adjustl(sNx))
       folder = trim(folder)//'_dt_'//trim(adjustl(sdt))
       folder = trim(folder)//'_lambda_'//trim(adjustl(slambda))
       folder = trim(folder)//'_As_'//trim(adjustl(sAs))
       folder = trim(folder)//'_Af_'//trim(adjustl(sAf))
       folder = trim(folder)//'_Dc_'//trim(adjustl(sDc))
       folder = trim(folder)//'_Dphi_'//trim(adjustl(sDu))
       folder = trim(folder)//'_sigma_'//trim(adjustl(ssig))
       folder = trim(folder)//'_lambda1_'//trim(adjustl(sl1))
       folder = trim(folder)//'_alpha_'//trim(adjustl(salpha))
       folder = trim(folder)//'_lambda2_'//trim(adjustl(sl2))
       folder = trim(folder)//'_Dc1_'//trim(adjustl(sDc1))

c --- create directory ---
c --- delete entire folder if exists ---
      cmd = 'rm -rf ' // trim(folder)
      call system(cmd)

c --- recreate fresh folder ---
      cmd = 'mkdir -p ' // trim(folder)
      call system(cmd)
      
      write(*,*) 'FINAL FOLDER:', folder

c***************************************************************


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
      allocate(uh_L(NxC, Ny, Nz))

      
C***************Initiation of fields*************************
c Initial phi field
      c_x=(Nx*1.d0+1)/2
      c_y=(Ny*1.d0+1)/2
      c_z=(Nz*1.d0+1)/2
      
      
c For elipsoidal only
      a_x = 25.d0    ! semi-axis length along x
      a_y = 15.d0    ! semi-axis length along y
      a_z = 15.d0    ! semi-axis length along z
   
       do 13 i=1,Nx
        do 14 j=1,Ny
         do 15 k=1,Nz

c    Cubic Initial Field 
! cube half width (distance from center to face)
c           a = 10         ! use r0 as cube half-length or rename as needed

c           dx = abs(dble(i) - c_x) 
c           dy = abs(dble(j) - c_y) 
c           dz = abs(dble(k) - c_z) 
           
c           if ((dx.le.a).and.(dy.le.a).and.(dz.le.a)) then
c             u(i,j,k)=1.d0
c           else
c             u(i,j,k)=-1.d0
c           endif


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

c**************************************************************
 15     continue
 14    continue
 13   continue
 
 

      call compute_grad_u(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
        call compute_curvature_real(Nx,Ny,Nz,NxC,uh,dx,dy,dz,
     1 tol,norm,ux,uy,uz,iku,Curve)
      eps = 0.2   ! thickness of surface shell
      band = 0.5   ! thickness of belly region

      do 23 i=1,Nx
       do 24 j=1,Ny
        do 25 k=1,Nz

          dx = (dble(i) - c_x)/a_x
          dy = (dble(j) - c_y)/a_y
          dz = (dble(k) - c_z)/a_z

          ellipsoid = dx*dx + dy*dy + dz*dz

c --- Only belly surface
          if (abs(ellipsoid - 1.0).le.eps.and.abs(dx).le.band) then
              c(i,j,k) = c_ini
          else
              c(i,j,k) = 0.0
          endif

 25     continue
 24    continue
 23    continue
       
     

c********Saving files*******     
       

       call fft_backward(Nx, Ny, Nz, uh, u)
       u=u*norm
       
       fname = trim(folder)//'/step_0.dat'
       open(8,file=fname)     
       do 41 i=1,Nx
        do 52 j=1,Ny
         do 70 k=1, Nz
          write(8,*)i, j, k, u(i,j,k), c(i,j,k), Curve(i,j,k), 
     1      sqrt(iku(i,j,k))
 70       continue
 52      continue
 41     continue
       close(8)
      
       call fft_forward(Nx, Ny, Nz, c, ck)
       call fft_forward(Nx, Ny, Nz, u, uh)

c *******************************************************

c*****************Time Marching********************************

      do 2 it=1, NT_f

      call fft_backward(Nx, Ny, Nz, ck, c)
      call fft_backward(Nx, Ny, Nz, uh, u)
      
      c=c*norm
      u=u*norm
             
       do 119 i=1,Nx
        do 120 j=1,Ny
         do 122 k=1, Nz

          fc(i,j,k)=Dc*(2*As*(((u(i,j,k)**2)-1)**2)*((c(i,j,k)-c_max)
     1    *c(i,j,k)*(2*c(i,j,k)-c_max))+Af*2*(u(i,j,k)**2)*c(i,j,k))
             
           fu(i,j,k)=Du*(4*As*u(i,j,k)*(u(i,j,k)**2-1)*c(i,j,k)**2
     1       *(c(i,j,k)-c_max)**2+Af*2*u(i,j,k)*c(i,j,k)**2)
     
           cc(i,j,k)=Curve(i,j,k)*c(i,j,k) 



     
 122      continue    
 120     continue    
 119    continue
 
       call fft_forward(Nx, Ny, Nz, fc, fch)
       call fft_forward(Nx, Ny, Nz, c, ck)
       
       call fft_forward(Nx, Ny, Nz, u, uh)
       call fft_forward(Nx, Ny, Nz, fu, fuh)

       call fft_forward(Nx, Ny, Nz, cc, cch)
 
 
       do 126 i=1, NxC
         do 127 j=1, Ny
          do 128 k=1, Nz      	    
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)


      	    nuc=Dc1*lambda*k2*2*As
     	    
            ck(i,j,k)=exp(-k2*nuc*dt)*(ck(i,j,k)-k2*fch(i,j,k)*dt
     1        +lambda2*(-k2)*cch(i,j,k)*dt)
     
c            if (i.eq.8.and.j.eq.8.and.k.eq.8) 
c     1       write(*,*) it, k2*fch(i,j,k),-k2*lambda2*cch(i,j,k)
            
            

            nuu=Du*2*sig*k2
            uh(i,j,k)=(exp(-nuu*k2*dt))*(uh(i,j,k)-k2*fuh(i,j,k)*dt)
 128      continue
 127     continue
 126    continue
 

 
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
 
      call compute_grad_u(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
     
      call compute_curvature_real(Nx,Ny,Nz,NxC,uh,dx,dy,dz,
     1 tol,norm,ux,uy,uz,iku,Curve)
         
c***********************************************
       
     

c********Saving files*******

      if (mod(it,snap_int).eq.0) then      
      
       
       call fft_backward(Nx, Ny, Nz, ck, c)
       call fft_backward(Nx, Ny, Nz, uh, u)
       c=c*norm
       u=u*norm
       
      
      write(fname,'(A,I0,A)') 'step_', it, '.dat'
      fname = trim(folder)//'/'//trim(fname)
      open(8,file=fname)   
      
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
      

 2    continue


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

       do 80 i = 2, Nx-1
        do 81 j = 2, Ny-1
         do 82 k = 2, Nz-1
          mag=sqrt(iku(i,j,k)) 
            dnx_dx = (n_x(i+1,j,k) - n_x(i-1,j,k)) / (2.d0*dx)
            dny_dy = (n_y(i,j+1,k) - n_y(i,j-1,k)) / (2.d0*dy)
            dnz_dz = (n_z(i,j,k+1) - n_z(i,j,k-1)) / (2.d0*dz)
            kappa = -(dnx_dx + dny_dy + dnz_dz)
            Curve(i,j,k)=kappa/2.d0*(tanh(10*(u(i,j,k)+0.6))
     1                   +tanh(10*(-u(i,j,k)+0.6)))

 82      continue
 81     continue
 80    continue

c-------------------3 point formula with boundary correction-----------------------------
c       do 83 i = 3, Nx-2
c       	do 84 j = 3, Ny-2
c       	 do 85 k = 3, Nz-2
c       	   mag=sqrt(iku(i,j,k))
c       	   if (mag.le.tol) then
c       	    Curve(i,j,k)=0
c       	   else
c================X-direction========================
c            valid_p=((sqrt(iku(i+1,j,k)).ge.tol))
c            valid_m=((sqrt(iku(i-1,j,k)).ge.tol))
                
c            if (valid_p.and.valid_m) then
c             dnx_dx = (n_x(i+1,j,k) - n_x(i-1,j,k)) / (2.d0*dx)
c            else if (valid_p.and.sqrt(iku(i+2,j,k)).ge.tol) then
c             dnx_dx = (-3*n_x(i,j,k) + 4*n_x(i+1,j,k) - n_x(i+2,j,k))
c     1        / (2.d0*dx)
c            else if (valid_m.and.sqrt(iku(i-2,j,k)).ge.tol) then
c             dnx_dx = (3*n_x(i,j,k) - 4*n_x(i-1,j,k) + n_x(i-2,j,k))
c     1       / (2.d0*dx)
c            else
c             dnx_dx=0.d0
c            endif
c===============Y-direction=======================
c            valid_p=((sqrt(iku(i,j+1,k)).ge.tol))
c            valid_m=((sqrt(iku(i,j-1,k)).ge.tol))
c                 
c            if (valid_p.and.valid_m) then
c             dny_dy = (n_y(i,j+1,k) - n_y(i,j-1,k)) / (2.d0*dy)
c            else if (valid_p.and.sqrt(iku(i,j+2,k)).ge.tol) then
c             dny_dy = (-3*n_y(i,j,k) + 4*n_y(i,j+1,k) - n_y(i,j+2,k))
c     1        / (2.d0*dy)
c            else if (valid_m.and.sqrt(iku(i,j-2,k)).ge.tol) then
c             dny_dy = (3*n_y(i,j,k) - 4*n_y(i,j-1,k) + n_y(i,j-2,k))
c     1       / (2.d0*dy)
c            else
c             dny_dy=0.d0
c            endif                  
c===============Z-direction=======================
c            valid_p=((sqrt(iku(i,j,k+1)).ge.tol))
c            valid_m=((sqrt(iku(i,j,k-1)).ge.tol))
c                 
c            if (valid_p.and.valid_m) then
c             dnz_dz = (n_z(i,j,k+1) - n_z(i,j,k-1)) / (2.d0*dz)
c            else if (valid_p.and.sqrt(iku(i,j,k+2)).ge.tol) then
c             dnz_dz = (-3*n_z(i,j,k) + 4*n_z(i,j,k+1) - n_z(i,j,k+2))
c     1        / (2.d0*dz)
c            else if (valid_m.and.sqrt(iku(i,j,k-2)).ge.tol) then
c             dnz_dz = (3*n_z(i,j,k) - 4*n_z(i,j,k-1) + n_z(i,j,k-2))
c     1       / (2.d0*dz)
c            else
c            dnz_dz=0.d0
c            endif
c            kappa=-(dnx_dx+dny_dy+dnz_dz)
c            Curve(i,j,k)=kappa
c            
c          endif
          
c 85      continue
c 84     continue
c 83    continue
c-------------------------------------------------------------------------------------
       call fft_forward(Nx, Ny, Nz, u, uh)
       end subroutine compute_curvature_real

c**************************************************************

      subroutine strip_spaces(str)
      character*(*) str
      integer i, j, n

      n = len(str)
      j = 1

      do 100 i = 1, n
         if (str(i:i) .ne. ' ') then
            str(j:j) = str(i:i)
            j = j + 1
         endif
 100   continue

      do 200 i = j, n
         str(i:i) = ' '
 200   continue
       end subroutine strip_spaces
       
c*****************************************************************
