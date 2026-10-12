c Replicating Hernandez 1 completely with curvature..... real space curvature calculated , time increment..... rotated ellipsoid eqn... surface tension dependent on concentration^2.... conserved noise.... Sriram curvature term...input file changed... Adams bashforth

      implicit double precision(a-h, o-z)
      include 'fftw3.f'


      double precision dx, dy, dz, k2, norm, kv, nuc
      double precision bx,by,bz, c_x, c_y,c_z, r, h
      double precision nuu, r_ellip,  F, a
      include 'params.inc'

      double precision, allocatable :: u(:,:,:), c(:,:,:), fc(:,:,:)
      double precision, allocatable :: kx(:), ky(:), kz(:),ux(:,:,:)
      double precision, allocatable :: uy(:,:,:), uz(:,:,:), ub(:,:,:)
      double precision, allocatable :: iku(:,:,:), fu(:,:,:),k2u(:,:,:)
      double precision, allocatable :: Curve(:,:,:), n_x(:,:,:)
      double precision, allocatable :: n_y(:,:,:), n_z(:,:,:)
      double precision, allocatable :: cc(:,:,:), ciku(:,:,:)
      double precision, allocatable :: k2ub(:,:,:), zeta_x(:,:,:)
      double precision, allocatable :: zeta_y(:,:,:), zeta_z(:,:,:)
      
      double precision, allocatable :: cx(:,:,:), cy(:,:,:), cz(:,:,:)
 
      complex*16, allocatable :: ck(:,:,:), fch(:,:,:), ikuhx(:,:,:)
      complex*16, allocatable :: ikuhy(:,:,:), ikuhz(:,:,:), uh(:,:,:)
      complex*16, allocatable :: fuh(:,:,:), Curveh(:,:,:)
      complex*16, allocatable :: cch(:,:,:), ikbcu(:,:,:), cikuh(:,:,:)
      complex*16, allocatable :: k2uh(:,:,:), k2ubh(:,:,:), ubh(:,:,:)
      complex*16, allocatable :: zeta_xh(:,:,:), zeta_yh(:,:,:)
      complex*16, allocatable :: zeta_zh(:,:,:), Fu_o(:,:,:)
      complex*16, allocatable :: Fu_n(:,:,:), Fc_o(:,:,:), Fc_n(:,:,:)
      

      integer*8 :: plan
      integer :: NxC, it, isnap, l
      complex*16, parameter :: iota = (0.d0, 1.d0)
      complex*16, parameter :: minus_iota = (0.d0, -1.d0)
      character(len=1024) :: fname
      character*10 run_name
      logical :: file_exists

      call get_command_argument(1, run_name)

c      call execute_command_line('rm -f step_*.dat')
c      call execute_command_line('rm -f step_*.vtu')
c      call execute_command_line('rm -f cutc_step_*.dat')
c      call execute_command_line('rm -f cutc_step_*.vtu')
c      call execute_command_line('rm -f cute_step_*.dat')
c     call execute_command_line('rm -f cute_step_*.vtu')  

c************Folder Creation*******************************
c --- master parameter file (CSV) ---
c --- CSV master parameter file ---
      fname = '../params_master_file.csv'

c --- check if file already exists ---
      inquire(file=fname,exist=file_exists)

c --- open file in append mode ---
      open(99,file=fname,position='append')

c --- write header only if file does not exist ---
      if (.not. file_exists) then
         write(99,'(A)') 
     1 'File name,Nx,axis1,axis2,axis3,angle,lambda,As,Af,Dc,Du,'//
     2 'sigma,lambda1,alpha,beta,Ab,eps,cmax,c_ini,kH'
      endif
      
c --- write parameter row ---
       write(99,'(A,",",I0,17(",",F12.6))') 
     1 trim(run_name), Nx, a_x, a_y, a_z, theta_rot, lambda,
     2 As, Af, Dc, Du, sig, lambda1,
     3 alpha, beta, Ab, eps, c_max, c_ini,kH

      close(99)


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
      allocate(k2u(Nx,Ny,Nz))
      allocate(k2ub(Nx,Ny,Nz))
      allocate(ub(Nx,Ny,Nz))
      allocate(cx(Nx,Ny,Nz))
      allocate(cy(Nx,Ny,Nz))
      allocate(cz(Nx,Ny,Nz))
      allocate(zeta_x(Nx,Ny,Nz))
      allocate(zeta_y(Nx,Ny,Nz))
      allocate(zeta_z(Nx,Ny,Nz))
              
      allocate(ck(NxC, Ny, Nz))
      allocate(fch(NxC, Ny, Nz))
      allocate(uh(NxC, Ny, Nz))
      allocate(fuh(NxC, Ny, Nz))
      allocate(Curveh(NxC, Ny, Nz))
      allocate(ikbcu(NxC, Ny, Nz))
      allocate(cch(NxC, Ny, Nz))
      allocate(cikuh(NxC, Ny, Nz))
      allocate(k2uh(NxC, Ny, Nz))
      allocate(k2ubh(NxC,Ny,Nz))
      allocate(ubh(NxC,Ny,Nz))
      allocate(zeta_xh(NxC,Ny,Nz))
      allocate(zeta_yh(NxC,Ny,Nz))
      allocate(zeta_zh(NxC,Ny,Nz))
      allocate(Fu_o(Nxc, Ny, Nz))
      allocate(Fu_n(Nxc, Ny, Nz))
      allocate(Fc_o(Nxc, Ny, Nz))
      allocate(Fc_n(Nxc, Ny, Nz))

      fname ='step_0.dat'
      
      open(unit=8, file=fname, status='old', action='read')
      do 41 i=1, Nx
      	do 52 j=1, Ny
      	  do 70 k=1, Nz
      		read(8, *) ii, jj, kk, u(i,j,k), c(i,j,k), 
     1        Curve(i,j,k), dummy
 70       continue
 52     continue
 41   continue
      
      close(8)
      
       call fft_forward(Nx, Ny, Nz, c, ck)
       call fft_forward(Nx, Ny, Nz, u, uh)

c *******************************************************

c*****************Time Marching********************************

      do 2 it=1, NT_f

      do 6 i=1, NxC
      	do 7 j=1, Ny
      	  do 8 k=1,Nz
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)
      	    k2uh(i,j,k)=-k2*uh(i,j,k)      	  
 8        continue
 7       continue
 6     continue
       
      call fft_backward(Nx, Ny, Nz, ck, c)
      call fft_backward(Nx, Ny, Nz, uh, u)
      call fft_backward(Nx, Ny, Nz, k2uh, k2u)
      
      c=c*norm
      u=u*norm
      k2u=k2u*norm
      
      do 9 i=1,Nx
      	do 10 j=1,Ny
      	  do 110 k=1,Nz
      	    ub(i,j,k)=-u(i,j,k)+u(i,j,k)**3-eps**2*k2u(i,j,k)
     1 	    +beta*c(i,j,k)**2*eps*(1-(u(i,j,k)**2))
 110        continue
 10     continue
 9     continue
 
      call fft_forward(Nx, Ny, Nz, ub, ubh)
 
      do 16 i=1, NxC
      	do 17 j=1, Ny
      	  do 18 k=1,Nz
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)
      	    k2ubh(i,j,k)=-k2*ubh(i,j,k) 
 18       continue
 17     continue
 16    continue
 
       call fft_backward(Nx, Ny, Nz, k2ubh, k2ub)
       k2ub=k2ub*norm
       
       call compute_grad_u2(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
     1  uy, uz, iku, norm)
     
       call fft_backward(Nx, Ny, Nz, uh, u)
       u=u*norm
       
       call compute_grad_uh(Nx,Ny,Nz,NxC,u,uh,kx,ky,kz
     1       ,ux,uy,uz,norm)
     
       call compute_grad_uh(Nx,Ny,Nz,NxC,c,ck,kx,ky,kz
     1       ,cx,cy,cz,norm)
     

       
       do 119 i=1,Nx
        do 120 j=1,Ny
         do 122 k=1, Nz
          
         fc(i,j,k)=Dc*(Ab*4*eps*beta*c(i,j,k)*(1-u(i,j,k)**2)*ub(i,j,k)
     1    +2*As*(((u(i,j,k)**2)-1)**2)*((c(i,j,k)-c_max)
     1    *c(i,j,k)*(2*c(i,j,k)-c_max))+Af*2*(u(i,j,k)**2)*c(i,j,k)
c     1     +2*sig*c(i,j,k)*iku(i,j,k)
     1      -kH*k2*u(i,j,k))
             
           fu(i,j,k)=Du*(Ab*(2*ub(i,j,k)*(-1+3*u(i,j,k)**2
     1      -2*eps*beta*c(i,j,k)**2*u(i,j,k))-2*eps**2*k2ub(i,j,k))
     1      +4*As*u(i,j,k)*(u(i,j,k)**2-1)*c(i,j,k)**2
     1       *(c(i,j,k)-c_max)**2+Af*2*u(i,j,k)*c(i,j,k)**2
c     1        -2*sig*(2*c(i,j,k)*(ux(i,j,k)*cx(i,j,k)+uy(i,j,k)
c     1        *cy(i,j,k)+uz(i,j,k)*cz(i,j,k))+c(i,j,k)**2*k2u(i,j,k))
     1         -kH*k2*c(i,j,k))
     
           cc(i,j,k)=Curve(i,j,k)*c(i,j,k) 
           
c noise           call box_muller(zeta_xi)
c noise           call box_muller(zeta_yj)
c noise           call box_muller(zeta_zk)
           
c noise           zeta_x(i,j,k)=zeta_xi
c noise           zeta_y(i,j,k)=zeta_yj
c noise           zeta_z(i,j,k)=zeta_zk

     
 122      continue    
 120     continue    
 119    continue
 
       call fft_forward(Nx, Ny, Nz, fc, fch)
       call fft_forward(Nx, Ny, Nz, c, ck)
       
       call fft_forward(Nx, Ny, Nz, u, uh)
       call fft_forward(Nx, Ny, Nz, fu, fuh)

       call fft_forward(Nx, Ny, Nz, cc, cch)
 
c noise       call fft_forward(Nx, Ny, Nz, zeta_x, zeta_xh)
c noise       call fft_forward(Nx, Ny, Nz, zeta_y, zeta_yh)
c noise       call fft_forward(Nx, Ny, Nz, zeta_z, zeta_zh)
       
       
       do 126 i=1, NxC
         do 127 j=1, Ny
          do 128 k=1, Nz      	    
      	    k2=(kx(i)**2+ky(j)**2+kz(k)**2)


      	    nuc=Dc*lambda*k2*2*As

c noise      	    eta=iota*str*(kx(i)*zeta_xh(i,j,k)+ky(j)*zeta_yh(i,j,k)
c noise     1 	       +kz(k)*zeta_zh(i,j,k))
     
     	    
     	    Fc_n(i,j,k)=-k2*fch(i,j,k)+lambda1*(-k2)*cch(i,j,k)       
     
            Fu_n(i,j,k)=-k2*fuh(i,j,k)
         
            nuu=Du*2*sig*k2
            
            if (it.eq.1) then
                ck(i,j,k) = exp(-k2*nuc*dt)*(ck(i,j,k)+dt*Fc_n(i,j,k))
                uh(i,j,k) = exp(-nuu*k2*dt)*(uh(i,j,k)+dt*Fu_n(i,j,k))

                Fc_o(i,j,k) = Fc_n(i,j,k)
                Fu_o(i,j,k) = Fu_n(i,j,k)
                
           else
           
            ck(i,j,k)=exp(-k2*nuc*dt)*(ck(i,j,k)+(1.5d0*Fc_n(i,j,k)
     1       -0.5d0*Fc_o(i,j,k))*dt)
                 
            uh(i,j,k)=(exp(-nuu*k2*dt))*(uh(i,j,k)+(1.5d0*Fu_n(i,j,k)
     1       -0.5d0*Fu_o(i,j,k))*dt)
     
            Fc_o(i,j,k)=Fc_n(i,j,k)            
            Fu_o(i,j,k)=Fu_n(i,j,k)
            
            endif
         
 128      continue
 127     continue
 126    continue
 

 
      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
 
      call compute_grad_u2(Nx, Ny, Nz, NxC, u, uh, kx, ky, kz, ux,
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
      subroutine compute_grad_u2(Nx,Ny,Nz,NxC,u,uh,kx,ky,kz
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
        do 520 j = 1, Ny
            do 53 k = 1, Nz
                uhx(i,j,k) = iota * kx(i) * uh(i,j,k)
                uhy(i,j,k) = iota * ky(j) * uh(i,j,k)
                uhz(i,j,k) = iota * kz(k) * uh(i,j,k)
 53         continue
 520     continue
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

      end subroutine compute_grad_u2






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

            Curve(i,j,k)=kappa/2.d0*(tanh(10.d0*(u(i,j,k)+0.6d0))
     1                   +tanh(10.d0*(-u(i,j,k)+0.6d0)))

 82      continue
 81     continue
 80    continue


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

c*****************************************************************      
      subroutine compute_grad_uh(Nx,Ny,Nz,NxC,u,uh,kx,ky,kz
     1       ,ux,uy,uz,norm)
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


      do 511 i = 1, NxC
        do 521 j = 1, Ny
            do 531 k = 1, Nz
                uhx(i,j,k) = iota * kx(i) * uh(i,j,k)
                uhy(i,j,k) = iota * ky(j) * uh(i,j,k)
                uhz(i,j,k) = iota * kz(k) * uh(i,j,k)
 531         continue
 521     continue
 511    continue 

      call fft_backward(Nx, Ny, Nz, uh, u)
      u=u*norm
      
      call fft_backward(Nx, Ny, Nz, uhx, ux)
      call fft_backward(Nx, Ny, Nz, uhy, uy)
      call fft_backward(Nx, Ny, Nz, uhz, uz)
      ux=ux*norm
      uy=uy*norm
      uz=uz*norm

      end subroutine compute_grad_uh
      
c*****************************************************************

      
      subroutine box_muller(z)
      implicit none

      real*8, intent(out) :: z
      real :: u1, u2
      real, parameter :: pi = ACOS(-1.d0)

      call random_number(u1)
      call random_number(u2)

    ! Prevent log(0)
      u1 = max(u1, tiny(1.0))

      z = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)

      end subroutine box_muller
      
c****************************************************************

