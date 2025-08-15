      implicit double precision(a-h, o-z)
      parameter(N=64)
      
      real, dimension(N,N,N) :: A
      
      c_x=(N*1.d0+1)/2
      c_y=(N*1.d0+1)/2
      c_z=(N*1.d0+1)/2
      rad_cell=24.d0
      layer=1.d0

      do 3 i=1,N
         do 4 j=1,N
            do 5 k=1,N
              xp=i*1.d0
              yp=j*1.d0
              zp=k*1.d0
              grid_d=grid_distance(xp,yp,zp,c_x,c_y,c_z)
              inner_b=(rad_cell-layer)
              outer_b=(rad_cell+layer)
              if (grid_d.le.inner_b) then
                 A(i,j,k)=1
              else if (grid_d.gt.inner_b.and.grid_d.le.outer_b) then
                 grid_c=grid_d-(rad_cell-layer)
                 A(i,j,k)=1-grid_c/(2*layer)
              else
                 A(i,j,k)=0
              endif
 5          continue
 4          continue
 3          continue









      open(100, file= 'Heat_Map_Cell.dat', status='replace')
c    Writing a matrix output to generate a heat map    
      do 1 i = 1, N
         do 2 j = 1, N
           do 6 k= 1,N
            write(100,'(I4,1x,I4,1x,I4,1x,F12.6)') i,j,k,A(i,j,k)
 6         continue
         write(100,*)  ! blank line between rows
 2       continue
       write(100,*)  ! blank line between rows   
 1     continue
      
      stop
      end
      
      
      double precision function grid_distance(x,y,z,xc,yc,zc)
      implicit double precision(a-h, o-z)
      grid_distance=sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
      return
      end
