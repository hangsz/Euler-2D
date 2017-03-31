!pay attention: add two dummy grids ,notice the index range
module Grid_handle   !grid handle module  
    use Control_para
    implicit none      
    
    real(8)::x(2:imax+1,2:jmax+1),y(2:imax+1,2:jmax+1)                       !the coordinates of every node
    real(8),dimension(2,2:imax+1,2:jmax+1)::vector1,vector2                  !the normal vector of every edge in i,j direction
    real(8),dimension(2:imax+1,2:jmax+1)::ds1,ds2                           !the length of every edge in i,j direction
    
    real(8)::ds(4,2:imax+1,2:jmax+1)                                         !the four distances of every node to its adjacent cell centres
    real(8)::vol(0:imax+2,0:jmax+2)                                          !the volume of every cell 
contains  

subroutine Grid
    implicit none
    call Read_grid                  !读取网格
    call Geometrical_para           !进行有限体积法中需要的向量、体积等几何参数计算
end subroutine

subroutine  Read_grid
    implicit none
    integer,parameter::fileid=10 
    character(len=20)::firstline
    character(len=20)::secondline
    character(len=80)::filename   
    integer::error                
    logical::alive                
    integer::i,j
    
    !write(*,*)  "Read_grid"
    
   
    !write(*,*)  "Filename:"         
    !read(*,"(A80)") filename      
    !filename="XY-C-Y.dat"
    filename="naca0012.grd"
    
    inquire(file=filename,exist=alive)    
    if(.not.alive) then
        write(*,*) trim(filename),"Does't exist!"   
        stop
    end if
    
    open(fileid,file=filename,action='read',status='old')
    read(fileid,*,iostat=error,err=10)  firstline
    read(fileid,*,iostat=error,err=10)  secondline
    
    !read(fileid,"(2X,F11.8,2X,F11.8)",iostat=error,err=10)  ((x(i,j),y(i,j),i=2,imax+1),j=2,jmax+1)
    read(fileid,"(1X,E16.9E3,1X,E16.3E3)",iostat=error,err=10)  ((x(i,j),y(i,j),i=2,imax+1),j=2,jmax+1)
10  if(error>0)  then
        write(*,*) "...Read error..."
        stop 
    end if   
    close(fileid)
end subroutine

subroutine  Geometrical_para
    implicit none
    integer::i,j
    real(8)::center(2,2:imax,2:jmax)
    write(*,*)  "Geometrical_para"
    !************************************************************************
    !calculate the coordinates of every cell centre
    do j=2,jmax
    do i=2,imax
       center(1,i,j)=0.5*(  0.5*(x(i+1,j+1)+x(i+1,j))+0.5* (x(i,j+1)+x(i,j))  )
       center(2,i,j)=0.5*(  0.5*(y(i+1,j+1)+y(i+1,j))+0.5* (y(i,j+1)+y(i,j))  )
    end do
    end do 
    !****************************
    !calculate the distance of every node to the four adjacent cell centres (for left bottom in counterclockwise direction）
    !inner points
    do j=3,jmax
    do i=3,imax
       ds(1,i,j)=sqrt( (x(i,j)-center(1,i-1,j-1))**2+ ( y(i,j)-center(2,i-1,j-1))**2 )
       ds(2,i,j)=sqrt( (x(i,j)-center(1,i,j-1))**2  + ( y(i,j)-center(2,i,j-1))**2 )
       ds(3,i,j)=sqrt( (x(i,j)-center(1,i,j))**2    + ( y(i,j)-center(2,i,j))**2 )
       ds(4,i,j)=sqrt( (x(i,j)-center(1,i-1,j))**2  + ( y(i,j)-center(2,i-1,j))**2 )
    end do
    end do
   
    !boundaries,excluding vertex
    !j=2
    do i=3,imax
        ds(3,i,2)=sqrt( (x(i,2) -center(1,i,2))**2  + (y(i,2)  -center(2,i,2))**2 )
        ds(4,i,2)=sqrt( (x(i,2) -center(1,i-1,2))**2+ (y(i,2)  -center(2,i-1,2))**2 )
    end do
    !j=jmax+1
    do i=3,imax
        ds(1,i,jmax+1)=sqrt( (x(i,jmax+1)-center(1,i-1,jmax))**2+ (y(i,jmax+1)  -center(2,i-1,jmax))**2 )
        ds(2,i,jmax+1)=sqrt( (x(i,jmax+1)-center(1,i,jmax))**2  + (y(i,jmax+1)  -center(2,i,jmax))**2 )
    end do
    
    !i=2
    do j=3,jmax
        ds(2,2,j)=sqrt( (x(2,j)  -center(1,2,j-1))**2+ (y(2,j)  -center(2,2,j-1))**2 )
        ds(3,2,j)=sqrt( (x(2,j)  -center(1,2,j))**2  + (y(2,j)  -center(2,2,j))**2 )
    end do
    !i=imax+1
    do j=3,jmax
        ds(1,imax+1,j)=sqrt( (x(imax+1,j)  -center(1,imax,j-1))**2+ (y(imax+1,j)  -center(2,imax,j-1))**2 )
        ds(4,imax+1,j)=sqrt( (x(imax+1,j)  -center(1,imax,j))**2  + (y(imax+1,j)  -center(2,imax,j))**2 )
    end do
     !four vertex 
    ds(3,2,2)=sqrt( (x(2,2)  -center(1,2,2))**2  + (y(2,2)  -center(2,2,2))**2 )
    ds(4,imax+1,2)=sqrt( (x(imax+1,2)  -center(1,imax,2))**2  + (y(imax+1,2)  -center(2,imax,2))**2 )
    ds(1,imax+1,jmax+1)=sqrt( (x(imax+1,jmax+1)-center(1,imax,jmax))**2+ (y(imax+1,jmax+1)  -center(2,imax,jmax))**2 )
    ds(2,2,jmax+1)=sqrt( (x(2,jmax+1)-center(1,2,jmax))**2  + (y(2,jmax+1)  -center(2,2,jmax))**2 )
  
    do j=2,jmax+1
    do i=2,imax+1
       ds(:,i,j)=ds(:,i,j)/sum(ds(:,i,j) )
    end do
    end do
    !********************************************************************
    !the volume of every volume
    do j=2,jmax
    do i=2,imax
       vol(i,j)=0.5*( (x(i,j)-x(i+1,j+1))*(y(i+1,j)-y(i,j+1))+ (x(i,j+1)-x(i+1,j))*(y(i,j)-y(i+1,j+1))  )
    end do
    end do
    !***************************************************************
    !along i increase direction,the normal vector is clockwise
    do j=2,jmax+1
    do i=2,imax
       vector1(1,i,j)=  y(i+1,j) - y(i,j)
       vector1(2,i,j)=- x(i+1,j) + x(i,j)
       ds1(i,j)=sqrt(vector1(1,i,j)**2+vector1(2,i,j)**2)
    end do
    end do
    !along j increase direction,the normal vector is clockwise
    do i=2,imax+1 
    do j=2,jmax
       vector2(1,i,j)= y(i,j+1) -y (i,j)
       vector2(2,i,j)=-x(i,j+1) + x(i,j)
       ds2(i,j)=sqrt( vector2(1,i,j)**2+vector2(2,i,j)**2 )
    end do
    end do
    
end subroutine
end module