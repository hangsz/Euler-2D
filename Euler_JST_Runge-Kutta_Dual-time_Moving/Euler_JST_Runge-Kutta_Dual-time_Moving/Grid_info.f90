module Grid_info
    
    implicit none
    integer::nnodes,ncells,nedges 
    integer,allocatable:: iedge(:,:),icell(:,:)
    real(8),allocatable:: xy(:,:),vol(:)
    real(8),allocatable:: vector(:,:),ds(:)  
    real(8),allocatable:: cell_center(:,:)
    
contains

subroutine Grid
    call Read_grid
    call Grid_data
end subroutine

subroutine  Read_grid
    implicit none
	integer::i	
		!nnodes     :   node numbers
		!ncells     :   cell numbers
		!nedges     :   edge numbers
		!iedge      :   matrix for edge
		                !iedge(1,i) : edge's start point 'a'
		                !iedge(2,i) : edge's end point 'b'
		                !iedge(3,i) : edge's left cell  (all positive)
		                !iedge(4,i) : positive for right cell
		                            !iedge(4,i) -1 for wall boundary
		                            !iedge(4,i) -2 for farfiled boundary
		!xy         :   cartesian coordinates  xy(1,i) for x  xy(2,i) for y
		!icell      :   triangle cell's three nodes index
		!vol        :   cell's volume(area in 2d)

    open(10,file='grid/naca0012.grd')
            
        read(10,*)  nnodes,nedges,ncells
            
        !allocate memory
        allocate(iedge(4,nedges))
        allocate(xy(2,nnodes))
        allocate(icell(3,ncells))
        allocate(vol(ncells))
            
        do i=1,nnodes
            read(10,*)   xy(1,i),xy(2,i)
        end do
            
        do i=1,nedges
            read(10,*) iedge(1,i),iedge(2,i),iedge(3,i),iedge(4,i)
        end do
        
        do i=1,ncells
            read(10,*) icell(1,i),icell(2,i),icell(3,i)
        end do
            
        do i=1,ncells
            read(10,*) vol(i)
        end do

 end subroutine
 
 subroutine  Grid_data
     implicit none
     integer::i
     integer::node(3)                        !store the three nodes of every cell
    
     !-------------------------------------------------
     !allocate memory
     allocate(vector(2,nedges))               !the normal vector  of every edge
     allocate(ds(nedges))                     !the  length of every edge
     allocate(cell_center(2,ncells))
     
     !the normal vector and the length of every edge
     vector(1,:)= xy(2,iedge(2,:)) - xy(2,iedge(1,:))
     vector(2,:)=-xy(1,iedge(2,:)) + xy(1,iedge(1,:))
     ds(:)=sqrt( vector(1,:)**2+vector(2,:)**2 )
     !----------------------------------------------------------------
     cell_center = ( xy(:,icell(1,:) ) + xy(:,icell(2,:) ) + xy(:,icell(3,:) ) ) /3.0
   
     
 end subroutine
 
 end module