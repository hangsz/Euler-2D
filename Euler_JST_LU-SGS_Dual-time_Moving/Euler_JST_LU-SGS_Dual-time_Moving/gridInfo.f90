module gridInfo
    !--------------------------------------------------------------
        !purpose: read the grid infomation and calculate the values which are necessary for flow caluclation
    !--------------------------------------------------------------
        
    implicit none
    
    integer::nnodes,ncells,nedges 
    integer,allocatable:: iedge(:,:),icell(:,:)
    real(8),allocatable:: xy(:,:),vol(:)
    
    real(8),allocatable:: vector(:,:),ds(:)             
    real(8),allocatable:: cell_center(:,:)
    real(8),allocatable:: rL(:,:),rR(:,:)
    real(8),allocatable:: rij(:,:),lij(:),tij(:,:)
    
    real(8),allocatable::d(:)
    
    !---------------------------------------------------------------
    
               !variable specification  
               
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
    !vector     :   the normal vector of each edge
    !ds         :   the length of each edge
    !cell_center:   the (x,y) coordinates of each cell
    !rL,rR      :   the vectors points from the edge center to the left and right cell center respectively
    !rij,lij,tij:   the vector points from left cell center to the right; the magnitude of the vector; the unitized vector
    !d          :   the nearest distance of each cell to the wall  
    
    !--------------------------------------------------------------------
    
contains


subroutine grid

    call readGrid
    call gridData
    
end subroutine

subroutine  readGrid
    !----------------------------------------------
         ! purpose : read the grid information
    !----------------------------------------------
         
    implicit none
	integer::i	
	
    integer::idGrid
    character(len = 30)::filename
    logical::alive
    
    !---------------------------------------------
    
    filename = 'grid/naca0012.grd'
    
    inquire(file =filename, exist =alive)
    if ( .NOT. alive)  then
        write(*,*)  filename,'not exist.'
        stop
    end if
    
    open(newunit=idGrid,file=filename)
            
        read(idGrid,*)  nnodes,nedges,ncells
            
        ! allocate memory
        allocate(iedge(4,nedges))
        allocate(xy(2,nnodes))
        allocate(icell(3,ncells))
        allocate(vol(ncells))
            
        do i=1,nnodes
            read(idGrid,*)   xy(1,i),xy(2,i)
        end do
            
        do i=1,nedges
            read(idGrid,*) iedge(1,i),iedge(2,i),iedge(3,i),iedge(4,i)
        end do
        
        do i=1,ncells
            read(idGrid,*) icell(1,i),icell(2,i),icell(3,i)
        end do
            
        do i=1,ncells
            read(idGrid,*) vol(i)
        end do
    close(idGrid)

end subroutine
 
subroutine  gridData
    !------------------------------------------------
        ! purpose: calculate related geometric variables
    !------------------------------------------------
         
    implicit none
    integer::i
     
    integer::start_node,end_node
    integer::ncl,ncr
     
    !------------------------------------------------
     
    allocate(vector(2,nedges))              
    allocate(ds(nedges))                    
    allocate(cell_center(2,ncells))
    allocate(rL(2,nedges))
    allocate(rR(2,nedges))
    allocate(rij(2,nedges))
    allocate(lij(nedges))
    allocate(tij(2,nedges))
     
     
    ! the normal vector and the length of every edge
    vector(1,:)= xy(2,iedge(2,:)) - xy(2,iedge(1,:))
    vector(2,:)=-xy(1,iedge(2,:)) + xy(1,iedge(1,:))
    ds(:)=sqrt( vector(1,:)**2+vector(2,:)**2 )
     
    cell_center = ( xy(:,icell(1,:) ) + xy(:,icell(2,:) ) + xy(:,icell(3,:) ) ) /3.0
     
end subroutine

 
end module