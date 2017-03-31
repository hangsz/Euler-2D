subroutine Cp_cal(filename)
    implicit none
    integer::num_b
    integer,allocatable::cell_b(:)
    real(8),allocatable::Cp_b(:,:),Cp_u(:,:),Cp_l(:,:)
    real(8)::Cn = 0.0
    real(8)::temp(2)
    integer,parameter::fileid=7
    integer::node1,node2,node3
    integer::i,j,k
    logical::alive
    character(len =30)::filename,filename1
    
    inquire(file="Cell_b.dat",exist = alive )
    
    if( alive )  then
        
        open(fileid,file = "Cell_b.dat")
            
        read(fileid,*)  num_b
        
        allocate( cell_b(num_b) )  
        
        do i=1,num_b
                read(fileid,*)  cell_b(i)
         end do
    
        close(fileid)    
    else
        num_b=0.0
        do i =1,nedges
            if(iedge(4,i) == -1) then 
                num_b = num_b + 1    
            end if  
        end do
    
        allocate( cell_b(num_b) )  
     
        j=1
        do i =1,nedges 
            if( iedge(4,i) .EQ. -1 ) then 
                cell_b(j) = iedge(3,i) 
                !write(*,*) j,cell_b(j) 
                j=j+1
            end if
        end do
        
        
        !output the numberof cells  on the wall
        open(fileid,file = "Cell_b.dat")
            
        write(fileid,*)  num_b
        
         do i=1,num_b
                write(fileid,*)  cell_b(i)
         end do
    
        close(fileid)
        
    end if
    
    
    allocate( Cp_b(3,num_b) ) 
    allocate( Cp_u(2,num_b/2) ) 
    allocate( Cp_l(2,num_b/2) )
    
    !the lateral coordinates and pressure coefficients  of cp
    !do i=1,num_b
     !  write(*,*)  num_b
        !node1=icell( 1,cell_b(i) ) 
        !node2=icell( 2,cell_b(i) ) 
        !node3=icell( 3,cell_b(i) ) 
        !write(*,*) node1,node2,node3
        !Cp_b(1,i)=( xy(1, node1) + xy(1, node2 ) + xy (1, node3)  ) /3.0
        Cp_b(1,:)=( xy(1, icell(1,cell_b) ) + xy(1, icell(2,cell_b) ) + xy (1, icell(3,cell_b) )  ) /3.0 
        Cp_b(2,:)=( xy(2, icell(1,cell_b) ) + xy(2, icell(2,cell_b) ) + xy (2, icell(3,cell_b) )  ) /3.0 
        
        Cp_b(3,:)=( U(5, cell_b )-p_inf ) / ( 1.0/2*rou_inf*(u_inf**2+v_inf**2) )
      
     
    where( Cp_b(2,:)>1.5)  
        Cp_b(2,:) = 1.5
    else where( Cp_b(2,:)<-1.5)
        Cp_b(2,:) = -1.5
    end where
    !end do
    
    j=1
    k=1 
   
    ! divide the Cp into Cp_up_boundary and Cp_lown_boundary
    do i=1,num_b
        if( Cp_b(2,i) .GT. 0.0 ) then
          
            Cp_u(1,j)= Cp_b(1,i)
            Cp_u(2,j)= Cp_b(3,i)
            j=j+1
        else
           
            Cp_l(1,k)= Cp_b(1,i)
            Cp_l(2,k)= Cp_b(3,i)
            k=k+1
        end if
        
    end do
    !sort 
    do i=1,num_b/2
    do j=num_b/2-1,i,-1
        if( Cp_u(1,j) .GT. Cp_u(1,j+1) ) then 
            temp=Cp_u(:,j+1)
            Cp_u(:,j+1)=Cp_u(:,j)
            Cp_u(:,j)=temp
        end if
        
        if( Cp_l(1,j).GT. Cp_l(1,j+1)) then 
            temp=Cp_l(:,j+1)
            Cp_l(:,j+1)=Cp_l(:,j)
            Cp_l(:,j)=temp
        end if
    end do
    end do
    
    
    
    i=1
    Cn = - Cp_u(1,1)*Cp_u(2,1)
    do i=2,num_b/2
        Cn = Cn- Cp_u(2,i)*( Cp_u(1,i)-Cp_u(1,i-1))
    end do
    
    i=1
    Cn = Cn + Cp_l(1,1)*Cp_l(2,1)
    do i=2,num_b/2
        Cn = Cn + Cp_l(2,i)*( Cp_l(1,i)-Cp_l(1,i-1))
    end do
    
    write(*,*)  Cn
    
        
    filename1 = trim(filename)//"-Cp_u.dat"
    open(fileid,file=filename1)
    do i=1,num_b/2 
        write(fileid,*)  Cp_u(1,i),Cp_u(2,i)
    end do
    
    close(fileid)
    
    filename1 = trim(filename) //"-Cp_l.dat"
    open(fileid,file=filename1)
    do i=1,num_b/2
        write(fileid,*) Cp_l(1,i),Cp_l(2,i)
    end do
    close(fileid)
    
    
end subroutine

            
    
        
    
    
    
    
    
    
    