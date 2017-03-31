module Print_data               !handle the flow information, print it          
    use FVM
    implicit none
    
    real(8)::Cp(2:imax+1,2:jmax+1)   !the pressure coefficient of nodes
       
contains

subroutine  Cp_cal                  !calculate the flow information of nodes from cell centres
    implicit none
    integer::i,j
    real(8)::Cp_center(0:imax+1,0:jmax+1) !the pressure coefficient of cell centres
    
    Cp_center=0.0
    Cp_center(2:imax,2:jmax)=(p(2:imax,2:jmax)/p_inf-1.0)*2.0/(Ma_inf**2*gamma)

    do j=2,jmax+1
    do i=2,imax+1
       Cp(i,j)=Cp_center(i-1,j-1)*ds(1,i,j)+Cp_center(i,j-1)*ds(2,i,j)+Cp_center(i,j)*ds(3,i,j)+Cp_center(i-1,j)*ds(4,i,j)
    end do
    end do
    
    
end subroutine

subroutine Output
    implicit none
    integer::i,j
    
    integer,parameter::fileid=7
    !write(*,*)   "Output"
    !**************************************
    
    open(fileid,file="xy.dat")
    write(fileid,99) ( (i,j, x(i,j),y(i,j),i=2,imax+1 ),j=2,jmax+1) 
99  format(I3,1X,I3,3X,F15.8,3X,F15.8)
    
    open(fileid,file="Cp.dat")
    write(fileid,*)  ' VARIABLES="X","Y","Cp" '
    write(fileid,"(A7,I3,1X,A4,I3)")  "zone I=",imax,", J=",jmax
    write(fileid,100) ( ( x(i,j),y(i,j),Cp(i,j),i=2,imax+1 ),j=2,jmax+1) 
    close(fileid)
    
    open(fileid,file="Cp-b.dat")
    write(fileid,*)  'VARIABLES="X","Y","Cp"'
    write(fileid,"(A7,I3,1X,A4,I3)")  "zone I=",imax,", J=",jmax
    write(fileid,100) (  x(i,2),y(i,2),Cp(i,2),i=2,imax+1 )
100  format(3X,F15.8,3X,F15.8,3X,F15.8)
    close(fileid)
    
    open(fileid,file="Flow_info.dat")
    write(fileid,110) ( ( i,j,rou(i,j),p(i,j),u(i,j),v(i,j),i=2,imax ),j=2,jmax) 
110 format(I3,1X,I3,3X,F15.8,3X,F15.8,3X,F15.8,3X,F15.8)
    
    close(fileid)
    
end subroutine

end module 