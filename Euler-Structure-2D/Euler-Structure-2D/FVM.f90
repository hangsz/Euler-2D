module FVM          !Finite volume method
    use Grid_handle
    implicit none
    real(8)::u_inf,v_inf,a_inf      
    real(8),dimension(0:imax+2,0:jmax+2)::rou,p,u,v,Ma
    real(8)::W(0:imax+2,0:jmax+2,5) 
  
contains
!include "Boundary_cond.f90"
include "Boundary_cond - 1.f90"
include "Convective_flux.f90"
include "Dissipation.f90"
include "Timestep.f90"

subroutine  Solver
    implicit none
   integer::iter        !iterative variable
    integer::flag        !judge wheather the average density converge
    
    !write(*,*)  "Solver"
    
    do iter=1,itermax
        call Runge_kutta
        call Converge(flag)
        if(flag==1)  then
            write(*,*)  "...Converge..."
            exit 
        else
            write(*,*)  iter,"...Misconverge..."
            write(*,*)  "            .Continue..."
        end if
        write(*,*)
    end do
       
end subroutine 

subroutine Flow_init           !initialize flow
    implicit none
    
    !write(*,*)  "Flow_init"
    
    !原始量量初始化
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cos(att/180.0*pi)
    v_inf=Ma_inf*a_inf*sin(att/180.0*pi)
    
    rou=rou_inf
    p=p_inf
    u=u_inf
    v=v_inf
    
    !守恒量W的初始化
    W(:,:,1)=rou_inf
    W(:,:,2)=rou_inf*u_inf    
    W(:,:,3)=rou_inf*v_inf
    W(:,:,5)=p_inf/(gamma-1)+rou_inf*(u_inf**2+v_inf**2)/2.0   
    
end subroutine

subroutine Runge_kutta
    implicit none
    integer::m,mm
    real(8),dimension(0:imax+2,0:jmax+2,5)::W0,Fc,Dissi,R
    real(8)::dt(0:imax+2,0:jmax+2)
     
    !write(*,*)   "Runge_kutta"
    
    W0=W
    do m=1,Stage  
      call Boundary_cond                 !handle the boundary's conditions          
      call Convective_flux(Fc)              
      if(m==1)  then
         call Dissipation(Dissi,dt)                     
      end if
      
      do mm=1,5                           !only calculate the dissipation and time step at the first step
        R(:,:,mm)=-( Fc(:,:,mm)-Dissi(:,:,mm))/vol       
        W(:,:,mm)=W0(:,:,mm)+alpham(m)*dt*R(:,:,mm)
      end do
    rou=W(:,:,1)
    u=W(:,:,2)/W(:,:,1)
    v=W(:,:,3)/W(:,:,1)
    p=(gamma-1)*(  W(:,:,5)-rou*(  u**2+v**2  )/2.0  )     
    end do 
    
end subroutine

subroutine Converge(flag)           !judge wheather the mean density converge
    implicit none
    integer::i,j
    integer::flag                       !flag,1:converge;0:disconverge
    real(8)::rou_ncell=0.0         !mean density of n+1 layer
    real(8),save::rou_mean=1.225   !mean density of n layer
    
    !write(*,*)  "Converge"

    rou_ncell=sum(rou(2:imax,2:jmax))/( (imax-1)*(jmax-1) )

    flag=0
    if (abs(rou_ncell-rou_mean)<eps)   flag=1
    
    write(*,*)  abs(rou_ncell-rou_mean)
    
    rou_mean=rou_ncell
    
end subroutine

end module
    
    
    
    
    
    
    
  