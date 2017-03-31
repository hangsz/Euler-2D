module FVM                                           !finite volume method 
    use Control_para
    use Grid_info
    
    implicit none 
    
    real(8)::u_inf,v_inf,a_inf                       !incoming flow's velocities and sound velocity
  
    real(8),allocatable::U(:,:),&
                         W(:,:),Wn(:,:),Wn1(:,:),Q(:,:)   , &                 !the conservative variables of every cell
                         W0(:,:)  , &                 !the zeroth step conservative variables
                         Rsi(:,:) , &                 !residual 
                         Fc(:,:)  , &                 !convective flux of every cell
                         alf(:)   , &                 !spectral radius of every edge
                         Dissi(:,:)  ,D_last(:,:) ,&              !artificial dissipation of every cell
                         dt(:), & !time step of every cell
                         U_av(:,:)
    
     !------------------------
     !rotate
     real(8),allocatable::U_Rot(:,:)
     real(8),allocatable::xy0(:,:),vector0(:,:)
     
     !----------------------------
     real(8)::oldRsi =1.0E-10
     
contains

include "Read_flow.f90"

include "Mean_edge.f90"
include "Con_flux.f90"
include "Art_dissipation.f90"

include "Rotation.f90"

include "Runge_Kutta.f90"

include "Output.f90"

include "outputFreeFlow.f90"


subroutine Allocate_memory    

    !allocate memory     
    allocate( U(5,ncells) )
    
    allocate( W(5,ncells) )
    allocate( Wn(5,ncells) )
    allocate( Wn1(5,ncells) )
    allocate( Q(5,ncells) )
    
    allocate( W0(5,ncells) )
    allocate( Rsi(5,ncells) )
    allocate( Fc(5,ncells) )
    allocate( alf(nedges) )
    allocate( Dissi(5,ncells) )
    allocate( D_last(5,ncells) )
    allocate( dt(ncells) )
    
    allocate( U_av(6,nedges) )
    
    !moving
    allocate( U_Rot(2,nnodes) )
    allocate( xy0(2,nnodes) )
    allocate( vector0(2,nedges) )
    
end subroutine

subroutine Flow_init
   
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
    
    W(1,:) = U(1,:)
    W(2,:) = U(1,:)*U(2,:)
    W(3,:) = U(1,:)*U(3,:)
    W(5,:)=  U(5,:) /(gamma-1.0) + U(1,:)*( U(2,:)**2 + U(3,:)**2 ) / 2.0 
    

end subroutine

subroutine Solver          !the Solver
    implicit none
    integer::i,j
    integer::count
    integer::iter          !iterative variable
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=30):: filename!= "flow_info-.dat"
    logical::alive
    
    real(8)::t_total = 0.0
    real(8)::omg
    real(8)::AoA,angular_rate
    
    write(*,*)  "Solver"
    
    
    call Grid
    call Allocate_memory
    call Read_flow
    
    call Flow_init
    
  
    !geometry set
    !translate the grid to where the rotational center's coordinates are zero
    
    xy0 = xy
    
    vector0 = vector
    
    do i=1,nnodes
        xy0(:,i) = xy0(:,i)- rot_cen(:)
    end do 
    !----------------------
    
    
    Wn =  W
    Wn1 = Wn
    
    omg =  2.0*kr*Ma_inf*a_inf
    
    do i =1,phase  
        do iter =1,itermax
            
            Wn1= Wn
            Wn = W
            do j=1,5
                Q(j,:) = 2.0/dt_r*vol*Wn(j,:) - 1.0/2/dt_r*vol*Wn1(j,:)    
            end do
            
            t_total = t_total + dt_r
            
            AoA = att + att_ampl* sin( omg * t_total )
            angular_rate = omg * ( att_ampl/180.0 *pi ) *cos( omg * t_total)
          
        
            call Rotation(AoA,angular_rate)  !renew the node speed and coordinates
           
            do count=1,iter_inner 
                write(*,*)  count,iter,i
                write(*,*) "t:",t_total,"dt:",dt_r
                write(*,*) "AoA:",AoA
                
                call Runge_Kutta
                call Converge(flag)
                
                if( flag==1)  exit
            end do
            
            oldRsi = sum(abs(Rsi(1,:)))/ncells
            
        end do
        
        write(filename,"(I2)") i
            
        filename = "flow-info-"//trim(filename)
        call Output(filename) 
        write(*,*) "Output" 
        
   end do
   
     
end subroutine 

!subroutine Converge(flag)          !verify wheather the flow converge
!    implicit none
!    integer::i,j
!    integer::flag                  !flag, 1:converge;0:disconverge
!   
!    real(8)::Rsi_old = 0.0,Rsi_new
!    !write(*,*)  "Converge"
!
!    
!    Rsi_new = sum ( Rsi(1,:) )/ncells
!    
!    flag = 0 
!    
!    if(abs(Rsi_new) .LE. eps ) flag = 1
!    
!    write(*,*) abs(Rsi_new) 
!           
!    Rsi_old = Rsi_new
!    
!end subroutine

subroutine Converge(flag)          !verify wheather the flow converge
    implicit none
    integer::i,j
    integer::flag                  !flag, 1:converge;0:disconverge
    real(8)::rou_ncell=0.0    !the mean density of n+1 layer
    real(8)::u_ncell=0.0    !the mean density of n+1 layer
    real(8)::v_ncell=0.0    !the mean density of n+1 layer
    real(8)::p_ncell=0.0    !the mean density of n+1 layer
      
    real(8),save::rou_mean = 1.225  !the mean density of n layer
    real(8),save::u_mean = 0.0  !the mean density of n layer
    real(8),save::v_mean = 0.0  !the mean density of n laye
    real(8),save::p_mean = 103150.0  !the mean density of n layer


    real(8)::newRsi =0.0
    !write(*,*)  "Converge"

    rou_ncell = sum(U(1,:))/ncells
    u_ncell = sum(U(2,:))/ncells
    v_ncell = sum(U(3,:))/ncells
    p_ncell = sum(U(5,:))/ncells
    
    flag = 0
     
    newRsi = sum(abs(Rsi(1,:)))/ncells
    
    if( newRsi/oldRsi<1.0E-2 )   flag=1
    
    
    write(*,*)  U(1,1),U(2,1),U(3,1),U(5,1)
    write(*,*)  abs(rou_ncell-rou_mean),abs(p_ncell-p_mean)
    write(*,*)  abs(u_ncell-u_mean),abs(v_ncell-v_mean)
    write(*,*)
       
    rou_mean = rou_ncell
    u_mean = u_ncell
    v_mean = v_ncell
    p_mean = p_ncell
    
end subroutine

end module
    