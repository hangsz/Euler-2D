module  Control_para    !control parameters
    implicit none
    !**********************************************
    !constant parameters
    real(8),parameter::pi=3.14159265
    
    !**********************************************
    !grid parameters
    integer,parameter::imax=193,jmax=33              !!the number of points  in x and y direciton respectively
    integer,parameter::airfoil_num=127               !the number of points on the aerofoil
    
    !**********************************************
    !incoming flow information
    real(8),parameter::Ma_inf=0.8               !the mach number
    real(8),parameter::att=1.25                 !attack angle,бу
    real(8),parameter::rou_inf=1.225            !density,kg/m^3
    real(8),parameter::p_inf=1.01325E5          !pressure,N/m^2
    !real(8),parameter::u_inf=10.0,v_inf=0.0    !teh velocities in x and y directions,m/s
    real(8),parameter::R=287.0                  !the ideal gas constant,J/(kg*K) 
    real(8),parameter::gamma=1.4                !specific heat ratio
    !**********************************************
    !Finite volume method,Jamson scheme
    integer,parameter::itermax=10000                  !the maximun iterative times
    real(8),parameter::k2=0.6,k4=1.0/64         !artifical dissipation,k2(1/2,1),k4(1/256,1/32)
    real(8),parameter::CFL=2.0                  !CFL number
    real(8),parameter::eps=1.0E-10               !the convergence precision  of the average density
    !adjustable parameters of Runge_kutta scheme
    integer,parameter::stage=4                       !Runge_kutta's step number
    real(8),parameter::alpham(stage)=(/ 1.0/4, 1.0/3, 1.0/2, 1.0 /)    !optimised stage coefficients
    
end module