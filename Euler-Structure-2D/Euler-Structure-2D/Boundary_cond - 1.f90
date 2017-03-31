subroutine Boundary_cond     !boundary conditons, determine the flow information in dummy cell
    implicit none
    integer::i,j
    real(8)::u_b,v_b,p_b,rou_b                  !the original variables on boundaries
    real(8)::W_b(5)                             !the consevative variables on boundaries
    real(8)::Z,Ma,rou0,c0
    
    !write(*,*)  "Boundary_cond"
    !**************************************************
    !bottom boundary
    !aerofoil surface
    do i=(imax-airfoil_num)/2+1,(imax+airfoil_num)/2+1  
        !there's no flow normal to the surface,v.n=0
        !atain the pressure and density by extrapolation from the interior
        p_b=0.5*(3*p(i,2)-p(i,3))   
        rou_b=0.5*(3*rou(i,2)-rou(i,3)) 
        !the first layer dummy cell
        p(i,1)=2*p_b-p(i,2) 
        rou(i,1)=2*rou_b-rou(i,2) 
        u(i,1)=-u(i,2)     !the real u and v are not necessary, it's enough to be compatible to solid wall conditon
        v(i,1)=-v(i,2) 
        !the second layer dummy cell 
        p(i,0)=2*p(i,1)-p(i,2)
        W(i,1,:)=2.0*W(i,2,:)-W(i,3,:)          
        W(i,0,:)=2.0*W(i,1,:)-W(i,2,:)  
    end do 
           
    !point isn't on aerofoil
    !the flow information in dummy cells is set with the value of the opposite real cell 
    do i=2,(imax-airfoil_num)/2
        !the first layer of dummy cells
        u(i,1)=u(imax-i+2,2)
        v(i,1)=v(imax-i+2,2)
        p(i,1)=p(imax-i+2,2)
        rou(i,1)=rou(imax-i+2,2)
        !the second layer of dummy cells
        p(i,0)=p(i,3)
        W(i,1,:)=W(imax-i+2,2,:)
        W(i,0,:)=W(imax-i+2,3,:)
      
        u(imax-i+2,1)=u(i,2)
        v(imax-i+2,1)=v(i,2)
        p(imax-i+2,1)=p(i,2)
        rou(imax-i+2,1)=rou(i,2)
      
        p(imax-i+2,0)=p(i,3)
        W(imax-i+2,1,:)=W(i,2,:)
        W(imax-i+2,0,:)=W(i,3,:)
    end do
  
    !*********************************************************** 
    !top boundary(j=jmax+1)
    do i=2,imax
        rou0=rou(i,jmax)
        c0=sqrt(p(i,jmax)*gamma/rou(i,jmax))
        Z=u(i,jmax)*vector1(1,i,jmax+1)+v(i,jmax)*vector1(2,i,jmax+1)
        Ma=sqrt(u(i,jmax)**2+v(i,jmax)**2)/c0
        if (Z .GE. 0.0)   then        !inflow
            if( Ma .GE. 1.0 )  then   !supersonic
                p_b=p_inf
                rou_b=rou_inf
                u_b=u_inf
                v_b=v_inf
            else                      !subsonic
                p_b=0.5*( p_inf+p(i,jmax)+rou0*c0/ds1(i,jmax+1)*((u_inf-u(i,jmax))*vector1(1,i,jmax+1)+(v_inf-v(i,jmax))*vector1(2,i,jmax+1) )   )
                rou_b=rou_inf+(p_b-p_inf)/c0**2
                u_b=u_inf+(p_inf-p_b)*vector1(1,i,jmax+1)/ds1(i,jmax+1)/(rou0*c0)
                v_b=v_inf+(p_inf-p_b)*vector1(2,i,jmax+1)/ds1(i,jmax+1)/(rou0*c0)  
            end if
            
        else                           !outflow
             if( Ma .GE. 1.0 )  then   !supersonic
                 p_b=p(i,jmax)
                 rou_b=rou(i,jmax)
                 u_b=u(i,jmax)
                 v_b=v(i,jmax)
             else                   !subsonic
                p_b=p_inf
                rou_b=rou(i,jmax)+(p_b-p(i,jmax))/c0**2
                U_b=u(i,jmax)-(p(i,jmax)-p_b)*vector1(1,i,jmax+1)/ds1(i,jmax+1)/(rou0*c0)
                V_b=v(i,jmax)-(p(i,jmax)-p_b)*vector1(2,i,jmax+1)/ds1(i,jmax+1)/(rou0*c0)
             end if
        end if 
      
        rou(i,jmax+1)=2*rou_b-rou(i,jmax) 
        p(i,jmax+1)=2*p_b-p(i,jmax)
        u(i,jmax+1)=2*U_b-u(i,jmax) 
        v(i,jmax+1)=2*V_b-v(i,jmax) 
        p(i,jmax+2)=2*p(i,jmax+1)-p(i,jmax)
        
        W_b(1)=rou_b
        W_b(2)=rou_b*U_b
        W_b(3)=rou_b*V_b
        W_b(5)=p_b/(gamma-1)+rou_b*(U_b**2+V_b**2)/2.0 
        
        W(i,jmax+1,:)=2*W_b-W(i,jmax,:)
        
        W(i,jmax+2,:)=2*W(i,jmax+1,:)-W(i,jmax,:)  
    end do
    
     !i=2
    do j=2,jmax
        rou0=rou(2,j)
        c0=sqrt(p(2,j)*gamma/rou(2,j))
        Z=u(2,j)*vector2(1,2,j)+v(2,j)*vector2(2,2,j)
        Ma=sqrt(u(2,j)**2+v(2,j)**2)/c0
        if(Z .LE. 0.0)  then          !outflow
            if( Ma .GE. 1.0)  then    !supersonic
                p_b=p(2,j)
                rou_b=rou(2,j)
                u_b=u(2,j)
                v_b=v(2,j)
            else                       !subsonic
                p_b=p_inf
                rou_b=rou(2,j)+(p_b-p(2,j))/c0**2
                u_b=u(2,j)-(p(2,j)-p_b)*vector2(1,2,j)/ds2(2,j)/(rou0*c0)
                v_b=v(2,j)-(p(2,j)-p_b)*vector2(2,2,j)/ds2(2,j)/(rou0*c0)
            end if
        else                          !inflow
            if( Ma .GE. 1.0)  then    !supersonic
                p_b=p_inf
                rou_b=rou_inf
                u_b=u_inf
                v_b=v_inf
            else                      !subsonic
                p_b=0.5*( p_inf+p(2,j)+rou0*c0/ds2(2,j)*((u_inf-u(2,j))*vector2(1,2,j)+(v_inf-v(2,j))*vector2(2,2,j) )   )
                rou_b=rou_inf+(p_b-p_inf)/c0**2
                u_b=u_inf+(p_inf-p_b)*vector2(1,2,j)/ds2(2,j)/(rou0*c0)
                v_b=v_inf+(p_inf-p_b)*vector2(2,2,j)/ds2(2,j)/(rou0*c0)
            end if
            
        end if      
        rou(1,j)=2*rou_b-rou(2,j) 
        p(1,j)=2*p_b-p(2,j)
        u(1,j)=2*U_b-u(2,j)
        v(1,j)=2*V_b-v(2,j) 
        p(0,j)=2*p(1,j)-p(2,j)
        
        W_b(1)=rou_b
        W_b(2)=rou_b*U_b
        W_b(3)=rou_b*V_b
        W_b(5)=p_b/(gamma-1)+rou_b*(U_b**2+V_b**2)/2.0 
        
        W(1,j,:)=2*W_b-W(2,j,:)
        W(0,j,:)=2*W(1,j,:)-W(2,j,:)
    end do
    
    ! i=imax+1
    do j=2,jmax
        rou0=rou(imax,j)
        c0=sqrt(p(imax,j)*gamma/rou(imax,j))
        Z=u(imax,j)*vector2(1,imax+1,j)+v(imax,j)*vector2(2,imax+1,j)
        Ma=sqrt(u(imax,j)**2+v(imax,j)**2)/c0
        if( Z .GE. 0.0)  then         !outflow
            if( Ma .GE. 1.0)  then    !supersonic
                p_b=p(imax,j)
                rou_b=rou(imax,j)
                u_b=u(imax,j)
                v_b=v(imax,j)
            else                     !subsonic
                p_b=p_inf           
                rou_b=rou(imax,j)+(p_b-p(imax,j))/c0**2
                u_b=u(imax,j)+(p(imax,j)-p_b)*vector2(1,imax+1,j)/ds2(imax+1,j)/(rou0*c0)
                v_b=v(imax,j)+(p(imax,j)-p_b)*vector2(2,imax+1,j)/ds2(imax+1,j)/(rou0*c0)
            end if
        else                           !inflow
            if( Ma .GE. 1.0 )  then    !supersoinc
                p_b=p_inf
                rou_b=rou_inf
                u_b=u_inf
                v_b=v_inf
            else                       !subsonic
                p_b=0.5*( p_inf+p(imax,j)-rou0*c0/ds2(imax+1,j)*(  (u_inf-u(imax,j))*vector2(1,imax,j)+(v_inf-v(imax,j))*vector2(2,imax+1,j) )   )
                rou_b=rou_inf+(p_b-p_inf)/c0**2
                u_b=u_inf - (p_inf-p_b)*vector2(1,imax+1,j)/ds2(imax+1,j)/(rou0*c0)
                v_b=v_inf - (p_inf-p_b)*vector2(2,imax+1,j)/ds2(imax+1,j)/(rou0*c0)   
            end if
        end if
        rou(imax+1,j)=2*rou_b -rou(imax,j)
        p(imax+1,j)=2*p_b-p(imax,j)
        u(imax+1,j)=2*U_b -u(imax,j)
        v(imax+1,j)=2*V_b-v(imax,j) 
        p(imax+2,j)=2*p(imax+1,j)-p(imax,j)
        
        W_b(1)=rou_b
        W_b(2)=rou_b*U_b
        W_b(3)=rou_b*V_b
        W_b(5)=p_b/(gamma-1)+rou_b*( U_b**2+V_b**2)/2.0 
        
        W(imax+1,j,:)=2*W_b-W(imax,j,:)
        W(imax+2,j,:)=2*W(imax+1,j,:)- W(imax,j,:)
    end do
end subroutine