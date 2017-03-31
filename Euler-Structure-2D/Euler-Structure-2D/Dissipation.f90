subroutine Dissipation(Dissi,dt)
    implicit none
    integer::i,j
    real(8),intent(out)::Dissi(0:imax+2,0:jmax+2,5),dt(0:imax+2,0:jmax+2)
    real(8)::u_av,v_av,rou_av,p_av           !original variables on edges
    real(8)::Nu1,Nu2,epsi2,epsi4             !pressure sensor,and coefficients of dissipation
    real(8)::ac,alpha,dis(5)
    
    !write(*,*) "Dissipation"
    Dissi=0.0   
    dt=0.0
    !in i direction
    do j=3,jmax
    do i=2,imax
        nu1=abs(p(i,j+1)-2*p(i,j)+p(i,j-1))/( p(i,j+1)+2*p(i,j)+p(i,j-1))
        nu2=abs(p(i,j)-2*p(i,j-1)+p(i,j-2))/( p(i,j)+2*p(i,j-1)+p(i,j-2))
        epsi2=k2*max(Nu1,Nu2)
        epsi4=max(0.0,k4-epsi2)
        
        u_av=(u(i,j)+u(i,j-1))/2
        v_av=(v(i,j)+v(i,j-1))/2
        p_av=(p(i,j)+p(i,j-1))/2
        rou_av=(rou(i,j)+rou(i,j-1))/2
        ac=sqrt(p_av*gamma/rou_av)
        alpha=abs(U_av*vector1(1,i,j)+V_av*vector1(2,i,j) ) + ac*ds1(i,j)
        
        dis=alpha*epsi2*(W(i,j,:)-W(i,j-1,:) )-alpha*epsi4*( W(i,j+1,:)-3*W(i,j,:) +3*W(i,j-1,:)-W(i,j-2,:)  )
        dissi(i,j-1,:)=dissi(i,j-1,:) + dis 
        dissi(i,j,:)= dissi(i,j,:) - dis    
        !the sum of every edge's spectral radius 
        dt(i,j)=dt(i,j)+alpha
        dt(i,j-1)=dt(i,j-1)+alpha
    end do
    end do
    !in j direction 
    do i=3,imax
    do j=2,jmax
        !nu=abs(p(i,j)-p(i-1,j))/abs(p(i,j)+p(i-1,j))
        nu1=abs(p(i+1,j)-2*p(i,j)+p(i-1,j))/( p(i+1,j)+2*p(i,j)+p(i-1,j))
        nu2=abs(p(i,j)-2*p(i-1,j)+p(i-2,j))/( p(i,j)+2*p(i-1,j)+p(i-2,j))
        epsi2=k2*max(nu1,nu2)
        epsi4=max(0.0,k4-epsi2)
        
        U_av=(u(i,j)+u(i-1,j))/2
        V_av=(v(i,j)+v(i-1,j))/2
        p_av=(p(i,j)+p(i-1,j))/2
        rou_av=(rou(i,j)+rou(i-1,j))/2
        ac=sqrt(p_av*gamma/rou_av)
        alpha=abs(U_av * vector2(1,i,j)+V_av * vector2(2,i,j) ) + ac*ds2(i,j)
        
        dis=alpha*epsi2*(W(i,j,:)-W(i-1,j,:) )-alpha*epsi4*( W(i+1,j,:)-3*W(i,j,:) +3*W(i-1,j,:)-W(i-2,j,:) )
        dissi(i-1,j,:)=dissi(i-1,j,:) + dis  
        dissi(i,j,:)=dissi(i,j,:) - dis  
        
        dt(i,j)=dt(i,j)+alpha
        dt(i-1,j)=dt(i,j)+alpha
    end do
    end do
     !time step
    dt=CFL * vol / dt
end subroutine
