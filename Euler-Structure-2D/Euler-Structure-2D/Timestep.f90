subroutine Timestep(dt)
    implicit none
    integer::i,j
    real(8)::dt(0:imax+2,0:jmax+2)
    real(8)::ac,alpha,vector_av(2),ds_av
    
    write(*,*)  "Timestep"
    
    dt=0.0      !归零
    !i方向扫略
    do j=2,jmax       
    do i=2,imax
      ac=sqrt(p(i,j)*gamma/rou(i,j))
      vector_av(:)=0.5*(vector1(:,i,j)/ds1(i,j)+vector1(:,i,j+1)/ds1(i,j+1) )
      ds_av=0.5*(ds1(i,j)+ds1(i,j+1))
      alpha=( abs(U(i,j)*vector_av(1)+V(i,j)*vector_av(2) ) + ac ) *ds_av   
      dt(i,j)=dt(i,j)+alpha
    end do
    end do
    !****************************************
    !j方向扫略
    do i=2,imax        
    do j=2,jmax
      ac=sqrt(p(i,j)*gamma/rou(i,j))
      vector_av(:)=0.5*(vector2(:,i,j)/ds2(i,j)+vector2(:,i+1,j)/ds2(i+1,j) )
      ds_av=0.5*(ds2(i,j)+ds2(i+1,j))
      alpha=( abs(u(i,j)*vector_av(1)+v(i,j)*vector_av(2) ) + ac )*ds_av   
      dt(i,j)=dt(i,j)+alpha
    end do
    end do
      dt=CFL * vol / dt
end subroutine

!subroutine Timestep(dt)   !按照非结构网格添加的时间步长算法，按边循环
   ! implicit none
   ! integer::i,j
   ! real(8)::dt(0:imax+2,0:jmax+2)
   ! real(8)::u_av,v_av,rou_av,p_av
   ! real(8)::ac,alpha
    
   ! write(*,*)  "Timestep"
    
    !dt=0.0      !归零
    !i方向扫略
    !do j=2,jmax+1        
    !do i=2,imax
      !U_av=(u(i,j)+u(i,j-1))/2
      !V_av=(v(i,j)+v(i,j-1))/2
      !p_av=(p(i,j)+p(i,j-1))/2
      !rou_av=(rou(i,j)+rou(i,j-1))/2
      !ac=sqrt(p_av*gamma/rou_av)
      !alpha=abs(U_av*vector1(1,i,j)+V_av*vector1(2,i,j) ) + ac*ds1(i,j)   
      !dt(i,j)=dt(i,j)+alpha
      !dt(i,j-1)=dt(i,j-1)+alpha
      
   ! end do
    !end do
    !****************************************
    !j方向扫略
    !do i=2,imax          !imax+1
    !do j=2,jmax
    !  U_av=(u(i,j)+u(i-1,j))/2
     ! V_av=(v(i,j)+v(i-1,j))/2
     ! p_av=(p(i,j)+p(i-1,j))/2
     ! rou_av=(rou(i,j)+rou(i-1,j))/2
     ! ac=sqrt(p_av*gamma/rou_av)
      
     ! alpha=abs(U_av*vector2(1,i,j)+V_av*vector2(2,i,j) ) + ac *ds2(i,j)   
     ! dt(i,j)=dt(i,j)+alpha
     ! dt(i-1,j)=dt(i,j)+alpha
      
   ! end do
   ! end do
    !  dt=CFL * vol / dt
!end subroutine



