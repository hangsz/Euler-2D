subroutine  Convective_flux(Fc)
    implicit none
    integer::i,j
    real(8),intent(out)::Fc(0:imax+2,0:jmax+2,5)
    real(8)::W1,W2,W3,W5              !original variables on edges
    real(8)::U_av,V_av,p_av           !conservative variables on edges
    real(8)::Z       
    
    ! write(*,*) "Convective_flux"
    
    Fc=0.0  
    !in i direction 
    do j=2,jmax+1
    do i=2,imax
        p_av=(p(i,j)+p(i,j-1))/2
        U_av=(u(i,j)+u(i,j-1))/2
        V_av=(v(i,j)+v(i,j-1))/2
        
        W1=(W(i,j,1)+W(i,j-1,1))/2   
        W2=(W(i,j,2)+W(i,j-1,2))/2   
        W3=(W(i,j,3)+W(i,j-1,3))/2
        W5=(W(i,j,5)+W(i,j-1,5))/2   
        W5=W5+p_av                   !rouH=rouE+p
        
        Z=U_av*vector1(1,i,j)+V_av*vector1(2,i,j)
        !the right cell subtract
        Fc(i,j-1,1)=Fc(i,j-1,1) - Z*W1
        Fc(i,j-1,2)=Fc(i,j-1,2) -(Z*W2 + p_av*vector1(1,i,j) )
        Fc(i,j-1,3)=Fc(i,j-1,3) -(Z*W3 + p_av*vector1(2,i,j) )
        Fc(i,j-1,5)=Fc(i,j-1,5) - Z*W5
        !the left cell puls
        Fc(i,j,1)=Fc(i,j,1)+ Z*W1
        Fc(i,j,2)=Fc(i,j,2)+ (Z*W2+ p_av*vector1(1,i,j) )
        Fc(i,j,3)=Fc(i,j,3)+ (Z*W3+ p_av*vector1(2,i,j) )
        Fc(i,j,5)=Fc(i,j,5)+ Z*W5
    end do
    end do
    
    !in j direction
    do i=2,imax+1
    do j=2,jmax
        p_av=(p(i,j)+p(i-1,j))/2
        U_av=(u(i,j)+u(i-1,j))/2
        V_av=(v(i,j)+v(i-1,j))/2
        
        W1=(W(i,j,1)+W(i-1,j,1))/2   
        W2=(W(i,j,2)+W(i-1,j,2))/2   
        W3=(W(i,j,3)+W(i-1,j,3))/2
        W5=(W(i,j,5)+W(i-1,j,5))/2   
        W5=W5+p_av 
        Z=U_av*vector2(1,i,j)+V_av*vector2(2,i,j)
        !the left cell plus
        Fc(i-1,j,1)=Fc(i-1,j,1)+ Z*W1
        Fc(i-1,j,2)=Fc(i-1,j,2)+ (Z*W2+ p_av*vector2(1,i,j) )
        Fc(i-1,j,3)=Fc(i-1,j,3)+ (Z*W3+ p_av*vector2(2,i,j) )
        Fc(i-1,j,5)=Fc(i-1,j,5)+ Z*W5
        !the rigt cell subtract
        Fc(i,j,1)=Fc(i,j,1)- Z*W1
        Fc(i,j,2)=Fc(i,j,2)- (Z*W2+ p_av*vector2(1,i,j) )
        Fc(i,j,3)=Fc(i,j,3)- (Z*W3+ p_av*vector2(2,i,j) )
        Fc(i,j,5)=Fc(i,j,5)- Z*W5
    end do
    end do
    
end subroutine
