subroutine rotation(AoA,angular_rate) 
    !------------------------------------------
       ! purpose:  renew all the geometric information with grid rotating. 
    !------------------------------------------
       
    implicit none 
    real(8)::AoA                     ! attact angle              
    real(8)::angular_rate            ! angular velocity 
    !------------------------------------------
    ! all the geometric information has to be renewed with grid rotating
    xy(1,:) =   xy0(1,:) * cosd(AoA-att) + xy0(2,:) * sind(AoA-att) 
    xy(2,:) = - xy0(1,:) * sind(AoA-att) + xy0(2,:) * cosd(AoA-att) 
         
    U_Rot(1,:) =  angular_rate *  xy(2,:)   
    U_Rot(2,:) = -angular_rate *  xy(1,:)      
              
    vector(1,:) =  vector0(1,:) * cosd( AoA-att ) + vector0(2,:) * sind( AoA-att )
    vector(2,:) = -vector0(1,:) * sind( AoA-att ) + vector0(2,:) * cosd( AoA-att )   
              
end subroutine

    