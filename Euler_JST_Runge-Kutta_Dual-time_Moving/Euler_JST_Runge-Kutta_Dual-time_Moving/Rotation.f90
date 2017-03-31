subroutine Rotation(AoA,angular_rate)  !  the next step's node speed and coordinates
    implicit none 
    real(8)::AoA   !angle of attack
    real(8)::angular_rate
     
    xy(1,:) =   xy0(1,:) * cosd(AoA-att) + xy0(2,:) * sind(AoA-att) 
    xy(2,:) = - xy0(1,:) * sind(AoA-att) + xy0(2,:) * cosd(AoA-att) 
         
    U_Rot(1,:) =   angular_rate *  xy(2,:)   
    U_Rot(2,:) = - angular_rate *  xy(1,:)      
              
  
    vector(1,:) =  vector0(1,:) * cosd( AoA-att ) + vector0(2,:) * sind( AoA-att )
    vector(2,:) = -vector0(1,:) * sind( AoA-att ) + vector0(2,:) * cosd( AoA-att )   
              
end subroutine

    