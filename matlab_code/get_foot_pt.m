function [footx, footy, nx, ny] = get_foot_pt( x, y )

    [c_x, c_y, r, rot] = get_domain();
    
    R = [ cos(rot), -sin(rot);
          sin(rot), cos(rot) ];
      
    x = x - c_x;
    y = y - c_y;
    
    xy = R * [x;y];
    x = xy( 1 ); 
    y = xy( 2 );
       
    if ~is_domain_disk()
        %------------- square domain
        if x~= 0
            slp = y / x;
            
            if ( abs(x) > r && abs(y) > r ) || ( abs( abs(slp) - 1 ) < 1e-10 )
                
                if x > 0
                    footx = r;
                    nx = 1 / sqrt( 2 );
                else
                    footx = -r;
                    nx = -1 / sqrt( 2 );
                end
                
                if y > 0 
                    footy = r;
                    ny = 1 / sqrt( 2 );
                else
                    footy = -r;
                    ny = -1 / sqrt( 2 );
                end
                
                %[ nx, ny ]
                                
            elseif abs( slp ) > 1
                
                footx = x;
                nx = 0;
                if y > 0 
                    ny = 1; 
                    footy = r;
                else 
                    ny = -1; 
                    footy = -r;
                end
                
            else
                footy = y;
                ny = 0;
                if x > 0 
                    nx = 1; 
                    footx = r;
                else 
                    nx = -1; 
                    footx = -r;
                end
            end

        else
            
            if y > 0
                footx = 0;
                footy = r;
                nx = 0;
                ny = 1;
            else
                footx = 0;
                footy = -r;
                nx = 0;
                ny = -1;
            end
            
        end
        %------------- square domain
    else
        %------------- disk domain
        nx = x;
        ny = y;
        len = sqrt( nx*nx + ny*ny );
        nx = nx / len;
        ny = ny / len;

        footx = nx * r;
        footy = ny * r;
        %------------- disk domain
    end
    
    foot = inv(R) * [ footx; footy ];
    footx = foot(1);
    footy = foot(2);
    
    n = inv(R) * [nx; ny];
    nx = n(1);
    ny = n(2);
    
    footx = footx + c_x;
    footy = footy + c_y;
    
end