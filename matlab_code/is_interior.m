function ret = is_interior( i, j, resol )
      
    [c_x, c_y, r, rot] = get_domain();
    
    R = [ cos(rot), -sin(rot);
          sin(rot), cos(rot) ];
    
    h = 1 / resol;
    x = i * h - c_x;
    y = j * h - c_y;
    
    xy = R * [x;y];
    x = xy(1);
    y = xy(2);
    
    if ~is_domain_disk()
        %------------- square domain
        if( abs( x ) < r && abs( y ) < r )
            ret = 1;
        else
            ret = 0;
        end
        %------------- square domain
    else
        %------------- disk domain
        assert( abs(x*x + y*y - r*r) > 1e-10 );
        
        if( x*x + y*y < r*r )
            ret = 1;
        else
            ret = 0;
        end
        %------------- disk domain
    end
    
end