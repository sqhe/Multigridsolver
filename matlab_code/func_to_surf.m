function A = func_to_surf( f, resol )
    % convert a vector to the 2d matrix form for drawing
    
    for i = 0 : resol
        for j = 0 : resol
            ind = ij2ind( i, j, resol );
            
            if f(ind) == 0 %abs( f(ind) ) < 1e-6;
                A( i+1, j+1 ) = inf;
            else
                A( i+1, j+1 ) = f(ind);
            end
        end
    end

end