function ret = is_real_interior( i, j, resol )
    
    if( ~is_interior( i, j, resol ) )
        ret = 0;
    else
        
        ngbr = get_direct_ngbr( i, j, resol );
        for ii = 1 : 4
            [ni, nj] = ind2ij( ngbr(ii), resol );
                        
            if( ~is_interior( ni, nj, resol ) )
                ret = 0;
                return;
            end
        end
        
        ret = 1;
    end
    
end