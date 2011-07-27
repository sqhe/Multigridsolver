function I_h_2h = construct_I_h_2h( resol, ind2dm_h )

    % I_h_2h: map the grid from size h to size 2h

    resol2h = resol / 2;
    data_num_h = ( resol + 1 )^2;
    data_num_2h = ( resol2h + 1 )^2;
  
    I_h_2h = sparse( data_num_2h, data_num_h );
    
    for i2h = 0 : resol2h
        for j2h = 0 : resol2h
            
            ind2h = ij2ind( i2h, j2h, resol2h );
            ih = i2h * 2;
            jh = j2h * 2;
            indh = ij2ind( ih, jh, resol );
            
            if is_bndy_exterior( i2h, j2h, resol2h )
                
                self = [];
                [ni, nj] = ind2ij( indh, resol );
                if is_bndy_exterior( ni, nj, resol )
                    self = indh;
                end
                
                dngbr_tmp = get_direct_ngbr( ih, jh, resol );
                cngbr_tmp = get_corner_ngbr( ih, jh, resol );
                               
                dngbr = [];
                for i = 1 : 4
                    [ni, nj] = ind2ij( dngbr_tmp(i), resol );
                    if is_bndy_exterior( ni, nj, resol )
                        dngbr = [ dngbr, dngbr_tmp(i) ];
                    end
                end
                
                cngbr = [];
                for i = 1 : 4
                    [ni, nj] = ind2ij( cngbr_tmp(i), resol );
                    if is_bndy_exterior( ni, nj, resol )
                        cngbr = [ cngbr, cngbr_tmp(i) ];
                    end
                end
                                
                nonempty = 0;
                [m,n]=size( dngbr ); if n > 0 nonempty = nonempty + 1; end
                [m,n]=size( cngbr ); if n > 0 nonempty = nonempty + 1; end
                [m,n]=size( self ); if n > 0 nonempty = nonempty + 1; end
                assert( nonempty > 0 );
                       
                [m,n] = size( self );
                if n>0 I_h_2h( ind2h, self ) = 1/nonempty; end
                
                [m,n] = size( dngbr );
                if n>0 I_h_2h( ind2h, dngbr ) = 1/nonempty/n; end
                
                [m,n] = size( cngbr );
                if n>0 I_h_2h( ind2h, cngbr ) = 1/nonempty/n; end
                
                assert( abs( sum(I_h_2h(ind2h,:)) - 1 ) < 1e-10 );
                
            elseif is_bndy_interior( i2h, j2h, resol2h )
                
                self = [];
                [ni, nj] = ind2ij( indh, resol );
                if is_bndy_interior( ni, nj, resol )
                    self = indh;
                end
                
                dngbr_tmp = get_direct_ngbr( ih, jh, resol );
                cngbr_tmp = get_corner_ngbr( ih, jh, resol );
                
                dngbr = [];
                for i = 1 : 4
                    [ni, nj] = ind2ij( dngbr_tmp(i), resol );
                    if is_bndy_interior( ni, nj, resol )
                        dngbr = [ dngbr, dngbr_tmp(i) ];
                    end
                end
                
                cngbr = [];
                for i = 1 : 4
                    [ni, nj] = ind2ij( cngbr_tmp(i), resol );
                    if is_bndy_interior( ni, nj, resol )
                        cngbr = [ cngbr, cngbr_tmp(i) ];
                    end
                end
                                
                nonempty = 0;
                [m,n]=size( dngbr ); if n > 0 nonempty = nonempty + 1; end
                [m,n]=size( cngbr ); if n > 0 nonempty = nonempty + 1; end
                [m,n]=size( self ); if n > 0 nonempty = nonempty + 1; end
%                 if nonempty > 2 || nonempty == 0
%                     nonempty
%                 end
                assert( nonempty > 0 );
                
                [m,n] = size( self );
                if n>0 I_h_2h( ind2h, self ) = 1/nonempty; end
                
                [m,n] = size( dngbr );
                if n>0 I_h_2h( ind2h, dngbr ) = 1/nonempty/n; end
                
                [m,n] = size( cngbr );
                if n>0 I_h_2h( ind2h, cngbr ) = 1/nonempty/n; end
             
                assert( abs( sum(I_h_2h(ind2h,:)) - 1 ) < 1e-10 );
                
            elseif is_real_interior( i2h, j2h, resol2h )
                
                I_h_2h( ind2h, indh ) = 1/4;
                assert( ind2dm_h( indh ) ~= 0 );

                ngbr = get_direct_ngbr( ih, jh, resol );
                I_h_2h( ind2h, ngbr ) = 1/8;
                assert( ind2dm_h( ngbr ) ~= 0 );

                ngbr = get_indirect_ngbr( ih, jh, resol ); %%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                I_h_2h( ind2h, ngbr ) = 1/16;
                assert( ind2dm_h( ngbr ) ~= 0 );
                
                assert( abs( sum(I_h_2h(ind2h,:)) - 1 ) < 1e-10 );
                
            end
        end
    end
    
end
