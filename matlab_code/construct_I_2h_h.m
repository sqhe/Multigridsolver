function I_2h_h = construct_I_2h_h( resol )

    % I_h_2h: map the grid from size 2h to size h

    resol2h = resol / 2;
    data_num_h = ( resol + 1 )^2;
    data_num_2h = ( resol2h + 1 )^2;

    %------------------------------------------------------
    % construct I_2h_h
    I_2h_h = sparse( data_num_h, data_num_2h );
      
    for ih = 0 : resol
        for jh = 0 : resol
        
            indh = ij2ind( ih, jh, resol );
            
            if( mod( ih, 2 ) == 0 && mod( jh, 2 ) == 0 )
                
                i2h = ih/2;
                j2h = jh/2;
                ind2h = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1;
                
            elseif( mod( ih, 2 ) == 0 )
                
                i2h = ih / 2;
                j2h = ( jh + 1 ) / 2;
                ind2h = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/2;
                
                j2h = ( jh - 1 ) / 2;
                ind2h = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/2;
                
            elseif( mod( jh, 2 ) == 0 )
                
                j2h = jh / 2;
                i2h = ( ih + 1 ) / 2;
                ind2h = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/2;
                
                i2h = ( ih - 1 ) / 2;
                ind2h = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/2;
                
            else
                
                i2h = ( ih + 1 ) / 2;
                j2h = ( jh + 1 ) / 2;
                ind2h( 1 ) = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/4;
                
                i2h = ( ih - 1 ) / 2;
                j2h = ( jh + 1 ) / 2;
                ind2h( 1 ) = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/4;
                
                i2h = ( ih + 1 ) / 2;
                j2h = ( jh - 1 ) / 2;
                ind2h( 1 ) = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/4;
                
                i2h = ( ih - 1 ) / 2;
                j2h = ( jh - 1 ) / 2;
                ind2h( 1 ) = ij2ind( i2h, j2h, resol2h );
                I_2h_h( indh, ind2h ) = 1/4;
            end
            
        end % jh = reosl
    end % ih = resol
    
end