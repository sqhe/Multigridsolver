function [ retA, retf ] = fill_bndy_exterior( i, j, resol, A, f, ind2dm, dm2ind )

    % first order boundary condition

    ind = ij2ind( i, j, resol );
    dm = ind2dm( ind );
    h = 1 / resol;
    scale = get_first_scale( resol );
            
    [c_x, c_y, r] = get_domain();
    
    pt_x = i * h;
    pt_y = j * h;
    
    [foot_x, foot_y, n_x, n_y] = get_foot_pt( pt_x, pt_y );
    
    dx = foot_x - pt_x;
    dy = foot_y - pt_y;
    
    if( dx > 0 )
        ngbr_i = i + 1;
    elseif( dx <= 0 )
        ngbr_i = i - 1;
    end
    
    if( dy > 0 )
        ngbr_j = j + 1;
    elseif( dy <= 0 )
        ngbr_j = j - 1;
    end
    
    ngbr_x = ngbr_i * h;
    ngbr_y = ngbr_j * h;
    
    dm_a = ind2dm( ij2ind( i, j, resol ) );
    dm_b = ind2dm( ij2ind( i, ngbr_j, resol ) );
    dm_c = ind2dm( ij2ind( ngbr_i, j, resol ) );
          
    %if ( ngbr_x - pt_x ) * n_x > 0
   %     A( dm, dm_a ) = A( dm, dm_a ) + (-1 ) / h * abs( n_x );
   %     A( dm, dm_c ) = A( dm, dm_c ) + ( 1 ) / h * abs( n_x );
    %else
    %    A( dm, dm_a ) = A( dm, dm_a ) + ( 1 ) / h * abs( n_x );
    %    A( dm, dm_c ) = A( dm, dm_c ) + (-1 ) / h * abs( n_x );
    %end
            
    %if ( ngbr_y - pt_y ) * n_y > 0        
    %    A( dm, dm_a ) = A( dm, dm_a ) + (-1 ) / h * abs( n_y );
    %    A( dm, dm_b ) = A( dm, dm_b ) + ( 1 ) / h * abs( n_y );
    %else
   %     A( dm, dm_a ) = A( dm, dm_a ) + ( 1 ) / h * abs( n_y );
    %    A( dm, dm_b ) = A( dm, dm_b ) + (-1 ) / h * abs( n_y );
    %end
    A(dm,dm)=1;
    
    %assert( abs( A(dm, dm_a) + A(dm, dm_b) + A(dm, dm_c) ) < 1e-6 );
    
    %f( dm ) = bd_con(  i, j, resol, ind2dm, dm2ind );
    f( dm ) = 0;
    %f( dm) = bound_con( i, j, resol, ind2dm, dm2ind );
    
    retA = A;
    retf = f;
    
end
