function [ retA, retf ] = fill_bndy_interior( i, j, resol, A, f, ind2dm, dm2ind )

    % second order boundary condition
    
    ind = ij2ind( i, j, resol );
    dm = ind2dm( ind );
    h = 1 / resol;
        
    pt_x = i * h;
    pt_y = j * h;
    
    [foot_x, foot_y, n_x, n_y] = get_foot_pt( pt_x, pt_y );
    
    dx = foot_x - pt_x;
    dy = foot_y - pt_y;
    
    if( dx > 0 )
        ngbr_i_inc = 1;
    elseif( dx <= 0 )
        ngbr_i_inc = -1;
    end
    
    if( dy > 0 )
        ngbr_j_inc = 1;
    elseif( dy <= 0 )
        ngbr_j_inc = -1;
    end
    
    i_arr( 1 ) = i;                 j_arr( 1 ) = j;
    i_arr( 2 ) = i;                 j_arr( 2 ) = j+1;
    i_arr( 3 ) = i;                 j_arr( 3 ) = j-1;
    i_arr( 4 ) = i-1;               j_arr( 4 ) = j;
    i_arr( 5 ) = i+1;               j_arr( 5 ) = j;
    i_arr( 6 ) = i + ngbr_i_inc;    j_arr( 6 ) = j + ngbr_j_inc;
    
    for ii = 1 : 6
        ind_arr = ij2ind( i_arr(ii), j_arr(ii), resol );
        dm_arr( ii ) = ind2dm( ind_arr );
        x = i_arr( ii ) - i;
        y = j_arr( ii ) - j;
        
        arr_x( ii ) = x;
        arr_y( ii ) = y;
        
        P( 1, ii ) = x*x;
        P( 2, ii ) = y*y;
        P( 3, ii ) = x*y;
        P( 4, ii ) = x;
        P( 5, ii ) = y;
        P( 6, ii ) = 1;
    end
   
%     det(P*P');
%     assert( abs( det(P*P') ) > 1e-6 );
    W = inv( P' ); %inv(P*P')*P;
    
    a = 2 * W( 1, : ) * n_x * n_x;
    b = 2 * W( 2, : ) * n_y * n_y;
    c = 2 * W( 3, : ) * n_x * n_y;
    w = -( a + b + c ) / (h^2);
        
    %temp=get_direct_ngbr(i,j,resol);
    
    %dm_arr=[dm];
    
    %for i=1:4
    %    tsd=ind2dm(temp(i));
    %    dm_arr=[dm_arr tsd];
    %end
    
    %w=[4 -1 -1 -1 -1];
    
    %A( dm, dm_arr ) = w;
    dbg_sum = sum( w );
    A(dm,dm)=1;
    
    %str = sprintf( 'dm2: %d, %f %f %f %f %f %f ', dm, w(1,1), w(1,2), w(1,3), w(1,4), w(1,5), w(1,6) );
    %disp( str );
    %assert( abs(dbg_sum) < 1e-8 );
    
    %f( dm ) = bd_con(  i, j, resol, ind2dm, dm2ind );
    f( dm ) = 0;
    %f( dm) = bound_con( i, j, resol, ind2dm, dm2ind );
    
    %A( dm, : ) = -A( dm, : );
    %f( dm ) = 1;
     
    retA = A;
    retf = f;
    
end
