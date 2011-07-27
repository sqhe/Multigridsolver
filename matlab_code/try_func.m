function [tsm] = try_func( A, x, b, resol, dm2ind )
    
    %ret_x = sor( A, x, b, 1, 5, 1e-10 );
    %return;

    %ret_x = A \ b;
    %return;

    [ m,n ] = size( A );
    
    in_id = 1;
    bd_id = 1;
    for dmi = 1 : n
        ind = dm2ind( dmi );
        [ i,j ] = ind2ij( ind, resol );
        
        if is_real_interior( i, j, resol )
            dm2in( dmi ) = in_id;
            in2dm( in_id ) = dmi;
            in_id = in_id + 1;
        else
            dm2bd( dmi ) = bd_id;
            bd2dm( bd_id ) = dmi;
            bd_id = bd_id + 1;
        end
    end
    
    n_in = in_id - 1;
    n_bd = bd_id - 1;

    % combined to separate
    DM_SEP = sparse( n, n );
    for dmi = 1 : n
        ind = dm2ind( dmi );
        [ i,j ] = ind2ij( ind, resol );
        
        if is_real_interior( i, j, resol )
            DM_SEP( dm2in( dmi ), dmi ) = 1;
        else
            DM_SEP( n_in + dm2bd(dmi), dmi ) = 1;
        end
        
    end
    
    tsm=zeros(n,n);
    
    for i=1:n
        for j=1:n
            tsm(i,j)=DM_SEP(i,j);
        end
    end
    
    SEP_DM = DM_SEP';
    assert( all( all( SEP_DM*DM_SEP == eye(n) ) ) );
        
    x = DM_SEP * x;
    A = DM_SEP * A * SEP_DM;
    b = DM_SEP * b;
    
    II = A( 1:n_in, 1:n_in );
    IB = A( 1:n_in, n_in+1:n );
    BI = A( n_in+1:n, 1:n_in );
    BB = A( n_in+1:n, n_in+1:n );
    
%     if all( eig(II) > 0 )
%         spd = 1;
%     else
%         assert( 0 );
%     end
    
    assert( all( all( II - II' == zeros(n_in) ) ) );
    
    bi = b( 1:n_in );
    bb = b( n_in+1:n );
    
    xi = x( 1:n_in );
    xb = x( n_in+1:n );
    
    for i = 1 : 20
        xb = BB \ (bb - BI * xi );
        xi = sor( II, xi, bi-IB*xb, 1, 1, 1e-10 );
    end
    
    %xi = II \ (bi-IB*xb);
    x = [ xi; xb ];
    
    ret_x = SEP_DM * x;

end