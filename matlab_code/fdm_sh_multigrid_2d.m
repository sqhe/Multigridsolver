function ret = fdm_sh_multigrid_2d( base_resol, depth, lambda )

    resol = base_resol * 2^(depth-1);
    resol
    
    %----------------------------------------------
    % construct A, f
    [ A, f, n, ind2dm, dm2ind ] = construct_Af( resol, lambda );
    [ DM_IND, IND_DM ] = construct_dm_ind( dm2ind, resol, n );
  
    %----------------------------------------------
    % initial guess
    v = zeros( n, 1 );
    
    [ A_base, f_base, n_base, ind2dm_base, dm2ind_base ] = construct_Af( base_resol, lambda );
    [ DM_IND_base, IND_DM_base ] = construct_dm_ind( dm2ind_base, base_resol, n_base );
    vbase = A_base \ f_base;
    
    % interpolate the coarse initial guess to the finest grid
    for i = 1 : depth - 1
        [ A_base, na0, n_base, ind2dm_base, dm2ind_base ] = construct_Af( base_resol * 2^i, lambda );
        I_2h_h = construct_I_2h_h( base_resol * 2^i );
        vbase = I_2h_h * DM_IND_base * vbase;
        [ DM_IND_base, IND_DM_base ] = construct_dm_ind( dm2ind_base, base_resol * 2^i, n_base );
        vbase = IND_DM_base * vbase;
    end
    
    % get the initial guess: v
    v = vbase;   
   
    % draw the initial guess on the finest grid
    sol_vbase = DM_IND_base * vbase;
    Surf_guess = func_to_surf( sol_vbase, resol );
    figure( 1 );
    surf( Surf_guess );
                     
    %----------------------------------------------
    % get ideal result by the direct method for comparison
       
    u = A \ f;
    sol_u = DM_IND * u; 
 
    %----------------------------------------------
    % do multigrid
    for i = 1 : 6
        i
        v = multigrid( A, v, f, DM_IND, IND_DM, ind2dm, dm2ind, resol, depth, lambda );
                
        error_curve( i ) = norm( u - v );
    end
   
    num_ite = i
    sol_v = DM_IND * v;
    
    figure( 2 );
    plot( error_curve );
       
    %----------------------------------------------
    % draw the output
    
    Surf_u = func_to_surf( sol_u, resol ); % ideal result
    Surf_v = func_to_surf( sol_v, resol ); % result by multigrid
        
    figure( 3 );
    surf( Surf_v );
    figure( 4 );
    surf( Surf_u );
    
    %----------------------------------------------
    % error statistics
        
    disp( '----------for mg and u' );
    str = sprintf( 'err_norm = %f, err_pc = %f %%, Linf_u = %f, Linf_v = %f, Linf_err = %f', ...
        norm( v - u ), ...
        norm( v - u ) / norm( u )*100, ...
        norm( u, inf ), ...
        norm( v, inf ), ...
        norm( v-u, inf ) );
    disp( str );
            
    return;

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [ A, f, n, ind2dm, dm2ind ] = construct_Af( resol, lambda )

    % A: the linear system
    % f: right hand side
    % n: number of nodes in the domain
    % ind2dm: a mapping for the node from grid ID to domain ID
    % dm2ind: inverse mapping of ind2dm

    num_data = ( resol + 1 )^2;
    ind2dm = zeros( 1, num_data );
       
    % get all the nodes in the domain and assign them with domain ID
    dmi = 1;
    for ii = 0 : resol 
        for jj = 0 : resol            
            if( is_real_interior( ii, jj, resol ) || ...
                is_bndy_interior( ii, jj, resol ) || ...
                is_bndy_exterior( ii, jj, resol ) )                
                ind = ij2ind( ii, jj, resol );
                dm2ind( dmi ) = ind;
                ind2dm( ind )= dmi;                                
                dmi = dmi + 1;
            end
        end
    end
   
    A = sparse( dmi-1, dmi-1 );
    f = zeros( dmi-1, 1 );
    n = dmi - 1;
  
    % construct the linear system and right hand side
    for ii = 0 : resol
        for jj = 0 : resol
            
            ind = ij2ind( ii, jj, resol );
            dmi = ind2dm( ind );
             
            if( is_real_interior( ii, jj, resol ) ) % interior
                [ A, f ] = fill_real_interior( ii, jj, resol, A, f, lambda, ind2dm, dm2ind );
            elseif( is_bndy_interior( ii, jj, resol ) ) % second order
                [ A, f ] = fill_bndy_interior( ii, jj, resol, A, f, ind2dm, dm2ind );
            elseif( is_bndy_exterior( ii, jj, resol ) ) % first order
                [ A, f ] = fill_bndy_exterior( ii, jj, resol, A, f, ind2dm, dm2ind );
            end
        end
    end

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = multigrid( A, v, f, DM_IND, IND_DM, ind2dm, dm2ind, resol, depth, lambda )
    
    if( depth ~= 1 )
        blanks = '';
        for i = 1 : depth - 1
            blanks = strcat( blanks, '------' );
        end
    
        u = A \ f;

        %----------------------------------------------
        % relax a1 times on Au = f on initial guess v
        abserr = abs_error( u, v );
        str = sprintf( '%d before relax = %f', depth, abserr );
        str = strcat( blanks, str );
        disp( str );
    
        %[v, error, iter, flag]  = sor( A, v, f, 1, 20, 1e-10 ); %3 steps of S.O.R. is the minimum for this problem?
        v = separate_iter_solver( A, v, f, resol, dm2ind );
   
        abserr = abs_error( u, v );
        str = sprintf( '%d before cg = %f', depth, abserr );
        str = strcat( blanks, str );
        disp( str );
    
        %----------------------------------------------
        % recursive do multigrid
        
        resol_2h = resol / 2;
         
        % construct the coarser matrix
        [ A_2h, f_nouse, n_2h, ind2dm_2h, dm2ind_2h ] = construct_Af( resol_2h, lambda );
         
        % construct I_h_2h, I_2h_h
        [ DM_IND_2h, IND_DM_2h ] = construct_dm_ind( dm2ind_2h, resol_2h, n_2h );
        I_h_2h = construct_I_h_2h( resol, ind2dm );
        I_2h_h = construct_I_2h_h( resol );

        % construct A_2h, v_2h, f_2h
        f_2h = IND_DM_2h * I_h_2h * DM_IND * ( f - A * v );
        v_2h = zeros( n_2h, 1 );
        
        
        Surf11 = func_to_surf( DM_IND*(f-A*v), resol );  
        Surf12 = func_to_surf( I_h_2h*DM_IND*(f-A*v), resol_2h );
                
        % sub-level 2 times
        for i = 1 : 2
            v_2h = multigrid( A_2h, v_2h, f_2h, ...
                DM_IND_2h, IND_DM_2h, ind2dm_2h, dm2ind_2h, ...
                resol_2h, depth - 1, lambda );
        end
   
        % correction
        v = v + IND_DM * I_2h_h * DM_IND_2h * v_2h;
   
        % relax a1 times on Au = f on initial guess v
        abserr = abs_error( u, v );
        str = sprintf( '%d after cg = %f', depth, abserr );
        str = strcat( blanks, str );
        disp( str );
    
        v = separate_iter_solver( A, v, f, resol, dm2ind );
                
        abserr = abs_error( u, v );
        str = sprintf( '%d after relax = %f', depth, abserr );
        str = strcat( blanks, str );
        disp( str );
 
        ret = v;
        
    else
        
        % at the coarses level, solve it by direct method
               
        v = A \ f;
        ret = v;
        
        str = sprintf( '%d after relax = %f', depth, 0 );
        disp( str );
        
    end
  
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [ DM_IND, IND_DM ] = construct_dm_ind( dm2ind, resol, n )

    % DM_IND: map a vector from domain to grid
    % IND_DM: map a vector from grid to domain
   
    DM_IND = sparse( (resol+1)^2, n );
    IND_DM = sparse( n, (resol+1)^2 );
    
    for i = 1 : n
        ind = dm2ind( i );
        DM_IND( ind, i ) = 1;
        IND_DM( i, ind ) = 1;
    end

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [ retA, retf ] = fill_real_interior( i, j, resol, A, f, lambda, ind2dm, dm2ind )

    ind = ij2ind( i, j, resol );
    dm = ind2dm( ind );
    h = 1 / resol;
    scale = get_interior_scale( resol );
    
    %----------------------------------------------
    % fill laplace^2
    %A( dm, dm ) = A( dm, dm ) + 20 / h^4;
    A( dm, dm ) = A( dm, dm ) + 20 * scale;
   
    ngbr = get_direct_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) - 8 / h^4;
    A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) - 8 * scale;
    
    ngbr = get_corner_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 2 / h^4;
    A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 2 * scale;
    
    ngbr = get_indirect_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 / h^4;
    A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 * scale;
    
    
    %----------------------------------------------
    %fill laplace
    %A( dm, dm ) = A( dm, dm ) - 4 / h^2 * ( -2 * lambda );
    A( dm, dm ) = A( dm, dm ) - 4 * h^2 * ( -2 * lambda ) * scale;

    ngbr = get_direct_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 / h^2 * ( -2 * lambda );
    A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 * h^2 * ( -2 * lambda ) * scale;
    
    %----------------------------------------------
    %fill
    %A( dm, dm ) = A( dm, dm ) + lambda * lambda;
    A( dm, dm ) = A( dm, dm ) + lambda * lambda * h^4 * scale;
   
    f( dm ) = 0;
    
    retA = A;
    retf = f;
    
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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
          
    if ( ngbr_x - pt_x ) * n_x > 0
        A( dm, dm_a ) = A( dm, dm_a ) + (-1 ) / h * abs( n_x );
        A( dm, dm_c ) = A( dm, dm_c ) + ( 1 ) / h * abs( n_x );
    else
        A( dm, dm_a ) = A( dm, dm_a ) + ( 1 ) / h * abs( n_x );
        A( dm, dm_c ) = A( dm, dm_c ) + (-1 ) / h * abs( n_x );
    end
            
    if ( ngbr_y - pt_y ) * n_y > 0        
        A( dm, dm_a ) = A( dm, dm_a ) + (-1 ) / h * abs( n_y );
        A( dm, dm_b ) = A( dm, dm_b ) + ( 1 ) / h * abs( n_y );
    else
        A( dm, dm_a ) = A( dm, dm_a ) + ( 1 ) / h * abs( n_y );
        A( dm, dm_b ) = A( dm, dm_b ) + (-1 ) / h * abs( n_y );
    end
    
    assert( abs( A(dm, dm_a) + A(dm, dm_b) + A(dm, dm_c) ) < 1e-6 );
    
    f( dm ) = 1e-10;
    
    retA = A;
    retf = f;
    
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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
        
    A( dm, dm_arr ) = w;
    dbg_sum = sum( w );
    
    str = sprintf( 'dm2: %d, %f %f %f %f %f %f ', dm, w(1,1), w(1,2), w(1,3), w(1,4), w(1,5), w(1,6) );
    %disp( str );
    assert( abs(dbg_sum) < 1e-8 );
    
    f( dm ) = -1;
    
    %A( dm, : ) = -A( dm, : );
    %f( dm ) = 1;
     
    retA = A;
    retf = f;
    
end


%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = ij2ind( i, j, resol )
    ind = i * (resol+1) + j + 1;
    ret = ind;
    
%     n = (resol+1)^2;
%     ret = n - ind + 1;
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [ i, j ] = ind2ij( ind, resol )
    
%     n = (resol+1)^2;
%     ind = n - ind + 1;

    j = mod( ind-1, resol+1 );
    i = ( ind - 1 - j ) / ( resol + 1 );
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = get_direct_ngbr( i, j, resol )
    ret( 1 ) = ij2ind( i+1, j, resol );
    ret( 2 ) = ij2ind( i-1, j, resol );
    ret( 3 ) = ij2ind( i, j+1, resol );
    ret( 4 ) = ij2ind( i, j-1, resol );        
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = get_corner_ngbr( i, j, resol )
    ret( 1 ) = ij2ind( i+1, j+1, resol );
    ret( 2 ) = ij2ind( i-1, j+1, resol );
    ret( 3 ) = ij2ind( i-1, j-1, resol );
    ret( 4 ) = ij2ind( i+1, j-1, resol );
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = get_indirect_ngbr( i, j, resol )
    ret( 1 ) = ij2ind( i+2, j, resol );
    ret( 2 ) = ij2ind( i-2, j, resol );
    ret( 3 ) = ij2ind( i, j+2, resol );
    ret( 4 ) = ij2ind( i, j-2, resol );
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = is_bndy_interior( i, j, resol )

    if( ~is_interior( i, j, resol ) )
        ret = 0;
    else
        
        ngbr = get_direct_ngbr( i, j, resol );
        for ii = 1 : 4
            [ni, nj] = ind2ij( ngbr(ii), resol );
            
            if( ~is_interior( ni, nj, resol ) )
                ret = 1;
                return;
            end
        end
        
        ret = 0;
    end
    
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = is_bndy_exterior( i, j, resol ) 

    if( is_interior( i, j, resol ) )
        ret = 0;
    else
        
        ngbr = get_direct_ngbr( i, j, resol );
        for ii = 1 : 4
            [ni, nj] = ind2ij( ngbr(ii), resol );
            
            if( is_interior( ni, nj, resol ) )
                ret = 1;
                return;
            end
        end
        
        ngbr = get_corner_ngbr( i, j, resol );
        for ii = 1 : 4
            [ni, nj] = ind2ij( ngbr(ii), resol );
            
            if( is_interior( ni, nj, resol ) )
                ret = 1;
                return;
            end
        end
        
        ret = 0;
    end
    
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret_x = separate_iter_solver( A, x, b, resol, dm2ind )
    
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = assert( b )
    if( ~b )
        if( nargin<2 || isempty(s) )
            s = ' '; 
        end
        error( s ); 
    end
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [ M, N, b ] = split( A, b, w, flag )
%
% function [ M, N, b ] = split( A, b, w, flag )
%
% split.m sets up the matrix splitting for the stationary
% iterative methods: jacobi and sor (gauss-seidel when w = 1.0 )
%
% input   A        DOUBLE PRECISION matrix
%         b        DOUBLE PRECISION right hand side vector (for SOR)
%         w        DOUBLE PRECISION relaxation scalar
%         flag     INTEGER flag for method: 1 = jacobi
%                                           2 = sor
%
% output  M        DOUBLE PRECISION matrix
%         N        DOUBLE PRECISION matrix such that A = M - N
%         b        DOUBLE PRECISION rhs vector ( altered for SOR )

  [m,n] = size( A );
       
  if ( flag == 1 ),                   % jacobi splitting

     M = diag(diag(A));
     N = diag(diag(A)) - A;

  elseif ( flag == 2 ),               % sor/gauss-seidel splitting

     b = w * b;
     M =  w * tril( A, -1 ) + diag(diag( A ));
     N = -w * triu( A,  1 ) + ( 1.0 - w ) * diag(diag( A ));

  end

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [x, error, iter, flag]  = sor(A, x, b, w, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag]  = sor(A, x, b, w, max_it, tol)
%
% sor.m solves the linear system Ax=b using the 
% Successive Over-Relaxation Method (Gauss-Seidel method when omega = 1 ).
%
% input   A        REAL matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         w        REAL relaxation scalar
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

  flag = 0;                                   % initialization
  iter = 0;

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  [ M, N, b ] = split( A, b, w, 2 );          % matrix splitting

  for iter = 1:max_it                         % begin iteration

     x_1 = x;
     x   = M \ ( N*x + b );                   % update approximation

     error = norm( x - x_1 ) / norm( x );     % compute error
     if ( error <= tol ), break, end          % check convergence

  end
  b = b / w;                                  % restore rhs

  if ( error > tol ) flag = 1; end;           % no convergence

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function error = est_error( A, x, b )
    bnrm2 = norm( b );
    r = b - A*x;
    error = norm( r ) / bnrm2;
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function error = abs_error( u, v )
    error = norm( u - v );
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [ c_x, c_y, r, rot ] = get_domain()
    % c_x, c_y: center of square / disk
    % r: radius of disk / half of the edge length of a square
    % rot: rotate the domain by rot degree
    
    c_x = 0.5;
    c_y = 0.5;
    r = 0.30001; 
    rot = 0 / 180 * pi;
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = is_domain_disk()
    % if the domain is disk?
    ret = 1;
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
function ret = get_first_scale( resol )
    ret = 1;
end

function ret = get_interior_scale( resol )
    ret = resol^4;
end
