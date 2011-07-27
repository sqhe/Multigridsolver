function ret = fdm_sh_multigrid_2d_new( base_resol, depth, lambda )

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
    surf(Surf_u );

    
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
