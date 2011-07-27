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
