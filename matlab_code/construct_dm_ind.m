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
