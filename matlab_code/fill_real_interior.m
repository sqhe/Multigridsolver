function [ retA, retf ] = fill_real_interior( i, j, resol, A, f, lambda, ind2dm, dm2ind )

    ind = ij2ind( i, j, resol );
    dm = ind2dm( ind );
    h = 1 / resol;
    scale = get_interior_scale( resol );
    
    %----------------------------------------------
    % fill laplace^2
    %A( dm, dm ) = A( dm, dm ) + 20 / h^4;
    A( dm, dm ) = A( dm, dm ) + 4 * scale;
   
    ngbr = get_direct_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) - 8 / h^4;
    A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) - 1*scale;
    
    %ngbr = get_corner_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 2 / h^4;
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 2 * scale;
    
    %ngbr = get_indirect_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 / h^4;
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 * scale;
    
    
    %----------------------------------------------
    %fill laplace
    %A( dm, dm ) = A( dm, dm ) - 4 / h^2 * ( -2 * lambda );
    %A( dm, dm ) = A( dm, dm ) - 4 * h^2 * ( -2 * lambda ) * scale;

    %ngbr = get_direct_ngbr( i, j, resol );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 / h^2 * ( -2 * lambda );
    %A( dm, ind2dm(ngbr) ) = A( dm, ind2dm(ngbr) ) + 1 * h^2 * ( -2 * lambda ) * scale;
    
    %----------------------------------------------
    %fill
    %A( dm, dm ) = A( dm, dm ) + lambda * lambda;
    %A( dm, dm ) = A( dm, dm ) + lambda * lambda * h^4 * scale;
   
    f( dm ) = bd_con(  i, j, resol, ind2dm, dm2ind );
    
    retA = A;
    retf = f;
    
end
