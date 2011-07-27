function ret = get_indirect_ngbr( i, j, resol )
    ret( 1 ) = ij2ind( i+2, j, resol );
    ret( 2 ) = ij2ind( i-2, j, resol );
    ret( 3 ) = ij2ind( i, j+2, resol );
    ret( 4 ) = ij2ind( i, j-2, resol );
end