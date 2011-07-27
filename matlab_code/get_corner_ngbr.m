function ret = get_corner_ngbr( i, j, resol )
    ret( 1 ) = ij2ind( i+1, j+1, resol );
    ret( 2 ) = ij2ind( i-1, j+1, resol );
    ret( 3 ) = ij2ind( i-1, j-1, resol );
    ret( 4 ) = ij2ind( i+1, j-1, resol );
end