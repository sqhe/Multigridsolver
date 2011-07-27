function error = est_error( A, x, b )
    bnrm2 = norm( b );
    r = b - A*x;
    error = norm( r ) / bnrm2;
end