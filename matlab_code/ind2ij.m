function [ i, j ] = ind2ij( ind, resol )
    
%     n = (resol+1)^2;
%     ind = n - ind + 1;

    j = mod( ind-1, resol+1 );
    i = ( ind - 1 - j ) / ( resol + 1 );
end