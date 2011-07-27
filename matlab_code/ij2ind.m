function ret = ij2ind( i, j, resol )
    ind = i * (resol+1) + j + 1;
    ret = ind;
    
%     n = (resol+1)^2;
%     ret = n - ind + 1;
end
