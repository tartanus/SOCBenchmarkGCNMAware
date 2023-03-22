function [newv,newfv] = vsort(v, fv)
    % sort so v(1,:) has the lowest function value
    [newfv,j] = sort(fv);
    newv = v ( : , j) ;
end