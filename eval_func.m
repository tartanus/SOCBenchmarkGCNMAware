function fv = eval_func(v)
global funfcn
    n = size(v,2);
    fv = zeros(1,n);
    % evaluate function values
    for j = 1 :n
    fv(j) = feval(funfcn,v(: ,j));
    end
end