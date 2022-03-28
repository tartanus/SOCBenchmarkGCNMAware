function [x,fval,xHist,fHist,exitflag,output] = NM4_s3d6k1n8(fun,x,options,varargin)
% NM4 is (basically) Matlab's FMINSEARCH with frame based
% convergence algorithm added to it and the NM step moved
% to a separate function.
%
% Other functions inside this one: eval_func, vsort
%
% Other functions called by this one:
% simplex
% frame
% NM_std
%
% define global variables so the NM step can see which function
% to use and update the number of function evaluations
global funfcn func_evals counter
funfcn= fun;
counter=1;
% set up Matlab's fminsearch environment
if nargin<3, options = []; end


n = prod(size(x));
numberOfVariables = n;

defaultopt = optimset('display' ,'final' ,'maxiter',...
200*numberOfVariables ,'maxfunevals', ...
200*numberOfVariables ,'TolX' ,1e-4,'TolFun' ,1e-4);

options = optimset(defaultopt,options);
printtype = optimget(options, 'display');
tolx = optimget(options,'tolx');
tolf = optimget(options, 'tolfun');
maxfun = optimget(options,'maxfuneval');
maxiter = optimget(options,'maxiter');
switch printtype
case {'none' ,'off'}
    prnt = 0;
case 'iter'
    prnt = 2;
case 'final'
    prnt = 1;
case 'simplex'
    prnt = 3;
otherwise
    prnt = 1;
end

header= 'Iteration Func-count min f(x)';
header= [header blanks(15) 'QMF       '];
header= [header 'Det Max length Status'] ;

x=x(:);
rho = 1; chi = 2; psi=0.5; sigma=0.5;
onesn= ones(1,n);
two2np1= 2 :n+1;
two2np2= 2:n+2;
one2n= 1:n;
one2np1= 1:n+1;
% set up parameters for new algorithm
h = 1; K = 1.0e+03; N = 100;
nu = 4.50;      % shrink factor for epsilon
kappa=0.25;     % shrink factor for h
delta=1e-9;  % min for det. before collapse signalled

% keep track of numbers of each step performed
reflect = 0;
expand = 0;
cont_outside = 0;
cont_inside = 0;
shrink = 0;
frames = 0;
reshapes = 0;
total_qmf = 0;

% setup initial simplex using Matlab's fminsearch's method
v = simplex('mat' ,x);
fv = eval_func(v);
func_evals = n+1;
mu = prod(frame('1' ,v));
[v fv] = vsort(v,fv);
N = (fv(end) - fv(1))/(N * n * h^nu);
epsilon = N * h^nu;
itercount = 1;
how = 'initial';
if prnt == 2
    disp(' ')
    disp(header)
    disp([sprintf(' %5.0f %10.0f %15.6g %s', ...
    itercount, func_evals, fv(l)), how])
elseif prnt == 3
    clc
    formatsave = get(O,{'format' ,'formatspacing'});
    format compact
    format short e
    disp(' ')
    disp(how)
    v
    fv
    func evals
end

exitflag = 1;

% Main algorithm
% Iterate until the diameter of the simplex is less than
% tolx AND the function values differ from the min by less
% than tolf, or the max function evaluations are exceeded.
% (Cannot use OR instead of AND.)
while func_evals < maxfun && itercount < maxiter
    if max(max(abs(v(: ,two2np1)-v(: ,onesn)))) <= tolx && ...
            max(abs(fv(1)-fv(two2np1))) <= tolf
        break
    end
    basis_status = 'ok';
    more_info = 0;
    qmf = 0;

    % generate new (sorted) simplex by NM step
%     [vnm fnm how] = nm_std(v, fv);

    
    xbar = sum(v(:,one2n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,end);
    x(:) = xr; 
    fxr = funfcn(x,varargin{:});
    func_evals = func_evals+1;
    
    if fxr < fv(:,1)
        % Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
        x(:) = xe; fxe = funfcn(x,varargin{:});
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,end) = xe;
            fv(:,end) = fxe;
            how = 'expand';
        else
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        if fxr < fv(:,n)
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(:,end)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                x(:) = xc; fxc = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(:,end) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,end);
                x(:) = xcc; fxcc = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxcc < fv(:,end)
                    v(:,end) = xcc;
                    fv(:,end) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j=two2np1
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    x(:) = v(:,j); fv(:,j) = funfcn(x,varargin{:});
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fv,j] = sort(fv);
    v = v(:,j);
    itercount = itercount + 1;
    if prnt == 3
        fprintf(' %5.0f        %5.0f     %12.6g         %s', itercount, func_evals, fv(1), how)
    elseif prnt == 4
        disp(' ')
        disp(how)
        v
        fv
        func_evals
    end
    
    fnm=fv;
    
    
    
    
    
    
    
    % new algorithm begins
    if fv(: ,end) - fnm(: ,end) >= epsilon
    % NM made sufficient progress, accept simplex
    v = vnm;
    fv = fnm;
    switch how
        case 'reflect'
            reflect = reflect + 1;
        case 'expand'
            expand = expand + 1;
            mu = mu * chi;
        case 'contract outside'
            cont_outside = cont_outside + 1;
            mu = mu * psi;
        case 'contract inside'
            cont_inside = cont_inside + 1;
            mu = mu * psi;
        case 'shrink'
            shrink = shrink + 1;
            mu = mu * sigma^n;
    end

    % given the simplex, complete the frame
    v_frame = frame('f' ,v);
    fv_frame = [fveval_func(v_frame(:,end))];
    func_evals = func_evals + 1;


    % until there is epsilon descent, reduce the frame
    while max([fv_frame(1) - fv_frame(two2np2)]) <= epsilon
        % have we been in the qmf too long?
        if func_evals > maxfun || itercount > maxiter
            basis_status = 'stuck in QMF' ;
            exi tflag = 0;
            break
        end
        if reshape_flag == 0
            % reshape the simplex
            v = simplex('qr' ,v_frame(: ,1:n+1));
            fv(two2np1) = eval_func(v(: ,two2np1));
            func_evals = func_evals + n;
            mu = prod(frame('l' ,v));
            % complete the frame
            v_frame = frame('f' ,v);
            fv_frame = [fv eval_func(v_frame(: ,end))];
            func_evals = func_evals + 1;
            reshape_flag = 1;
            reshapes = reshapes + 1;
            how= [how' + reshape'];
        else
            qmf = qmf + 1;
            kappa = -kappa;
            h = abs(h * kappa);
            epsilon = N * h^nu;
        end
        % shrink the frame
        v_frame = frame('sh', v_frame, kappa);
        mu = mu * abs(kappa)^n;
        temp = v_frame ( : ,two2np1) - v_frame (: ,onesn)
        simplex_diam = max(max(abs(temp)));
        % has the frame collapsed?
        if simplex_diam == 0
            fv_frame(two2np2) = fv_frame(l);
            basis_status = [basis_status' + frame zero'] ;
            break
        else
            fv_frame(two2np2) = eval_func(v_frame(: ,two2np2));
            func evals func evals + n+1;
        end
        % has required tolerance been reached?
        if simplex_diam <= tolx && ...
            max(abs(fv_frame(l) - fv_frame(two2np1))) <= tolf
            break
        end
    end



    total_qmf=total_qmf+qmf;

    % reduce frame to simplex again
    v= v_frame(: ,one2np1);
    fv= fv_frame(: ,one2np1);
    % check if new frame point is the lowest
    if fv_frame(end) < fv(l)
        v(: ,1) = v_frame(: ,end);
        fv(l) = fv_frame(end);
        if strcmp(how, 'reshape')
            mu=mu * (n+1);
        else
            mu = mu * chi;

        end
        how=[how' + swap'];

    end
% 
%     fvHist(:,counter)=fv;
%     vHist(:,counter)=v;
    fv
    v

        [v fv] = vsort(v, fv);
        % original algorithm continues on from here
        itercount = itercount + 1;
            if prnt == 2
                how = [how blanks(22 - length(how))];
                if more_info
                    disp([sprintf('%5.0f %12.0f %16.6g %25s %3.0f %12.4g %9.3g %s', ...
                    itercount, func_evals, fv(1), how, qmf,...
                    basis_det, max_length, basis_status) ])
                else
                    disp([sprintf(' %5.0f %5.0f %12.6g %18s',...
                    itercount, func_evals, fv(1), how) ])
                end
            elseif prnt == 3
                disp(' ')
                disp(how)
                v
                fv
                func_evals

            end
            
        
    end % while


    x=v(:,1);
    
    xHist(:,counter)=x;
    fHist(counter)=min(fv);
    counter=counter+1;
    if prnt == 3
        % reset format
        set(O,{'format' ,'formatspacing'},formatsave);
    end
    output.iterations = itercount;
    output.funcCount = func_evals;
    output.algorithm = 'Frame based convergence algorithm';

    fval = min(fv);
    if func_evals >= maxfun
        if prnt > 0
        disp(' ')
        disp('Exiting: ')
        %18s', ...
        disp('Maximum number of function evaluations has been exceeded')
        disp(' - increase MaxFunEvals option.')
        msg = sprintf('Current function value: %f \n', fval);
        disp(msg)
        end

        exitflag = 0;

    elseif itercount >= maxiter
        if prnt > 0
            disp(' ')
            disp('Exiting: Maximum number of iterations has been exceeded')
            disp(' - increase MaxIter option.')
            msg = sprintf('                  Current function value: %f \n', fval);
            fprintf(msg)
        end
        exitflag = 0;
    else
        if prnt > 0
            convmsg1 = sprintf([ ...
            '\nOptimization terminated successfully:\n', ...
            'the current x satisfies the termination criteria\n',...
            'using OPTIONS.TolX of %e \n',...
            'and F(X) satisfies the convergence criteria\n',...
            'using OPTIONS.TolFun of %e \n'],...
            options.TolX, options.TolFun);
            disp(convmsg1)
            exitflag = 1;
        end
  end    
end

output.Reflect = reflect;
output.Expand = expand;
output.Outside = cont_outside;
output.Inside = cont_inside;
output.Frames = frames;
output.Total_QMF = total_qmf;
output.Reshapes = reshapes;










