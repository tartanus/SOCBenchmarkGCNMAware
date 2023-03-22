function v = simplex(type,x,b)
% SIMPLEX generates a matrix whose columns are the vertices
% of an orthogonal simplex with root-vertex x.
%
% v = simplex(type, x, [b])
%
% Where type is either:
%
% mckinnon = proper starting simplex for mckinnon's example
% mat Matlab's method based on initial vertex x
% qr QR decomposition of basis b about initial vertex x
% ijk orth. reg. simplex about x with side lengths in b
% hh orth. simplex about x which is the
% 'orthogonal complement' of -b

switch type
case 'mckinnon'
    SQR33 = sqrt(33);
    lambda1 = (1 + SQR33)/8;
    lambda2 = (1 - SQR33)/8;
    v = [0 lambda1 1; 0 lambda2 1];
otherwise
    if nargin < 2
    error('Insufficient information supplied')
    end
    % initialise
    n = size(x,1);
    if n == 1
        x=x(:);
        n = length(x);
    end
    v = zeros(n,n+1);
    v(:,1) = x(:,1);
    switch type
        case 'mat'
            usual_delta = 0.05;
            zero_term_delta = 0.00025;
        for j = 1:n
            y = x;
            if y(j)~=0
                y(j)=(1 + usual_delta)*y(j);
            else
                y(j) = zero_term_delta;
            end
                v(:,j+1)=y;
        end
    
        case 'ijk'
            if nargin ~= 3
                error('Side length information has not been supplied'),
            end
            if length(b) == 1
                % create a regular orthogonal simplex
                b = ones(1,n)*b;
                % create an orthogonal simplex using side_lengths
                % check there are no zero lengths
            end
            if any(b == 0)
                error('Side lengths must be non-zero')
            elseif length(b) < n
                error('Not enough side lengths given')
            else
                for j = 1:n
                    y = x;
                    y(j) = y(j) + b(j);
                    v(: ,j+1) = y;
                end
            end
        case 'hh'
            % create simplex about xC: ,1) with orthogonal
            % decomposition from vector -b
            % sum of orthogonal decomposition vectors = -b
            % given simplex, find the longest basis vector
            [basis basis_lengths] = frame('pb' ,x);
            j = find(basis_Iengths(1:n) == max(basis_Iengths(1:n)));
            j = j(1);
            bl = basis(: ,j);
            % calculate the orthogonal decomposition, vectors sum to -b
            basis_orth = hh(bl);
            % create the new simplex
            for j = 1 :n
                v(: ,j+1) = v(: ,1) + basis_orth(: ,j);
            end
       case 'qr'
        % create simplex about xC: ,1) using QR decomposition
        % of basis vectors for simplex x
        % get the basis vectors for the current simplex
        [basis basis_lengths]=frame('pb' ,x);
        % order basis vectors according to length of first
        % n basis vectors
        [sorted_lengths, j] = sort(basis_lengths(1:n));
        % get in descending order
        j = fliplr(j);
        basis = basis(: ,j);
        % find QR decomposition of the ordered basis vectors
        [Q,R] = qr(basis);
        % setup new length criteria
        d = diag(R);
        davg = sum(abs(d)) / n;
        sign_d = sign(0.5 + sign(d));
        d_new = sign_d .* max (abs (d) , davg/10);
        D = diag(d_new);
        % calculate new basis vectors
        basis = Q*D;
        % create new simplex about xC; ,1)
        % new simplex is x(:,1) and x(:,1) + basis(: ,j) for j=1 .. n
        for j = 1:n
            v(:,j+1)=v(:,1) + basis(:,j);
        end
    
        
     otherwise
        error('An unknown simplex type has been used')
     end
end   
    