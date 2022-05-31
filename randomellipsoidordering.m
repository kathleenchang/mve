%% Random Ellipsoidal Ordering
%% Solving Dual Problem 
X = X3;
volumes = [];

for counter = 1:50

    X = X3;
    [n, p] = size(X);
    
    % select vectors to span R^n
    randidx = randi([1, length(X)], n, 1);

    while length(unique(randidx)) < length(randidx)
        randidx = randi([1, length(X)], n, 1);
    end

    % store basis vectors
    Xbasis = [];
    for i=1:length(randidx)
        Xbasis = [Xbasis X(:, randidx(i))]; 
    end
    
    
    U = dual_function(X, p);
    [A, b] = mve(X, n);
    volumes = [ volumes 1/det(A) ]; 

    mhb_dist = []; 

    for i=1:size(X, 2)
        xk = X(:, i);
        d = xk' * inv(X*U*X') * xk;
        mhb_dist = [mhb_dist d];

    end

    [val, mhb_idx] = sort(mhb_dist); % sort by Mahabolis distance
    
    k = n;
    h = 8; % number of points to remove
    
    for i = 1:length(mhb_idx)    
        if ~any(ismember(mhb_idx(i), randidx)) % if the current index is not already in the basis        
            k = k + 1;
            Xbasis = [Xbasis X(:, mhb_idx(i))];
        end

        if k == p-h % stop when you have removed h points          
            break
        end
        
    end
 
    [n, p] = size(Xbasis);
    [A, b] = mve(Xbasis, n);
    volumes = [ volumes 1/det(A) ]; % recompute volume 

end

%% Percent reduction in volume histogram

vol = reshape(volumes, [2, 50]);

pdiff = (vol(1,:)-vol(2,:))./vol(1,:);

histogram(pdiff, 5)
title('R^{5}')
xlabel('Percent Reduction')
ylabel('Counts')

mean(pdiff)


%% Functions

function U = dual_function(X, m)
   
    cvx_begin

        variable U(m, m) diagonal;
        variable u(m) 

        maximize (det_rootn(X * U * X'))
        subject to
           diag(U) == u;
           sum(u) == 1;
           u >= 0;

    cvx_end
end


function [A, b] = mve(X, n)
    cvx_begin
        variable A(n, n) symmetric
        variable b(n)
        dual variable v
        maximize (det_rootn(A))
        subject to
            v : norms(A*X + b*ones(1,size(X,2))) <= 1
	cvx_end

end


function [] = plotEllipsoid(A, b, X, figno)
    
    noangles = 200;
    angles   = linspace( 0, 2 * pi, noangles );
    ellipse  = A \ [ cos(angles) - b(1) ; sin(angles) - b(2) ];
    figure(figno)
    hold on
    plot( X(1,:), X(2,:), 'ro', ellipse(1,:), ellipse(2,:), 'b-' );
    title('Ellipsoid 3')

end
