%% Ellipsoidal Ordering
%% Solving Dual Problem 
X = X1;
[n, p] = size(X);
volumes = [];

U = dual_function(X, p);
[A, b] = mve(X);
volumes = [ volumes 1/det(A) ]; 
plotEllipsoid(A, b, X, 1)

%% Obtain Mahalanobis Distance
mhb_dist = []; 

for i=1:size(X, 2)
    xk = X(:, i);
    d = xk' * inv(X*U*X') * xk;
    mhb_dist = [mhb_dist d];
end

[val, idx] = sort(mhb_dist);

%% Select points to keep based on Mahalanobis distance

h = 8;
keep = p - h;

X_valid = [];
for i=1:keep   
    X_valid = [X_valid X(:, idx(i))];
end

removed = [];
for i=keep+1:p
    removed = [removed X(:, idx(i))];
end

[n, p] = size(X_valid);
U = dual_function(X_valid, p);
[A, b] = mve(X_valid);
volumes = [ volumes 1/det(A) ]; 
plotEllipsoid(A, b, X_valid, 1)



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


function [A, b] = mve(X)
    cvx_begin
        variable A(2,2) symmetric
        variable b(2)
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
    title('Ellipsoid 1')

end


