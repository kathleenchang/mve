%% Ellipsoid Peeling
% Boyd, S., & Vandenberghe, L. (2022). Additional Exercises for Convex Optimization. 
h = 8; % number of points to peel
X = X3;
[n, p] = size(X);

volumes = [];
removed = [];

figure 
for i = 1:h

    % fit ellipsoid
    cvx_begin
        variable A(n,n) symmetric 
        variable b(n)
        dual variable v
        maximize (det_rootn(A))
        subject to
            v : norms(A*X + b*ones(1,size(X,2))) <= 1
    cvx_end


    % detect outliers
    [vm idx] = max(v);
    removed = [ removed X(:,idx) ];
    X(:,idx) = [];
    volumes = [ volumes 1/det(A) ]; 

    % plot ellipsoid
    noangles = 200;
    angles   = linspace( 0, 2 * pi, noangles );
    ellipse  = A \ [ cos(angles) - b(1) ; sin(angles) - b(2) ];
    plot( X(1,:), X(2,:), 'ro', ellipse(1,:), ellipse(2,:), 'b-' );
    hold on


end

title('Ellipsoid 1')

% show points that were removed
figure  
plot(X(1,:), X(2,:), 'bx'); % normal points
hold on
plot(removed(1,:), removed(2,:), 'ro'); % outliers   
legend('original points', 'removed') 
title('Ellipsoid 1')

% volume of the ellipsoid as points are removed
figure
hold on
semilogy(volumes)
xlabel('number of removed points');
ylabel('ellipsoid volume');
title('Ellipsoid 1 Volume')
set(gca,'XTick', 1:length(volumes));