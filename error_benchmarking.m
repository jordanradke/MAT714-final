%% error testing for ACM.m against analytic solution

pp = 12;
err = zeros(pp,1);
h   = zeros(pp,1);

for p = 1:pp
    p
    % grid refinement + time steps
    R = 1;                          % Reynolds number
    N   = 2^p; 
    dx1 = 1/N;
    dx2 = 1/N;

    a = 0;
    b = 1;
    c = 0;
    d = 1;

    x1 = a:dx1:b;
    x2 = c:dx2:d;

    [X1,X2] = meshgrid(x1,x2);
    
    % numerical solution
    fprintf('numerical solution...')
    [u1,u2] = ACM(N);
    
    % analytic solution sampled on grid (only for f=0 init. data)
    fprintf('analytic solution...')
    u1_sol = 4*X1.*(1-X1);    % just the long-time solution.
    u2_sol = zeros(size(X1));

    
    err(p) = max(max(abs(u1(:,:,3)-u1_sol))) + ...
             max(max(abs(u2(:,:,3)-u2_sol)));
    h(p) = dx2;
    
end

% plot loglog error 
figure(1); clf();
loglog(h,err,'o-', 'LineWidth', 2)
hold on;
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error','FontSize',24);
xlabel('$h$','Interpreter','latex','FontSize',24)
ylabel('relative $\ell^\infty$ error', 'Interpreter','latex','FontSize',24)
    
    

    
    