%% Chorin's algorithm for simple test case
function [u1,u2] = ACM(N)
% parameters:
% N: grid size
%
% returns:
% a steady state velocity field

% create grid
f = 1;
R = 1;                          % Reynolds number
delta = .006*f^2;               % artificial compressibility
%N   = 50; 
dx1 = 1/N;
dx2 = 1/N;
dt  = .6*dx1*R*sqrt(delta);     % according to stability condition


a = 0;
b = 1;
c = 0;
d = 1;

x1 = a:dx1:b;
x2 = c:dx2:d;

[X1,X2] = meshgrid(x1,x2);

% initialize solution (u1,u2,rho). pressure = rho/delta
u1 = zeros([size(X1),3]);
u2 = zeros([size(X1),3]);       % velocity field

rho = zeros([size(X1), 3]);     % fluid density

% must initialize u(t=2dt) since this is a two-step algorithm. just
% satisify boundary conditions for now, maybe can optimize this for
% faster convergence later
b1         = zeros(size(X1));
b1(1,:)    = 4*x1.*(1-x1);    % just the long-time solution.
b1(N+1,:)  = b1(1,:);
b2         = zeros(size(X2));
rho0       = zeros(size(X1));

u1(:,:,1) = b1;
u1(:,:,2) = b1;

u2(:,:,1) = b2;
u2(:,:,2) = b2;

rho(:,:,1) = rho0;
rho(:,:,2) = rho0;

w = (1+2*dt*(1/(dx1)^2 + 1/(dx2)^2))^-1;
t = 3; % this is stupid

s = 1;

err = 10^6;

while err > 10^-5   
    if mod(s,10000) == 0
        s
    end
    
    for i = 1:N+1
        %% bdy updates: left + right boundaries
        if i == 1
             u1(1,:,t)    = b1(1,:);
             u2(1,:,t)    = b2(1,:);
             rho(1,1,t)   = -2*dt/dx1*(u1(2,1,t-1)-u1(1,1,t-1)) ...
                            - 2*dt/dx2*(u2(1,2,t-1)-u2(1,1,t-1)) ...
                            + rho(1,1,t-2);
             rho(1,2:N,t) = -2*dt/dx1*(u1(2,2:N,t-1)-u1(1,2:N,t-1)) ...
                             - dt/dx2*(u2(1,3:N+1,t-1)-u2(1,1:N-1,t-1)) ...
                             + rho(1,2:N,t-2);
             rho(1,N+1,t) = -2*dt/dx1*(u1(2,N+1,t-1)-u1(1,N+1,t-1)) ...
                            - 2*dt/dx2*(u2(1,N+1,t-1)-u2(1,N,t-1)) ...
                            + rho(1,N+1,t-2);
        elseif i == N+1
             u1(N+1,:,t)    = b1(N+1,:);
             u2(N+1,:,t)    = b2(N+1,:);
             rho(N+1,1,t)   = -2*dt/dx1*(u1(N+1,1,t-1)-u1(N,1,t-1)) ...
                              - 2*dt/dx2*(u2(N+1,2,t-1)-u2(N+1,1,t-1)) ...
                              + rho(N+1,1,t-2);
             rho(N+1,2:N,t) = -2*dt/dx1*(u1(N+1,2:N,t-1)-u1(N,2:N,t-1)) ...
                              - dt/dx2*(u2(N+1,3:N+1,t-1)-u2(N+1,1:N-1,t-1)) ...
                              + rho(1,2:N,t-2);
             rho(N+1,N+1,t) = -2*dt/dx1*(u1(N+1,N+1,t-1)-u1(N,N+1,t-1)) ...
                              -2*dt/dx2*(u2(N+1,N+1,t-1)-u2(N+1,N,t-1)) ...
                              + rho(N+1,N+1,t-2);
        else     
            for j = 1:N+1
                %% bdy updates: these should be split and cleaned up
                % for u, simply enforce boundary conditions of steady-state
                % problem
                % for rho, use first-order one-sided difference scheme. this
                % should not affect overall stability, but it is possible to
                % use second order if you're having troubles near the boundary
                if j == 1
                    u1(:,1,t)    = b1(:,1);
                    u2(:,1,t)    = b2(:,1);
                    rho(2:N,1,t) = -2*dt/dx2*(u2(2:N,2,t-1)-u2(2:N,1,t-1)) ...
                                   - dt/dx1*(u1(3:N+1,1,t-1)-u1(1:N-1,1,t-1)) ...
                                   + rho(2:N,1,t-2);
                elseif j == N+1
                    u1(:,N+1,t)    = b1(:,N+1);
                    u2(:,N+1,t)    = b2(:,N+1);
                    rho(2:N,N+1,t) = -2*dt/dx2*(u2(2:N,N+1,t-1)-u2(2:N,N,t-1)) ...
                                     - dt/dx1*(u1(3:N+1,N+1,t-1)-u1(1:N-1,N+1,t-1)) ...
                                     + rho(2:N,N+1,t-2); 
                %% interior updates
                else 
                    u1(i,j,t) = w*(-R*dt/dx1*(u1(i+1,j,t-1)^2-u1(i-1,j,t-1)^2) ...
                                   - R*dt/dx2*(u1(i,j+1,t-1)*u2(i,j+1,t-1) - u1(i,j-1,t-1)*u2(i,j-1,t-1)) ...
                                   + 2*dt/(dx1^2)*(u1(i+1,j,t-1)+u1(i-1,j,t-1) - u1(i,j,t-2)) ...
                                   + 2*dt/(dx2^2)*(u1(i,j+1,t-1)+u1(i,j-1,t-1) - u1(i,j,t-2)) ...
                                   - dt/dx1*1/delta*(rho(i+1,j,t-1) - rho(i-1,j,t-1)) ...
                                   + u1(i,j,t-2));
                    u2(i,j,t) = w*(-R*dt/dx1*(u2(i+1,j,t-1)*u1(i+1,j,t-1) - u2(i-1,j,t-1)*u1(i-1,j,t-1)) ...
                                   - R*dt/dx2*(u2(i,j+1,t-1)^2 - u2(i,j-1,t-1)^2) ...
                                   + 2*dt/(dx1^2)*(u2(i+1,j,t-1) + u2(i-1,j,t-1) - u2(i,j,t-2)) ...
                                   + 2*dt/(dx2^2)*(u2(i,j+1,t-1) + u2(i,j-1,t-1) - u2(i,j,t-2)) ...
                                   - dt/dx2*1/delta*(rho(i,j+1,t-1) - rho(i,j-1,t-1)) ...
                                   + u2(i,j,t-2));

                    rho(i,j,t) = -dt/dx1*(u1(i+1,j,t-1)-u1(i-1,j,t-1)) ...
                                 -dt/dx2*(u2(i,j+1,t-1)-u2(i,j-1,t-1)) ...
                                 + rho(i,j,t-2);
                end
            end
        end
    end
    
    % update values for next iteration
    u1(:,:,1) = u1(:,:,2);
    u1(:,:,2) = u1(:,:,3);
    
    u2(:,:,1) = u2(:,:,2);
    u2(:,:,2) = u2(:,:,3);
    
    rho(:,:,1) = rho(:,:,2);
    rho(:,:,2) = rho(:,:,3);
    
    s = s+1;
    err = max(max(abs(u1(:,:,3)-u1(:,:,1)))) + ...
          max(max(abs(u2(:,:,3)-u2(:,:,1))));

end

end
           


%surf(X1,X2,u1(:,:,3))


