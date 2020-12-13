%% Chorin's algorithm for simple test case. can generalize as needed.

% create grid
R = f;                          % Reynolds number
delta = .006*f^2;                 % artificial compressibility
N   = 50; 
T   = .2;
dx1 = 1/N;
dx2 = 1/N;
dt  = .6*dx1*R*sqrt(delta);   % what does stability analysis say we need here?


a = 0;
b = 1;
c = 0;
d = 1;

x1 = a:dx1:b;
x2 = c:dx2:d;

[X2,X1] = meshgrid(x2,x1); s = size(X2);
M = s(1)-1;
N = s(2)-1;

% initialize solution (u,v,rho). pressure = rho/delta
u = zeros([size(X1),3]);
v = zeros([size(X1),3]);       % velocity field

rho = zeros([size(X1), 3]);     % fluid density

% must initialize u(t=2dt) since this is a two-step algorithm. just
% satisify boundary conditions for now, maybe can optimize this for
% faster convergence later
bdy1         = zeros(size(X1));
bdy1(1,:)    = 4*x2.*(1-x2);    % just the long-time solution.
bdy1(N+1,:)  = bdy1(1,:);
bdy2         = zeros(size(X2));
rho0         = zeros(size(X1));

u(:,:,1) = bdy1;
u(:,:,2) = bdy1;

v(:,:,1) = bdy2;
v(:,:,2) = bdy2;

rho(:,:,1) = rho0;
rho(:,:,2) = rho0;

w = (1+2*dt*(1/(dx1)^2 + 1/(dx2)^2))^-1;
t = 3; % this is stupid

for s = 3:floor(T/dt)
    if mod(s,10000) == 0
        s
    end
    
    for i = 1:N+1
        %% bdy updates: left + right boundaries periodic
        if i == 1
             u(1,:,t)    = bdy1(1,:);
             v(1,:,t)    = bdy2(1,:);
             rho(1,1,t)   = -2*dt/dx1*(u(2,1,t-1)-u(1,1,t-1)) ...
                            - 2*dt/dx2*(v(1,2,t-1)-v(1,1,t-1)) ...
                            + rho(1,1,t-2);
             rho(1,2:N,t) = -2*dt/dx1*(u(2,2:N,t-1)-u(1,2:N,t-1)) ...
                             - dt/dx2*(v(1,3:N+1,t-1)-v(1,1:N-1,t-1)) ...
                             + rho(1,2:N,t-2);
             rho(1,N+1,t) = -2*dt/dx1*(u(2,N+1,t-1)-u(1,N+1,t-1)) ...
                            - 2*dt/dx2*(v(1,N+1,t-1)-v(1,N,t-1)) ...
                            + rho(1,N+1,t-2);
        elseif i == N+1
             u(N+1,:,t)    = bdy1(N+1,:);
             v(N+1,:,t)    = bdy2(N+1,:);
             rho(N+1,1,t)   = -2*dt/dx1*(u(N+1,1,t-1)-u(N,1,t-1)) ...
                              - 2*dt/dx2*(v(N+1,2,t-1)-v(N+1,1,t-1)) ...
                              + rho(N+1,1,t-2);
             rho(N+1,2:N,t) = -2*dt/dx1*(u(N+1,2:N,t-1)-u(N,2:N,t-1)) ...
                              - dt/dx2*(v(N+1,3:N+1,t-1)-v(N+1,1:N-1,t-1)) ...
                              + rho(1,2:N,t-2);
             rho(N+1,N+1,t) = -2*dt/dx1*(u(N+1,N+1,t-1)-u(N,N+1,t-1)) ...
                              -2*dt/dx2*(v(N+1,N+1,t-1)-v(N+1,N,t-1)) ...
                              + rho(N+1,N+1,t-2);
        else     
            for j = 1:N+1
                %% bdy updates: top and bottom
                % for u, simply enforce boundary conditions of steady-state
                % problem
                % for rho, use first-order one-sided difference scheme
                if j == 1
                    u(:,1,t)    = bdy1(:,1);
                    v(:,1,t)    = bdy2(:,1);
                    rho(2:N,1,t) = -2*dt/dx2*(v(2:N,2,t-1)-v(2:N,1,t-1)) ...
                                   - dt/dx1*(u(3:N+1,1,t-1)-u(1:N-1,1,t-1)) ...
                                   + rho(2:N,1,t-2);
                elseif j == N+1
                    u(:,N+1,t)    = bdy1(:,N+1);
                    v(:,N+1,t)    = bdy2(:,N+1);
                    rho(2:N,N+1,t) = -2*dt/dx2*(v(2:N,N+1,t-1)-v(2:N,N,t-1)) ...
                                     - dt/dx1*(u(3:N+1,N+1,t-1)-u(1:N-1,N+1,t-1)) ...
                                     + rho(2:N,N+1,t-2); 
                %% interior updates
                else 
                    u(i,j,t) = w*(-R*dt/dx1*(u(i+1,j,t-1)^2-u(i-1,j,t-1)^2) ...
                                   - R*dt/dx2*(u(i,j+1,t-1)*v(i,j+1,t-1) - u(i,j-1,t-1)*v(i,j-1,t-1)) ...
                                   + 2*dt/(dx1^2)*(u(i+1,j,t-1)+u(i-1,j,t-1) - u(i,j,t-2)) ...
                                   + 2*dt/(dx2^2)*(u(i,j+1,t-1)+u(i,j-1,t-1) - u(i,j,t-2)) ...
                                   - dt/dx1*1/delta*(rho(i+1,j,t-1) - rho(i-1,j,t-1)) ...
                                   + u(i,j,t-2));
                    v(i,j,t) = w*(-R*dt/dx1*(v(i+1,j,t-1)*u(i+1,j,t-1) - v(i-1,j,t-1)*u(i-1,j,t-1)) ...
                                   - R*dt/dx2*(v(i,j+1,t-1)^2 - v(i,j-1,t-1)^2) ...
                                   + 2*dt/(dx1^2)*(v(i+1,j,t-1) + v(i-1,j,t-1) - v(i,j,t-2)) ...
                                   + 2*dt/(dx2^2)*(v(i,j+1,t-1) + v(i,j-1,t-1) - v(i,j,t-2)) ...
                                   - dt/dx2*1/delta*(rho(i,j+1,t-1) - rho(i,j-1,t-1)) ...
                                   + v(i,j,t-2));

                    rho(i,j,t) = -dt/dx1*(u(i+1,j,t-1)-u(i-1,j,t-1)) ...
                                 -dt/dx2*(v(i,j+1,t-1)-v(i,j-1,t-1)) ...
                                 + rho(i,j,t-2);
                end
            end
        end
    end
    
    % update values for next iteration
    u(:,:,1) = u(:,:,2);
    u(:,:,2) = u(:,:,3);
    
    v(:,:,1) = v(:,:,2);
    v(:,:,2) = v(:,:,3);
    
    rho(:,:,1) = rho(:,:,2);
    rho(:,:,2) = rho(:,:,3);

end
           


surf(X1,X2,u(:,:,3))


