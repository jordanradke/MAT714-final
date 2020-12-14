%% benard's experiment: fluid heated from bottom
function [u,v,temp,rho,dt,x2,x3,dx2,dx3,a0,sigma,s] = benard(f,NN)
% parameters:
% f: the Rayleigh number as a fraction of R_c = 1708
% NN: the grid size
%
% returns:
% steady-state velocity and temperature profile, and various parameters 
% needed to compute Nusselt number in nusselt.m

% physical parameters
T0 = 1;
T1 = 0;

R = 1708*f; % Rayleigh number
sigma = 1;  % Prandtl number


q = 1000*f*(T0 - T1);     % strength of heat flux. really a derived quantity from T0, T1

% create grid
delta = 0.001*f;
%NN   = 40; 
T   = 1;
dx2 = 1/NN;
dx3 = 1/NN;
c = R*(1/(sigma*q*delta))^(1/2);
dt  = .6*dx2/c;   % what does stability analysis say we need here?
TT  = ceil(T/dt);

a0 = 2*pi/3.117;  % period of first unstable mode
a = 0;
b = 1;
c = 0;
d = 2*pi/a0;

x2 = c:dx2:a0; 
x3 = a:dx3:b;

[X3,X2] = meshgrid(x3,x2); s = size(X2);
M = s(1)-1;  % length
N = s(2)-1;  % height

%% initialize//boundary conditions
% initialize solution (u,v,rho,T). pressure = rho/delta
u   = zeros([size(X2),3]);
v   = zeros([size(X2),3]);      % velocity field

rho  = zeros([size(X2), 3]);     % fluid density

temp = zeros([size(X2), 3]);     % temperature

% must initialize u(t=2dt) since this is a two-step algorithm. just
% satisify boundary conditions for now, maybe can optimize this for
% faster convergence later
bdy1         = zeros(size(X2));
%bdy1(1,:)    = 4*x2.*(1-x2);    % just the long-time solution.
%bdy1(N+1,:)  = bdy1(1,:);
bdy2         = zeros(size(X2));

rho0         = zeros(size(X2));

% interpolates from given boundary values
[x3_samp, x2_samp] = meshgrid([a,b], x2);
temp0_samp = [T0*ones(size(x2)).', T1*ones(size(x2)).'];
temp0 = interp2(x3_samp, x2_samp, temp0_samp, X3,X2);

% add random noise to initial distribution
eps = .01;
temp0 = temp0 + eps.*[zeros(M+1,1), rand(M+1,N-1), zeros(M+1,1)]; 

u(:,:,1) = bdy1;
u(:,:,2) = bdy1;
u(:,:,3) = bdy1;

v(:,:,1) = bdy2;
v(:,:,2) = bdy2;
v(:,:,3) = bdy2;

rho(:,:,1) = rho0;
rho(:,:,2) = rho0;
rho(:,:,3) = rho0;

temp(:,:,1) = temp0;
temp(:,:,2) = temp0;
temp(:,:,3) = temp0;


%% main loop
w      = (1+2*dt*(1/(dx2).^2 + 1/(dx3).^2))^-1;
wsigma = (1+2*dt/sigma*(1/(dx2).^2 + 1/(dx3).^2))^-1;
t=3;

s = 1;
err = 10^6;

while err > 10^-5

    if mod(s,10000) == 0
        s
    end
    
    for i = 1:M+1
        %% bdy updates: left + right boundaries
        % these boundaries are periodic 
        if i == 1
             % velocity (momentum conservation)
            u(1,2:N,t) = w*(-dt/dx2*(u(2,2:N,t-1).^2-u(M+1,2:N,t-1).^2) ...
                          - dt/dx3*(u(1,3:N+1,t-1).*v(1,3:N+1,t-1) - u(1,1:N-1,t-1).*v(1,1:N-1,t-1)) ...
                          + 2*dt/(dx2^2)*(u(2,2:N,t-1)+u(M+1,2:N,t-1) - u(1,2:N,t-2)) ...
                          + 2*dt/(dx3^2)*(u(1,3:N+1,t-1)+u(1,1:N-1,t-1) - u(1,2:N,t-2)) ...
                          - dt/dx2*R/(sigma*q*delta)*(rho(2,2:N,t-1) - rho(M+1,2:N,t-1)) ...
                          + u(1,2:N,t-2));
            v(1,2:N,t) = w*(-dt/dx2*(v(2,2:N,t-1).*u(2,2:N,t-1) - v(M+1,2:N,t-1).*u(M+1,2:N,t-1)) ...
                          - dt/dx3*(v(1,3:N+1,t-1).^2 - v(1,1:N-1,t-1).^2) ...
                          + 2*dt/(dx2.^2)*(v(2,2:N,t-1) + v(M+1,2:N,t-1) - v(1,2:N,t-2)) ...
                          + 2*dt/(dx3.^2)*(v(1,3:N+1,t-1) + v(1,1:N-1,t-1) - v(1,2:N,t-2)) ...
                          - 2*dt*R/(sigma*q)*(1-q*(temp(1,2:N,t-1)-1)) ...
                          - dt/dx3*R/(sigma*q*delta)*(rho(1,3:N+1,t-1) - rho(1,1:N-1,t-1)) ...
                          + v(1,2:N,t-2));

            % deMsity (mass conservation)
            rho(1,2:N,t) = -dt/dx2*(u(2,2:N,t-1)-u(M+1,2:N,t-1)) ...
                           -dt/dx3*(v(1,3:N+1,t-1)-v(1,1:N-1,t-1)) ...
                           + rho(1,2:N,t-2);

            % temperature (advection-diffusion)
            temp(1,2:N,t) = wsigma*(-2*dt/dx2*(temp(2,2:N,t-1).*u(2,2:N,t-1)-temp(M+1,2:N,t-1).*u(M+1,2:N,t-1)) ...
                                  -2*dt/dx3*(temp(1,3:N+1,t-1).*v(1,3:N+1,t-1)-temp(1,1:N-1,t-1).*v(1,1:N-1,t-1)) ...
                                  +2*dt/(dx2.^2)*(temp(2,2:N,t-1)+temp(M+1,2:N,t-1)-temp(1,2:N,t-2)) ...
                                  +2*dt/(dx3.^2)*(temp(1,3:N+1,t-1)+temp(1,1:N-1,t-1)-temp(1,2:N,t-2)) ...
                          + temp(1,2:N,t-2));
        elseif i == M+1
             % velocity (momentum conservation)
            u(M+1,2:N,t) = w*(-dt/dx2*(u(1,2:N,t-1).^2-u(M,2:N,t-1).^2) ...
                          - dt/dx3*(u(M+1,3:N+1,t-1).*v(M+1,3:N+1,t-1) - u(M+1,1:N-1,t-1).*v(M+1,1:N-1,t-1)) ...
                          + 2*dt/(dx2.^2)*(u(1,2:N,t-1)+u(M,2:N,t-1) - u(M+1,2:N,t-2)) ...
                          + 2*dt/(dx3.^2)*(u(M+1,3:N+1,t-1)+u(M+1,1:N-1,t-1) - u(M+1,2:N,t-2)) ...
                          - dt/dx2*R/(sigma*q*delta)*(rho(1,2:N,t-1) - rho(M,2:N,t-1)) ...
                          + u(M+1,2:N,t-2));
            v(M+1,2:N,t) = w*(-dt/dx2*(v(1,2:N,t-1).*u(1,2:N,t-1) - v(M,2:N,t-1).*u(M,2:N,t-1)) ...
                          - dt/dx3*(v(M+1,3:N+1,t-1).^2 - v(M+1,1:N-1,t-1).^2) ...
                          + 2*dt/(dx2.^2)*(v(1,2:N,t-1) + v(M,2:N,t-1) - v(M+1,2:N,t-2)) ...
                          + 2*dt/(dx3.^2)*(v(M+1,3:N+1,t-1) + v(M+1,1:N-1,t-1) - v(M+1,2:N,t-2)) ...
                          - 2*dt*R/(sigma*q)*(1-q*(temp(M+1,2:N,t)-1)) ...
                          - dt/dx3*R/(sigma*q*delta)*(rho(M+1,3:N+1,t-1) - rho(M+1,1:N-1,t-1)) ...
                          + v(M+1,2:N,t-2));

            % density (mass conservation)
            rho(M+1,2:N,t) = -dt/dx2*(u(1,2:N,t-1)-u(M,2:N,t-1)) ...
                             -dt/dx3*(v(M+1,3:N+1,t-1)-v(M+1,1:N-1,t-1)) ...
                             + rho(M+1,2:N,t-2);

            % temperature (advection-diffusion)
            temp(M+1,2:N,t) = wsigma*(-2*dt/dx2*(temp(1,2:N,t-1).*u(1,2:N,t-1)-temp(M,2:N,t-1).*u(M,2:N,t-1)) ...
                                  -2*dt/dx3*(temp(M+1,3:N+1,t-1).*v(M+1,3:N+1,t-1)-temp(M+1,1:N-1,t-1).*v(M+1,1:N-1,t-1)) ...
                                  +2*dt/(dx2.^2)*(temp(1,2:N,t-1)+temp(M,2:N,t-1)-temp(M+1,2:N,t-2)) ...
                                  +2*dt/(dx3.^2)*(temp(M+1,3:N+1,t-1)+temp(M+1,1:N-1,t-1)-temp(M+1,2:N,t-2)) ...
                            + temp(M+1,2:N,t-2));
        else     
            for j = 1:N+1
                %% bdy updates: top + bottom (including corners)
                % for (u,v) and temp, simply enforce boundary conditions
                % for rho, use first-order one-sided difference scheme in x3
                % derivatives, periodicity in x2 derivatives. One-sided
                % scheme should not affect overall stability, but try to
                % use second order if you're having troubles near the boundary
                if j == 1
                    u(:,1,t)    = bdy1(:,1);
                    v(:,1,t)    = bdy2(:,1);
                    
                    rho(1,1,t)   = -dt/dx2*(u(2,1,t-1)-u(M+1,1,t-1)) ...
                                   - 2*dt/dx3*(v(1,2,t-1)-v(1,1,t-1)) ...
                                   + rho(1,1,t-2);
                    rho(2:M,1,t) = -dt/dx2*(u(3:M+1,1,t-1)-u(1:M-1,1,t-1)) ...
                                   - 2*dt/dx3*(v(2:M,2,t-1)-v(2:M,1,t-1)) ...
                                   + rho(2:M,1,t-2);
                    rho(M+1,1,t) = - dt/dx2*(u(1,1,t-1)-u(M,1,t-1)) ...
                                   - 2*dt/dx3*(v(M+1,2,t-1)-v(M+1,1,t-1)) ...
                                   + rho(M+1,1,t-2);           
                               
                    temp(:,1,t) = temp0(:,1);
                elseif j == N+1
                    u(:,N+1,t)    = bdy1(:,N+1);
                    v(:,N+1,t)    = bdy2(:,N+1);
                    
                    rho(1,N+1,t)   = -dt/dx2*(u(2,N+1,t-1)-u(M+1,N+1,t-1)) ...
                                     - 2*dt/dx3*(v(1,N+1,t-1)-v(1,N,t-1)) ...
                                     + rho(1,N+1,t-2);
                    rho(2:M,N+1,t) = - dt/dx2*(u(3:M+1,N+1,t-1)-u(1:M-1,N+1,t-1)) ...
                                     - 2*dt/dx3*(v(2:M,N+1,t-1)-v(2:M,N,t-1)) ...
                                     + rho(2:M,N+1,t-2); 
                    rho(M+1,N+1,t) = -dt/dx2*(u(1,N+1,t-1)-u(M,N+1,t-1)) ...
                                     -2*dt/dx3*(v(M+1,N+1,t-1)-v(M+1,N,t-1)) ...
                                     + rho(M+1,N+1,t-2);
                       
                    
                    temp(:,N+1,t) = temp0(:,N+1);
                    
                %% interior updates
                else 
                    % velocity (momentum conservation)
                    % todo: fix u veloctiy shouldn't contain bouyancy term
                    u(i,j,t) = w*(-dt/dx2*(u(i+1,j,t-1).^2-u(i-1,j,t-1).^2) ...
                                  - dt/dx3*(u(i,j+1,t-1).*v(i,j+1,t-1) - u(i,j-1,t-1).*v(i,j-1,t-1)) ...
                                  + 2*dt/(dx2.^2)*(u(i+1,j,t-1)+u(i-1,j,t-1) - u(i,j,t-2)) ...
                                  + 2*dt/(dx3.^2)*(u(i,j+1,t-1)+u(i,j-1,t-1) - u(i,j,t-2)) ...
                                  - dt/dx2*R/(sigma*q*delta)*(rho(i+1,j,t-1) - rho(i-1,j,t-1)) ...
                                  + u(i,j,t-2));
                    v(i,j,t) = w*(-dt/dx2*(v(i+1,j,t-1).*u(i+1,j,t-1) - v(i-1,j,t-1).*u(i-1,j,t-1)) ...
                                  - dt/dx3*(v(i,j+1,t-1).^2 - v(i,j-1,t-1).^2) ...
                                  + 2*dt/(dx2.^2)*(v(i+1,j,t-1) + v(i-1,j,t-1) - v(i,j,t-2)) ...
                                  + 2*dt/(dx3.^2)*(v(i,j+1,t-1) + v(i,j-1,t-1) - v(i,j,t-2)) ...
                                  - 2*dt*R/(sigma*q)*(1-q*(temp(i,j,t)-1)) ...
                                  - dt/dx3*R/(sigma*q*delta)*(rho(i,j+1,t-1) - rho(i,j-1,t-1)) ...
                                  + v(i,j,t-2));

                    % density (mass conservation)
                    rho(i,j,t) = -dt/dx2*(u(i+1,j,t-1)-u(i-1,j,t-1)) ...
                                 -dt/dx3*(v(i,j+1,t-1)-v(i,j-1,t-1)) ...
                                 + rho(i,j,t-2);
                             
                    % temperature (advection-diffusion) 
                    temp(i,j,t) = wsigma*(-2*dt/dx2*(temp(i+1,j,t-1).*u(i+1,j,t-1)-temp(i-1,j,t-1).*u(i-1,j,t-1)) ...
                                          -2*dt/dx3*(temp(i,j+1,t-1).*v(i,j+1,t-1)-temp(i,j-1,t-1).*v(i,j-1,t-1)) ...
                                          +2*dt/(dx2.^2)*(temp(i+1,j,t-1)+temp(i-1,j,t-1)-temp(i,j,t-2)) ...
                                          +2*dt/(dx3.^2)*(temp(i,j+1,t-1)+temp(i,j-1,t-1)-temp(i,j,t-2)) ...
                                  + temp(i,j,t-2));
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
    
    temp(:,:,1) = temp(:,:,2);
    temp(:,:,2) = temp(:,:,3);

    s = s+1;
    err = max(max(abs(u(:,:,3)-u(:,:,1)))) + ...
          max(max(abs(v(:,:,3)-v(:,:,1)))) + ...
          max(max(abs(temp(:,:,3)-temp(:,:,1))));

end
    


end

%surf(x3,x2,temp(:,:,3))




