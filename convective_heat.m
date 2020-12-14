% script to compute nusselt number given the Navier-Stokes-Fourier dynamics
% at a given Rayleigh number. Uses a grid size of 40x80. Extra resolution
% around [.5,1] where onset of stability occurs


f = [1/1708,.3,.5,.6,.7,.8,.9,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10];
F = length(f);
nu = zeros(F,1);
ra   = zeros(F,1);

for j = 1:F
    f(j)

    [u,v,temp,rho,dt,x2,x3,dx2,dx3,a0,sigma,s] = benard(f(j),40);
    
    ra(j) = f(j);
    nu(j) = nusselt(1,temp,v,a0,x2,dx3,sigma);
end

plot(ra,nu)

figure(1); clf();
plot(ra(1:14),nu(1:14),'o-', 'LineWidth', 1)
hold on;
plot(ra(1:14), nu(1:14), 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Nu as a function of $\frac{R}{R*}$', 'Interpreter', 'latex', 'FontSize',24);
xlabel('$\frac{R}{R_c}$','Interpreter','latex','FontSize',24)
ylabel('Nusselt number', 'Interpreter','latex','FontSize',24)

