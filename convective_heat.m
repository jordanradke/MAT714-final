F = 10;
nu = zeros(F+1,1);
ra   = zeros(F+1,1);

for f = 0:F
    f
    if f == 0
        [u,v,temp,rho,dt,x2,x3,dx2,dx3,a0,sigma] = benard(f/1750);
    else
        [u,v,temp,rho,dt,x2,x3,dx2,dx3,a0,sigma] = benard(f/2);
    end
    
    ra(f+1) = f/2;
    nu(f+1) = nusselt(1,temp,v,a0,x2,dx3,sigma);
end

plot(ra,nu)