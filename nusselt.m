% compute nusselt number given veloctiy/temperature profile
function [nu] = nusselt(L,temp,v,a0,x2,dx3,sigma)

% 1st order finite difference in x3 direction 
dT    = 1/dx3*onesidediff(temp,L);

nu = 1/a0*trapz(x2, sigma*v(:,L,3).*temp(:,L,3) - dT);

end 



    
    
    
    
    
    