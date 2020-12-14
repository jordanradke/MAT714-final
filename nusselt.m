% compute nusselt number given veloctiy/temperature profile
function [nu] = nusselt(L,temp,v,a0,x2,dx3,sigma)
% parameters:
% L: an index, the cross section of box heat integrated over
% temp: temp profile, matrix
% v: vertical velocity profile, matrix
% a0: width of channel
% x2: horizontal grid
% dx3: vertical grid size
% sigma: Prandtl number
%
% returns:
% nu: dimensionless Nusselt number


% 1st order finite difference in x3 direction 
dT    = 1/dx3*onesidediff(temp,L);

nu = 1/a0*trapz(x2, sigma*v(:,L,3).*temp(:,L,3) - dT);

end 



    
    
    
    
    
    
