
function [du] = onesidediff(u,L)

du = .5*(3*u(:,L+2)-4*u(:,L+1)+u(:,L));

end
