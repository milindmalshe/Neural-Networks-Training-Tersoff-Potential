%%%retuns dg(theta)/dc
function [dgdc] = dg_dc(i,j,k,r,c,d,h,R,D)

[fC] = fc(i,k,r,R,D);

rij=r(i,j);
rik=r(i,k);
rjk=r(j,k);

cosTh = (rij^2 + rik^2 - rjk^2)/(2 * rij * rik);

dgdc = fC*(2*c/d^2 - 2*c/(d^2+(h-cosTh)^2));