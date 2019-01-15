%%% returns dg(theta)/dd
function [dgdd] = dg_dd(i, j, k, r, c, d, h,R,D)

[fC] = fc(i,k,r,R,D);

rij=r(i,j);
rik=r(i,k);
rjk=r(j,k);

cosTh = (rij^2 + rik^2 - rjk^2)/(2 * rij * rik);

% dgdd = fC*(1- 2*c/d^3 - 2*d*c^2/(d^2+(h-cosTh)^2)^2);
dgdd = fC*(-2*c^2/d^3+2*c^2*d/(d^2+(h-cosTh)^2)^2);