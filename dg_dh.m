%%% returns dg(theta)/dh
function [dgdh] = dg_dh(i, j, k, r, c, d, h,R,D)

[fC] = fc(i,k,r,R,D);

rij=r(i,j);
rik=r(i,k);
rjk=r(j,k);

cosTh = (rij^2 + rik^2 - rjk^2)/(2 * rij * rik);

dgdh = fC*( 2*(h-cosTh)*c^2/(d^2+(h-cosTh)^2)^2);