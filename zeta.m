function [zetaij] = zeta(i, j, k, r,c,d,h,R,D)

[fC] = fc(i,k,r,R,D);

rij=r(i,j);
rik=r(i,k);
rjk=r(j,k);

[g_ijk] = g(rij,rik,rjk,c,d,h);

zetaij = fC*g_ijk;



