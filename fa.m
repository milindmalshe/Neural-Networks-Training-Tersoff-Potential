%%%returns fC
function [fA] = fa(i,j,r,B,lambda2)

rij=r(i,j);
fA = -1.*B*exp(-1.*lambda2.*rij);