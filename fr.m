%%%returns fR
function [fR] = fr(i,j,r,A,lambda1)

rij=r(i,j);
fR = A*exp(-1.*lambda1.*rij);



