%%%returns fC
function [fC] = fc(i,j,r, R,D)

rij=r(i,j);
if (rij < (R-D))
	
    fC = 1;
	
elseif (rij >= (R-D) && rij <= (R+D))
	                     
	fC = 0.5+0.5.*cos(pi.*(rij-(R-D))./(2.*D));

elseif (rij > (R+D))
    
	fC=0;
    
end