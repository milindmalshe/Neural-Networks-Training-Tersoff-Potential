function [bij] = b(zetaij, beta_eta, zeta_eta, eta)


bij=(1+ beta_eta * zeta_eta)^(-1/(2*eta));

if(zetaij == 0)
	zetaij = 1.0e-10; % to avoid infinity if zeta is 0 then in the calculation of derivative of the bij index of beta becomes -ve, therefore it becomes infinity
end
