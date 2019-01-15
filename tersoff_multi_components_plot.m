
%Milind Malshe' code for Multicomponent system
function [Vij] = tersoff_multi_components_plot(X,rs,Q,clust_size,type,V)%GA
R=1.815;
D=0.335;
% param(2,2,2,1)=X(5);% type_1, type_2, type_3, A
% param(2,2,2,2)=X(6);% type_1, type_2, type_3, B
param(2,2,2,3)=X(5);% type_1, type_2, type_3, lamda1
param(2,2,2,4)=X(6);% type_1, type_2, type_3, lamda2
param(2,2,2,5)=0.0;% type_1, type_2, type_3, lamda3
param(2,2,2,6)=0.0;% type_1, type_2, type_3, alpha
param(2,2,2,7)=X(7);% type_1, type_2, type_3, beta
param(2,2,2,8)=X(8);% type_1, type_2, type_3, eta
param(2,2,2,9)=X(9);% type_1, type_2, type_3, c
param(2,2,2,10)=X(10);% type_1, type_2, type_3, d
param(2,2,2,11)=X(11);% type_1, type_2, type_3, h
param(2,2,2,12)=R;% type_1, type_2, type_3, R
param(2,2,2,13)=D;% type_1, type_2, type_3, D


for iQ=1:1:Q
    
    Vij=0;
    for i=1:1:clust_size
         for j=1:1:clust_size
             if i~=j
                 r(i,j)=rs(i,j,iQ);
             end
         end
    end
    for i=1:clust_size

        iMov = i;%GA
         for j=1:clust_size%list(config,1)+1 

             jBond = j;
            if(jBond == iMov)
                continue;
            end
            %rij= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
            rij= r(i,j);

            %passing params
            A=X(1)*rij+X(3);
            B=X(2)*rij+X(4);
            lamda1=param(type(i),type(j),type(j),3);
            lamda2=param(type(i),type(j),type(j),4);
            lamda3=param(type(i),type(j),type(j),5);
            alpha=param(type(i),type(j),type(j),6);
            beta=param(type(i),type(j),type(j),7);
            eta=param(type(i),type(j),type(j),8);
            c=param(type(i),type(j),type(j),9);
            d=param(type(i),type(j),type(j),10);
            h=param(type(i),type(j),type(j),11);
            R=param(type(i),type(j),type(j),12);
            D=param(type(i),type(j),type(j),13);




            [fR] = fr(rij,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);
            [fA] = fa(rij,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);
            % 			[fC,Dfcij] = fc(rij,iMov,jBond);


            zetaij=0;

            % 			Dzetak= zeros(total,3);

            %START 3-body parameters calculation
    % 		for k=2:list(iMov,1)+1
    % 			k3Body=list(iMov,k);
            for k=1:clust_size%list(config,1)+1
                k3Body = k;
                if((k3Body == iMov) | (k3Body == jBond))
                    continue
                end
    % 			rik= sqrt( (x(iMov)-x(k3Body)).^2 + (y(iMov)-y(k3Body)).^2 + (z(iMov)-z(k3Body)).^2 );
                rik= r(i,k);
    % 			rjk= sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
                rjk= r(j,k);

                %passing params
                A=X(1)*rik+X(3);
                B=X(2)*rik+X(4);
                lamda1=param(type(i),type(j),type(k),3);
                lamda2=param(type(i),type(j),type(k),4);
                lamda3=param(type(i),type(j),type(k),5);
                alpha=param(type(i),type(j),type(k),6);
                beta=param(type(i),type(j),type(k),7);
                eta=param(type(i),type(j),type(k),8);
                c=param(type(i),type(j),type(k),9);
                d=param(type(i),type(j),type(k),10);
                h=param(type(i),type(j),type(k),11);
                R=param(type(i),type(j),type(k),12);
                D=param(type(i),type(j),type(k),13);					

                [zetaij] = zeta(k3Body,rij,rik,rjk,zetaij,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);		

            end

            %passing params
            A=X(1)*rij+X(3);
            B=X(2)*rij+X(4);
            lamda1=param(type(i),type(j),type(j),3);
            lamda2=param(type(i),type(j),type(j),4);
            lamda3=param(type(i),type(j),type(j),5);
            alpha=param(type(i),type(j),type(j),6);
            beta=param(type(i),type(j),type(j),7);
            eta=param(type(i),type(j),type(j),8);
            c=param(type(i),type(j),type(j),9);
            d=param(type(i),type(j),type(j),10);
            h=param(type(i),type(j),type(j),11);
            R=param(type(i),type(j),type(j),12);
            D=param(type(i),type(j),type(j),13);

            [bij] = b(iMov,jBond,zetaij,A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);		

            [fC] = fc(rij,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);

            Vij = Vij + fC*(fR + bij * fA); 
    % 		Vij(config) = Vij(config) + fC*(fR + bij * fA); 	



        end   %end for jBond

    end   %end for iMov
    
    % end   %end for config

    Vhat(iQ) = Vij/2;
end
% Vij(config) = Vij(config)/2;%WAS WRITTEN FOR GA

figure;
for i=1:1:Q
    rij(i)=rs(1,2,i);
end
plot(rij,Vhat);
hold on;
plot(rij,V);


%**************************************************************************

function [fR] = fr(r,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D)
% global x y z;
% global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;


fR = A.*exp(-1.*lamda1.*r);


% DrijDxi= (x(iMov)-x(jBond))./r;
% DrijDyi= (y(iMov)-y(jBond))./r;
% DrijDzi= (z(iMov)-z(jBond))./r;


%-------------------------------------------------------

function [fA] = fa(r,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D)
% global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;


fA = -1.*B.*exp(-1.*lamda2.*r);

%-------------------------------------------------------

function [fC] = fc(r,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D)
% global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;


if (r < (R-D))
	fC = 1;
	% end
elseif (r >= (R-D) & r <= (R+D))
	%                     fC(j,k,i)=0.5-0.5.*sin(pi./2.*(r(j,k,i)-R)./D);
	fC = 0.5+0.5.*cos(pi.*(r-(R-D))./(2.*D));
	% end
elseif (r > (R+D))
	fC=0;
end



%--------------------------------------------------------------------------

function [g_ijk] = g(rij,rik,rjk,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D)
% global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;



cosTh = (rij^2 + rik^2 - rjk^2)/(2 * rij * rik);

th=acos(cosTh);
sinTh = sin(th);

g_ijk = 1+ c^2/d^2 - c^2/(d^2+(h-cosTh)^2);



%-------------------------------------------------------------------

function [zetaij] = zeta(k3Body,rij,rik,rjk,zetaij,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D)

% global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;



	[fC] = fc(rik,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);
	
	[g_ijk] = g(rij,rik,rjk,   A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D);
	
	zetaij = zetaij + fC*g_ijk;
		
%--------------------------------------------------------------------------

function [bij] = b(iMov,jBond,zetaij, A,B,lamda1,lamda2,lamda3,alpha,beta,eta,c,d,h,R,D)

% global A  B lamda1 lamda2 lamda3 alpha beta eta c d h R D;

bij=(1+ (beta^eta * zetaij^eta))^(-1/(2*eta));

if(zetaij == 0)
	zetaij = 1.0e-10; % to avoid infinity if zeta is 0 then in the calculation of derivative of the bij index of beta becomes -ve, therefore it becomes infinity
end
		

