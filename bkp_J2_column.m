%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% This will return Nth column matrix for J2  %%%%%%%%%%%%%%

function [dVdlambda1, dVdlambda2, dVdbeta, dVdeta, dVdc, dVdd, dVdh] = J2_column(param, r, clust_size, type, energy)

dVdlambda1=0.0;
dVdlambda2=0.0;
dVdbeta=0.0;
dVdeta=0.0;
dVdc=0.0;
dVdd=0.0;
dVdh=0.0;

for i=1:1:clust_size
    for j=1:1:clust_size
        if i==j
            continue;
        end
        A=param(type(i),type(j),type(j),1);
        B=param(type(i),type(j),type(j),2);
        lambda1=param(type(i),type(j),type(j),3);
        lambda2=param(type(i),type(j),type(j),4);
%         lambda3=param(type(i),type(j),type(j),5);
%         alpha=param(type(i),type(j),type(j),6);
        beta=param(type(i),type(j),type(j),7);
        eta=param(type(i),type(j),type(j),8);
%         c=param(type(i),type(j),type(j),9);
%         d=param(type(i),type(j),type(j),10);
%         h=param(type(i),type(j),type(j),11);
        R=param(type(i),type(j),type(j),12);
        D=param(type(i),type(j),type(j),13);        
        [fC]=fc(i,j,r,R,D);
        [fR]=fr(i,j,r,A,lambda1);
        [fB]=fb(i,j,r,B,lambda2);
        %dV/dlambda1
        dVdlambda1=dVdlambda1-r(i,j)*fC*fR;
        zetaij=0.0;
        dgdc=0.0;
        dgdd=0.0;
        dgdh=0.0;
        for k=1:1:clust_size
            if i==j || i==k
                continue;
            end
%             A=param(type(i),type(j),type(k),1);
%             B=param(type(i),type(j),type(k),2);
%             lambda1=param(type(i),type(j),type(k),3);
%             lambda2=param(type(i),type(j),type(k),4);
%             lambda3=param(type(i),type(j),type(k),5);
%             alpha=param(type(i),type(j),type(k),6);
%             beta=param(type(i),type(j),type(k),7);
%             eta=param(type(i),type(j),type(k),8);
            c=param(type(i),type(j),type(k),9);
            d=param(type(i),type(j),type(k),10);
            h=param(type(i),type(j),type(k),11);
            [dgdc]   = dgdc+dg_dc(i,j,k,r,c,d,h);
            [dgdd]   = dgdd+dg_dd(i,j,k,r,c,d,h);
            [dgdh]   = dgdh+dg_dh(i,j,k,r,c,d,h);
            [zetaij] = zetaij+zeta(i, j, k, r, c, d, h);		
        end
        
        beta_eta=beta^eta;
        zeta_eta=zetaij^eta;
        [bij]=b(zetaij,beta_n, zeta_n);
        
        %dV/dlambda2
        dVdlambda2=dVdlambda2-fC*bij*r(i,j)*fB;
        
        %dV/dbeta
        dVdbeta=dVdbeta-0.5*fB*fC*bij*zeta_eta*(beta^(eta-1))/(1+zeta_eta*beta_eta);
        
        %dV/deta
        dVdeta=dVdeta+fC*fB*bij*0.5*(log(1+beta_n*zeta_n)/eta^2-beta_eta*zeta_eta(log(beta)+log(zeta))/(n*(1+beta_eta*zeta_eta)));
        
        %common factor
        cf=0.5*fB*fC*bij*beta_eta*(zeta^(eta-1))/(1+zeta_eta*beta_eta);
        
        %dV/dc
        dVdc=dVdc-cf*dgdc;
        
        %dV/dd
        dVdd=dVdd-cf*dgdd;
        
        %dV/dh
        dVdh=dVdh-cf*dgdh;

        
    end
  
end
% dVdlambda1=energy-dVdlambda1/2.0;
% dVdlambda2=energy-dVdlambda2/2.0;
% dVdbeta=energy-dVdbeta/2.0;
% dVdeta=energy-dVdeta/2.0;
% dVdc=energy-dVdc/2.0;
% dVdd=energy-dVdd/2.0;
% dVdh=energy-dVdh/2.0;









