%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% dV/dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%K  LOOP NEEDED
function [dVdB] = dV_dB(I, J, param, r, clust_size, type)

dVdB=0.0;
for i=I:1:I
    
    for j=J:1:J
        if i==j
            continue;
        end
%         A=param(type(i),type(j),type(j),1);
%         B=param(type(i),type(j),type(j),2);
%         lambda1=param(type(i),type(j),type(j),3);
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
        [dfAdB]=dfa_db(i,j,r,lambda2);
        [fC] = fc(i,j,r, R, D);
        zetaij=0.0;
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
            R=param(type(i),type(j),type(k),12);
            D=param(type(i),type(j),type(k),13);  
            [zetaij] = zetaij+zeta(i, j, k, r, c, d, h, R, D);		
        end
        beta_eta=beta^eta;
        zeta_eta=zetaij^eta;
        [bij]=b(zetaij,beta_eta, zeta_eta, eta);
        dVdB=dVdB+bij*fC*dfAdB;
    end
    
end

% dVdB=dVdB/2.0;

function [dfAdB] = dfa_db(i,j,r,lambda2)

dfAdB = -1.0*exp(-1.*lambda2.*r(i,j));

