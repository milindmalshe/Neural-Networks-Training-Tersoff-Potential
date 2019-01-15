%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% dV/dA


function [dVdA] = dV_dA(I, J, param, r, type)

dVdA=0.0;
for i=I:1:I
    
    for j=J:1:J
        if i==j
            continue;
        end
%         A=param(type(i),type(j),type(j),1);
%         B=param(type(i),type(j),type(j),2);
        lambda1=param(type(i),type(j),type(j),3);
%         lambda2=param(type(i),type(j),type(j),4);
%         lambda3=param(type(i),type(j),type(j),5);
%         alpha=param(type(i),type(j),type(j),6);
%         beta=param(type(i),type(j),type(j),7);
%         eta=param(type(i),type(j),type(j),8);
%         c=param(type(i),type(j),type(j),9);
%         d=param(type(i),type(j),type(j),10);
%         h=param(type(i),type(j),type(j),11);
        R=param(type(i),type(j),type(j),12);
        D=param(type(i),type(j),type(j),13);        
        [dfRdA]=dfr_da(i,j,r,lambda1);
        [fC] = fc(i,j,r, R, D);
        dVdA=dVdA+fC*dfRdA;
    end
    
end

% dVdA=dVdA/2.0;

function [dfRdA] = dfr_da(i,j,r,lambda1)

dfRdA = exp(-1.*lambda1.*r(i,j));


