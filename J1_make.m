function [dVdw1, dVdw2, dVdb1,dVdb2] = J1_make(net,param, r, clust_size, type)

    for i=1:1:3
        for j=1:1:3
            if(j>i)
                dVdA(i,j)=dV_dA(i, j, param, r, type);
                dVdB(i,j)=dV_dB(i, j, param, r, clust_size,type);
            end
        end
    end

%         dVdw1= 2*(dVdA(1,2)*(-r(1,2))+dVdA(1,3)*(-r(1,3))+dVdA(2,3)*(-r(2,3))); %Multiplied by 2, to take into account A(i,j) as well as A(j,i) terms
%         dVdw2= 2*(dVdB(1,2)*(-r(1,2))+dVdB(1,3)*(-r(1,3))+dVdB(2,3)*(-r(2,3)));
%         dVdb1= 2*(dVdA(1,2)*(-1)+dVdA(1,3)*(-1)+dVdA(2,3)*(-1));
%         dVdb2= 2*(dVdB(1,2)*(-1)+dVdB(1,3)*(-1)+dVdB(2,3)*(-1));
    
    
    if strcmp(net.layers{1}.transferFcn,'purelin')
        dVdw1 = 1*(dVdA(1,2)*(-1* purelin('dn',r(1,2))* r(1,2) ) + dVdA(1,3)*(-1* purelin('dn',r(1,3))* r(1,3) )+dVdA(2,3)*(-1* purelin('dn',r(2,3))*r(2,3) )); %Multiplied by 2, to take into account A(i,j) as well as A(j,i) terms
        dVdw2 = 1*(dVdB(1,2)*(-1* purelin('dn',r(1,2))* r(1,2) ) + dVdB(1,3)*(-1* purelin('dn',r(1,3))* r(1,3) )+dVdB(2,3)*(-1* purelin('dn',r(2,3))*r(2,3) ));

        dVdb1 = 1*(dVdA(1,2)*(-1* purelin('dn',r(1,2)) ) + dVdA(1,3)*(-1* purelin('dn',r(1,3)) ) + dVdA(2,3)*(-1* purelin('dn',r(2,3)) ));
        dVdb2 = 1*(dVdB(1,2)*(-1* purelin('dn',r(1,2)) ) + dVdB(1,3)*(-1* purelin('dn',r(1,3)) ) + dVdB(2,3)*(-1* purelin('dn',r(2,3)) ));

    elseif strcmp(net.layers{1}.transferFcn,'tansig')
        
        n12 = net.IW{1,1}.*r(1,2) + net.b{1,1};
        n13 = net.IW{1,1}.*r(1,3) + net.b{1,1};
        n23 = net.IW{1,1}.*r(2,3) + net.b{1,1};
        
        a12 = tansig(n12);
        a13 = tansig(n13);
        a23 = tansig(n23);
        
        Dn_tansig_a12 = (1-a12.^2);
        Dn_tansig_a13 = (1-a13.^2);
        Dn_tansig_a23 = (1-a23.^2);
        
%         dVdw1= 2*(dVdA(1,2)*(-1* tansig('dn',r(1,2))* r(1,2) ) + dVdA(1,3)*(-1* tansig('dn',r(1,3))* r(1,3) )+dVdA(2,3)*(-1* tansig('dn',r(2,3))*r(2,3) )); %Multiplied by 2, to take into account A(i,j) as well as A(j,i) terms
%         dVdw2= 2*(dVdB(1,2)*(-1* tansig('dn',r(1,2))* r(1,2) ) + dVdB(1,3)*(-1* tansig('dn',r(1,3))* r(1,3) )+dVdB(2,3)*(-1* tansig('dn',r(2,3))*r(2,3) ));
% 
%         dVdb1= 2*(dVdA(1,2)*(-1* tansig('dn',r(1,2)) ) + dVdA(1,3)*(-1* tansig('dn',r(1,3)) ) + dVdA(2,3)*(-1* tansig('dn',r(2,3)) ));
%         dVdb2= 2*(dVdB(1,2)*(-1* tansig('dn',r(1,2)) ) + dVdB(1,3)*(-1* tansig('dn',r(1,3)) ) + dVdB(2,3)*(-1* tansig('dn',r(2,3)) ));


%%DELETE FOLLOWING CODE
%         dVdw1= 1*(dVdA(1,2)*(-1* Dn_tansig_a12* r(1,2) ) + dVdA(1,3)*(-1* Dn_tansig_a13* r(1,3) )+dVdA(2,3)*(-1* Dn_tansig_a23*r(2,3) )); %Multiplied by 2, to take into account A(i,j) as well as A(j,i) terms
%         dVdw2= 1*(dVdB(1,2)*(-1* Dn_tansig_a12* r(1,2) ) + dVdB(1,3)*(-1* Dn_tansig_a13* r(1,3) )+dVdB(2,3)*(-1* Dn_tansig_a23*r(2,3) ));
% 
%         dVdb1= 1*(dVdA(1,2)*(-1* Dn_tansig_a12 ) + dVdA(1,3)*(-1* Dn_tansig_a13 ) + dVdA(2,3)*(-1* Dn_tansig_a23 ));
%         dVdb2= 1*(dVdB(1,2)*(-1* Dn_tansig_a12 ) + dVdB(1,3)*(-1* Dn_tansig_a13 ) + dVdB(2,3)*(-1* Dn_tansig_a23 ));
%%DELETE ABOVE CODE

        
        dVdw1= 1*(dVdA(1,2)*(-1* Dn_tansig_a12(1)* r(1,2) ) + dVdA(1,3)*(-1* Dn_tansig_a13(1)* r(1,3) )+dVdA(2,3)*(-1* Dn_tansig_a23(1)*r(2,3) ));
        dVdw2= 1*(dVdB(1,2)*(-1* Dn_tansig_a12(2)* r(1,2) ) + dVdB(1,3)*(-1* Dn_tansig_a13(2)* r(1,3) )+dVdB(2,3)*(-1* Dn_tansig_a23(2)*r(2,3) ));
        
        dVdb1= 1*(dVdA(1,2)*(-1* Dn_tansig_a12(1) ) + dVdA(1,3)*(-1* Dn_tansig_a13(1) ) + dVdA(2,3)*(-1* Dn_tansig_a23(1) ));
        dVdb2= 1*(dVdB(1,2)*(-1* Dn_tansig_a12(2) ) + dVdB(1,3)*(-1* Dn_tansig_a13(2) ) + dVdB(2,3)*(-1* Dn_tansig_a23(2) ));
    end

    
        
    
    %dummy line
    i;