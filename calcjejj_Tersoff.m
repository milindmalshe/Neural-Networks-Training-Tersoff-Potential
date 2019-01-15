function [JE,JtJ,normJE]=calcjejj_Tersoff(net,rs,param,clust_size,Q,type,Ex)

J1=zeros(Q,11);
J2=zeros(Q,11);

%%% Q=no. of rows
for iQ=1:1:Q
    
    %converting rs to my format
    for i=1:1:clust_size
         for j=1:1:clust_size
             if i~=j
                 r(i,j)=rs(i,j,iQ);
             end
         end
    end
    
    %%%J1
    [dVdw1, dVdw2, dVdb1,dVdb2] = J1_make(net,param, r, clust_size, type);
    
    
    %making iQth row
    J1(iQ,1)=dVdw1;
    J1(iQ,2)=dVdw2;
    J1(iQ,3)=dVdb1;
    J1(iQ,4)=dVdb2;
    
    %%%J2
    [dVdlambda1, dVdlambda2, dVdbeta, dVdeta, dVdc, dVdd, dVdh]=J2_make(net,param, r, clust_size, type);
    
    %making iQth row
    J2(iQ,5)=dVdlambda1;
    J2(iQ,6)=dVdlambda2;
    J2(iQ,7)=dVdbeta;
    J2(iQ,8)=dVdeta;
    J2(iQ,9)=dVdc;
    J2(iQ,10)=dVdd;
    J2(iQ,11)=dVdh;
  
end
%%% forming J matrix

J=J1+J2;
JE=J'*Ex;
JtJ=J'*J;
normJE=sqrt(JE'*JE);
