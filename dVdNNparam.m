function [dVdW1, dVdW2, dVdb1, dVdb2, dVdA,dVdB]=dVdNNparam(net,param,rs,Q,clust_size)


for iQ=1:1:iQ
    dVdA=dV_dA();
    dVdB=dV_dB();
    dVdW1(iQ)=dVdA(1)*(-rs(1,2,iQ))+dVdA(2)*(-rs(1,3,iQ))+dVdA(3)*(-rs(2,3,iQ));
    dVdW2(iQ)=dVdB(1)*(-rs(1,2,iQ))+dVdB(2)*(-rs(1,3,iQ))+dVdB(3)*(-rs(2,3,iQ));
    dVdb1(iQ)=-1.0;
    dVdb2(iQ)=-1.0;
    
end