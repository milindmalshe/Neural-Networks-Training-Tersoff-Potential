
function [rs,V]=read_data()

data=xlsread('H3_trial.xls');

%no of data pts
Q=length(data);

for iQ=1:1:Q
    V(iQ)=data(iQ,4);
    rs(1,2,iQ)=data(iQ,1);
    rs(2,1,iQ)=data(iQ,1);
    rs(1,3,iQ)=data(iQ,2);
    rs(3,1,iQ)=data(iQ,2);
    rs(2,3,iQ)=data(iQ,3);
    rs(3,2,iQ)=data(iQ,3);
end
   

