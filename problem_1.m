
clc();
clear();
data=xlsread('O3_1point.xls');

%no of data pts
Q=length(data);
Q = 1;
for iQ=1:1:Q
    %dummy Targets
    T(1,iQ)=data(iQ,4);
    T(2,iQ)=data(iQ,4);
    %dummy input
    P(1,iQ)=data(iQ,1);
    %P(2,iQ)=data(iQ,2);
    %P(3,iQ)=data(iQ,3);
end

net = newff([0.5 3.5],[2],{'purelin'},'trainlm_Tersoff');
net.trainParam.epochs=5000;
% 
% net.IW{1} = [0.0006;122.3017;];
% net.b{1} =  [1.0288;113.5513;];



net = train(net,P,T);
