

epochs = 10000;

goal = 0;

min_grad = 1e-10;

mu = 0.001;
mu_dec = 0.1;
mu_inc = 10;
mu_max = 1e10;

net = newff([1.2 4.7],[2],{'purelin'},'trainlm');
% net.IW{1} = [1.8851;5.0087];
net.b{1} =  [1.00;1.00];
net.IW{1} = [1.00;1.00];

data=xlsread('O3_UMP4_corrected_qudraticFit.xls');
clust_size=3;
type=[2,2,2];
%no of data pts
Q=length(data);
% Q=1;
for iQ=1:1:Q
    V(iQ)=data(iQ,4);
    rs(1,2,iQ)=data(iQ,1);
    rs(2,1,iQ)=data(iQ,1);
    rs(1,3,iQ)=data(iQ,2);
    rs(3,1,iQ)=data(iQ,2);
    rs(2,3,iQ)=data(iQ,3);
    rs(3,2,iQ)=data(iQ,3);
end


R=1.815;
D=0.335;
%%O-O-O
% param(2,2,2,1)=3778.5548;
% param(2,2,2,2)= 905.0261;

param(2,2,2,3)= -1.5292;
param(2,2,2,4)= -1.1087;
param(2,2,2,5)=0.0;% type_1, type_2, type_3, lamda3
param(2,2,2,6)=0.0;% type_1, type_2, type_3, alpha
param(2,2,2,7)= 0.3845;
param(2,2,2,8)= 9.2704;
param(2,2,2,9)= 100000.0023;
param(2,2,2,10)= 999.534;  
param(2,2,2,11)= 9.9291;
param(2,2,2,12)=R;% type_1, type_2, type_3, R
param(2,2,2,13)=D;% type_1, type_2, type_3, D

% param(2,2,2,3)=4.7858;
% param(2,2,2,4)=3.3228;
% param(2,2,2,5)=0.0;% type_1, type_2, type_3, lamda3
% param(2,2,2,6)=0.0;% type_1, type_2, type_3, alpha
% param(2,2,2,7)=1e0;
% param(2,2,2,8)=1.e0;
% param(2,2,2,9)=1e5;
% param(2,2,2,10)=1.00e3;  
% param(2,2,2,11)= -1.0;
% param(2,2,2,12)=R;% type_1, type_2, type_3, R
% param(2,2,2,13)=D;% type_1, type_2, type_3, D



X = getx(net);
%%%%xtra
for i=1:1:7
    if i<3
        X(i+4)=param(2,2,2,i+2);
    else
        X(i+4)=param(2,2,2,i+4);
    end
        
end
numParameters = length(X);
ii = sparse(1:numParameters,1:numParameters,ones(1,numParameters));%%%
[perf,Ex, Vhat] = calcperf_NN_Tersoff(net,rs,V,Q,param,clust_size,type);% [perf,El,Ac,N,Zb,Zi,Zl] = calcperf(net,X,Pd,Tl,Ai,Q,TS);


[gXt,jjt,normgXt]=calcjejj_Tersoff(net,rs,param,clust_size,Q,type,Ex);%%%%





for epoch=0:epochs

  % Jacobian
  [je,jj,normgX]=calcjejj_Tersoff(net,rs,param,clust_size,Q,type,Ex);%%[je,jj,normgX]=calcjejj(net,Pd,Zb,Zi,Zl,N,Ac,El,Q,TS,mem_reduc);
  normgX;
  % Training Record
  epochPlus1 = epoch+1;
  tr.perf(epochPlus1) = perf;
  tr.mu(epochPlus1) = mu;
  tr.gradient(epochPlus1) = normgX;  
  
  % Stopping Criteria
%   currentTime = etime(clock,startTime);
  if (perf <= goal)
    stop = 'Performance goal met.';
  elseif (epoch == epochs)
    stop = 'Maximum epoch reached, performance goal was not met.';
%   elseif (currentTime > time)
%     stop = 'Maximum time elapsed, performance goal was not met.';
  elseif (normgX < min_grad)
    stop = 'Minimum gradient reached, performance goal was not met.';
  elseif (mu > mu_max)
    stop = 'Maximum MU reached, performance goal was not met.';
%   elseif (doValidation) & (VV.numFail > max_fail)
%     stop = 'Validation stop.';
%   elseif flag_stop
%     stop = 'User stop.';
  end
  
  % Progress
%   if isfinite(show) & (~rem(epoch,show) | length(stop))
%     fprintf('%s%s%s',this,'-',gradientFcn);
  if isfinite(epochs) fprintf(', Epoch %g/%g',epoch, epochs); end
%   if isfinite(time) fprintf(', Time %4.1f%%',currentTime/time*100); end
  if isfinite(goal) fprintf(', %s %g/%g',upper(net.performFcn),perf,goal); end
  if isfinite(min_grad) fprintf(', Gradient %g/%g',normgX,min_grad); end
  fprintf('\n')
%   flag_stop=plotperf(tr,goal,this,epoch);
%     if length(stop) fprintf('%s, %s\n\n',this,stop); end
%   end
 
  % Stop when criteria indicate its time
%   if length(stop)
%     if (doValidation)
%     net = VV.net;
%   end
%     break
%   end
  
  % Levenberg Marquardt
  while (mu <= mu_max)
    % CHECK FOR SINGULAR MATRIX
    [msgstr,msgid] = lastwarn;
    lastwarn('MATLAB:nothing','MATLAB:nothing')
    warnstate = warning('off','all');
    dX = -(jj+ii*mu) \ je;
    [msgstr1,msgid1] = lastwarn;
    flag_inv = isequal(msgid1,'MATLAB:nothing');
    if flag_inv, lastwarn(msgstr,msgid); end;
    warning(warnstate)
    X2 = X + dX;
    
%     X2(1:2,1) = 0;%%%%%*********************Setting NN weights to 0
    
    net2 = setx(net,X2);
    for ijk=1:1:7
        if ijk<3
            param(2,2,2,ijk+2)=X2(ijk+4);
        else
            param(2,2,2,ijk+4)=X2(ijk+4);
        end
    end
    [perf2,Ex, Vhat] = calcperf_NN_Tersoff(net2,rs,V,Q,param,clust_size,type);%[perf2,El2,Ac2,N2,Zb2,Zi2,Zl2] = calcperf(net2,X2,Pd,Tl,Ai,Q,TS);
    V-Vhat;
    if (perf2 < perf) && flag_inv
      X = X2; net = net2; %Zb = Zb2; Zi = Zi2; Zl = Zl2;
      %N = N2; Ac = Ac2; El = El2; 
      perf = perf2;
      mu = mu * mu_dec;
      if (mu < 1e-20)
        mu = 1e-20;
      end
      break   % Must be after the IF
    end
    mu = mu * mu_inc;
  end
  
  % Validation
%   if (doValidation)
%       [vperf,Ex, Vhat] = calcperf_NN_Tersoff(net2,rs,V,Q,param,clust_size,type);%vperf = calcperf(net,X,VV.Pd,VV.Tl,VV.Ai,VV.Q,VV.TS);
%       if (vperf < VV.perf)
%           VV.perf = vperf; VV.net = net; VV.numFail = 0;
%       elseif (vperf > VV.perf)
%           VV.numFail = VV.numFail + 1;
%       end
%   end
end
