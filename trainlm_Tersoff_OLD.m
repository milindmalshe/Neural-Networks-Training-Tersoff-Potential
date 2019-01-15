function [net,tr,Ac,El,v5,v6,v7,v8] = ...
  trainlm_Tersoff_OLD(net,Pd,Tl,Ai,Q,TS,VV,TV,v9,v10,v11,v12)

% data=xlsread('O3_1point.xls');
data=xlsread('O3_UMP4_form_smooth_corrected.xls');
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

param(2,2,2,3)=-1.01154;%[-1.1884;];%4.7858;
param(2,2,2,4)=-0.512073;%[0.2502;];% 3.3228;
param(2,2,2,5)=0.0;% type_1, type_2, type_3, lamda3
param(2,2,2,6)=0.0;% type_1, type_2, type_3, alpha
param(2,2,2,7)= 33.053892;%[33.0331;];%1e1;
param(2,2,2,8)=10.104213;%[10.1063;];%1.e1;
param(2,2,2,9)=[100000.008496;];%1e5;
param(2,2,2,10)=[998.296592;];%1.00e3;  
param(2,2,2,11)= -0.867249;%-1;%[41.9552;];%-1.0;
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



%TRAINLM Levenberg-Marquardt backpropagation.
%
%  Syntax
%  
%    [net,tr] = trainlm(net,Pd,Tl,Ai,Q,TS,VV,TV)
%    info = trainlm(code)
%
%  Description
%
%    TRAINLM is a network training function that updates weight and
%    bias values according to Levenberg-Marquardt optimization.
%
%    TRAINLM(NET,Pd,Tl,Ai,Q,TS,VV,TV) takes these inputs,
%      NET - Neural network.
%      Pd  - Delayed input vectors.
%      Tl  - Layer target vectors.
%      Ai  - Initial input delay conditions.
%      Q   - Batch size.
%      TS  - Time steps.
%      VV  - Either empty matrix [] or structure of validation vectors.
%      TV  - Either empty matrix [] or structure of test vectors.
%    and returns,
%      NET - Trained network.
%      TR  - Training record of various values over each epoch:
%            TR.epoch - Epoch number.
%            TR.perf  - Training performance.
%            TR.vperf - Validation performance.
%            TR.tperf - Test performance.
%            TR.mu    - Adaptive mu value.
%
%    Training occurs according to the TRAINLM's training parameters
%    shown here with their default values:
%      net.trainParam.epochs     100  Maximum number of epochs to train
%      net.trainParam.goal         0  Performance goal
%      net.trainParam.max_fail     5  Maximum validation failures
%      net.trainParam.mem_reduc    1  Factor to use for memory/speed trade off.
%      net.trainParam.min_grad 1e-10  Minimum performance gradient
%      net.trainParam.mu       0.001  Initial Mu
%      net.trainParam.mu_dec     0.1  Mu decrease factor
%      net.trainParam.mu_inc      10  Mu increase factor
%      net.trainParam.mu_max    1e10  Maximum Mu
%      net.trainParam.show        25  Epochs between displays (NaN for no displays)
%      net.trainParam.time       inf  Maximum time to train in seconds
%
%    Dimensions for these variables are:
%      Pd - NoxNixTS cell array, each element P{i,j,ts} is a DijxQ matrix.
%      Tl - NlxTS cell array, each element P{i,ts} is a VixQ matrix.
%    Ai - NlxLD cell array, each element Ai{i,k} is an SixQ matrix.
%    Where
%      Ni = net.numInputs
%    Nl = net.numLayers
%    LD = net.numLayerDelays
%      Ri = net.inputs{i}.size
%      Si = net.layers{i}.size
%      Vi = net.targets{i}.size
%      Dij = Ri * length(net.inputWeights{i,j}.delays)
%
%    If VV or TV is not [], it must be a structure of vectors:
%      VV.PD, TV.PD - Validation/test delayed inputs.
%      VV.Tl, TV.Tl - Validation/test layer targets.
%      VV.Ai, TV.Ai - Validation/test initial input conditions.
%      VV.Q,  TV.Q  - Validation/test batch size.
%      VV.TS, TV.TS - Validation/test time steps.
%    Validation vectors are used to stop training early if the network
%    performance on the validation vectors fails to improve or remains
%    the same for MAX_FAIL epochs in a row.  Test vectors are used as
%    a further check that the network is generalizing well, but do not
%    have any effect on training.
%
%    TRAINLM(CODE) return useful information for each CODE string:
%      'pnames'    - Names of training parameters.
%      'pdefaults' - Default training parameters.
%
%  Network Use
%
%    You can create a standard network that uses TRAINLM with
%    NEWFF, NEWCF, or NEWELM.
%
%    To prepare a custom network to be trained with TRAINLM:
%    1) Set NET.trainFcn to 'trainlm'.
%       This will set NET.trainParam to TRAINLM's default parameters.
%    2) Set NET.trainParam properties to desired values.
%
%    In either case, calling TRAIN with the resulting network will
%    train the network with TRAINLM.
%
%    See NEWFF, NEWCF, and NEWELM for examples.
%
%  Algorithm
%
%    TRAINLM can train any network as long as its weight, net input,
%    and transfer functions have derivative functions.
%
%    Backpropagation is used to calculate the Jacobian jX of performance
%    PERF with respect to the weight and bias variables X.  Each
%    variable is adjusted according to Levenberg-Marquardt,
%
%      jj = jX * jX
%      je = jX * E
%      dX = -(jj+I*mu) \ je
%
%    where E is all errors and I is the identity matrix.
%
%    The adaptive value MU is increased by MU_INC until the change above
%    results in a reduced performance value.  The change is then made to
%    the network and mu is decreased by MU_DEC.
%
%    The parameter MEM_REDUC indicates how to use memory and speed to
%    calculate the Jacobian jX.  If MEM_REDUC is 1, then TRAINLM runs
%    the fastest, but can require a lot of memory. Increasing MEM_REDUC
%    to 2, cuts some of the memory required by a factor of two, but
%    slows TRAINLM somewhat.  Higher values continue to decrease the
%    amount of memory needed and increase training times.
%
%    Training stops when any of these conditions occurs:
%    1) The maximum number of EPOCHS (repetitions) is reached.
%    2) The maximum amount of TIME has been exceeded.
%    3) Performance has been minimized to the GOAL.
%    4) The performance gradient falls below MINGRAD.
%    5) MU exceeds MU_MAX.
%    6) Validation performance has increased more than MAX_FAIL times
%       since the last time it decreased (when using validation).
%
%  See also NEWFF, NEWCF, TRAINGD, TRAINGDM, TRAINGDA, TRAINGDX.

% Mark Beale, 11-31-97, ODJ 11/20/98
% Updated by Orlando De Jesús, Martin Hagan, Dynamic Training 7-20-05
% Copyright 1992-2005 The MathWorks, Inc.
% $Revision: 1.1.6.1 $ $Date: 2005/11/15 01:18:19 $

% **[ NNT2 Support ]**
if ~isa(net,'struct') & ~isa(net,'char')
  nntobsu('trainlm','Use NNT2FF and TRAIN to update and train your network.')
  switch(nargin)
  case 5, [net,tr,Ac,El] = tlm1(net,Pd,Tl,Ai,Q); return
  case 6, [net,tr,Ac,El] = tlm1(net,Pd,Tl,Ai,Q,TS); return
  case 8, [net,tr,Ac,El,v5,v6] = tlm2(net,Pd,Tl,Ai,Q,TS,VV,TV); return
  case 9, [net,tr,Ac,El,v5,v6] = tlm2(net,Pd,Tl,Ai,Q,TS,VV,TV,v9); return
  case 11, [net,tr,Ac,El,v5,v6,v7,v8] = tlm3(net,Pd,Tl,Ai,Q,TS,VV,TV,v9,v10,v11); return
  case 12, [net,tr,Ac,El,v5,v6,v7,v8] = tlm3(net,Pd,Tl,Ai,Q,TS,VV,TV,v9,v10,v11,v12); return
  end
end

% FUNCTION INFO
% =============

if isstr(net)
  switch (net)
    case 'pnames',
    net = fieldnames(trainlm('pdefaults'));
    case 'pdefaults',
    trainParam.epochs = 100;
    trainParam.goal = 0;
    trainParam.max_fail = 5;
    trainParam.mem_reduc = 1;
    trainParam.min_grad = 1e-10;
    trainParam.mu = 0.001;
    trainParam.mu_dec = 0.1;
    trainParam.mu_inc = 10;
    trainParam.mu_max = 1e10;
    trainParam.show = 25;
    trainParam.time = inf;
    net = trainParam;
    % Command to get default gradient function
    case 'gdefaults',
       % Pd contains information about a dynamic (~=0) or static (==0) network
       if Pd ==0
          net='calcjx';
       else
          net='calcjxfp';
       end
    otherwise,
    error('Unrecognized code.')
  end
  return
end

% CALCULATION
% ===========

% Parameters
epochs = net.trainParam.epochs;
goal = net.trainParam.goal;
max_fail = net.trainParam.max_fail;
mem_reduc = net.trainParam.mem_reduc;
min_grad = net.trainParam.min_grad;
mu = net.trainParam.mu;
mu_inc = net.trainParam.mu_inc;
mu_dec = net.trainParam.mu_dec;
mu_max = net.trainParam.mu_max;
show = net.trainParam.show;
time = net.trainParam.time;
gradientFcn = net.gradientFcn;

% Parameter Checking
if (~isa(epochs,'double')) | (~isreal(epochs)) | (any(size(epochs)) ~= 1) | ...
  (epochs < 1) | (round(epochs) ~= epochs)
  error('Epochs is not a positive integer.')
end
if (~isa(goal,'double')) | (~isreal(goal)) | (any(size(goal)) ~= 1) | ...
  (goal < 0)
  error('Goal is not zero or a positive real value.')
end
if (~isa(max_fail,'double')) | (~isreal(max_fail)) | (any(size(max_fail)) ~= 1) | ...
  (max_fail < 1) | (round(max_fail) ~= max_fail)
  error('Max_fail is not a positive integer.')
end
if (~isa(mem_reduc,'double')) | (~isreal(mem_reduc)) | (any(size(mem_reduc)) ~= 1) | ...
  (mem_reduc < 1) | (round(mem_reduc) ~= mem_reduc)
  error('Mem_reduc is not a positive integer.')
end
if (~isa(min_grad,'double')) | (~isreal(min_grad)) | (any(size(min_grad)) ~= 1) | ...
  (min_grad < 0)
  error('Min_grad is not zero or a positive real value.')
end
if (~isa(mu,'double')) | (~isreal(mu)) | (any(size(mu)) ~= 1) | ...
  (mu <= 0)
  error('Mu is not a positive real value.')
end
if (~isa(mu_dec,'double')) | (~isreal(mu_dec)) | (any(size(mu_dec)) ~= 1) | ...
  (mu_dec < 0) | (mu_dec > 1)
  error('Mu_dec is not a real value between 0 and 1.')
end
if (~isa(mu_inc,'double')) | (~isreal(mu_inc)) | (any(size(mu_inc)) ~= 1) | ...
  (mu_inc < 1)
  error('Mu_inc is not a real value greater than 1.')
end
if (~isa(mu_max,'double')) | (~isreal(mu_max)) | (any(size(mu_max)) ~= 1) | ...
  (mu_max <= 0)
  error('Mu_max is not a positive real value.')
end
if (mu > mu_max)
  error('Mu is greater than Mu_max.')
end
if (~isa(show,'double')) | (~isreal(show)) | (any(size(show)) ~= 1) | ...
  (isfinite(show) & ((show < 1) | (round(show) ~= show)))
  error('Show is not ''NaN'' or a positive integer.')
end
if (~isa(time,'double')) | (~isreal(time)) | (any(size(time)) ~= 1) | ...
  (time < 0)
  error('Time is not zero or a positive real value.')
end

% Constants
this = 'TRAINLM';
doValidation = ~isempty(VV);
doTest = ~isempty(TV);

% Initialize
flag_stop=0;
stop = '';
startTime = clock;
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
[perf,Ex, Vhat] = calcperf_Tersoff(net,rs,V,Q,param,clust_size,type);% [perf,El,Ac,N,Zb,Zi,Zl] = calcperf(net,X,Pd,Tl,Ai,Q,TS);
if (doValidation)
  VV.net = net;
  [vperf,Ex] = calcperf_Tersoff(net,rs,V,Q,param,clust_size,type);%changes made
  VV.perf = vperf;
  VV.numFail = 0;
end
tr = newtr(epochs,'perf','vperf','tperf','mu','gradient');  

% Testing gradient calculation for accuracy with numerical gradient
if 0 % If this flag is zero, no test will be done
  [gXt,jjt,normgXt]=calcjejj_Tersoff(net,rs,param,clust_size,Q,type,Ex);%%%%
  [gX11] = approxGrad(net,Pd,Tl,Ai,Q,TS,1e-6);
 
  c = iscell(El);
  if c
    ee = cell2mat(El);
  end
  numElementsA = prod(size(ee));
  flag_test = 1;
  if isequal(net.performFcn,'msereg'),
    gXt = 2*gXt*net.performParam.ratio/numElementsA + 2*(1-net.performParam.ratio)*X/length(X);
  elseif isequal(net.performFcn,'mse'),
    gXt = 2*gXt/numElementsA;
  elseif  isequal(net.performFcn,'mse'),
    gXt = 2*gXt;
  else
    flag_test = 0;
  end
  if flag_test,
    sseg = sumsqr(gXt-gX11);
    gXzero = gX11==0;
    den_perc = max(abs(gX11));
    if den_perc~=0,
      gXperc = 100*abs((gXt-gX11))./den_perc;
    else
      den_perc2 = max(abs(gXt));
      if den_perc2~=0,
        gXperc = 100*abs((gXt-gX11))./den_perc2;
      else
        gXperc = zeros(size(gXt));
      end
    end
  
    rmseg = sqrt(sseg/length(gXperc));
    if(any(gXperc>1)&(rmseg>1e-4))
      fprintf(['error in jacobian'  '\n'])
      zzz=clock;
      fname = cat(2,'jac_err',num2str(zzz(6)));
      fname = strrep(fname,'.','_');
      fprintf(['file name for saved data is ' fname '\n\n'])
      save(fname)
    end
  end
end
%end gradient test

% Train
for epoch=0:epochs

  % Jacobian
  [je,jj,normgX]=calcjejj_Tersoff(net,rs,param,clust_size,Q,type,Ex);%%[je,jj,normgX]=calcjejj(net,Pd,Zb,Zi,Zl,N,Ac,El,Q,TS,mem_reduc);
  normgX;
  % Training Record
  epochPlus1 = epoch+1;
  tr.perf(epochPlus1) = perf;
  tr.mu(epochPlus1) = mu;
  tr.gradient(epochPlus1) = normgX;  
  if (doValidation)
    tr.vperf(epochPlus1) = VV.perf;
  end
  if (doTest)
    [tr.tperf(epochPlus1),Ex] = calcperf_Tersoff(net,rs,V,Q,param,clust_size,type);%Changes made
  end
  
  % Stopping Criteria
  currentTime = etime(clock,startTime);
  if (perf <= goal)
    stop = 'Performance goal met.';
  elseif (epoch == epochs)
    stop = 'Maximum epoch reached, performance goal was not met.';
  elseif (currentTime > time)
    stop = 'Maximum time elapsed, performance goal was not met.';
  elseif (normgX < min_grad)
    stop = 'Minimum gradient reached, performance goal was not met.';
  elseif (mu > mu_max)
    stop = 'Maximum MU reached, performance goal was not met.';
  elseif (doValidation) & (VV.numFail > max_fail)
    stop = 'Validation stop.';
  elseif flag_stop
    stop = 'User stop.';
  end
  
  % Progress
  if isfinite(show) & (~rem(epoch,show) | length(stop))
    fprintf('%s%s%s',this,'-',gradientFcn);
  if isfinite(epochs) fprintf(', Epoch %g/%g',epoch, epochs); end
  if isfinite(time) fprintf(', Time %4.1f%%',currentTime/time*100); end
  if isfinite(goal) fprintf(', %s %g/%g',upper(net.performFcn),perf,goal); end
  if isfinite(min_grad) fprintf(', Gradient %g/%g',normgX,min_grad); end
  fprintf('\n')
  flag_stop=plotperf(tr,goal,this,epoch);
    if length(stop) fprintf('%s, %s\n\n',this,stop); end
  end
 
  % Stop when criteria indicate its time
  if length(stop)
    if (doValidation)
    net = VV.net;
  end
    break
  end
  
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
    [perf2,Ex, Vhat] = calcperf_Tersoff(net2,rs,V,Q,param,clust_size,type);%[perf2,El2,Ac2,N2,Zb2,Zi2,Zl2] = calcperf(net2,X2,Pd,Tl,Ai,Q,TS);
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
  if (doValidation)
    [vperf,Ex, Vhat] = calcperf_Tersoff(net2,rs,V,Q,param,clust_size,type);%vperf = calcperf(net,X,VV.Pd,VV.Tl,VV.Ai,VV.Q,VV.TS);
  if (vperf < VV.perf)
    VV.perf = vperf; VV.net = net; VV.numFail = 0;
  elseif (vperf > VV.perf)
      VV.numFail = VV.numFail + 1;
  end
  end
end

% Finish
tr = cliptr(tr,epoch);
