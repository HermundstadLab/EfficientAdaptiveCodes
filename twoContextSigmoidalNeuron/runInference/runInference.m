function res = runInference(varargin)

% RUNINFERENCE performs online inference with compression
%
%   RUNINFERENCE performs inference of a single distributional parameter
%   theta, and returns an estimate thetahat. By default, theta specifies the
%   variance of an incoming signal x (to specify inference of the mean, choose 
%   the option 'estimator', 'mean')
%
%   RUNINFERENCE assumes that theta is drawn from a environment that randomly
%   switches between values with a hazard rate h.  To specify a two-state
%   environment (in which theta switches between thetaL and thetaH), choose
%   the option 'environment', 'twoState'.  To specify a multi-state
%   environment (in which theta switches between values drawn uniformly on
%   the interval thetaL<=theta<=thetaH), choose the option 'environment',
%   'multiState'.  By default, the environment is a twoState environment that
%   switches every 1/h timesteps.  To instead specify that the environment switches
%   with a fixed probability h on each timestep, specify the option 'probe',false 
%   within the function call 'theta = generateEnvironment'
%
%   RUNINFERENCE performs inference on a signal x that is generated from the
%   specified values of theta.  x is drawn from a gaussian distribution
%   P(x|theta) whose variance or mean is specified by theta.  An additional
%   (optional) parameter 'param' specifies the mean or variance, respectively.
%   If not supplied by the user, 'param' is set to 0 (if theta=variance) or
%   1 (if theta=mean).  RUNINFERENCE allows the user to unput a predefined 
%   environment theta and/or signal x by choosing the options 'signal', x,
%   'theta',theta. If not provided as an input, theta and x will be generated 
%   within RUNINFERENCE.  
%
%   To perform adaptive compression, specify a type of GLM (see
%   options below; specify 'GLM', type).  The corresponding NL operation will be
%   performed at each timestep i by mapping a signal x(i) onto a compressed 
%   value y(i). The user can also specify the variable on which the NL is 
%   optimized by setting the option 'GLMvar', type.  See options below.
%
%   thetahatX = RUNINFERENCE(h,thetaL,thetaH) performs online inference
%   (without compression) of a parameter theta. Theta is generated from a 
%   two-state process that switches between thetaL and thetaH every 1/h 
%   timesteps.  By default, theta specifies the variance of P(x|theta), 
%   and the mean of P(x|theta) is set to 0.  
%
%   thetahatX = RUNINFERENCE(h,thetaL,thetaH,param) defines the optional
%   parameter 'param' of P(x|theta).  By default, 'param' is set to 0 and
%   specifies the mean of P(x|theta).  If the user specifies the option
%   'estimator','mean', then 'param' is set to 1 and specifies the variance
%   of P(x|theta)
%
%   thetahatX = RUNINFERENCE(h,thetaL,thetaH,param,nT) defines the total
%   duration of the inference process. By default, nT = 100/h (which
%   generates ~100 switching events).
%
%   [thetahatX,thetahatY] = RUNINFERENCE(h,thetaL,thetaH...) returns the estimate
%   of theta generated from a compressed signal y. By default, this compression
%   is performed using an optimal quantizer that minimizes the squared error 
%   between abs(x) and abs(y)
%
%   [thetahatX,thetahatY,theta,x] = RUNINFERENCE(h,thetaL,thetaH...) returns the
%   trace of theta and the original signal x.
%
%   [thetahatX,thetahatY,theta,x,y,quantCounts,quantVals,quantEdges] = RUNINFERENCE(h,thetaL,thetaH...) 
%   returns the compressed signal y, the number of times that the signal was 
%   compressed to each quantization level (quantCounts), the quantization
%   levels at each timestep of the inference (quantVals), and the corresponding
%   binning on signal space (quantEdges).
%
%   Options:
%   'signal'     : input signal (array of length nT; empty by default)
%   'estimator'  : 'mean' (default), 'variance'


%
%   See also: GENERATEENVIRONMENT,GENERATESIGNAL,COMPRESSSIGNAL


%-------------------- parse inputs ------------------------%

defaultGLMtype = 'adaptability';
validGLMtypes  = {'adaptability','nonlocal','fixed'};
checkGLMtype   = @(x) any(validatestring(x,validGLMtypes));

% estimator
inputExist = find(cellfun(@(x) strcmpi(x, 'estimator') , varargin));
if inputExist
    estimator = varargin{inputExist+1};
else
    estimator = 'mean';
end
defaultEstimator = 'mean';
validEstimators  = {'mean','variance'};
checkEstimator   = @(x) any(validatestring(x,validEstimators));

% environment 
defaultH = .01;
if strcmp(estimator,'mean')==1
    defaultParam  =  1;
    defaultThetaL = -2;
    defaultThetaH =  2;
else
    defaultParam  = 0;
    defaultThetaL = 1;
    defaultThetaH = 4;
end 
thetaExist = find(cellfun(@(x) strcmpi(x, 'theta') , varargin));
if thetaExist
    theta = varargin{thetaExist+1};
else
    theta = [];
end

% stimulus 
signalExist = find(cellfun(@(x) strcmpi(x, 'signal') , varargin));
if signalExist
    signal = varargin{signalExist+1};
    defaultNC = numel(signal)./(2./defaultH);
else
    signal = [];
    defaultNC = 1;
end


defaultLevels = 8;
defaultEps    = 0;


p = inputParser();
p.addOptional('h',defaultH,@isnumeric)
p.addOptional('thetaL',defaultThetaL,@isnumeric)
p.addOptional('thetaH',defaultThetaH,@isnumeric)
p.addOptional('param',defaultParam,@isnumeric)
p.addOptional('nCycles',defaultNC,@isnumeric)
p.addParameter('estimator',defaultEstimator,checkEstimator)
p.addParameter('GLMtype',defaultGLMtype,checkGLMtype)
p.addParameter('nLevels',defaultLevels,@isnumeric)
p.addParameter('eps',defaultEps,@isnumeric)
p.addParameter('signal',signal,@isnumeric)
p.addParameter('theta',theta,@isnumeric)

p.parse(varargin{:});
x           = p.Results.signal;
h           = p.Results.h;
nCycles     = p.Results.nCycles;
theta       = p.Results.theta;
thetaL      = p.Results.thetaL;
thetaH      = p.Results.thetaH;
param       = p.Results.param;
estimator   = validatestring(p.Results.estimator,validEstimators);
GLMtype     = validatestring(p.Results.GLMtype,validGLMtypes);
nLevels     = p.Results.nLevels;
eps         = p.Results.eps;
%----------------------------------------------------------%




paramObj = loadParams('estimator',estimator,'GLMtype',GLMtype,'nLevels',nLevels,'h',h,...
    'eps',eps,'thetaL',thetaL,'thetaH',thetaH);

%divide data into streams
nTperCycle = 2*floor(1./h);
nStreams = 10;
nCperStream = floor(nCycles/nStreams);
nTperStream = nTperCycle*(nCperStream+1); %add+remove one cycle beginning
nT = nTperStream*nStreams;

%generate probe theta
theta = generateEnvironment(h,thetaL,thetaH,nT,'probe', true);

%generate signal
rng('default');
rng(1);
x = generateSignal(theta,param,'estimator',estimator);
x = reshape(x',[nTperStream,nStreams]);
theta = reshape(theta',[nTperStream,nStreams]);


%initialize vectors
y           = zeros(nTperCycle*nCperStream,nStreams);
tt          = zeros(nTperCycle*nCperStream,nStreams);
xx          = zeros(nTperCycle*nCperStream,nStreams);
xhat        = zeros(nTperCycle*nCperStream,nStreams);
thetahatX   = zeros(nTperCycle*nCperStream,nStreams);
thetahatY   = zeros(nTperCycle*nCperStream,nStreams);
encodingParams1 = zeros(nTperCycle*nCperStream,nStreams);  
encodingParams2 = zeros(nTperCycle*nCperStream,nStreams);  
decodingParams1 = zeros(nTperCycle*nCperStream,nStreams);  
decodingParams2 = zeros(nTperCycle*nCperStream,nStreams);  


parfor stream = 1:nStreams
    
    thetas = theta(:,stream)
    xs = x(:,stream)
    ys = zeros(nTperStream,1);
    xhats = zeros(nTperStream,1);
    thetahatXs = zeros(nTperStream,1);
    thetahatYs = zeros(nTperStream,1);
    encodingParams1s = zeros(nTperStream,1);
    encodingParams2s = zeros(nTperStream,1);
    decodingParams1s = zeros(nTperStream,1);
    decodingParams2s = zeros(nTperStream,1);
    
    %initialize posterior
    posteriorX  = initializePosterior();
    posteriorY  = initializePosterior();

    %encode / decode signal
    [ys(1),params] = encodeSignal(xs(1),posteriorY,paramObj,'nLevels',nLevels);
    encodingParams1s(1) = params(1);
    encodingParams2s(1) = params(1);

    [xhats(1),params] = decodeSignal(ys(1),posteriorY,paramObj,'nLevels',nLevels);
    decodingParams1s(1) = params(1);
    decodingParams2s(1) = params(2);

    %initialize posterior, thetahat with no compression
    posteriorX    = updatePosterior(posteriorX,xs(1),h,thetaL,thetaH,param,'estimator',estimator);
    thetahatXs(1) = updateEstimate( posteriorX,thetaL,thetaH);

    posteriorY   = updatePosterior(posteriorY,xhats(1),h,thetaL,thetaH,param,'estimator',estimator);
    thetahatYs(1) = updateEstimate( posteriorY,thetaL,thetaH);



    %run inference
    for i=2:nTperStream

        %update estimates with compression
        [ys(i),params] = encodeSignal(xs(i),posteriorY,paramObj,'nLevels',nLevels);
        encodingParams1s(i) = params(1);
        encodingParams2s(i) = params(2);

        [xhats(i),params] = decodeSignal(ys(i),posteriorY,paramObj,'nLevels',nLevels);
        decodingParams1s(i) = params(1);
        decodingParams2s(i) = params(2);


        %update estimates
        posteriorX   = updatePosterior(posteriorX,xs(i),h,thetaL,thetaH,param,'estimator',estimator);
        thetahatXs(i) = updateEstimate( posteriorX,thetaL,thetaH);

        posteriorY   = updatePosterior(posteriorY,xhats(i),h,thetaL,thetaH,param,'estimator',estimator);
        thetahatYs(i) = updateEstimate( posteriorY,thetaL,thetaH);

    end
    
    %assign to matrices
    inds = nTperCycle+1:nTperStream;
    y(:,stream)    = ys(inds);
    xx(:,stream)   = xs(inds);
    tt(:,stream)   = thetas(inds);
    xhat(:,stream) = xhats(inds);
    thetahatX(:,stream) = thetahatXs(inds);
    thetahatY(:,stream) = thetahatYs(inds);
    encodingParams1(:,stream) = encodingParams1s(inds);
    encodingParams2(:,stream) = encodingParams2s(inds);
    decodingParams1(:,stream) = decodingParams1s(inds);
    decodingParams2(:,stream) = decodingParams2s(inds);
    
end

%collapse variables
x = xx(:);
y = y(:);
xhat  = xhat(:);
theta = tt(:);
thetahatX = thetahatX(:);
thetahatY = thetahatY(:);
encodingParams = [encodingParams1(:),encodingParams2(:)];
decodingParams = [decodingParams1(:),decodingParams2(:)];


errRec = (x-xhat).^2;
errInf = (thetahatX - thetahatY).^2;

res.full.theta     = theta;
res.full.thetahatX = thetahatX;
res.full.thetahatY = thetahatY;
res.full.x         = x;
res.full.y         = y;
res.full.xhat      = xhat;
res.full.errRec    = errRec;
res.full.errInf    = errInf;
res.full.encodingParams = encodingParams;
res.full.decodingParams = decodingParams;


nT = nStreams*nTperCycle*nCperStream; 
T = nTperCycle;
nCycles = nT/T;

% compute mutual information
mi = MI(x,y,nLevels,h);


res.avg.end.theta     = mean(reshape(theta',[T,nCycles]),2);
res.avg.end.thetahatX = mean(reshape(thetahatX,[T,nCycles]),2);
res.avg.end.thetahatY = mean(reshape(thetahatY,[T,nCycles]),2);
res.avg.end.x         = mean(reshape(x,[T,nCycles]),2);
res.avg.end.y         = mean(reshape(y,[T,nCycles]),2);
res.avg.end.xhat      = mean(reshape(xhat,[T,nCycles]),2);
res.avg.end.errRec    = mean(reshape(errRec,[T,nCycles]),2);
res.avg.end.errInf    = mean(reshape(errInf,[T,nCycles]),2);
res.avg.end.encodingParams = [mean(reshape(encodingParams(:,1),[T,nCycles]),2),mean(reshape(encodingParams(:,2),[T,nCycles]),2)];
res.avg.end.decodingParams = [mean(reshape(decodingParams(:,1),[T,nCycles]),2),mean(reshape(decodingParams(:,2),[T,nCycles]),2)];
res.avg.end.MI = mi;

res.std.end.theta     = std(reshape(theta',[T,nCycles]),[],2);
res.std.end.thetahatX = std(reshape(thetahatX,[T,nCycles]),[],2);
res.std.end.thetahatY = std(reshape(thetahatY,[T,nCycles]),[],2);
res.std.end.x         = std(reshape(x,[T,nCycles]),[],2);
res.std.end.y         = std(reshape(y,[T,nCycles]),[],2);
res.std.end.xhat      = std(reshape(xhat,[T,nCycles]),[],2);
res.std.end.errRec    = std(reshape(errRec,[T,nCycles]),[],2);
res.std.end.errInf    = std(reshape(errInf,[T,nCycles]),[],2);
res.std.end.encodingParams = [std(reshape(encodingParams(:,1),[T,nCycles]),[],2),std(reshape(encodingParams(:,2),[T,nCycles]),[],2)];
res.std.end.decodingParams = [std(reshape(decodingParams(:,1),[T,nCycles]),[],2),std(reshape(decodingParams(:,2),[T,nCycles]),[],2)];


%store mid to mid results
inds = [T/4+1:T,1:T/4];

res.avg.mid.theta       = res.avg.end.theta(inds);
res.avg.mid.thetahatX   = res.avg.end.thetahatX(inds);
res.avg.mid.thetahatY   = res.avg.end.thetahatY(inds);
res.avg.mid.x           = res.avg.end.x(inds);
res.avg.mid.y           = res.avg.end.y(inds);
res.avg.mid.xhat        = res.avg.end.xhat(inds);
res.avg.mid.errRec      = res.avg.end.errRec(inds);
res.avg.mid.errInf      = res.avg.end.errInf(inds);
res.avg.mid.encodingParams = res.avg.end.encodingParams(inds,:);
res.avg.mid.decodingParams = res.avg.end.decodingParams(inds,:);
res.avg.mid.MI          = res.avg.end.MI(inds);

res.std.mid.theta       = res.std.end.theta(inds);
res.std.mid.thetahatX   = res.std.end.thetahatX(inds);
res.std.mid.thetahatY   = res.std.end.thetahatY(inds);
res.std.mid.x           = res.std.end.x(inds);
res.std.mid.y           = res.std.end.y(inds);
res.std.mid.xhat        = res.std.end.xhat(inds);
res.std.mid.errRec      = res.std.end.errRec(inds);
res.std.mid.errInf      = res.std.end.errInf(inds);
res.std.mid.encodingParams = res.std.end.encodingParams(inds,:);
res.std.mid.decodingParams = res.std.end.decodingParams(inds,:);

folder = 'outputFolder/'; 
if strcmp(GLMtype,'nonlocal')== 1
    filename = [folder,'resNL_',estimator,'_sep',num2str(thetaH-thetaL),'_',num2str(nLevels),'levels','_eps',num2str(eps),'_h',num2str(100*h),'revised.mat'];
else
    filename = [folder,'res_',estimator,'_sep',num2str(thetaH-thetaL),'_',num2str(nLevels),'levels','_eps',num2str(eps),'_h',num2str(100*h),'.mat'];
end
save(filename,'res')

end
