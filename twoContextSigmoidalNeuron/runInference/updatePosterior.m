function posterior = updatePosterior(posterior,x,h,thetaL,thetaH,varargin)

% UPDATEPOSTERIOR updates the posterior belief P(theta|x), given the
% current posterior and the input signal x
%
%   The input 'posterior' is a function (that can be initialized by calling
%   the function initializePosterior.m) that depends on the environment 
%   (two-state ot multi-state).  If the environment is not specified by the
%   user, UPDATEPOSTERIOR uses the form of the posterior to automatically
%   determine the correct environment
%  
%   The correct posterior depends on whether the theta represents the mean
%   or the variance.  By default, UPDATEPOSTERIOR assumes that
%   theta=variance.  To specify theta=mean, call the option
%   'estimator','mean'
%
%   posterior = UPDATEPOSTERIOR(posterior,x,h,thetaL,thetaH) returns the
%   updated posterior based on the previous posterior, the input signal x,
%   the environmental parameters (hazard rate h, thetaL, thetaH).
%
%   Options:
%   'estimator'  : 'variance' (default),'mean'
%   'environment': 'twoState','multiState' (default set by form of
%                  posterior)
%
%   See also: INITIALIZEPOSTERIOR,RUNINFERENCE


%-------------------- parse inputs ------------------------%
inputExist = find(cellfun(@(x) strcmpi(x, 'estimator') , varargin));
if inputExist
    estimator = varargin{inputExist+1};
else
    estimator = 'mean';
end
if strcmp(estimator,'mean')==1
    defaultParam = 1;
else
    defaultParam = 0;
end   
defaultEstimator = 'mean';
validEstimators  = {'mean','variance'};
checkEstimator   = @(x) any(validatestring(x,validEstimators));


checkPosterior = @(x)  validateattributes(x, {'function_handle'}, {'scalar'});

p = inputParser();
p.addRequired('posterior',checkPosterior);
p.addRequired('x',@isnumeric)
p.addRequired('h',@isnumeric)
p.addRequired('thetaL',@isnumeric)
p.addRequired('thetaH',@isnumeric)
p.addOptional('param',defaultParam,@isnumeric)
p.addParameter('estimator',defaultEstimator,checkEstimator)

p.parse(posterior,x,h,thetaL,thetaH,varargin{:});
h         = p.Results.h;
x         = p.Results.x;
thetaL    = p.Results.thetaL;
thetaH    = p.Results.thetaH;
param     = p.Results.param;
posterior = p.Results.posterior;
estimator = validatestring(p.Results.estimator,validEstimators);
%----------------------------------------------------------%


f = posterior()*(1-h) + h.*(1-posterior());

if strcmp(estimator,'mean')==1
     
    norm = f.*gaussian(x,thetaL,param) + (1-f).*gaussian(x,thetaH,param);
    posterior = @(sig) f.*gaussian(x,thetaL,param)./norm;

else
    
    norm = f.*gaussian(x,param,thetaL) + (1-f).*gaussian(x,param,thetaH);
    posterior = @(sig) f.*gaussian(x,param,thetaL)./norm;
end


end
