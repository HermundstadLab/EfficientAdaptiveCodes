function signal  = generateSignal(theta,varargin)

% GENERATESIGNAL generates a signal x from the input trace of theta
%
%   The output signal is drawn from a gaussian distribution P(x|theta)
%   whose variance or mean is specified by theta.  By default, theta is 
%   taken to be the variance of P(x|theta).  To change this, specify the 
%   option 'estimator','mean' (see options below).  An additional
%   (optional) parameter 'param' specifies the mean (if theta=variance) or 
%   variance (if theta=mean).  If not supplied by the user, 'param' is set
%   to 0 (if theta=variance) or 1 (if theta=mean).
%
%   signal  = GENERATESIGNAL(theta) generates a signal x drawn at each
%   timestep from the distribution P(x|theta).  P(x|theta) is a Gaussian
%   distribution.  By default, the mean of P is 0, and the variance is 
%   specified by theta. 
%
%   signal  = GENERATESIGNAL(theta,param) defines the optional
%   parameter 'param' of P(x|theta).  By default, 'param' is set to 0 and
%   specifies the mean of P(x|theta).  If the user specifies the option
%   'estimator','mean', then 'param' is set to 1 and specifies the variance
%   of P(x|theta)
%
%   Options:
%   'estimator'  : 'variance' (default),'mean'
%
%   See also: GENERATEENVIRONMENT,RUNINFERENCE

%-------------------- parse inputs ------------------------%
inputExist = find(cellfun(@(x) strcmpi(x, 'estimator') , varargin));
if inputExist
    estimator = varargin{inputExist+1};
else
    estimator = 'variance';
end
if strcmp(estimator,'mean')==1
    defaultParam = 1;
else
    defaultParam = 0;
end 

defaultEstimator = 'mean';
validEstimators  = {'mean','variance'};
checkEstimator   = @(x) any(validatestring(x,validEstimators));

p = inputParser();
p.addRequired('theta',@isnumeric)
p.addOptional('param',defaultParam,@isnumeric)
p.addParameter('estimator',defaultEstimator,checkEstimator)

p.parse(theta,varargin{:});
theta       = p.Results.theta;
param       = p.Results.param;
estimator   = p.Results.estimator;
%----------------------------------------------------------%


if ~isempty(p.UsingDefaults)
   disp('Using default settings for signal generation: ')
   str = [];
   for i=1:numel(p.UsingDefaults)
       str = [str,p.UsingDefaults{i},': ',num2str(p.Results.(p.UsingDefaults{i})),'; '];
   end
   disp(str)
end

signal = zeros(1,numel(theta));
[theta,~,inds] = unique(theta);

if strcmp(estimator,'mean')
    for i=1:numel(theta)
        rng(0);
        signal(1,inds==i) = normrnd(theta(i),param,[1,numel(find(inds==i))]);
    end
elseif strcmp(estimator,'variance')
    for i=1:numel(theta)
        signal(1,inds==i) = normrnd(param,theta(i),[1,numel(find(inds==i))]);
    end
end