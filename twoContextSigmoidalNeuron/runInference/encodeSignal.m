function [y,encodingParams] = encodeSignal(x,posterior,paramObj,varargin)

% COMPRESSSIGNAL performs compression of a sample x
%
%   By default, COMPRESSSIGNAL implements an optimal compression of the 
%   signal x that minimizes (abs(x)-abs(y))^2.  To perform a different
%   compression, specify a type of quantization (see options below; specify
%   'quantizer', type).  The quantization will be performed by mapping a 
%   signal x onto a compressed value y (if 'type'='none', y=x).
%   The compressed value y is one of n quantization levels; to
%   choose n, specify the option 'nLevels', n (by default, n=4).  The user 
%   can also specify the variable on which the quantization is performed 
%   by setting the option 'quantVar', type.  See options below.
%
%   COMPRESSSIGNAL uses precomputed quantizers that are indexed by the 
%   posterior probability P(theta|x).  By default, these quantizers are
%   loaded from within COMPRESSSIGNAL.  The user can alternatively supply
%   a structure quantObj (which stores precomputed quantizers) by specifying
%   the option 'quantObj',quantObj.  
%
%   If performing 'optimal' compression on the variable 'Thetahat' (by
%   specifying the options 'quantizer','optimal','quantVar','Thetahat'),
%   COMPRESSSIGNAL will compress the signal in a manner than minimizes
%   (theta(x)-theta(y))^2.  The computation of theta(x), theta(y) requires
%   several additional inputs, including the type of environment
%   ('twoState' or 'multiState'), the related environmental parameters h,
%   thetaL, thetaH, param, and the type of estimator ('mean' or 'variance').
%   If not specified by the user, these options will be set to default
%   values: h = .01; thetaL = 1; thetaH = 2; param = 0; 'environment','twoState',
%   'estimator','mean'.
%
%   y = COMPRESSSIGNAL(x,posterior) returns the compressed signal y,
%   computed using an optimal quantization on x given the current posterior 
%   p(theta|x).
%
%   y = COMPRESSSIGNAL(x,posterior,h,thetaL,thetaH,param) defines the 
%   environmental parameters h (hazard rate), thetaL/thetaH (low/high
%   values of theta), and param (additional parameter of P(x|theta).
%   
%   [y,ind,binCenters,binEdges,binProb] = COMPRESSSIGNAL(x,posterior,...) 
%   returns the properties of the quantizer: the quantization level (ind) 
%   to which the original signal was mapped, the set of quantization levels 
%   (binCenters), the edges of the quantization levels (binEdges), and the
%   total signal probability p(x|thetahat) in each level
%
%   Options:
%   'estimator'  : 'variance' (default),'mean'
%   'environment': 'twoState' (default),'multiState'
%   'quantizer'  : 'optimal' (default),'uniform','random','fixed'
%   'quantVar'   : for 'optimal' quantizer: 'X' (default),'Thetahat', 'XThetahat'
%                : for 'uniform'/'random' quantizers: 'X' (default), 'P','LogP'
%   'quantObj'   : structure containing precomputed quantizers
%   'quantParams': for compositional quantizer (quantVar = XThetahat)
%   'nLevels'    : number of quantization levels (integer; 4 by default)
%
%   See also: LOADQUANTIZER,RUNINFERENCE

%-------------------- parse inputs ------------------------%
  
nLevelsExist = find(cellfun(@(x) strcmpi(x, 'nLevels'), varargin));
if nLevelsExist
    defaultLevels = varargin{nLevelsExist+1};
else
    if strcmp(GLMtype,'optimal')==1
        defaultLevels = 1;
    else
        defaultLevels = 8;
    end
end  

checkPosterior = @(x)  validateattributes(x, {'function_handle'}, {'scalar'});

p = inputParser();
p.addRequired('x',@isnumeric)
p.addRequired('posterior',checkPosterior);
p.addParameter('nLevels',defaultLevels,@isnumeric)

p.parse(x,posterior,varargin{:});
x           = p.Results.x;
posterior   = p.Results.posterior;
nLevels     = p.Results.nLevels;
%----------------------------------------------------------%

encodingParams = getParams(posterior,nLevels,paramObj);
y  = encodeWithDiscreteGLM(x, encodingParams(1), encodingParams(2) ,encodingParams(3), encodingParams(4), nLevels); 

end


function y = encodeWithDiscreteGLM(x, k, x0, phi, sigmaNoise, nLevels)

%pass through nonlinearity
y = 1 ./ (1 + exp(-k * (phi .* x - x0) ) );

%add Gaussian noise
y = y + randn(size(y)) * sigmaNoise;
y(y<0)=0;
y(y>1)=1;

%discretize
bins = linspace(0,1,nLevels+1);
[~,~,binCounts] = histcounts(y,bins);
y = binCounts-1;
end


function [params] = getParams(posterior,nLevels,paramObj)

dVec = (paramObj.pVec - posterior()).^2;
[~,ind] = min(dVec(2:end-1));
ind = ind+1;



indBins = find(paramObj.nBinsV == nLevels);
if isempty(indBins)
    error('Wrong bin number');
end

params(1) = paramObj.k(indBins,ind);
params(2) = paramObj.x0(indBins,ind);
params(3) = paramObj.phi;
params(4) = paramObj.sigmaNoise;

end

