function paramObj = loadParams(varargin)

% LOADGLMPARAMS loads precomputed parameters of GLM
%
%   quanparamObjtObj = LOADGLMPARAMS() returns a structure
%   (paramObj) that contains precomputed GLM parameters.  It is loaded from 
%   a file based on the on the type of GLM (GLMtype) and the variable being
%   optimized (GLMvar).  By default, LOADGLMPARAMS will load the set of
%   GLM parameters that are optimized to minimize error in X.  To specify a 
%   different scheme, set the options below.
%
%   Options:
%   'estimator'  : 'variance' (default), 'mean'
%   'GLMtype'    : 'adpatability' (default),'nonlocal','fixed'
%   'GLMvar'     : for 'optimal' GLM: 'X' (default),'Thetahat', 'XThetahat'
%                : for 'fixed' GLM: 'X' (default), 'Thetahat'
%   'GLMparams'  : for compositional GLM (GLMvar = XThetahat)
%
%   See also: INITIALIZEPOSTERIOR,RUNINFERENCE

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
if strcmp(estimator,'mean')==1
    defaultThetaL = -2;
    defaultThetaH =  2;
else
    defaultThetaL = 1;
    defaultThetaH = 4;
end 

defaultEps = 0;
defaultH = .01;
defaultLevels = 8;

p = inputParser();
p.addParameter('estimator',defaultEstimator,checkEstimator)
p.addParameter('GLMtype',defaultGLMtype,checkGLMtype)
p.addParameter('eps',defaultEps,@isnumeric)
p.addParameter('h',defaultH,@isnumeric)
p.addOptional('thetaL',defaultThetaL,@isnumeric)
p.addOptional('thetaH',defaultThetaH,@isnumeric)
p.addParameter('nLevels',defaultLevels,@isnumeric)

p.parse(varargin{:});
estimator = validatestring(p.Results.estimator,validEstimators);
GLMtype   = validatestring(p.Results.GLMtype,validGLMtypes);
eps       = p.Results.eps;
h         = p.Results.h;
nLevels   = p.Results.nLevels;
thetaL    = p.Results.thetaL;
thetaH    = p.Results.thetaH;

%--------------------------------------------------------%
dir = '/Users/hermundstada/Dropbox/adaptiveInferenceSimulations4/GLMparams/';

%modified 3/13/2020 to multiple h by 1000, rather than 100 (default)
separation = thetaH - thetaL;
if strcmp(GLMtype,'nonlocal')==1
    %filename = [dir,'discreteNonlocalGLM',estimator,num2str(nLevels),'levels',num2str(separation),'sep',num2str(1000*h),'h.mat'];
    filename = [dir,'discreteNonlocalGLM',estimator,num2str(nLevels),'levels',num2str(separation),'sep',num2str(100*h),'h_revised.mat'];
else
    %filename = [dir,'discreteGLM',estimator,num2str(nLevels),'levels',num2str(separation),'sep',num2str(1000*h),'h.mat'];
    filename = [dir,'discreteGLM',estimator,num2str(nLevels),'levels',num2str(separation),'sep',num2str(100*h),'h_revised.mat'];
end

rs = load(filename);
paramObj = rs.paramObj.(['eps',num2str(eps)]);

end