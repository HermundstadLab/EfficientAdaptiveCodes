function theta = generateEnvironment(h,thetaL,thetaH,varargin)

% GENERATEENVIRONMENT generates a switching environment governed by a
% single parameter theta.  
%   
%   GENERATE ENVIRONMENT assumes that theta is drawn from a environment that randomly
%   switches between values with a hazard rate h.  To specify a two-state
%   environment (in which theta switches between thetaL and thetaH), choose
%   the option 'environment', 'twoState'.  To specify a multi-state
%   environment (in which theta switches between values drawn uniformly on
%   the interval thetaL<=theta<=thetaH), choose the option 'environment',
%   'multiState'.  The environment is two-state by default.
%
%   If using a two-state environment, the default 'probe' environment switches
%   every 1/h timesteps.  To instead specify that the environment switches
%   with a fixed probability h on each timestep, turn off the 'probe'
%   by setting the option 'probe',false.
%
%   theta = GENERATEENVIRONMENT(h,thetaL,thetaH) generates a trace of theta
%   that switches between thetaL and thetaH every 1/h timesteps.
%
%   theta = GENERATEENVIRONMENT(h,thetaL,thetaH,nT) defines the total
%   duration of the inference process. By default, nT = 100/h (which
%   generates ~100 switching events).
%   
%   Options:
%   'probe'      : true (default), false
%   'environment': 'twoState' (default),'multiState'
% 
%   See also: GENERATESIGNAL,RUNINFERENCE

%-------------------- parse inputs ------------------------%

defaultT = 100./h;
p = inputParser();
p.addRequired('h',@isnumeric)
p.addRequired('thetaL',@isnumeric)
p.addRequired('thetaH',@isnumeric)
p.addOptional('nT',defaultT,@isnumeric)
p.addParameter('probe',true,@islogical)

p.parse(h,thetaL,thetaH, varargin{:});
h      = p.Results.h;
nT     = p.Results.nT;
thetaL = p.Results.thetaL;
thetaH = p.Results.thetaH;
probe  = p.Results.probe;
%----------------------------------------------------------%

theta = twoStateEnvironment(h,thetaL,thetaH,nT,'probe',probe);



end

function theta = twoStateEnvironment(h,thetaL,thetaH,nT,varargin)

%-------------------- parse inputs ------------------------%
p = inputParser();
p.addRequired('h',@isnumeric)
p.addRequired('thetaL',@isnumeric)
p.addRequired('thetaH',@isnumeric)
p.addRequired('nT',@isnumeric)
p.addParameter('probe',true,@islogical)
p.parse(h,thetaL,thetaH,nT,varargin{:});
h      = p.Results.h;
nT     = p.Results.nT;
thetaL = p.Results.thetaL;
thetaH = p.Results.thetaH;
probe  = p.Results.probe;
%----------------------------------------------------------%


if probe
    %ts = round(1./h);
    ts = floor(1./h);
    tswitch = ts:ts:nT;
else
    switches = poissrnd(h,[1,nT]);
    tswitch  = find(switches==1);
end
stateTrace = zeros(1,nT);
tstart = 1;
state  = 0;
%state = 1;
for i=1:numel(tswitch)
    stateTrace(tstart:tswitch(i)) = state;
    tstart = tswitch(i)+1;
    state  = 1-state;
end
stateTrace(tstart:end) = state;

theta = zeros(1,nT);
theta(stateTrace<.5) = thetaL;
theta(stateTrace>.5) = thetaH;

end
