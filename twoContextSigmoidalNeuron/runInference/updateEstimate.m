function thetahat = updateEstimate(posterior,thetaL,thetaH)

% UPDATEESTIMATE updates the estimate thetahat, based on the posterior
% belief P(theta|x)
%
%   The input 'posterior' is a function (that can be initialized by calling
%   the function iinitializePosterior.m) that depends on the environment 
%   (two-state ot multi-state).  If the environment is not specified by the
%   user, UPDATEESTIMATE uses the form of the posterior to automatically
%   determine the correct environment
%
%   thetahat = UPDATEESTIMATE(posterior,thetaL,thetaH) returns the
%   updated estimate thetahat based on the current posterior and
%   the environmental parameters thetaL, thetaH.
%
%   Options:
%   'environment': 'twoState','multiState' (default set by form of
%                  posterior)
%
%   See also: UPDATEPOSTERIOR,RUNINFERENCE

%-------------------- parse inputs ------------------------%

checkPosterior = @(x)  validateattributes(x, {'function_handle'}, {'scalar'});

p = inputParser();
p.addRequired('posterior',checkPosterior);
p.addRequired('thetaL',@isnumeric)
p.addRequired('thetaH',@isnumeric)

p.parse(posterior,thetaL,thetaH);
thetaL    = p.Results.thetaL;
thetaH    = p.Results.thetaH;
posterior = p.Results.posterior;
%----------------------------------------------------------%

thetahat = posterior()*thetaL + (1-posterior())*thetaH;

end
