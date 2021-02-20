function posterior = initializePosterior()

% INITIALIZEPOSTERIOR initializes the posterior belief P(theta|x), given
% the environment parameters thetaL and thetaH 
%
%   The output 'posterior' is a function that depends on the environment 
%   (two-state ot multi-state).  By default, the environment is assumed to
%   be two-state.  To specify a mutli-state environment, call the option
%   'environment','multiState'
%
%   posterior = INITIALIZEPOSTERIOR(thetaL,thetaH) returns the
%   initialized posterior based on the environmental parameters thetaL
%   and thetaH.
%
%   Options:
%   'environment': 'twoState','multiState' (default set by form of
%                  posterior)
%
%   See also: UPDATEPOSTERIOR,RUNINFERENCE

posterior = @(x).5;

end
