function z = truncnorm(mu, sig2, a, b)
%TRUNCNORM Truncated normal random number generator.
%   z = TRUNCNORM(mu, sig2, a, b), returns a value for a random variable
%   distributed according to the distribution p(x) proportional to N(x, mu,
%   sig2) for x in [a, b] and p(x) = 0 for x < a or x > b, where N(x, mu,
%   sig2) is the pdf of a normal random variable with mean mu and variance
%   sig2.
%
% Inputs:
%           mu scalar. Mean of the normal random variable to truncate.
%
%           sig2 scalar. Variance of the normal random variable to 
%                   truncate.
%
%           a scalar. Left endpoint of truncation interval.
%
%           b scalar. Right endpoint of truncation interval.
%
% Outputs:
%           z scalar. a draw from a truncated normal distribution with
%                   range [a, b].

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

s=sqrt(sig2);
u=rand;
up=normcdf((b-mu)/s);
low=normcdf((a-mu)/s);
x=low+u*(up-low);
z=mu+s*norminv(x);