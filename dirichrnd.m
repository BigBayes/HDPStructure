function theta=dirichrnd(alpha)
%DIRICHRND Dirichlet random number generator.
%   theta = DIRICHRND(alpha), returns a value for a random variable
%   distributed according to the law Dirichlet(alpha(1) ... alpha(D)).
%
% Inputs:
%           alpha 1xD parameters for the Dirchlet law. alpha(d) > 0.
%
% Outputs:
%           theta 1xD a draw from Dirichlet(alpha(1) ... alpha(D)) such
%                   that sum(theta)=1 and theta(d) > 0.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

assert(all(alpha(:)>0));
g=gamrnd(alpha,1);
theta=g/sum(g);