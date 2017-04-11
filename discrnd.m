function x=discrnd(lambda)
%DISCRND Discrete random number generator.
%   x=DISCRND(lambda), returns a random discrete variable x with support 1
%       ... D where D is the length of lambda such that
%       Pr(x=k) = lambda(k).
%
% Inputs:
%           lambda 1xD. Probability vector for the random discrete variable
%                   x.
%
% Outputs:
%           x scalar. The random draw for x.
%
%   Copyright (c) 2017, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.
%
%
%   Change log:
%
%       March 2017: Added warning for numerical issues in normalization.

u = rand();
cdf = lambda(1);
if cdf >= u
    x = 1;
    return;
end

for index = 2:numel(lambda)
    cdf = cdf + lambda(index);
    if cdf >= u
        x = index;
        return;
    end
end

warning('Numerical issue in normalization (%g != 1)', sum(lambda(:)));
[~, x] = nanmax(lambda(:)); % return the MAP