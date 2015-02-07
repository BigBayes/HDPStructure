function p = chi2tst(probs, freqs)

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

denom = sum(freqs);
N = length(probs);
[~, p] = chi2gof(1:N, 'Ctrs', 1:N, 'Expected', denom * probs, 'Frequency', freqs);