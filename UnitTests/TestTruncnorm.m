classdef TestTruncnorm < matlab.unittest.TestCase

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

    properties
        
    end
    
    methods (Test)
        function test1(testCase)
            function p = doTestTruncnorm(mu, sig2, a, b)
                iters = 100;
                s = nan(1, iters);
                norm = makedist('normal');
                for j = 1:iters
                    sample = truncnorm(mu, sig2, a, b);
                    testCase.verifyGreaterThan(sample, a);
                    testCase.verifyLessThan(sample, b);
                    s(j) = sample;
                end
                
                x = sort(s');
                sig = sqrt(sig2);
                denom = cdf(norm, (b - mu)/sig) - cdf(norm, (a - mu)/sig);
                c = (cdf(norm, (x - mu)/sig) - cdf(norm, (a - mu)/sig)) / denom;
                [~, p, ~] = kstest(s, 'CDF', [x c]);
            end
            
            p = doTestTruncnorm(0.0, 1.0, -1.0, 1.0);
            testCase.verifyGreaterThan(p, 1e-5);
            p = doTestTruncnorm(1.0, 1.0, 0, 2.0);
            testCase.verifyGreaterThan(p, 1e-5);
            p = doTestTruncnorm(0.0, 0.1, -1.0, 1.0);
            testCase.verifyGreaterThan(p, 1e-5);
        end
    end
end