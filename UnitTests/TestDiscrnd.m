classdef TestDiscrnd < matlab.unittest.TestCase

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

    properties
        
    end
    
    methods (Test)
        function test2(testCase)
            N = 2;
            lambdas = [1/3 2/3];
            freqs = zeros(1, N);
            trials = 100;
            for trial = 1:trials
                result = discrnd(lambdas);
                freqs(result) = freqs(result) + 1;
            end
            
            p = chi2tst(lambdas, freqs);
            testCase.verifyGreaterThan(p, 1e-5);
        end
        
        function test3(testCase)
            N = 3;
            lambdas = [1/4 1/2 1/4];
            freqs = zeros(1, N);
            trials = 100;
            for trial = 1:trials
                result = discrnd(lambdas);
                freqs(result) = freqs(result) + 1;
            end
            
            p = chi2tst(lambdas, freqs);
            testCase.verifyGreaterThan(p, 1e-5);
        end
    end
end