classdef TestDirichrnd < matlab.unittest.TestCase

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

    properties
        
    end
    
    methods (Test)
        function test1(testCase)
            alpha=[1.0,2.0,3.0,4.0];
            D=numel(alpha);
            iters=100;
            samples=nan(iters,D);
            for j=1:iters
                sample=dirichrnd(alpha);
                testCase.verifyEqual(sum(sample),1.0,'AbsTol',1e-12);
                for d=1:D
                    testCase.verifyGreaterThan(sample(d),0.0);
                end
                samples(j,:)=sample;
            end
            
            z=sum(alpha);
            for d=1:D
                a=alpha(d);
                b=z-a;
                [~,p,~]=kstest(samples(:, d)','CDF',makedist('beta',a,b));
                testCase.verifyGreaterThan(p,1e-5);
            end
        end
    end
end