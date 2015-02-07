classdef TestClusinitial < matlab.unittest.TestCase

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

    properties
        
    end
    
    methods (Test)
        
        %%%% With 20 rows of 1s, 20 rows of 0s, need >90% consistency.
        function test2(testCase)
            X=[ones(20,20);zeros(20,20)];
            [initclus,~,~]=clusinitial(X);
            denom=0;
            numer=0;
            for i1=1:20
                for i2=1:(i1-1)
                    denom=denom+1;
                    if initclus(i1)==initclus(i2)
                        numer=numer+1;
                    end
                    
                    if initclus(i1+20)==initclus(i2+20);
                        numer=numer+1;
                    end
                end
            end
            
            for i1=1:20
                for i2=1:20
                    if initclus(i1)~=initclus(i2)
                        numer=numer+1;
                    end
                end
            end
            
            testCase.verifyGreaterThan(numer/denom,.90);
        end
        
        %%%% Test dimensions and types of outputs.
        function test1(testCase)
            N=100;
            M=20;
            X=rand(N,M);
            [initclus,C,mm]=clusinitial(X);
            testCase.verifyEqual(size(initclus,1),N);
            testCase.verifyEqual(size(initclus,2),1);
            testCase.verifyEqual(isscalar(mm),true);
            testCase.verifyLessThan(mm,0.0);
            testCase.verifyEqual(isscalar(C),true);
            testCase.verifyEqual(C,round(C),'AbsTol',1e-12);
            testCase.verifyGreaterThan(C,0);
        end
    end
end