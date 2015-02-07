% RUN_TESTS Execute all unit tests.
%   RUN_TESTS, executes all of the unit tests in the UnitTests directory.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

%%%% Make sure that the unit tests can see the HDP code.
[script_path,~,~]=fileparts(mfilename('fullpath'));
addpath(script_path);
addpath([script_path '/UnitTests']);

%%%% Run the unit tests.
import matlab.unittest.TestSuite;
totalResults=run(TestSuite.fromFolder('UnitTests'));