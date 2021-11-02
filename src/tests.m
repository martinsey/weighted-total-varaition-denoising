%% Main function to generate tests
function tests = exampleTest
    tests = functiontests(localfunctions);
end

%% Test Functions
function testReadImage(testCase)
    y = readImage("../data/parrotgray.png");
    y_test = readmatrix("test.txt");
    diff = abs(y_test - y);
    assert(max(diff(:)) < 5e-16)
end



function testD1p(testCase)
    [d1p, d2p] = diffOperator.fd(4);
    d1pTest = [4, 0, 0,; -4, 4, 0; 0, -4, 4; 0, 0, -4];
    d2pTest = [4, -4, 0, 0; 0, 4, -4, 0; 0, 0, 4, -4];
    
    diffdp1 = d1p - d1pTest;
    diffdp2 = d2p - d2pTest;
    
    assert(max(diffdp1(:)) < 5e-16)
    assert(max(diffdp2(:)) < 5e-16)
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end