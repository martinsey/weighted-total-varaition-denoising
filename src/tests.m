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
    d1pTest = [-4, 0, 0,; 4, -4, 0; 0, 4, -4; 0, 0, 4];
    d2pTest = [-4, 4, 0, 0; 0, -4, 4, 0; 0, 0, -4, 4];
    
    diffdp1 = abs(d1p - d1pTest);
    diffdp2 = abs(d2p - d2pTest);
    
    assert(max(diffdp1(:)) < 5e-16);
    assert(max(diffdp2(:)) < 5e-16);
end

function testD1m(testCase)
    [d1m, d2m] = diffOperator.bd(4);
    d1mTest = [-4, 0, 0,; 4, -4, 0; 0, 4, -4; 0, 0, 4];
    d2mTest = [-4, 4, 0, 0; 0, -4, 4, 0; 0, 0, -4, 4];
    
    diffdm1 = abs(d1m - d1mTest);
    diffdm2 = abs(d2m - d2mTest);
    
    assert(max(diffdm1(:)) < 5e-16);
    assert(max(diffdm2(:)) < 5e-16);
end

function testGrad(testCase)
    u = [1, 2, 3; 4, 5, 6; 7, 8, 9];
    ux_test = [1, 1; 1, 1; 1, 1] * 3;
    uy_test = [3, 3, 3; 3, 3, 3] * 3;
    [ux, uy] = diffOperator.grad(u, 3);
    
    diff_ux = abs(ux - ux_test);
    diff_uy = abs(uy - uy_test);
    
    assert(max(diff_ux(:)) < 5e-16);
    assert(max(diff_uy(:)) < 5e-16);
end

function testDiv(testCase)
    p1 = [1, 2, 3; 4, 5, 6];
    p2 = [1, 2; 3, 4; 5, 6];
    
    div_test = [1, 1; 1, 1] * 3 + [3, 3; 3, 3] * 3;
    div = diffOperator.div(p1, p2);
end

function testD1mNeumann(testCase)
    [d1m, d2m] = diffOperator.bdNeumann(4);
    d1mTest = [0, -4, 0, 0; 0, 4, -4, 0; 0, 0, 4, -4; 0, 0, 0, 4];
    d2mTest = [0, 0, 0, 0; -4, 4, 0, 0; 0, -4, 4, 0; 0, 0, -4, 4];
    
    diffdm1 = abs(d1m - d1mTest);
    diffdm2 = abs(d2m - d2mTest);
    
    assert(max(diffdm1(:)) < 5e-16);
    assert(max(diffdm2(:)) < 5e-16);
end

function testD1pNeumann(testCase)
    [d1p, d2p] = diffOperator.fdNeumann(4);
    d1pTest = [-4, 0, 0, 0; 4, -4, 0, 0; 0, 4, -4, 0; 0, 0, 4, 0];
    d2pTest = [-4, 4, 0, 0; 0, -4, 4, 0; 0, 0, -4, 4; 0, 0, 0, 0];
    
    diffdp1 = abs(d1p - d1pTest);
    diffdp2 = abs(d2p - d2pTest);
    
    assert(max(diffdp1(:)) < 5e-16);
    assert(max(diffdp2(:)) < 5e-16);
end

function testLaplaceNeumann(testCase)
    n = 10;
    f = zeros(n, n);
    for x = 1:n
        for y = 1:n
            f(x, y) = 1;
        end
    end
    [d1p, d2p] = diffOperator.fdNeumann(n);
    [d1m, d2m] = diffOperator.bdNeumann(n);
    laplace_1 = d2p * d2m;
    laplce_2 = d1p * d1m;
    laplace_1
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