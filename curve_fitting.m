function a = curve_fitting(x, y, f, df, aStart, minResDiff = 1E-4)
  % CURVE_FITTING fits a generic function to some measured data points by using 
  % damped gauss-newton to minimize the least-squares-distance
  %
  % BE AWARE: There is no error handling or input checking whatsoever. Make sure
  % the params have the right sizes.
  %
  % a = curve_fitting(x, y, f, df, aStart, minResDiff = 1E-4)
  % IN
  %    x            - n*1 column vector with the measure points
  %    y            - n*1 column vector with the corresponding measure values
  %    f = @(a, x)  - anonymus function (m*1, n*1) -> n*1. Function that should 
  %                   be fitted to the measurment. 'a' are the parameters to be 
  %                   optimized, x is the vector with the measure points.
  %    df = @(a, x) - anonymus function (m*1, n*1) -> n*m. The Jacobian matrix - 
  %                   meaning first column has df/da_1(x), second column has 
  %                   df/da_2(x) and so on. Use https://www.wolframalpha.com/ 
  %                   if you have trouble getting the derivates.
  %    aStart       - m*1 column vector with some start values for the function 
  %                   parameters
  %    minResDiff   - scalar value. Minimum distance the residue should decrease
  %                   between iterations before the process stops
  % OUT
  %    a            - m*1 column vector with the last function params
  %
  % Examples:
  %   linear fitting
  %     a = curve_fitting(
  %           [10;  15;  -5;   2;   12; 0.1; 12;  12; -3; 6; 9; -2; 3], ...
  %           [2.4; -15; 10.6; 2.3; -5; 0; -10; 5;  20; 4; 3; 15; 12], ...
  %           @(a, x) a(1) .* x + a(2), ...
  %           @(a, x) [ ...
  %             x, ...
  %             ones(size(x)) ...
  %           ], ...
  %           [1; 0]
  %         );
  %   quadratic fitting
  %     a = curve_fitting(
  %           [10;  15;  -5;   2;   12; 0.1; 12;  12; -3; 6; 9; -2; 3], ...
  %           [2.4; -15; 10.6; 2.3; -5; 0; -10; 5;  20; 4; 3; 15; 12], ...
  %           @(a, x) a(1) .* (x - a(2)).^2 + a(3), ...
  %           @(a, x) [ ...
  %             (x - a(2)).^2, ...
  %             -2*a(1)*(x-a(2)), ... 
  %             ones(size(x)) ...
  %           ], ...
  %           [1; 0; 0]
  %         );
  %   exponential fitting
  %     a = curve_fitting(
  %           [10;  15;  -5;   2;   12; 0.1; 12;  12; -3; 6; 9; -2; 3], ...
  %           [2.4; -15; 10.6; 2.3; -5; 0; -10; 5;  20; 4; 3; 15; 12], ...
  %           @(a, x) a(1)*exp(a(2).*x), ...
  %           @(a, x) [ ...
  %             exp(a(2).*x), ...
  %             a(1)*x.*exp(a(2).*x) ...
  %           ], ...
  %           [1; 1]
  %         );
  %   power fitting
  %     a = curve_fitting(
  %           [10;  15;  -5;   2;   12; 0.1; 12;  12; -3; 6; 9; -2; 3], ...
  %           [2.4; -15; 10.6; 2.3; -5; 0; -10; 5;  20; 4; 3; 15; 12], ...
  %           @(a, x) a(1).*x.^a(2), ...
  %           @(a, x) [ ...
  %             x.^a(2), ...
  %             a(1).*x.^a(2).*log(x) ...
  %           ], ...
  %           [1; 1]
  %         );
  %
  % See also example.m for a usage example.
  %
  % Author: Sebastian Moehler, 2021
  % License: MIT. Do what you want, but dont blame me :-D
  
 

  % for better readability, define some functions for later use
  res = @(a) 1/2 * norm(f(a, x) - y)^2;                             % residue
  reg = @(a) sum((f(a, x) - mean(y)).^2) / sum((y - mean(y)).^2);   % coefficient of determination
  B = @(a) df(a, x);                                                % Jacobian
  d = @(a) y - f(a, x);                                             % distance between measurement and calculation

  % let's go 
  a = aStart;
  oldRes = res(a);

  i = 0;
  fprintf('start-residue: %e\n\n', oldRes);

  % We want at least five iterations, after that we check if the residue converges.
  % there's no point in trying to push the residue to zero. With least squares, 
  % there will always be some distance between the 
  % measurement and the calculation.
  while ( i < 5 || abs(oldRes - res(a)) > minResDiff) && i < 1000
    i++;
    oldRes = res(a);

    % one gauss-newton-step. We could implement a qr-solver for this, but why bother 
    % if we can use matlab functions?
    v = B(a) \ d(a);
    
    % damping... factors being 1, 1/2, 1/4, 1/8, ...
    minRes = res(a + v);
    bestFactor = 1;
    
    for k = 1 : 4 % damping factors < 2^-5 (~0.03) are not useful, you have no 
                  % movement there and will get nowhere
      aktRes = res(a + (2^-k) * v);
      if aktRes < minRes
        bestFactor = 2^-k;
        minRes = aktRes;
      endif
    end  
    a = a + bestFactor * v;
    fprintf('Step %03d: res = %e, R^2 = %f, damping = %f\n', i, res(a), reg(a), bestFactor);
  end

 end;
