% example file to show the usage of curve_fitting
% Author: Sebastian Moehler, 2021
% License: MIT. Do what you want, but dont blame me :-D

% measurements. For now, we take some random data points.
x = [10;  15;  -5;   2;   12; 0.1; 12;  12; -3; 6; 9; -2; 3];
y = [2.4; -15; 10.6; 2.3; -5; 0; -10; 5;  20; 4; 3; 15; 12];

% If you want, you can introduce some noise - e.g. ten percent
%noise = @(y) y + rand(size(y))*max(abs(y))*0.10;
%y = noise(y);

% so you can experiment a little, here are some possible fitting functions
f_lin = @(a, x) a(1) .* x + a(2);
f_qua  = @(a, x) a(1) .* (x - a(2)).^2 + a(3);
f_exp = @(a, x) a(1)*exp(a(2).*x);
f_pow = @(a, x) a(1).*x.^a(2);        % Don't use this if x contains negative numbers.
                                      % The derivate contains log(x)...

% the necessary derivates
df_lin = @(a, x) [ ...
  x, ...                % df/da(1)
  ones(size(x)) ...     % df/da(2). Note that we can't just write '1' because we
                        % need a column vector of the right size
];

df_qua = @(a, x) [ ...
  (x - a(2)).^2, ...
  -2*a(1)*(x-a(2)), ... 
  ones(size(x)) ...
];

df_exp = @(a, x) [ ...
  exp(a(2).*x), ...                     
  a(1)*x.*exp(a(2).*x) ...             
 ];

df_pow = @(a, x) [ ...
  x.^a(2), ...
  a(1).*x.^a(2).*log(x) ...
];

 
% we need a useful start vector. Note that a start vector which produces 0*x^0
% is not useful ;-) 
% Also: the more 'outlandish' your fitting function is, the more the results 
% may depend on the start vector. With linear fitting, almost every start vector 
% will produce the same params; but try f_pow with [-1; 1] and [1; 1] for example
% this is because gauss-newton searches for a local minimum - and there may be 
% more than 1...
a_lin = [1; 0];
a_qua = [1; 0; 0];
a_exp = [1; 1];
a_pow = [1; 1];

% now it's time to choose. 
f = f_qua;
df = df_qua;
a = a_qua;

% if you want to test f_pow, make sure you dont have negative numbers in x.
% also make sure x has no zeros. The derivate contains log(x)...
%x = abs(x);

% doing the fitting
a = curve_fitting(x, y, f, df, a)

% now let's check what we got
close all;

% first, some decent intervals to plot the resulting fitting function
plotX = linspace(min(x) - (max(x) - min(x))*0.2, max(x) + (max(x) - min(x))*0.2, 100);
plotY = f(a, plotX);

% and printing the original data ...
figure;
scatter(x, y);
hold on;
% and the fitting
plot(plotX, plotY, 'b--');
axis([min(x) - (max(x) - min(x))*0.2, max(x) + (max(x) - min(x))*0.2, min(y) - (max(y) - min(y))*0.2, max(y) + (max(y) - min(y))*0.2]);
