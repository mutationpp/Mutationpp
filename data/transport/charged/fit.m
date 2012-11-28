clear; clc;

% Start by loading in the the data that we want fitted
M = importdata('integrals.dat');

% Loop over each column to be fitted and fit to a form
% f(x) = exp(A + B*ln(x) + C*ln(x)^2 + D*ln(x)^3 + E*ln(x)^4)
cols = length(M.data(1,:));
lnT = log(M.data(:,1));

for i = 2:cols
    lnX = log(M.data(:,i));
    coeffs = polyfit(lnT, lnX, 4);
    res    = sqrt(sum(((exp(polyval(coeffs, lnT)) - M.data(:,i)).^2)./M.data(:,i)));
    fprintf('%-20s %15.7e %15.7e %15.7e %15.7e %15.7e %5.3f\n', ...
        M.colheaders{i}, coeffs(5), coeffs(4), coeffs(3), coeffs(2), coeffs(1), res);
    
    plot(lnT, polyval(coeffs, lnT), lnT, log(M.data(:,i)), 'o');
    hold on;
end

