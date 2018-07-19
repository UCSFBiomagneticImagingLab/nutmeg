function [H,p]=nut_normtest(test)
% NUT_NORMTEST  tests if sample is normally distributed with a Kolmogorov-Smirnov
%               goodness-of-fit hypothesis test.
%
% [H,p] = normtest(test_sample)
%
% H     false if normally distributed !
% p     significance of Kolmogorov-Smirnov goodness-of-fit hypothesis test

test = test(:);
test = test(isfinite(test));

CDF=[test normcdf(test,mean(test),std(test))];
[H,p]=kstest(test,CDF);
