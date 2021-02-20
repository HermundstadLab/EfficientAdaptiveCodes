function [rV,rV0] = encodeWithDiscrete(x, phi, k, x0, params)

sigmaNoise = params(1);
nLevels    = params(2);

%response vector
rV = 1 ./ (1 + exp(-k * (phi .* x - x0) ) );
rV0 = rV;

%Gaussian noise
rV = rV + randn(size(rV)) * sigmaNoise;

rV(rV<0)=0;
rV(rV>1)=1;

bins = linspace(0,1,nLevels+1);

[~,~,binCounts] = histcounts(rV,bins);
rV = binCounts-1;
    