function [errTotal, pR] = discNonlinearityLogisticDecodingManySpikes(paramsV, L, envData, nSeq, priorV, noiseSigma)

if nargin < 6
    noiseSigma = 0;
end

x0 = paramsV(1);
k = abs(paramsV(2));

%env data is data across environments
[nEnv, nSmp] = size(envData);

%prior data
if nargin < 5
    priorV = ones(nEnv, 1) ./ nEnv;
end

%compute responses to stimuli
bns = 0:L;

%number of training sequences
nHist = 1e4;

%probability distributions of responses across environments
pR = zeros(nEnv, length(bns));
rV = zeros(nEnv, nSmp);

rng('default');
for n = 1:nEnv
        tmpS = L ./ (1 + exp(-k .* (envData(n, :) + randn(size(envData(n, :))) .* noiseSigma - x0)));
    
        rV(n, :) = poissrnd(tmpS);
        pR(n, :) = histc(rV(n, :), bns);
        
        %add epsilon to avoid zeros
        pR(n, :) = pR(n, :) + 1e-3;
        
        %p(r | ENV)
        pR(n, :) = pR(n, :) ./ sum(pR(n, :));
end

%saturate
rV(rV > L) = L;

%decode
errTotal = 0;
for n = 1:nEnv
    for t = 1:nHist
        tInds = rV(n, randi(nSmp, 1, nSeq)) + 1;
        postV = prod(pR(:, tInds), 2);
        postV = postV .* priorV;
       
        [~, dec] = max(postV);
        errTotal = errTotal + double(dec == n);
    end
end
errTotal = errTotal ./ (nHist * nEnv);

