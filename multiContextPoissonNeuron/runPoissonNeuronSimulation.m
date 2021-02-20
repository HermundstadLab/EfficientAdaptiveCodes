%% simulation parameters
clc;
addpath(genpath('.'));

%number of time steps
p.nT = 2e3; 

%hazard rate
p.h = 1e-2; 

%number of alpha codes
p.nAlpha = 3; 
p.alphaVec = linspace(0, .4, p.nAlpha);

%report flag
saveFig = 0;

%type of code biasing 
typeOfCode = 1; %NONLOCALITY
%typeOfCode = 2; %INFERENCE


%precomputed nonlinearity parameters
p.paramPath = 'POISSON_NEURON_SIMULATION_PARAMS.mat'; 

if typeOfCode == 1
    doAdaptability = 0; doLocality = 1;
elseif typeOfCode == 2
    doAdaptability = 1; doLocality = 0;
end

load(p.paramPath);
p.L = L;
p.nEnv = nVariances^2;
p.nMemory = nHist;
p.noiseSigma = noiseSigma;
fprintf('Parameters loaded\n');

spikingPInf = bsxfun(@rdivide, spikingPInf', sum(spikingPInf'))';

if doAdaptability == 1
	fprintf('ADAPTABILITY\n'); pause(2);
	repPath = sprintf('CORRECT_LONG_SIMULATION_ADAPTABILITY_LMAX_%d_H_%.3f.pdf', L, p.h);
end

if doLocality == 1
	fprintf('NONLOCALITY\n'); pause(2);
	repPath = sprintf('CORRECT_LONG_SIMULATION_NONLOCALITY_LMAX_%d_H_%.3f.pdf', L, p.h);
end

%save path
if doAdaptability == 1
	p.savePath = sprintf('SIMULATION_ADAPTABILITY_LMAX_%d_H_%.3f.mat', L, p.h);
end

if doLocality == 1
	p.savePath = sprintf('SIMULATION_NONLOCALITY_LMAX_%d_H_%.3f.mat', L, p.h);
end

%% allocate results

%spiking responses
res.spkMAlpha = zeros(p.nT, p.nAlpha);
res.spkMInf = zeros(p.nT, 1);
res.spkMStatic = zeros(p.nT, 1);
res.spkMOracle = zeros(p.nT, 1);

%decoded stimulus values
res.xHatMAlpha = zeros(p.nT, p.nAlpha);
res.xHatMInf = zeros(p.nT, 1);
res.xHatMStatic = zeros(p.nT, 1);
res.xHatMOracle = zeros(p.nT, 1);

%posteriors
res.posteriorMAlpha = zeros(p.nT, p.nEnv, p.nAlpha);
res.posteriorMInf = zeros(p.nT, p.nEnv);

%estimated environments
res.estEnvMAlpha = zeros(p.nT, p.nAlpha);
res.estEnvMInf = zeros(p.nT, 1);

%nonlinearity parameters
res.nonlinParamsMAlpha = zeros(p.nT, 2, p.nAlpha);

%decoding vector
res.decodingVec = linspace(-16, 16, 2^14);


%% precompute conditional distributions for all alphas

%number of samples used to estimate conditional distributions of context
%given response statistics
nSmp = 1e4;

%number of environments x number of alphas
condProbs = cell(p.nEnv, p.nAlpha);

%training data
envData = zeros(p.nEnv, nSmp);
for n = 1:nVariances
    for m = 1:nVariances
        tInd = (n - 1) * nVariances + m;
        envData(tInd, :) = randn(1, nSmp) .* sqrt(varThetaV(n)) + muThetaV(m);
        
        fprintf('%d out of %d\n', (n-1)*nVariances+m, p.nEnv);
    end
end

%absoulute value
paramM(:, 1) = abs(paramM(:, 1));
paramInf(1) = abs(paramInf(1));

%figure;
for n = 1:p.nEnv
    for m = 1:p.nAlpha
        if doAdaptability
            if sum(ismember(envIndsToModify, n)) > 0
                kA = (1 - p.alphaVec(m)) * paramM(n, 1) + p.alphaVec(m) * paramInf(1);
                xA = (1 - p.alphaVec(m)) * paramM(n, 2) + p.alphaVec(m) * paramInf(2);
            else
                kA = paramM(n, 1);
                xA = paramM(n, 2);
            end
        end
        
        if doLocality
                    
                kA = (1 - p.alphaVec(m)) * paramM(n, 1) + p.alphaVec(m) * paramMismatched(n, 1);
                xA = (1 - p.alphaVec(m)) * paramM(n, 2) + p.alphaVec(m) * paramMismatched(n, 2);
        end
            
            
        [~, pR] = discNonlinearityLogisticDecodingManySpikes([xA, kA], p.L, envData, p.nMemory, p.noiseSigma);      
        condProbs{n, m} = pR;
        
    end
end


%% generate data

clc;
%variance indices and values
res.thetaLabelV = zeros(p.nT, 1);
res.thetaV = zeros(p.nT, 1);
res.muV = zeros(p.nT, 1);

%initialize and sample
tmpVar = find(logical(mnrnd(1, ones(nVariances, 1) ./ nVariances)));
tmpMu = find(logical(mnrnd(1, ones(nVariances, 1) ./ nVariances)));
tInd = (tmpVar - 1) * nVariances + tmpMu;

res.thetaLabelV(1) = tInd;
res.thetaV(1) = varThetaV(tmpVar);
res.muV(1) = muThetaV(tmpMu);

for n = 2:p.nT
    if rand() < p.h 
        tmpVar = find(logical(mnrnd(1, ones(nVariances, 1) ./ nVariances)));
        tmpMu = find(logical(mnrnd(1, ones(nVariances, 1) ./ nVariances)));

        tInd = (tmpVar - 1) * nVariances + tmpMu;
        
        while tInd == res.thetaLabelV(n - 1)
            tmpVar = find(logical(mnrnd(1, ones(nVariances, 1) ./ nVariances)));
            tmpMu = find(logical(mnrnd(1, ones(nVariances, 1) ./ nVariances)));

            tInd = (tmpVar - 1) * nVariances + tmpMu;  
        end
        
        res.thetaLabelV(n) = tInd;
        res.thetaV(n) = varThetaV(tmpVar);
        res.muV(n) = muThetaV(tmpMu);
    else
         res.thetaV(n) = res.thetaV(n - 1);
         res.muV(n) = res.muV(n - 1);

         res.thetaLabelV(n) = res.thetaLabelV(n - 1);
    end
end


%stimulus vector
res.xV = randn(p.nT, 1) .* sqrt(res.thetaV) + res.muV;


%% run the main loop simpler

%initialize
res.estEnvMAlpha(1, :) = res.thetaLabelV(1);
res.estEnvMInf(1) = res.thetaLabelV(1);

%select encoders
for n = 1:p.nAlpha
    %here predicted index is equal to the current index - trivial
    %prediction
    predEnvInd = res.estEnvMAlpha(1, n);
    
    if doAdaptability
        if sum(ismember(envIndsToModify, predEnvInd)) > 0
            res.nonlinParamsMAlpha(2, 1, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 1) + p.alphaVec(n) * paramInf(1);
            res.nonlinParamsMAlpha(2, 2, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 2) + p.alphaVec(n) * paramInf(2);
        else
            res.nonlinParamsMAlpha(2, 1, n) = paramM(predEnvInd, 1);
            res.nonlinParamsMAlpha(2, 2, n) = paramM(predEnvInd, 2);
        end
    end
    
   
    if doLocality        
        res.nonlinParamsMAlpha(2, 1, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 1) + p.alphaVec(n) * paramMismatched(predEnvInd, 1);
        res.nonlinParamsMAlpha(2, 2, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 2) + p.alphaVec(n) * paramMismatched(predEnvInd, 2); 
    end
end

tic;
for t = 2:p.nT
    
    if ~mod(t, 100)
        dT = toc;
        fprintf('Step %d out of %d --- t = %.2f[sec] --- proj. total T. = %.2f [min] ETA = %.2f [min]\n', t, p.nT, dT, dT * (p.nT / 100) / 60, dT * (p.nT / 100) / 60 * (1 - t / p.nT));
        tic;
    end

    
    % ----- generate spiking responses -----  
    
    for n = 1:p.nAlpha
        sT = p.L / (1 + exp(-res.nonlinParamsMAlpha(t, 1, n) * (res.xV(t) + randn() * noiseSigma - res.nonlinParamsMAlpha(t, 2, n))));
        res.spkMAlpha(t, n) = poissrnd(sT, 1, 1);
    end
    
    sT = p.L / (1 + exp(-paramInf(1) * (res.xV(t) + randn() * p.noiseSigma - paramInf(2))));
    res.spkMInf(t) = poissrnd(sT, 1, 1);
    
    sT = p.L / (1 + exp(-paramStatic(1) * (res.xV(t) + randn() * p.noiseSigma - paramStatic(2))));
    res.spkMStatic(t) = poissrnd(sT, 1, 1);
    
    sT = p.L / (1 + exp(-paramM(res.thetaLabelV(t), 1) * (res.xV(t) + randn() * p.noiseSigma - paramM(res.thetaLabelV(t), 2))));
    res.spkMOracle(t) = poissrnd(sT, 1, 1);
      
    
    % ----- decode the spiking response -----
    
    for n = 1:p.nAlpha
        sLV = p.L ./ (1 + exp(-res.nonlinParamsMAlpha(t, 1, n) * (res.decodingVec - res.nonlinParamsMAlpha(t, 2, n))));
        logLV = log(poisspdf(res.spkMAlpha(t, n), sLV));
        res.xHatMAlpha(t, n) = res.decodingVec(find(logLV == max(logLV), 1, 'first'));
    end

    sLV = p.L ./ (1 + exp(-paramInf(1) * (res.decodingVec - paramInf(2))));
    logLV = log(poisspdf(res.spkMInf(t), sLV));
    res.xHatMInf(t) = res.decodingVec(find(logLV == max(logLV), 1, 'first'));    
    
    
    sLV = p.L ./ (1 + exp(-paramStatic(1) * (res.decodingVec - paramStatic(2))));
    logLV = log(poisspdf(res.spkMStatic(t), sLV));
    res.xHatMStatic(t) = res.decodingVec(find(logLV == max(logLV), 1, 'first'));    
    
    
    sLV = p.L ./ (1 + exp(-paramM(res.thetaLabelV(t), 1) * (res.decodingVec - paramM(res.thetaLabelV(t), 2))));
    logLV = log(poisspdf(res.spkMOracle(t), sLV));
    res.xHatMOracle(t) = res.decodingVec(find(logLV == max(logLV), 1, 'first'));    
    
    
    % ----- update the posterior - ML decoding with p.nHist past samples -----
    
    nHist = min([t - 1, p.nMemory]);
    
    for n = 1:p.nAlpha
        pR = condProbs{res.estEnvMAlpha(t - 1, n), n};
        
        tInds = res.spkMAlpha(t-nHist:t, n);
        tInds(tInds > p.L) = p.L;
        tInds = tInds + 1;
            
        postV = prod(pR(:, tInds), 2);
        postV = postV ./ sum(postV);
        res.posteriorMAlpha(t, :, n) = postV';      
        
       [~, res.estEnvMAlpha(t, n)] = max(res.posteriorMAlpha(t, :, n));
    end
    
    tInds = res.spkMInf(t-nHist:t);
    tInds(tInds > p.L) = p.L;
    tInds = tInds + 1;
    
    postV = prod(spikingPInf(:, tInds), 2);
    postV = postV ./ sum(postV);
    res.posteriorMInf(t, :) = postV';
    
    [~, res.estEnvMInf(t)] = max(res.posteriorMInf(t, :));
    
    
    
    % ----- predict and set-up the encoders -----
   
    if t < p.nT
        for n = 1:p.nAlpha
            %here predicted index is equal to the current index - trivial
            %prediction
            predEnvInd = res.estEnvMAlpha(t, n);

            if doAdaptability
                if sum(ismember(envIndsToModify, predEnvInd)) > 0
                    res.nonlinParamsMAlpha(t + 1, 1, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 1) + p.alphaVec(n) * paramInf(1);
                    res.nonlinParamsMAlpha(t + 1, 2, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 2) + p.alphaVec(n) * paramInf(2);
                else
                    res.nonlinParamsMAlpha(t + 1, 1, n) = paramM(predEnvInd, 1);
                    res.nonlinParamsMAlpha(t + 1, 2, n) = paramM(predEnvInd, 2);
                end
            end
            
            if doLocality        
                res.nonlinParamsMAlpha(t + 1, 1, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 1) + p.alphaVec(n) * paramMismatched(predEnvInd, 1);
                res.nonlinParamsMAlpha(t + 1, 2, n) = (1 - p.alphaVec(n)) * paramM(predEnvInd, 2) + p.alphaVec(n) * paramMismatched(predEnvInd, 2);
            end
        end
    end
      
end

if saveFig
    fprintf('Saving...\n');
    save(p.savePath, 'p', 'res');
    fprintf('Saved\n');
end



%% excess vs matched error analysis

errAlpha = zeros(2, p.nAlpha);
errInf = zeros(2, 1);

nMistmatchedT = zeros(1, p.nAlpha);
mismatchedErrV = zeros(1, p.nAlpha);
matchedErrV = zeros(1, p.nAlpha);
maxMismatchedErrV = zeros(1, p.nAlpha);
avgMismatchedErrV = zeros(1, p.nAlpha);

errTrsT = 1;
for n = 1:p.nAlpha
    indsMatched = res.thetaLabelV == res.estEnvMAlpha(:, n);
    
    errAlpha(1, n) = sum((res.xHatMAlpha(indsMatched, n) - res.xV(indsMatched)).^2);
    errAlpha(2, n) = sum((res.xHatMAlpha(~indsMatched, n) - res.xV(~indsMatched)).^2);
    
    nMistmatchedT(n) = sum(~indsMatched) ./ p.nT;
    
    mismatchedErrV(n) = sum((res.xHatMAlpha(~indsMatched, n) - res.xV(~indsMatched)).^2) ./ p.nT;
    matchedErrV(n) = sum((res.xHatMAlpha(indsMatched, n) - res.xV(indsMatched)).^2) ./ p.nT;
    
    maxMismatchedErrV(n) = mean((res.xHatMAlpha(~indsMatched, n) - res.xV(~indsMatched)).^2 > errTrsT);
    avgMismatchedErrV(n) = mean((res.xHatMAlpha(~indsMatched, n) - res.xV(~indsMatched)).^2);
end

indsMatched = res.thetaLabelV == res.estEnvMInf;
errInf(1) = sum((res.xHatMInf(indsMatched) - res.xV(indsMatched)).^2);
errInf(2) = sum((res.xHatMInf(~indsMatched) - res.xV(~indsMatched)).^2);

errAlpha = bsxfun(@rdivide, errAlpha, p.nT);
errInf = bsxfun(@rdivide, errInf, p.nT);

errStatic = [0, mean((res.xV - res.xHatMStatic).^2)];
errOracle = [0, mean((res.xV - res.xHatMOracle).^2)];

figure;
subplot(1, 6, 1:4); bar([errStatic; errAlpha'; errOracle], 'stacked'); title(sprintf('ALL ENVS h=%.2f, mem=%.2f', p.h, p.nMemory), 'FontSize', 22);
lbl = cell(p.nAlpha + 2, 1); lbl{1}='STATIC'; lbl{p.nAlpha+2}='ORAC'; for n=1:p.nAlpha; lbl{n+1}=sprintf('A=%.2f',p.alphaVec(n)); end;    set(gca, 'XTickLabel', lbl);
subplot(1, 6, 5); cc = cell(p.nEnv, 1); for n = 1:p.nEnv; cc{n} = envData(n, :); end; nhist(cc); xlim([-12, 12]); ylim([0, 0.6]); 
ccmap = othercolor('Reds9', p.nAlpha);
subplot(1, 6, 6); scatter(nMistmatchedT, mismatchedErrV, 100, ccmap, 'filled'); xlabel('Time of mismatch'); ylabel('Mismatch error');
set(gcf, 'Position', [24, 142, 1275, 260]); drawnow; pause(.1);

if saveFig
    export_fig(repPath, '-append');
end


%rescaled errors
errA = sum(errAlpha);
errA = (errA - errOracle(2)) / (errStatic(2) - errOracle(2));

figure;
bar(errA);
xlim([0, p.nAlpha + 1]); 

ylabel('RELATIVE ORACLE TO STATIC ERROR');
set(gca, 'XTick', 1:p.nAlpha, 'XTickLabel', lbl(2:end-1));
set(gcf, 'Position', [55, 1335, 961, 450]);

title('RELATIVE ERROR')

if saveFig
    export_fig(repPath, '-append');
end



%% decomposition of mismatched error

figure; 

xT = linspace(0, 1, 100);
xE = linspace(min(avgMismatchedErrV) * 0.9, max(avgMismatchedErrV) * 1.1, 100);
errSurf = xT' * xE;

imagesc(xE, xT, errSurf); axis xy; hold on; colormap(othercolor('Greys3'));

ccmap = othercolor('Reds9', p.nAlpha);
scatter(avgMismatchedErrV, nMistmatchedT, 100, ccmap, 'filled', 'MarkerEdgeColor', [.5, .5, .5]); 
xlabel('MISMATCHED ERROR RATE'); ylabel('TIME OF MISMATCHED'); box;


if saveFig
    export_fig(repPath, '-append');
end


%% matched vs mismatched tradeoff

figure; 
ccmap = othercolor('Reds9', p.nAlpha);
scatter(errAlpha(1, :), errAlpha(2, :), 100, ccmap, 'filled'); 
xlabel('MATCHED ERROR'); ylabel('MISMATCHED ERROR');

if saveFig
    export_fig(repPath, '-append');
end


