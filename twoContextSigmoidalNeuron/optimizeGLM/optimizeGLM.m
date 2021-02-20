function [res,paramObjTemp] = optimizeGLM(estimator,nLevels,sig,h,saveOut)

if nargin<5
    saveOut = false;
end



%discretization of prior probability (pLow)
nP = 100;
pVec = linspace(0, 1, nP+1);
pVec = pVec(2:end-1);

%environment parameters
if strcmp(estimator,'mean')==1
    thetaL = -sig;
    thetaH = sig;
    infError = @dAveragePosteriorTwoStateMean;
elseif strcmp(estimator,'variance')==1
    thetaL = 1;
    thetaH = sig;
    infError = @dAveragePosteriorTwoState;
else
    error('unrecognized estimator');
end


%encoding noise
sigmaNoise = 1e-2;

%encoding functions
encode = @encodeWithDiscrete;
paramsEncode = [sigmaNoise,nLevels];

%auxiliary parameter
phi = 1;

%parameters to be optimized
N = 201;


if strcmp(estimator,'mean')==1    
    k = linspace(0, 2, N);          %slope
    x0 = linspace(-2.5, 2.5, N);    %shift
elseif strcmp(estimator,'variance')==1
    
    k = linspace(0, 1.5, N);         %slope
    x0 = linspace(-1, 0, N);         %shift 
    error('unrecognized estimator');
end

%training dataset
xCore = randn(1e5, 1);


nn = numel(pVec);
errMInf = zeros(N, N, nn);
errMRec = zeros(N, N, nn);
resM = zeros(N, N, nn);
p1 = zeros(N, N, nn);
p0 = zeros(N, N, nn);

errMInf_flip = zeros(N, N, nn);
errMRec_flip = zeros(N, N, nn);
resM_flip = zeros(N, N, nn);
p1_flip = zeros(N, N, nn);
p0_flip = zeros(N, N, nn);

parfor i = 1:numel(pVec)

    %set random number generator (depends on Matlab version)
    rng('default');
    
    %prior parameters
    pLowPrev = pVec(i);
    if strcmp(estimator,'mean')==1
        priorMu = (pLowPrev * thetaL) + (1 - pLowPrev) * thetaH;
        priorSigma = 1;
    else
        priorMu = 0;
        priorSigma = (pLowPrev * thetaL) + (1 - pLowPrev) * thetaH;
    end

    %training data
    xM = xCore * priorSigma + priorMu;

    
    %optimize
    for n = 1:N
        for m = 1:N
            %disp([m,n])
            rV = encode(xM, phi, k(n), x0(m), paramsEncode);
            [xHat,paramsOLE] = decodeOLEFromDiscreteGLM(rV, xM);

            errMRec(n, m, i) = mean((xM - xHat).^2);
            errMInf(n, m, i) = mean(infError(xM, xHat, pLowPrev, thetaL, thetaH, h));
            p1(n, m, i) = paramsOLE(1);
            p0(n, m, i) = paramsOLE(2);
            resM(n, m, i) = mean(rV);
            
            if strcmp(estimator,'variance')==1
                rV = encode(xM, phi, k(n), -x0(m), paramsEncode);
                [xHat,paramsOLE] = decodeOLEFromDiscreteGLM(rV, xM);

                errMRec_flip(n, m, i) = mean((xM - xHat).^2);
                errMInf_flip(n, m, i) = mean(infError(xM, xHat, pLowPrev, thetaL, thetaH, h));
                p1_flip(n, m, i) = paramsOLE(1);
                p0_flip(n, m, i) = paramsOLE(2);
                resM_flip(n, m, i) = mean(rV);
            end
        end  
    end
    
    fprintf('iteration %d\n', i);

end

res.errMInf = errMInf;
res.errMRec = errMRec;
res.resM = resM;
res.p1 = p1;
res.p0 = p0;

if strcmp(estimator,'variance')==1
    res.errMInf_flip = errMInf_flip;
    res.errMRec_flip = errMRec_flip;
    res.resM_flip = resM_flip;
    res.p1_flip = p1_flip;
    res.p0_flip = p0_flip;
end

if strcmp(estimator,'mean')==1
    priorSigma = 1;
    estimatorParam = priorSigma;
else
    priorMu = 0;
    estimatorParam = priorMu;
end


Nn = length(pVec);

nE = 100;
eps = linspace(0,1,nE+1);
x0V = zeros(Nn, 2+nE+1);
kV  = zeros(Nn, 2+nE+1);
p0  = zeros(Nn, 2+nE+1);
p1  = zeros(Nn, 2+nE+1);


nSmooth = 12; 

for n = 1:Nn
 
    %---- reconstruction + inference ----%
    if strcmp(estimator,'mean')==1
        %enforce symmetry
        LagM  = .5*(squeeze(res.errMRec(:, :, n)) + fliplr(squeeze(res.errMRec(:, :, Nn-n+1))));
        LagMI = .5*(squeeze(res.errMInf(:, :, n)) + fliplr(squeeze(res.errMInf(:, :, Nn-n+1))));
    elseif strcmp(estimator,'variance')==1
        LagM  = .5*(squeeze(res.errMRec(:, :, n)) + squeeze(res.errMRec_flip(:, :, n)));
        LagMI = .5*(squeeze(res.errMInf(:, :, n)) + squeeze(res.errMInf_flip(:, :, n)));
    else
        error('unrecognized estimator')
    end

    
    %recale error to be between 0 and 1
    LagM  = (LagM-min(LagM(:)))./(max(LagM(:))-min(LagM(:)));
    LagMI = (LagMI-min(LagMI(:)))./(max(LagMI(:))-min(LagMI(:)));
    
    LagM  = smooth2a(LagM, nSmooth, nSmooth);
    LagMI = smooth2a(LagMI, nSmooth, nSmooth);
    
    %use eps to weight contributions froms lagM and lagMI
    for i=1:nE+1
        %sum
        LagE = eps(i).*LagM + (1-eps(i)).*LagMI;
    
        %find parameters minimizing the inference erorr
        minE = min(LagE(:));
        [xx, yy] = find(LagE == minE); xx = xx(1); yy = yy(1);

        %save optimal parameters 
        x0V(n, 2+i) = x0(yy);
        kV(n, 2+i)  = k(xx);

        %save optimal decoding parameters
        p1(n, 2+i) = res.p1(xx,yy,n);
        p0(n, 2+i) = res.p0(xx,yy,n);
    end

end


folder = '/outputFolder/'; 
filename = ['discreteGLM',estimator,num2str(nLevels),'levels',num2str(thetaH-thetaL),'sep',num2str(1000*h),'h.mat'];
paramObjTemp = struct();
for i=1:nE+1
    paramObjTemp = struct();    
    paramObjTemp.pVec = pVec;
    paramObjTemp.phi = 1;
    paramObjTemp.sigmaNoise = 1e-2;
    paramObjTemp.estimator = estimator;
    paramObjTemp.thetaL = thetaL;
    paramObjTemp.thetaH = thetaH;
    paramObjTemp.h = h;
    paramObjTemp.param  = estimatorParam;
    paramObjTemp.nBinsV = nLevels;
    paramObjTemp.k  = kV(:,2+i)' ;
    paramObjTemp.x0 = x0V(:,2+i)';
    paramObjTemp.p0 = p0(:,2+i)' ;
    paramObjTemp.p1 = p1(:,2+i)' ;
    paramObj.(['eps',num2str(round(100*eps(i)))]) = paramObjTemp;
end

if saveOut
    save([folder,filename],'paramObj');
end


end
