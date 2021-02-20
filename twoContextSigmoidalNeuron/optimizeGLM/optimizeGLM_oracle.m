function [res,paramObj] = optimizeGLM_oracle(estimator,nLevels,sig,h,saveOut)

if nargin<5
    saveOut = false;
end

%training dataset
xCore = randn(1e5, 1);

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
else
    error('unrecognized estimator');
end
    
for i=1:2

    %environment parameters
    if strcmp(estimator,'mean')==1
        thetaL = -sig;
        thetaH = sig;
        sigma = 1;
        if i<1.5
            xM = xCore * sigma + thetaL;
        elseif i>1.5
            xM = xCore * sigma + thetaH;
        end
    elseif strcmp(estimator,'variance')==1
        thetaL = 1;
        thetaH = sig;
        mu = 0;
        if i<1.5
            xM = xCore * thetaL + mu;
        elseif i>1.5
            xM = xCore * thetaH + mu;
        end
    else
        error('unrecognized estimator');
    end

    errMRec = zeros(N, N);
    resM = zeros(N, N);
    p1 = zeros(N, N);
    p0 = zeros(N, N);

    errMRec_flip = zeros(N, N);
    resM_flip = zeros(N, N);
    p1_flip = zeros(N, N);
    p0_flip = zeros(N, N);

    %set random number generator (depends on Matlab version)
    rng('default');


    %optimize
    for n = 1:N
        for m = 1:N
            %disp([m,n])
            rV = encode(xM, phi, k(n), x0(m), paramsEncode);
            [xHat,paramsOLE] = decodeOLEFromDiscreteGLM(rV, xM);

            errMRec(n, m) = mean((xM - xHat).^2);
            p1(n, m) = paramsOLE(1);
            p0(n, m) = paramsOLE(2);
            resM(n, m) = mean(rV);

            if strcmp(estimator,'variance')==1
                rV = encode(xM, phi, k(n), -x0(m), paramsEncode);
                [xHat,paramsOLE] = decodeOLEFromDiscreteGLM(rV, xM);

                errMRec_flip(n, m) = mean((xM - xHat).^2);
                p1_flip(n, m) = paramsOLE(1);
                p0_flip(n, m) = paramsOLE(2);
                resM_flip(n, m) = mean(rV);
            end
        end  
    end


    res.errMRec = errMRec;
    res.resM = resM;
    res.p1 = p1;
    res.p0 = p0;

    if strcmp(estimator,'variance')==1
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



    nSmooth = 12; 
    LagM  = res.errMRec;

    %recale error to be between 0 and 1
    LagM  = (LagM-min(LagM(:)))./(max(LagM(:))-min(LagM(:)));
    LagM  = smooth2a(LagM, nSmooth, nSmooth);


    %find parameters minimizing the inference erorr
    minE = min(LagM(:));
    [xx, yy] = find(LagM == minE); xx = xx(1); yy = yy(1);

    %save optimal parameters 
    x0V(i) = x0(yy);
    kV(i)  = k(xx);

    %save optimal decoding parameters
    p1V(i) = res.p1(xx,yy);
    p0V(i) = res.p0(xx,yy);
end


%toc;

folder = '/outputFolder/'; 
filename = ['discreteOracleGLM',estimator,num2str(nLevels),'levels',num2str(thetaH-thetaL),'sep',num2str(1000*h),'h.mat'];

paramObj = struct();    
paramObj.phi = 1;
paramObj.sigmaNoise = 1e-2;
paramObj.estimator = estimator;
paramObj.thetaL = thetaL;
paramObj.thetaH = thetaH;
paramObj.h = h;
paramObj.param  = estimatorParam;
paramObj.nBinsV = nLevels;
paramObj.k  = kV;
paramObj.x0 = x0V;
paramObj.p0 = p0V;
paramObj.p1 = p1V;


if saveOut
    save([folder,filename],'paramObj');
end


end
