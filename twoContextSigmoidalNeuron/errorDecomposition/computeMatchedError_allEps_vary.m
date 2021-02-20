function res = computeMatchedError_allEps_varyh(h)


nLevels = 8;
estimators = {'mean','variance'};
sep = [4,2];

eps = 0:.02:1;
inds = [1:5:46,47:51];
eps = eps(inds);

for j=1:numel(estimators)
    [err,excess,errF,errO,FRprops] = matchedError(estimators{j},sep(j),nLevels,h,100*eps);

    res.(estimators{j}).eps     = eps;
    res.(estimators{j}).nLevels = nLevels;

    res.(estimators{j}).totalError(1,1:numel(eps))   = err(1,:);
    res.(estimators{j}).matchedError(1,1:numel(eps)) = err(2,:);
    res.(estimators{j}).excessError(1,1:numel(eps))  = err(3,:);
    
    res.(estimators{j}).excessErrorAvg(1,1:numel(eps)) = excess(2,:);
    res.(estimators{j}).excessTime(1,1:numel(eps))     = excess(1,:);
    
    res.(estimators{j}).dtSwitchUp   = excess(3,:);      
    res.(estimators{j}).dtSwitchDown = excess(4,:); 
    
    res.(estimators{j}).fixedError  = errF;
    res.(estimators{j}).oracleError = errO;
    
    res.(estimators{j}).blHigh   = FRprops(1,:);      
    res.(estimators{j}).blLow    = FRprops(2,:);
    res.(estimators{j}).peakHigh = FRprops(3,:);      
    res.(estimators{j}).peakLow  = FRprops(4,:);
            
    
end

save(['matchedError_allEps_h',num2str(1000*h),'temp.mat'],'res');


end

function [errRec,excess,errF,errO,FRprops] = matchedError(estimator,sep,nLevels,h,eps)


if strcmp(estimator,'mean')==1
    thetaH = sep/2;
    thetaL = -sep/2;
    thresh = .5;
else
    thetaH = sep+1;
    thetaL = 1;
    thresh = .1;
end 




%--------------- compute task error for eps strategies -------------------%

folder = '/outputFolder/'; 

for i=1:numel(eps)
    filename = ['resNL_',estimator,'_sep',num2str(sep),'_',num2str(nLevels),'levels','_eps',num2str(eps(i)),'_h',num2str(100*h),'revised.mat'];
    rI = load([folder,filename]);   
    res.(['eps',num2str(eps(i))]) = rI.res;
end

dt = floor(1./h)/2;




for i=1:numel(eps)
    %all error
    resTemp = res.(['eps',num2str(eps(i))]);
    errRec(1,i) = mean(resTemp.full.errRec);
    pLow = (thetaH-resTemp.avg.mid.thetahatY)./(thetaH - thetaL);
    pHigh = 1-pLow;
    
    
    %switch down
    t0 = dt; 
    indsU = (t0+1):(2*dt);
    
    
    
    tU  = [t0, t0 + (1:find( pHigh(indsU)>(1-thresh),1,'first'))];
    
    
    t0 = 3*dt; 
    indsD = (t0+1):(4*dt); %start with previous timestep  
    pLow = (thetaH-resTemp.avg.mid.thetahatY)./(thetaH - thetaL);
    
    tD = [t0, t0 + (1:find( pLow(indsD)>(1-thresh),1,'first'))];
    
    tTransient = [tU,tD];
    tStable = 1:4*dt;
    tStable(tTransient) = []; 
 
    %matched-state error
    errRec(2,i) = sum(resTemp.avg.mid.errRec(tStable))./(4*dt);
    
    %excess error
    errRec(3,i) = sum(resTemp.avg.mid.errRec(tTransient))./(4*dt);  
    
    pLow = (thetaH-resTemp.avg.end.thetahatY(1:2*dt))./(thetaH - thetaL);
    pHigh = 1-(thetaH-resTemp.avg.end.thetahatY((2*dt+1):4*dt))./(thetaH - thetaL);
    
    ddt = floor(dt/2);

    boundLow  = mean(pLow((end-ddt):end));
    boundHigh = mean(pHigh((end-ddt):end));
    
    sdLow  = std(pLow((end-ddt):end));
    sdHigh = std(pHigh((end-ddt):end));

    
    n=5; 
    dtU = [t0, t0 + (1: (find( pHigh<(boundHigh - n*sdHigh) | pHigh>(boundHigh + n*sdHigh) ,1,'last')+1) )];
    dtD = [t0, t0 + (1: (find( pLow <(boundLow -  n*sdLow ) | pLow >(boundLow  + n*sdLow ) ,1,'last')+1) )];
    
       
    excess(1,i)   = numel(tTransient)./(4*dt);
    excess(2,i)   = mean(resTemp.avg.mid.errRec(tTransient));
    excess(3,i)   = numel(dtU);
    excess(4,i)   = numel(dtD);
    
    
    
    %load entropy rates
    if h~=0.01
        if ismember(i,[1,9,15])
            r = load(['results/entropyCodes_',estimator,'_eps',num2str(eps(i)),'_h',num2str(100*h),'.mat']);
            r = r.rate;
            r = mean(r,1)';

            frHigh = r((dt+1):3*dt);
            frLow  = [r((3*dt+1):4*dt);r(2:dt)];

            frH = mean(frHigh((end-ddt):end));
            frL = mean(frLow((end-ddt):end));
            [~,ii] = max(abs(frHigh-frH));
            frPeakH = frHigh(ii);
            [~,ii] = max(abs(frLow-frL));
            frPeakL = frLow(ii);

            FRprops(1,i) = frH;
            FRprops(2,i) = frL;
            FRprops(3,i) = frPeakH;
            FRprops(4,i) = frPeakL;
        end
    else
        r = load(['results/entropyCodes_',estimator,'_eps',num2str(eps(i)),'_h',num2str(100*h),'.mat']);
        r = r.rate;
        r = mean(r,1)';

        frHigh = r((dt+1):3*dt);
        frLow  = [r((3*dt+1):4*dt);r(2:dt)];

        frH = mean(frHigh((end-ddt):end));
        frL = mean(frLow((end-ddt):end));
        [~,ii] = max(abs(frHigh-frH));
        frPeakH = frHigh(ii);
        [~,ii] = max(abs(frLow-frL));
        frPeakL = frLow(ii);

        FRprops(1,i) = frH;
        FRprops(2,i) = frL;
        FRprops(3,i) = frPeakH;
        FRprops(4,i) = frPeakL;
    end
end
    











folder = '/Users/hermundstada/Dropbox/adaptiveInferenceSimulations4/GLMparams/'; 

filename = ['discreteFixedGLM',estimator,num2str(nLevels),'levels',num2str(sep),'sep',num2str(1000*h),'h.mat'];
pobj = load([folder,filename]);   
paramObjFixed = pobj.paramObj;

filename = ['discreteOracleGLM',estimator,num2str(nLevels),'levels',num2str(sep),'sep',num2str(1000*h),'h.mat'];
pobj = load([folder,filename]);   
paramObjOracle = pobj.paramObj;

x = res.eps0.full.x;
theta = res.eps0.full.theta;


%----------------- compute task error for fixed code----------------------%
[r,~] = encodeWithDiscrete(x, 1, paramObjFixed.k, paramObjFixed.x0, [.01,nLevels]);
xHat   = r * paramObjFixed.p1 + paramObjFixed.p0;

errF = mean((x-xHat).^2);


%---------------- compute task error for oracle code ---------------------%
xL = x(theta==min(theta));
[rL,~] = encodeWithDiscrete(xL, 1, paramObjOracle.k(1), paramObjOracle.x0(1), [.01,nLevels]);
xHatL   = rL * paramObjOracle.p1(1) + paramObjOracle.p0(1);

xH = x(theta==max(theta));
[rH,~] = encodeWithDiscrete(xH, 1, paramObjOracle.k(2), paramObjOracle.x0(2), [.01,nLevels]);
xHatH   = rH * paramObjOracle.p1(2) + paramObjOracle.p0(2);

xHat = [xHatL;xHatH];
xx   = [xL;xH];

errO = mean((xx-xHat).^2);






end

