function rate = entropyCoding(estimator,nLevels,eps,h,sep)

folder = 'outputFolder/'; 
filename = ['res_',estimator,'_sep',num2str(sep),'_',num2str(nLevels),'levels','_eps',num2str(eps),'_h',num2str(100*h),'.mat'];
res = load([folder,filename]);   
res = res.res;

%------------------------- full --------------------------%

nTrials = 5000;
nT = 2*floor(1./h);
inds = [nT/4+1:nT,1:nT/4];
xM = reshape(res.full.x,[nT,nTrials])';xM = xM(:,inds);
thetaM = res.avg.end.theta';thetaM = thetaM(:,inds);
thetahatM = reshape(res.full.thetahatY,[nT,nTrials])';thetahatM = thetahatM(:,inds);
x0M = reshape(res.full.encodingParams(:,2),[nT,nTrials])';x0M = x0M(:,inds);
kM  = reshape(res.full.encodingParams(:,1),[nT,nTrials])';kM = kM(:,inds);


thetaH = max(thetaM(:));
thetaL = min(thetaM(:));
pLow = linspace(0,1,101);
thetahatVec = pLow.*thetaL + (1-pLow).*thetaH;

nS = 50000;
dists = zeros(numel(thetahatVec),nS);
for i=1:numel(thetahatVec)
    if strcmp(estimator,'mean')
        xhat = normrnd(thetahatVec(i),1,[1,nS]);
    else
        xhat = normrnd(0,thetahatVec(i),[1,nS]);
    end
    dists(i,:) = xhat;
end


phi = 1;
params = [.01,nLevels];
bins = (0:(nLevels))-.5;
rate = zeros(nTrials,nT);
spikes = 0:(nLevels-1);
parfor i=1:nTrials
    for j=2:nT
        k  = kM(i,j);
        x0 = x0M(i,j);

        thetahat = thetahatM(i,j-1);
        
        [~,ind] = min(abs(thetahat-thetahatVec));
        xhat = dists(ind,:);
        yhat = encodeWithDiscrete(xhat, phi, k, x0, params);
        chat = histcounts(yhat,bins)./nS;
        [~,inds] = sort(chat,'descend');
        s = zeros(1,nLevels);
        s(inds) = spikes;

        y = encodeWithDiscrete(xM(i,j), phi, k, x0, params);
        iBin = find(histcounts(y,bins));
        rate(i,j) = s(iBin);
    end
end

filename = [folder,'entropyCodes_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'.mat'];
save(filename,'rate')





