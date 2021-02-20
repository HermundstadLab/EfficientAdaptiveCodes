function M = plotNonlinearities(estimator,nLevels,eps,h)

if strcmp(estimator,'variance')==1
    thetaH = 3;
    thetaL = 1;
    x = -6:.01:6;
else
    thetaH = 2;
    thetaL = -2;
    x = -5:.01:5;
end
sep = thetaH-thetaL;

folder = 'outputFolder/'; 
filename = ['discreteGLM',estimator,num2str(nLevels),'levels',num2str(sep),'sep',num2str(100*h),'h.mat'];
res = load([folder,filename]);   
paramObj = res.paramObj.(['eps',num2str(eps)]);

phi = 1;
params = [.01,nLevels];




R = [];
for i=1:20
    rate = entropyCoding_posterior(estimator,nLevels,eps,h);
    R = cat(3,R,rate);
end
rate = mean(R,3);
M = zeros(numel(paramObj.pVec),numel(x));
ybin = linspace(0,1,nLevels+1);
ybin = ybin(2:end-1);

for i=1:numel(paramObj.pVec)
    k  = paramObj.k(i);
    x0 = paramObj.x0(i);

    
    [~,r] = encodeWithDiscrete(x, phi, k, x0, params);
    jmin = 1;
    for j=1:nLevels-1
        [~,ii] = min((ybin(j)-r).^2);
        M(i,jmin:ii) = rate(i,j);
        jmin=ii+1;
    end
    M(i,jmin:end) = rate(i,nLevels);
end

figure;imagesc(x,paramObj.pVec,M);colormap(flipud(gray));