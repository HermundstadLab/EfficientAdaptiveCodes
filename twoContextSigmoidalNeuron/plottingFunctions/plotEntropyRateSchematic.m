function plotEntropyRateSchematic()

folder = 'outputFolder/';
paramObj = load([folder,'GLMparams/discreteGLMmean8levels4sep1h.mat']);paramObj = paramObj.paramObj.eps100;
rateM100 = load([folder,'runInference/results/entropyCodes_mean_eps100_h1.mat']);rateM100 = mean(rateM100.rate);rateM100(1) = rateM100(end);
rateM80  = load([folder,'runInference/results/entropyCodes_mean_eps80_h1.mat']); rateM80  = mean(rateM80.rate);rateM80(1) = rateM80(end);
rateM0   = load([folder,'runInference/results/entropyCodes_mean_eps0_h1.mat']);  rateM0   = mean(rateM0.rate);rateM0(1) = rateM0(end);

rateV100 = load([folder,'runInference/results/entropyCodes_variance_eps100_h1.mat']);rateV100 = mean(rateV100.rate);rateV100(1) = rateV100(end);
rateV80  = load([folder,'runInference/results/entropyCodes_variance_eps80_h1.mat']); rateV80  = mean(rateV80.rate);rateV80(1) = rateV80(end);
rateV0   = load([folder,'runInference/results/entropyCodes_variance_eps0_h1.mat']);  rateV0   = mean(rateV0.rate);rateV0(1) = rateV0(end);



rM100 = load([folder,'runInference/results/res_mean_sep4_8levels_eps100_h1.mat']);rM100 = rM100.res.avg.mid.y;
rM80  = load([folder,'runInference/results/res_mean_sep4_8levels_eps80_h1.mat']); rM80  = rM80.res.avg.mid.y;
rM0   = load([folder,'runInference/results/res_mean_sep4_8levels_eps0_h1.mat']);  rM0   = rM0.res.avg.mid.y;

rV100 = load([folder,'runInference/results/res_variance_sep2_8levels_eps100_h1.mat']);rV100 = rV100.res.avg.mid.y;
rV80  = load([folder,'runInference/results/res_variance_sep2_8levels_eps80_h1.mat']); rV80  = rV80.res.avg.mid.y;
rV0   = load([folder,'runInference/results/res_variance_sep2_8levels_eps0_h1.mat']);  rV0   = rV0.res.avg.mid.y;


f = @(x,k,x0) ( 1 ./ (1 + exp(-k * (x - x0) ) ));
ybins = linspace(0,1,9);

x = -6:.01:6;
ind = 80;
mu = paramObj.pVec(ind)*(-2) + (1-paramObj.pVec(ind))*2;
k = paramObj.k(ind);x0 = paramObj.x0(ind);


xs = normrnd(mu,1,[1,100000]);
ys = f(xs,k,x0);
spk = histcounts(ys,ybins);
spk = spk./sum(spk);
[~,ii] = sort(spk,'descend');
kk = 0:7;

xs = normrnd(-2,1,[1,100000]);
ys = f(xs,k,x0);
spkT = histcounts(ys,ybins);
spkT = spkT./sum(spkT);

figure;
subplot(4,3,1);hold on;
plot(x,gaussian(x,mu,1));plot(x,gaussian(x,-2,1));set(gca,'fontsize',16);
subplot(4,3,4);plot(x,f(x,k,x0));set(gca,'fontsize',16);
subplot(4,3,7);bar(ii,kk);set(gca,'fontsize',16);

subplot(4,3,2);bar(spk);set(gca,'fontsize',16);
subplot(4,3,5);bar(spk(ii));set(gca,'fontsize',16);
subplot(4,3,8);bar(spkT);set(gca,'fontsize',16);
subplot(4,3,11);bar(spkT(ii));set(gca,'fontsize',16);

subplot(4,3,3);hold on;plot(rM100);plot(rM80);plot(rM0);ylim([0,7]);set(gca,'fontsize',16);
subplot(4,3,6);hold on;plot(rateM100);plot(rateM80);plot(rateM0);ylim([0,7]);set(gca,'fontsize',16);
subplot(4,3,9);hold on;plot(rV100);plot(rV80);plot(rV0);ylim([0,7]);set(gca,'fontsize',16);
subplot(4,3,12);hold on;plot(rateV100);plot(rateV80);plot(rateV0);ylim([0,7]);set(gca,'fontsize',16);


set(gcf,'Position',[200 200 1200 800])