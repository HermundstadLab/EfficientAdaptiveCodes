function fitFiringRates(estimator)

eps = [0:10:90,92:2:100];
h = [.005,.01,.015,.02,.025];

cmap = viridis(126);
cmap = cmap(1:101,:);
cmapR = cmap(round(100-eps+1),:);

if strcmp(estimator,'variance')==1
    ind = 9;
else
    ind = 15;
end

figure;hold on;
for j=1:15
    tauU = zeros(1,5);
    tt   = zeros(1,5);
    for i=1:5
        [tauU(i),~,tt(i)] = fitFR(estimator,eps(j),h(i));
    end
    tbase = mean(tt);
    p = polyfit(1./h,(tauU-tbase),1);
    m(j) = p(1);
    b(j) = p(2);
    if j==ind
        subplot(1,4,1);hold on;
        plot(1./h,tauU-tbase,'-o','linewidth',2,'MarkerSize',10,'color',cmapR(j,:));
        plot(1./h,p(1)./h+p(2),'--','linewidth',2,'color',cmapR(j,:));xlim([40,200])
    end
    if ismember(j,[1,9,15])
        subplot(1,4,2);hold on;
        plot(1./h,p(1)./h+p(2),'-','linewidth',2,'color',cmapR(j,:));xlim([40,200])
    end
end
subplot(1,4,3);hold on;plot(1-eps./100,m,'-k')
subplot(1,4,4);hold on;plot(1-eps./100,b,'-k')
for i=1:15
    subplot(1,4,3);plot(1-eps(i)./100,m(i),'o','markerfacecolor',cmapR(i,:),'markeredgecolor','none','markersize',15)
    subplot(1,4,4);plot(1-eps(i)./100,b(i),'o','markerfacecolor',cmapR(i,:),'markeredgecolor','none','markersize',15)
end

subplot(1,4,1);set(gca,'fontsize',16);%xlim([40,200]);
subplot(1,4,2);set(gca,'fontsize',16)
subplot(1,4,3);set(gca,'fontsize',16)
subplot(1,4,4);set(gca,'fontsize',16)
set(gcf,'Position',[200 200 2000 400])




end

function [tauU,x,t1] = fitFR(estimator,eps,h)
r = load(['outputFolder/entropyCodes_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'.mat']);
r = r.rate;
r = mean(r,1);

dt = floor(1./h)/2;

baseline = mean(r((3*dt-5):(3*dt)));
r1 = r(dt+1:3*dt);
t1 = find(r1>(baseline+.5*(max(r1)-baseline)),1,'first');
r1 = r1(t1:40);

t  = (1:numel(r1))-1;

options = optimoptions(@fmincon,'Display','none');
fun = @(x)sum(((x(1).*exp(-(t)/x(2))+x(3))-r1).^2);

[x,~]  = fmincon(fun,[max(r1),1,baseline],[],[],[],[],[0,0,0],[100,100,100],[],options);
tauU = x(2)+t1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%























% % f = @rateErr;
% % alpha = 5;
% % A = [];
% % b = [];
% % Aeq = [];
% % beq = [];
% % lb = [3,-5,-5,3,-5,-5];
% % ub = [8,5,5,8,5,5];
% % 
% % copts=[];errs=[];
% % c0  = [alpha,.5,.5,alpha,.1,.1];
% % [c_opt,errors]  = fmincon(f,c0,A,b,Aeq,beq,lb,ub,[],options);



% % r1 = r(dt+1:3*dt);
% % [~,ii] = max(r1);
% % r1 = r1(2:end)-mean(r(floor(3*dt-5):(3*dt)));
% % r1 = r1./max(r1);
% % t  = (1:numel(r1))-1;
% % 
% % options = optimoptions(@fminunc,'Display','none');
% % fun = @(x)sum((exp(-t/x)-r1).^2);
% % [x,~] = fminunc(fun,1,options);
% % tauU = x+ii;
% % 
% % plot(t,r1,t,f1(t,x))
%plot(t,f1(t,tauU))

% % if strcmp(estimator,'variance')==1
% %     %[~,ii] = max(r1);
% %     r1 = r1(1:40)-mean(r(floor(3*dt-5):(3*dt))); %orig
% % else
% %     r1 = r1(1:40)-mean(r(floor(3*dt-5):(3*dt)));
% % end
% % r1 = r1./max(r1);
% % 
% % t  = (1:numel(r1))-1;
% % sp = 1;
% % 
% % params = fit(t',r1',g,'StartPoint',sp);
% % tauU = params.a;



% Switch down
% r1 = r([3*dt+1:4*dt,2:dt]); 
% r1 = r1(2:39)-mean(r(floor(dt-5):(dt)));
% r1 = r1./min(r1);
% 
% t  = (1:numel(r1))-1; 
% sp = 1;



%plot(t,r1,t,fun(t,params.a,params.b))
%plot(t,fun(t,params.a,params.b))
%plot(t,r1)
%plot(t,r1,t,fun(t,params.a))
%plot(t,fun(t,params.a))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f1 = @(x,a)(exp(-(x)/a));
% g   = fittype('exp(-(x)/a)');

% if strcmp(estimator,'variance')==1
%     [~,ii] = max(r1);
%     r1 = r1(2:40)-mean(r(floor(3*dt-5):(3*dt))); %orig
%     %r1 = r1(ii:40)-mean(r(floor(3*dt-5):(3*dt))); 
% else
%     r1 = r1(1:40)-mean(r(floor(3*dt-5):(3*dt)));
% end
% r1 = r1./max(r1);
% 
% t  = (1:numel(r1))-1;
% sp = 1;

%r1 = r1-mean(r(floor(3*dt-dt/4):(3*dt)));
%r1 = r1(1:40)-mean(r(floor(3*dt-5):(3*dt))); %orig

%r = load(['results/entropyCodes_avg/entropyCodes_avg_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'.mat']);

%fun = @(x,a,b,c)(a-b*exp(-x/c));
%g = fittype('a-b*exp(-x/c)');

%fun = @(x,a,b)(a*exp(-(x)/b));
%g   = fittype('a*exp(-(x)/b)');

%sp = [min(r1);1];

%r1 = r1-mean(r(floor(dt-dt/4):dt)); r1 =
%r1(1:39)-mean(r(floor(dt-5):dt));

% [~,ii] = max(r1);
% r1 = r1(ii:end);
% t = (1:numel(r1))-1;

%sp = [max(r1);1];




%params = fit(t',r1',g,'StartPoint',sp); 
%tauD = params.a;
%plot(t,r1)

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %times = ((dt+1):(dt+floor(dt/2)))';
% % % 
% % % %ddt = 1;
% % % %times = ((dt+ddt):3*dt-1)';
% % % times = ((dt+ddt):(dt+39))';
% % % t = (1:numel(times))';
% % % f = r(times)';
% % % 
% % % params = fit(t,f,g,'StartPoint',[[ones(size(t)), -exp(-t)]\f; 1]);
% % % tauU = params.c;
% % % %times = ((3*dt+1):(3*dt+floor(dt/2)))';
% % % 
% % % %subplot(1,2,1);hold on;plot(t+ddt,f,t+ddt,fun(t,params.a,params.b,params.c))
% % % %pause();
% % % 
% % % %ddt = 1;
% % % r = r([3*dt+1:4*dt,2:3*dt]);
% % % times = (ddt:39)';
% % % %times = [(3*dt+ddt):4*dt,2:dt]';
% % % t = (1:numel(times))';
% % % f = r(times)';
% % % 
% % % params = fit(t,f,g,'StartPoint',[[ones(size(t)), -exp(-t)]\f; 1]);
% % % tauD = params.c;
% % % 
% % % %subplot(1,2,2);hold on;plot(t+ddt,f,t+ddt,fun(t,params.a,params.b,params.c))
% % % %pause();

























% % r1 = r(dt+1:3*dt);
% % t = 1:numel(r1);
% % 
% % bins = linspace(1,numel(t),31);
% % [~,~,binid] = histcounts(t,bins);
% % 
% % f    = accumarray(binid',r1',[],@(x) mean(x));
% % tavg = accumarray(binid',t', [],@(x) mean(x));
% % 
% % params = fit(tavg,f,g,'StartPoint',[[ones(size(tavg)), -exp(-tavg)]\f; 1]);
% % tauU = params.c;
% % 
% % 
% % 
% % r1 = r([3*dt+1:4*dt,2:dt]);
% % t = 1:numel(r1);
% % 
% % bins = linspace(1,numel(t),31);
% % [~,~,binid] = histcounts(t,bins);
% % 
% % f    = accumarray(binid',r1',[],@(x) mean(x));
% % tavg = accumarray(binid',t', [],@(x) mean(x));
% % 
% % params = fit(tavg,f,g,'StartPoint',[[ones(size(tavg)), -exp(-tavg)]\f; 1]);
% % tauD = params.c;









% %,peakU,blU,peakD,blD,r
% 
% r = load(['results/entropyCodes_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'.mat']);
% r = r.rate;
% r = mean(r,1)';
% 
% 
% % tswitchU = dt+1;
% % tswitchD = 3*dt+1;
% dt = floor(1./h)/2;
% 
% blU = mean(r((3*dt-10):3*dt));
% blD = mean(r((dt-10):dt));
% 
% inds = (dt+1):2*dt;
% [~,ii] = max(abs(r(inds)-blU));
% peakU = r(inds(ii));
% 
% inds = (3*dt+1):4*dt;
% [~,ii] = max(abs(r((3*dt+1):4*dt)-blD));
% peakD = r(inds(ii));
% 
% tauU = find(r((dt+1):2*dt)<blU+(peakU-blU)/5,1,'first');
% tauD = find(r((3*dt+1):4*dt)>blD+(peakD-blD)/5,1,'first');


% g = fittype('a-b*exp(-x/c)');
% fun = @(x,a,b,c)(a-b*exp(-x/c));
% 
% r = load(['results/entropyCodes_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'.mat']);
% r = r.rate;
% r1 = load(['results/entropyCodes_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'_v2.mat']);
% r1 = r1.rate;
% r2 = load(['results/entropyCodes_',estimator,'_eps',num2str(eps),'_h',num2str(100*h),'_v3.mat']);
% r2 = r2.rate;
% % kk = randperm(5000);
% % kk = kk(1:2500);
% % r = mean(r(kk,:),1)';
% r = mean([r;r1;r2],1)';
% 
% dt = floor(1./h)/2;
% t0 = 1:4*dt;
% 
% if strcmp(estimator,'variance')==1
%     tswitchU = dt+3;
%     tswitchD = 3*dt+4;
% else
%     tswitchU = dt+2;
%     tswitchD = 3*dt+2;
% end
% 
% % tswitchU = dt+1;
% % tswitchD = 3*dt+1;
% 
% 
% delta = 2*dt-10;
% %delta = 20;
% t = (1:(delta+1))';
% 
% %switch up
% %find peak
% inds = (tswitchU:(tswitchU+floor(delta/2)))';
% [~,iiU] = max(abs(r(inds)-r(inds(end))));
% 
% tprobeU = (inds(iiU):inds(iiU)+delta)';
% frU = r(tprobeU);
% fitParamsUp = fit(t,frU,g,'StartPoint',[[ones(size(t)), -exp(-t)]\frU; 1]);
% 
% %switch down
% inds = (tswitchD:(tswitchD+floor(delta/2)))';
% [~,iiD] = max(abs(r(inds)-r(inds(end))));
% 
% if inds(iiD)+delta>4*dt
%     tprobeD = [inds(iiD):4*dt,2:(2+ delta-(4*dt-inds(iiD)+1) )]';
%     tplot = 1:numel(inds(iiD):4*dt);
% else
%     tprobeD = (inds(iiD):inds(iiD)+delta)';
%     tplot = 1:numel(tprobeD);
% end
% frD = r(tprobeD);
% fitParamsDown = fit(t,frD,g,'StartPoint',[[ones(size(t)), -exp(-t)]\frD; 1]);
% 
% if strcmp(estimator,'variance')==1
%     tauU = fitParamsUp.c+2;
%     tauD = fitParamsDown.c+3;
% else
%     tauU = fitParamsUp.c+1;
%     tauD = fitParamsDown.c+1;
% end
% 
% figure;hold on;
% plot(t0,r);
% plot(tprobeU,frU,tprobeU,fun(t,fitParamsUp.a,fitParamsUp.b,fitParamsUp.c));
% plot(tprobeD(tplot),frD(tplot),tprobeD(tplot),fun(t(tplot),fitParamsDown.a,fitParamsDown.b,fitParamsDown.c));






















%upward switches
%inds = dt+1:2*dt;
%plot(1:numel(inds),r(inds),tprobeU-numel(inds),fun(t,fitParamsUp.a,fitParamsUp.b,fitParamsUp.c));

% indsU = (dt+1):(2*dt);
% indsD = (3*dt+1):(4*dt);
% tauU = find(r(indsU) < fitParamsUp.a-.1*fitParamsUp.b,1,'first');
% tauD = find(r(indsD) > fitParamsDown.a-.1*fitParamsDown.b,1,'first');

% figure;hold on;
% plot(t0,r);
% plot(tprobeU,frU,tprobeU,fun(t,fitParamsUp.a,fitParamsUp.b,fitParamsUp.c),indsU(tauU),r(indsU(tauU)),'xk');
% plot(tprobeD(tplot),frD(tplot),tprobeD(tplot),fun(t(tplot),fitParamsDown.a,fitParamsDown.b,fitParamsDown.c),indsD(tauD),r(indsD(tauD)),'xk');