function hA = plotMatchedError()



res = load('matchedError_allEps_h5.mat');  res05 = res.res;
res = load('matchedError_allEps_h10.mat'); res10 = res.res;
res = load('matchedError_allEps_h15.mat'); res15 = res.res;
res = load('matchedError_allEps_h20.mat'); res20 = res.res;
res = load('matchedError_allEps_h25.mat'); res25 = res.res;

eps = 1-res05.mean.eps;
cmap = viridis(126);
cmap = cmap(1:101,:);
cmapR = cmap(round(100*eps+1),:);


%all five hazard rates
iF = 5;
res = {res05,res15,res20,res25,res10};
c = ones(5,3);
m = [.6,.5,.4,.3,0]'+.2;
m = repmat(m,[1,3]);
c = c.*m;
lw = .5*ones(1,5);lw(iF) = 1.5;


estimators = {'variance','mean'};


figure;hold on;
r = res{5};

for k=1:numel(estimators)

    subplot(2,3,3*(k-1)+1);hold on;
    r0 = res{5};
    e0 = r0.(estimators{k}).excessErrorAvg(end)*r0.(estimators{k}).excessTime(end);
    if k==1
        x = .35:.005:.85;
        y = .13:.0005:.18;
        eMax = ceil(1.2*x(end)*y(end));
        eMin = floor(.8*x(1)*y(1));
        de = (eMax-eMin)./70;
        ll = [fliplr(e0:-de:eMin),(e0+de):de:eMax];
    else
        x = 0:.03:3;
        y = .02:.0005:.08;
        eMax = ceil(1.2*x(end)*y(end));
        eMin = floor(.8*x(1)*y(1));
        de = (eMax-eMin)./50;
        ll = [fliplr(e0:-de:eMin),(e0+de):de:eMax];
    end

    for j=1:numel(ll)
        plot(x,ll(j)./x,'linewidth',.5,'color',[.8,.8,.8])
    end

    plot(r.(estimators{k}).excessErrorAvg,r.(estimators{k}).excessTime,'-','color',c(iF,:),'linewidth',lw(iF));

    for j=1:numel(eps)
        plot(r.(estimators{k}).excessErrorAvg(j),r.(estimators{k}).excessTime(j),'o','markeredgecolor','none','markerfacecolor',cmapR(j,:),'markersize',15)
    end

    set(gca,'fontsize',16)
    xlabel('excess error')
    ylabel('excess time')
    if k==1
        xlim([.35,.85])
        ylim([.13,.18])
    else
        xlim([0,3])
        ylim([.02,.08])
    end


    subplot(2,3,3*(k-1)+2);hold on;
    r0 = res{5};
    e0 = r0.(estimators{k}).matchedError(end)+r0.(estimators{k}).excessError(end);
    if k==1
        ll = [fliplr(e0:-.02:0),(e0+.02):.02:.6];
    else
        ll = [fliplr(e0:-.025:0),(e0+.025):.025:.5];
    end

    for j=1:numel(ll)
        plot([0,ll(j)],[ll(j),0],'linewidth',.5,'color',[.8,.8,.8])
    end

    plot(r.(estimators{k}).matchedError,r.(estimators{k}).excessError,'-','color',c(iF,:),'linewidth',lw(iF));

    for j=1:numel(eps)
        plot(r.(estimators{k}).matchedError(j),r.(estimators{k}).excessError(j),'o','markeredgecolor','none','markerfacecolor',cmapR(j,:),'markersize',15)
    end
    [~,ii] = min(r.(estimators{k}).totalError);
    plot(r.(estimators{k}).matchedError(ii),r.(estimators{k}).excessError(ii),'xk','markersize',12)
    set(gca,'fontsize',16)
    xlabel('matched error')
    ylabel('mismatched error')
    if k==1
        xlim([.18,.3])
        ylim([.05,.14])
    else
        xlim([.04,.09])
        ylim([0,.2])
    end


    subplot(2,3,3*(k-1)+3);hold on;
    plot(eps,r.(estimators{k}).totalError,'-','color',c(iF,:),'linewidth',lw(iF));

    for j=1:numel(eps)
        plot(eps(j),r.(estimators{k}).totalError(j),'o','markeredgecolor','none','markerfacecolor',cmapR(j,:),'markersize',15)
    end

    [~,ii] = min(r.(estimators{k}).totalError);
    plot(eps(ii),r.(estimators{k}).totalError(ii),'xk','markersize',12)
    set(gca,'fontsize',16)
    xlabel('bias')
    ylabel('total error')
end
set(gcf,'Position',[200 200 1400 800])




cmap = plasma(10);
cmap = cmap([1,3,5,7,9],:);
lw = 1.5*ones(1,5);
res = {res05,res10,res15,res20,res25};
figure;hold on;
for i=1:numel(res)
    r = res{i};
    
    for k=1:numel(estimators)
        
        subplot(2,3,3*(k-1)+1);hold on;
        
        plot(r.(estimators{k}).excessErrorAvg,r.(estimators{k}).excessTime,'-o','color',cmap(i,:),'linewidth',lw(i),'markeredgecolor','none','markerfacecolor',cmap(i,:),'markersize',15);
        set(gca,'fontsize',16)
        xlabel('excess error')
        ylabel('excess time')
        
        
        subplot(2,3,3*(k-1)+2);hold on;
        
        plot(r.(estimators{k}).matchedError,r.(estimators{k}).excessError,'-o','color',cmap(i,:),'linewidth',lw(i),'markeredgecolor','none','markerfacecolor',cmap(i,:),'markersize',15);
        [~,ii] = min(r.(estimators{k}).totalError);
        plot(r.(estimators{k}).matchedError(ii),r.(estimators{k}).excessError(ii),'xk','markersize',12)
        set(gca,'fontsize',16)
        xlabel('matched error')
        ylabel('mismatched error')
        
       
        subplot(2,3,3*(k-1)+3);hold on;
        plot(eps,r.(estimators{k}).totalError,'-o','color',cmap(i,:),'linewidth',lw(i),'markeredgecolor','none','markerfacecolor',cmap(i,:),'markersize',15);
        [~,ii] = min(r.(estimators{k}).totalError);
        plot(eps(ii),r.(estimators{k}).totalError(ii),'xk','markersize',12)
        set(gca,'fontsize',16)
        xlabel('bias')
        ylabel('total error')
    end
    
    
    
end
set(gcf,'Position',[200 200 1400 800])

h = [.005,.01,.015,.02,.025];
dt = floor(1./h);

figure;hold on;
for i=1:numel(res)
    r = res{i};
    
    for k=1:numel(estimators)
        
        eF = r.(estimators{k}).fixedError;
        eO = r.(estimators{k}).oracleError;
        
        subplot(2,3,3*(k-1)+1);hold on;           
        x = r.(estimators{k}).excessErrorAvg;
        y = r.(estimators{k}).excessTime;
        
        x = (x-x(end))./(eF-eO);
        y = (y-y(end));   

        plot(x,y,'-o','color',cmap(i,:),'linewidth',lw(i),'markeredgecolor','none','markerfacecolor',cmap(i,:),'markersize',15);
        set(gca,'fontsize',16)
        xlabel('excess error')
        ylabel('excess time')
        
        
        subplot(2,3,3*(k-1)+2);hold on;        
        x = r.(estimators{k}).matchedError;
        y = r.(estimators{k}).excessError;
        x = (x-x(end))./(eF-eO);
        y = (y-y(end))./(eF-eO);
        
        plot(x,y,'-o','color',cmap(i,:),'linewidth',lw(i),'markeredgecolor','none','markerfacecolor',cmap(i,:),'markersize',15);
        set(gca,'fontsize',16)
        xlabel('matched error')
        ylabel('mismatched error')
        
       
        subplot(2,3,3*(k-1)+3);hold on;
        x = eps;
        y = r.(estimators{k}).totalError;
        y = (y-y(end))./(eF-eO);
        
        plot(x,y,'-o','color',cmap(i,:),'linewidth',lw(i),'markeredgecolor','none','markerfacecolor',cmap(i,:),'markersize',15);
        set(gca,'fontsize',16)
        xlabel('bias')
        ylabel('total error')
    end
    
    
    
end
set(gcf,'Position',[200 200 1400 800])
