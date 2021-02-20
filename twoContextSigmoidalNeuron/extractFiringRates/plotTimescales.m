function plotTimescales(estimator)

if strcmp(estimator,'variance')
    thetaH = 3;
    thetaL = 1; 
    sep    = 2;
    ddt0   = 0;
    ddt1   = 20;
else
    thetaH = 2;
    thetaL = -2; 
    sep    = 4;
    ddt0   = 0;
    ddt1   = 10;
end

figure;hold on;
eps = [100,80,0];
h = [.005,.01,.015,.02,.025];
dt = floor(1./h)/2;
for j=1:numel(eps)
    

    for i=1:numel(h)
        rs = load(['outputFolder/res_',estimator,'_sep',num2str(sep),'_8levels_eps',num2str(eps(j)),'_h',num2str(100*h(i)),'.mat']);
        res{i} = rs.res;
        tUp{i}   = (dt(i)-ddt0):(dt(i)+ddt1);
        tDown{i} = (3*dt(i)-ddt0):(3*dt(i)+ddt1);
    end


    cmap = plasma(10);
    cmap = cmap([1,3,5,7,9],:);
    lw = 1.5*ones(1,5);

    
    for i=1:5
        rs = res{i};
        pL = (thetaH - rs.avg.mid.thetahatY)./(thetaH-thetaL);
        pH = 1 - pL;
        
        tLinterp = 0:.01:(numel(tDown{i})-1);
        tHinterp = 0:.01:(numel(tUp{i})-1);
        pLinterp = interp1(0:(numel(tDown{i})-1),pL(tDown{i}), tLinterp);
        pHinterp = interp1(0:(numel(tUp{i})-1),  pH(tUp{i}),   tHinterp);
        
        subplot(2,numel(eps)+1,j);hold on;
        plot(0:(numel(tUp{i})-1),pH(tUp{i}),'--','color',cmap(i,:),'linewidth',lw(i))
        xlim([0,(numel(tUp{i})-1)])
        xlabel('time after switch')
        if j==1
            ylabel('posterior')
        end
        set(gca,'fontsize',16)
        title(['switch up, \alpha = ',num2str((100-eps(j))./100)])

        plot(0:(numel(tDown{i})-1),pL(tDown{i}),'color',cmap(i,:),'linewidth',lw(i))
        xlim([0,(numel(tUp{i})-1)])
        set(gca,'fontsize',16)
        title(['\alpha = ',num2str((100-eps(j))./100)])
        
        iD = find(pLinterp>.5,1,'first');
        iH = find(pHinterp>.5,1,'first');
        dtDown(j,i) = tLinterp(iD);
        dtUp(j,i)   = tHinterp(iH);
    end
end

rs = res{i};
alpha = [0    0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000    0.9000    0.9200    0.9400    0.9600    0.9800    1.0000];
cmapR = viridis(126);
cmapR = cmapR(1:101,:);
cmapR = cmapR(round(100*alpha+1),:);
cmapR = cmapR([1,6,15],:);

ms = [10,13,16];
dt = floor(1./h);
for j=1:numel(eps)
    subplot(2,numel(eps)+1,numel(eps)+1+j);hold on;
    plot(dt,dtDown(j,:),'-v',dt,dtUp(j,:),'--^','markersize',ms(j),'linewidth',lw(1),'color',cmapR(j,:),'markerfacecolor','w')
    xlim([40,200])
    xlabel('period')
    if j==1
        ylabel('time for posterior to cross 1/2')
    end
    set(gca,'fontsize',16)
    
    subplot(2,numel(eps)+1,2*(numel(eps)+1));hold on;
    plot(dt,dtUp(j,:),'--^','markersize',ms(j),'linewidth',lw(1),'color',cmapR(j,:),'markerfacecolor','w')
    xlim([40,200])
    xlabel('period')
    title('upward switch')
    set(gca,'fontsize',16)
    
    subplot(2,numel(eps)+1,numel(eps)+1);hold on;
    plot(dt,dtDown(j,:),'-v','markersize',ms(j),'linewidth',lw(1),'color',cmapR(j,:),'markerfacecolor','w')
    xlim([40,200])
    xlabel('period')
    title('downward switch')
    set(gca,'fontsize',16)
    
end

set(gcf,'Position',[200 200 1400 600])

eps = 100*([0    0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000    0.9000    0.9200    0.9400    0.9600    0.9800    1.0000]);
h = [.005,.01,.015,.02,.025];
dt = floor(1./h)/2;
for j=1:numel(eps)
    

    for i=1:numel(h)
        rs = load(['outputFolder/res_',estimator,'_sep',num2str(sep),'_8levels_eps',num2str(eps(j)),'_h',num2str(100*h(i)),'.mat']);
        rs = rs.res;
        tUp   = (dt(i)-ddt0):(dt(i)+ddt1);
        tDown = (3*dt(i)-ddt0):(3*dt(i)+ddt1);

        pL = (thetaH - rs.avg.mid.thetahatY)./(thetaH-thetaL);
        pH = 1 - pL;
        
        tLinterp = 0:.01:(numel(tDown)-1);
        tHinterp = 0:.01:(numel(tUp)-1);
        pLinterp = interp1(0:(numel(tDown)-1),pL(tDown), tLinterp);
        pHinterp = interp1(0:(numel(tUp)-1),  pH(tUp),   tHinterp);
        
        iD = find(pLinterp>.5,1,'first');
        iH = find(pHinterp>.5,1,'first');
        dtDown(j,i) = tLinterp(iD);
        dtUp(j,i)   = tHinterp(iH);
    end
end

M = (dtDown./dtUp)';
figure;imagesc((dtDown./dtUp)');caxis([-.2,2.2]);colormap(redblue)
figure;hold on;
for i=1:5
    plot(1-eps./100,M(i,:),'-o','color',cmap(i,:),'markerfacecolor',cmap(i,:),'linewidth',lw(i),'markersize',12)
end