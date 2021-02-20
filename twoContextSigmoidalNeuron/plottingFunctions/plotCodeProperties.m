function [tA,eSS,eMax] = plotCodeProperties(thetaL,thetaH,h,estimator,figNums,saveFig)

if nargin<6
    saveFig = false;
end

sep = thetaH-thetaL;
cmap = viridis(100);
cR = cmap(10,:);
cM = cmap(50,:);
cI = cmap(80,:);

lw = 2;
nLevels = 8;
    
%define aux functions
g = @(x,mu,sig) (1./sqrt(2.*pi.*sig.^2).*exp(-(x-mu).^2./(2.*sig.^2)));
f = @(x,k,x0) (1 ./ (1 + exp(-k * (x - x0) ) ));
x = -6:.01:6;

%load data
saveFolder = ['figs/',estimator,'/'];
folder = 'outputFolder/'; 

filename = ['res_',estimator,'_sep',num2str(sep),'_',num2str(nLevels),'levels','_eps',num2str(0),'_h',num2str(100*h),'.mat'];
rI = load([folder,filename]);   
resI = rI.res;

filename = ['res_',estimator,'_sep',num2str(sep),'_',num2str(nLevels),'levels','_eps',num2str(80),'_h',num2str(100*h),'.mat'];
rM = load([folder,filename]);  
resM = rM.res;

filename = ['res_',estimator,'_sep',num2str(sep),'_',num2str(nLevels),'levels','_eps',num2str(100),'_h',num2str(100*h),'.mat'];
rR = load([folder,filename]);  
resR = rR.res;



if strcmp(estimator,'mean')==1
    mu  = [thetaL,thetaH];
    sig = [1,1];
    
    inds = [51:200,1:50];
    mi = MI(resI.full.x,resI.full.y,8,h);
    resI.avg.mid.MI = mi(inds);
    resI.avg.end.MI = mi;

    mi = MI(resR.full.x,resR.full.y,8,h);
    resR.avg.mid.MI = mi(inds);
    resR.avg.end.MI = mi;

    mi = MI(resM.full.x,resM.full.y,8,h);
    resM.avg.mid.MI = mi(inds);
    resM.avg.end.MI = mi;
    
elseif strcmp(estimator,'variance')==1
    sig = [thetaL,thetaH];
    mu  = [0,0];
    
    inds = [51:200,1:50];
    mi = MI(resI.full.x,resI.full.y,8,h);
    resI.avg.mid.MI = mi(inds);
    resI.avg.end.MI = mi;

    mi = MI(resR.full.x,resR.full.y,8,h);
    resR.avg.mid.MI = mi(inds);
    resR.avg.end.MI = mi;

    mi = MI(resM.full.x,resM.full.y,8,h);
    resM.avg.mid.MI = mi(inds);
    resM.avg.end.MI = mi;
else
    error('unrecognized estimator')
end
    
%------------------------- FIGURE 1 --------------------------------------%
if ismember(1,figNums)  
    
    figure;hold on;
    if strcmp(estimator,'mean')==1
        xlims = [48,60];
        ylimInfo = [0,3];
        ylimRec  = [0,5];
    else
        xlims = [45,75];
        ylimInfo = [1,3];
        ylimRec  = [0,3];
    end 

    %----------------plot stimulus distribution---------------------------%
    subplot(3,3,1);hold on;
    plot([mu(1),mu(1)],[0,.5],'--','color',[.6,.6,.6],'linewidth',lw)
    plot([mu(2),mu(2)],[0,.5],'--','color',[.6,.6,.6],'linewidth',lw)
    plot(x,g(x,mu(1),sig(1)),'color',[.6,.6,.6],'linewidth',lw)
    plot(x,g(x,mu(2),sig(2)),'color',[.6,.6,.6],'linewidth',lw)
    if strcmp(estimator,'mean')==1
        mup = mu(1) + sep/4;
        plot([mup, mup],[0,.5],'--','color',cI,'linewidth',lw)
        plot(x,g(x,mup,sig(1)),'color',cI,'linewidth',lw)
    else
        plot(x,g(x,mu(1),sig(1)+sep/2),'color',cI,'linewidth',lw)
    end
    set(gca,'fontsize',16)
    
    %----------------------plot xhat (upward switch)----------------------%
    m = 3;
    inds = ((m-1)*200+51):(m*200+50);
    subplot(3,3,4);hold on;
    plot(resI.full.x(inds),'color',[.6,.6,.6],'linewidth',lw);
    plot(resI.full.xhat(inds),'color',cI,'linewidth',lw);
    plot(resR.full.xhat(inds),'color',cR,'linewidth',lw);
    plot(resM.full.xhat(inds),'color',cM,'linewidth',lw);
    xlim(xlims);
    ylim([-6,6])
    set(gca,'fontsize',16)
    
    %---------------------plot xhat (downward switch)---------------------%
    subplot(3,3,7);hold on;
    plot(resI.full.x(inds),'color',[.6,.6,.6],'linewidth',lw);
    plot(resI.full.xhat(inds),'color',cI,'linewidth',lw);
    plot(resR.full.xhat(inds),'color',cR,'linewidth',lw);
    plot(resM.full.xhat(inds),'color',cM,'linewidth',lw);
    xlim(xlims+100);
    ylim([-6,6])
    set(gca,'fontsize',16)
    

    %--------------------plot thetahat: full cycle------------------------%
    subplot(3,3,2);hold on;
    plot(resI.avg.mid.theta,'color',[.6,.6,.6],'linewidth',lw);
    plot(resI.avg.mid.thetahatY,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.thetahatY,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.thetahatY,'color',cM,'linewidth',lw);
    set(gca,'fontsize',16)
    
    %-------------------plot thetahat: upward switch----------------------%
    subplot(3,3,5);hold on;
    plot(resI.avg.mid.theta,'color',[.6,.6,.6],'linewidth',lw);
    plot(resI.avg.mid.thetahatY,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.thetahatY,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.thetahatY,'color',cM,'linewidth',lw);
    xlim(xlims);
    set(gca,'fontsize',16)
    
    %-------------------plot thetahat: downward switch--------------------%
    subplot(3,3,8);hold on;
    plot(resI.avg.mid.theta,'color',[.6,.6,.6],'linewidth',lw);
    plot(resI.avg.mid.thetahatY,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.thetahatY,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.thetahatY,'color',cM,'linewidth',lw);
    xlim(xlims+100);
    set(gca,'fontsize',16)

    

    %--------------------plot firing rate: full cycle---------------------%
    subplot(3,3,6);hold on;
    plot(resI.avg.mid.y,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.y,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.y,'color',cM,'linewidth',lw);
    ylim([0,7])
    set(gca,'fontsize',16)
    
    %-------------------plot firing rate: upward switch-------------------%
    subplot(3,3,6);hold on;
    plot(resI.avg.mid.y,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.y,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.y,'color',cM,'linewidth',lw);
    xlim(xlims)
    ylim([0,7])
    set(gca,'fontsize',16)
    
    %------------------plot firing rate: downward switch------------------%
    subplot(3,3,9);hold on;
    plot(resI.avg.mid.y,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.y,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.y,'color',cM,'linewidth',lw);
    xlim(xlims+100)
    ylim([0,7])
    set(gca,'fontsize',16)


    set(gcf,'Position',[200 200 1600 1200])
    if saveFig
        saveas(gcf,[saveFolder,'firingRates.eps'],'epsc')
    end
end
%-------------------------------------------------------------------------%



%------------------------- FIGURE 2 --------------------------------------%
if ismember(2,figNums)  
    figure;hold on;
    if strcmp(estimator,'mean')==1
        xlims = [48,60];
        ylimInfo = [0,3];
        ylimRec  = [0,5];
    else
        xlims = [45,75];
        ylimInfo = [1,3];
        ylimRec  = [0,3];
    end 

    %------------------plot MI rates: full cycle--------------------------%
    subplot(2,3,1);hold on;
    plot(resI.avg.mid.MI,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.MI,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.MI,'color',cM,'linewidth',lw);
    ylim(ylimInfo)
    set(gca,'fontsize',16)
    
    %------------------plot MI rates: upward switch-----------------------%
    subplot(2,3,2);hold on;
    plot(resI.avg.mid.MI,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.MI,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.MI,'color',cM,'linewidth',lw);
    ylim(ylimInfo)
    xlim(xlims)
    set(gca,'fontsize',16)
    
    %------------------plot MI rates: downward switch---------------------%
    subplot(2,3,3);hold on;
    plot(resI.avg.mid.MI,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.MI,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.MI,'color',cM,'linewidth',lw);
    ylim(ylimInfo)
    xlim(xlims+100)
    set(gca,'fontsize',16)
    
    
    %------------------plot rec error: full cycle-------------------------%
    subplot(2,3,4);hold on;
    plot(resI.avg.mid.errRec,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.errRec,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.errRec,'color',cM,'linewidth',lw);
    ylim(ylimRec)
    set(gca,'fontsize',16)

    %------------------plot rec error: upward switch----------------------%
    subplot(2,3,5);hold on;
    plot(resI.avg.mid.errRec,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.errRec,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.errRec,'color',cM,'linewidth',lw);
    ylim(ylimRec)
    xlim(xlims)
    set(gca,'fontsize',16)
    
    %------------------plot rec error: downward switch--------------------%
    subplot(2,3,6);hold on;
    plot(resI.avg.mid.errRec,'color',cI,'linewidth',lw);
    plot(resR.avg.mid.errRec,'color',cR,'linewidth',lw);
    plot(resM.avg.mid.errRec,'color',cM,'linewidth',lw);
    ylim(ylimRec)
    xlim(xlims+100)
    set(gca,'fontsize',16)
    

    set(gcf,'Position',[200 200 1600 800])
    if saveFig
        saveas(gcf,[saveFolder,'MIrates.eps'],'epsc')
    end
end
%-------------------------------------------------------------------------%




%------------------------- FIGURE 3 --------------------------------------%
if ismember(3,figNums)  
    
    
    r = load([folder,'entropyCodes_',estimator,'_eps',num2str(0),'.mat']);
    rs.adapt.r = r.rate;
    
    r = load([folder,'entropyCodes_',estimator,'_eps',num2str(80),'.mat']);
    rs.mixed.r = r.rate;
    
    r = load([folder,'entropyCodes_',estimator,'_eps',num2str(100),'.mat']);
    rs.task.r = r.rate;

    rs.task.r( :,1) = rs.task.r( :,end);
    rs.mixed.r(:,1) = rs.mixed.r(:,end);
    rs.adapt.r(:,1) = rs.adapt.r(:,end);
    
    rTavg = [mean(resR.full.y),mean(rs.task.r( :))];
    rMavg = [mean(resM.full.y),mean(rs.mixed.r(:))];
    rIavg = [mean(resI.full.y),mean(rs.adapt.r(:))];
    
    figure;hold on;

    subplot(1,3,1);hold on;
    plot(mean(rs.task.r ),'color',cR,'linewidth',lw)
    plot(mean(rs.mixed.r),'color',cM,'linewidth',lw)
    plot(mean(rs.adapt.r),'color',cI,'linewidth',lw)
    ylim([0,7])
    set(gca,'fontsize',16)

    subplot(1,3,2);hold on;
    bar([rTavg;rMavg;rIavg]);
    ylim([0,4])
    set(gca,'fontsize',16)

    set(gcf,'Position',[200 200 1600 400])
    
    if saveFig
        saveas(gcf,[saveFolder,'entropyCodes.eps'],'epsc')
    end
end
%-------------------------------------------------------------------------%


end
