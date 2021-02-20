function [M1,M2]=plotNonlinearityLines(estimator,nLevels,eps,h)

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
paramObj0 = res.paramObj.eps0;
paramObj100 = res.paramObj.eps100;
alpha = eps./100;

phi = 1;
params = [.01,nLevels];

ybins = linspace(0,1,nLevels+1);
M1 = zeros(99,nLevels-1);
M2 = zeros(99,nLevels-1);
x = -20:.01:20;
for i=1:99
    [~,y1] = encodeWithDiscrete(x, phi, paramObj.k(i), paramObj.x0(i), params);
    for j=2:nLevels
        [~,jj] = min((y1-ybins(j)).^2);
        M1(i,j-1) = x(jj);
    end
    
    [~,y2] = encodeWithDiscrete(x, phi, (1-alpha).*paramObj0.k(i)+alpha.*paramObj100.k(i), (1-alpha).*paramObj0.x0(i)+alpha.*paramObj100.x0(i), params);
    for j=2:nLevels
        [~,jj] = min((y2-ybins(j)).^2);
        M2(i,j-1) = x(jj);
    end
end

figure;
subplot(1,2,1);plot(M1,.01:.01:.99,'-k','linewidth',2);set(gca,'fontsize',16)
subplot(1,2,2);plot(M2,.01:.01:.99,'-k','linewidth',2);set(gca,'fontsize',16)