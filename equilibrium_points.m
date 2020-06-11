%% SIR affected distancing
%% Find the zero crossing points of the function on the left of eq. 10
% no distancing
ll = 0:0.01:1; ii = 1;
f = zeros(length(ll),4);
beta = 2.5; delta = 1;
for x = ll
    f(ii,1) = x - exp(-beta*(1-x));
    ii = ii +1;
end

%% lower bound on R(infty) --- upper bound bound on S
for k = 1:1:3
    ii = 1;
for x = ll
    f(ii,3+k+1) = x - power((1+k*beta*(1-x)/delta),-1/k);
    ii = ii +1;
end
end
%% lower bound on R(infty) --- upper bound on S
for k = 1:1:3
    ii = 1;
for x = ll
    f(ii,k+1) = x - exp(-beta*x^k*(1-x)/delta);
    ii = ii +1;
end
end

%% plot Figure 1
figure 

hold on
lll = xlabel('$S(\infty)$', 'FontSize',18);
set(lll,'Interpreter','Latex');
lll=ylabel('$S-(1+\alpha_1 \mathcal{R}_1  (1-S))^{-1/\alpha_1}$','FontSize',18);

scatter(0.11,0,100,'kd','Filled', 'LineWidth',2)
scatter(0.36,0,100,'rd','Filled', 'LineWidth',2)
scatter(0.61,0,100,'bd','Filled', 'LineWidth',2)

% scatter(0.11,0,'ko','Filled', 'LineWidth',2)
scatter(0.4,0,100,'ro','Filled', 'LineWidth',2)
scatter(0.65,0,100,'bo','Filled', 'LineWidth',2)

% lll=ylabel('$S-\exp(-\mathcal{R}_1 S^{\alpha_1} (1-S))$','FontSize',18);
set(lll,'Interpreter','Latex');
% lll = legend('$k=0$','$k=1$','$k=2$', '$k=3$')
lll = legend('${S}(\infty): \alpha_1=0$','${S}(\infty): \alpha_1=1$', '${S}(\infty): \alpha_1=3$','$\hat{S}(\infty): \alpha_1=1$', '$\hat{S}(\infty): \alpha_1=3$')
set(lll,'Interpreter','Latex');

legend boxoff

set(gca,'DefaultLegendAutoUpdate','off')

plot(ll,f(:,1),'k','LineWidth', 1.5)
plot(ll,f(:,5),'r-','LineWidth', 1.5)
% plot(ll,f(:,6),'b--','LineWidth', 1.5)
plot(ll,f(:,7),'b-','LineWidth', 1.5)

% plot(ll,f(:,2),'r--','LineWidth', 1.5)
% plot(ll,f(:,3),'b--','LineWidth', 1.5)
% plot(ll,f(:,4),'b--','LineWidth', 1.5)

plot(ll,zeros(length(ll),1),'k--')

% scatter(delta/beta,0,'ro','Filled', 'LineWidth',2)


% scatter(delta/beta,0,'ro','Filled', 'LineWidth',2)
% legend('Recovered')
% set(gca,'XTick',[0 200]);
% set(gca,'YTick',[0:0.08:0.4]);
% ax = gca;
% ax.XTickLabel = {'0','\infty','\infty'};
% ax.YTickLabel = {'0','8%','16%','24%','32%', '40%'};
% ylim([0 0.4])

set(gca,'FontSize',18)



%%
%% upper bound on R(infty) importing locality
beta = 2.5; ll = 0:0.01:1; ii = 1;
g = zeros(length(ll),6);
omega = 1/2;
r_1 = [1-0.11 1-0.4 1-0.65];
for mm =1:3 % origin outbreak size counter
    zz = 1;
    for k = [1 3] % strength of response at the locality 2
        ii = 1;
        beta_hat = beta * (1-(1-omega)*r_1(mm))^k;
        for x = ll
            g(ii,(mm-1)*2+zz) = x - exp(-beta_hat*(1-x));
            ii = ii +1;
        end
        zz= zz+1;
    end
end
%%
h = zeros(length(ll),6);
omega = 1/4;
r_1 = [1-0.11 1-0.4 1-0.65]; % origin outbreak size
for mm =1:3 % origin outbreak size counter
    zz = 1;
    for k = [1 3] % strength of response at the locality 2
        ii = 1;
        beta_hat = beta * (1-(1-omega)*r_1(mm))^k;
        for x = ll
            h(ii,(mm-1)*2+zz) = x - exp(-beta_hat*(1-x));
            ii = ii +1;
        end
        zz= zz+1;
    end
end

%%
set(groot,'defaultLegendAutoUpdate','off')
figure 
% scatter(0.7,0,50,'kd', 'LineWidth',2)
% scatter(0.997,0,50,'ks', 'LineWidth',2)
scatter(0.5,0,100,'rd', 'Filled','LineWidth',2)
hold on
scatter(1,0,100,'bd', 'Filled','LineWidth',2)
% scatter(0.36,0,100,'bo', 'LineWidth',2)
% scatter(0.83,0,100,'bx','LineWidth',2)

% scatter(0.51,0,50,'ko','Filled', 'LineWidth',2)
% scatter(1,0,50,'k*','Filled', 'LineWidth',2)
scatter(0.29,0,100,'ro','Filled', 'LineWidth',2)
scatter(1,0,100,'bo','Filled', 'LineWidth',2)
% scatter(0.19,0,100,'bo','Filled', 'LineWidth',2)
% scatter(0.49,0,100,'bh','Filled', 'LineWidth',2)




lll = xlabel('$S(\infty)$', 'FontSize',18);
set(lll,'Interpreter','Latex');
lll=ylabel('$S-\exp(-\tilde\mathcal{R}_2 (1-S))$','FontSize',18);
set(lll,'Interpreter','Latex');


lll = legend('${S}_2(\infty): \alpha_2=1$','${S}_2(\infty): \alpha_2=3$', '$\tilde{S}_2(\infty): \alpha_2=1$','$\tilde{S}_2(\infty): \alpha_2=3$')
set(lll,'Interpreter','Latex');
legend boxoff

set(gca,'FontSize',18)
% plot(ll,g(:,1),'k-','LineWidth', 1.5)
% plot(ll,g(:,2),'k--','LineWidth', 1.5)
plot(ll,g(:,3),'r-','LineWidth', 1.5)
plot(ll,g(:,4),'b-','LineWidth', 1.5)
% plot(ll,g(:,5),'b-','LineWidth', 1.5)
% plot(ll,g(:,6),'b--','LineWidth', 1.5)
plot(ll,zeros(length(ll),1),'k--')


%%
figure 
% scatter(0.7,0,50,'kd', 'LineWidth',2)
% scatter(0.997,0,50,'ks', 'LineWidth',2)
% scatter(0.5,0,100,'ro', 'LineWidth',2)
hold on
% scatter(1,0,100,'rx', 'LineWidth',2)
scatter(0.36,0,100,'rd','Filled', 'LineWidth',2)
scatter(0.83,0,100,'bd','Filled','LineWidth',2)

% scatter(0.51,0,50,'ko','Filled', 'LineWidth',2)
% scatter(1,0,50,'k*','Filled', 'LineWidth',2)
% scatter(0.29,0,100,'ro','Filled', 'LineWidth',2)
% scatter(1,0,100,'rh','Filled', 'LineWidth',2)
scatter(0.19,0,100,'ro','Filled', 'LineWidth',2)
scatter(0.49,0,100,'bo','Filled', 'LineWidth',2)

lll = xlabel('$S(\infty)$', 'FontSize',18);
set(lll,'Interpreter','Latex');
lll=ylabel('$S-\exp(-\tilde\mathcal{R}_2 (1-S))$','FontSize',18);
set(lll,'Interpreter','Latex');

lll = legend('${S}_2(\infty): \alpha_2=1$','${S}_2(\infty): \alpha_2=3$', '$\tilde{S}_2(\infty): \alpha_2=1$','$\tilde{S}_2(\infty): \alpha_2=3$')
% lll = legend('${S}_2(\infty): \alpha_2=1, \alpha_1=1$','${S}_2(\infty): \alpha_2=3, \alpha_1=1$', '${S}_2(\infty): \alpha_2=1, \alpha_1=3$', '${S}_2(\infty): \alpha_2=3, \alpha_1=3$','$\tilde{S}_2(\infty): \alpha_2=1, \alpha_1=1$','$\tilde{S}_2(\infty): \alpha_2=3, \alpha_1=1$', '$\tilde{S}_2(\infty): \alpha_2=1, \alpha_1=3$', '$\tilde{S}_2(\infty): \alpha_2=3, \alpha_1=3$')
set(lll,'Interpreter','Latex');
legend boxoff

set(gca,'FontSize',18)
% plot(ll,g(:,1),'k-','LineWidth', 1.5)
% plot(ll,g(:,2),'k--','LineWidth', 1.5)
% plot(ll,g(:,3),'r-','LineWidth', 1.5)
% plot(ll,g(:,4),'r--','LineWidth', 1.5)
plot(ll,g(:,5),'r-','LineWidth', 1.5)
plot(ll,g(:,6),'b-','LineWidth', 1.5)
plot(ll,zeros(length(ll),1),'k--')

%%
figure 
% scatter(0.7,0,50,'kd', 'LineWidth',2)
% scatter(0.997,0,50,'ks', 'LineWidth',2)
scatter(0.7,0,100,'rd', 'Filled','LineWidth',2)
hold on
scatter(1,0,100,'bd', 'Filled','LineWidth',2)
% scatter(0.36,0,100,'bo', 'LineWidth',2)
% scatter(0.83,0,100,'bx','LineWidth',2)

% scatter(0.51,0,50,'ko','Filled', 'LineWidth',2)
% scatter(1,0,50,'k*','Filled', 'LineWidth',2)
scatter(0.51,0,100,'ro','Filled', 'LineWidth',2)
scatter(1,0,100,'bo','Filled', 'LineWidth',2)
% scatter(0.19,0,100,'bo','Filled', 'LineWidth',2)
% scatter(0.49,0,100,'bh','Filled', 'LineWidth',2)

lll = xlabel('$S(\infty)$', 'FontSize',18);
set(lll,'Interpreter','Latex');
lll=ylabel('$S-\exp(-\tilde\mathcal{R}_2 (1-S))$','FontSize',18);
set(lll,'Interpreter','Latex');


lll = legend('${S}_2(\infty): \alpha_2=1$','${S}_2(\infty): \alpha_2=3$', '$\tilde{S}_2(\infty): \alpha_2=1$','$\tilde{S}_2(\infty): \alpha_2=3$')
set(lll,'Interpreter','Latex');
legend boxoff

set(gca,'FontSize',18)
% plot(ll,g(:,1),'k-','LineWidth', 1.5)
% plot(ll,g(:,2),'k--','LineWidth', 1.5)
plot(ll,h(:,3),'r-','LineWidth', 1.5)
plot(ll,h(:,4),'b-','LineWidth', 1.5)
% plot(ll,g(:,5),'b-','LineWidth', 1.5)
% plot(ll,g(:,6),'b--','LineWidth', 1.5)
plot(ll,zeros(length(ll),1),'k--')


%%
figure 
% scatter(0.7,0,50,'kd', 'LineWidth',2)
% scatter(0.997,0,50,'ks', 'LineWidth',2)
% scatter(0.5,0,100,'ro', 'LineWidth',2)
hold on
% scatter(1,0,100,'rx', 'LineWidth',2)
scatter(0.38,0,100,'rd','Filled', 'LineWidth',2)
scatter(0.998,0,100,'bd','Filled','LineWidth',2)

% scatter(0.51,0,50,'ko','Filled', 'LineWidth',2)
% scatter(1,0,50,'k*','Filled', 'LineWidth',2)
% scatter(0.29,0,100,'ro','Filled', 'LineWidth',2)
% scatter(1,0,100,'rh','Filled', 'LineWidth',2)
scatter(0.25,0,100,'ro','Filled', 'LineWidth',2)
scatter(1,0,100,'bo','Filled', 'LineWidth',2)

lll = xlabel('$S(\infty)$', 'FontSize',18);
set(lll,'Interpreter','Latex');
lll=ylabel('$S-\exp(-\tilde\mathcal{R}_2 (1-S))$','FontSize',18);
set(lll,'Interpreter','Latex');

lll = legend('${S}_2(\infty): \alpha_2=1$','${S}_2(\infty): \alpha_2=3$', '$\tilde{S}_2(\infty): \alpha_2=1$','$\tilde{S}_2(\infty): \alpha_2=3$')
% lll = legend('${S}_2(\infty): \alpha_2=1, \alpha_1=1$','${S}_2(\infty): \alpha_2=3, \alpha_1=1$', '${S}_2(\infty): \alpha_2=1, \alpha_1=3$', '${S}_2(\infty): \alpha_2=3, \alpha_1=3$','$\tilde{S}_2(\infty): \alpha_2=1, \alpha_1=1$','$\tilde{S}_2(\infty): \alpha_2=3, \alpha_1=1$', '$\tilde{S}_2(\infty): \alpha_2=1, \alpha_1=3$', '$\tilde{S}_2(\infty): \alpha_2=3, \alpha_1=3$')
set(lll,'Interpreter','Latex');
legend boxoff

set(gca,'FontSize',18)
% plot(ll,g(:,1),'k-','LineWidth', 1.5)
% plot(ll,g(:,2),'k--','LineWidth', 1.5)
% plot(ll,g(:,3),'r-','LineWidth', 1.5)
% plot(ll,g(:,4),'r--','LineWidth', 1.5)
plot(ll,h(:,5),'r-','LineWidth', 1.5)
plot(ll,h(:,6),'b-','LineWidth', 1.5)
plot(ll,zeros(length(ll),1),'k--')
% plot(ll,g(:,1),'k-','LineWidth', 1.5)
% plot(ll,g(:,2),'k--','LineWidth', 1.5)
