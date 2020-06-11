N=2; T= 200000; dt = 0.001;

Peaks = zeros(N,11,2);
T_Peak = zeros(N,11,2);
Final_Size = zeros(N,11,2);
ii = 1; jj=1;

for resp_1 = [1 3] % 1
    lambda_ijs = 1/100000;
    %     for lambda_ijs = [1/1000000 1/10000] %1/1000000
    lambda_ij = (ones(N,N)-diag(ones(N,1)))/100000; %lambda_ijs % baseline immigration constant
    for omega_2 =  flip(0:0.1:1) %0
        P_level = 1*ones(N,T);
        % generate network
        A_contact = ones(N,N)-diag(ones(N,1));
        
        %% migration
        % need to tie migration to network structure when you include the network
        lambda_in = zeros(N,N,T);
        lambda_out = zeros(N,N,T);
        % 0.005 gives ~7-11 day separation between peaks
        % 0.05 gives ~1-2 day separation between peaks
        
        %% disease state
        Susc = zeros(N,T); % fraction of population susceptible
        Exp = zeros(N,T); % fraction of population exposed
        Inf = zeros(N,T); % fraction of population infected
        beta_dist = zeros(N,T);
        
        %% Set model constants
        
        beta = 0.5; % baseline infectivity rate (S to E)
        mu = 0.5; %
        delta = 0.33; % healing rate (I to R)
        
        omega = ones(N,1)* 1; %[0.5,1]
        omega(1) = 1;
        omega(2) = omega_2;
        kappa = ones(N,1)*0; % vary between 0.5-5 % greater than one values give you overreaction
        kappa(1) = resp_1; % strength of response at Locality 1
        kappa(2) = 1; % strength of response at Locality 2
        alpha_ij = ones(N,1)*1; % vary between 1-5 % NOT USED ANYMORE
        
        k_ij = 0;
        %% Store data
        % state_init = struct('S',S_init,'I',I_init, 'E', E_init, 'mu',mu_init);
        % state_post = struct('S',S_init,'I',I_init, 'E', E_init, 'mu',mu_init);
        
        %% Initialize disease
        for i = 1:N
            Susc(i,1) = 1;
        end
        j = 1; %randi(N) ;
        Susc(j,1) = (P_level(j,1)-0.001)/P_level(j,1);
        Exp(j,1) = 0.001;
        %%
        for t = 1 : T-1
            % Compute in and out migration
            for j = 1 : N
                N_j = P_level(j,t);
                contact_neighs = find(A_contact(j,:) == 1);
                for i = 1:N
                    awareness_at_i = 1-(Susc(i,t)+Exp(i,t));
                    awareness_at_j = 1-(Susc(j,t)+Exp(j,t));
                    lambda_in(i,j,t) = lambda_ij(i,j) * ((1+awareness_at_i)/(1+awareness_at_j))^k_ij;
                    lambda_in(j,i,t) = lambda_ij(j,i) * ((1+awareness_at_j)/(1+awareness_at_i))^k_ij;
                end
                migration_in = zeros(N,2);
                migration_out = zeros(2,1);
                migration_in(:,1) = lambda_in(:,j,t).*Susc(:,t)./(Susc(:,t)+Exp(:,t)); % Susc emigration
                migration_in(:,2) = lambda_in(:,j,t).*Exp(:,t)./(Susc(:,t)+Exp(:,t)); % Exposed emigration
                migration_out(1) = sum(lambda_in(j,:,t)*Susc(j,t)/(Susc(j,t)+Exp(j,t))); %  Susc migration
                migration_out(2) = sum(lambda_in(j,:,t)*Exp(j,t)/(Susc(j,t)+Exp(j,t))); %  Exposed migration
                % propagate dynamics for j
                % x = [S, E, I,beta,delta,mu]
                awareness_at_neighbors = 1-(Susc(contact_neighs,t)+Exp(contact_neighs,t)); % infected and recovered at neighbors
                
                x = [Susc(j,t), Exp(j,t), Inf(j,t), beta, delta, mu, alpha_ij(j),kappa(j),omega(j)];
                x_integrated = SEIR_j(N_j,x,sum(migration_in), migration_out,j,awareness_at_neighbors,dt);
                %        x_integrated = struct('S',Susc(2),'E', Exp(2), 'I',Infec(2),'beta_dist',beta);
                Susc(j,t+1) = x_integrated.S;
                Exp(j,t+1) = x_integrated.E;
                Inf(j,t+1) = x_integrated.I;
                beta_dist(j,t+1) = x_integrated.beta_dist;
                P_level(j,t+1) = N_j + N_j*(sum(lambda_in(:,j,t))-sum(lambda_in(j,:,t)))*dt;
            end
        end
        Susc(1,end)
        Susc(2,end)
        Peaks(1,ii,jj) = max(Inf(1,:));
        Peaks(2,ii,jj) = max(Inf(2,:));
        T_Peak(1,ii,jj) = find(max(Inf(1,:))==Inf(1,:))*dt;
        T_Peak(2,ii,jj) = find(max(Inf(2,:))==Inf(2,:))*dt;
        Final_Size(1,ii,jj) = 1-Susc(1,end);
        Final_Size(2,ii,jj) = 1-Susc(2,end)
        ii = ii+1;
    end
    ii = 1;
    jj=jj+1;
end

% end
%%
beta = 2.5; ll = 0:0.01:1; 
g = zeros(length(ll),10,4);
r_1 = [1-0.4 1-0.65]; % outbreak size at the origin with weak and strong response
for mm =[1 2] % origin outbreak size counter
    zz = 1;
    for k = [1 3] % strength of response at the locality 2
        jj= 1;
        for omega = flip(0:0.1:1)
            beta_hat = beta * (1-(1-omega)*r_1(mm))^k;
            ii = 1;
            for x = ll
                g(ii,jj,(mm-1)*2+zz) = x - exp(-beta_hat*(1-x));
                ii = ii +1;
            end
            jj = jj +1;
        end
        zz= zz+1;
    end
end

upper_outbreak_val = zeros(11,4); 
%ii = 1 weak Locality 1 - weak Locality 2; = 2 low Locality 1 - strong Locality 2
%ii = 3 strong Locality 1 - weak Locality 2; = 4 strong Locality 1 - strong Locality 2
upper_outbreak_size = zeros(11,4);
for ii = 1:4
    for jj = 1:11
        upper_outbreak_size(jj,ii) = 1-(max(find(g(:,jj,ii)<0)))/100;
        upper_outbreak_val(jj,ii) = g(max(find(g(:,jj,ii)<0)),jj,ii);
    end
end
%%
% h = zeros(length(ll),6);
% omega = 1/4;
% r_1 = [1-0.4 1-0.65]; % origin outbreak size
% for mm =[1 2] % origin outbreak size counter
%     zz = 1;
%     for k = [1 3] % strength of response at the locality 2
%         ii = 1;
%         beta_hat = beta * (1-(1-omega)*r_1(mm))^k;
%         for x = ll
%             h(ii,(mm-1)*2+zz) = x - exp(-beta_hat*(1-x));
%             ii = ii +1;
%         end
%         zz= zz+1;
%     end
% end


%%
figure
hold on
plot(Susc(1,:),'LineWidth',2)
plot(Exp(1,:),'LineWidth',2)
plot(Inf(1,:),'LineWidth',2)
plot(Susc(2,:),'LineWidth',2)
plot(Exp(2,:),'LineWidth',2)
plot(Inf(2,:),'LineWidth',2)
lll = legend('S-Locality 1','E-Locality 1','I-Locality 1','S-Locality 2','E-Locality 2','I-Locality 2');
set(lll,'FontSize',16)
% title('S-E-I City 1 vs 2')
yyy = ylabel('Ratio','FontSize',16);
xxx = xlabel('Time','FontSize',16);
ax = gca;
set(ax,'FontSize',16)
ax.XTickLabel = {'0','20','40','60','80','100'};
legend boxoff
% %%
% %%
% figure
% hold on
% plot(Susc(2,:),'LineWidth',2)
% plot(Exp(2,:),'LineWidth',2)
% plot(Inf(2,:),'LineWidth',2)
% lll = legend('S','E','I');
% set(lll,'FontSize',16)
% title('City 2')
% yyy = ylabel('Ratio','FontSize',16);
% xxx = xlabel('Time','FontSize',16);
%%
% figure
% hold on
% plot(P_level(1,:),'LineWidth',2)
% plot(P_level(2,:),'LineWidth',2)
% lll = legend('Locality 1','Locality 2');
% set(lll,'FontSize',16)
% yyy = ylabel('Population','FontSize',16);
% xxx = xlabel('Time','FontSize',16);
%%
figure
hold on
plot(beta_dist(1,2:T)/beta,'LineWidth',2)
plot(beta_dist(2,2:T)/beta,'LineWidth',2)
lll = legend('Locality 1','Locality 2');
set(lll,'FontSize',16)
yyy = ylabel('Distancing effort','FontSize',16);
xxx = xlabel('Time','FontSize',16);
ax = gca;
set(ax,'FontSize',16)

%%
figure
hold on
plot(0:.1:1,Final_Size(2,:,1),'k-o','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,Final_Size(2,:,2),'b-o','LineWidth',2,'MarkerSize',8)
% change indices 1 and 3 to 2 and 4 when response at Locality 2 is strong
plot(0:.1:1,upper_outbreak_size(:,1),'k--','LineWidth',2,'MarkerSize',8) 
plot(0:.1:1,upper_outbreak_size(:,3),'b--','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Outbreak Size at Locality 2','FontSize',16);
xxx = xlabel('Adopted awareness at Locality 2 ($\omega_{21}$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Weak response at Locality 1','Strong response at Locality 1',...
    'Upper bound weak response','Upper bound strong response');
legend boxoff
set(lll,'Interpreter','latex')
set(lll,'FontSize',16);

axis([0 1 0.0 1]);
set(gca,'XTick',0:0.2:1);
% ax.XTickLabel = {'1','0.8','0.6','0.4','0.2','0'};
%
% plot(flip(Final_Size(2,:,1,1)),'r-','LineWidth',2)
% plot(flip(Final_Size(2,:,1,2)),'b-','LineWidth',2)
% plot(flip(Final_Size(2,:,2,1)),'r-','LineWidth',2)
% plot(flip(Final_Size(2,:,2,2)),'b-','LineWidth',2)
%% lower bound on the outbreak size at the origin
% ll = 0:0.01:1; ii = 1;
% f = zeros(length(ll),4);
% beta = 2.5; delta = 1;
% for k = 1:1:3
%     ii = 1;
% for x = ll
%     f(ii,3+k+1) = x - power((1+k*beta*(1-x)/delta),-1/k);
%     ii = ii +1;
% end
% end

%%
figure
hold on
plot(0:.1:1,(Final_Size(1,:,1)+Final_Size(2,:,1))/2,'k-o','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,(Final_Size(1,:,2)+Final_Size(2,:,2))/2,'b-o','LineWidth',2,'MarkerSize',8)

r_1 = [1-0.4 1-0.65];
plot(0:.1:1,(r_1(1)+upper_outbreak_size(:,2))/2,'k--','LineWidth',2,'MarkerSize',8) 
plot(0:.1:1,(r_1(2)+upper_outbreak_size(:,4))/2,'b--','LineWidth',2,'MarkerSize',8)

% plot(0:.1:1,(Final_Size(1,:,2,1)+Final_Size(2,:,2,1))/2,'k-+','LineWidth',2,'MarkerSize',8)
% plot(0:.1:1,(Final_Size(1,:,2,2)+Final_Size(2,:,2,2))/2,'b-+','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Total Outbreak Size','FontSize',16);
xxx = xlabel('Adopted awareness at Locality 2 ($\omega_{21}$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Weak response at Locality 1 ($\alpha_1=1$)','Strong response at Locality 1 ($\alpha_1=3$)',...
    'Approximate ($\hat{R}_1(\infty) + \tilde{R}_2(\infty)$) ($\alpha_1=1$)','Approximate ($\hat{R}_1(\infty) + \tilde{R}_2(\infty)$) ($\alpha_1=3$)');
legend boxoff
set(lll,'Interpreter','latex')
set(lll,'FontSize',16);

axis([0 1 0.0 0.8]);
set(gca,'XTick',0:0.2:1);