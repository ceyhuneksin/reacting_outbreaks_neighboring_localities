N=2; T= 100000; dt = 0.001;
P_level = 1*ones(N,T);
% generate network
A_contact = ones(N,N)-diag(ones(N,1));

%% migration 
% need to tie migration to network structure when you include the network
lambda_in = zeros(N,N,T);
lambda_out = zeros(N,N,T);
lambda_ij = (ones(N,N)-diag(ones(N,1)))/10000; 
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
omega(2) = 1; 
kappa = ones(N,1)*0; % vary between 0.5-5 % greater than one values give you overreaction
kappa(1) = 0; % distancing coefficient
kappa(2) = 0; % distancing coefficient
alpha_ij = ones(N,1); % vary between 1-5 % NOT USED ANYMORE

c_ij = 0;
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
new_cases = zeros(2,1);
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
        
        x = [Susc(j,t), Exp(j,t), Inf(j,t), beta, delta, mu, alpha_ij(j),kappa(j),omega(j),new_cases(j)];
        x_integrated = SEIR_j(N_j,x,sum(migration_in), migration_out,j,awareness_at_neighbors,dt);
%        x_integrated = struct('S',Susc(2),'E', Exp(2), 'I',Infec(2),'beta_dist',beta);
        Susc(j,t+1) = x_integrated.S;
        Exp(j,t+1) = x_integrated.E;
        Inf(j,t+1) = x_integrated.I;
        beta_dist(j,t+1) = x_integrated.beta_dist;
        new_cases(j) = max((Inf(j,t+1)-Inf(j,t))/dt,0);
        P_level(j,t+1) = N_j + N_j*(sum(lambda_in(:,j,t))-sum(lambda_in(j,:,t)))*dt;
    end
end
%%
Peak_1 = max(Inf(1,:))
Peak_2 = max(Inf(2,:))
E_Peak_1 = max(Exp(1,:))
E_Peak_2 = max(Exp(2,:))

T_Peak_1 = find(max(Inf(1,:))==Inf(1,:))*dt
T_Peak_2 = find(max(Inf(2,:))==Inf(2,:))*dt
Size_1 = 1-Susc(1,end)
Size_2 = 1-Susc(2,end)
%%
time_horizon = (1:T)*dt;
figure
hold on
plot(Susc(1,:),'b','LineWidth',2)
% plot(Exp(1,:),'LineWidth',2)
plot(Inf(1,:),'r','LineWidth',2)
plot(Susc(2,:),'b-.','LineWidth',2)
% plot(Exp(2,:),'LineWidth',2)
plot(Inf(2,:),'r-.','LineWidth',2)
lll = legend('$S_1(t)$','$I_1(t)$','$S_2(t)$','$I_2(t)$');
set(lll,'Interpreter','latex')
set(lll,'FontSize',16)
% title('S-E-I City 1 vs 2')
yyy = ylabel('Ratio of population','FontSize',16);
xxx = xlabel('Time (days)','FontSize',16);
ax = gca;
% set(gca,'XTick',time_horizon);
set(ax,'FontSize',16)
ax.XTickLabel = {'0','20','40','60','80','100'};
legend boxoff
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
