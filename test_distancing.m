N=2; T= 100000; dt = 0.001;

Peaks = zeros(N,11,2);
T_Peak = zeros(N,11,2);
Final_Size = zeros(N,11,2);
ii = 1; jj=1;

for resp_1 = [0:0.5:5]
    for lambda_ijs = [1/10000]
        lambda_ij = (ones(N,N)-diag(ones(N,1)))*lambda_ijs; % baseline immigration constant
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
        small_positive = 0.0001;
        
        omega = ones(N,1)* 1; %[0.5,1]
        omega(1) = 1;
        omega(2) = 1;
        kappa = ones(N,1)*0; % vary between 0.5-5 % greater than one values give you overreaction
        kappa(1) = resp_1;
        kappa(2) = resp_1;
        alpha_ij = ones(N,1)*1; % vary between 1-5 % NOT USED ANYMORE
        
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
        Peaks(1,ii,jj) = max(Inf(1,:));
        Peaks(2,ii,jj) = max(Inf(2,:));
        T_Peak(1,ii,jj) = find(max(Inf(1,:))==Inf(1,:))*dt;
        T_Peak(2,ii,jj) = find(max(Inf(2,:))==Inf(2,:))*dt;
        Final_Size(1,ii,jj) = 1-Susc(1,end);
        Final_Size(2,ii,jj) = 1-Susc(2,end)
        jj = jj+1;
    end
    ii = ii+1;
    jj=1;
end

%%

time_horizon = (1:T)*dt;
figure
hold on
plot(Susc(1,:),'LineWidth',2)
% plot(Exp(1,:),'LineWidth',2)
plot(Inf(1,:),'LineWidth',2)
plot(Susc(2,:),'LineWidth',2)
% plot(Exp(2,:),'LineWidth',2)
plot(Inf(2,:),'LineWidth',2)
lll = legend('S-Locality 1','I-Locality 1','S-Locality 2','I-Locality 2');
set(lll,'FontSize',16)
% title('S-E-I City 1 vs 2')
yyy = ylabel('Ratio','FontSize',16);
xxx = xlabel('Time (days)','FontSize',16);
ax = gca;
% set(gca,'XTick',time_horizon);
set(ax,'FontSize',16)
ax.XTickLabel = {'0','10','20','30','40','50','60','70','80','90','100'};

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
percentage_reduction = zeros(11,1);

max_outbreak = Final_Size(2,1,1);
percentage_reduction = 100*(max_outbreak-Final_Size(2,:,1))/max_outbreak;

percentage_reduction_peak = zeros(11,1);
max_peak = Peaks(2,1,1);
percentage_reduction_peak = 100*(max_peak-Peaks(2,:,1))/max_peak;
%%
figure
hold on
plot(0:0.5:5, percentage_reduction,'k-','LineWidth',2,'MarkerSize',8)
plot(0:0.5:5, percentage_reduction_peak,'b-','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Percentage reduction','FontSize',16);
xxx = xlabel('Social distancing ($\alpha_i$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Outbreak size','Peak infected');
legend boxoff
% set(lll,'Interpreter','latex')
set(lll,'FontSize',16);

 axis([0 5 0 100]);
% set(gca,'XTick',0:0.5:5);
%  ax.XTickLabel = {'0','1','2','3','4','5'};

%  ax.XTickLabel = {'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'};

%%
figure
hold on
plot(0:.5:5,(T_Peak(2,:,1)),'k-o','LineWidth',2,'MarkerSize',8)
plot(0:.5:5,(T_Peak(2,:,2)),'b-o','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Time of peak','FontSize',16);
xxx = xlabel('Social distancing ($\alpha_i$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Low mobility','High mobility');
legend boxoff
% set(lll,'Interpreter','latex')
set(lll,'FontSize',16);

% axis([0 1 0 0.2]);
% set(gca,'XTick',0:0.2:1);
% ax.XTickLabel = {'1','0.8','0.6','0.4','0.2','0'};

%%
figure
hold on
plot(0:.5:5,(Peaks(2,:,1)),'k-o','LineWidth',2,'MarkerSize',8)
plot(0:.5:5,(Peaks(2,:,2)),'b-o','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Ratio infected at peak','FontSize',16);
xxx = xlabel('Social distancing ($\alpha_i$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Low mobility','High mobility');
legend boxoff
% set(lll,'Interpreter','latex')
set(lll,'FontSize',16);