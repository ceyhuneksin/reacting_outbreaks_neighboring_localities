%% generates Figures 2 and 3
%% Current values are for Figures 2 and 3 (right)

N=2; T= 100000; dt = 0.001;

Peaks = zeros(N,11,2,2);
T_Peak = zeros(N,11,2,2);
Final_Size = zeros(N,11,2,2);
ii = 1; jj=1; kk = 1;

for resp_1 = [1 3] % 1
    for lambda_ijs = [1/1000000 1/10000] %1/1000000
        lambda_ij = (ones(N,N)-diag(ones(N,1)))*lambda_ijs; % baseline immigration constant
        for omega_2 =  flip(0:0.1:1) %0
            P_level = 1*ones(N,T);
            % generate network
            A_contact = ones(N,N)-diag(ones(N,1));
            
            %% migration
            % need to tie migration to network structure when you include the network
            lambda_in = zeros(N,N,T);
            lambda_out = zeros(N,N,T);
            
            %% disease state
            Susc = zeros(N,T); % fraction of population susceptible
            Exp = zeros(N,T); % fraction of population exposed
            Inf = zeros(N,T); % fraction of population infected
            beta_dist = zeros(N,T);
            
            %% Set model constants
            
            beta = 0.5; % baseline infectivity rate (S to E)
            mu = 1/2; % incubation period, can vary between 1/2-1/4
            delta = 1/3; % healing rate (I to R) can vary between 2-4
            
            omega = ones(N,1)* 1; %[0.5,1]
            omega(1) = 1;% origin's weight on local outbreak size $\omega_11$
            omega(2) = omega_2;% importing locality's weight on local outbreak size $\omega_22$
            kappa = ones(N,1)*0; % vary between 0.5-5 % greater than one values give you overreaction
            kappa(1) = resp_1; % value of social distancing at origin
            kappa(2) = 3; % linear (1) or strong (3) response at locality 2
            
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
                    lambda_in(i,j,t) = lambda_ij(i,j); % in case mobility is time-varying
                    lambda_in(j,i,t) = lambda_ij(j,i); % in case mobility is time-varying
                    migration_in = zeros(N,2);
                    migration_out = zeros(2,1);
                    migration_in(:,1) = lambda_in(:,j,t).*Susc(:,t)./(Susc(:,t)+Exp(:,t)); % Susc emigration
                    migration_in(:,2) = lambda_in(:,j,t).*Exp(:,t)./(Susc(:,t)+Exp(:,t)); % Exposed emigration
                    migration_out(1) = sum(lambda_in(j,:,t)*Susc(j,t)/(Susc(j,t)+Exp(j,t))); %  Susc migration
                    migration_out(2) = sum(lambda_in(j,:,t)*Exp(j,t)/(Susc(j,t)+Exp(j,t))); %  Exposed migration
                    % propagate dynamics for j
                    % x = [S, E, I,beta,delta,mu]
                    awareness_at_neighbors = 1-(Susc(contact_neighs,t)+Exp(contact_neighs,t)); % infected and recovered at neighbors
                    
                    x = [Susc(j,t), Exp(j,t), Inf(j,t), beta, delta, mu, kappa(j),omega(j)];
                    x_integrated = SEIR_j(N_j,x,sum(migration_in), migration_out,j,awareness_at_neighbors,dt);
                    %        x_integrated = struct('S',Susc(2),'E', Exp(2), 'I',Infec(2),'beta_dist',beta);
                    Susc(j,t+1) = x_integrated.S;
                    Exp(j,t+1) = x_integrated.E;
                    Inf(j,t+1) = x_integrated.I;
                    beta_dist(j,t+1) = x_integrated.beta_dist;
                    P_level(j,t+1) = N_j + N_j*(sum(lambda_in(:,j,t))-sum(lambda_in(j,:,t)))*dt;
                end
            end
            Peaks(1,ii,jj,kk) = max(Inf(1,:));
            Peaks(2,ii,jj,kk) = max(Inf(2,:));
            T_Peak(1,ii,jj,kk) = find(max(Inf(1,:))==Inf(1,:))*dt;
            T_Peak(2,ii,jj,kk) = find(max(Inf(2,:))==Inf(2,:))*dt;
            Final_Size(1,ii,jj,kk) = 1-Susc(1,end);
            Final_Size(2,ii,jj,kk) = 1-Susc(2,end)
            ii = ii+1;
        end
        ii = 1;
        jj = jj+1;
    end
    jj=1;
    kk=kk+1;
end
%%

%% Susceptible and infected trajectories at both localities
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
%% shows the effective reduction as ratio of contacts due to distancing for the last run
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

%% Figure 2 (left or right): Outbreak size at Locality 2 as a function of adopted awareness
figure
hold on
plot(0:.1:1,Final_Size(2,:,1,1),'k-o','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,Final_Size(2,:,1,2),'b-o','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,Final_Size(2,:,2,1),'k-+','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,Final_Size(2,:,2,2),'b-+','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Outbreak Size at Locality 2','FontSize',16);
xxx = xlabel('Adopted awareness at Locality 2 ($1-\omega_{22}$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Low mobility, weak response at Locality 1','Low mobility, strong response at Locality 1',...
    'High mobility, weak response at Locality 1','High mobility, strong response at Locality 1');
legend boxoff
set(lll,'Interpreter','latex')
set(lll,'FontSize',16);

axis([0 1 0.0 0.8]);
set(gca,'XTick',0:0.2:1);
% ax.XTickLabel = {'1','0.8','0.6','0.4','0.2','0'};
%
% plot(flip(Final_Size(2,:,1,1)),'r-','LineWidth',2)
% plot(flip(Final_Size(2,:,1,2)),'b-','LineWidth',2)
% plot(flip(Final_Size(2,:,2,1)),'r-','LineWidth',2)
% plot(flip(Final_Size(2,:,2,2)),'b-','LineWidth',2)

%% Figure 3 (left or right): Total Outbreak size as a function of adopted awareness
figure
hold on
plot(0:.1:1,(Final_Size(1,:,1,1)+Final_Size(2,:,1,1))/2,'k-o','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,(Final_Size(1,:,1,2)+Final_Size(2,:,1,2))/2,'b-o','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,(Final_Size(1,:,2,1)+Final_Size(2,:,2,1))/2,'k-+','LineWidth',2,'MarkerSize',8)
plot(0:.1:1,(Final_Size(1,:,2,2)+Final_Size(2,:,2,2))/2,'b-+','LineWidth',2,'MarkerSize',8)

yyy = ylabel('Total Outbreak Size','FontSize',16);
xxx = xlabel('Adopted awareness at Locality 2 ($1-\omega_{22}$)','FontSize',16);
set(xxx,'Interpreter','latex')
ax = gca;
set(ax,'FontSize',16);
lll = legend('Low mobility, weak response at Locality 1','Low mobility, strong response at Locality 1',...
    'High mobility, weak response at Locality 1','High mobility, strong response at Locality 1');
legend boxoff
set(lll,'Interpreter','latex')
set(lll,'FontSize',16);

axis([0 1 0.0 0.7]);
set(gca,'XTick',0:0.2:1);