%% seir dynamics with social distancing response
function x_integrated = SEIR_j(N,x,migration_in, migration_out,j,awareness_at_neighbors,dt)

% x = [S, E, I,beta,delta,mu, alpha, omega]
horizon = 2;
Susc = zeros(horizon);
Exp = zeros(horizon);
Infec = zeros(horizon);

Susc(1)=x(1); Exp(1)=x(2); Infec(1)=x(3); beta=x(4); delta=x(5); mu=x(6); kappa = x(7); omega = x(8); %new_cases = x(10);
Recov = 1-Susc(1)-Exp(1)-Infec(1);

% (1-omega*(Infec(1)+Recov) - (1-omega)*awareness_at_neighbors^alpha)^kappa
beta_dist = beta*((1-(omega*(Infec(1)+Recov) + (1-omega)*awareness_at_neighbors))^kappa);%
% beta_dist = beta*(Susc(1)+Exp(1))^kappa;

Susc(2) = Susc(1) + dt*(-beta_dist*Susc(1)*(Infec(1)+Exp(1))- migration_out(1)+migration_in(1));
Exp(2) = Exp(1) + dt*(beta_dist*Susc(1)*(Infec(1)+Exp(1))-mu*Exp(1)- migration_out(2)+migration_in(2));
Infec(2) = Infec(1) + dt*((mu*Exp(1)/N) - delta*Infec(1)/N);

x_integrated = struct('S',Susc(2),'E', Exp(2), 'I',Infec(2),'beta_dist',beta_dist);
