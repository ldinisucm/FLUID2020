function [lambda,xend] = SimulationNumeriqueLDAverageSlope(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar)
%function [lambda] = SimulationNumerique(pi1, pi2, kappa1, kappa2, k1, k2, n_env, T_max)
% SimulationNumerique return the lyapunov exponent for a 2 states stochastic
% model
% pi1, pi2, kappa1, kappa2, kA1.. sont les paramètres de la simulation
% T_max est le temps que durera la simulation
% The function also returns x, population fractions, to check distributions

%% Variables
x = zeros(1,2); % Vecteur pourcentage de population
x(:) = 0.5; % Pourcentages initiaux
%x(1) = 1; % Pourcentages initiaux
%x(2) = 0; % Pourcentages initiaux

A = zeros(2,2,2); % Matrice de croissance/transition
%%A(1,:,:) = [k1-pi1, pi2; pi1, -pi2+0.2]; % Dans l'environnement 1
%%A(2,:,:) = [-2.-pi1, pi2; pi1, k2-pi2]; % Dans 2
A(1,:,:) = [kA1-pi1, pi2; pi1, -pi2+kB1]; % Dans l'environnement 1
A(2,:,:) = [-pi1+kA2, pi2; pi1, kB2-pi2]; % Dans 2
%Luis Dinis. Diagonalize outside of the loop, "without the t"
B_1 = reshape(A(1,:,:),[2,2]);
B_2 = reshape(A(2,:,:),[2,2]);
[V1,D1]= eig(B_1);
[V2,D2]= eig(B_2);
V(:,:,1)=V1;
V(:,:,2)=V2;
D(:,:,1)=D1;
D(:,:,2)=D2;
%Using the last index  for S avoids the use of reshape, I think, because
%D(:,:,1) is already a 2x2 matrix

lim = 500/max(kA1,kB2); % Limite en Ts d'overflow de 500/max(k1,k2), valeur
% de 500 est empirique

%% Code
% Tirage des temps que dureront chacun des environnements 1 et 2
%Ts_liste1 = exprnd(1/kappa1,1,floor(n_env/2));
%Ts_liste2 = exprnd(1/kappa2,1,floor(n_env/2));
%Ts_liste1=load('Ts1.dat');
%Ts_liste2=load('Ts2.dat');
%Ts_liste = zeros(1,n_env);
% On alterne les environnements 1 et 2 de façon périodique
%Ts_liste(1:2:end) = Ts_liste1;
%Ts_liste(2:2:end) = Ts_liste2;

%Modify the construction of Ts_liste, so that it always ends in the next change after T_max. We can then get rid of n_env
Ts_liste=[];
Ts_liste(1)=exprnd(1/kappa1);
Ts_liste(2)=exprnd(1/kappa2);
t_env=cumsum(Ts_liste);
j=2;
while t_env(j)<T_max
  Ts_liste(j+1)=exprnd(1/kappa1);
  t_env(j+1) = t_env(j) + Ts_liste(j+1);
  Ts_liste(j+2)=exprnd(1/kappa2);
  t_env(j+2)= t_env(j+1) + Ts_liste(j+2);
%  t_env(j+2)
  j=j+2;
end
%t_env(end)

%Alternative, just to check

%Dt=1/max(kappa1,kappa2)/100;
%t=[0:Dt:T_max];
%S=1;
%kappa=kappa1;
%t_env=[];
%for i=1:length(t)
%	r=rand();
%	if r<kappa*Dt
%		%annotate the time of change
%		t_env=[t_env t(i)];
%		%change environment
%		if S==1
%			kappa=kappa2;
%		else
%			kappa=kappa1;
%		end
%                S=3-S;
%	end
%end
%%t_env(1:10)
%%var(t_env)
%index=length(t_env);
		


% On calcule la somme cumulée des temps, puis l'on coupe les listes à T_max
%n_env
%T_max
index = find(t_env>T_max);

%t_env(end)
%T_max
% n_max sera le tour de boucle nécessaire pour atteindre T_max
n_max = index(1);
t_env = t_env(1:n_max);
Ts_liste = Ts_liste(1:n_max);
t_env(end) = T_max;
Ts_liste(end) = t_env(end)-t_env(end-1);
x_tot = zeros(1,n_max); % Liste de la populatino totale

S = 1; % Index de l'environnement actuel
for n=1:n_max
    Ts = Ts_liste(n);
    %B_tmp = reshape(A(S,:,:),[2,2]);
    % Les lignes 43 et 44 permettent de calculer expm(B_tmp) de façon moins
    % couteuse que la fonction expm car les vecteurs et valeurs propres de
    % B_tmp sont bien définies
    
    %[V,D] = eig(B_tmp);
    %B = V*diag(exp(diag(D*Ts)))/V;
    %Luis Dinis
    B = V(:,:,S)*diag(exp(diag(D(:,:,S)*Ts)))/V(:,:,S);
    % En cas d'overflow on découpe le temps Ts trop grand en petit temps
    % Ts2 qui ne provoquent pas d'overflow
    if isinf(B(1,1)) || isinf(B(2,2)) || isinf(B(1,2)) || isinf(B(2,1))
        ratio = fix(Ts/lim) + 1; % On divise Ts
        Ts2 = Ts/ratio;
        x_tmp = 0; % Variable temporaire de log(x_tot)
        for i=1:ratio % On calcule l'évolution sur des petits Ts
	    %Luis Dinis, same trick here
            %B_tmp = reshape(A(S,:,:),[2,2]);
            %[V,D] = eig(B_tmp);
            %B = V*diag(exp(diag(D*Ts2)))/V;
	    B = V(:,:,S)*diag(exp(diag(D(:,:,S)*Ts2)))/V(:,:,S);
	    disp('Hola!') %Uncomment to check if we pass through here.  With RunDL parameters there is no overflow
            x(:) = B*x(:);
            x_tmp = x_tmp + log(sum(x(:)));
            x(:) = x(:)/sum(x(:));
        end
        x_tot(n+1) = x_tot(n) + x_tmp;
        S = 3-S; % On inverse l'environnement
        continue
    end
    x(:) = B*x(:); % Évolution de la population après Ts
    x_tot(n+1) = x_tot(n) + log(sum(x(:))); % Évolution du log(population)
    x(:) = x(:)/sum(x(:)); % Renormalisation
    S = 3-S;
end

lambda = (x_tot(end)-x_tot(1))/(t_env(end)-t_env(1));
% Pour des raisons de longueur de liste, on enlève la valeur initiale
x_tot = x_tot(2:end); 
xend = x;
%% Extraction des données
%%fit = polyfit(t_env,x_tot,1);
%%lambda = fit(1); % Exposant de Lyapunov
%% Hey!!! The fit gives a higher variance!!
if pintar
  %figure
  plot(t_env,x_tot)
  xlabel('t')
  ylabel('log(N(t))')
  hold on
  %plot(t_env,fit(1).*t_env+fit(2))
end
end
