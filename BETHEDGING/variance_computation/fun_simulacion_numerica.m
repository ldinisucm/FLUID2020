function [lambda_sim,Var_simul]=fun_simulacion_numerica(kappa_pi_set,ks,T_max,num_simulations)

%% Paramètres
% kappa1, kappa2, kA1.. sont les paramètres de la simulation
% n_sim est le nombre de simulations
% T_max est le temps que durera la simulation
% pia, pib sont les extrema de la carte pour pi
% n_points est le nombre de points pour un axe, la carte aura donc une
% taille de n_points*n_points
% NumReal number of realizations for the distribution P histogram 

%% extract parameters from input
kappa1=kappa_pi_set(1);
kappa2=kappa_pi_set(2);
pi1=kappa_pi_set(3);
pi2=kappa_pi_set(4);
kA1=ks(1);
kB1=ks(2);
kA2=ks(3);
kB2=ks(4);

%n_env = 500;

%Growth rate by simulation. build a histogram
disp('Running simulations...')
lambda_simulation=zeros(1,num_simulations);
pintar=1;%Do not plot simulations
for i=1:num_simulations
 [lambda_aux,xend]= SimulationNumeriqueLDAverageSlope(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar);
 lambda_simulation(i)=lambda_aux;
end 

lambda_simulation_ave=mean(lambda_simulation);
lambda_simulation_var=var(lambda_simulation);

T_times_var_simulation=T_max*lambda_simulation_var;

%output values
Var_simul=T_times_var_simulation;
lambda_sim=lambda_simulation_ave;
