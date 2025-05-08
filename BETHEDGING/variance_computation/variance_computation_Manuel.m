
%Compute variance for different times and parameters from simulations and compare with Jeremie expression.

varteor=[];
varsim=[];
lambdateor=[];
lambdasim=[];
times=[];


set='a'

kappasandpis={'fig2a','fig2b','fig2c'}
kappasandpis={'MinExt'}
% kappasandpis{1} 

for i=1:length(kappasandpis)
kappa_pi_set=kappasandpis{i}

switch kappa_pi_set
case {'MinExt'}
  kappa1=1.0;kappa2=1.0;
  pi1=0.7;
  pi2=0.1;
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
        pi1=0.4;pi2=0.4;
case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
  pi1=0.24;
  pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
  pi1=0.064;pi2=0.064;
case {'DL'}
kappa1 = 0.1; kappa2 = 0.33;
pi1 = kappa1; pi2 = kappa2;
end

switch set
case {'a'}
%set a) Hufton
kA1 = 2; kB2 = -0.2;
kA2=-2; kB1=0.2;
case {'b'}
%set b) Hufton
kA1 = 0.5;
kB1=0.0001;
kA2=0.0001;
kB2=0.3250;
otherwise
        kA1=2;
        kB2=0.2;
        kA2=0;
        kB1=0;
end

%call the function
%pi1=0.3841
%pi2=0.2777
kappa_pi_set=[kappa1 kappa2 pi1 pi2];
ks=[kA1 kB1 kA2 kB2];

T_max_vec=[10 11 13 15:5:20 30:15:100 110:30:400];
numtimes=length(T_max_vec);

%numtimes=1;

%Theory computation
[varteor(i,1),lambdateor(i,1)]=fun_variance_computation(kappa_pi_set,ks)
varteor(i,1:numtimes)=varteor(i,1)*ones(1,numtimes);
lambdateor(i,1:numtimes)=lambdateor(i,1)*ones(1,numtimes);

%Simulations
num_simulations_max=3200; %Number of realizations of simulations for each T_max
num_simulations_min=800;
for j=1:numtimes
 % T_max=j/numtimes*80/kappa1
%  n_env=T_max*10
T_max=T_max_vec(j)
times(j)=T_max;
i,j
num_simulations=num_simulations_min+j/numtimes*(num_simulations_max-num_simulations_min)
  %Simulations
  [lambdasim(i,j),varsim(i,j)]=fun_simulacion_numerica(kappa_pi_set,ks,T_max,num_simulations)
  %[lambdasim(i,j),varsim(i,j)]=fun_simulacion_numerica(kappa_pi_set,ks,n_env,T_max,num_simulations)
end 

figure
rojo_fuerte=[128 0 0]/255;
plot(times,varteor(i,:),'--',times,varsim(i,:),'o','LineWidth',1.2,'MarkerEdge',rojo_fuerte,'MarkerFace','r')
xlabel('t (simulation time)')
ylabel('Growth rate variance')
%title(['kappa&pi=' kappasandpis{i} ' , set=' set])
%legend('Var J','VarSimul*Tmax')

end

figure
plot(times,varteor,'--',times,varsim,'o')
xlabel('t (simulation time)')
ylabel('Growth rate variance')
%legend('VarTheorJ','VarTheorJ','VarSimPar1*Tmax','SimPar2*Tmax')
%title(['set=' set]);

figure
plot(times,lambdateor,'--',times,lambdasim,'o')
xlabel('Tmax (simulation time)')
ylabel('lambda')
title(['kappa&pi=' kappasandpis{i} ' , set=' set])



relative_error_percent_growthrate=abs(lambdasim-lambdateor)./lambdateor*100
average_relerror_growthrate=mean(mean(relative_error_percent_growthrate))
relative_error_percent_variance=abs(varsim-varteor)./varteor*100

average_relerror_variance=mean(mean(relative_error_percent_variance))
