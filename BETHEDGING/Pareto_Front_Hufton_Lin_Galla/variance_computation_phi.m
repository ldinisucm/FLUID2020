
%Compute variance for different times and parameters from simulations and compare with direct expression

varteor=[];
varsim=[];
lambdateor=[];
lamdasim=[];
times=[];


set='a'

kappasandpis={'fig2c'}
%kappasandpis={'other'}
 kappasandpis{1} 

for i=1:1
kappa_pi_set=kappasandpis{i}

switch kappa_pi_set
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
case {'other'}
  kappa1=9; kappa2=12;
  pi1=0.4;pi2=0.4;
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
kappa_pi_set=[kappa1 kappa2 pi1 pi2];
ks=[kA1 kB1 kA2 kB2];

numtimes=1;
for j=1:numtimes
T_max=1600/kappa1;
T_max=20/kappa1;
 % T_max=j/numtimes*1600/kappa1;
%times(j)=T_max;
numpuntos(j)=j/numtimes*80000;
i,j
  [varteor(i,j),varsim(i,j),lambdateor(i,j),lambdasim(i,j),varlambda_naive(i,j)]=fun_variance_computation_phi(kappa_pi_set,ks,T_max,numpuntos(j))
end 

figure
%plot(times,varteor(i,:),'--',times,varsim(i,:),'o')
plot(numpuntos,varteor(i,:),'*',numpuntos,varsim(i,:),'o')
%xlabel('Tmax (simulation time)')
xlabel('number of points')
ylabel('Variance')
title(['kappa&pi=' kappasandpis{i} ' , set=' set])
legend('Var J','VarSimul*Tmax')

end

%figure
%plot(times,varteor,'--',times,varsim,'o')
%xlabel('Tmax (simulation time)')
%ylabel('Variance')
%legend('VarTheorJ','VarTheorJ','VarSimPar1*Tmax','SimPar2*Tmax')
%title(['set=' set]);
%
%figure
%plot(times,lambdateor,'--',times,lambdasim,'o')
%xlabel('Tmax (simulation time)')
%ylabel('lambda')

%figure
%plot(times,varteor,'--',times,varsim,'o',times,varlambda_naive,'*')


relative_error_percent_growthrate=abs(lambdasim-lambdateor)./lambdateor*100
average_relerror_growthrate=mean(mean(relative_error_percent_growthrate))
relative_error_percent_variance=abs(varsim-varteor)./varteor*100

average_relerror_variance=mean(mean(relative_error_percent_variance))
