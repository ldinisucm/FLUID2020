%Parameter definitions
set='a'


kappa_pi_set='fig2b'

switch kappa_pi_set
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
  %pi1=0.24;
  %pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
  %pi1=0.064;pi2=0.064;
case {'DLa'}
kappa1 = 0.01; kappa2 = 0.033;
case {'DLb'}
kappa1 = 0.1; kappa2 = 0.33;
%pi1 = kappa1; pi2 = kappa2;
case {'DLc'}
kappa1=1; kappa2 = 3.3;  
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

thres_ext_vector=[0:-0.5:-8];
thres_ext_vector=[0:-0.01:-5];
%thres_ext_vector=[0 -5 -10]
numpoints=length(thres_ext_vector); %number of different thresholds to simulate
pext1_vs_thres=zeros(size(thres_ext_vector));
pext2_vs_thres=zeros(size(thres_ext_vector));
r=pext2_vs_thres;



num_realizations=8000
%optimal for DLb according to Pareto_front.m
%pi1=3.6959879e-02
%pi2=1.9757696e-01
%optimal for DLa according to Pareto_front.m
%pi1=5.5491742e-03
%pi2=3.2720220e-02
%optimal for fig2b according to Pareto_front.m
pi1=2.4715947e-01
pi2=2.1983367e-01
%optimal for fig2b for paper
pi1=2.6278436e-01
pi2=2.4628221e-01
T_max=500;
%T_max=50;
pintar=1;
verde=[0 0.5 0];
color=verde
%color='b'
xend1=zeros(1,num_realizations);
figure
hold on
extinciones1=zeros(1,numpoints);
lambda1_sim=zeros(1,num_realizations);
N1=zeros(numpoints,num_realizations);
N1ave_survival=zeros(1,numpoints);
for i=1:num_realizations
[lambda1,x_tot] = SimNum_forthreshold(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar,color);
i
x_tot=[0 x_tot];%First point is removed in simulation program
lambda1_sim(i)=lambda1;
xend1(i)=x_tot(end);
for j=1:numpoints
	thres_ext=thres_ext_vector(j);
if any(x_tot<thres_ext)
        extinciones1(j)=extinciones1(j)+1;
else
 %This one survived. Count the final population
 N1(j,i)=exp(x_tot(end));
end
end
N1ave_survival=mean(N1,2);

end
lambda1_ave=mean(lambda1_sim)
pext1_vs_thres=extinciones1/num_realizations;
%Sub optimal for DLb according to Pareto_front.m
%pi1=8.3447644e-02
%pi2=1.3805260e-01
%pi1=7.3894776e-02
%pi2=1.9334449e-01
%Sub optimal for DLa according to Pareto_front.m
%pi1=2.4409171e-02
%pi2=1.9169336e-02
%slightly suboptimal for fig2b according to Pareto_front.m
%pi1=3.8022976e-01
%pi2=2.4027120e-01
%very small growth and variance
%pi1=2.3433165e+00
%pi2=5.6036966e-03

%sub-optimal for fig2b for paper
pi1=3.4686696e-01
pi2=2.5471285e-01
color='r'
%%Small growth and variance "maroon/null"
%pi1=6.8948419e+00
%pi2=1.2863575e-03
%color=[0.5 0 0]
% "middle/magenta"
%pi1=5.6682983e-01
%pi2=2.0706331e-01
%color='m'
xend2=zeros(1,num_realizations);
extinciones2=zeros(1,numpoints);
lambda2_sim=zeros(1,num_realizations);
N2=zeros(length(thres_ext_vector),num_realizations);
N2ave_survival=zeros(1,numpoints);
for i=1:num_realizations
	i
[lambda2,x_tot] = SimNum_forthreshold(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar,color);
x_tot=[0 x_tot];%First point is removed in simulation program
xend2(i)=x_tot(end);
lambda2_sim(i)=lambda2;
for j=1:numpoints
	thres_ext=thres_ext_vector(j);
        if any(x_tot<thres_ext)
	  extinciones2(j)=extinciones2(j)+1;
        else
        N2(j,i)=exp(x_tot(end));

end
	
end
N2ave_survival=mean(N2,2);

end
lambda2_ave=mean(lambda2_sim)
pext2_vs_thres=extinciones2/num_realizations;
%Compute probabilities
%threshold=linspace(min([xend1 xend2])*0.99,max([xend1 xend2])*1.1,100);

%for i=1:100
%	thres=threshold(i);
%fraction1(i)=sum(xend1>thres)/num_realizations;
%
%	fraction2(i)=sum(xend2>thres)/num_realizations;
%end
%figure

%plot(threshold,fraction1,'g',threshold,fraction2,'r')
%analyze distribution
std1=std(xend1)
std2=std(xend2)
lambda1_ave
lambda2_ave

%figure
%h1=histogram(xend1,20);
%hold on
%h2=histogram(xend2,20);
%legend('green','red')
%figure
%plot(h1.BinEdges(2:end),cumsum(h1.Values),'g')
%hold on
%plot(h2.BinEdges(2:end),cumsum(h2.Values),'r')
xlabel('t')
ylabel('log(N)')
box on

r=pext1_vs_thres./pext2_vs_thres;
figure
%plot(thres_ext_vector,pext1_vs_thres,'MarkerEdge',verde,thres_ext_vector,pext2_vs_thres,'r')
plot(thres_ext_vector,pext1_vs_thres,'-','color',verde)
hold on
plot(thres_ext_vector,pext2_vs_thres,'-r')
xlabel('Threshold for extinction')
ylabel('Prob. of ext.')
figure
plot(thres_ext_vector,(r-1)*100)
xlabel('Threshold for extinction')
ylabel('Percent increase in prob of extinction')
%figure
%plot(thres_ext_vector,1./r)
%figure
%semilogy(thres_ext_vector,pext1_vs_thres,'g',thres_ext_vector,pext2_vs_thres,'r')
%
%xlabel('Threshold for extinction')
%ylabel('Prob. of ext.')
%
%figure
%plot(thres_ext_vector,1-pext1_vs_thres,'g',thres_ext_vector,1-pext2_vs_thres,'r')
%xlabel('Threshold for extinction')
%ylabel('Survival probability')
%figure
%semilogx(abs(thres_ext_vector),1-pext1_vs_thres,'g',abs(thres_ext_vector),1-pext2_vs_thres,'r')
%xlabel('Threshold for extinction')
%ylabel('Survival probability')
%figure
%plot(thres_ext_vector,(1-pext2_vs_thres)./(1-pext1_vs_thres),'o')
%xlabel('Threshold for extinction')
%ylabel('Survival probability ratio Sred/Sgreen')
%

%Compute average population including extinction.
%Try with the first extinction threshold for the moment
%Sum population over realizations that survived.
figure
plot(thres_ext_vector,N1ave_survival,'g',thres_ext_vector,N2ave_survival,'r')


%
%%Extinction. It's the same except you don't recover once you're extinct
