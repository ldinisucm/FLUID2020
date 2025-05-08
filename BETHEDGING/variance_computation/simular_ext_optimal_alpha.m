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
thres_ext_vector=-0.5;
%thres_ext_vector=[0 -5 -10]
numpoints=length(thres_ext_vector); %number of different thresholds to simulate
pext1_vs_thres=zeros(size(thres_ext_vector));
pext2_vs_thres=zeros(size(thres_ext_vector));
r=pext2_vs_thres;



verde=[0 0.5 0];
num_realizations=5000

M=load('pareto_fig2b_maspuntos_01.dat');
alpha_vec=M(:,1);
pi1_vec=M(:,2);
pi2_vec=M(:,3);

%alpha_vec=alpha_vec(1:2)

pext_vs_alpha=zeros(size(alpha_vec));

%optimal for fig2b for paper
%pi1=2.6278436e-01
%pi2=2.4628221e-01
%color=verde
%sub-optimal for fig2b for paper
%pi1=3.4686696e-01
%pi2=2.5471285e-01
%color='r'
%%Small growth and variance "maroon/null"
%pi1=6.8948419e+00
%pi2=1.2863575e-03
%color=[0.5 0 0]
% "middle/magenta"
%pi1=5.6682983e-01
%pi2=2.0706331e-01
%color='m'



for k=1:length(alpha_vec)
pi1=pi1_vec(k);
pi2=pi2_vec(k);
color='k'


T_max=500;
%T_max=50;
pintar=0;
xend1=zeros(1,num_realizations);
%figure(1)
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
pext_vs_alpha(k)=pext1_vs_thres
end

%figure(2)
%hold on
%plot(thres_ext_vector,pext1_vs_thres,'-','color',color)
%xlabel('E')
%ylabel('Prob. of ext.')

figure
[a,indices]=sort(alpha_vec)
pext=pext_vs_alpha(indices);
%pext=pext(1:2:end)
%alpha_pintar=a(1:2:end)
alpha_pintar=a;
plot(alpha_pintar,pext,'o')
