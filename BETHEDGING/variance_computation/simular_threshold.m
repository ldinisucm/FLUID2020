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

num_realizations=4000
%optimal for DLb according to Pareto_front.m
%pi1=3.6959879e-02
%pi2=1.9757696e-01
%optimal for DLa according to Pareto_front.m
pi1=5.5491742e-03
pi2=3.2720220e-02
%optimal for fig2b according to Pareto_front.m
pi1=2.4715947e-01
pi2=2.1983367e-01
T_max=100;
pintar=1;
color='g'
xend1=zeros(1,num_realizations);
thres_ext=-5;
figure
hold on
extinciones1=0;
for i=1:num_realizations
[lambda1,x_tot] = SimNum_forthreshold(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar,color);
i
xend1(i)=x_tot(end);
if any(x_tot<thres_ext)
        extinciones1=extinciones1+1;
end

end
pext1=extinciones1/num_realizations;
%Sub optimal for DLb according to Pareto_front.m
%pi1=8.3447644e-02
%pi2=1.3805260e-01
%pi1=7.3894776e-02
%pi2=1.9334449e-01
%Sub optimal for DLa according to Pareto_front.m
pi1=2.4409171e-02
pi2=1.9169336e-02
%optimal for fig2b according to Pareto_front.m
pi1=3.8022976e-01
pi2=2.4027120e-01
color='r'
xend2=zeros(1,num_realizations);
extinciones2=0;
for i=1:num_realizations
[lambda1,x_tot] = SimNum_forthreshold(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar,color);
xend2(i)=x_tot(end);
if any(x_tot<thres_ext)
	extinciones2=extinciones2+1;
end

end
pext2=extinciones2/num_realizations
%Compute probabilities
threshold=linspace(min([xend1 xend2])*0.99,max([xend1 xend2])*1.1,100);

for i=1:100
	thres=threshold(i);
fraction1(i)=sum(xend1>thres)/num_realizations;

	fraction2(i)=sum(xend2>thres)/num_realizations;
end
figure
plot(threshold,fraction1,'g',threshold,fraction2,'r')
%analyze distribution
std1=std(xend1)
std2=std(xend2)
figure
h1=histogram(xend1,20);
hold on
h2=histogram(xend2,20);
legend('green','red')
figure
plot(h1.BinEdges(2:end),cumsum(h1.Values),'g')
hold on
plot(h2.BinEdges(2:end),cumsum(h2.Values),'r')

%Extinction. It's the same except you don't recover once you're extinct
pext1
pext2
r=pext1/pext2
