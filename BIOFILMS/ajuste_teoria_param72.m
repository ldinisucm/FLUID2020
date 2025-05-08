%cargar los datos experimentales
tiempo=4;%12 24 48 72 96
datos_exp_fun('Static',tiempo)

%omega=[1:10 10:5:100 100:10:200 200:50:300];
omega=[0.01 0.1 1:10 10:5:100 100:10:200];
lambda=.00020;

%lambda=.12;
eta=0.02;
x=1.1;
wc=40;
A=0.0001;
B=0.000005
C=6
%beta=0.25;
%Gprima
%gp=A*(lambda^2*omega.^(x-1).*exp(-omega.^2/wc^2));
%gp_alt=A*(lambda^2*omega.^(x-1).*exp(-omega/wc))
%gp_alt2=A*(lambda^2*omega.^(x-1).*exp(-beta*(omega-wc).^2));
X=omega./sinh(pi*omega/wc);
%gp_alt2=A*lambda^2*omega.^(x-1).*besseli(0,X);
%gp_alt2=A*(lambda^2*omega.^(x-1).*exp(omega-wc).^(-1.5));
gp=A*omega.^(x-1).*besseli(0,X)


%Gprimaprima
%gpp=A*(eta*omega+lambda*omega.^(x-1).*exp(-omega.^2/wc^2));
%gpp=A*(eta*omega+lambda*omega.^(x-1).*besseli(0,X));
gpp=C*(eta*omega+B*omega.^(x-1).*besseli(0,X));
%figure
%loglog(omega,gp,omega,gpp,omega,gp_alt,omega,gp_alt2)
loglog(omega,gp,omega,gpp)
%xlabel('\omega')
%ylabel('G')
%legend('G''','G"')

%x=1.05
%wc=10
%tau0=1/wc;
%%tau0=0.01;
%a=0;
%a_alt=0;
%b=0;
%N=1e6;
%omega=10.^[-2:.1:2];
%for i=1:length(omega)
%  w=omega(i);
%  %u=[tau0*w:Du:10000*tau0*w];
%  u=linspace(tau0*w,100000*tau0*w,N);
%  Du=(max(u)-min(u))/N;
%int1=Du*sum(u.^(2-x)./(1+u.^2));
%int2=Du*sum(u.^(1-x)./(1+u.^2));
%G0=1;
%Z=1/G0;
%a(i)=int1*w^(x-1)/Z;
%a_alt(i)=int1*w^(x-1)/Z.*exp(-w^2/wc^2);
%b(i)=int2*w^(x-1)/Z;
%%a(i)=int
%end
%%a=a/a(1);
%datos_exp_fun('Static',tiempo)
%hold on
%loglog(omega,a,omega,b,omega,a_alt)
%%%loglog(omega,a)
%%ratio_at_w0_num=b(1)/a(1)
%%ratio_at_w0_teor=pi/2*(x-1)
