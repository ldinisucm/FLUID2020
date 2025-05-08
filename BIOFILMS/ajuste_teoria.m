%cargar los datos experimentales
tiempo=3;%12 24 48 72 96
datos_exp_fun('Static',tiempo)
hold on

%omega=[1:10 10:5:100 100:10:200 200:50:300];
omega=[0.01 0.1 1:10 10:5:100 100:10:200];
%lambda=.00020;

lambda=.00000012;
eta=0.06;
x=1.05;
fc=2.7
wc=2*pi*fc;
A=2000000000;
Gp0=4;
Gpp0=.5;
%B=0.00001
%C=1.2
%F=4
%beta=0.25;
%Gprima
%gp=A*(lambda^2*omega.^(x-1).*exp(-omega.^2/wc^2));
gp=Gp0*omega.^(x-1).*exp(-omega.^2/wc^2);
wc_exp=2*pi*1
gp_alt=Gp0*omega.^(x-1).*exp(-omega/wc_exp);
%gp_alt2=A*(lambda^2*omega.^(x-1).*exp(-beta*(omega-wc).^2));
%X=omega./sinh(pi*omega/wc);
%gp_alt2=A*lambda^2*omega.^(x-1).*besseli(0,X);
%gp_alt2=A*(lambda^2*omega.^(x-1).*exp(omega-wc).^(-1.5));
%gp=A*omega.^(x-1).*besseli(0,F*X)


%Gprimaprima
gpp=(eta*omega+Gpp0*omega.^(x-1).*exp(-omega.^2/wc^2));
%gpp=A*(eta*omega+lambda*omega.^(x-1).*besseli(0,X));
%gpp=C*(eta*omega+B*omega.^(x-1).*besseli(0,F*X));
%figure
%loglog(omega,gp,omega,gpp,omega,gp_alt,omega,gp_alt2)
loglog(omega,gp,omega,gpp,omega,gp_alt)
%xlabel('\omega')
%ylabel('G')
%legend('G''','G"')

%
x=1.05
wc=2*pi*2.7
tau0=1;
%%tau0=0.01;
a=0;
a_alt=0;
b=0;
N=1e6;
G0=0.5;
omega=10.^[-2:.1:2];
for i=1:length(omega)
  w=omega(i);
  %u=[tau0*w:Du:10000*tau0*w];
  u=linspace(tau0*w,100000*tau0*w,N);
  Du=(max(u)-min(u))/N;
int1=Du*sum(u.^(2-x)./(1+u.^2));
int2=Du*sum(u.^(1-x)./(1+u.^2));
Z=1/G0;
a(i)=int1*w^(x-1)/Z;
a_alt(i)=int1*w^(x-1)/Z.*exp(-w^2/wc^2);
b(i)=int2*w^(x-1)/Z+eta*w;
%%a(i)=int
end
%%a=a/a(1);
datos_exp_fun('Static',tiempo)
hold on
loglog(omega,a,omega,b,omega,a_alt)
%%loglog(omega,a)
%%ratio_at_w0_num=b(1)/a(1)
%%ratio_at_w0_teor=pi/2*(x-1)
