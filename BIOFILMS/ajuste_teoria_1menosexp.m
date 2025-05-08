%cargar los datos experimentales
tiempo=4;%12 24 48 72 96
datos_exp_fun('Static',tiempo)

%omega=[1:10 10:5:100 100:10:200 200:50:300];
omega=[0.01 0.1 1:10 10:5:100 100:10:200];


x=1.2
wc=1
tau0=10
a=0;
atc=00;
aheavi=0;
tc=1;
N=1e5;
G0=5e-3;
G0=1;
omega=10.^[-2:.02:2];
for i=1:length(omega)
  w=omega(i);
  %u=[tau0*w:Du:10000*tau0*w];
  tau=linspace(tau0,200000*tau0,N);
  u=linspace(tau0*w,200000*tau0*w,N);
  Du=(max(u)-min(u))/N;
  Dtau=(max(tau)-min(tau))/N;
  %int_tau=Dtau*sum((w^2*tau.^(2-x))./(1+w^2*tau.^2));
  int_tau=Dtau*sum((w^2*tau.^(4-x).*(1-2*tc./tau))./((tau+tc*tau.^2*w^2-tc).^2+(2*tc*tau*w+w*tau.^2).^2));
  int_heavi=Dtau*sum(exp(-tc./tau).*(w^2*tau.^(2-x).*cos(w*tc)+w*tau.^(1-x)*sin(w*tc))./(1+w^2*tau.^2));
int1=Du*sum(u.^(2-x)./(1+u.^2));
%int2=Du*sum(u.^(1-x)./(1+u.^2));
Z=1/G0;
a(i)=int1*w^(x-1)/Z;
atc(i)=int_tau/Z;
aheavi(i)=int_heavi/Z;
%a_alt(i)=int1*w^(x-1)/Z.*exp(-w^2/wc^2);
%b(i)=int2*w^(x-1)/Z;
%%a(i)=int
end
loglog(omega,a,omega,atc,omega,aheavi)
%%a=a/a(1);
%datos_exp_fun('Static',tiempo)
%hold on
%loglog(omega,a,omega,b,omega,a_alt)
%%%loglog(omega,a)
%%ratio_at_w0_num=b(1)/a(1)
%%ratio_at_w0_teor=pi/2*(x-1)
