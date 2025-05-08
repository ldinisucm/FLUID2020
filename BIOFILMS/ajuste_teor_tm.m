%cargar los datos experimentales
tiempo=5;%12 24 48 72 96
datos_exp_fun('Static',tiempo)

%omega=[1:10 10:5:100 100:10:200 200:50:300];
omega=[0.01 0.1 1:10 10:5:100 100:10:200];


x=1.01
wc=1
tau0=10
a=0;
atc=00;
a_dft=0;
aheavi=0;
tc=.01;
b=0.2
wc=0
N=1000;
G0=.01;
Gprima_num=0;
dft_num=0;
omega=10.^[-2:.01:2];
m=1.5
for i=1:length(omega)
  w=omega(i);
  %u=[tau0*w:Du:10000*tau0*w];
  tau=linspace(tau0,100*tau0,N);
  u=linspace(tau0*w,100*tau0*w,N);
  Du=(max(u)-min(u))/N;
  Dtau=(max(tau)-min(tau))/N;
  %int_tau=Dtau*sum((w^2*tau.^(2-x))./(1+w^2*tau.^2));
  %int_tau=Dtau*sum((w^2*tau.^(2-x))./(1+w^2*tau.^2));
  %int_tau=Dtau*sum(tau.^(-x).*real(j*w*j^(2*m)*(-1)^m*(1./tau+j*w).^(-m-1)));
%  int_tau=Dtau*sum(tau.^(-x).*real(w*tc*tau./(tc*tau*w+j*(tc+tau))));
%  int_tau=Dtau*sum(tau.^(-x).*real(j*w*tau.^2./(-tc*tau.^2*w^2+2*j*tc*tau*w+tc+j*tau.^2*w+tau)));
%  int_tau=Dtau*sum(tau.^(-x).*real(w*(tau./(j+tau*w)-tc*b*tau./(tc*tau*w+j*(tc+tau)))));
%o  int_tau=Dtau*sum(tau.^(-x).*real(j*w*tau.*(-j*tc*tau*w+tc+b*tau+tau)./(tau*w+j)./(tc*(tau*w+j)+j*tau)));
int_tau=Dtau*sum(tau.^(-x*b*(w-wc)^m/(b*(w-wc)^m+x)).*real(j*w*tau./(1+j*w*tau)));
int_tau=Dtau*sum(tau.^(-x-b*(w-wc)^m).*real(j*w*tau./(1+j*w*tau)));
 % int_heavi=Dtau*sum(exp(-tc./tau).*(w^2*tau.^(2-x).*cos(w*tc)+w*tau.^(1-x)*sin(w*tc))./(1+w^2*tau.^2));
int1=Du*sum(u.^(2-x)./(1+u.^2));
%int2=Du*sum(u.^(1-x)./(1+u.^2));
Z=1/G0;
a(i)=int1*w^(x-1)/Z;
atc(i)=int_tau/Z;
%aheavi(i)=int_heavi/Z;
end
loglog(omega,a,omega,atc)
