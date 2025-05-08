%cargar los datos experimentales
tiempo=4;%12 24 48 72 96
[omega_exp,Gp_exp,Gpp_exp]=datos_exp_fun_minimizar('Static',tiempo)
%Comprobamos si hay algún dato negativo y lo quitamos
Gp_prune=0;
wp=0;
Gpp_prune=0;
wpp=0;
jp=1;
jpp=1;
display('Quitando valores negativos')
for i=1:length(Gp_exp)
	if Gp_exp(i)>=0
		Gp_prune(jp)=Gp_exp(i);
		wp(jp)=omega_exp(i);
		jp=jp+1;
	end
end
for i=1:length(Gpp_exp)
        if Gpp_exp(i)>=0
		Gpp_prune(jpp)=Gpp_exp(i);
		wpp(jpp)=omega_exp(i);
		jpp=jpp+1;
	end
end
Gp_exp=log(Gp_prune)'
omega_exp_Gp=(wp)'
Gpp_exp=log(Gpp_prune)'
omega_exp_Gpp=(wpp)'



%Guess inicial de los parámetros
switch tiempo

case 1

%Bueno par 12h
eta=.25;
x=1.24;
fc=90
Gp0=0.00002;
Gpp0=.4;
F=1.5;

case 2
%Buenos para 24h
eta=.6;
x=1.1;
fc=5
Gp0=0.00002;
Gpp0=.4;
F=4;
eta=.9;
x=1.14;
fc=20
Gp0=0.00002;
Gpp0=.17;
F=3.4;

case 3
%Buenos para 48h
eta=0.8;
x=1.09;
fc=4.1
Gp0=0.00002;
Gpp0=.4;
F=10;

eta=0.6;
x=1.11;
fc=9
Gp0=0.00002;
Gpp0=.1;
F=8;
	
case 4
%Para 72h
eta=0.4;
x=1.07;
fc=3.9
Gp0=0.00002;
Gpp0=.24;
F=1.8;

case 5
%Bueno par 96h
eta=0.01;
x=1.1;
fc=1.2
Gp0=0.00002;
Gpp0=1.2;
F=4;
end

parametros=[x fc Gp0  eta Gpp0 F];

tolerancia=1e-8
iter=500000
options=optimset('TolFun',tolerancia,'TolX',tolerancia,'MaxFunEvals',iter,'MaxIter',iter);
options=optimset
parmFit=fminsearch('fun_minimizar_logaritmico_3',parametros,options,omega_exp_Gp,omega_exp_Gpp,Gp_exp,Gpp_exp)

x=parmFit(1);
fc=parmFit(2);
wc=2*pi*fc;
Gp0=parmFit(3);
eta=parmFit(4);
Gpp0=parmFit(5);
F=parmFit(6);

omega=[0.01 0.1 1:10 10:5:100 100:10:200];
X=omega./sinh(pi*omega./wc);
I0=besseli(0,F*X);
Gp_teor=Gp0*omega.^(x-1).*I0.^(-x);

X=omega./sinh(pi*omega./wc);
I0=besseli(0,F*X);
Gpp_teor=Gpp0*(eta*omega+Gp0*omega.^(x-1).*I0.^(-x));


loglog(omega,Gp_teor,'DisplayName','G'' Teor')
loglog(omega,Gpp_teor,'DisplayName','G'''' Teor')
legend('Location','SouthWest')

