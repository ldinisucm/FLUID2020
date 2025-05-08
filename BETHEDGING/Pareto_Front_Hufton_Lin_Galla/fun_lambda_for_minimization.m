function [lambda_theor]=fun_lambda(pis,extra_param)

%% Paramètres
% kappa1, kappa2, kA1.. sont les paramètres de la simulation

%% extract parameters from input
kappa1=extra_param(1);
kappa2=extra_param(2);
pi1=pis(1);
pi2=pis(2);
kA1=extra_param(3);
kB1=extra_param(4);
kA2=extra_param(5);
kB2=extra_param(6);
numpuntos_ave=extra_param(7);

%Compute theoretical distribution according to Jérémie's notes
%compute phi plus minus
%Delta1=kA1-kA2; %wrong, but anyway kA2=0,kB1=0 then nothing changed.
%Delta2=kB1-kB2;
Delta1=kA1-kB1; %Correct
Delta2=kA2-kB2;
phiv1plus=(Delta1-pi1-pi2+sqrt((Delta1-pi1-pi2)^2+4*pi2*Delta1))/(2*Delta1);
phiv1minus=(Delta1-pi1-pi2-sqrt((Delta1-pi1-pi2)^2+4*pi2*Delta1))/(2*Delta1);
phiv2plus=(Delta2-pi1-pi2+sqrt((Delta2-pi1-pi2)^2+4*pi2*Delta2))/(2*Delta2);
phiv2minus=(Delta2-pi1-pi2-sqrt((Delta2-pi1-pi2)^2+4*pi2*Delta2))/(2*Delta2);

phi1plus=max(phiv1plus,phiv1minus);
phi2plus=min(phiv2plus,phiv2minus);

phi1minus=min(phiv1plus,phiv1minus);
phi2minus=max(phiv2plus,phiv2minus);

g=kappa1/Delta1/(phi1plus-phi1minus);
h=-kappa2/Delta2/(phi2minus-phi2plus);
Q1=kappa2/(kappa1+kappa2);
Q2=kappa1/(kappa1+kappa2);

%Normalization of P1, definition of p11plus, the regular part
%numpuntos=4000; %number of points for integration
%numpuntos_ave=1e6; %number of points for integration
x=linspace(0,(phi1plus-phi2plus)^g-eps,numpuntos_ave);
Dx=x(2)-x(1);
integrand1=(phi1plus-x.^(1/g)-phi2plus).^(h).*(phi1plus-phi1minus-x.^(1/g)).^(-g-1).*(phi2minus+x.^(1/g)-phi1plus).^(-h);
%figure
%plot(x,(phi1plus-x.^(1/g)-phi2plus),x,(phi1plus-phi1minus-x.^(1/g)),x,(phi2minus+x.^(1/g)-phi1plus))
int1=Dx*sum(integrand1);
N1=g*Delta1*Q1/int1;
p11plus=N1/Delta1*integrand1;
P1=x.^(-1/g+1).*p11plus;
%figure
%plot(x,p11plus)
%title('p11plus')
Norm_P1=Dx*sum(p11plus)/g;
%figure
%plot(x,P1)
phi=phi1plus-x.^(1/g);
%plot(phi,P1)
%hold on


%Normalization of P2, definition of p22plus the regular part
%numpuntos=1000; %number of points for integration
x=linspace(0,(phi1plus-phi2plus)^h,numpuntos_ave);
Dx=x(2)-x(1);
integrand2=(phi1plus-x.^(1/h)-phi2plus).^g.*(phi2plus-phi1minus+x.^(1/h)).^(-g).*(phi2minus-x.^(1/h)-phi2plus).^(-h-1);
int2=Dx*sum(integrand2);
N2=h*abs(Delta2)*Q2/int2;
p22plus=N2/abs(Delta2)*(phi1plus-x.^(1/h)-phi2plus).^g.*(phi2plus-phi1minus+x.^(1/h)).^(-g).*(phi2minus-x.^(1/h)-phi2plus).^(-h-1);
P2=p22plus.*x.^(-1/h+1);
Norm_P2=Dx*sum(p22plus)/h;
Norm_P=Norm_P1+Norm_P2;
%plot(x,P2);
phi=phi2plus+x.^(1/h);
%plot(phi,P2)
%legend('P1','P2')
%title('P1 and P2')

%figure
%plot(x,p22plus)
%title('p22plus')
%Average growth computation
%first part
x_bar=(phi1plus-phi2plus)^g;
x=linspace(0,x_bar-eps,numpuntos_ave);
Dx=x(2)-x(1);
p11plus=N1/Delta1*(phi1plus-x.^(1/g)-phi2plus).^h.*(phi1plus-phi1minus-x.^(1/g)).^(-g-1).*(phi2minus+x.^(1/g)-phi1plus).^(-h);
phi=phi1plus-x.^(1/g);
lambda1=1/g*sum(p11plus.*(kA1*phi+kB1*(1-phi)))*Dx;

%second part
x_bar=(phi1plus-phi2plus)^h;
x=linspace(0,x_bar,numpuntos_ave);
Dx=x(2)-x(1);

p22plus=N2/abs(Delta2)*(phi1plus-x.^(1/h)-phi2plus).^g.*(phi2plus-phi1minus+x.^(1/h)).^(-g).*(phi2minus-x.^(1/h)-phi2plus).^(-h-1);
phi=phi2plus+x.^(1/h);
lambda2=1/h*sum(p22plus.*(kA2*phi+kB2*(1-phi)))*Dx;

lambda=lambda1+lambda2

lambda_theor=-lambda;
end
