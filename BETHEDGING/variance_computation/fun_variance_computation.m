function [Var_theor,lambda_theor]=fun_variance_computation(kappa_pi_set,ks)

%% Paramètres
% kappa1, kappa2, kA1.. sont les paramètres de la simulation
% n_sim est le nombre de simulations
% pia, pib sont les extrema de la carte pour pi
% n_points est le nombre de points pour un axe, la carte aura donc une
% taille de n_points*n_points
% NumReal number of realizations for the distribution P histogram 

%% extract parameters from input
kappa1=kappa_pi_set(1);
kappa2=kappa_pi_set(2);
pi1=kappa_pi_set(3);
pi2=kappa_pi_set(4);
kA1=ks(1);
kB1=ks(2);
kA2=ks(3);
kB2=ks(4);

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
numpuntos_ave=1e6; %number of points for integration
numpuntos_var=2000; %number of points for integration
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
Norm_P1=Dx*sum(p11plus)/g
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
Norm_P2=Dx*sum(p22plus)/h
Norm_P=Norm_P1+Norm_P2
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

%Variance computation

%first part
phi_bar=0.5*(phi2plus+phi1plus);
x_bar=(phi_bar-phi2plus)^h;
y=linspace(0,x_bar,numpuntos_var);
Dy=y(2)-y(1);

p12plus=N1/Delta1*(phi1plus-phi2plus-y.^(1/h)).^(g-1).*(y.^(1/h)+phi2plus-phi1minus).^(-g-1).*(phi2minus-y.^(1/h)-phi2plus).^(-h);
p22plus=N2/abs(Delta2)*(phi1plus-y.^(1/h)-phi2plus).^g.*(phi2plus-phi1minus+y.^(1/h)).^(-g).*(phi2minus-y.^(1/h)-phi2plus).^(-h-1);
phi=phi2plus+y.^(1/h);
phi1=phi;
v1=Delta1*phi.*(1-phi)-pi1*phi+pi2*(1-phi);
Gofy=(kappa1*y.^(1/h).*p12plus+kappa2*p22plus)./(v1.*p12plus).^2;
%Dt=0.001;
Dt=1/numpuntos_var;
t=0:Dt:1;
%figure
%plot(y,p12plus,y,p22plus)
%legend('p12plus','p22plus')
for i=1:length(y);
  yprime=y(i);
  z=t*yprime;
  p12plus=N1/Delta1*(phi1plus-phi2plus-z.^(1/h)).^(g-1).*(z.^(1/h)+phi2plus-phi1minus).^(-g-1).*(phi2minus-z.^(1/h)-phi2plus).^(-h);
  p22plus=N2/abs(Delta2)*(phi1plus-z.^(1/h)-phi2plus).^g.*(phi2plus-phi1minus+z.^(1/h)).^(-g).*(phi2minus-z.^(1/h)-phi2plus).^(-h-1);
  f1=kA1*(z.^(1/h)+phi2plus)+kB1*(1-z.^(1/h)-phi2plus)-lambda;
  f2=kA2*(z.^(1/h)+phi2plus)+kB2*(1-z.^(1/h)-phi2plus)-lambda;
	Fofz=z.^(1/h).*p12plus.*f1+p22plus.*f2;
  Jofy(i)=(sum(Fofz)*Dt)^2;
end
V1=1/h^3*sum(Gofy.*Jofy)*Dy
%suma=sum(Gofy.*Jofy)*Dy


%second part
x_bar=((phi1plus-phi2plus)*0.5)^g;
x=linspace(0,x_bar-eps,numpuntos_var);
Dx=x(2)-x(1);
p11plus=N1/Delta1*(phi1plus-x.^(1/g)-phi2plus).^h.*(phi1plus-phi1minus-x.^(1/g)).^(-g-1).*(phi2minus+x.^(1/g)-phi1plus).^(-h);
p21plus=N2/abs(Delta2)*(phi1plus-x.^(1/g)-phi2plus).^(h-1).*(phi1plus-x.^(1/g)-phi1minus).^(-g).*(phi2minus+x.^(1/g)-phi1plus).^(-h-1);
phi=phi1plus-x.^(1/g); %change of variables
phi2=phi;
%v1=Delta1*phi.*(1-phi)-pi1*phi+pi2*(1-phi);
v1=Delta1*(phi1plus-x.^(1/g)).*(1-(phi1plus-x.^(1/g)))-pi1*(phi1plus-x.^(1/g))+pi2*(1-(phi1plus-x.^(1/g)));
Gofx=x.^(2/g).*( kappa1*p11plus+kappa2*p21plus.*x.^(1/g)  )./(v1.*p11plus).^2;
%Alternative Gofx, using v1=Delta1*x^(1/g)(phi1plu-x^(1/g))
%Gofx=( kappa1*p11plus+kappa2*p21plus.*x.^(1/g) )./(Delta1.*(phi1plus-x.^(1/g)-phi1minus).*p11plus).^2;
v2=Delta2*phi.*(1-phi)-pi1*phi+pi2*(1-phi);
%Jeremie's version, gets the same result as "my Alternative Gofx"
Gofx=( kappa1*p11plus+kappa2*p21plus.*x.^(1/g) )./(v2.*p21plus).^2;
%Dt=0.001;
Dt=1/numpuntos_var;
t=0:Dt:1;
%figure
%plot(y,p11plus,y,p21plus)
%legend('p11plus','p21plus')
for i=1:length(x);
xprime=x(i);
z=t*xprime;
p11plus=N1/Delta1*(phi1plus-z.^(1/g)-phi2plus).^h.*(phi1plus-phi1minus-z.^(1/g)).^(-g-1).*(phi2minus+z.^(1/g)-phi1plus).^(-h);
p21plus=N2/abs(Delta2)*(phi1plus-z.^(1/g)-phi2plus).^(h-1).*(phi1plus-z.^(1/g)-phi1minus).^(-g).*(phi2minus+z.^(1/g)-phi1plus).^(-h-1);
f1=kA1*(phi1plus-z.^(1/g))+kB1*(1-phi1plus+z.^(1/g))-lambda;
f2=kA2*(phi1plus-z.^(1/g))+kB2*(1-phi1plus+z.^(1/g))-lambda;
Fofz2=p11plus.*f1+p21plus.*f2.*z.^(1/g);
Jofx(i)=(sum(Fofz2)*Dt)^2;
end
V2=1/g^3*sum(Gofx.*Jofx)*Dx
%figure
%plot(phi1,Jofy)
%hold on 
%plot(phi2,Jofx)
%legend('J1','J2')
%figure
%plot(Fofz)
%hold on
%plot(Fofz2)
%legend('F1','F2')

%figure
%plot(phi1,Gofy)
%hold on
%plot(phi2,Gofx)
%%plot(phi2,Gofxv2p2)
%legend('G1','G2')
%%legend('G1','G2','G2v2p2')
%
%figure
%plot(phi1,Jofy.*Gofy)
%hold on
%plot(phi2,Gofx.*Jofx)
%legend('G1*J1','G2*J2')



%%KL approx
%t2 = 0.1e1 / (kappa1 + kappa2);
%t3 = t2 * kappa2;
%t4 = (kA1 ^ 2);
%t11 = (kB1 ^ 2);
%t16 = (pi1 ^ 2);
%t18 = 2 * pi2 * pi1;
%t19 = pi2 ^ 2;
%t21 = sqrt((-2 * kB1 * kA1 - 2 * kA1 * pi1 + 2 * pi2 * kA1 + 2 * kB1 * pi1 - 2 * kB1 * pi2 + t11 + t16 + t18 + t19 + t4));
%t25 = (kA2 ^ 2);
%t32 = (kB2 ^ 2);
%t38 = sqrt((-2 * kB2 * kA2 - 2 * kA2 * pi1 + 2 * pi2 * kA2 + 2 * kB2 * pi1 - 2 * kB2 * pi2 + t16 + t18 + t19 + t25 + t32));
%t41 = kB2 - pi2 - kA2 + pi1 + t38;
%t44 = kB1 - pi2 - kA1 + pi1 + t21;
%t52 = sqrt(0.1e1 + 0.4e1 / t44 ^ 2 * t19);
%t59 = sqrt(0.1e1 + 0.4e1 / t41 ^ 2 * t19);
%t62 = log(0.1e1 / t59 / t52 * (0.1e1 + 0.4e1 / t44 / t41 * t19));
%t65 = (kB1 - pi2 + kA1 - pi1 + t21) * t3 / 0.2e1 + (kB2 - pi2 + kA2 - pi1 + t38) * t2 * kappa1 / 0.2e1 + t3 * kappa1 * t62;
%lambda_KL=t65


%disp('Results')
%lambda_KL

%output values
Var_theor=V1+V2;
lambda_theor=lambda;
