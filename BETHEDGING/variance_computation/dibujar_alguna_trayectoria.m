
set='a'

kappasandpis={'fig2a','fig2b','fig2c'}

kappa_pi_set=kappasandpis{2}

switch kappa_pi_set
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
        pi1=0.4;pi2=0.4;
case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
  pi1=0.24;
  pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
  pi1=0.064;pi2=0.064;
case {'DL'}
kappa1 = 0.1; kappa2 = 0.33;
pi1 = kappa1; pi2 = kappa2;
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

%call the function
pi1=0.263
pi2=0.246
%pi1=1
%pi2=0.1
pi1=6
pi2=0.001
T_max=1000;
pintar=1;
SimulationNumeriqueLDAverageSlope(pi1, pi2, kappa1, kappa2, kA1, kB1, kA2, kB2, T_max,pintar);
