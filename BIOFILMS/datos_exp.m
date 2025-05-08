%condicion='Static' %Static o Shaking
condicion='Shaking' %Static o Shaking
archivo_Elast=[condicion '_Elast_Monroy.csv'];
archivo_Loss=[condicion '_Loss_Monroy.csv'];
A=load(archivo_Elast);
B=load(archivo_Loss);

%A=load('Static_Elast_Monroy.csv');
%B=load('Static_Loss_Monroy.csv');
%A=load('StaticElast.txt');
%B=load('StaticLoss.txt');
%A=load('ShakingElast.txt');
%B=load('ShakingLoss.txt');
figure
welast=A(:,1);
%gprime_96=10.^A(:,6);
tiempo=6; %2:12h 3:24h 4:48h 5:72h 6:96h
gprime_tiempo=A(:,tiempo); %/5
loglog(welast,gprime_tiempo,'-o')
hold on
wloss=B(:,1);
gpp_tiempo=B(:,tiempo);
%gpp_96=10.^B(:,);
loglog(wloss,gpp_tiempo,'-o')
legend('Gp','Gpp')
tiempos=[12 24:24:96];
cadena=[condicion ' ' num2str(tiempos(tiempo-1)) 'h'];
title(cadena)
