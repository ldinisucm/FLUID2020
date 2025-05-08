function color=choose_color(T,Tlimits)
%decides the color for a T
colores=[
0 0 1 %azul
.8500 0.3250 0.0980 %naranj
0.9290 0.6940 0.1250 %ama
0.4940 0.1840 0.5560 %mag
.4660 0.6740 0.1880 %verde
];
[f,c]=size(colores);
if length(Tlimits<=f)
%colores=['rgbcmyk'];
else
	disp('Demasiados límites')
	cuantos=floor(Tlimits/f)+1;
	for j=1:cuantos
		colores=[colores;colores];
	end
end

color=colores(1,:);
for i=1:length(Tlimits)
	if T>Tlimits(i)
		color=colores(i+1,:);
	end
end
