logC0=0;
tray=[0;1]
num_pasos=1
p=.6;
P=[1-p;p] %prob 0, prob 1

odds=[2.3 2]';
b=P;
b=[0.44 1-0.44]
delta_0=log(b(1)*odds(1));
delta_1=log(b(2)*odds(2));


%delta_plus=0.6;
%delta_minus=-0.6;
threshold=-0.001;
if logC0+delta_0<threshold
	P(1)=0;
end
if logC0+delta_1<threshold
	P(2)=0;
end
P

for i=2:num_pasos
	aux=[tray;tray];
	col_nueva=[zeros(length(tray),1);ones(length(tray),1)];
	tray=[col_nueva aux]
	[f,c]=size(tray);
	Pold=P
	P=zeros(f,1);
	for j=1:f
		tray(j,:)	
		logC=logC0+sum(tray(j,:))*delta_1+(c-sum(tray(j,:)))*delta_0
 		if logC>threshold
			index=dec(tray(j,:))
			previous_index=dec(tray(j,1:i-1))
			P(index)=(p*tray(j,i)+(1-p)*(1-tray(j,i)))*Pold(previous_index);
			P(index)
		end
	end
	P
end
Pext=1-sum(P)
