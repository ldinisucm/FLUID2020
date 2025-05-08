
k1=0.5
k2=4
m=1;
A=1.1;

gamma=0.1;
w0=sqrt((k2+k1)/m);
a=k2*A/m;
tau=m/gamma;
w=[0.01:.01:10*w0];
x10=a./sqrt((w0^2-w.^2).^2+w.^2/tau^2)
phi=atan(-w/tau./(w0^2-w.^2))
Tmax=0;
for i=1:length(w)
	t=[0:w(i)*2*pi/1000:2*pi/w(i)];
Tmax(i)=k2*max(A*cos(w(i)*t)-x10(i)*cos(w(i)*t-phi(i)));
end
plot(w,Tmax)
