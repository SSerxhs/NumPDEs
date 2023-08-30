n=6;
x=zeros(n,1);
y0=zeros(n,1);
y1=zeros(n,1);
y2=zeros(n,1);
y3=zeros(n,1);
y4=zeros(n,1);
for i=1:n
    x(i)=1/(n-1)*(i-1);
end
for i=1:n
	y0(i)=sin(pi*1*x(i));
	y1(i)=sin(pi*2*x(i));
	y2(i)=sin(pi*3*x(i));
	y3(i)=sin(pi*4*x(i));
	y4(i)=sin(pi*5*x(i));
end
plot(x,y0,'r',x,y1,'g',x,y2,'b',x,y3,'y',x,y4,'k');