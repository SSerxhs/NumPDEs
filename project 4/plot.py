from math import sin,cos,exp,acos,pi
import matplotlib.pyplot as plt
n=201
x=[0]*n
y=[0]*n
def phi(x):
	return max(0,1-20*abs(x-0.5))
def A(k):
	k*=pi
	return 40/(k**2)*(-sin(9/20*k)+2*sin(1/2*k)-sin(11/20*k))
def u(x,t):
	r=0
	for k in range(1,100):
		r+=A(k)*exp(-k*k*pi*pi*t)*sin(k*pi*x)
	return r
def f(x):
	return exp(-20*(x-19)**2)+exp(-(x-22)**2)
for i in range(n):
	x[i]=i*10/200-2
	y[i]=f(x[i]+17)
filename="advection initial condition"
plt.xlim(-2,8)
plt.ylim(-0.4,1)
plt.plot(x,y)
plt.subplots_adjust(left=0.15,right=0.98,top=0.92,bottom=0.05)
plt.suptitle(filename, fontsize = 12)
plt.plot(x,y,color='k')
plt.savefig(filename+'.jpg')
plt.clf()
for i in range(n):
	x[i]=i/200
for i in range(n):
	y[i]=phi(x[i])
filename="heat initial condition"
plt.xlim(0,1)
plt.plot(x,y)
plt.subplots_adjust(left=0.15,right=0.98,top=0.92,bottom=0.05)
plt.suptitle(filename, fontsize = 12)
plt.plot(x,y,color='k')
plt.savefig(filename+'.jpg')
plt.clf()
for r in [1,2,0.5,10]:
	k=r/400
	for t in [1,2,10]:
		filename="exact solution(r="+str(r)+")(t="+str(t)+"k)"
		plt.xlim(0,1)
		plt.suptitle(filename, fontsize = 12)
		for i in range(n):
			y[i]=u(x[i],t*k)
		if (max(y)<=0.5+1e-9 and min(y)>=-1e-9 and max(y)>0.01):
			plt.ylim(0,0.5)
		plt.plot(x,y,color='k')
		plt.subplots_adjust(left=0.15,right=0.98,top=0.92,bottom=0.05)
		plt.savefig(filename+'.jpg')
		plt.clf()
for filename in ["CN(r=1)","CN(r=2)","BTCS(r=1)","BTCS(r=10)","FTCS(r=1)","FTCS(r=0.5)","RKM2(r=1)","RKM2(r=10)","GLRKM1(r=1)","GLRKM1(r=10)"]:
	fi=open(filename+'.txt',"r")
	x=list(map(float,fi.readline().split()))
	for k in [1,2,10]:
		plt.xlim(0,1)
		plt.suptitle(filename+'(t='+str(k)+'k)', fontsize = 12)
		y=list(map(float,fi.readline().split()))
		n=len(y)
		if (max(y)<=0.5+1e-9 and min(y)>=-1e-9 and max(y)>0.01):
			plt.ylim(0,0.5)
		plt.plot(x,y,color='k')
		plt.subplots_adjust(left=0.15,right=0.98,top=0.92,bottom=0.05)
		plt.savefig(filename+'(t='+str(k)+'k).jpg')
		plt.clf()
for filename in ["leapfrog(k=0.8h)","LF(k=0.8h)","LW(k=0.8h)","upwind(k=0.8h)","BW(k=0.8h)","leapfrog(k=h)","LW(k=h)"]:
	fi=open(filename+'.txt',"r")
	plt.suptitle(filename, fontsize = 12)
	x=list(map(float,fi.readline().split()))
	plt.xlim(x[0],x[-1])
	plt.ylim(-0.4,1)
	y=list(map(float,fi.readline().split()))
	n=len(y)
	plt.plot(x,y,color='k')
	for i in range(n):
		y[i]=f(x[i])
	plt.plot(x,y,color='r')
	plt.subplots_adjust(left=0.15,right=0.98,top=0.92,bottom=0.05)
	plt.savefig(filename+'.jpg')
	plt.clf()