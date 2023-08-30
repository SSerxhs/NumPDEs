from math import sin,cos,exp
import matplotlib.pyplot as plt
for file in ["Euler","classical RK","Dormand-Prince embedded RK","tmp"]:
    fi=open(file+".in","r")
    n=int(fi.readline())
    x,y=[0]*n,[0]*n
    for i in range(n):
        x[i],y[i]=map(float,fi.readline().split())
    plt.plot(x,y)
    plt.savefig(file+'.jpg')
    plt.clf()