name=["Adams-Bashforth(1)","Adams-Bashforth(2)","Adams-Bashforth(3)","Adams-Bashforth(4)","Adams-Moulton(2)","Adams-Moulton(3)","Adams-Moulton(4)","Adams-Moulton(5)","BDF(1)","BDF(2)","BDF(3)","BDF(4)","classical RK","ESDIRK","Gauss-Legendre RK(1)","Gauss-Legendre RK(2)","Gauss-Legendre RK(3)","Fehlberg embedded RK","Dormand-Prince embedded RK"]
fo=open("bsearch.in","w")
def tran(s):
	t=""
	for c in s:
		if (c=='\''):
			c='"'
		t+=c
	return t
ban=["Adams-Bashforth(1)","BDF(1)","Adams-Bashforth(2)","BDF(2)","Adams-Moulton(2)"]
for s in name:
	if s not in ban:
		fo.write(tran(str([s,0,0,1]))+'\n')
for s in name:
	if s not in ban:
		fo.write(tran(str([s,1,0,1]))+'\n')
fo=open("error1.in","w")
for s in name:
	step=[500000,1000000]
	if (s in ["Gauss-Legendre RK(3)","Dormand-Prince embedded RK"]):
		step=[20000,40000]
	if (s in ["Adams-Moulton(5)"]):
		step=[100000,200000]
	if (s in ["Fehlberg embedded RK"]):
		step=[200000,400000]
	for x in step:
		fo.write(tran(str([s,0,x,1]))+'\n')
fo=open("error2.in","w")
for s in name:
	step=[50000,100000]
	if (s in ["Adams-Moulton(4)","BDF(3)","BDF(4)","classical RK","ESDIRK","Gauss-Legendre RK(1)","Gauss-Legendre RK(2)","Gauss-Legendre RK(3)","Fehlberg embedded RK","Dormand-Prince embedded RK"]):
		step=[5000,10000]
	if (s in ["Gauss-Legendre RK(3)","Dormand-Prince embedded RK"]):
		step=[2000,4000]
	if (s in ["Adams-Moulton(5)"]):
		step=[1000,2000]
	if (s in ["Fehlberg embedded RK"]):
		step=[2000,4000]
	for x in step:
		fo.write(tran(str([s,1,x,1]))+'\n')