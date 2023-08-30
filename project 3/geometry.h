namespace geometry//不要用 int！
{
#define tmpl template<typename T>
	typedef long long ll;
	typedef double db;
	const db eps=1e-8;
#define all(x) (x).begin(),(x).end()
	inline int sgn(const ll &x)
	{
		if (x<0) return -1;
		return x>0;
	}
	inline int sgn(const db &x)
	{
		if (fabs(x)<eps) return 0;
		return x>0?1:-1;
	}
	tmpl struct point//* 为叉乘，& 为点乘，只允许使用 double 和 ll
	{
		T x,y;
		point() {}
		point(T a,T b) :x(a),y(b) {}
		operator point<ll>() const { return point<ll>(x,y); }
		operator point<db>() const { return point<db>(x,y); }
		point<T> operator+(const point<T> &o) const { return point(x+o.x,y+o.y); }
		point<T> operator-(const point<T> &o) const { return point(x-o.x,y-o.y); }
		point<T> operator*(const T &k) const { return point(x*k,y*k); }
		point<T> operator/(const T &k) const { return point(x/k,y/k); }
		T operator*(const point<T> &o) const { return x*o.y-y*o.x; }
		T operator&(const point<T> &o) const { return x*o.x+y*o.y; }
		void operator+=(const point<T> &o) { x+=o.x; y+=o.y; }
		void operator-=(const point<T> &o) { x+=o.x; y+=o.y; }
		void operator*=(const T &k) { x*=k; y*=k; }
		void operator/=(const T &k) { x/=k; y/=k; }
		bool operator==(const point<T> &o) const { return x==o.x&&y==o.y; }
		bool operator!=(const point<T> &o) const { return x!=o.x||y!=o.y; }
		db len() const { return sqrt(len2()); }//模长
		T len2() const { return x*x+y*y; }
	};
	const point<db> npos=point<db>(514e194,9810e191),apos=point<db>(145e174,999e180);
	const int DS[4]={1,2,4,3};
	tmpl int quad(const point<T> &o)//坐标轴归右上象限，返回值 [1,4]
	{
		return DS[(sgn(o.y)<0)*2+(sgn(o.x)<0)];
	}
	tmpl bool angle_cmp(const point<T> &a,const point<T> &b)
	{
		int c=quad(a),d=quad(b);
		if (c!=d) return c<d;
		return a*b>0;
	}
	tmpl db dis(const point<T> &a,const point<T> &b) { return (a-b).len(); }
	tmpl T dis2(const point<T> &a,const point<T> &b) { return (a-b).len2(); }
	tmpl point<T> operator*(const T &k,const point<T> &o) { return point<T>(k*o.x,k*o.y); }
	tmpl bool operator<(const point<T> &a,const point<T> &b)
	{
		int s=sgn(a*b);
		return s>0||s==0&&sgn(a.len2()-b.len2())<0;
	}
	tmpl istream &operator>>(istream &cin,point<T> &o) { return cin>>o.x>>o.y; }
	tmpl ostream &operator<<(ostream &cout,const point<T> &o)
	{
		if ((point<db>)o==apos) return cout<<"all position";
		if ((point<db>)o==npos) return cout<<"no position";
		return cout<<'('<<o.x<<','<<o.y<<')';
	}
	tmpl struct line
	{
		point<T> o,d;
		line() {}
		line(const point<T> &a,const point<T> &b,int twopoint);
		bool operator!=(const line<T> &m) { return !(*this==m); }
	};
	template<> line<ll>::line(const point<ll> &a,const point<ll> &b,int twopoint)
	{
		o=a;
		d=twopoint?b-a:b;
		ll tmp=gcd(d.x,d.y);
		assert(tmp);
		if (d.x<0||d.x==0&&d.y<0) tmp=-tmp;
		d.x/=tmp; d.y/=tmp;
	}
	template<> line<db>::line(const point<db> &a,const point<db> &b,int twopoint)
	{
		o=a;
		d=twopoint?b-a:b;
		int s=sgn(d.x);
		if (s<0||!s&&d.y<0) d.x=-d.x,d.y=-d.y;
	}
	tmpl line<T> rotate_90(const line<T> &m) { return line(m.o,point(m.d.y,-m.d.x),0); }
	tmpl line<db> rotate(const line<T> &m,db angle)
	{
		return {(point<db>)m.o,{m.d.x*cos(angle)-m.d.y*sin(angle),m.d.x*sin(angle)+m.d.y*cos(angle)},0};
	}
	tmpl db get_angle(const line<T> &m,const line<T> &n) { return asin((m.d*n.d)/(m.d.len()*n.d.len())); }
	tmpl bool operator<(const line<T> &m,const line<T> &n)
	{
		int s=sgn(m.d*n.d);
		return s?s>0:m.d*m.o<n.d*n.o;
	}
	bool operator==(const line<ll> &m,const line<ll> &n) { return m.d==n.d&&(m.o-n.o)*m.d==0; }
	bool operator==(const line<db> &m,const line<db> &n) { return fabs(m.d*n.d)<eps&&fabs((n.o-m.o)*m.d)<eps; }
	tmpl ostream &operator<<(ostream &cout,const line<T> &o) { return cout<<'('<<o.d.x<<" k + "<<o.o.x<<" , "<<o.d.y<<" k + "<<o.o.y<<")"; }
	tmpl point<db> intersect(const line<T> &m,const line<T> &n)
	{
		if (!sgn(m.d*n.d))
		{
			if (!sgn(m.d*(n.o-m.o))) return apos;
			return npos;
		}
		return (point<db>)m.o+(n.o-m.o)*n.d/(db)(m.d*n.d)*(point<db>)m.d;
	}
	tmpl db dis(const line<T> &m,const point<T> &o) { return abs(m.d*(o-m.o)/m.d.len()); }
	tmpl db dis(const point<T> &o,const line<T> &m) { return abs(m.d*(o-m.o)/m.d.len()); }
	struct circle
	{
		point<db> o;
		db r;
		circle() {}
		circle(const point<db> &O,const db &R=0) :o(point<db>((db)O.x,(db)O.y)),r(R) {}//圆心半径构造
		circle(const point<db> &a,const point<db> &b)//直径构造
		{
			o=(a+b)*0.5;
			r=dis(b,o);
		}
		circle(const point<db> &a,const point<db> &b,const point<db> &c)//三点构造外接圆（非最小圆）
		{
			auto A=(b+c)*0.5,B=(a+c)*0.5;
			o=intersect(rotate_90(line(A,c,1)),rotate_90(line(B,c,1)));
			r=dis(o,c);
		}
		circle(vector<point<db>> a)
		{
			int n=a.size(),i,j,k;
			mt19937 rnd(75643);
			shuffle(all(a),rnd);
			*this=circle(a[0]);
			for (i=1; i<n; i++) if (!cover(a[i]))
			{
				*this=circle(a[i]);
				for (j=0; j<i; j++) if (!cover(a[j]))
				{
					*this=circle(a[i],a[j]);
					for (k=0; k<j; k++) if (!cover(a[k])) *this=circle(a[i],a[j],a[k]);
				}
			}
		}
		circle(const vector<point<ll>> &b)
		{
			vector<point<db>> a(b.size());
			int n=a.size(),i,j,k;
			for (i=0; i<a.size(); i++) a[i]=(point<db>)b[i];
			*this=circle(a);
		}
		tmpl bool cover(const point<T> &a) { return sgn(dis((point<db>)a,o)-r)<=0; }
	};
	tmpl struct segment
	{
		point<T> a,b;
		segment() {}
		segment(point<T> o,point<T> p)
		{
			int s=sgn(o.x-p.x);
			if (s>0||!s&&o.y>p.y) swap(o,p);
			a=o; b=p;
		}
	};
	tmpl bool intersect(const segment<T> &m,const segment<T> &n)
	{
		auto a=n.b-n.a,b=m.b-m.a;
		auto d=n.a-m.a;
		if (sgn(n.b.x-m.a.x)<0||sgn(m.b.x-n.a.x)<0) return 0;
		if (sgn(max(n.a.y,n.b.y)-min(m.a.y,m.b.y))<0||sgn(max(m.a.y,m.b.y)-min(n.a.y,n.b.y))<0) return 0;
		return sgn(b*d)*sgn((n.b-m.a)*b)>=0&&sgn(a*d)*sgn((m.b-n.a)*a)<=0;
	}
	tmpl struct convex
	{
		vector<point<T>> p;
		convex(vector<point<T>> a);
		db peri()//周长
		{
			int i,n=p.size();
			db C=(p[n-1]-p[0]).len();
			for (i=1; i<n; i++) C+=(p[i-1]-p[i]).len();
			return C;
		}
		db area() { return area2()*0.5; }//面积
		T area2()//两倍面积
		{
			int i,n=p.size();
			T S=p[n-1]*p[0];
			for (i=1; i<n; i++) S+=p[i-1]*p[i];
			return abs(S);
		}
		db diam() { return sqrt(diam2()); }
		T diam2()//直径平方
		{
			T r=0;
			int n=p.size(),i,j;
			if (n<=2)
			{
				for (i=0; i<n; i++) for (j=i+1; j<n; j++) r=max(r,dis2(p[i],p[j]));
				return r;
			}
			p.push_back(p[0]);
			for (i=0,j=1; i<n; i++)
			{
				while ((p[i+1]-p[i])*(p[j]-p[i])<=(p[i+1]-p[i])*(p[j+1]-p[i])) if (++j==n) j=0;
				r=max({r,dis2(p[i],p[j]),dis2(p[i+1],p[j])});
			}
			p.pop_back();
			return r;
		}
		bool cover(const point<T> &o) const//点是否在凸包内
		{
			if (o.x<p[0].x||o.x==p[0].x&&o.y<p[0].y) return 0;
			if (o==p[0]) return 1;
			if (p.size()==1) return 0;
			ll tmp=(o-p[0])*(p.back()-p[0]);
			if (tmp==0) return dis2(o,p[0])<=dis2(p.back(),p[0]);
			if (tmp<0||p.size()==2) return 0;
			int x=upper_bound(1+all(p),o,[&](const point<T> &a,const point<T> &b) { return (a-p[0])*(b-p[0])>0; })-p.begin()-1;
			return (o-p[x])*(p[x+1]-p[x])<=0;
		}
		convex<T> operator+(const convex<T> &A) const
		{
			int n=p.size(),m=A.p.size(),i,j;
			vector<point<T>> c;
			if (min(n,m)<=2)
			{
				c.reserve(n*m);
				for (i=0; i<n; i++) for (j=0; j<m; j++) c.push_back(p[i]+A.p[j]);
				return convex<T>(c);
			}
			point<T> a[n],b[m];
			for (i=0; i+1<n; i++) a[i]=p[i+1]-p[i];
			a[n-1]=p[0]-p[n-1];
			for (i=0; i+1<m; i++) b[i]=A.p[i+1]-A.p[i];
			b[m-1]=A.p[0]-A.p[m-1];
			c.reserve(n+m);
			c.push_back(p[0]+A.p[0]);
			for (i=j=0; i<n&&j<m;) c.push_back(c.back()+(a[i]*b[j]>0?a[i++]:b[j++]));
			while (i<n-1) c.push_back(c.back()+a[i++]);
			while (j<m-1) c.push_back(c.back()+b[j++]);
			return convex<T>(c);
		}
		void operator+=(const convex &a) { *this=*this+a; }
	};
	tmpl convex<T>::convex(vector<point<T>> a)
	{
		int n=a.size(),i;
		if (!n) return;
		p=a;
		for (i=1; i<n; i++) if (p[i].x<p[0].x||p[i].x==p[0].x&&p[i].y<p[0].y) swap(p[0],p[i]);
		a.resize(0); a.reserve(n);
		for (i=1; i<n; i++) if (p[i]!=p[0]) a.push_back(p[i]-p[0]);
		sort(all(a));
		for (i=0; i<a.size(); i++) a[i]+=p[0];
		point<T> *st=p.data()-1;
		int tp=1;
		for (auto &v:a)
		{
			while (tp>1&&sgn((st[tp]-st[tp-1])*(v-st[tp-1]))<=0) --tp;
			st[++tp]=v;
		}
		p.resize(tp);
	}
	template<> bool convex<db>::cover(const point<db> &o) const//点是否在凸包内
	{
		if (o.x<p[0].x||o.x==p[0].x&&o.y<p[0].y) return 0;
		if (o==p[0]) return 1;
		if (p.size()==1) return 0;
		ll tmp=(o-p[0])*(p.back()-p[0]);
		if (tmp==0) return dis2(o,p[0])<=dis2(p.back(),p[0]);
		if (tmp<0||p.size()==2) return 0;
		int x=upper_bound(1+all(p),o,[&](const point<db> &a,const point<db> &b) { return (a-p[0])*(b-p[0])>eps; })-p.begin()-1;
		return (o-p[x])*(p[x+1]-p[x])<=0;
	}
	tmpl struct half_plane//默认左侧
	{
		point<T> o,d;
		operator half_plane<ll>() const { return {(point<ll>)o,(point<ll>)d,0}; }
		operator half_plane<db>() const { return {(point<db>)o,(point<db>)d,0}; }
		half_plane() {}
		half_plane(const point<T> &a,const point<T> &b,bool twopoint)
		{
			o=a;
			d=twopoint?b-a:b;
		}
		bool operator<(const half_plane<T> &a) const
		{
			int p=quad(d),q=quad(a.d);
			if (p!=q) return p<q;
			p=sgn(d*a.d);
			if (p) return p>0;
			return sgn(d*(a.o-o))>0;
		}
	};
	tmpl ostream &operator<<(ostream &cout,half_plane<T> &m) { return cout<<m.o<<" | "<<m.d; }
	tmpl point<db> intersect(const half_plane<T> &m,const half_plane<T> &n)
	{
		if (!sgn(m.d*n.d))
		{
			if (!sgn(m.d*(n.o-m.o))) return apos;
			return npos;
		}
		return (point<db>)m.o+(n.o-m.o)*n.d/(db)(m.d*n.d)*(point<db>)m.d;
	}
	const db inf=1e9;
	tmpl convex<db> intersect(vector<half_plane<T>> a)
	{
		T I=inf;
		a.push_back({{-I,-I},{I,-I},1});
		a.push_back({{I,-I},{I,I},1});
		a.push_back({{I,I},{-I,I},1});
		a.push_back({{-I,I},{-I,-I},1});
		sort(all(a));
		int n=a.size(),i,h=0,t=-1;
		half_plane<db> q[n];
		point<db> p[n];
		vector<point<db>> r;
		for (i=0; i<n; i++) if (i==n-1||sgn(a[i].d*a[i+1].d))
		{
			auto x=(half_plane<db>)a[i];
			while (h<t&&sgn((p[t-1]-x.o)*x.d)>=0) --t;
			while (h<t&&sgn((p[h]-x.o)*x.d)>=0) ++h;
			q[++t]=x;
			if (h<t) p[t-1]=intersect(q[t-1],q[t]);
		}
		while (h<t&&sgn((p[t-1]-q[h].o)*q[h].d)>=0) --t;
		if (h==t) return convex<db>(vector<point<db>>(0));
		p[t]=intersect(q[h],q[t]);
		return convex<db>(vector<point<db>>(p+h,p+t+1));
	}
	tmpl db dis(const point<db> &o,const segment<T> &l)
	{
		if ((l.b-l.a&o-l.a)<0||(l.a-l.b&o-l.b)<0) return min(dis(o,l.a),dis(o,l.b));
		return dis(o,line(l.a,l.b,1));
	}
	tmpl db dis(const segment<T> &l,const point<db> &o)
	{
		if ((l.b-l.a&o-l.a)<0||(l.a-l.b&o-l.b)<0) return min(dis(o,l.a),dis(o,l.b));
		return dis(o,line(l.a,l.b,1));
	}
#undef tmpl
}