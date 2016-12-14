// Linear algebra with C++
// Real vectors
//
// EB 93 94 95

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>

#include <alp.h>

void
Vector::Allocate()
    { if (n>0) { x=new double[n];s=1;shared=0; } else shared=1; }

//
// Constructors
// ============
Vector::Vector(int d) : n(d)
    { Allocate(); }
Vector::Vector(double x1,double x2) : n(2)
    { Allocate();Set(x1,x2); }
Vector::Vector(double x1,double x2,double x3) : n(3)
    { Allocate();Set(x1,x2,x3); }
Vector::Vector(double x1,double x2,double x3,double x4) : n(4)
    { Allocate();Set(x1,x2,x3,x4); }
Vector::Vector(const Vector &u) : n(u.n)
    { Allocate();Copy(u,*this); }

//
// Destructor
// ==========
Vector::~Vector()
    { if (!shared) delete [] x; }

//
// Initializations and structure access
// ====================================
void
Vector::Zero()
    {
    int i,ii=0;
    for (i=n;i>0;i--) { x[ii]=0.0;ii+=s; }
    }
void
Vector::Set(double x1,double x2)
    { (*this)(1)=x1;(*this)(2)=x2; }
void
Vector::Set(double x1,double x2,double x3)
    { (*this)(1)=x1;(*this)(2)=x2;(*this)(3)=x3; }
void
Vector::Set(double x1,double x2,double x3,double x4)
    { (*this)(1)=x1;(*this)(2)=x2;(*this)(3)=x3;(*this)(4)=x4; }
void
Vector::Get(double &x1,double &x2) const
    { x1=(*this)(1);x2=(*this)(2); }
void
Vector::Get(double &x1,double &x2,double &x3) const
    { x1=(*this)(1);x2=(*this)(2);x3=(*this)(3); }
void
Vector::Get(double &x1,double &x2,double &x3,double &x4) const
    { x1=(*this)(1);x2=(*this)(2);x3=(*this)(3);x4=(*this)(4); }

SharedVector
Vector::SubVector(int first,int size)
   { return SharedVector(x+(first-1)*s,size,s); }

//const 
//SharedVector
//Vector::SubVector(int first,int size) 
// const
//  {  return SharedVector(x+(first-1)*s,size,s); }

//
// Bounds operators
// ================
int
MinAbsIndex(const Vector& a)
    { return a.MinAbsIndex(); }
double
MinAbs(const Vector& a)
    { return a.MinAbs(); }
void
MinAbsIndexAndValue(const Vector& a,int& index,double& value)
    { a.MinAbsIndexAndValue(index,value); }
int
MaxAbsIndex(const Vector& a)
    { return a.MaxAbsIndex(); }
double
MaxAbs(const Vector& a)
    { return a.MaxAbs(); }
void
MaxAbsIndexAndValue(const Vector& a,int& index,double& value)
    { a.MaxAbsIndexAndValue(index,value); }
int
MinIndex(const Vector& a)
    { return a.MinIndex(); }
double
Min(const Vector& a)
    { return a.Min(); }
void
MinIndexAndValue(const Vector& a,int& index,double& value)
    { a.MinIndexAndValue(index,value); }
int
MaxIndex(const Vector& a)
    { return a.MaxIndex(); }
double
Max(const Vector& a)
    { return a.Max(); }
void
MaxIndexAndValue(const Vector& a,int& index,double& value)
    { a.MaxIndexAndValue(index,value); }

//
// Bounds member functions
// =======================
int
Vector::MinAbsIndex() const
    { int i;double v;MinAbsIndexAndValue(i,v);return i; }
double
Vector::MinAbs() const
    { int i;double v;MinAbsIndexAndValue(i,v);return v; }
void
Vector::MinAbsIndexAndValue(int& index,double& value) const
    {
    int i,ii=s,bi=0;
    double h,bh=fabs(x[0]);
    for (i=1;i<n;i++,ii+=s)
	if ((h=fabs(x[ii]))<bh) { bh=h;bi=i; }
    index=bi+1;value=bh;
    }
int
Vector::MaxAbsIndex() const
    { int i;double v;MaxAbsIndexAndValue(i,v);return i; }
double
Vector::MaxAbs() const
    { int i;double v;MaxAbsIndexAndValue(i,v);return v; }
void
Vector::MaxAbsIndexAndValue(int& index,double& value) const
    {
    int i,ii=s,bi=0;
    double h,bh=fabs(x[0]);
    for (i=1;i<n;i++,ii+=s)
	if ((h=fabs(x[ii]))>bh) { bh=h;bi=i; }
    index=bi+1;value=bh;
    }
int
Vector::MinIndex() const
    { int i;double v;MinIndexAndValue(i,v);return i; }
double
Vector::Min() const
    { int i;double v;MinIndexAndValue(i,v);return v; }
void
Vector::MinIndexAndValue(int& index,double& value) const
    {
    int i,ii=s,bi=0;
    double h,bh=x[0];
    for (i=1;i<n;i++,ii+=s)
	if ((h=x[ii])<bh) { bh=h;bi=i; }
    index=bi+1;value=bh;
    }
int
Vector::MaxIndex() const
    { int i;double v;MaxIndexAndValue(i,v);return i; }
double
Vector::Max() const
    { int i;double v;MaxIndexAndValue(i,v);return v; }
void
Vector::MaxIndexAndValue(int& index,double& value) const
    {
    int i,ii=s,bi=0;
    double h,bh=x[0];
    for (i=1;i<n;i++,ii+=s)
	if ((h=x[ii])>bh) { bh=h;bi=i; }
    index=bi+1;value=bh;
    }

//
// Memory operators
// ================
Vector &
Vector::operator = (const Vector &p)
    {
    Copy(p,*this);
    return *this;
    }
void
Copy(const Vector &a,Vector &b)
    {
    int i,i1=0,i2=0;
    for (i=b.n;i>0;i--) { b.x[i2]=a.x[i1];i1+=a.s;i2+=b.s; }
    }
void
Swap(Vector &a,Vector &b)
    {
    int i,i1=0,i2=0;
    double h;
    for (i=b.n;i>0;i--)
	{ h=b.x[i2];b.x[i2]=a.x[i1];a.x[i1]=h;i1+=a.s;i2+=b.s; }
    }

//
// Norms, distances and dot product operators
// ==========================================
double
AbsSum(const Vector &a)
    { return a.AbsSum(); }
double
Sum(const Vector &a)
    { return a.Sum(); }
double
Norm2(const Vector &a)
    { return a.Norm2(); }
double
Norm(const Vector &a)
    { return a.Norm(); }
double
Dist2(const Vector &a,const Vector &b)
    {
    int i,i1=0,i2=0;
    double h,s=0.0;
    for (i=a.n;i>0;i--)
	{ h=a.x[i1]-b.x[i2];s+=h*h;i1+=a.s;i2+=b.s; }
    return s;
    }
double
Dist(const Vector &a,const Vector &b)
    { return sqrt(Dist2(a,b)); }
double
Dot(const Vector &a,const Vector &b)
    {
    int i,i1=0,i2=0;
    double s=0.0;
    for (i=a.n;i>0;i--)
	{ s+=a.x[i1]*b.x[i2];i1+=a.s;i2+=b.s; }
    return s;
    }
void
Normalize(const Vector &u,Vector &v)
    { Div(Norm(u),u,v); }

//
// Norms member functions
// ======================
double
Vector::AbsSum() const
    {
    int i,ii=0;
    double h=0.0;
    for (i=n;i>0;i--) { h+=fabs(x[ii]);ii+=s; }
    return h;
    }
double
Vector::Sum() const
    {
    int i,ii=0;
    double h=0.0;
    for (i=n;i>0;i--) { h+=x[ii];ii+=s; }
    return h;
    }
double
Vector::Norm2() const
    {
    int i,ii=0;
    double h=0.0;
    for (i=n;i>0;i--) { h+=x[ii]*x[ii];ii+=s; }
    return h;
    }
double
Vector::Norm() const
    { return sqrt(Norm2()); }
void
Vector::Normalize()
    { Div(Norm()); }

//
// Vector space operators
// ======================
void
Mul(const double k,const Vector &u,Vector &v)
    {
    int i,i1=0,i2=0;
    for (i=v.n;i>0;i--)
	{ v.x[i2]=k*u.x[i1];i1+=u.s;i2+=v.s; }
    }
void
Div(const double k,const Vector &u,Vector &v)
    { Mul(1.0/k,u,v); }
void
Add(const Vector &u,const Vector &v,Vector &w)
    {
    int i,i1=0,i2=0,i3=0;
    for (i=w.n;i>0;i--)
	{ w.x[i3]=u.x[i1]+v.x[i2];i1+=u.s;i2+=v.s;i3+=w.s; }
    }
void
Sub(const Vector &u,const Vector &v,Vector &w)
    {
    int i,i1=0,i2=0,i3=0;
    for (i=w.n;i>0;i--)
	{ w.x[i3]=u.x[i1]-v.x[i2];i1+=u.s;i2+=v.s;i3+=w.s; }
    }
void
Add(const Vector &u,double k,const Vector &v,Vector &w)
    {
    int i,i1=0,i2=0,i3=0;
    for (i=w.n;i>0;i--)
	{ w.x[i3]=u.x[i1]+k*v.x[i2];i1+=u.s;i2+=v.s;i3+=w.s; }
    }
void
Combine(double k1,const Vector &u1,
	double k2,const Vector &u2,Vector &v)
    {
    int i,i1=0,i2=0,i3=0;
    for (i=v.n;i>0;i--)
	{ v.x[i3]=k1*u1.x[i1]+k2*u2.x[i2];i1+=u1.s;i2+=u2.s;i3+=v.s; }
    }
void
Combine(double k1,const Vector &u1,
	double k2,const Vector &u2,
	double k3,const Vector &u3,Vector &v)
    {
    int i,i1=0,i2=0,i3=0,i4=0;
    for (i=v.n;i>0;i--)
	{
	v.x[i4]=k1*u1.x[i1]+k2*u2.x[i2]+k3*u3.x[i3];
	i1+=u1.s;i2+=u2.s;i3+=u3.s;i4+=v.s;
	}
    }
void
Middle(const Vector &u,const Vector &v,Vector &w)
    {
    int i,i1=0,i2=0,i3=0;
    for (i=w.n;i>0;i--)
	{ w.x[i3]=0.5*(u.x[i1]+v.x[i2]);i1+=u.s;i2+=v.s;i3+=w.s; }
    }
void
Barycenter(double k1,const Vector& a1,
	   double k2,const Vector& a2,Vector& a)
    {
    double h=1.0/(k1+k2);
    Combine(k1*h,a1,k2*h,a2,a);
    }
void
Barycenter(double k1,const Vector& a1,
	   double k2,const Vector& a2,
	   double k3,const Vector& a3,Vector& a)
    {
    double h=1.0/(k1+k2+k3);
    Combine(k1*h,a1,k2*h,a2,k3*h,a3,a);
    }

//
// Vector space member functions
// =============================
Vector &
Vector::operator += (const Vector &u)
    {
    int i,i1=0,i2=0;
    for (i=n;i>0;i--) { x[i1]+=u.x[i2];i1+=s;i2+=u.s; }
    return *this;
    }
Vector &
Vector::operator -= (const Vector &u)
    {
    int i,i1=0,i2=0;
    for (i=n;i>0;i--) { x[i1]-=u.x[i2];i1+=s;i2+=u.s; }
    return *this;
    }
Vector &
Vector::operator *= (double k)
    {
    int i,i1=0;
    for (i=n;i>0;i--) { x[i1]*=k;i1+=s; }
    return *this;
    }
Vector &
Vector::operator /= (double k)
    {
    int i,i1=0;
    double l=1.0/k;
    for (i=n;i>0;i--) { x[i1]*=l;i1+=s; }
    return *this;
    }
void
Vector::Add(const Vector& u)
    { (*this)+=u; }
void
Vector::Sub(const Vector& u)
    { (*this)-=u; }
void
Vector::Mul(double k)
    { (*this)*=k; }
void
Vector::Div(double k)
    { (*this)/=k; }
void
Vector::Add(double k,const Vector& u)
    {
    int i,i1=0,i2=0;
    for (i=n;i>0;i--) { x[i1]+=k*u.x[i2];i1+=s;i2+=u.s; }
    }
void
Vector::Combine(double k,double l,const Vector& u)
    {
    int i,i1=0,i2=0;
    for (i=n;i>0;i--) { x[i1]=k*x[i1]+l*u.x[i2];i1+=s;i2+=u.s; }
    }

//
// Vector order (added Aug 95)
// ============
int
IsSup(const Vector& p,const Vector& ref)
    {
    int is_sup=1,i,c;
    for (i=p.N();i>0;i--)
	{
	c=AlpCompare(p(i),ref(i));
	if (c<0) return -1;	// We can stop now, it is not superior
	if (c==0) is_sup=0;
	}
    return is_sup;
    }

int
IsInf(const Vector& p,const Vector& ref)
    {
    int is_inf=1,i,c;
    for (i=p.N();i>0;i--)
	{
	c=AlpCompare(p(i),ref(i));
	if (c>0) return -1;	// We can stop now, it is not inferior
	if (c==0) is_inf=0;
	}
    return is_inf;
    }

//
// Dimension specific operators
// ============================
double
Det(const Vector &u,const Vector &v)
    { return u(1)*v(2)-u(2)*v(1); }
double
Det(const Vector &u,const Vector &v,const Vector &w)
    {
    return (u(1)*v(2)*w(3)+
	    u(2)*v(3)*w(1)+
	    u(3)*v(1)*w(2)-
	    u(1)*v(3)*w(2)-
	    u(2)*v(1)*w(3)-
	    u(3)*v(2)*w(1));
    }
void
Cross(const Vector &u,const Vector &v,Vector &w)
    {
    Vector h(3);
    h(1)=u(2)*v(3)-u(3)*v(2);
    h(2)=u(3)*v(1)-u(1)*v(3);
    h(3)=u(1)*v(2)-u(2)*v(1);
    Copy(h,w);
    }

//
// IO operators
// ============
/*
std::ostream &
operator << (std::ostream &o,const Vector &u)
    {
    int i,n=u.N();
    o << '<';
    for (i=1;i<=n;i++)
	{
	o << u(i);
	if (i<n) o << ','; else o << '>';
	}
    return o;
    }

std::ostream &
TeX(std::ostream &o,const Vector &u)
    {
    int i,n=u.N();
    
    // The output is a vertical array of reals between ()
    o << "\\left( \\begin{array}{c} ";
    for (i=1;i<=n;i++)
	{
	o << u(i);
	if (i<n) o << "\\\\ ";
	else o << "\\end{array} \\right)";
	}
    return o;
    }
*/
