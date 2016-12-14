// Linear algebra with C++
// Real matrices
//
// $Header: /home/bainvil/Modules/alp/RCS/rmatbase.C,v 2.2 1995/06/10 13:29:42 bainvil Exp bainvil $

#include <alp.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

#define LOOPL(u,v) for (register int u=1;u<=n1;u++) for (register int v=1;v<=n2;v++)
#define LOOPA(a,u,v) int n1=a.R(),n2=a.C();LOOPL(u,v)
#define min(a,b) ((a)<(b))?(a):(b)

// using namespace std;

void Matrix::Allocate()
    {
    x=new double [n1*n2];s1=n2;s2=1;
//    std::cerr<<"Matrix::Allocate : n1="<<n1<<" n2="<<n2<<" s1="<<s1<<" s2="<<s2<<std::endl;
    }

//
// Constructors
// ============
Matrix::Matrix(const int rows,const int columns) : n1(rows), n2(columns)
    { Allocate(); }
Matrix::Matrix(const int n) : n1(n), n2(n)
    { Allocate(); }
Matrix::Matrix(const Matrix &a) : n1(a.n1), n2(a.n2)
    { Allocate();Copy(a,*this); }
Matrix::Matrix(const Vector& u) : n1(u.N()), n2(1)
    { SetColumn(1,u); }

//
// Destructeur
// ===========
Matrix::~Matrix()
    {
//    std::cerr<<"Matrix::~Matrix"<<std::endl;
    delete [] x;
    }

//
// Diagonal, rows and columns operators
// ====================================
void GetDiag(const Matrix& a,Vector &u)
    { a.GetDiag(u); }
void GetRow(const Matrix& a,int i,Vector &u)
    { a.GetRow(i,u); }
void GetColumn(const Matrix& a,int j,Vector &u)
    { a.GetColumn(j,u); }
void SetDiag(Matrix& a,const Vector &u)
    { a.SetDiag(u); }
void SetRow(Matrix& a,int i,const Vector &u)
    { a.SetRow(i,u); }
void SetColumn(Matrix& a,int j,const Vector &u)
    { a.SetColumn(j,u); }

//
// Diagonal, rows and columns member functions
// ===========================================
void Matrix::GetDiag(Vector &u) const
    { Copy(ShareDiag(),u); }
void Matrix::GetRow(int i,Vector &u) const
    { Copy(ShareRow(i),u); }
void Matrix::GetColumn(int j,Vector &u) const
    { Copy(ShareColumn(j),u); }
void Matrix::SetDiag(const Vector &u)
    { Copy(u,ShareDiag()); }
void Matrix::SetRow(int i,const Vector &u)
    { Copy(u,ShareRow(i)); }
void Matrix::SetColumn(int j,const Vector &u)
    { Copy(u,ShareColumn(j)); }

void Matrix::Zero()
    {
    memset((char *)x,0,n1*n2*sizeof(x[0]));
    }
void Matrix::Id()
    {
    int i,p=min(n1,n2);
    memset((char *)x,0,n1*n2*sizeof(x[0]));
    for (i=1;i<=p;i++) (*this)(i,i)=1.0;
    }
void Matrix::Scaling(const double k)
    {
    int i,p=min(n1,n2);
    memset((char *)x,0,n1*n2*sizeof(x[0]));
    for (i=1;i<=p;i++) (*this)(i,i)=k;
    }
void Matrix::Rotation2D(const double a)
    {
    register double sn,cs;
    sn=sin(a);cs=cos(a);	// sincos ?
    (*this)(1,1)=(*this)(2,2)=cs;
    (*this)(1,2)= -((*this)(2,1)=sn);
    }
void Matrix::Rotation3D(const Vector &u,const double a)
    {
    double lambda,mu,x,y,z;
    u.Get(x,y,z);
    mu=1-cos(a);lambda=sin(a);
    (*this)(1,1)=1.0-mu*(y*y+z*z);
    (*this)(2,2)=1.0-mu*(x*x+z*z);
    (*this)(3,3)=1.0-mu*(x*x+y*y);
    (*this)(1,2)=-lambda*z+mu*x*y;
    (*this)(1,3)=lambda*y+mu*x*z;
    (*this)(2,1)=lambda*z+mu*x*y;
    (*this)(2,3)=-lambda*x+mu*y*z;
    (*this)(3,1)=-lambda*y+mu*x*z;
    (*this)(3,2)=lambda*x+mu*y*z;
    }
/* LD 
void Matrix::Rotation3D(const Vector &rot)
    {
    Vector v(3);
    double h;
    h= ::Norm(rot);  // Sans le ::, g++ 2.5.8 se trompe de Norm() :-)
    if (h<1.0e-12) { Id();return; } // Evite un gag !
    ::Mul(1.0/h,rot,v);		    // v = rot normalise
    Rotation3D(v,h);
    }
*/
void Matrix::Rotation3DFromUnitQuaternion(const Vector &q)
    {
    double a,b,c,d;
    q.Get(a,b,c,d);
    (*this)(1,1)=a*a+b*b-c*c-d*d;
    (*this)(1,2)=2.0*(b*c-a*d);
    (*this)(1,3)=2.0*(b*d+a*c);
    (*this)(2,1)=2.0*(b*c+a*d);
    (*this)(2,2)=a*a-b*b+c*c-d*d;
    (*this)(2,3)=2.0*(c*d-a*b);
    (*this)(3,1)=2.0*(b*d-a*c);
    (*this)(3,2)=2.0*(c*d+a*b);
    (*this)(3,3)=a*a-b*b-c*c+d*d;
    }

//
// Memory operators
// ================
Matrix & Matrix::operator = (const Matrix &a)
    { Copy(a,*this);return *this; }
void Copy(const Matrix &a,Matrix &b)
    {
    if (&a==&b) return;
    memcpy((char *)b.x,(char *)a.x,b.n1*b.n2*sizeof(double));b.s1=a.s1;b.s2=a.s2;
    }
void Swap(Matrix &a,Matrix &b)
    {
    register double k;
    LOOPA(a,i,j) { k=a(i,j);a(i,j)=b(i,j);b(i,j)=k; }
    }

//
// Bounds operators
// ================
void MinAbsIndex(const Matrix &a,int &i,int &j)
    { a.MinAbsIndex(i,j); }
void MaxAbsIndex(const Matrix &a,int &i,int &j)
    { a.MaxAbsIndex(i,j); }
void MinIndex(const Matrix &a,int &i,int &j)
    { a.MinIndex(i,j); }
void MaxIndex(const Matrix &a,int &i,int &j)
    { a.MaxIndex(i,j); }
double MinAbs(const Matrix &a)
    { return a.MinAbs(); }
double MaxAbs(const Matrix &a)
    { return a.MaxAbs(); }
double Min(const Matrix &a)
    { return a.Min(); }
double Max(const Matrix &a)
    { return a.Max(); }

//
// Bounds member functions
// =======================
void Matrix::MinAbsIndex(int &i,int &j) const
    {
    double s=fabs((*this)(1,1)),t;
    i=j=1;LOOPL(u,v) if ((t=fabs((*this)(u,v)))<s) { s=t;i=u;j=v; }
    }
void Matrix::MaxAbsIndex(int &i,int &j) const
    {
    double s=fabs((*this)(1,1)),t;
    i=j=1;LOOPL(u,v) if ((t=fabs((*this)(u,v)))>s) { s=t;i=u;j=v; }
    }
void Matrix::MinIndex(int &i,int &j) const
    {
    double s=(*this)(1,1),t;
    i=j=1;LOOPL(u,v) if ((t=(*this)(u,v))<s) { s=t;i=u;j=v; }
    }
void Matrix::MaxIndex(int &i,int &j) const
    {
    double s=(*this)(1,1),t;
    i=j=1;LOOPL(u,v) if ((t=(*this)(u,v))>s) { s=t;i=u;j=v; }
    }
double Matrix::MinAbs() const
    {
    double s=fabs((*this)(1,1)),t;
    LOOPL(u,v) if ((t=fabs((*this)(u,v)))<s) s=t;
    return s;
    }
double Matrix::MaxAbs() const
    {
    double s=fabs((*this)(1,1)),t;
    LOOPL(u,v) if ((t=fabs((*this)(u,v)))>s) s=t;
    return s;
    }
double Matrix::Min() const
    {
    double s=(*this)(1,1),t;
    LOOPL(u,v) if ((t=(*this)(u,v))<s) s=t;
    return s;
    }
double Matrix::Max() const
    {
    double s=(*this)(1,1),t;
    LOOPL(u,v) if ((t=(*this)(u,v))>s) s=t;
    return s;
    }

//
// Norms and dot product operators
// ===============================
double AbsSum(const Matrix &a)
    { return a.AbsSum(); }
double Sum(const Matrix &a)
    { return a.Sum(); }
double Norm2(const Matrix &a)
    { return a.Norm2(); }
double Norm(const Matrix &a)
    { return a.Norm(); }
double Dot(const Matrix &a,const Matrix &b)
    {
    register double s=0.0;
    LOOPA(a,i,j) s+=a(i,j)*b(i,j);
    return s;
    }
double Trace(const Matrix& a)
    { return a.Trace(); }

//
// Norms and dot product member functions
// ======================================
double Matrix::AbsSum() const
    { return ShareAsVector().AbsSum(); }
double Matrix::Sum() const
    { return ShareAsVector().Sum(); }
double Matrix::Norm2() const
    { return ShareAsVector().Norm2(); }
double Matrix::Norm() const
    { return ShareAsVector().Norm(); }
double Matrix::Trace() const
    { return ShareDiag().Sum(); }

//
// Arithmetic matrix-vector operators
// ==================================
void Image(const Vector &u,const Matrix &a,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { Image(Vector(u),a,v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=Dot(a.ShareRow(i),u);
    }
void Image(const Vector &u,
	   double alpha,const Matrix &a,
	   double beta,const Vector &b,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { Image(Vector(u),alpha,a,beta,b,v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=alpha*Dot(a.ShareRow(i),u)+beta*b(i);
    }
void TransImage(const Vector &u,const Matrix &a,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { TransImage(Vector(u),a,v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=Dot(a.ShareColumn(i),u);
    }
void TransImage(const Vector &u,
		double alpha,const Matrix &a,
		double beta,const Vector &b,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { TransImage(Vector(u),alpha,a,beta,b,v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=alpha*Dot(a.ShareColumn(i),u)+beta*b(i);
    }
void Update(const Matrix &a,double alpha,const Vector &x,const Vector &y,Matrix &b)
    { LOOPA(b,i,j) b(i,j)=a(i,j)+alpha*x(i)*y(j); }
double Quad(const Matrix &a,const Vector& x,const Vector& y)
    {
    int i;
    double s=0.0;
    for (i=1;i<=a.R();i++) s+=x(i)*Dot(a.ShareRow(i),y);
    return s;
    }
//
// Arithmetic matrix-vector member functions
// =========================================
void Matrix::Image(const Vector& u,Vector& v)
    { // Plus lent si &u == &v
    if (&u==&v) { Image(Vector(u),v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=Dot(ShareRow(i),u);
    }
void Matrix::Image(const Vector &u,double alpha,
	   double beta,const Vector &b,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { Image(Vector(u),alpha,beta,b,v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=alpha*Dot(ShareRow(i),u)+beta*b(i);
    }
void Matrix::TransImage(const Vector &u,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { TransImage(Vector(u),v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=Dot(ShareColumn(i),u);
    }
void Matrix::TransImage(const Vector &u,double alpha,
		double beta,const Vector &b,Vector &v)
    { // Plus lent si &u == &v
    if (&u==&v) { TransImage(Vector(u),alpha,beta,b,v);return; }
    for (register int i=1;i<=v.N();i++) v(i)=alpha*Dot(ShareColumn(i),u)+beta*b(i);
    }
void Matrix::Update(double alpha,const Vector &x,const Vector &y)
    { LOOPL(i,j) (*this)(i,j)+=alpha*x(i)*y(j); }
double Matrix::Quad(const Vector& x,const Vector& y)
    {
    int i;
    double s=0.0;
    for (i=1;i<=n1;i++) s+=x(i)*Dot(ShareRow(i),y);
    return s;
    }

//
// Transposition (takes the time of swapping 2 integers : fast)
// =============
void Transpose(const Matrix &a,Matrix &b)
    { Copy(a,b);b.Transpose(); }
void Matrix::Transpose()
    { int u=n1;n1=n2;n2=u;u=s1;s1=s2;s2=u; }

//
// Matrix additions operators
// ==========================
void Add(const Matrix &a,const Matrix &b,Matrix &c)
    { LOOPA(c,i,j) c(i,j)=a(i,j)+b(i,j); }
void Sub(const Matrix &a,const Matrix &b,Matrix &c)
    { LOOPA(c,i,j) c(i,j)=a(i,j)-b(i,j); }
void Add(double alpha,const Matrix &a,
	 double beta,const Matrix &b,Matrix &c)
    { LOOPA(c,i,j) c(i,j)=alpha*a(i,j)+beta*b(i,j); }
void TransAdd(const Matrix &a,const Matrix &b,Matrix &c)
    { // Plus lent si &b == &c
    if (&b==&c) { TransAdd(a,Matrix(b),c);return; }
    LOOPA(c,i,j) c(i,j)=a(i,j)+b(j,i);
    }
void TransAdd(double alpha,const Matrix &a,
	      double beta,const Matrix &b,Matrix &c)
    { // Plus lent si &b == &c
    if (&b==&c) { TransAdd(a,Matrix(b),c);return; }
    LOOPA(c,i,j) c(i,j)=alpha*a(i,j)+beta*b(j,i);
    }

//
// Matrix additions member functions
// =================================
Matrix& Matrix::operator += (const Matrix &a)
    { LOOPL(i,j) (*this)(i,j)+=a(i,j); return *this; }
Matrix& Matrix::operator -= (const Matrix &a)
    { LOOPL(i,j) (*this)(i,j)-=a(i,j); return *this; }
void Matrix::Add(const Matrix &a)
    { LOOPL(i,j) (*this)(i,j)+=a(i,j); }
void Matrix::Sub(const Matrix &a)
    { LOOPL(i,j) (*this)(i,j)-=a(i,j); }
void Matrix::Add(double k,const Matrix &a)
    { LOOPL(i,j) (*this)(i,j)+=k*a(i,j); }
void Matrix::Combine(double k,double l,const Matrix &a)
    { LOOPL(i,j) (*this)(i,j)=k*(*this)(i,j)+l*a(i,j); }

//
// Matrix-scalar products operators
// ================================
void Mul(double k,const Matrix &a,Matrix &b)
    { Copy(a,b);b.Mul(k); }
void Div(double k,const Matrix &a,Matrix &b)
    { Copy(a,b);b.Div(k); }
//
// Matrix-scalar products member functions
// =======================================
Matrix& Matrix::operator *= (double k)
    { Mul(k);return *this; }
Matrix& Matrix::operator /= (double k)
    { Div(k);return *this; }
void Matrix::Mul(double k)
    { LOOPL(i,j) (*this)(i,j)*=k; }
void Matrix::Div(double k)
    { double l=1.0/k;LOOPL(i,j) (*this)(i,j)*=l; }

//
// Matrix-matrix multiplication operators
// ======================================
void Mul(const Matrix &a,const Matrix &b,Matrix &c)
    {
    if (&a==&c || &b==&c) { Matrix d(c.n1,c.n2);Mul(a,b,d);Copy(d,c);return; }
    LOOPA(c,i,j) c(i,j)=Dot(a.ShareRow(i),b.ShareColumn(j));
    }
void Mul(double alpha,const Matrix &a,const Matrix &b,
	 double beta,const Matrix &c,Matrix &m)
    {
    if (&a==&m || &b==&m) { Matrix d(m.n1,m.n2);Mul(alpha,a,b,beta,c,d);Copy(d,m);return; }
    LOOPA(m,i,j) m(i,j)=beta*c(i,j)+alpha*Dot(a.ShareRow(i),b.ShareColumn(j));
    }
void RightTransMul(const Matrix &a,const Matrix &b,Matrix &c)
    {
    if (&a==&c || &b==&c) { Matrix d(c.n1,c.n2);RightTransMul(a,b,d);Copy(d,c);return; }
    LOOPA(c,i,j) c(i,j)=Dot(a.ShareRow(i),b.ShareRow(j));
    }
void LeftTransMul(const Matrix &a,const Matrix &b,Matrix &c)
    {
    if (&a==&c || &b==&c) { Matrix d(c.n1,c.n2);LeftTransMul(a,b,d);Copy(d,c);return; }
    LOOPA(c,i,j) c(i,j)=Dot(a.ShareColumn(i),b.ShareColumn(j));
    }
void RightTransMul(double alpha,const Matrix &a,const Matrix &b,
		   double beta,const Matrix &c,Matrix &m)
    {
    if (&a==&m || &b==&m)
	{ Matrix d(m.n1,m.n2);RightTransMul(alpha,a,b,beta,c,d);Copy(d,m);return; }
    LOOPA(m,i,j) m(i,j)=beta*c(i,j)+alpha*Dot(a.ShareRow(i),b.ShareRow(j));
    }
void LeftTransMul(double alpha,const Matrix &a,const Matrix &b,
		  double beta,const Matrix &c,Matrix &m)
    {
    if (&a==&m || &b==&m)
	{ Matrix d(m.n1,m.n2);LeftTransMul(alpha,a,b,beta,c,d);Copy(d,m);return; }
    LOOPA(m,i,j) m(i,j)=beta*c(i,j)+alpha*Dot(a.ShareColumn(i),b.ShareColumn(j));
    }

//
// Matrix-matrix multiplication member functions
// =============================================
void Matrix::RightMul(const Matrix &a)
    { ::Mul(*this,a,*this); }
void Matrix::LeftMul(const Matrix &a)
    { ::Mul(a,*this,*this); }
void Matrix::RightMul(double alpha,const Matrix &a,double beta,const Matrix &b)
    { ::Mul(alpha,*this,a,beta,b,*this); }
void Matrix::LeftMul(double alpha,const Matrix &a,double beta,const Matrix &b)
    { ::Mul(alpha,a,*this,beta,b,*this); }

//
// IO operators
// ============
std::ostream & operator << (std::ostream &o,const Matrix &a)
    {
    int i,j;
    for (i=1;i<=a.R();i++)
    for (j=1;j<=a.C();j++)
	{
	o << a(i,j);
	if (j==a.C()) o << std::endl; else o << "  ";
	}
    return o;
    }

std::ostream & TeX(std::ostream &o,const Matrix &a)
    {
    int i,j;
    
    // The output is a a.R rows by a.C columns of numbers
    // between ()
    o << "\\left( \\begin{array}{";
    for (i=0;i<a.C();i++) o << 'c';
    o << "} ";
    for (i=1;i<=a.R();i++)
    for (j=1;j<=a.C();j++)
	{
	o << a(i,j);
	if (j==a.C()) o << "\\\\ "; else o << " & ";
	}
    o << "\\end{array} \\right)";
    return o;
    }
