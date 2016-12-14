// Isothetic boxes
//
// EB Aug 95

#include <alp.h>

//
// Construction - destruction
// ==========================
//
IBox::IBox(const Vector& a,const Vector& b) : n(a.N()), m(a.N()), M(a.N())
    { Set(a,b); }

IBox::IBox(int d) : n(d), m(d), M(d), is_empty(1)
    {
    }

//
// Initialization
// ==============

void
IBox::Empty()
    { is_empty=1; }

void
IBox::Unit()
    {
    int i;
    is_empty=0;
    for (i=n;i>0;i--) { m(i)=-1;M(i)=1; }
    }

void
IBox::Space()
    {
    int i;
    is_empty=0;
    for (i=n;i>0;i--) { m(i)=-AlpInfty;M(i)=AlpInfty; }
    }

void
IBox::Set(int i,double a,double b)
    {
    is_empty=0;
    if (a<b) { m(i)=a;M(i)=b; } else { m(i)=b;M(i)=a; }
    }

void
IBox::Set(const Vector& a,const Vector& b)
    {
    int i;
    is_empty=0;
    for (i=n;i>0;i--) if (a(i)<=b(i)) { m(i)=a(i);M(i)=b(i); }
    else { M(i)=a(i);m(i)=b(i); }
    }

//
// Queries
// =======
int
IBox::Compare(int i,double x) const
    { return AlpCompare(x,m(i))+AlpCompare(x,M(i)); }

int
Compare(int i,const IBox& a,const IBox& b)
    { return AlpCompare(a.M(i),b.m(i))-AlpCompare(b.M(i),a.m(i)); }

int
IBox::Contains(const Vector& p) const
    {
    int cm,cM;
    cm=IsSup(p,m);
    if (cm<0) return -1;	// one of the p(i) < m(i) -> outside
    cM=IsInf(p,M);
    if (cM<0) return -1;	// one of the p(i) > M(i) -> outside
    if (cm==0 || cM==0) return 0; // Boundary
    return 1;			  // In other cases, -> inside
    }

int
IBox::Contains(const IBox& a) const
    {
    int cont=0,i,cm,cM;
    for (i=n;i>0;i--)
	{
	cm=AlpCompare(m(i),a.m(i));
	if (cm>0) return -1;	// a is not included in this
	if (cm<0) cont=1;	// Inclusion is strict
	cM=AlpCompare(M(i),a.M(i));
	if (cM<0) return -1;	// a is not included in this
	if (cM>0) cont=1;	// Inclusion is strict
	}
    return cont;
    }

int
IBox::Intersects(const IBox& a) const
    {
    int inter=1,i,c;
    for (i=n;i>0;i--)
	{
	c=AlpCompare(M(i),a.m(i));
	if (c<0) return -1;	// No intersection
	if (c==0) inter=0;	// No interior intersection is possible now
	c=AlpCompare(a.M(i),m(i));
	if (c<0) return -1;	// No intersection
	if (c==0) inter=0;	// No interior intersection is possible now
	}
    return inter;
    }

int
Intersects(const IBox& a,const IBox& b)
    { return a.Intersects(b); }

int
IBox::Intersects(const Vector& A,const Vector& U)
   {
   double kmin=-AlpEpsilon,kmax=AlpInfty,k1,k2;
   int i,sgn;
   for (i=1;i<=n;i++) if ((sgn=AlpCompare(U(i),0.0))!=0)
      {
      if (sgn==1)		// U(i)>0
	 { k1=(m(i)-A(i))/U(i);k2=(M(i)-A(i))/U(i); }
      else
	 { k1=(M(i)-A(i))/U(i);k2=(m(i)-A(i))/U(i); }
      // We have k1<k2, the values of k when the ray
      // enters in, and leaves, the I-th slab.
      if (kmax<=k1 || k2<=kmin) return 0;
      if (k1>kmin) kmin=k1;
      if (k2<kmax) kmax=k2;
      }
   return 1;			// If we haven't returned before,
				// then we intersect the box.
   }

//
// Operations
// ==========
void
IBox::Point(const Vector& x,Vector& a) const
    {
    int i;
    for (i=n;i>0;i--) a(i)=0.5*((1.0+x(i))*M(i)+(1.0-x(i))*m(i));
    }

void
IBox::Center(Vector& a) const
    {
    int i;
    for (i=n;i>0;i--) a(i)=0.5*(M(i)+m(i));
    }

void
IBox::ChildBox(int child_id,IBox& b) const
   {
   int i,bit;
   for (i=bit=1;i<=n;i++,bit<<=1)
      if (child_id & bit)
	 { b.m(i)=m(i);b.M(i)=0.5*(M(i)+m(i)); }
      else
	 { b.M(i)=M(i);b.m(i)=0.5*(M(i)+m(i)); }
   }

int
IBox::ChildId(const Vector& a) const
   {
   int i,bit,id=0;
   for (i=bit=1;i<=n;i++,bit<<=1)
      if (2.0*a(i)<m(i)+M(i)) id|=bit;
   return id;
   }

void
Intersection(const IBox& a,const IBox& b,IBox& r)
    {
    int i;
    // Eliminate cases where a or b is empty
    if (a.is_empty || b.is_empty)
	{ r.is_empty=1;return; } // Intersection is empty : stop now !
    for (i=r.N();i>0;i--)
	{
	r.m(i)=(a.m(i)>b.m(i))?a.m(i):b.m(i); // g++ max operator
	r.M(i)=(a.M(i)<b.M(i))?a.M(i):b.M(i); // g++ min operator
	if (r.m(i)>r.M(i)+AlpEpsilon)
	    { r.is_empty=1;return; } // Intersection is empty : stop now !
	}
    r.is_empty=0;
    }

void
Union(const IBox& a,const IBox& b,IBox& r)
    {
    int i;
    // Eliminate cases where a or b is empty
    if (a.is_empty) { r=b;return; }
    if (b.is_empty) { r=a;return; }
    for (i=r.N();i>0;i--)
	{
	r.m(i)=(a.m(i)<b.m(i))?a.m(i):b.m(i); // g++ min operator
	r.M(i)=(a.M(i)>b.M(i))?a.M(i):b.M(i); // g++ max operator
	}
    }

IBox&
IBox::operator *= (const IBox& a)
    {
    int i;
    // Eliminate cases where a or *this is empty
    if (a.is_empty || is_empty)
	{ is_empty=1;return *this; } // Intersection is empty : stop now !
    for (i=n;i>0;i--)
	{
	m(i)=(m(i)>a.m(i))?m(i):a.m(i); // g++ max operator
	M(i)=(M(i)<a.M(i))?M(i):a.M(i); // g++ min operator
	if (m(i)>M(i)+AlpEpsilon)
	    { is_empty=1;return *this; } // Intersection is empty : stop now !
	}
    return *this;
    }

IBox&
IBox::operator += (const IBox& a)
    {
    int i;
    // Eliminate cases where a or *this is empty
    if (a.is_empty) return *this;
    if (is_empty) { *this=a;return *this; }
    for (i=n;i>0;i--)
	{
	m(i)=(m(i)<a.m(i))?m(i):a.m(i); // g++ min operator
	M(i)=(M(i)>a.M(i))?M(i):a.M(i); // g++ max operator
	}
    return *this;
    }

IBox&
IBox::operator += (const Vector& p)
    {
    int i;
    if (is_empty)
	{ is_empty=0;Copy(p,m);Copy(p,M);return *this; } // null set -> {p}
    for (i=n;i>0;i--)
	{
	if (p(i)>M(i)) M(i)=p(i);
	else if (p(i)<m(i)) m(i)=p(i);
	}
    return *this;
    }

/*
std::ostream& operator << (std::ostream& o,const IBox& b)
   {
   return o<<"Box("<<b.m<<","<<b.M<<")";
   }
*/
