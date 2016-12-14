// Linear algebre with C++
// Real tensors
// $Header: /home/bainvil/Modules/alp/RCS/rtensbase.C,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#include <alp.h>

// BUGS : Share est foireux, je m'en occuperai une autre fois :-)

//
// Coordinates <-> array index
// ===========================
int Tensor::CtoI(const int *c) const
    {
    int j,i=0;
    for (j=0;j<d;j++) i+=(c[j]-1)*s[j];
    return i;
    }
void Tensor::ItoC(int i,int *c) const
    {
    int j;
    for (j=0;j<d;j++) c[j]=1+((i/s[j])%n[j]);
    }
//
// Memory allocation
// =================
void Tensor::AllocD()
    { shared=0;n=new int[d];s=new int[d]; }
void Tensor::AllocX()
    {
    int j;
    s[0]=1;for (j=1;j<d;j++) s[j]=n[j-1]*s[j-1];
    size=s[d-1]*n[d-1];
    x=new double[size];
#if 1
    cerr<<"Tensor::AllocX()"<<endl;
    cerr<<"  order = "<<d<<endl;
    cerr<<"  size = "<<size<<endl;
    cerr<<"  increments =";
    for (j=0;j<d;j++) cerr<<" "<<s[j];
    cerr<<endl<<"  dimensions =";
    for (j=0;j<d;j++) cerr<<" "<<n[j];
    cerr<<endl;
#endif
    }
//
// Protected constructor (for shared tensors)
// ==========================================
Tensor::Tensor(int order,char a) : d(order)
    { AllocD();shared=1; }
//
// Constructors
// ============
Tensor::Tensor(int n1) : d(1)
    { AllocD();n[0]=n1;AllocX(); }
Tensor::Tensor(int n1,int n2) : d(2)
    { AllocD();n[0]=n1;n[1]=n2;AllocX(); }
Tensor::Tensor(int n1,int n2,int n3) : d(3)
    { AllocD();n[0]=n1;n[1]=n2;n[2]=n3;AllocX(); }
Tensor::Tensor(int n1,int n2,int n3,int n4) : d(4)
    { AllocD();n[0]=n1;n[1]=n2;n[2]=n3;n[3]=n4;AllocX(); }
Tensor::Tensor(int n1,int n2,int n3,int n4,int n5) : d(5)
    { AllocD();n[0]=n1;n[1]=n2;n[2]=n3;n[3]=n4;n[4]=n5;AllocX(); }
Tensor::Tensor(int n1,int n2,int n3,int n4,int n5,int n6) : d(6)
    { AllocD();n[0]=n1;n[1]=n2;n[2]=n3;n[3]=n4;n[4]=n5;n[5]=n6;AllocX(); }
Tensor::Tensor(const Tensor& t) : d(t.d)
    {
    int j;
    AllocD();
    for (j=0;j<d;j++) n[j]=t.n[j];
    AllocX();
    Copy(t,*this);
    }
//
// Destructor
// ==========
Tensor::~Tensor()
    {
    if (!shared) delete [] x;
    delete [] n;delete [] s;
    }
//
// Memory operations
// =================
void Copy(const Tensor& a,Tensor& b)
    {
    if (b.d==1)
	{
	int i,pa=0,pb=0,sa=a.s[0],sb=b.s[0];
	for (i=b.n[0];i>0;i--) { b.x[pb]=a.x[pa];pa+=sa;pb+=sb; }
	}
    else
	{
	int i;
	for (i=b.n[0];i>0;i--) Copy(a.Share1(i),b.Share1(i));
	}
    }
void Swap(Tensor& a,Tensor& b)
    { Tensor u(a);Copy(b,a);Copy(u,b); }
Tensor& Tensor::operator = (const Tensor& a)
    { Copy(a,*this);return *this; }
//
// Initialisations
// ===============
void Tensor::Zero()
    { for (int i=0;i<size;i++) x[i]=double(i); }
//
// Sharing parts of the tensor
// ===========================
SharedTensor Tensor::Share(int k,int *si,int *c)
    {
    int j,*dd=new int[k],*ds=new int[k];
    SharedTensor *u;
    for (j=0;j<k;j++) { dd[j]=si[j];ds[j]=s[dd[j]];c[dd[j]]=0; }
    u=new SharedTensor(k,x+CtoI(c),dd,ds);
    delete [] dd;
    delete [] ds;
    return *u;
    }
const SharedTensor Tensor::Share(int k,int *si,int *c) const
    {
    int j,*dd=new int[k],*ds=new int[k];
    SharedTensor *u;
    for (j=0;j<k;j++) { dd[j]=si[j];ds[j]=s[dd[j]];c[dd[j]]=0; }
    u=new SharedTensor(k,x+CtoI(c),dd,ds);
    delete [] dd;
    delete [] ds;
    return *u;
    }
SharedTensor Tensor::Share1(int i)
    { return SharedTensor(d-1,x+(s[0]*(i-1)),n+1,s+1); }
const SharedTensor Tensor::Share1(int i) const
    { return SharedTensor(d-1,x+(s[0]*(i-1)),n+1,s+1); }
//
// Total aggregation
// =================
double Tensor::Aggregate(AggregationFunction &F) const
    {
    int j,*i=new int[d],index;
    cerr<<(*this)<<endl;
    F.Initialize();
    for (j=0;j<d;j++) i[j]=1;
    index=0;
    do
	{
	F.Aggregate(x[index]);
	j=d-1;
	do
	    {
	    index+=s[j];
	    if ((++i[j])>n[j]) { index-=n[j]*s[j];i[j--]=1; }
	    else break; // On n'a pas fini cette dimension, on continue
	    }
	while (j>=0);
	}
    while (j>=0);
    delete [] i;
    return F.Result();
    }
//
// Partial aggregation
// ===================
void Tensor::Aggregate(int k,int *si,int *di,AggregationFunction &F,Tensor& a) const
    {
    int *i=new int[d],j,aj,index;

    for (j=0;j<d;j++) i[j]=1; // Initialise le compteur
    for (j=0;j<k;j++) di[si[j]]=-1; // Marque les images des indices de la somme
    index=0;
    do
	{
	a.x[index]=Share(k,si,i).Aggregate(F);
	j=d-1;
	do
	    {
	    if ((aj=di[j])<0) j--;
	    else
		{
		index+=a.s[aj];
		if ((++i[j])>n[j]) { index-=a.s[aj]*a.n[aj];i[j--]=1; }
		else break;
		}
	    }
	while (j>=0);
	}
    while (j>=0);
    delete [] i;
    }
//
// Text output
// ===========
ostream& operator << (ostream& o,const Tensor& a)
    {
    int j;
    o<<endl<<"Tensor["<<a.d;
    for (j=0;j<a.d;j++) cout<<((j==0)?';':',')<<a.n[j];
    for (j=0;j<a.d;j++) cout<<((j==0)?';':',')<<a.s[j];
    o<<']';
    if (a.d==1)
	{ for (j=0;j<a.n[0];j++) o<<((j==0)?'(':',')<<a.x[j*a.s[0]]; o<<')'; }
    else
	{ for (j=1;j<=a.n[0];j++) o<<((j==1)?'(':',')<<a.Share1(j); o<<')'; }
    return o;
    }
