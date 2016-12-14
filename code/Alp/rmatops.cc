// Linear algebra with C++
// Real matrixes : high level operators
//
// $Header: /home/bainvil/Modules/alp/RCS/rmatops.C,v 1.5 1995/03/16 18:15:21 bainvil Exp bainvil $

#include <alp.h>

#if 0
void SingularValueDecomposition(const Matrix &a,Matrix &u,Vector &sigma,Matrix &vt)
    {
    cerr<<"SVD not implemented"<<endl;
    }
#endif

double Det(const Matrix &a)
    {
    int n=a.R();
    if (n==1) return a(1,1);
    else if (n==2) return a(1,1)*a(2,2)-a(1,2)*a(2,1);
    else if (n==3) return (a(1,1)*a(2,2)*a(3,3)
			   +a(2,1)*a(3,2)*a(1,3)
			   +a(3,1)*a(1,2)*a(2,3)
			   -a(1,1)*a(3,2)*a(2,3)
			   -a(2,1)*a(1,2)*a(3,3)
			   -a(3,1)*a(2,2)*a(1,3));
    else
	{
	Matrix lu(n);
	Permutation p(n),q(n);
	int i;
	double s;

	LU(a,p,q,lu);
	for (i=1,s=1.0;i<=n;i++) s*=lu(i,i);
	return s;
	}
    }

int Solve(const Matrix& a,const Vector& b,Vector& x)
    {
    int n=a.R(),i,j;
    Permutation p(n),q(n);
    Matrix lu(n);
    Vector v(n);
    double s;

    (void)LU(a,p,q,lu);   // decompose a en p.a.q = l.u
    // Alors, a.x = b s'ecrit p.a.q.q^-1.x = p.b,
    // soit l.u.q^-1.x = p.b soit encore
    // d'abord l.v = p.b et ensuite u.q^-1(x) = v
    // Resout l.v = p.b
    for (i=1;i<=n;i++)
	{
	for (j=1,s=b(p(i));j<i;j++) s-=lu(i,j)*v(j);
	v(i)=s;  // Puisque l(i,i)=1
	}
    // Resout u.q^-1.x = v
    for (i=n;i>=1;i--)
	{
	for (j=n,s=v(i);j>i;j--) s-=lu(i,j)*x(q(j));
	if (fabs(lu(i,i))>1.0e-8) x(q(i))=s/lu(i,i);
	else if (fabs(s)<1.0e-8) x(q(i))=0.0;  // Non determine, on met 0
	else return 0;  // Pas de solution
	}
    return 1;
    }

int LU(const Matrix& a,Permutation& p,Permutation& q,Matrix& lu)
    {
    int n=a.R(),i,j,r,c,pivr,pivc,sing=0;
    double g,h,k,piv,pivh=0.0;
    Matrix s(a);  // On travaille sur une copie de a, permet &a==&lu

    // Initialise les permutations sur les lignes et les colonnes
    p.Id();q.Id();
    // Calcule un ordre de grandeur des coefficients de la matrice
    g=Norm(s)/double(n*n);

    // Boucle sur les lignes
    for (i=1;i<=n;i++)
	{
	// Cherche le meilleur pivot dans la matrice restante
	pivr=pivc=-1;
	for (c=i;c<=n;c++)
	for (r=i;r<=n;r++)
	    {
	    piv=s(p(r),q(c));
	    if (fabs(piv)>1.0e-8*g)
		{
		h=fabs(piv-g);		       // Cherche le plus proche de g
		if ((pivr<0) || (h<pivh)) { pivh=h;pivr=r;pivc=c; }
		}
	    }
	// Si pas de pivot non nul, la matrice n'est pas reguliere, on passe...
	if (pivr<0) { sing=1;continue; }
	// Echange (virtuellement) les lignes et les colonnes pour positionner
	// le pivot en (i,i)
	p.Exch(i,pivr);pivr=p(i);
	q.Exch(i,pivc);pivc=q(i);
	// Annule la colonne du pivot dans les lignes restantes
	h=1.0/s(pivr,pivc);      // Inverse du pivot
	for (r=i+1;r<=n;r++)
	    {
	    j=p(r);
	    k=h*s(j,pivc);
	    s(j,pivc)=k;    // Construit la matrice L...
	    if (fabs(k)>1.0e-8*g)
		{
		for (c=i+1;c<=n;c++) s(j,q(c))-=k*s(pivr,q(c));  // Et un peu la matrice U...
		}
	    }
	}

    // Calcule le resultat non permute dans lu
    for (r=1;r<=n;r++)
    for (c=1;c<=n;c++)
	{ lu(r,c)=s(p(r),q(c)); }

    return sing;
    }

int Invert(const Matrix& a,Matrix& b)
    {
    int n=a.R(),i,j,r,c,pivr,pivc;
    double g,h,k,piv,pivh=0.0;
    Permutation p(n),q(n);
    Matrix s(a),t(n);  // On travaille sur une copie de a

    // Initialise les permutations sur les lignes et les colonnes
    p.Id();q.Id();
    // Calcule un ordre de grandeur des coefficients de la matrice
    g=Norm(s)/double(n*n);
    // Initialise le resultat permute t
    t.Id();

    // Boucle sur les lignes
    for (i=1;i<=n;i++)
	{
	// Cherche le meilleur pivot dans la matrice restante
	pivr=pivc=-1;
	for (c=i;c<=n;c++)
	for (r=i;r<=n;r++)
	    {
	    piv=s(p(r),q(c));
	    if (fabs(piv)>1.0e-8*g)
		{
		h=fabs(piv-g);		       // Cherche le plus proche de g
		if ((pivr<0) || (h<pivh)) { pivh=h;pivr=r;pivc=c; }
		}
	    }
	// Si pas de pivot non nul, la matrice n'est pas reguliere, on arrete
	if (pivr<0) return 1;
	// Echange (virtuellement) les lignes et les colonnes pour positionner
	// le pivot en (i,i)
	p.Exch(i,pivr);pivr=p(i);
	q.Exch(i,pivc);pivc=q(i);
	// Divise la ligne du pivot par le pivot
	h=1.0/s(pivr,pivc);      // Inverse du pivot
	for (c=i+1;c<=n;c++) s(pivr,q(c))*=h;
	for (c=1;c<=n;c++) t(pivr,q(c))*=h;
	// Annule la colonne du pivot dans toutes les lignes (sauf celle du pivot)
	for (r=1;r<=n;r++) if (r!=i)
	    {
	    j=p(r);
	    k=s(j,pivc);
	    if (fabs(k)>1.0e-8*g)
		{
		for (c=i+1;c<=n;c++) s(j,q(c))-=k*s(pivr,q(c));
		for (c=1;c<=n;c++) t(j,c)-=k*t(pivr,c);
		}
	    }
	}

    // Calcule le resultat non permute dans b
    for (r=1;r<=n;r++)
    for (c=1;c<=n;c++)
	{ b(q(r),c)=t(p(r),c); }

    return 0;
    }

// Fonction interne qui effectue une rotation de la methode de Jacobi
// en annulant deux termes non diagonaux de la matrice symetrique a
void _JacobiIteration_AlpPrivate(Matrix& a,Matrix& omega,int p,int q)
    {
    double api,aqi,oip,oiq;
    double t,v,c,s;
    int i,n=a.R();

    v=(a(q,q)-a(p,p))/(2.0*a(p,q));  // Cotangente de l'angle de rotation
    if (v<0.0) t=-v-sqrt(1.0+v*v); else t=-v+sqrt(1.0+v*v);  // Tangente de l'angle
    c=1.0/sqrt(1.0+t*t);  // Cosinus de l'angle
    s=t*c;  // Sinus de l'angle

    // Met a jour la matrice a (matrice symetrique a diagonaliser)
    for (i=1;i<=n;i++)
	{
	if (i!=p && i!=q)
	    {
	    api=a(p,i);aqi=a(q,i);
	    a(p,i)=a(i,p)=c*api-s*aqi;
	    a(q,i)=a(i,q)=c*aqi+s*api;
	    }
	}
    v=t*a(p,q);
    a(p,p)-=v;
    a(q,q)+=v;
    a(p,q)=a(q,p)=0.0;

    // Met a jour la matrice omega (cumul des rotations)
    for (i=1;i<=n;i++)
	{
	oip=omega(i,p);oiq=omega(i,q);
	omega(i,p)=c*oip-s*oiq;
	omega(i,q)=s*oip+c*oiq;
	}
    }

void SolveEigenProblem(const Matrix& a,Vector& eigenvalues,Matrix& eigenvectors)
    {
    int i,j,p=0,q=0,n=a.R();
    double seuil,mx;
    Matrix s(a);  // On travaille sur une copie de a
    eigenvectors.Id();  // Initialise la base de vecteurs propres
    // Calcule le seuil a atteindre en fonction de la norme de la matrice
    seuil=1.0e-10*Norm(s)/double(n*n);
    while (1)
	{
	// Cherche (p,q) dans la moitie superieure de la matrice s
	mx=-1.0;
	for (i=1;i<n;i++)
	for (j=i+1;j<=n;j++)
	    {
	    if (fabs(s(i,j))>mx) { mx=fabs(s(i,j));p=i;q=j; }
	    }
	if (mx<seuil) break;  // Sort de la boucle si tout vaut "zero"
	_JacobiIteration_AlpPrivate(s,eigenvectors,p,q);  // Effectue la rotation adequate
	}
    s.GetDiag(eigenvalues);  // Recupere les n valeurs propres dans un vecteur
    }

#if 0
void
LargestEigenValue(const Matrix& a,double& eigenvalue,Vector& eigenvector)
   {
   cerr<<"LargestEigenValue not implemented"<<endl;
   }
#endif
