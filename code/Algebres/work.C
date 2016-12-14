
//========================================================================
void SLSCGreg3dpg(SparseColMat *a,Vector **x,Vector **b, 
		int dim1, int dim2, int dim3, Scalaire tau,
                int *out, int proj,
		Scalaire arret,int itermin,int itermax,int & iter,
		Vector **r,int verbose)
     // FOR PARALLEL GEOMETRIES !!!
     // Sparse Least Squares Conjugate Gradient with Regularisation
     // allows us to solve a sparse linear least squares 
     // Min_{x} || A P x - b || +tau x^t PDP x
     // where D is the 3d Laplacian discretisation
     // x is suppose to be a 3D discretisation stored in the following way:
     // (*(x[k]))(i*dim2+j+1)=x(i,j,k)
     // Dx=6x(i,j,k)-(x(i+1,j,k)+x(i-1,j,k)+ .....)
     // P is a projector onto support information 
     //       (equal to Identity if noproj=1)
     //
     // *a is a SparseColMat of dimension DIM1A x DIM2A
     // **x are dim2  Vector of dimension DIM2
     // *out is an integer array of dimension dim2*DIM2 for support 
     //              information 
     //     introduction: if proj is equal to 0 then *out is suppose 
     //     to be zero. Else out[j*dim2+i-1]=1 means that (*(x[j]))(i-1) is nul
     //            (*out) is the complementary of the projector P
     // *b are dim2  Vector of dimension DIM1
     // *r are dim2  Vector of dimension DIM2 [residual (A^t A +tau D)x - A^t b ]
     // arret : STOP criterion
     //         if ||*r(k)|| <= arret or if ||*r(k)||/||*r(0)||
     //           where *r(k) is the residual at iteration k.
     // verbose ==1 then the program is verbose...
{
  Scalaire rkrk,rzrz, ri, norme, alpha, beta;
  int i,j,n;
  int fini; /* 0=false, 1=true */

  n=a->DIM2();
  cout <<"  DIMENSION n "<< n <<endl;
  if( (dim1*dim2) != n)
    {
      cout << " Err in SLSCGreg3d dim1,dim2,dim3: " << endl;
      cout << dim1 << " "  << dim2 << " " << dim2 << endl;
      exit(1);
    }

  Vector **d, **atad, **atb,**vinter;
  d = new Vector* [dim3];
  atad= new Vector* [dim3];
  atb= new Vector* [dim3];
  vinter= new Vector* [dim3];
  for(j=0;j<dim3;j++) {
    d[j]= new Vector(n);
    atad[j]=new Vector(n);
    atb[j]=new Vector(n);
    vinter[j]=new Vector(a->DIM1());
  }

// atb = A^t b
  for(j=0;j<dim3;j++) 
    ProdtransSparseColMatVector(a,b[j],atb[j]);
  if(proj)
    {
        for(j=0;j<dim3;j++) 
	  for(i=1;i<=n;i++)
	    if(out[j*n+i-1]) (*(atb[j]))(i)=0;
    }
      
// 
  iter=0;
  fini=0;

// r = P(A^t A+tau D)P x - P A^t b
  if(proj)
    {
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++)
	  if(out[j*n+i-1]) (*(x[j]))(i)=0;
    }
      
  for(j=0;j<dim3;j++) 
    ProdSparseColMatVector(a,x[j],vinter[j]);
  for(j=0;j<dim3;j++) 
    ProdtransSparseColMatVector(a,vinter[j],r[j]);
  AddLapla3Dgp(x,dim1,dim2,dim3, tau,r);
  if(proj)
    {
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++)
	  if(out[j*n+i-1]) (*(r[j]))(i)=0;
    }
      

  for(j=0;j<dim3;j++) 
    for(i=1;i<=n;i++)
      (*(r[j]))(i)=(*(r[j]))(i)-(*(atb[j]))(i);
// d = -r; rkrk=||r||^2
  rkrk=(double)0.;
  for(j=0;j<dim3;j++) 
    for(i=1;i<=n;i++){
      ri=(*(r[j]))(i);
      (*(d[j]))(i)=-ri;
      rkrk=rkrk+ ri * ri;
    }
  if((rkrk/(n*dim3))==0) fini=TRUE;
  if(verbose) cout << " Energie du residu (rkrk/n) : "<< rkrk/n << endl;
  rzrz=rkrk;

  while(fini==FALSE)
    {
      iter=iter+1;
      /* matadk */
      if(proj)
	{
	  for(j=0;j<dim3;j++) 
	    for(i=1;i<=n;i++)
	      if(out[j*n+i-1]) (*(d[j]))(i)=0;
	}
      for(j=0;j<dim3;j++) 
	ProdSparseColMatVector(a,d[j],vinter[j]);
      for(j=0;j<dim3;j++) 
	ProdtransSparseColMatVector(a,vinter[j],atad[j]);
      AddLapla3Dgp(d,dim1,dim2,dim3, tau, atad);
      if(proj)
	{
	  for(j=0;j<dim3;j++) 
	    for(i=1;i<=n;i++)
	      if(out[j*n+i-1]) (*(atad[j]))(i)=0;
	}
      /* ak=rktrk/dktadk */
      norme=(Scalaire) 0.;
      /* <ad,d> */
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++)
	  {
	    norme=norme + (*(d[j]))(i) * (*(atad[j]))(i);
	  }
      if(norme==0.)
        {
          fini=TRUE;
          return;
        }
      else
        {
          alpha=rkrk/norme;
        }
      /* xk+1=xk+akdk */
      /*      x(i)=x(i)+alpha*d(i); */
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++)
	  (*(x[j]))(i)=(*(x[j]))(i) + alpha * (*(d[j]))(i);
      /* rk+1=rk+akadk */
      /*      r(i)=r(i)+alpha*ad(i); */
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++)
	  (*(r[j]))(i)=(*(r[j]))(i) + alpha * (*(atad[j]))(i);
      
      /* bk=rk+1trk+1/rkrk */
      /* norme=norme+r(i)**2 */
      norme=(Scalaire) 0.0;
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++){
	  ri = (*(r[j]))(i);
	  norme=norme + ri * ri ;
	}
      beta=norme/rkrk;
      /* remise a jour de rkrk */
      rkrk=norme;
      /* dk+1=-rk+1+bk dk */
      /* d(i)=-r(i)+beta*d(i) */
      for(j=0;j<dim3;j++) 
	for(i=1;i<=n;i++)
	  (*(d[j]))(i) = -(*(r[j]))(i)+ beta * (*(d[j]))(i) ;
      if((rkrk/n)<=(arret*arret))fini=TRUE;
      if((rkrk/rzrz)<=arret)fini=TRUE;
      if(iter<itermin)fini=FALSE;
      if(iter>itermax)fini=TRUE;
      if(verbose) cout << " iteration no "<< iter;
      if(verbose) cout << " Energie du residu (rkrk/n) : "<< rkrk/n << endl;
    }
  for(j=0;j<dim3;j++) {
    delete d[j];
    delete atad[j];
    delete atb[j];
    delete vinter[j];
  }
  delete [] d;
  delete [] atad;
  delete [] atb;
  delete [] vinter;
}


void AddLapla3Dgp(Vector **x, int dim1,int dim2,int dim3, 
		Scalaire tau,Vector **y)
     // FOR PARALLEL GEOMETRIES !!!
     //  *y=*y+tau D *x
{
  int i,j,k;
  for(k=1;k<dim3-1;k++)
    for(i=1;i<dim1-1;i++)
      for(j=1;j<dim2-1;j++)
	{
	  (*(y[k]))(k*dim2+i+1)=(*(y[k]))(k*dim2+i+1)
	    + tau * (6* (*(x[k]))(i*dim2+j+1)
		     - (*(x[k-1]))(i*dim2+j+1)
		     - (*(x[k+1]))(i*dim2+j+1)
		     - (*(x[k]))(i*dim2+j+1+1)
		     - (*(x[k]))(i*dim2+j+1-1)
		     - (*(x[k]))((i+1)*dim2+j+1)
		     - (*(x[k]))((i-1)*dim2+j+1) );
	}

  for(j=0;j<dim2;j++)
    for(k=0;k<dim3;k++)
      {
	(*(y[k]))(j+1)=(*(y[k]))(j+1)
	  + tau * (*(x[k]))(j+1);
	(*(y[k]))((dim1-1)*dim2+j+1)=
	  (*(y[k]))((dim1-1)*dim2+j+1)
	    + tau * (*(x[k]))((dim1-1)*dim2+j+1);
      }

  for(i=0;i<dim1;i++)
    for(k=0;k<dim3;k++)
      {
	(*(y[k]))(i*dim2+1)=(*(y[k]))(i*dim2+1)
	  + tau *  (*(x[k]))(i*dim2+1);
	(*(y[k]))(i*dim2+dim2-1+1)= 
	  (*(y[(k]))(i*dim2+dim2-1+1)
	  + tau *  (*(x[k]))(i*dim2+dim2-1+1);
      }

  for(j=0;j<dim2;j++)
    for(i=0;i<dim1;i++)
      {
	(*(y[0]))(i*dim2+j+1)=(*(y[0]))(i*dim2+j+1)
	  + tau *  (*(x[0]))(i*dim2+j+1);
	(*(y[dim3-1]))(i*dim2+j+1)=
	  (*(y[dim3-1]))(i*dim2+j+1)
	  + tau *  (*(x[dim3-1]))(i*dim2+j+1);
      }
}
