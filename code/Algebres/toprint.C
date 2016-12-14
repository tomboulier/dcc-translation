void SLSCGreg3dpg(SparseColMat *a,Vector **x,Vector **b, 
		int dim1, int dim2, int dim3, Scalaire tau,
                int *out, int proj,
		Scalaire arret,int itermin,int itermax,int & iter,
		  Vector **r,int verbose, Scalaire ParallelRatio)
     // FOR PARALLEL GEOMETRIES !!!
     // Sparse Least Squares Conjugate Gradient with Regularisation
     // allows us to solve a sparse linear least squares 
     // Min_{x} || A P x - b || +tau x^t PDP x
     // where D is the 3d Laplacian discretisation
     // x is suppose to be a 3D discretisation stored in the following way:
     // (*(x[j]))(k*dim1+i+1)=x(i,j,k)
     // Dx=6x(i,j,k)-(x(i+1,j,k)+x(i-1,j,k)+ .....)
     // P is a projector onto support information 
     //       (equal to Identity if noproj=1)
     //
     // *a is a SparseColMat of dimension DIM1A x DIM2A 
     // **x are dim1  Vectors of dimension DIM2A=dim3*dim2
     // *out is an integer array of dimension dim1*DIM2A (dim1*dim2*dim3)
     //     for support  information 
     //     introduction: if proj is equal to 0 then *out is suppose 
     //     to be zero. Else out[j*dim2+i-1]=1 means that (*(x[j]))(i-1) is nul
     //            (*out) is the complementary of the projector P
      //
     // Thu Oct  9 19:04:38 CEST 2008
     // *b are NbSclice  Vector of dimension DIM1A
     // *r are dim1  Vector of dimension DIM1A [residual (A^t A +tau D)x - A^t b ]
     //        where NbSlice=dim1*ParallelRatio (see below) [NbSlice is a local variable]
     // arret : STOP criterion
     //         if ||*r(k)|| <= arret or if ||*r(k)||/||*r(0)||
     //           where *r(k) is the residual at iteration k.
     // verbose ==1 then the program is verbose...
     //
     // Thu Oct  9 19:04:38 CEST 2008
     // ParallelRatio is the sampling ratio between NbSlice/dim1


{
  Scalaire rkrk,rzrz, ri, norme, alpha, beta;
  int i,j,k,n;
  int fini; /* 0=false, 1=true */

  int NbSlice=(int) rint(dim1*ParallelRatio);
  std::cout << " dim1= "<<dim1<<";  NbSlice= "<<NbSlice<<std::endl; 
  if((ParallelRatio<=0)||(ParallelRatio>1)){
    std::cout <<" ParallelRatio= "<<ParallelRatio << 
      "  Should be between 0 and 1 !!!"<<std::endl;
  }
  Scalaire InvParallelRatio=1./ParallelRatio;

  n=a->DIM2();
  cout <<"  DIMENSION n "<< n <<endl;
  if( (dim3*dim2) != n)
    {
      cout << " Err in SLSCGreg3d dim1,dim2,dim3: " << endl;
      cout << dim1 << " "  << dim2 << " " << dim3 << endl;
      exit(1);
    }

  Vector **d, **atad, **atb,**vinter;
  atb= new Vector* [dim1];
  d = new Vector* [dim1];
  atad= new Vector* [dim1];
  vinter= new Vector* [NbSlice];
  for(k=0;k<NbSlice;k++) {
    vinter[k]=new Vector(a->DIM1());
  }
  for(j=0;j<dim1;j++) {
    atb[j]=new Vector(n);
    d[j]= new Vector(n);
    atad[j]=new Vector(n);
    for(i=1;i<=n;i++){
      (*(atb[j]))(i)=0;
      (*(d[j]))(i)=0;
      (*(atad[j]))(i)=0;
    }
  }


// atb = A^t b
  for(k=0;k<NbSlice;k++) {
    j=(int) round( ((2*k+1)*InvParallelRatio -1)/2.);
    ProdtransSparseColMatVector(a,b[k],atb[j]);
  }
  if(proj){
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++)
	if(out[j*n+i-1]) (*(atb[j]))(i)=0;
  }
      
// 
  iter=0;
  fini=0;

// r = P(A^t A+tau D)P x - P A^t b
  if(proj)    {
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++)
	if(out[j*n+i-1]) (*(x[j]))(i)=0;
  }
      
  for(k=0;k<NbSlice;k++) {
    j=(int) round( ((2*k+1)*InvParallelRatio -1)/2.);
    ProdSparseColMatVector(a,x[j],vinter[k]);
  }
  for(k=0;k<NbSlice;k++) {
    j=(int) round( ((2*k+1)*InvParallelRatio -1)/2.);
    ProdtransSparseColMatVector(a,vinter[k],r[j]);
  }
  AddLapla3Dgp(x,dim1,dim2,dim3, tau,r);
  if(proj)
    {
      for(j=0;j<dim1;j++) 
	for(i=1;i<=n;i++)
	  if(out[j*n+i-1]) (*(r[j]))(i)=0;
    }
      

  for(j=0;j<dim1;j++) 
    for(i=1;i<=n;i++)
      (*(r[j]))(i)=(*(r[j]))(i)-(*(atb[j]))(i);
// d = -r; rkrk=||r||^2
  rkrk=(double)0.;
  for(j=0;j<dim1;j++) 
    for(i=1;i<=n;i++){
      ri=(*(r[j]))(i);
      (*(d[j]))(i)=-ri;
      rkrk=rkrk+ ri * ri;
    }
  if((rkrk/(n*dim1))==0) fini=TRUE;
  if(verbose) cout << " Energie du residu (rkrk/(n*dim1)) : "<< rkrk/(n*dim1) << endl;
  if(verbose) std::cout << " n : "<< n << endl;
  if(verbose) std::cout << " dim1 : "<< dim1 << endl;
  rzrz=rkrk;

  while(fini==FALSE){
    iter=iter+1;
    /* matadk */
    if(proj){
      for(j=0;j<dim1;j++) 
	for(i=1;i<=n;i++)
	  if(out[j*n+i-1]) (*(d[j]))(i)=0;
    }
    for(k=0;k<NbSlice;k++) {
      j=(int) round( ((2*k+1)*InvParallelRatio -1)/2.);
      ProdSparseColMatVector(a,d[j],vinter[k]);
    }
    for(k=0;k<NbSlice;k++) {
      j=(int) round( ((2*k+1)*InvParallelRatio -1)/2.);
      ProdtransSparseColMatVector(a,vinter[k],atad[j]);
    }
    //ICICI
    AddLapla3Dgp(d,dim1,dim2,dim3, tau, atad);
    if(proj){
      for(j=0;j<dim1;j++) 
	for(i=1;i<=n;i++)
	  if(out[j*n+i-1]) (*(atad[j]))(i)=0;
    }
    /* ak=rktrk/dktadk */
    norme=(Scalaire) 0.;
    /* <ad,d> */
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++)  {
	norme=norme + (*(d[j]))(i) * (*(atad[j]))(i);
      }
    if(norme==0.)   {
      fini=TRUE;
      return;
    }
    else{
      alpha=rkrk/norme;
    }
    /* xk+1=xk+akdk */
    /*      x(i)=x(i)+alpha*d(i); */
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++)
	(*(x[j]))(i)=(*(x[j]))(i) + alpha * (*(d[j]))(i);
    /* rk+1=rk+akadk */
    /*      r(i)=r(i)+alpha*ad(i); */
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++)
	(*(r[j]))(i)=(*(r[j]))(i) + alpha * (*(atad[j]))(i);
      
    /* bk=rk+1trk+1/rkrk */
    /* norme=norme+r(i)**2 */
    norme=(Scalaire) 0.0;
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++){
	ri = (*(r[j]))(i);
	norme=norme + ri * ri ;
      }
    beta=norme/rkrk;
    /* remise a jour de rkrk */
    rkrk=norme;
    /* dk+1=-rk+1+bk dk */
    /* d(i)=-r(i)+beta*d(i) */
    for(j=0;j<dim1;j++) 
      for(i=1;i<=n;i++)
	(*(d[j]))(i) = -(*(r[j]))(i)+ beta * (*(d[j]))(i) ;
    if((rkrk/(n*dim1))<=(arret*arret))fini=TRUE;
    if((rkrk/rzrz)<=arret)fini=TRUE;
    if(iter<itermin)fini=FALSE;
    if(iter>=itermax)fini=TRUE;
    if(verbose) cout << " iteration no "<< iter;
    if(verbose) cout << " Energie du residu (rkrk/(n*dim1)) : "
		     << rkrk/(n*dim1) << endl;
  }

  for(j=0;j<dim1;j++) {
    delete d[j];
    delete atad[j];
    delete atb[j];
  }
  for(k=0;k<NbSlice;k++) 
    delete vinter[k];

  delete [] d;
  delete [] atad;
  delete [] atb;
  delete [] vinter;
}

