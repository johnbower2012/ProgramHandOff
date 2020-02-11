#include "coshfunc.h"

#define PI 4.0*atan(1.0)

//////////////////////////
//    CDistribution     //
//////////////////////////

CDistribution::CDistribution(){
}
void CDistribution::GenX(int ab,double &x,double &weight){
}

//////////////////////////
//         CCOSH        //
//////////////////////////

CCosh::CCosh(){
  this->nmax=-1;
}
CCosh::CCosh(int nmaxset){
  this->Set_NMax(nmaxset);
}
void CCosh::Set_NMax(int nmaxset){
  this->nmax=nmaxset;
  Eigen::VectorXd x,y,ytest;
  int n,nprime,m,mprime;
  double norm;
  I.resize(2*nmax+2);
  I(0)=0.5*PI;
  I(1)=1.0;

  for(n=2;n<=2*nmax+1;n++){
    I(n)=I(n-2)*double(n-1)/double(n);
  }
  
  A = Eigen::MatrixXd::Zero(nmax+1,nmax+1);
  Z = Eigen::VectorXd::Zero(nmax+1);
  Eigen::MatrixXd K;

  A(0,0)=1.0;
  Z(0) = A(0,0)*I(0);
  for(n=1;n<=nmax;n++){
    A(n,n)=1.0;
    K = Eigen::MatrixXd::Zero(n,n);
    x = Eigen::VectorXd::Zero(n);
    y = Eigen::VectorXd::Zero(n);
    ytest = Eigen::VectorXd::Zero(n);
    for(nprime=0;nprime<n;nprime++){
      for(m=0;m<n;m++){
	K(nprime,m)=0.0;
	for(mprime=0;mprime<=nprime;mprime++){
	  K(nprime,m)+=I(m+mprime+1)*A(nprime,mprime);
	}
      }
    }
    for(nprime=0;nprime<n;nprime++){
      y(nprime)=0.0;
      for(mprime=0;mprime<=nprime;mprime++)
	y(nprime)+=I(n+mprime+1)*A(nprime,mprime);
    }
    x = K.colPivHouseholderQr().solve(y);
    for(m=0;m<n;m++)
      A(n,m)=-x(m);
    norm=0.0;
    for(m=0;m<=n;m++){
      for(mprime=0;mprime<=n;mprime++){
	norm+=A(n,m)*A(n,mprime)*I(m+mprime+1);
      }
    }
    if(norm<0.0){
      printf("norm(%d)=%g\n",n,norm);
      printf("negative norm, try making nmax smaller\n");
      exit(1);
    }
    x=x/sqrt(norm);
    y=y/sqrt(norm);
    //ytest=K*x-y;                                                                                                              
    //ytest.print("ytest:");                                                                                                    
    for(m=0;m<=n;m++){
      A(n,m)=A(n,m)/sqrt(norm);
      Z(n) += A(n,m)*I(m);
    }
   }

  seed = 123456;
  rng = std::mt19937(seed);
  uniform_dist = std::uniform_real_distribution<double>(0,1);
}
double CCosh::GetG(int n,double x){
  double g,coshx;
  int m;
  coshx=cosh(x);
  g=0.0;
  for(m=0;m<n+1;m++){
    g+=this->A(n,m)*pow(coshx,-m-1);
  }
  return g;
}
double CCosh::GetRandomX(){
  //random number from (2/pi)/cosh(x)
  return asinh(tan(0.5*PI*uniform_dist(rng)));
  //  return acosh(1.0/uniform_dist(rng));
}
void CCosh::TestOverlap(){
  int n,nprime,nmax=this->nmax;
  for(n=0;n<nmax+1;n++){
    for(nprime=0;nprime<nmax+1;nprime++){
      this->CalcOverlap(n,nprime);
      this->CalcOverlapNumerical(n,nprime);
      printf("----------------------\n");
    }
  }
}
void CCosh::CalcOverlap(int n,int nprime){
  double overlap=0.0;
  int m,mprime;
  for(m=0;m<=n;m++){
    for(mprime=0;mprime<=nprime;mprime++){
      overlap+=A(n,m)*I(m+mprime+1)*A(nprime,mprime);
    }
  }
  printf("overlap(%d,%d)=%g\n",n,nprime,overlap);
}
void CCosh::CalcOverlapNumerical(int n,int nprime){
  double overlap=0.0,g,gprime;
  double x,dx=0.0002,xmax=35.0;
  for(x=0.5*dx;x<xmax;x+=dx){
    g=this->GetG(n,x);
    gprime=this->GetG(nprime,x);
    overlap+=g*gprime*dx;
  }
  if(abs(overlap)<5.0E-9)
    overlap=0.0;
  printf("numerical: overlap(%d,%d)=%g\n",n,nprime,overlap);
}
void CCosh::Integrate(double x){
  printf("Testing getting random no. from dist (2/pi)/cosh(x)\n");
  double answer=0,dx=0.0001,xx;
  for(xx=0.5*dx;xx<x;xx+=dx){
    answer+=dx/cosh(xx);
  }
  printf("%g =? %g\n",answer,atan(sinh(x)));

  double y,dy=0.01;
  int ibin,nbins=1000;
  double *prob=new double[nbins];
  int n,nmax=1000000;
  for(n=0;n<nmax;n++){
    y=GetRandomX();
    ibin=lrint(floor(y/dy));
    if(ibin<nbins)
      prob[ibin]+=1.0/(nmax*dy);
  }
  for(ibin=0;ibin<nbins;ibin++){
    y=(ibin+0.5)*dy;
    printf("%5d %6.3f %10.7f =? %10.7f\n",ibin,y,prob[ibin],(2.0/PI)/cosh(y));
  }
}

//////////////////////////
//       CDISTCOSH      //
//////////////////////////

CDistCosh::CDistCosh() : coshcalc(){
}
CDistCosh::CDistCosh(Eigen::VectorXd Width, Eigen::MatrixXd G) : coshcalc(){
  int ng=G.cols();
  this->width = Width;
  this->g = G;

  coshcalc.Set_NMax(ng+1);
}
void CDistCosh::Set_WG(Eigen::VectorXd Width, Eigen::MatrixXd G){
  this->width = Width;
  this->g = G;
  int ng=G.cols();
  if(coshcalc.nmax!=ng){
    coshcalc.Set_NMax(ng);
  }
}
void CDistCosh::GenX(int ab, double &x, double &weight){
  int n,ng=this->coshcalc.nmax;
  x=coshcalc.GetRandomX();
  //  double Z=0.0;
  weight=0.0;
  for(n=0;n<ng+1;n++){
    weight+=coshcalc.GetG(n,x/width(ab))*g(ab,n);
    //    Z += coshcalc.Z(n)*g(ab,n);
  }
  weight *= 0.5*PI/(coshcalc.GetG(0,x));//*g(ab,0));
  x*=width(ab);
  //  x*=PI/2.0/1.83193118835; //should return 1.0 on avg
  //2.0*Catalan: 1.83193118835
}
Eigen::MatrixXd CDistCosh::GenG0(int ab, Eigen::MatrixXd G){
  int nmax=G.cols()/ab-1,
    samples=G.rows();
  this->coshcalc.Set_NMax(nmax+1);
  Eigen::MatrixXd newG = Eigen::MatrixXd::Zero(samples,ab*(nmax+2));
  for(int i=0;i<samples;i++){
    for(int j=0;j<ab;j++){
      newG(i,j*(nmax+2)) = G(i,j*(nmax+1));
      newG(i,j*(nmax+2)+1) = 1.0/G(i,j*(nmax+1));
      for(int k=1;k<nmax+1;k++){
	newG(i,j*(nmax+2)+k+1) = G(i,j*(nmax+1)+k);
	newG(i,j*(nmax+2)+1) -= G(i,j*(nmax+1)+k)*this->coshcalc.Z(k);
      }
      newG(i,j*(nmax+2)+1) /= (this->coshcalc.Z(0));
    }
  }
  return newG;
}
Eigen::MatrixXd CDistCosh::FixG0(int ab, Eigen::MatrixXd G){
  int nmax=G.cols()/ab-2,
    samples=G.rows();
  this->coshcalc.Set_NMax(nmax+1);
  Eigen::MatrixXd newG = Eigen::MatrixXd::Zero(samples,ab*(nmax+2));
  for(int i=0;i<samples;i++){
    for(int j=0;j<ab;j++){
      newG(i,j*(nmax+2)) = G(i,j*(nmax+2));
      newG(i,j*(nmax+2)+1) = 1.0/G(i,j*(nmax+2));
      for(int k=1;k<nmax+1;k++){
	newG(i,j*(nmax+2)+k+1) = G(i,j*(nmax+2)+k);
	newG(i,j*(nmax+2)+1) -= G(i,j*(nmax+2)+k)*this->coshcalc.Z(k);
      }
      newG(i,j*(nmax+2)+1) /= (this->coshcalc.Z(0));
    }
  }
  return newG;
}
Eigen::MatrixXd CDistCosh::Functions(int n, double x){
  int ab = this->width.size();
  int nmax = this->g.cols();
  double dx=x/double(n), xx=0.0;
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n+1,ab+1);
  for(int xi=0;xi<n+1;xi++){
    M(xi,0) = xx;
    for(int species=0;species<ab;species++){
      for(int m=0;m<nmax;m++){
	M(xi,species+1) += coshcalc.GetG(m,xx/width(species))*g(species,m);
      }
    }
    xx += dx;
  }
  return M;
}
Eigen::MatrixXd CDistCosh::FunctionSet(int grid, double x, int samples, int ab, int nmax, Eigen::MatrixXd G){
  int Grid=grid+1,
    Cols=samples*ab+1;
  Eigen::MatrixXd Func = Eigen::MatrixXd::Zero(Grid,Cols);
  Eigen::VectorXd 
    width(ab);
  Eigen::MatrixXd 
    g(ab,nmax+1),
    m;
  this->coshcalc.Set_NMax(nmax+1);
  for(int i=0;i<samples;i++){
    for(int j=0;j<ab;j++){
      width(j) = G(i,j*(nmax+2));
      for(int k=0;k<nmax+1;k++){
	g(j,k) = G(i,j*(nmax+2)+k+1);
      }
    }
    this->Set_WG(width,g);
    m=this->Functions(grid,x);

    for(int j=0;j<grid+1;j++){
      for(int k=0;k<ab;k++){
	Func(j,1+i*ab+k) = m(j,k+1);
      }
      Func(j,0) = m(j,0);
    }
  }
  return Func;
}
double CDistCosh::Integrate(int ab){
  int m;
  int nmax = coshcalc.nmax;
  double sum=0.0;
  for(m=0;m<nmax+1;m++){
    sum += coshcalc.Z(m)*g(ab,m);
  }
  sum *= width(ab);
  return sum;
}
double CDistCosh::IntegrateNumerical(int ab){
  double
    x=100.0,
    xx=0.0,
    dx=0.0001,
    sum=0.0;
  int
    m;
  int nmax = coshcalc.nmax;
  for(m=0;m<nmax+1;m++){
    for(xx=0.0;xx<x;xx+=dx){
      sum += this->coshcalc.GetG(m,xx)*g(ab,m);
    }
  }
  sum *= dx;
  sum *= width(ab);
  return sum;
}


