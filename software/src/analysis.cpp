#include "analysis.h"

void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();

  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

  matrix.conservativeResize(numRows,numCols);
}

void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;

  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

  matrix.conservativeResize(numRows,numCols);
}

CAnalysis::CAnalysis(){
}
void CAnalysis::ComputeMean(){
  int rows=this->Data.rows(),
    cols=this->Data.cols();
  this->Mean = Eigen::VectorXd::Zero(rows);
  for(int row=0;row<rows;row++){
    for(int col=0;col<cols;col++){
      this->Mean(row) += this->Data(row,col);
    }
    this->Mean(row) /= (double) cols;
  }
}
void CAnalysis::SumErrorInQuadrature(Eigen::MatrixXd Error2){
  int 
    rows = this->Error.rows(),
    cols = this->Error.cols();
  if(rows==0 && cols==0){
    this->Error = Error2;
  }
  else if(Error2.rows() != rows || Error2.cols() != cols){
      printf("SumInQuadrature Failed. Either A.rows != B.rows or A.cols != B.cols\n");
      printf("A.rows=%d. A.cols=%d. B.rows=%d. B.cols=%d\n",rows,cols,(int)Error2.rows(),(int)Error2.cols());
      exit(1);
  }
  else{
    Eigen::MatrixXd sum = Eigen::MatrixXd::Zero(rows,cols);
    for(int i=0;i<rows;i++){
      for(int j=0;j<cols;j++){
	  sum(i,j) = sqrt(this->Error(i,j)*this->Error(i,j) + Error2(i,j)*Error2(i,j));
	  if(sum(i,j) == 0.0) sum(i,j) = 10e-5;
	}
    }
    this->Error=sum;
  }
}
void CAnalysis::SumErrorInQuadrature(Eigen::VectorXd error2){
  int 
    rows = this->Error.rows(),
    cols = this->Error.cols();
  if(rows==0 && cols==0){
    this->Error = Eigen::MatrixXd::Zero(error2.size(),1);
    this->Error.col(0) = error2;
  }
  else if(error2.size() != rows){
      printf("SumInQuadrature Failed. A.rows != b.size\n");
      printf("A.rows=%d. A.cols=%d. b.size=%d\n",rows,cols,(int) error2.size());
      exit(1);
  }
  else{
    Eigen::MatrixXd sum = Eigen::MatrixXd::Zero(rows,cols);
    for(int i=0;i<rows;i++){
      for(int j=0;j<cols;j++){
	sum(i,j) = sqrt(this->Error(i,j)*this->Error(i,j) + error2(i)*error2(i));
	if(sum(i,j) == 0.0) sum(i,j) = 10e-5;
      }
    }
    this->Error = sum;
  }
}
void CAnalysis::ComputeTilde(){
  int rows=this->Data.rows(),
    cols=this->Data.cols();
  this->Tilde = Eigen::MatrixXd::Zero(rows,cols);
  for(int row=0;row<rows;row++){
    for(int col=0;col<cols;col++){
      this->Tilde(row,col) = (this->Data(row,col) - this->Mean(row))/this->Error(row,col);
    }
  }
}
void CAnalysis::InvertTilde(){
  int rows=this->Data.rows(),
    cols=this->Data.cols();
  this->Tilde = Eigen::MatrixXd::Zero(rows,cols);
  for(int row=0;row<rows;row++){
    for(int col=0;col<cols;col++){
      this->Data(row,col) = this->Tilde(row,col)*this->Error(row,col) + this->Mean(row);
    }
  }
}
void CAnalysis::ComputeCovariance(){
  int rows=this->Tilde.rows(),
    cols=this->Tilde.cols();
  this->Covariance = Eigen::MatrixXd::Zero(rows,rows);
  for(int row1=0;row1<rows;row1++){
    for(int row2=0;row2<rows;row2++){
      for(int col=0;col<cols;col++){
	this->Covariance(row1,row2) += this->Tilde(row1,col)*this->Tilde(row2,col);
      }
      this->Covariance(row1,row2) /= (double) cols;
    }
  }
}
void CAnalysis::EigenSolve(){
  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(this->Covariance);
  if(eigensolver.info() != Eigen::Success) abort();
  this->EigenValues = eigensolver.eigenvalues().real();
  this->EigenVectors = eigensolver.eigenvectors().real();
}
void CAnalysis::EigenSort(){
  Eigen::VectorXd eigval = this->EigenValues;
  Eigen::MatrixXd eigvec = this->EigenVectors;
  int index,
    vals = eigval.size();
  double max;
  bool tick=false;
  Eigen::MatrixXd eigsort = Eigen::MatrixXd::Zero(vals,vals),
    eigval_matrix = Eigen::MatrixXd::Zero(vals,2);
  for(int val=0;val<vals;val++){
    eigval_matrix(val,0) = 0;
    eigval_matrix(val,1) = -1;
  } 
  for(int val=0;val<vals;val++){
    index=0;
    for(int search2=0;search2<val;search2++){
      if(eigval_matrix(search2,1)==index){
	index++;
	search2=-1;
      }
    }
    max = eigval(index,0);
    for(int search1=index;search1<vals;search1++){
      tick=false;
      for(int search2=0;search2<val;search2++){
	if(eigval_matrix(search2,1)==search1){
	  tick=true;
	  break;
	}
      }
      if(tick==true){
	continue;
      }
      if(eigval(search1)>max){
	max = eigval(search1);
	index = search1;
      }
      eigval_matrix(val,1) = index;
      eigval_matrix(val,0) = max;
    }
  }
  for(int col=0;col<vals;col++){
    eigval(col) = eigval_matrix(col,0);
    index=eigval_matrix(col,1);
    for(int row=0;row<vals;row++){
      eigsort(row,col) = eigvec(row,index);
    }
  }
  eigvec = eigsort;
  this->EigenValues = eigval;
  this->EigenVectors = eigvec;
}
void CAnalysis::ComputeZ(){
  this->Z = this->Tilde.transpose()*this->EigenVectors;
}

/*
void ZerothMoment(Eigen::MatrixXd Function, double &Moment, double &MomentError){
  int points=Function.rows()-1;
  double f=0,dx=0,error=0;
  Moment = MomentError = 0;
  for(int point=0;point<points;point++){
    error = (Function(point+1,2) + Function(point,2))/2.0;
    f = (Function(point+1,1) + Function(point,1))/2.0;
    dx = Function(point+1,0) - Function(point,0);
    Moment += f*dx;
    MomentError += error*dx;
  }
}
void FirstMoment(Eigen::MatrixXd Function, double &Moment, double &MomentError){
  int points=Function.rows()-1;
  double f=0,x=0,dx=0,error=0;  //zero=0,zeroerror=0;
  Moment = MomentError = 0;
  for(int point=0;point<points;point++){
    error = (Function(point+1,2) + Function(point,2))/2.0;
    f = (Function(point+1,1) + Function(point,1))/2.0;
    x = (Function(point+1,0) + Function(point,0))/2.0;    
    dx = Function(point+1,0) - Function(point,0);
    MomentError += error*x*dx;
    Moment += f*x*dx;
  }
  //ZerothMoment(Function,zero,zeroerror);
  //  Moment /= zero;
}
void SecondMoment(Eigen::MatrixXd Function, double &Moment, double &MomentError){
  int points=Function.rows()-1;
  double f=0,x=0,dx=0,error=0,
    zero=0,zeroerror=0,
    first=0,firsterror=0,
    m2=0,s0=0,
    sp1=0,sp0=0;
  //ZerothMoment(Function,zero,zeroerror);
  FirstMoment(Function,first,firsterror);
  if(first != 0.0)
    {
      sp1 = firsterror/first;
      sp1 *= sp1;
    }
  else
    {
      sp1 = 1.0;
    }
  Moment = MomentError = 0;
  for(int point=0;point<points;point++){
    error = (Function(point+1,2) + Function(point,2))/2.0;
    f = (Function(point+1,1) + Function(point,1))/2.0;
    x = (Function(point+1,0) + Function(point,0))/2.0;    
    dx = Function(point+1,0) - Function(point,0);
    m2 = f*(x-first)*(x-first)*dx;
    Moment += m2;
    if(f != 0.0)
      {
	sp0 = error/f;
	sp0 *= sp0;
      }
    else
      {
	sp0 = 0.0;
      }
    MomentError += m2*sqrt(2.0*(sp1) + sp0);
  }
  //  Moment /= zero;
}
void MatrixMoments(std::vector<Eigen::MatrixXd> matrix, std::vector<Eigen::MatrixXd> matrixerror, Eigen::VectorXd DelX, Eigen::MatrixXd &Obs, Eigen::MatrixXd &ObsError){
  int files = matrix.size(),
    points = matrix[0].rows(),
    obs_file = 3,
    obs = obs_file*matrix.size(),
    runs = matrix[0].cols();
  double zero, first, second,
    zeroerror, firsterror, seconderror;
  Eigen::MatrixXd function = Eigen::MatrixXd::Zero(points,3);
  Obs = Eigen::MatrixXd::Zero(obs,runs);
  ObsError = Eigen::MatrixXd::Zero(obs,runs);
  for(int file=0;file<files;file++){
    for(int run=0;run<runs;run++){
      for(int point=0;point<points;point++){
	function(point,0) = DelX(point);
	function(point,1) = matrix[file](point,run);
	function(point,2) = matrixerror[file](point,run);
      }
      ZerothMoment(function,zero,zeroerror);
      FirstMoment(function,first,firsterror);
      SecondMoment(function,second,seconderror);

      Obs(file*obs_file, run) = zero;
      Obs(file*obs_file +1, run) = first;
      Obs(file*obs_file +2, run) = second;

      ObsError(file*obs_file, run) = zeroerror;
      ObsError(file*obs_file +1, run) = firsterror;
      ObsError(file*obs_file +2, run) = seconderror;
    }
  }    
}
*/
