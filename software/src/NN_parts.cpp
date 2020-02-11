#include "NN_parts.h"

/**ACTIVATIONS--
******************/

Eigen::MatrixXd CRelu::Function(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      if(Z(irow,icol) < 0.0) temp(irow,icol) = 0.0;
      else temp(irow,icol) = Z(irow,icol);
    }
  }
  return temp;
}
Eigen::MatrixXd CRelu::Derivative(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      if(Z(irow,icol) < 0.0) temp(irow,icol) = 0.0;
      else Z(irow,icol) = 1.0;
    }
  }
  return temp;
}
Eigen::MatrixXd CSigmoid::Function(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      temp(irow,icol) = 1.0/(1.0+exp(-Z(irow,icol)));
    }
  }
  return temp;
}
Eigen::MatrixXd CSigmoid::Derivative(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = this->Function(Z);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      temp(irow,icol) = temp(irow,icol)*(1.0 - temp(irow,icol));
    }
  }
  return temp;
}
Eigen::MatrixXd CTanh::Function(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(rows,cols);
  double expp, expn;
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      expp = exp(Z(irow,icol));
      expn = exp(-Z(irow,icol));
      temp(irow,icol) = (expp - expn)/(expp+expn);
    }
  }
  return temp;
}
Eigen::MatrixXd CTanh::Derivative(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = this->Function(Z);
  double tanh;
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      tanh = temp(irow,icol);
      temp(irow,icol) = 1.0 - tanh*tanh;
    }
  }
  return temp;
}
Eigen::MatrixXd CSoftMax::Function(Eigen::MatrixXd Z){
  int rows=Z.rows(),
    cols=Z.cols();
  double sum=0.0;
  Eigen::MatrixXd softMax = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    sum=0.0;
    for(int icol=0;icol<cols;icol++){
      softMax(irow,icol) = exp(Z(irow,icol));
      sum += softMax(irow,icol);
    }
    softMax.row(irow) /= sum;
  }
  return softMax;
}
Eigen::MatrixXd CSoftMax::Derivative(Eigen::MatrixXd Z){ //NEED TO CREATE // PLACEHOLDER
  return Z;
}
Eigen::MatrixXd CIdentity::Function(Eigen::MatrixXd Z){
  return Z;
}
Eigen::MatrixXd CIdentity::Derivative(Eigen::MatrixXd Z){
  int rows=Z.rows(),
    cols=Z.cols();
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Zero(rows,cols);
  for(int icol=0;icol<cols;icol++){
    for(int irow=0;irow<rows;irow++){
      Ones(irow,icol) = 1.0;
    }
  }

  return Ones;
}


/**LOSS--
******************/

void CLoss::Construct(std::string SActivation, double reg_param){
  this->Regularization_Parameter = reg_param;
  if(this->Activation != nullptr){
    delete this->Activation;
  }

  if(SActivation=="SOFTMAX"){
    this->Activation = new CSoftMax;
  }else if(SActivation=="SIGMOID"){
    this->Activation = new CSigmoid;
  }else if(SActivation=="TANH"){
    this->Activation = new CTanh;
  }else if(SActivation=="IDENTITY"){
    this->Activation = new CIdentity;
  }else if(SActivation=="RELU"){
    this->Activation = new CRelu;
  }else{
    this->Activation = nullptr;
    printf("%s is not a valid option for CActivation in CLoss.\n",SActivation.c_str());
    exit(1);
  }
}
double CL2Loss::Function(Eigen::MatrixXd ZL,Eigen::MatrixXd Y){
  int rows=ZL.rows(),
    cols=ZL.cols();
  double loss=0.0;
  Eigen::MatrixXd AL = this->Activation->Function(ZL);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      loss += (AL(i,j) - Y(i,j))*(AL(i,j) - Y(i,j));
    }
  }
  return 0.5*loss/(double) AL.rows();
}
Eigen::MatrixXd CL2Loss::Derivative(Eigen::MatrixXd ZL, Eigen::MatrixXd Y){
  int rows=ZL.rows(),
    cols=ZL.cols();
  Eigen::MatrixXd Loss = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd AL = this->Activation->Function(ZL);
  Eigen::MatrixXd act_der = this->Activation->Derivative(ZL);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      Loss(i,j) = (AL(i,j) - Y(i,j))*act_der(i,j);
    }
  }
  return Loss;
}
double CEntropyLoss::Function(Eigen::MatrixXd ZL,Eigen::MatrixXd Y){
  int rows=ZL.rows(),
    cols=ZL.cols();
  double loss=0.0;
  Eigen::MatrixXd AL = this->Activation->Function(ZL);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      loss += -(double)Y(i,j)*log(AL(i,j) + 1e-6) - (1.0 - (double)Y(i,j))*log(1.0 - AL(i,j) + 1e-6); //add 1e-6 for stability
    }
  }

  return loss/(double) rows;
}
Eigen::MatrixXd CEntropyLoss::Derivative(Eigen::MatrixXd ZL, Eigen::MatrixXd Y){ //Assumes SoftMax
  int rows=ZL.rows(),
    cols=ZL.cols();
  Eigen::MatrixXd Loss = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd AL = this->Activation->Function(ZL);
  for(int j=0;j<cols;j++){
    for(int i=0;i<rows;i++){
      Loss(i,j) = AL(i,j)-Y(i,j); 
    }
  }
  return Loss;
}

/**SOLVERS--
******************/

void CSolver::Construct(double LR, double RP){
  this->learning_rate = LR;
  this->reg_param = RP;
}
void CSolver::ConstructSGD(double Momentum, std::vector<int> Layers){
  this->momentum = Momentum;
  int layers = Layers.size();
  for(int ilay=0;ilay<layers-1;ilay++){ //Update coeffs
    int layer=Layers[ilay];
    int n_layer=Layers[ilay+1];
    this->old_w.push_back(Eigen::MatrixXd::Zero(layer,n_layer));
    this->old_b.push_back(Eigen::MatrixXd::Zero(1,n_layer));
  }
}
void CSGD::Solve(std::vector<Eigen::MatrixXd> &Weight, std::vector<Eigen::MatrixXd> &Bias, std::vector<Eigen::MatrixXd> delta_w, std::vector<Eigen::MatrixXd> delta_b){
  int layers=Weight.size();
  double OneMinusBeta = (1.0 - this->momentum)*this->learning_rate; //takes into account learning rate & introduction of new delta all in one step
  for(int ilay=0;ilay<layers;ilay++){ //Update coeffs
    old_w[ilay] = this->momentum*old_w[ilay] + OneMinusBeta*delta_w[ilay];
    old_b[ilay] = this->momentum*old_b[ilay] + OneMinusBeta*delta_b[ilay];
    Weight[ilay] -= old_w[ilay];
    Bias[ilay] -= old_b[ilay];
  }
}
void CSolver::ConstructAdam(double beta1, double beta2, std::vector<int> Layers){
  this->t = 0;
  this->epsilon = 1e-8;
  this->beta1 = beta1;
  this->beta2 = beta2;
  this->moment1_w.resize(0);
  this->moment1_b.resize(0);
  this->moment2_w.resize(0);
  this->moment2_b.resize(0);
  int layers=Layers.size();
  for(int ilay=0;ilay<layers-1;ilay++){ //Update coeffs
    int layer=Layers[ilay];
    int n_layer=Layers[ilay+1];
    this->moment1_w.push_back(Eigen::MatrixXd::Zero(layer,n_layer));
    this->moment1_b.push_back(Eigen::MatrixXd::Zero(1,n_layer));
    this->moment2_w.push_back(Eigen::MatrixXd::Zero(layer,n_layer));
    this->moment2_b.push_back(Eigen::MatrixXd::Zero(1,n_layer));
  }
}
void CAdam::Solve(std::vector<Eigen::MatrixXd> &Weight, std::vector<Eigen::MatrixXd> &Bias, std::vector<Eigen::MatrixXd> delta_w, std::vector<Eigen::MatrixXd> delta_b){
  int layers=Weight.size();
  int layer=0;
  int n_layer=0;
  double OneMinusBeta1 = 1.0 - this->beta1;
  double OneMinusBeta2 = 1.0 - this->beta2;
  double OneMinusBeta1T = 1.0 - pow(this->beta1,this->t);
  double OneMinusBeta2T = 1.0 - pow(this->beta2,this->t);
  double alpha_t = this->learning_rate*sqrt(OneMinusBeta2T)/OneMinusBeta1T;
  for(int ilay=0;ilay<layers;ilay++){ //Update coeffs
    layer=this->moment1_w[ilay].rows();
    n_layer=this->moment1_w[ilay].cols();

    //updated first-moment based on previous time step
    this->moment1_w[ilay] *= this->beta1;
    this->moment1_b[ilay] *= this->beta1;
    //updated second-moment based on previous time step
    this->moment2_w[ilay] *= this->beta2;
    this->moment2_b[ilay] *= this->beta2;

    //update first-moment based on new gradients
    this->moment1_w[ilay] += OneMinusBeta1*delta_w[ilay];
    this->moment1_b[ilay] += OneMinusBeta1*delta_b[ilay];
    //update second-moment based on new gradients
    for(int j=0;j<n_layer;j++){
      this->moment2_b[ilay](0,j) += OneMinusBeta2*delta_b[ilay](0,j)*delta_b[ilay](0,j);
      for(int i=0;i<layer;i++){
	this->moment2_w[ilay](i,j) += OneMinusBeta2*delta_w[ilay](i,j)*delta_w[ilay](i,j);
      }
    }

    /*
    //updated first-moment
    this->moment1_w[ilay] /= OneMinusBeta1T;
    this->moment1_b[ilay] /= OneMinusBeta1T;
    //updated second-moment
    this->moment2_w[ilay] /= OneMinusBeta2T;
    this->moment2_b[ilay] /= OneMinusBeta2T;
    for(int j=0;j<n_layer;j++){
      Bias[ilay](j,0) -= this->learning_rate*this->moment1_b[ilay](j,0)/(sqrt(this->moment2_b[ilay](j,0)) + this->epsilon);
      for(int i=0;i<layer;i++){
	Weight[ilay](i,j) -= this->learning_rate*this->moment1_w[ilay](i,j)/(sqrt(this->moment2_w[ilay](i,j)) + this->epsilon);
      }
    }
    */
    for(int j=0;j<n_layer;j++){
      Bias[ilay](0,j) -= alpha_t*this->moment1_b[ilay](0,j)/(sqrt(this->moment2_b[ilay](0,j)) + this->epsilon);
      for(int i=0;i<layer;i++){
	Weight[ilay](i,j) -= alpha_t*this->moment1_w[ilay](i,j)/(sqrt(this->moment2_w[ilay](i,j)) + this->epsilon);
      }
    }
  }
}
