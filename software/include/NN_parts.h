#ifndef __NN_PARTS_H__
#define __NN_PARTS_H__

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>

/**ACTIVATIONS--
******************/

class CActivation{
 public:
  CActivation(){};
  virtual ~CActivation(){};

  virtual Eigen::MatrixXd Function(Eigen::MatrixXd Z) =0;
  virtual Eigen::MatrixXd Derivative(Eigen::MatrixXd Z) =0;
};
class CRelu : public CActivation{
 public:
  CRelu(){};
  ~CRelu(){};

  Eigen::MatrixXd Function(Eigen::MatrixXd Z);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd Z);
};
class CSigmoid : public CActivation{
 public:
  CSigmoid(){};
  ~CSigmoid(){};

  Eigen::MatrixXd Function(Eigen::MatrixXd Z);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd Z);
};
class CTanh : public CActivation{
 public:
  CTanh(){};
  ~CTanh(){};

  Eigen::MatrixXd Function(Eigen::MatrixXd Z);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd Z);
};
class CSoftMax : public CActivation{
 public:
  CSoftMax(){};
  ~CSoftMax(){};

  Eigen::MatrixXd Function(Eigen::MatrixXd Z);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd Z); //NOT DEFINED -- PLACEHOLDER
};
class CIdentity : public CActivation{
 public:
  CIdentity(){};
  ~CIdentity(){};

  Eigen::MatrixXd Function(Eigen::MatrixXd Z);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd Z);
};

/**LOSS--
******************/

class CLoss{
 public:
  double Regularization_Parameter;
  CActivation *Activation;

  CLoss(){
    this->Activation = nullptr;
    this->Regularization_Parameter = 0.0;
  };
  virtual ~CLoss(){
    if(this->Activation != nullptr) delete this->Activation;
  };
  
  void Construct(std::string SActivation, double reg_param);

  virtual double Function(Eigen::MatrixXd ZL, Eigen::MatrixXd Y) =0;
  virtual Eigen::MatrixXd Derivative(Eigen::MatrixXd ZL, Eigen::MatrixXd Y) =0;
};
class CL2Loss : public CLoss{
 public:
  CL2Loss(){};
  ~CL2Loss(){};

  double Function(Eigen::MatrixXd ZL, Eigen::MatrixXd Y);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd ZL, Eigen::MatrixXd Y);
};
class CEntropyLoss : public CLoss{ //Assumes SoftMax
 public:
  CEntropyLoss(){};
  ~CEntropyLoss(){};

  double Function(Eigen::MatrixXd ZL, Eigen::MatrixXd Y);
  Eigen::MatrixXd Derivative(Eigen::MatrixXd ZL, Eigen::MatrixXd Y);//Assumes SoftMax
};

/**SOLVERS--
******************/

class CSolver{
 public:
  CSolver(){};
  virtual ~CSolver(){};
  double learning_rate;
  double reg_param;
  int t;
  double momentum;
  double beta1;
  double beta2;
  double epsilon;
  std::vector<Eigen::MatrixXd> old_w;
  std::vector<Eigen::MatrixXd> old_b;
  std::vector<Eigen::MatrixXd> moment1_w;
  std::vector<Eigen::MatrixXd> moment1_b;
  std::vector<Eigen::MatrixXd> moment2_w;
  std::vector<Eigen::MatrixXd> moment2_b;

  void Construct(double LR, double RP);
  void ConstructAdam(double beta1, double beta2, std::vector<int> Layers);
  void ConstructSGD(double Momentum, std::vector<int> Layers);
  virtual void Solve(std::vector<Eigen::MatrixXd> &Weight, std::vector<Eigen::MatrixXd> &Bias, std::vector<Eigen::MatrixXd> delta_w, std::vector<Eigen::MatrixXd> delta_b) =0;
};
class CSGD : public CSolver{
 public:
  CSGD(){};
  ~CSGD(){};
  
  void Solve(std::vector<Eigen::MatrixXd> &Weight, std::vector<Eigen::MatrixXd> &Bias, std::vector<Eigen::MatrixXd> delta_w, std::vector<Eigen::MatrixXd> delta_b);
};
class CAdam : public CSolver{
 public:
  CAdam(){};
  ~CAdam(){};

  void Solve(std::vector<Eigen::MatrixXd> &Weight, std::vector<Eigen::MatrixXd> &Bias, std::vector<Eigen::MatrixXd> delta_w, std::vector<Eigen::MatrixXd> delta_b);
};

#endif 
