#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include<Eigen/Dense>
#include<chrono>
#include<random>
#include<cmath>
#include<iostream>
#include "parametermap.h"
#include "NN_parts.h"

void AddOnesColumn(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix, int column);

class CEmulator{
 public:
  virtual ~CEmulator(){}
  virtual void Construct(Eigen::MatrixXd newX, Eigen::MatrixXd newY, CParameterMap MAP) =0;
  virtual Eigen::MatrixXd Emulate(Eigen::MatrixXd testX) =0;
};

/**GAUSSIAN PROCESS--
******************/

class CGaussianProcess : public CEmulator{
 public:
  int trainPoints;
  int paramCount;
  int obsCount;
  int hyperparamCount; 
  
  double Epsilon; //Numerical stability factor for matrix inversion, 1e-8 by default
  double SigmaF; //Variance of func to be emulated
  double CharacLength; //Length overwhich func varies
  double SigmaNoise; //Noise of func
  Eigen::VectorXd Noise;

  Eigen::MatrixXd X;
  Eigen::MatrixXd Y;
  Eigen::MatrixXd Hyperparameters;
  Eigen::MatrixXd Beta;

  std::vector<Eigen::MatrixXd> Kernel;
  std::vector<Eigen::MatrixXd> KernelInv;
  std::vector<Eigen::MatrixXd> HMatrix;
  std::vector<Eigen::MatrixXd> HMatrix_KernelInv;
  std::vector<Eigen::MatrixXd> KernelInv_Y;
  std::vector<Eigen::MatrixXd> betaMatrix;

  CGaussianProcess();
  ~CGaussianProcess(){}
  CGaussianProcess(Eigen::MatrixXd newX, Eigen::MatrixXd newY, CParameterMap MAP);
  CGaussianProcess(Eigen::MatrixXd newX, Eigen::MatrixXd newY, std::string filename);
  void Construct(Eigen::MatrixXd newX, Eigen::MatrixXd newY, CParameterMap MAP);
  void Construct(Eigen::MatrixXd newX, Eigen::MatrixXd newY, std::string filename);

  void ConstructHyperparameters();  
  void ConstructBeta();
  void RegressionLinearFunction();
  
  Eigen::MatrixXd KernelFunction(Eigen::MatrixXd A, Eigen::MatrixXd B, int obsIndex);
  Eigen::MatrixXd RegressionLinearFunction(Eigen::MatrixXd testX, int obsIndex);

  Eigen::MatrixXd Emulate(Eigen::MatrixXd testX);
  Eigen::MatrixXd Emulate_NR(Eigen::MatrixXd testX);
};

/**NEURAL NETWORK--
******************/

class CNeuralNet : public CEmulator{
public:
  CParameterMap MAP;
  int layers;
  double learning_rate;
  double regular_param;
  double momentum;
  double beta1;
  double beta2;
  int mini_batch;
  int epochs;

  std::string SActivation;
  std::string SLoss;
  std::string SFinalActivation;
  std::string SSolver;

  Eigen::MatrixXd X;
  Eigen::MatrixXd Y;

  CActivation *Activation;
  CLoss *Loss;
  CSolver *Solver;

  std::vector<int> Layers;
  
  std::vector<Eigen::MatrixXd> Weight;
  std::vector<Eigen::MatrixXd> Bias;

  std::vector<Eigen::MatrixXd> activations;
  std::vector<Eigen::MatrixXd> zs;

  std::vector<Eigen::MatrixXd> delta;
  std::vector<Eigen::MatrixXd> delta_w;
  std::vector<Eigen::MatrixXd> delta_b;

  CNeuralNet(){
    this->Activation = nullptr;
    this->Loss = nullptr;
    this->Solver = nullptr;
  };
  ~CNeuralNet(){
    if(this->Activation != nullptr) delete this->Activation;
    if(this->Loss != nullptr) delete this->Loss;
    if(this->Solver != nullptr) delete this->Solver;
  };
  
  Eigen::MatrixXd Emulate(Eigen::MatrixXd X); //virtual from Emulator
  void Construct(Eigen::MatrixXd X_train, Eigen::MatrixXd Y_train, CParameterMap Map); //virtual from Emulator

  void FeedForward(Eigen::MatrixXd x);
  void BackPropagation(Eigen::MatrixXd x, Eigen::MatrixXd y);
  void Train(int Epochs);

  Eigen::MatrixXd Max(Eigen::MatrixXd Z);

  double Accuracy();
  double Accuracy(Eigen::MatrixXd x, Eigen::MatrixXd y);
};



#endif
