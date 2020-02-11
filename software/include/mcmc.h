#ifndef __MCMC_H__
#define __MCMC_H__

#include<Eigen/Dense>
#include<random>
#include<chrono>
#include<math.h>
#include<string>
#include "emulator.h"

class CRandom{
 public:
  unsigned seed;
  std::default_random_engine generator;
  std::normal_distribution<double> normal_dist;
  std::uniform_real_distribution<double> uniform_dist;

  CRandom();
  CRandom(unsigned Seed);

  void setSeed(unsigned Seed);
  void setSeedClock();

  double normal();
  double uniform();
};

class CLikelihood{
 public:
  CLikelihood(){};
  virtual ~CLikelihood(){};
  virtual double getLikelihood(Eigen::MatrixXd Z1, Eigen::MatrixXd Z2, Eigen::MatrixXd Error) =0;
};

class CLHGaussian : public CLikelihood{
  ~CLHGaussian(){};
  double getLikelihood(Eigen::MatrixXd Z1, Eigen::MatrixXd Z2, Eigen::MatrixXd Error);  
};

class CLHLorentzian : public CLikelihood{
  ~CLHLorentzian(){};
  double getLikelihood(Eigen::MatrixXd Z1, Eigen::MatrixXd Z2, Eigen::MatrixXd Error);  
};
class CLHCosh : public CLikelihood{
  ~CLHCosh(){};
  double getLikelihood(Eigen::MatrixXd Z1, Eigen::MatrixXd Z2, Eigen::MatrixXd Error);  
};

class CMCMC{
 public:
  double maxLogLikelihood;
  double Likelihood;
  CRandom random;

  int paramCount;
  int NSamples;

  std::string SLikelihood;
  CLikelihood *FLikelihood;

  Eigen::MatrixXd Range;
  Eigen::VectorXd Widths;

  Eigen::MatrixXd Position;
  Eigen::MatrixXd testPosition;

  CMCMC();
  CMCMC(Eigen::MatrixXd newRange, Eigen::VectorXd newWidths);
  CMCMC(std::string filename);
  ~CMCMC(){
    if(this->FLikelihood != nullptr)
      delete this->FLikelihood;
    this->FLikelihood = nullptr;
  }
  void Construct(Eigen::MatrixXd newRange, Eigen::VectorXd newWidths, std::string slikelihood);
  void setSeed(unsigned Seed);
  void setPosition();
  void setPosition(Eigen::MatrixXd newPosition);
  void setRange(Eigen::MatrixXd newRange);
  void setWidths(Eigen::VectorXd newWidths);

  void step();
  //void step(CCosh dist, Eigen::MatrixXd MinMax);
  //double getLikelihood(Eigen::MatrixXd Z1, Eigen::MatrixXd Z2, Eigen::MatrixXd Error);
  bool decide(double likelihood);
  //double getLogLikelihoodGaussian(Eigen::MatrixXd Z);
  //bool decideGaussian(Eigen::MatrixXd Z);
  //double getLogLikelihoodLorentzian(Eigen::MatrixXd Z);
  //bool decideLorentzian(Eigen::MatrixXd Z);

  Eigen::MatrixXd Run(CEmulator *Emulator, Eigen::MatrixXd Target, Eigen::MatrixXd Error, int NSamp);
};

#endif
