#ifndef __COSHFUNC_H__
#define __COSHFUNC_H__

#include<iostream>
#include<cmath>
#include<fstream>
#include<Eigen/Dense>
#include<random>
#include <cstdlib>
#include <cstdio>
#include <complex>

#define PI 4.0*atan(1.0)

class CDistribution{
 public:
  CDistribution();

  virtual void GenX(int ab,double &x,double &weight);
};

class CCosh{
 public:
  int nmax;
  Eigen::MatrixXd A;
  Eigen::VectorXd I;
  Eigen::VectorXd Z;

  CCosh();
  CCosh(int nmaxset);

  void Set_NMax(int nmaxset);
  double GetG(int n,double x);
  double GetRandomX();
  void TestOverlap();
  void CalcOverlap(int n,int nprime);
  void CalcOverlapNumerical(int n,int nprime);
  void Integrate(double x);
  
  int seed;
  std::mt19937 rng;
  std::uniform_real_distribution<double> uniform_dist;
};

class CDistCosh : public CDistribution{
 public:
  //Function set
  CCosh coshcalc; 

  //widths & coeffs of invcosh functions
  Eigen::VectorXd width;
  Eigen::MatrixXd g;

  CDistCosh();
  CDistCosh(Eigen::VectorXd Width, Eigen::MatrixXd G);

  void Set_WG(Eigen::VectorXd Width, Eigen::MatrixXd G);
  void GenX(int ab,double &x,double &weight);
  Eigen::MatrixXd GenG0(int ab, Eigen::MatrixXd G);
  Eigen::MatrixXd FixG0(int ab, Eigen::MatrixXd G);
  Eigen::MatrixXd FunctionSet(int grid, double x, int samples, int ab, int nmax, Eigen::MatrixXd G);
  double Integrate(int ab);
  double IntegrateNumerical(int ab);
};


#endif
