#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include<Eigen/Dense>
#include<vector>
#include<iostream>
#include "system.h"

class CAnalysis{
 public:
  Eigen::MatrixXd Data;
  Eigen::MatrixXd Tilde;
  Eigen::MatrixXd Z;

  Eigen::VectorXd Mean;
  Eigen::MatrixXd Error;

  Eigen::MatrixXd Covariance;

  Eigen::VectorXd EigenValues;
  Eigen::MatrixXd EigenVectors;  

  CAnalysis();

  void ComputeMean();
  void SumErrorInQuadrature(Eigen::MatrixXd Error2);
  void SumErrorInQuadrature(Eigen::VectorXd error2);
  void ComputeTilde();
  void InvertTilde();
  void ComputeCovariance();
  void EigenSolve();
  void EigenSort();
  void ComputeZ();
};

void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

void AddOnesRow(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix);
void AddOnesColumn(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix);

void ScaleMatrixColumns(Eigen::MatrixXd Matrix,Eigen::VectorXd &Mean, Eigen::VectorXd &Std, Eigen::MatrixXd &Scaled);
void ScaleMatrixRows(Eigen::MatrixXd Matrix,Eigen::VectorXd &Mean, Eigen::VectorXd &Stdd, Eigen::MatrixXd &Scaled);
void ScaleMatrixColumnsUniform(Eigen::MatrixXd Matrix,Eigen::MatrixXd &MinMax, Eigen::MatrixXd &Scaled);
void ScaleMatrixRowsUniform(Eigen::MatrixXd Matrix,Eigen::MatrixXd &MinMax, Eigen::MatrixXd &Scaled);

void AverageColumns(Eigen::VectorXd &average, Eigen::MatrixXd matrix);
void AverageRows(Eigen::VectorXd &average, Eigen::MatrixXd matrix);

void SumInQuadrature(Eigen::MatrixXd &sum, Eigen::MatrixXd A, Eigen::MatrixXd B);
void SumInQuadrature(Eigen::MatrixXd &sum, Eigen::MatrixXd A, Eigen::VectorXd b);
void SumInQuadrature(Eigen::VectorXd &sum, Eigen::VectorXd a, Eigen::VectorXd b);
/*

void ZerothMoment(Eigen::MatrixXd Function, double &Moment, double &MomentError);
void FirstMoment(Eigen::MatrixXd Function, double &Moment, double &MomentError);
void SecondMoment(Eigen::MatrixXd Function, double &Moment, double &MomentError);
void MatrixMoments(std::vector<Eigen::MatrixXd> Matrix, Eigen::VectorXd DelX, Eigen::MatrixXd &ObsX, Eigen::MatrixXd &ObsError);

*/

#endif
