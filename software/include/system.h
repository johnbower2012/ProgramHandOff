#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<Eigen/Dense>
#include<string>
#include<vector>
//#include "analysis.h"


bool BothAreSpaces(char lhs, char rhs);
void RemoveSpaces(std::string& str);
void LHCSampling(Eigen::MatrixXd &hypercube, int samples, Eigen::MatrixXd range);

class CSystem{
 public:
  std::string delimiter;

  CSystem();
  CSystem(std::string Delimiter);
  void setDelimiter(std::string Delimiter);

  void Mkdir(std::string folder);
  void MkdirLoop(std::string folder, int start, int finish);
  void Touch(std::string file);

  void PrintFile(std::string filename);
  void PrintFormattedFile(std::string filename);

  Eigen::MatrixXd LoadFile(std::string filename);
  std::vector<Eigen::MatrixXd> LoadFiles(std::string folder, std::string filename, int start, int finish);

  void WriteFile(std::string filename, Eigen::MatrixXd &Matrix);
  void WriteFile(std::string filename, Eigen::VectorXd &Vector);
  void WriteCSVFile(std::string filename, std::vector<std::string> header, Eigen::MatrixXd &Matrix);

/*************************************************************
   loads RangeName file to create parameter files 
      from Start(0) to Finish(1000) in Foldername/run%04d/Filename
   This function also hands a copy of the parameters
      in the Parameters matrix
   AB denotes how many different quark pairs, uu, ud, us, ss = 4
*************************************************************/

  void LoadParamFile(std::string filename, std::vector<std::string> &Distribution, std::vector<std::string> &Names, Eigen::MatrixXd &Matrix);
  void WriteParamFile(std::string fileName, std::vector<std::string> &header, Eigen::MatrixXd &file);
  void WriteParamFileLoop(std::string filename, std::string folder, int start, std::vector<std::string> &header, Eigen::MatrixXd &matrix);

};



#endif
