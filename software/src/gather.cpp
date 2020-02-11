#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<stdlib.h>
#include "system.h"
#include "analysis.h"
#include "coshfunc.h"
#include "emulator.h"
#include "mcmc.h"

int main(int argc, char* argv[]){
  //Check that usage is proper: enter the model output to use
  if(argc != 3)
    {
      printf("Usage: ./main start end\n");
      return 1;
    }
  std::string 
    foldername = "model_output",
    delimiter=" ";
  int
    start=atoi(argv[1]),
    finish=atoi(argv[2]),
    ab=4;
  Eigen::MatrixXd Parameters;

  /*************************************
   *   Load Data for analysis         *
  *************************************/
  //Model BF filenames
  std::vector<std::string> modelfilenames={"I211_J211.dat",
					   "I2212_J2212.dat",
					   "I321_J2212.dat",
					   "I321_J321.dat"};
  //EXP quark filenames
  std::vector<std::string> expfilenames={"star_pipi.dat",
					 "star_ppbar.dat",
					 "star_pK.dat",
					 "star_KK.dat"};
  //Number of model && exp files
  //data column to be selected from model/exp files
  std::vector<Eigen::MatrixXd> 
    ModelMatrix,
    ExpMatrix;
  Eigen::VectorXd
    modeldy,
    expdy;

  LoadMEDataFiles(modelfilenames, expfilenames,
		     ModelMatrix, ExpMatrix,
		     modeldy, expdy,
		     foldername, delimiter,
		     start, finish);
  Eigen::MatrixXd holder;
  for(int i=0;i<4;i++){
    AddOnesColumn(ModelMatrix[i],holder);
    holder.col(0) = modeldy;
    WriteFile(modelfilenames[i],holder,delimiter);
  }

  return 0;
}
