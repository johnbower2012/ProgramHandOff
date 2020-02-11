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

int main(int argc, char* argv[])
{
  if(argc!=3)
    {
      printf("Usage: ./parameters start finish\n");
      return 1;
    }

  std::vector<std::string> 
    Name;
  std::string 
    foldername="model_output",
    filename="posterior.dat",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    infilename=foldername+"/moments_parameters.dat",
    delimiter=" ";
  //start, finish for model runs to use
  int
    start=atoi(argv[1]),
    finish=atoi(argv[2]),
    ab=4;
  Eigen::MatrixXd 
    matrix,
    range,
    mean,
    std;
  LoadParamFile(rangename,Name,range,delimiter);
  WritePosteriorParameterFiles(foldername,filename,Name,writefilename,delimiter,start,finish);
  std::string cmd = "for((i="+std::to_string(start)+";i<"+std::to_string(finish)+";i++)); do fn=$(printf '"+foldername+"/run%04d/' $i); cp -v "+foldername+"/fixed_parameters.dat $fn; done";
  system(cmd.c_str());

}
