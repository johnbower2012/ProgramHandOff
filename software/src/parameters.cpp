#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<stdlib.h>
#include "parametermap.h"
#include "system.h"
#include "coshfunc.h"

int main(int argc, char* argv[])
{
  if(argc!=3)
    {
      printf("Usage: ./parameters start finish\n");
      return 1;
    }

  CParameterMap MAP;
  CSystem SYSTEM;
  CDistCosh DIST;

  MAP.ReadParsFromFile("write_parameters.dat");

  std::string 
    model_foldername=MAP.getS("MODEL_FOLDER","model_output"),
    param_filename=MAP.getS("PARAM_FILE","parameters.dat"),
    fixed_filename=MAP.getS("PARAM_FIXED_FILE","fixed_parameters.dat"),
    prior_filename=MAP.getS("PARAM_PRIOR_FILE","parameter_priors.dat"),
    lhc_filename=MAP.getS("PARAMS_FILE","lhc_parameters.dat"),
    gab_filename=MAP.getS("GAB_FILE","gabfunctions.dat"),
    delimiter=" ";
  //start, finish for model runs to use
  int
    start=atoi(argv[1]),
    finish=atoi(argv[2]),
    ab=MAP.getI("QUARK_PAIRS",4),
    GRID=MAP.getI("GRID",100);
  Eigen::VectorXd CHI(ab);
  CHI(0)=MAP.getD("CHI0",7.2728e+02);
  CHI(1)=MAP.getD("CHI1",-2.2654e+02);
  CHI(2)=MAP.getD("CHI2",-8.3627e+01);
  CHI(3)=MAP.getD("CHI3",3.0778e+02);
  double
    ETA_MAX=MAP.getD("ETA_MAX",8);
  bool
    FIX_G0=MAP.getB("FIX_G0",true),
    COMMON_WIDTH=MAP.getB("COMMON_WIDTH",false);

  Eigen::MatrixXd Parameters;
  Eigen::MatrixXd Priors;
  Eigen::MatrixXd GAB;

  std::vector<std::string> Distributions, Names;

  SYSTEM.LoadParamFile(model_foldername+"/"+prior_filename, Distributions, Names, Priors);
  LHCSampling(Parameters,finish-start,Priors);
  if(COMMON_WIDTH){
    int ParamsPerPair=Parameters.cols()/ab;
    for(int pair=1;pair<ab;pair++){
      Parameters.col(pair*ParamsPerPair) = Parameters.col(0);
    }
  }
  if(FIX_G0){
    Parameters = DIST.FixG0(ab,Parameters);
  }
  SYSTEM.WriteParamFileLoop(param_filename,model_foldername,start,Names,Parameters);
  SYSTEM.WriteFile(model_foldername+"/"+lhc_filename,Parameters);
  GAB=DIST.FunctionSet(GRID,ETA_MAX,finish-start,ab,Parameters.cols()/ab-2,Parameters);
  for(int samples=0;samples<finish-start;samples++){
    for(int pair=0;pair<ab;pair++){
      GAB.col(samples*ab+pair+1) *= CHI(pair);
    }
  }
  SYSTEM.WriteFile(model_foldername+"/"+gab_filename,GAB);

  std::string cmd = "for((i="+std::to_string(start)+";i<"+std::to_string(finish)+";i++)); do fn=$(printf '"+model_foldername+"/run%04d/' $i); cp -v "+model_foldername+"/"+fixed_filename+" $fn; done";
  system(cmd.c_str());
}
