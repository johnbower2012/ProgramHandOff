#ifndef __BALANCEMODEL_H__
#define __BALANCEMODEL_H__

#include<string>
#include<Eigen/Dense>
#include "parametermap.h"
#include "system.h"
#include "analysis.h"
#include "emulator.h"
#include "mcmc.h"
#include "coshfunc.h"

class CBalanceModel:public CSystem{
 public:
  CParameterMap MAP;
  std::vector<CAnalysis> PCA;
  CAnalysis PCA_RD;
  CEmulator* Emulator;
  CMCMC MCMC;

  int START;
  int FINISH;
  int QParameters; //Number of Quark Parameters
  int QuarkPairs; //Number of quark pairs
  int BFCount; //Number of diff balance functions
  int BFSamples; //Number of model samples
  int PrinComp; //Number of principal components to keep per BF
  int NSAMPLES; //Number of MCMC Samples
  int testZ; //Line from ModelZ to target, -1 = ExpZ

  bool CUT_MODEL;
  bool CUT_EXP;
  bool CUT_WIDTH;
  bool COMMON_WIDTH;
  bool CUT_G0;
  bool FIX_G0;
  bool SCALE_X;

  //DATA PROCESSING INFO
  double PCA_CUTOFF;
  bool MEAN_ERROR;
  bool USE_MEAN_AS_ERROR;
  double MEAN_AS_ERROR;
  bool REDUCE_DIM;
  bool REDUCE_DIM_SEPARATELY;

  //GAUSSIAN PROCESS INFO
  double EPSILON;
  double SIGMA_F;
  double CHARAC_LENGTH;
  double SIGMA_NOISE;

  //NEURAL NETWORK INFO
  double LEARNING_RATE;
  double REGULAR_PARAM;
  int MINI_BATCH_SIZE;
  double MOMENTUM;
  double BETA1;
  double BETA2;
  double LAYERS;
  std::string SACTIVATION;
  std::string SLOSS;
  std::string SSOLVER;

  //MCMC INFO
  int NTRACE;
  double MCMC_WIDTH;
  double MCMC_MIN;
  double MCMC_MAX;
  int MCMC_POST;
  double MCMC_PERCENT_ERROR;

  //GAB FUNCTIONS INFO
  int GRID;
  double ETA_MAX;

  //EMULATOR & MCMC CHOICES
  std::string EMULATOR_CHOICE;
  std::string MCMC_CHOICE;

  //GAB CHI
  Eigen::VectorXd CHI;

  //MODEL DATA INFO
  std::string MODEL_FOLDER;
  std::vector<std::string> MODEL_NAMES;

  //EXP DATA INFO
  std::string EXP_FOLDER;
  std::vector<std::string> EXP_NAMES;

  //OUTPUT AND FILE INFO
  std::string OUTPUT_FOLDER;
  std::string RUN_FOLDER;
  std::string TRAINPOINTS_FILE;
  std::string EXPZ_FILE;
  std::string PARAMS_FILE;
  std::string MCMCTRACE_FILE;
  std::string MCMCHISTORY_FILE;
  std::string GAB_FILE;
  std::string FULLG_FILE;
  std::string MINMAX_FILE;
  std::string POST_EXT;
  std::string FILE_EXT;
  std::string CSV_EXT;

  //BF VECTOR & MATRIX NAMES
  Eigen::VectorXd ModelRapidity;
  std::vector<Eigen::MatrixXd> Model;
  std::vector<Eigen::MatrixXd> ModelError;
  Eigen::VectorXd ExperimentRapidity;
  std::vector<Eigen::MatrixXd> Experiment;
  std::vector<Eigen::VectorXd> ExperimentalError;

  //ERROR VECTOR FOR USE IN MCMC SEARCH
  Eigen::MatrixXd McmcError;

  //MODEL PARAM MATRIX NAMES
  Eigen::MatrixXd WidthParameters;
  Eigen::MatrixXd Parameters;
  Eigen::MatrixXd UnscaledParameters;  
  Eigen::MatrixXd MinMax;

  //Z MATRIX NAMES
  Eigen::MatrixXd ModelZ;
  Eigen::MatrixXd ExpZ;

  //STORE MCMCHISTORY OF STEPS
  Eigen::MatrixXd MCMCHistory;

  //CLASS TO MANIPULATE BF FUNCTION INFO
  CBalanceModel(std::string filename, int testZ); 

  //LOAD AND MANIPULATE
  void LoadData();
  void CutData();
  void WriteBF();

  //PROCESS
  void ScaleParameters();
  void ReduceDimensionsSeparately();
  void ReduceDimensions();

  //CALCULATE Zs
  void WriteEigen();
  void WriteZ();
  void WriteReconBF();

  //CREATE & DELETE EMULATOR CHOICE
  void CreateEmulator();
  void DeleteEmulator();

  //CREATE AND RUN MCMC CHOICE
  void CreateMCMC();
  void RunMCMC();

  //WRITE MCMC & GAB FILES
  void WriteMCMC();
  void WriteCoshFunctions();
  
  //SCALE GAB BY CHI
  void ScaleByChi(Eigen::MatrixXd &Functions);

  //EXTRACT POSTERIOR (20) FROM MCMCTRACE
  Eigen::MatrixXd ExtractPosterior(Eigen::MatrixXd MCMC);

  //UNSCALE PARAM FROM [0,1]
  Eigen::MatrixXd UnscaleParameters(Eigen::MatrixXd ScaledParameters);
};


#endif
