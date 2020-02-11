#include "balancemodel.h"

CBalanceModel::CBalanceModel(std::string filename,int testZ){
  this->MAP.ReadParsFromFile(filename); //Load all parameters from file
  this->testZ = testZ; //Which model to reproduce; -1 is EXP

  this->LoadData(); //Load Model & Exp Files
  this->CutData(); //Cut out 0.05 rapidity, match exp rapidity range to model rap range
  this->WriteBF(); //Write out BF,BF_error

  if(this->SCALE_X){
    this->ScaleParameters(); //Scale to [0,1]
  }
  if(REDUCE_DIM_SEPARATELY){
    this->ReduceDimensionsSeparately(); //Dimensional reduction via PCA
  }else{
    this->ReduceDimensions(); //Dimensional reduction via PCA
  }
  this->WriteZ(); //Write out ModelZ and ExpZ

  this->CreateEmulator(); //Create Emulator Choice
  this->CreateMCMC(); //Create MCMC Choice
  this->RunMCMC(); //Run MCMC

  this->WriteMCMC(); //Write out MCMC trace & history
  this->WriteCoshFunctions(); //Write GAB Prior and Posterior

  this->DeleteEmulator(); //Delete Emulator memory allocation
}

///////////////////////////////////////////////////////////////////////////////////////

void CBalanceModel::LoadData(){
  CSystem SYSTEM;
  char buffer[30];
  string name;

  this->START=this->MAP.getI("START",0);
  this->FINISH=this->MAP.getI("FINISH",500);
  this->BFSamples=FINISH-START;
  this->PrinComp=this->MAP.getI("PRIN_COMP",4);
  this->QuarkPairs=this->MAP.getI("QUARK_PAIRS",4);
  this->NSAMPLES=this->MAP.getI("NSAMPLES",100000);
  this->NTRACE=this->MAP.getI("NTRACE",(int)(this->NSAMPLES*0.1));

  this->EMULATOR_CHOICE=this->MAP.getS("EMULATOR_CHOICE","GAUSSIAN_PROCESS");
  this->MCMC_CHOICE=this->MAP.getS("MCMC_CHOICE","GAUSSIAN");

  this->MODEL_FOLDER=this->MAP.getS("MODEL_FOLDER","model_output");
  this->EXP_FOLDER=this->MAP.getS("EXP_FOLDER","model_output");
  this->OUTPUT_FOLDER=this->MAP.getS("OUTPUT_FOLDER","stat_output");
  this->RUN_FOLDER=this->MAP.getS("RUN_FOLDER","run");
  if(this->testZ!=-1){
    sprintf(buffer,"%s%04d",this->RUN_FOLDER.c_str(),this->testZ);
    this->RUN_FOLDER = buffer;
    SYSTEM.Mkdir(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER);
  }else{
    sprintf(buffer,"%s_exp",this->RUN_FOLDER.c_str());
    this->RUN_FOLDER = buffer;
    SYSTEM.Mkdir(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER);
  }

  this->TRAINPOINTS_FILE=this->MAP.getS("TRAINPOINTS_FILE","trainplot");
  this->EXPZ_FILE=this->MAP.getS("EXPZ_FILE","expz");
  this->PARAMS_FILE=this->MAP.getS("PARAMS_FILE","parameters_lhc.dat");
  this->MCMCTRACE_FILE = this->MAP.getS("MCMCTRACE_FILE","mcmctrace");
  this->MCMCHISTORY_FILE = this->MAP.getS("MCMCHISTORY_FILE","mcmchistory");
  this->MINMAX_FILE = this->MAP.getS("MINMAX_FILE","minmax");
  this->GAB_FILE = this->MAP.getS("GAB_FILE","gabfunctions");
  this->FULLG_FILE = this->MAP.getS("FULLG_FILE","fullg");
  this->POST_EXT = this->MAP.getS("POST_EXT","posterior");
  this->FILE_EXT = this->MAP.getS("FILE_EXT","dat");
  this->CSV_EXT = this->MAP.getS("CSV_EXT","csv");

  this->CUT_MODEL=this->MAP.getB("CUT_MODEL",true);
  this->CUT_EXP=this->MAP.getB("CUT_EXP",true);
  this->COMMON_WIDTH=this->MAP.getB("COMMON_WIDTH",false);
  this->FIX_G0=this->MAP.getB("FIX_G0",true);
  this->SCALE_X = this->MAP.getB("SCALE_X",true);

  this->PCA_CUTOFF=MAP.getD("PCA_CUTOFF",0.5);
  this->MEAN_ERROR=MAP.getB("MEAN_ERROR",false);
  this->USE_MEAN_AS_ERROR=MAP.getB("USE_MEAN_AS_ERROR",false);
  this->MEAN_AS_ERROR=MAP.getD("MEAN_AS_ERROR",0.05);
  this->REDUCE_DIM=MAP.getB("REDUCE_DIM",true);
  this->REDUCE_DIM_SEPARATELY=MAP.getB("REDUCE_DIM_SEPARATELY",false);

  //GAUSSIAN PROCESS INFO
  this->EPSILON = MAP.getD("EPSILON",1e-8);
  this->SIGMA_F = MAP.getD("SIGMA_F",0.5);
  this->CHARAC_LENGTH = MAP.getD("CHARAC_LENGTH",0.45);
  this->SIGMA_NOISE = MAP.getD("SIGMA_NOISE",0.05);

  //NEURAL NETWORK INFO
  this->LEARNING_RATE = this->MAP.getD("LEARNING_RATE",0.001);
  this->REGULAR_PARAM = this->MAP.getD("REGULAR_PARAM",0.0);
  this->MINI_BATCH_SIZE = this->MAP.getD("MINI_BATCH_SIZE",10);
  this->MOMENTUM = this->MAP.getD("MOMENTUM",0.9);
  this->BETA1 = this->MAP.getD("BETA1",0.9);
  this->BETA2 = this->MAP.getD("BETA2",0.999);
  this->LAYERS = this->MAP.getD("LAYERS",1);
  this->SACTIVATION = this->MAP.getS("ACTIVATION","RELU");
  this->SLOSS = this->MAP.getS("LOSS","ENTROPY");
  this->SSOLVER = this->MAP.getS("SOLVER","SGD");

  this->MCMC_WIDTH=this->MAP.getD("MCMC_WIDTH",0.005);
  this->MCMC_MIN=this->MAP.getD("MCMC_MIN",0.0);
  this->MCMC_MAX=this->MAP.getD("MCMC_MAX",1.0);
  this->MCMC_POST=this->MAP.getI("MCMC_POSTERIOR",20);
  this->MCMC_PERCENT_ERROR=this->MAP.getD("MCMC_PERCENT_ERROR",0.06);

  this->GRID=this->MAP.getI("GAB_FUNC_GRID",100);
  this->ETA_MAX=this->MAP.getD("GAB_FUNC_CUTOFF",8);

  this->CHI = Eigen::VectorXd::Zero(4);
  this->CHI(0)=this->MAP.getD("CHI0",7.2728e+02);
  this->CHI(1)=this->MAP.getD("CHI1",-2.2654e+02);
  this->CHI(2)=this->MAP.getD("CHI2",-8.3627e+01);
  this->CHI(3)=this->MAP.getD("CHI3",3.0778e+02);

  Eigen::MatrixXd params=SYSTEM.LoadFile(MODEL_FOLDER+"/"+PARAMS_FILE);
  this->QParameters=params.cols();
  this->Parameters=params.block(START,0,this->BFSamples,this->QParameters);

  int i=0;
  sprintf(buffer,"MODEL_NAME_%i",i);
  name=buffer;
  while(this->MAP.getS(name,"NULL")!="NULL"){
    this->MODEL_NAMES.push_back(this->MAP.getS(name,"NULL"));
    i++;  
    sprintf(buffer,"MODEL_NAME_%i",i);
    name=buffer;
  }
  if(i==0){
    printf("Usage: Enter at least one MODEL_NAME_X, e.g. MODEL_NAME_0.\n");
    exit(1);
  }
   
  i=0;
  sprintf(buffer,"EXP_NAME_%i",i);
  name=buffer;
  while(this->MAP.getS(name,"NULL")!="NULL"){
    this->EXP_NAMES.push_back(this->MAP.getS(name,"NULL"));
    i++;  
    sprintf(buffer,"EXP_NAME_%i",i);
    name=buffer;
  }
  if(i==0){
    printf("Usage: Enter at least one EXP_NAME_X, e.g. EXP_NAME_0.\n");
    exit(1);
  }

  int MODEL_COUNT=this->MODEL_NAMES.size();
  int EXP_COUNT=this->EXP_NAMES.size();
  if(MODEL_COUNT!=EXP_COUNT){
    printf("Usage: Enter equal model and experimental balance functions.\n");
    exit(1);
  }else{
    this->BFCount=MODEL_COUNT;
  }
  
  this->Model=std::vector<Eigen::MatrixXd>(MODEL_COUNT);
  this->ModelError=std::vector<Eigen::MatrixXd>(MODEL_COUNT);

  for(int i=0;i<MODEL_COUNT;i++){
    printf("Loading Model output file %s/%s from run%04d to run%04d...\n",this->MODEL_FOLDER.c_str(),this->MODEL_NAMES[i].c_str(),this->START,this->FINISH-1); 
    std::vector<Eigen::MatrixXd> matrix = SYSTEM.LoadFiles(this->MODEL_FOLDER,this->MODEL_NAMES[i],this->START,this->FINISH);
    int YN = matrix[i].rows();
    this->Model[i] = Eigen::MatrixXd(YN,BFSamples);
    this->ModelError[i] = Eigen::MatrixXd(YN,BFSamples);
    this->ModelRapidity = matrix[0].col(0);
    for(int run=0;run<BFSamples;run++){
      this->Model[i].col(run) = matrix[run].col(1);
      this->ModelError[i].col(run) = matrix[run].col(2);
    }
  }

  this->Experiment=std::vector<Eigen::MatrixXd>(MODEL_COUNT);
  this->ExperimentalError=std::vector<Eigen::VectorXd>(MODEL_COUNT);
  for(int i=0;i<EXP_COUNT;i++){
    printf("Loading Data file %s/%s...\n",this->EXP_FOLDER.c_str(),this->EXP_NAMES[i].c_str()); 
    Eigen::MatrixXd matrix = SYSTEM.LoadFile(EXP_FOLDER+"/"+EXP_NAMES[i]);
    this->ExperimentRapidity = matrix.col(0);
    this->Experiment[i] = matrix.col(1);
    this->ExperimentalError[i] = matrix.col(2);
  }
}
void CBalanceModel::CutData(){
  if(this->CUT_MODEL){ //Remove 0.05 rapidity from Model
    Eigen::MatrixXd MR = Eigen::MatrixXd::Zero(this->ModelRapidity.size(),1);
    MR.col(0) = this->ModelRapidity;
    RemoveRow(MR,0);
    this->ModelRapidity = MR.col(0);
    for(int i=0;i<this->BFCount;i++){
      RemoveRow(this->Model[i],0);
      RemoveRow(this->ModelError[i],0);
    }
  }
  if(this->CUT_EXP){ //Remove 0.05 rapidity from Exp
    Eigen::MatrixXd ER = Eigen::MatrixXd::Zero(this->ExperimentRapidity.size(),1);
    ER.col(0) = this->ExperimentRapidity;
    RemoveRow(ER,0);
    this->ExperimentRapidity = ER.col(0);
    for(int i=0;i<this->BFCount;i++){
      ER = this->Experiment[i];
      RemoveRow(ER,0);
      this->Experiment[i] = ER;
      ER = this->ExperimentalError[i];
      RemoveRow(ER,0);
      this->ExperimentalError[i] = ER;
    }
  }
  if(this->FIX_G0){ //Remove all g0 from Parameters
    int ParamPerPair = this->QParameters/this->QuarkPairs;
    for(int i=this->QuarkPairs-1;i>-1;i--){
      RemoveColumn(this->Parameters,i*ParamPerPair+1);
    }
    this->QParameters = this->Parameters.cols();
  }
  /*
  if(this->CUT_WIDTH){ //Remove all widths from Parameters
    this->WidthParameters = Eigen::MatrixXd::Zero(this->BFSamples,this->QuarkPairs);
    int ParamPerPair = this->QParameters/this->QuarkPairs;
    for(int i=this->QuarkPairs-1;i>-1;i--){
      this->WidthParameters.col(i) = this->Parameters.col(i*ParamPerPair);
      RemoveColumn(this->Parameters,i*ParamPerPair);
    }
    this->QParameters = this->Parameters.cols();
  }
  */
  if(this->COMMON_WIDTH){ //Remove all but one width from Parameters
    this->WidthParameters = Eigen::MatrixXd::Zero(this->BFSamples,this->QuarkPairs);
    int ParamPerPair = this->QParameters/this->QuarkPairs;
    for(int i=this->QuarkPairs-1;i>-1;i--){
      this->WidthParameters.col(i) = this->Parameters.col(i*ParamPerPair);
      if(i!=0) RemoveColumn(this->Parameters,i*ParamPerPair);
    }
    this->QParameters = this->Parameters.cols();
  }

  int model = ModelRapidity.size();
  int exp = ExperimentRapidity.size();
  if(exp > model){ //Cut Exp BF rapidity range down to match Model rapidity range
    Eigen::MatrixXd ER = Eigen::MatrixXd::Zero(model,1),
      temp =  Eigen::MatrixXd::Zero(exp,1);
    temp.col(0) = this->ExperimentRapidity;
    ER = temp.block(0,0,model,1);
    this->ExperimentRapidity = ER.col(0);
    for(int i=0;i<this->BFCount;i++){
      temp = this->Experiment[i];
      this->Experiment[i] = temp.block(0,0,model,1);
      temp = this->ExperimentalError[i];
      this->ExperimentalError[i] = temp.block(0,0,model,1);
    }
  }
}
void CBalanceModel::WriteBF(){
  printf("Writing BF files to %s/%s...\n",this->OUTPUT_FOLDER.c_str(),this->RUN_FOLDER.c_str());

  int etaCount = this->ModelRapidity.size();
  Eigen::MatrixXd print = Eigen::MatrixXd::Zero(etaCount,this->BFSamples+1);
  print.block(0,0,etaCount,1) = this->ModelRapidity;
  for(int i=0;i<this->BFCount;i++){
    print.block(0,1,etaCount,this->BFSamples) = this->Model[i];
    this->WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->MODEL_NAMES[i],print);
    print.block(0,1,etaCount,this->BFSamples) = this->ModelError[i];
    this->WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+"error_"+this->MODEL_NAMES[i],print);
  }
}
void CBalanceModel::ScaleParameters(){
  int rows = this->Parameters.rows(),
    cols = this->Parameters.cols();
  this->MinMax = Eigen::MatrixXd::Zero(cols,2);
  this->UnscaledParameters = this->Parameters;
  for(int col=0;col<cols;col++){
      this->MinMax(col,0) = this->UnscaledParameters.col(col).minCoeff();
      this->MinMax(col,1) = this->UnscaledParameters.col(col).maxCoeff();
      for(int row=0;row<rows;row++){
	  this->Parameters(row,col) = (this->UnscaledParameters(row,col) - this->MinMax(col,0))/(this->MinMax(col,1) - this->MinMax(col,0));
      }
  }

}
void CBalanceModel::ReduceDimensionsSeparately(){
  int index,size;
  CAnalysis Mean;
  this->PCA=std::vector<CAnalysis>(this->BFCount);
  this->McmcError = Eigen::MatrixXd::Zero(1,0);
  this->ModelZ = Eigen::MatrixXd::Zero(this->BFSamples,0);
  this->ExpZ = Eigen::MatrixXd::Zero(1,0);
  if(this->MEAN_ERROR){
      printf("Using Mean of Model Error\n");
  }
  for(int i=0;i<this->BFCount;i++){
    this->PCA[i].Data = this->Model[i];
    this->PCA[i].ComputeMean();
    this->PCA[i].Error = Eigen::MatrixXd::Zero(this->Model[i].rows(),this->Model[i].cols());
    if(this->MEAN_ERROR){
      Mean.Data = this->ModelError[i];
      Mean.ComputeMean();
      this->PCA[i].SumErrorInQuadrature(Mean.Mean);
    }else{
      this->PCA[i].SumErrorInQuadrature(this->ModelError[i]);
    }
    if(this->USE_MEAN_AS_ERROR){
      Eigen::VectorXd mean = this->PCA[i].Mean*this->MEAN_AS_ERROR;
      this->PCA[i].SumErrorInQuadrature(mean);
    }

    this->PCA[i].SumErrorInQuadrature(this->ExperimentalError[i]);
    this->PCA[i].ComputeTilde();
    this->PCA[i].ComputeCovariance();
    if(this->REDUCE_DIM){
      this->PCA[i].EigenSolve();
      this->PCA[i].EigenSort();
      this->PCA[i].ComputeZ();
      index=0;
      for(int j=0;j<this->PCA[i].EigenValues.size();j++){
	if(this->PCA[i].EigenValues(j)<this->PCA_CUTOFF) break;
	index++;
      }
      this->WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+"eigenvectors"+std::to_string(i)+"."+this->FILE_EXT,this->PCA[i].EigenVectors);
      this->WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+"eigenvalues"+std::to_string(i)+"."+this->FILE_EXT,this->PCA[i].EigenValues);
      std::cout << "Keeping " << index << " of " << this->PCA[i].EigenValues.size()<<  " PC: " <<  this->PCA[i].EigenValues.transpose().block(0,0,1,index) << std::endl;

      size=this->ModelZ.cols();
      this->ModelZ.conservativeResize(this->BFSamples,size+index);
      this->ModelZ.block(0,size,this->BFSamples,index) = this->PCA[i].Z.block(0,0,this->BFSamples,index);

      this->McmcError.conservativeResize(1,size+index);
      for(int j=0;j<index;j++){
	this->McmcError(0,size+j) = sqrt(1+this->MCMC_PERCENT_ERROR*this->MCMC_PERCENT_ERROR*this->PCA[i].EigenValues(j));
      }
    }else{
      size=this->PCA[i].Tilde.rows()*i;
      index=this->PCA[i].Tilde.rows();
      this->ModelZ.conservativeResize(this->BFSamples,index+size);
      this->ModelZ.block(0,size,this->BFSamples,index) = this->PCA[i].Tilde.transpose();
      this->McmcError.conservativeResize(1,size+index);
      for(int j=0;j<index;j++){
	this->McmcError(0,size+j) = 1;
      }
    }

    PCA[i].Data = this->Experiment[i];
    PCA[i].ComputeTilde();
    PCA[i].ComputeZ();
    if(this->REDUCE_DIM){
      this->ExpZ.conservativeResize(1,size+index);
      this->ExpZ.block(0,size,1,index) = this->PCA[i].Z.block(0,0,1,index);
    }else{
      this->ExpZ.conservativeResize(1,size+index);
      this->ExpZ.block(0,size,1,index) = this->PCA[i].Tilde.transpose();
    }
  }
  std::cout << "MCMC Error: " << this->McmcError << std::endl;
}
void CBalanceModel::ReduceDimensions(){
  int index,size;
  CAnalysis Mean;
  size = this->Model[0].rows();
  Eigen::VectorXd ExpError = Eigen::VectorXd::Zero(this->BFCount*size);
  Eigen::MatrixXd ModError = Eigen::MatrixXd::Zero(this->BFCount*size,this->BFSamples);  
  this->PCA_RD.Data = Eigen::MatrixXd::Zero(this->BFCount*size,this->BFSamples);
  for(int i=0;i<this->BFCount;i++){
    this->PCA_RD.Data.block(size*i,0,size,this->BFSamples) = this->Model[i];
    ModError.block(size*i,0,size,this->BFSamples) = this->ModelError[i];
    ExpError.segment(size*i,size) = this->ExperimentalError[i];
  }

  this->PCA_RD.ComputeMean();
  this->PCA_RD.Error = Eigen::MatrixXd::Zero(this->BFCount*size,this->BFSamples);
  this->PCA_RD.SumErrorInQuadrature(ExpError);
  if(this->MEAN_ERROR){
    Mean.Data = ModError;
    Mean.ComputeMean();
    this->PCA_RD.SumErrorInQuadrature(Mean.Mean);
    printf("Using Mean of Model Error\n");
  }else{
    this->PCA_RD.SumErrorInQuadrature(ModError);
  }
  if(this->USE_MEAN_AS_ERROR){
    Eigen::VectorXd mean = this->PCA_RD.Mean*this->MEAN_AS_ERROR;
    printf("Using Mean of Model As Error, fraction: %f\n",this->MEAN_AS_ERROR);
    this->PCA_RD.SumErrorInQuadrature(mean);
  }

  this->PCA_RD.ComputeTilde();
  if(this->REDUCE_DIM){
    this->PCA_RD.ComputeCovariance();
    this->PCA_RD.EigenSolve();
    this->PCA_RD.EigenSort();
    this->PCA_RD.ComputeZ();
    index=0;
    for(int j=0;j<this->PCA_RD.EigenValues.size();j++){
      if(this->PCA_RD.EigenValues(j)<this->PCA_CUTOFF) break;
      index++;
    }
    this->WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+"eigenvectors."+this->FILE_EXT,this->PCA_RD.EigenVectors);
    this->WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+"eigenvalues."+this->FILE_EXT,this->PCA_RD.EigenValues);
    std::cout << "Keeping " << index << " of " << this->PCA_RD.EigenValues.size() << " PC: " << this->PCA_RD.EigenValues.transpose().block(0,0,1,index) << std::endl;

    this->ModelZ = this->PCA_RD.Z.block(0,0,this->BFSamples,index);
    this->McmcError = Eigen::MatrixXd::Zero(1,index);
    for(int j=0;j<index;j++){
      this->McmcError(0,j) = sqrt(1+this->MCMC_PERCENT_ERROR*this->MCMC_PERCENT_ERROR*this->PCA_RD.EigenValues(j));
    }
    for(int i=0;i<this->BFCount;i++){
      this->PCA_RD.Data.block(size*i,0,size,1) = this->Experiment[i];
    }
    PCA_RD.ComputeTilde();
    PCA_RD.ComputeZ();
    this->ExpZ = this->PCA_RD.Z.block(0,0,1,index);
  }else{
    index=this->PCA_RD.Tilde.transpose().cols();
    this->ModelZ = this->PCA_RD.Tilde.transpose();
    this->McmcError = Eigen::MatrixXd::Zero(1,index);
    for(int j=0;j<index;j++){
      this->McmcError(0,j) = 1;
    }
    for(int i=0;i<this->BFCount;i++){
      this->PCA_RD.Data.block(size*i,0,size,1) = this->Experiment[i];
    }
    PCA_RD.ComputeTilde();
    this->ExpZ = this->PCA_RD.Tilde.transpose();
  }
  
  std::cout << "MCMC Error: " << this->McmcError << std::endl;
}

void CBalanceModel::WriteZ(){
  CSystem SYSTEM;
  int Zs = this->ModelZ.cols();
  
  Eigen::MatrixXd plot = Eigen::MatrixXd::Zero(this->BFSamples,this->QParameters+Zs);
  plot.block(0,0,this->BFSamples,this->QParameters) = this->Parameters;
  plot.block(0,this->QParameters,this->BFSamples,Zs) = this->ModelZ;
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->TRAINPOINTS_FILE+"."+this->FILE_EXT,plot);

  plot = ExpZ;
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->EXPZ_FILE+"."+this->FILE_EXT,plot);
}
void CBalanceModel::WriteReconBF(){
  CSystem SYSTEM;
  int Zs = this->ModelZ.cols();
  
  Eigen::MatrixXd plot = Eigen::MatrixXd::Zero(this->BFSamples,this->QParameters+Zs);
  plot.block(0,0,this->BFSamples,this->QParameters) = this->Parameters;
  plot.block(0,this->QParameters,this->BFSamples,Zs) = this->ModelZ;
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->TRAINPOINTS_FILE+"."+this->FILE_EXT,plot);

  plot = ExpZ;
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->EXPZ_FILE+"."+this->FILE_EXT,plot);
}
void CBalanceModel::CreateEmulator(){
  printf("Creating %s emulator...\n",this->EMULATOR_CHOICE.c_str());
  if(this->EMULATOR_CHOICE=="GAUSSIAN_PROCESS"){
    this->Emulator = new CGaussianProcess;
    this->Emulator->Construct(this->Parameters,this->ModelZ,this->MAP);
  }else if(this->EMULATOR_CHOICE=="NEURAL_NET"){
    this->Emulator = new CNeuralNet;
    this->Emulator->Construct(this->Parameters,this->ModelZ,this->MAP);
  }else{
    printf("Error. Emulator choice %s not recognized. Terminating program.\n",this->EMULATOR_CHOICE.c_str());
    exit(1);
  }
}
void CBalanceModel::DeleteEmulator(){
  delete this->Emulator;
}
void CBalanceModel::CreateMCMC(){
  printf("Creating %s MCMC Likelihood Function...\n",this->MCMC_CHOICE.c_str());
  Eigen::MatrixXd 
    newRange = Eigen::MatrixXd::Zero(this->QParameters,2);
  Eigen::VectorXd
    newWidth = Eigen::VectorXd::Zero(this->QParameters);
  for(int i=0;i<this->QParameters;i++){
    newWidth(i) = this->MCMC_WIDTH;
    newRange(i,0) = this->MCMC_MIN;
    newRange(i,1) = this->MCMC_MAX;
  }
  this->MCMC.Construct(newRange,newWidth, this->MCMC_CHOICE);
}
void CBalanceModel::RunMCMC(){
  printf("Running MCMC with STEP %f, RUNS %d, recording history after %d\n",this->MCMC_WIDTH,this->NSAMPLES,this->NTRACE);
  Eigen::MatrixXd TargetZ;
  Eigen::MatrixXd print;

  CSystem SYSTEM;
  if(this->testZ == -1){
    TargetZ = this->ExpZ;
    std::cout << "ExpZ: " << TargetZ << std::endl;
    print = Eigen::MatrixXd::Zero(1,TargetZ.cols());
    print.block(0,0,1,TargetZ.cols()) = TargetZ;
  }else{
    TargetZ = this->ModelZ.row(this->testZ);
    print = Eigen::MatrixXd::Zero(1,this->QParameters + TargetZ.cols());
    print.block(0,0,1,this->QParameters) = this->Parameters.row(this->testZ);
    print.block(0,this->QParameters,1,TargetZ.cols()) = TargetZ;

    std::cout << "Param" << this->testZ << ": " << this->Parameters.row(this->testZ) << std::endl;
    std::cout << "ModZ" << this->testZ << ": " << TargetZ << std::endl;
  }
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/targetz."+this->FILE_EXT,print);

  this->MCMCHistory = this->MCMC.Run(this->Emulator,TargetZ,this->McmcError,this->NSAMPLES);
  
  //Eigen::MatrixXd toKeep = this->MCMCHistory.block(this->NTRACE/5,0,this->MCMCHistory.rows()-this->NTRACE/5,this->MCMCHistory.cols());
  //this->MCMCHistory = toKeep;
}
void CBalanceModel::WriteMCMC(){
  printf("Writing MCMC files to %s/%s...\n",this->OUTPUT_FOLDER.c_str(),this->RUN_FOLDER.c_str());

  CSystem SYSTEM;
  int printInt = this->MCMCHistory.rows();
  //int testZ = this->MAP.getI("MCMC_TEST_Z",-1);
  Eigen::MatrixXd print = this->MCMCHistory.block(0,1,printInt,this->QParameters);
  std::vector<std::string> header(this->QParameters);
  Eigen::VectorXd MaxLLH(1);
  MaxLLH(0) = this->MCMC.maxLogLikelihood;
  for(int i=0;i<this->QParameters;i++){
    header[i] = "QParameter_"+std::to_string(i);
  }

  SYSTEM.WriteCSVFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+MCMCTRACE_FILE+"."+this->CSV_EXT,header,print);
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+MCMCHISTORY_FILE+"."+this->FILE_EXT,this->MCMCHistory);
  SYSTEM.WriteFile(this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+"maxloglikelihood"+"."+this->FILE_EXT,MaxLLH);
}
void CBalanceModel::WriteCoshFunctions(){
  printf("Writing GAB COSH files to %s/%s...\n",this->OUTPUT_FOLDER.c_str(),this->RUN_FOLDER.c_str());
  CSystem SYSTEM;
  std::string filename;
  int ab = this->QuarkPairs;
  int nmax = this->QParameters/ab - 1;
  CDistCosh dist;

  Eigen::MatrixXd 
    FullG,
    temp;
  if(this->SCALE_X){
    temp = this->UnscaledParameters;
  }else{
    temp = this->Parameters;
  }
  /*
  if(this->CUT_WIDTH){
    int ParamPerPair = this->QParameters/this->QuarkPairs + 1;
    for(int i=0;i<this->QuarkPairs;i++){
      AddOnesColumn(temp, FullG, i*ParamPerPair);
      //FullG.col(i*ParamPerPair) = this->WidthParameters.col(i);
      FullG.col(i*ParamPerPair) *= this->WidthParameters(0,i);
      temp = FullG;
    }
  }else{
    nmax--;
    FullG = temp;
  }
  */

  if(this->COMMON_WIDTH){
    int ParamPerPair = (this->QParameters-1)/this->QuarkPairs + 1;
    for(int i=1;i<this->QuarkPairs;i++){
      AddOnesColumn(temp, FullG, i*ParamPerPair);
      //FullG.col(i*ParamPerPair) = this->WidthParameters.col(0);
      FullG.col(i*ParamPerPair) = FullG.col(0);
      temp = FullG;
    }
  }else{
    nmax--;
    FullG = temp;
  }

  if(this->FIX_G0){ //Add G0 from unscaled parameters
    FullG = dist.GenG0(this->QuarkPairs,temp);
    nmax++;
  }else{
    FullG = temp;
  }

  /*
  if(this->CUT_G0 && this->SCALE_X){ //Add G0 from unscaled parameters
    FullG = dist.GenG0(this->QuarkPairs,this->UnscaledParameters);
  }else if(this->CUT_G0){ //Add G0 back in to create gabfunctions 
    FullG = dist.GenG0(this->QuarkPairs,this->Parameters);
  }else if(this->SCALE_X){ //Use unscaled X
    FullG = this->UnscaledParameters;
  }else{ // No alterations, simply use parameters
    FullG = this->Parameters;
  }
  */
  filename = this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->FULLG_FILE+"."+this->FILE_EXT;
  SYSTEM.WriteFile(filename,FullG); //Write prior FullG
  Eigen::MatrixXd //Create prior gabfunctions
    Functions = dist.FunctionSet(this->GRID,this->ETA_MAX,
				 this->BFSamples,
				 this->QuarkPairs,
				 nmax,
				 FullG);
  this->ScaleByChi(Functions); //Change unit area functions to Chi area

  filename = this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->GAB_FILE+"."+this->FILE_EXT;
  SYSTEM.WriteFile(filename,Functions); //Write prior gabfunctions


  Eigen::MatrixXd //extract mcmctrace from mcmchistory
    mcmc = this->MCMCHistory.block(0,1,this->MCMCHistory.rows(),this->QParameters),
    posterior = Eigen::MatrixXd::Zero(MCMC_POST,this->QParameters);
  posterior = ExtractPosterior(mcmc); //Extract posterior samples from mcmctrace
  if(this->SCALE_X){ //Check if Parameters were scaled for mcmc
    filename = this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->MINMAX_FILE+"."+this->FILE_EXT;
    SYSTEM.WriteFile(filename,this->MinMax); //Write MinMax from scaling Parameters
    posterior = this->UnscaleParameters(posterior); //Unscale posterior trace
  }
  /*
  if(this->CUT_WIDTH){
    int ParamPerPair = this->QParameters/this->QuarkPairs + 1;
    for(int i=0;i<this->QuarkPairs;i++){
      AddOnesColumn(posterior, FullG, i*ParamPerPair);
      FullG.col(i*ParamPerPair) *= this->WidthParameters(0,i);
      posterior = FullG;
    }
  }else{
    FullG = posterior;    
  }
  */

  if(this->COMMON_WIDTH){
    int ParamPerPair = (this->QParameters-1)/this->QuarkPairs + 1;
    for(int i=1;i<this->QuarkPairs;i++){
      AddOnesColumn(posterior, FullG, i*ParamPerPair);
      FullG.col(i*ParamPerPair) = FullG.col(0);
      posterior = FullG;
    }
  }else{
    FullG = posterior;    
  }

  if(this->FIX_G0){ //Add G0 back in to create gabfunctions
    FullG = dist.GenG0(this->QuarkPairs,posterior);
    posterior = FullG;
  }else{
    FullG = posterior;
  }

  filename = this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->FULLG_FILE+this->POST_EXT+"."+this->FILE_EXT;
  SYSTEM.WriteFile(filename,FullG); //Write posterior FullG

  Functions = dist.FunctionSet(this->GRID,this->ETA_MAX, //Create posterior gabfunctions
			       this->MCMC_POST,
			       this->QuarkPairs,
			       nmax,
			       FullG);
  this->ScaleByChi(Functions); //Change unit area functions to Chi area
  

  filename = this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->MCMCTRACE_FILE+this->POST_EXT+"."+this->FILE_EXT;
  SYSTEM.WriteFile(filename,FullG); //Write posterior trace

  filename =this->OUTPUT_FOLDER+"/"+this->RUN_FOLDER+"/"+this->GAB_FILE+this->POST_EXT+"."+this->FILE_EXT;
  SYSTEM.WriteFile(filename,Functions); //Write posterior gabfunctions
}
Eigen::MatrixXd CBalanceModel::UnscaleParameters(Eigen::MatrixXd ScaledParameters){
  int rows = ScaledParameters.rows(),
    cols = ScaledParameters.cols();
  Eigen::MatrixXd UnscaledParameters = Eigen::MatrixXd::Zero(rows,cols);
  for(int col=0;col<cols;col++){
      for(int row=0;row<rows;row++){
	UnscaledParameters(row,col) = ScaledParameters(row,col)*(this->MinMax(col,1) - this->MinMax(col,0)) + this->MinMax(col,0);
      }
  }
  return UnscaledParameters;
}
void CBalanceModel::ScaleByChi(Eigen::MatrixXd &Functions){
  int cols=(Functions.cols()-1)/this->QuarkPairs;
  for(int i=0;i<this->QuarkPairs;i++){
    for(int j=0;j<cols;j++){
      Functions.col(1+j*this->QuarkPairs+i) = this->CHI(i)*Functions.col(1+j*this->QuarkPairs+i);
    }
  }
}
Eigen::MatrixXd CBalanceModel::ExtractPosterior(Eigen::MatrixXd mcmc){
  Eigen::MatrixXd 
    posterior = Eigen::MatrixXd::Zero(this->MCMC_POST,this->QParameters);
  int incr = mcmc.rows()/MCMC_POST;
  for(int i=0;i<this->MCMC_POST;i++){
    posterior.row(i) = mcmc.row(incr*i);
  }
  return posterior;
}
