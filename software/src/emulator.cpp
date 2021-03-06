#include "emulator.h"

void AddOnesColumn(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix, int column){
  int rows = matrix.rows(),
    cols = matrix.cols();
  outMatrix = Eigen::MatrixXd::Zero(rows,cols+1);
  for(int row=0;row<rows;row++){
    outMatrix(row,column) = 1.0;
    for(int col=0;col<cols;col++){
      if(col < column){
	outMatrix(row,col) = matrix(row,col);
      }else{
	outMatrix(row,col+1) = matrix(row,col);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////

/**GAUSSIAN PROCESS--
******************/

///////////////////////////////////////////////////////////////////////////////////////


CGaussianProcess::CGaussianProcess(){
}
CGaussianProcess::CGaussianProcess(Eigen::MatrixXd newX, Eigen::MatrixXd newY, CParameterMap MAP){
  this->Epsilon = MAP.getD("EPSILON",1e-8);
  this->SigmaF = MAP.getD("SIGMA_F",0.5);
  this->CharacLength = MAP.getD("CHARAC_LENGTH",0.45);
  this->SigmaNoise = MAP.getD("SIGMA_NOISE",0.05);

  this->X = newX;
  this->Y = newY;

  this->trainPoints = this->X.rows();
  this->paramCount = this->X.cols();
  this->obsCount = this->Y.cols();

  this->ConstructHyperparameters();
  this->ConstructBeta();

  this->Noise = this->Hyperparameters.col(this->hyperparamCount - 1);

  Eigen::MatrixXd trainI = Eigen::MatrixXd::Identity(trainPoints, trainPoints);

  this->Kernel = std::vector<Eigen::MatrixXd>(obsCount);
  this->KernelInv = std::vector<Eigen::MatrixXd>(obsCount);
  this->HMatrix = std::vector<Eigen::MatrixXd>(obsCount);
  for(int ob=0;ob<obsCount;ob++){
    this->Kernel[ob] = this->KernelFunction(this->X,this->X,ob);
    this->Kernel[ob] += trainI*(this->Noise(ob)*this->Noise(ob) + this->Epsilon);
    this->KernelInv[ob] = this->Kernel[ob].inverse();
    this->HMatrix[ob] = RegressionLinearFunction(this->X, ob);
  }
}
CGaussianProcess::CGaussianProcess(Eigen::MatrixXd newX, Eigen::MatrixXd newY, std::string filename){
  CParameterMap MAP;
  MAP.ReadParsFromFile(filename);

  this->Epsilon = MAP.getD("EPSILON",1e-8);
  this->SigmaF = MAP.getD("SIGMA_F",0.5);
  this->CharacLength = MAP.getD("CHARAC_LENGTH",0.45);
  this->SigmaNoise = MAP.getD("SIGMA_NOISE",0.05);

  this->X = newX;
  this->Y = newY;

  this->trainPoints = this->X.rows();
  this->paramCount = this->X.cols();
  this->obsCount = this->Y.cols();

  this->ConstructHyperparameters();
  this->ConstructBeta();
  this->Noise = this->Hyperparameters.col(this->hyperparamCount - 1);

  Eigen::MatrixXd trainI = Eigen::MatrixXd::Identity(trainPoints, trainPoints);

  this->Kernel = std::vector<Eigen::MatrixXd>(obsCount);
  this->KernelInv = std::vector<Eigen::MatrixXd>(obsCount);
  this->HMatrix = std::vector<Eigen::MatrixXd>(obsCount);
  for(int ob=0;ob<obsCount;ob++){
    this->Kernel[ob] = this->KernelFunction(this->X,this->X,ob);
    this->Kernel[ob] += trainI*(this->Noise(ob)*this->Noise(ob) + this->Epsilon);
    this->KernelInv[ob] = this->Kernel[ob].inverse();
    this->HMatrix[ob] = RegressionLinearFunction(this->X, ob);
  }
}
///////////////////////////////////////////////////////////////////////////////////////

void CGaussianProcess::Construct(Eigen::MatrixXd newX, Eigen::MatrixXd newY, CParameterMap MAP){
  this->Epsilon = MAP.getD("EPSILON",1e-8);
  this->SigmaF = MAP.getD("SIGMA_F",0.5);
  this->CharacLength = MAP.getD("CHARAC_LENGTH",0.45);
  this->SigmaNoise = MAP.getD("SIGMA_NOISE",0.05);

  this->X = newX;
  this->Y = newY;

  this->trainPoints = this->X.rows();
  this->paramCount = this->X.cols();
  this->obsCount = this->Y.cols();

  this->ConstructHyperparameters();
  this->ConstructBeta();
  this->Noise = this->Hyperparameters.col(this->hyperparamCount - 1);

  Eigen::MatrixXd trainI = Eigen::MatrixXd::Identity(trainPoints, trainPoints);
  Eigen::MatrixXd tempPP;

  this->Kernel = std::vector<Eigen::MatrixXd>(obsCount);
  this->KernelInv = std::vector<Eigen::MatrixXd>(obsCount);
  this->HMatrix = std::vector<Eigen::MatrixXd>(obsCount);
  this->HMatrix_KernelInv = std::vector<Eigen::MatrixXd>(obsCount);
  this->KernelInv_Y = std::vector<Eigen::MatrixXd>(obsCount);
  this->betaMatrix = std::vector<Eigen::MatrixXd>(obsCount);

  for(int ob=0;ob<obsCount;ob++){
    this->Kernel[ob] = this->KernelFunction(this->X,this->X,ob);
    this->Kernel[ob] += trainI*(this->Noise(ob)*this->Noise(ob) + this->Epsilon);
    this->KernelInv[ob] = this->Kernel[ob].inverse();
    this->HMatrix[ob] = RegressionLinearFunction(this->X, ob);
    this->HMatrix_KernelInv[ob] = this->HMatrix[ob]*this->KernelInv[ob];
    this->KernelInv_Y[ob] = this->KernelInv[ob]*this->Y.col(ob);
    tempPP = this->HMatrix_KernelInv[ob]*this->HMatrix[ob].transpose();
    this->betaMatrix[ob] = tempPP.inverse()*this->HMatrix_KernelInv[ob]*this->Y.col(ob);
  }
}
void CGaussianProcess::Construct(Eigen::MatrixXd newX, Eigen::MatrixXd newY, std::string filename){
  CParameterMap MAP;
  MAP.ReadParsFromFile(filename);

  this->Epsilon = MAP.getD("EPSILON",1e-8);
  this->SigmaF = MAP.getD("SIGMA_F",0.5);
  this->CharacLength = MAP.getD("CHARAC_LENGTH",0.45);
  this->SigmaNoise = MAP.getD("SIGMA_NOISE",0.05);

  this->X = newX;
  this->Y = newY;

  this->trainPoints = this->X.rows();
  this->paramCount = this->X.cols();
  this->obsCount = this->Y.cols();

  this->ConstructHyperparameters();
  this->ConstructBeta();
  this->Noise = this->Hyperparameters.col(this->hyperparamCount - 1);

  Eigen::MatrixXd trainI = Eigen::MatrixXd::Identity(trainPoints, trainPoints);

  this->Kernel = std::vector<Eigen::MatrixXd>(obsCount);
  this->KernelInv = std::vector<Eigen::MatrixXd>(obsCount);
  this->HMatrix = std::vector<Eigen::MatrixXd>(obsCount);
  for(int ob=0;ob<obsCount;ob++){
    this->Kernel[ob] = this->KernelFunction(this->X,this->X,ob);
    this->Kernel[ob] += trainI*(this->Noise(ob)*this->Noise(ob) + this->Epsilon);
    this->KernelInv[ob] = this->Kernel[ob].inverse();
    this->HMatrix[ob] = RegressionLinearFunction(this->X, ob);
  }
}
void CGaussianProcess::ConstructHyperparameters(){
  int cols = this->Y.cols();
  this->Hyperparameters = Eigen::MatrixXd::Zero(cols,3);
  for(int col=0;col<cols;col++){
      Hyperparameters(col,0) = this->SigmaF*(this->Y.col(col).maxCoeff() - this->Y.col(col).minCoeff());
      Hyperparameters(col,1) = this->CharacLength;
      Hyperparameters(col,2) = this->SigmaNoise*Hyperparameters(col,0);
  }
  this->hyperparamCount = this->Hyperparameters.cols();
}
void CGaussianProcess::ConstructBeta(){
  this->Beta = Eigen::MatrixXd::Zero(this->obsCount,this->paramCount+1); 

  Eigen::MatrixXd
    temp = Eigen::MatrixXd::Zero(this->paramCount+1,this->paramCount+1),
    X_;
  AddOnesColumn(this->X,X_,0);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(this->paramCount);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(this->trainPoints);

  temp = X_.transpose()*X_;
  temp = temp.inverse();
  for(int j=0;j<this->obsCount;j++){
    y = Y.col(j);
    beta = temp*X_.transpose()*y;
    this->Beta.row(j) = beta;
  }
}
Eigen::MatrixXd CGaussianProcess::RegressionLinearFunction(Eigen::MatrixXd testX, int obsIndex){
  int points = testX.rows();
  Eigen::MatrixXd HMat = Eigen::MatrixXd::Zero(1,points);
  for(int i=0;i<points;i++){
    HMat(0,i) = this->Beta(obsIndex,0);
    for(int j=0;j<this->paramCount;j++){
      HMat(0,i) += this->Beta(obsIndex,j+1)*testX(i,j);
    }
  }
  return HMat;
}
Eigen::MatrixXd CGaussianProcess::KernelFunction(Eigen::MatrixXd A, Eigen::MatrixXd B, int obsIndex){
  double xi, xj, sum, diff,
    theta1 = this->Hyperparameters(obsIndex,0)*this->Hyperparameters(obsIndex,0),
    theta2 = 0.5/(this->Hyperparameters(obsIndex,1)*this->Hyperparameters(obsIndex,1));
  int ai = A.rows(),
    bi = B.rows();
  Eigen::MatrixXd KernelMatrix = Eigen::MatrixXd::Zero(ai,bi);

  for(int a=0;a<ai;a++){
    for(int b=0;b<bi;b++){
      sum = 0.0;
      for(int pc=0;pc<this->paramCount;pc++){
	xi = A(a,pc);
	xj = B(b,pc);
	diff = xi - xj;
	sum += diff*diff;
      }
      KernelMatrix(a,b) = theta1*exp(-sum*theta2);
    }
  }
  return KernelMatrix;
}
Eigen::MatrixXd CGaussianProcess::Emulate(Eigen::MatrixXd testX){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0.0,1.0);

  int testPoints = testX.rows();
  //matrices and vectors for computational use
  Eigen::MatrixXd KernelS = Eigen::MatrixXd::Zero(this->trainPoints,testPoints),
    KernelSS = Eigen::MatrixXd::Zero(testPoints,testPoints),
    
    meanMatrix = Eigen::MatrixXd::Zero(testPoints,1),
    //varMatrix = Eigen::MatrixXd::Zero(testPoints, testPoints),
    //L = Eigen::MatrixXd::Zero(testPoints, testPoints),

    HS = Eigen::MatrixXd::Zero(1, testPoints),
    R = Eigen::MatrixXd::Zero(this->trainPoints, testPoints),
    //betaMatrix = Eigen::MatrixXd::Zero(this->trainPoints, 1),

    //tempPP = Eigen::MatrixXd::Zero(1,1),
    //tempPPInv = Eigen::MatrixXd::Zero(1,1),

    postFunc = Eigen::MatrixXd::Zero(testPoints, 1),
    randomSample = Eigen::MatrixXd::Zero(testPoints, 1),

    testI = Eigen::MatrixXd::Identity(testPoints,testPoints);
  
  //outMatrix = Eigen::MatrixXd::Zero(testPoints,this->paramCount+this->obsCount);
  Eigen::MatrixXd outMatrix;
    outMatrix = Eigen::MatrixXd::Zero(testPoints,this->obsCount);

  for(int ob=0;ob<this->obsCount;ob++){
    KernelS = this->KernelFunction(this->X, testX, ob);
    KernelSS = this->KernelFunction(testX, testX, ob);
    KernelSS += testI*this->Epsilon;

    HS = this->RegressionLinearFunction(testX, ob);
    R = HS - this->HMatrix_KernelInv[ob]*KernelS;

    meanMatrix = KernelS.transpose()*this->KernelInv_Y[ob] + R.transpose()*this->betaMatrix[ob];

    /*
    HS = this->RegressionLinearFunction(testX, ob);
    R = HS - this->HMatrix[ob]*this->KernelInv[ob]*KernelS;
    tempPP = this->HMatrix[ob]*KernelInv[ob]*this->HMatrix[ob].transpose();
    tempPPInv = tempPP.inverse();
    betaMatrix = tempPPInv*this->HMatrix[ob]*this->KernelInv[ob]*this->Y.col(ob);

    meanMatrix = KernelS.transpose()*this->KernelInv[ob]*Y.col(ob) + R.transpose()*betaMatrix;
    */
    //varMatrix = KernelSS - KernelS.transpose()*this->KernelInv[ob]*KernelS + R.transpose()*tempPPInv*R;
    //L = varMatrix.llt().matrixL();
    
    outMatrix.col(ob) = meanMatrix.col(0);
  }
  return outMatrix;
}
Eigen::MatrixXd CGaussianProcess::Emulate_NR(Eigen::MatrixXd testX){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0.0,1.0);

  int testPoints = testX.rows();
  
  //matrices and vectors for computational use
  Eigen::MatrixXd KernelS = Eigen::MatrixXd::Zero(this->trainPoints,testPoints),
    KernelSS = Eigen::MatrixXd::Zero(testPoints,testPoints),
    
    meanMatrix = Eigen::MatrixXd::Zero(testPoints,1),
    varMatrix = Eigen::MatrixXd::Zero(testPoints, testPoints),
    L = Eigen::MatrixXd::Zero(testPoints, testPoints),

    //postFunc = Eigen::MatrixXd::Zero(testPoints, 1),
    //randomSample = Eigen::MatrixXd::Zero(testPoints, 1),

    testI = Eigen::MatrixXd::Identity(testPoints,testPoints);

  Eigen::MatrixXd outMatrix;
    outMatrix = Eigen::MatrixXd::Zero(testPoints,this->obsCount);

  for(int ob=0;ob<obsCount;ob++){
    KernelS = this->KernelFunction(this->X, testX, ob);
    KernelSS = this->KernelFunction(testX, testX, ob);
    KernelSS += testI*this->Epsilon;

    meanMatrix = KernelS.transpose()*this->KernelInv[ob]*Y.col(ob);

    //varMatrix = KernelSS - KernelS.transpose()*this->KernelInv[ob]*KernelS;
    //L = varMatrix.llt().matrixL();
    outMatrix.col(ob) = meanMatrix.col(0);
  }
  return outMatrix;
}

///////////////////////////////////////////////////////////////////////////////////////

/**NEURAL NETWORK--
******************/

///////////////////////////////////////////////////////////////////////////////////////


void CNeuralNet::Construct(Eigen::MatrixXd X_train, Eigen::MatrixXd Y_train, CParameterMap Map){
  this->X = X_train;
  this->Y = Y_train;
  this->MAP = Map;

  int Params = this->X.cols();
  int Obs = this->Y.cols();

  this->epochs = this->MAP.getI("EPOCHS",30);
  this->learning_rate = this->MAP.getD("LEARNING_RATE",0.001);
  this->regular_param = this->MAP.getD("REGULAR_PARAM",0.0);
  this->mini_batch = this->MAP.getD("MINI_BATCH_SIZE",10);
  this->momentum = this->MAP.getD("MOMENTUM",0.9);
  this->beta1 = this->MAP.getD("BETA1",0.9);
  this->beta2 = this->MAP.getD("BETA2",0.999);
  
  this->layers = this->MAP.getD("LAYERS",1);
  this->Layers.push_back(this->MAP.getI("LAYER0_SIZE",100));
  for(int i=1;i<layers;i++){
    std::string layer_name = "LAYER"+std::to_string(i)+"_SIZE";
    this->Layers.push_back(this->MAP.getI(layer_name,0));
  }
  this->Layers.insert(this->Layers.begin(),Params);
  this->Layers.push_back(Obs);
  this->layers = this->Layers.size();

  //Set up Activation Functions
  this->SActivation = this->MAP.getS("ACTIVATION","RELU");
  if(this->Activation != nullptr) delete this->Activation;
  if(this->SActivation=="SOFTMAX"){
    printf("Creating %s for Activation Function\n",this->SActivation.c_str());
    this->Activation = new CSoftMax;
  }else if(this->SActivation=="SIGMOID"){
    printf("Creating %s for Activation Function\n",this->SActivation.c_str());
    this->Activation = new CSigmoid;
  }else if(this->SActivation=="TANH"){
    printf("Creating %s for Activation Function\n",this->SActivation.c_str());
    this->Activation = new CTanh;
  }else if(this->SActivation=="IDENTITY"){
    printf("Creating %s for Activation Function\n",this->SActivation.c_str());
    this->Activation = new CIdentity;
  }else if(this->SActivation=="RELU"){
    printf("Creating %s for Activation Function\n",this->SActivation.c_str());
    this->Activation = new CRelu;
  }else{
    this->Activation = nullptr;
    printf("%s is not a valid option for CActivation in CLoss.\n",SActivation.c_str());
    exit(1);
  }

  //Set up Loss & Final Act Functions
  this->SLoss = this->MAP.getS("LOSS","ENTROPY");
  this->SFinalActivation = this->MAP.getS("FINAL_ACTIVATION","SOFTMAX");
  if(this->Loss != nullptr) delete this->Loss;
  if(this->SLoss=="ENTROPY"){
    printf("Creating %s for Loss Function\n",this->SLoss.c_str());
    this->Loss = new CEntropyLoss;
    printf("Creating %s for Final Activation Function\n","SOFTMAX");
    this->Loss->Construct("SOFTMAX",this->regular_param);
  }else if(this->SLoss=="L2"){
    printf("Creating %s for Loss Function\n",this->SLoss.c_str());
    this->Loss = new CL2Loss;
    printf("Creating %s for Final Activation Function\n",this->SFinalActivation.c_str());
    this->Loss->Construct(this->SFinalActivation,this->regular_param);
  }else{
    this->Loss = nullptr;
    printf("%s is not a valid option for CLoss.\n",this->SLoss.c_str());
    exit(1);
  }

  //Set up Solver
  this->SSolver = this->MAP.getS("SOLVER","SGD");
  if(this->Solver != nullptr) delete this->Solver;
  if(this->SSolver=="SGD"){
    printf("Creating %s for Solver\n",this->SSolver.c_str());
    this->Solver = new CSGD;
    this->Solver->Construct(this->learning_rate, this->regular_param);
    this->Solver->ConstructSGD(this->momentum, this->Layers);
  }else if(this->SSolver=="ADAM"){
    printf("Creating %s for Solver\n",this->SSolver.c_str());
    this->Solver = new CAdam;
    this->Solver->Construct(this->learning_rate, this->regular_param);
    this->Solver->ConstructAdam(this->beta1,this->beta2,this->Layers);
  }else{
    this->Loss = nullptr;
    printf("%s is not a valid option for CSolver.\n",this->SSolver.c_str());
    exit(1);
  }

  Eigen::MatrixXd weight;
  Eigen::MatrixXd bias;
  Eigen::MatrixXd del;
  int layer, next_layer;

  this->activations.push_back(Eigen::MatrixXd::Zero(mini_batch,Params));
  this->zs.push_back(Eigen::MatrixXd::Zero(mini_batch,Params));

  for(int i=0;i<this->layers-1;i++){
    layer=this->Layers[i];
    next_layer=this->Layers[i+1];

    weight = Eigen::MatrixXd::Zero(layer,next_layer);
    this->Weight.push_back(weight);
    this->delta_w.push_back(weight);
    this->Weight[i].setRandom();
    if(this->regular_param != 0.0){
      //this->Weight[i] /= sqrt((double) this->Weight[i].size());
    }

    bias = Eigen::MatrixXd::Zero(1,next_layer);
    this->Bias.push_back(bias);
    this->delta_b.push_back(bias);
    this->Bias[i].setRandom();
    if(this->regular_param != 0.0){
      //this->Bias[i] /= sqrt((double) this->Bias[i].size());
    }

    del = Eigen::MatrixXd::Zero(mini_batch,next_layer);
    this->delta.push_back(del);

    this->activations.push_back(del);
    this->zs.push_back(del);
  }

  this->Train(this->epochs);
}
void CNeuralNet::FeedForward(Eigen::MatrixXd x){
  int samples=x.rows();
  this->zs[0] = x;
  this->activations[0] = x;
  for(int ilayer=0;ilayer<this->layers-1;ilayer++){
    this->zs[ilayer+1] = this->activations[ilayer]*this->Weight[ilayer];
    for(int i=0;i<samples;i++){
      this->zs[ilayer+1].row(i) += this->Bias[ilayer];
    }
    this->activations[ilayer+1] = this->Activation->Function(this->zs[ilayer+1]);
  }
}
void CNeuralNet::BackPropagation(Eigen::MatrixXd x, Eigen::MatrixXd y){
  int layer, n_layer,final_layer;
  int samples = x.rows();
  final_layer = this->layers-1;
  layer = this->Layers[final_layer-1];
  n_layer = this->Layers[final_layer];

  Eigen::MatrixXd yrun = this->Emulate(x);
  Eigen::MatrixXd cost_der = this->Loss->Derivative(this->zs[final_layer],y);
  Eigen::MatrixXd act_der;// = this->ActivationDerivative(zs[final_layer]);

  this->delta[final_layer-1] = Eigen::MatrixXd::Zero(samples,n_layer);

  for(int i=0;i<n_layer;i++){
    this->delta_b[final_layer-1](0,i) = 0.0; // rezero d_b
    for(int samp=0;samp<samples;samp++){
      this->delta[final_layer-1](samp,i) = cost_der(samp,i);//*act_der(i,samp); //set delta
      this->delta_b[final_layer-1](0,i) += this->delta[final_layer-1](samp,i); //sum d_b across all samples
    }
    this->delta_b[final_layer-1](0,i) /= (double) samples;
    for(int j=0;j<layer;j++){
      this->delta_w[final_layer-1](j,i) = 0.0; // rezero d_w
      for(int samp=0;samp<samples;samp++){
	this->delta_w[final_layer-1](j,i) += this->delta[final_layer-1](samp,i)*activations[final_layer-1](samp,j); // sum d_w across all samples
      }
      this->delta_w[final_layer-1](j,i) /= (double) samples;
    }
  }
  
  double del=0.0;
  for(int ilay=this->layers-2;ilay>0;ilay--){
    layer=this->Layers[ilay-1];
    n_layer=this->Layers[ilay];
    this->delta[ilay-1] = Eigen::MatrixXd::Zero(samples,n_layer); //rezero delta
    act_der = this->Activation->Derivative(zs[ilay]);
    this->delta[ilay-1] = this->delta[ilay]*this->Weight[ilay].transpose();
    for(int i=0;i<n_layer;i++){ //calculate delta & d_b
      this->delta_b[ilay-1](0,i) = 0.0; //rezero d_b
      for(int samp=0;samp<samples;samp++){
	del=this->delta[ilay-1](samp,i);
	this->delta[ilay-1](samp,i) = del*act_der(samp,i);//coeff wise product d(l+1) = d(l+1)*f'(a(l))
	this->delta_b[ilay-1](0,i) += this->delta[ilay-1](samp,i); //sum d_b over all samples
      }
      this->delta_b[ilay-1](0,i) /= (double) samples;

      for(int j=0;j<layer;j++){
	this->delta_w[ilay-1](j,i) = 0.0; //rezero d_w
	for(int samp=0;samp<samples;samp++){ //d_w_ji = a(j)*d_l(i)
	  this->delta_w[ilay-1](j,i) += this->delta[ilay-1](samp,i)*activations[ilay-1](samp,j);
	}
	this->delta_w[ilay-1](j,i) /= (double) samples;
      }
    }
  }
}
void CNeuralNet::Train(int Epochs){
  printf("Training NEURAL_NET emulator...\n");
  int samples=this->X.rows();
  int runs=ceil((double)samples/(double)this->mini_batch);
  int params = this->X.cols();
  int obs = this->Y.cols();
  int start, finish, batch = this->mini_batch,
    pos;
  int print_int=Epochs/10;
  std::vector<int> rand;
  for(int i=0;i<samples;i++){
    rand.push_back(i);
  }
  Eigen::MatrixXd x = Eigen::MatrixXd::Zero(batch,params);
  Eigen::MatrixXd y = Eigen::MatrixXd::Zero(batch,obs);
  printf("Starting Epoch %d of %d...",0,Epochs);
  printf("  Loss: %f\n",this->Loss->Function(this->Emulate(X),Y));
  for(int i=0;i<Epochs;i++){
    if((i+1)%print_int==0){
      printf("Starting Epoch %d of %d...",i+1,Epochs);
    }
    std::random_shuffle(rand.begin(),rand.end());
    for(int j=0;j<runs;j++){
      this->Solver->t += 1;
      start = j*mini_batch;
      finish = start + batch;
      if(start > samples-1) break;
      if(finish > samples){
	batch = samples - start;
	x = Eigen::MatrixXd::Zero(batch,params);
	y = Eigen::MatrixXd::Zero(batch,obs);
      }
      for(int k=0;k<batch;k++){
	pos = rand[j*this->mini_batch+k];
	x.row(k) = X.row(pos);
	y.row(k) = Y.row(pos);
      }
      /*
      start = rand[j]*mini_batch;
      finish = start + batch;
      if(start > samples-1) break;
      if(finish > samples) batch = samples - start;
      x = X.block(0,start,params,batch);
      y = Y.block(0,start,obs,batch);
      */
      if(this->regular_param != 0.0){
	for(int ilay=0;ilay<this->layers-1;ilay++){ //Update coeffs    
	  Weight[ilay] *= (1.0 - this->regular_param);
	}
      }
      this->FeedForward(x);
      this->BackPropagation(x,y);
      this->Solver->Solve(this->Weight,this->Bias,this->delta_w,this->delta_b);
    }
    if((i+1)%print_int==0){
      printf("  Loss: %f\n",this->Loss->Function(this->Emulate(X),Y));
    }
  }
}
Eigen::MatrixXd CNeuralNet::Max(Eigen::MatrixXd Z){
  int rows=Z.rows(),
    cols=Z.cols(),
    index=0;
  double max=0.0;
  Eigen::MatrixXd Max = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    max=Z(irow,0);
    index=0;
    for(int icol=0;icol<cols;icol++){
      if(Z(irow,icol) > max){
	max = Z(irow,icol);
	index = icol;
      }
    }
    Max(irow,index) = 1.0;
  }
  return Max;
}
Eigen::MatrixXd CNeuralNet::Emulate(Eigen::MatrixXd X){
  Eigen::MatrixXd 
    output=X,
    temp;
  int samples=X.rows();
  for(int ilayer=0;ilayer<this->layers-2;ilayer++){
    temp = output;
    output = temp*this->Weight[ilayer];
    for(int i=0;i<samples;i++){
      output.row(i) += this->Bias[ilayer];
    }
    temp = output;
    output = this->Activation->Function(temp);
  }
  temp = output;
  output = temp*this->Weight[this->layers-2];
  for(int i=0;i<samples;i++){
    output.row(i) += this->Bias[this->layers-2];
  }
  temp = this->Loss->Activation->Function(output);
  return temp;
}
double CNeuralNet::Accuracy(){
  Eigen::MatrixXd ytest = this->Max(this->Emulate(this->X));
  double acc=this->X.rows();
  for(int i=0;i<this->Y.rows();i++){
    for(int j=0;j<this->Y.cols();j++){
      acc -= abs(this->Y(i,j) - ytest(i,j))/2.0;
    }
  }
  return acc;
}
double CNeuralNet::Accuracy(Eigen::MatrixXd x, Eigen::MatrixXd y){
  Eigen::MatrixXd ytest = this->Max(this->Emulate(x));
  double acc=x.rows();
  for(int i=0;i<y.rows();i++){
    for(int j=0;j<y.cols();j++){
      acc -= abs(y(i,j) - ytest(i,j))/2.0;
    }
  }
  return acc;
}
