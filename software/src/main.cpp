#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<stdlib.h>
#include "balancemodel.h"


int main(int argc, char* argv[]){
  if(argc == 2){
    int testZ = atoi(argv[1]);
    printf("Running CBalanceModel on line %d from ModelZ.",testZ);
    std::cout << std::endl;
    CBalanceModel Run("run.dat",testZ);
  }else{
    int testZ = -1;
    printf("Running CBalanceModel on EXP data.\n");
    CBalanceModel Run("run.dat",testZ);
  }
  return 0;
}
