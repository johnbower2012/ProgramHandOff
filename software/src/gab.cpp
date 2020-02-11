#include<iostream>
#include<string>
#include "system.h"

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      printf("Usage: ./gab filename\n");
      exit(1);
    }
  std::string 
    folder="model_output",
    filename=argv[1],
    delimiter=" ";
  int ab=4;

  WriteGABFunctions(folder+"/"+filename,delimiter,ab);
}
