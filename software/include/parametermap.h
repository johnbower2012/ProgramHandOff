#ifndef __CParametermap_H__
#define __CParametermap_H__

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
//#include "random.h"

using namespace std;

//---------------------------------------------------------------
// This code helps one parse the map that contains configuration 
// parameters.
//
// If "nameOfParameter is not a key in the map, -1 is returned since that
// is the 3rd argment of the getI function.
//
//MH 22 jun04
//---------------------------------------------------------------

//This code only works with a map of the type below.  The type def is
//to make it easy to remember.

//These functions are all in the namespace parameter.
class CParameterMap : public map<string,string> {
 public:
  bool   getB(string ,bool);
  int    getI(string ,int);
  string getS(string ,string);
  double getD(string ,double);
  vector< double > getV(string, string);
  vector< string > getVS(string, string);
  vector< vector< double > > getM(string, double);
  void set(string, double);
  void set(string, int);
  void set(string, bool);
  void set(string, string);
  void set(string, char*);
  void set(string, vector< double >);
  void set(string, vector< string >);
  void set(string, vector< vector< double > >);
  void ReadParsFromFile(const char *filename);
  void ReadParsFromFile(string filename);
  void PrintPars();
};

#endif
