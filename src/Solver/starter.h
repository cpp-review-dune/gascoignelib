#ifndef  __Starter_h
#define  __Starter_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Starter

////
////
/////////////////////////////////////////////

#include  <string>
#include  "paramfile.h"

class Starter
{
private:

  bool   _erase;
  Gascoigne::ParamFile* _paramfile;

public:

//
////  Con(De)structor 
//

  Starter(int argc, char** argv, const std::string& paramfile);

  ~Starter();

  const Gascoigne::ParamFile* GetParamFile() const {return _paramfile;}
};


#endif
