#ifndef  __Starter_h
#define  __Starter_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Starter

////
////
/////////////////////////////////////////////

#include <string>

class Starter
{
private:

  bool   _erase;
  std::string _paramfile;

public:

//
////  Con(De)structor 
//

  Starter(int argc, char** argv, const std::string& paramfile);

  ~Starter();

  std::string GetParamFile() const {return _paramfile;}
};


#endif
