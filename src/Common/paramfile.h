#ifndef  __ParamFile_h
#define  __ParamFile_h


/////////////////////////////////////////////
////
////@brief
////  ... comments ParamFile

////
////
/////////////////////////////////////////////

#include  <string>

namespace Gascoigne
{

  class ParamFile : public std::string
{
private:


protected:


public:


//
////  Con(De)structor 
//
  
  ParamFile() : std::string() {}
  ParamFile(const std::string& name) : std::string(name) {}
  std::string GetName() const {return *this;}
  void SetName(const std::string& name) { std::string::operator=(name);} 
};
}

#endif
