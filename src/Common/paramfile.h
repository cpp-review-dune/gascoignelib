#ifndef  __ParamFile_h
#define  __ParamFile_h


#include  <string>

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments ParamFile

////
////
/////////////////////////////////////////////

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
