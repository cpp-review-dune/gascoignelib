#ifndef  __VectorInterface_h
#define  __VectorInterface_h

#include  <string>
#include  <iostream>


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments VectorInterface
////
////
/////////////////////////////////////////////

class VectorInterface : public std::string
{
protected:

  std::string& GetName() {return *this;}

public:

//
////  Con(De)structor 
//

  //VectorInterface()  {}
  VectorInterface(const std::string& name) : std::string(name)  {}
  VectorInterface(const VectorInterface& v) {
    SetName(v.GetName());
  }
  virtual ~VectorInterface() {}

  void SetName(const std::string& name) {GetName()=name;}

  const std::string& GetName() const {return *this;}

  friend std::ostream& operator<<(std::ostream& os, const VectorInterface& g) {
    os << "Name:\t" << g.GetName() << std::endl;
    return os;
  }

};
}

#endif
