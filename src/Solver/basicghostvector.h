#ifndef  __BasicGhostVector_h
#define  __BasicGhostVector_h

#include  <string>
#include  <iostream>


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments BasicGhostVector

////
////
/////////////////////////////////////////////

class BasicGhostVector : public std::string
{
private:

  std::string _type;
  mutable int _ncomp;

protected:

  std::string& GetName() {return *this;}

public:

//
////  Con(De)structor 
//

  BasicGhostVector() : _type("node"),_ncomp(-1) {}
  BasicGhostVector(const std::string& name) : std::string(name), _type("node"),_ncomp(-1) {}
  BasicGhostVector(const std::string& name, const std::string& type) : std::string(name), _type(type),_ncomp(-1) {}
  BasicGhostVector(const BasicGhostVector& v) {
    SetName(v.GetName());
    SetType(v.GetType());
    _ncomp = v.GetNcomp();
  }
  virtual ~BasicGhostVector() {}

  void SetName(const std::string& name) {GetName()=name;}
  void SetType(const std::string& type) {_type=type;}
  void SetNcomp(int n) const { _ncomp=n; }

  const std::string& GetName() const {return *this;}
  const std::string& GetType() const {return _type;}
  int GetNcomp() const { return _ncomp;}

  friend std::ostream& operator<<(std::ostream& os, const BasicGhostVector& g) {
    os << "Name:\t" << g.GetName() << std::endl;
    os << "Type:\t" << g.GetType() << std::endl;
    return os;
  }

};
}

#endif
