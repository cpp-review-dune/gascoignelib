#ifndef  __BasicGhostVector_h
#define  __BasicGhostVector_h


/////////////////////////////////////////////
////
////@brief
////  ... comments BasicGhostVector

////
////
/////////////////////////////////////////////

#include  <string>
#include  <iostream>

using namespace std;

class BasicGhostVector : public string
{
private:

  string _type;

protected:

  string& GetName() {return *this;}

public:


//
////  Con(De)structor 
//

  BasicGhostVector() : _type("node") {}
  BasicGhostVector(const string& name) : string(name), _type("node") {}
  BasicGhostVector(const string& name, const string& type) : string(name), _type(type) {}
  BasicGhostVector(const BasicGhostVector& v) {
    SetName(v.GetName());
    SetType(v.GetType());
  }
  virtual ~BasicGhostVector() {}

  void SetName(const string& name) {GetName()=name;}
  void SetType(const string& type) {_type=type;}

  const string& GetName() const {return *this;}
  const string& GetType() const {return _type;}

  friend ostream& operator<<(ostream& os, const BasicGhostVector& g) {
    os << "Name:\t" << g.GetName() << endl;
    os << "Type:\t" << g.GetType() << endl;
    return os;
  }

};


#endif
