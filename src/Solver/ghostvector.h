#ifndef __GhostVector_h
#define __GhostVector_h

#include  <string>
#include  <cassert>
#include  <iostream>

/*---------------------------------------*/

class GhostVector : public std::string
{
 private:

 public:

  GhostVector(std::string s):std::string(s){}

  GhostVector& operator=(const GhostVector& v) {
    std::cerr << "GhostVector::operator=" << std::endl;
    assert(0);
  }
  void set(const GhostVector& v) {
    std::string::operator=(v);
  }
};

  
#endif  
