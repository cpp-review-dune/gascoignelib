#ifndef __curvedshapes_h
#define __curvedshapes_h

#include <map>
#include <set>
#include <vector>
#include  <string>
#include "boundaryfunction.h"
#include "boundaryline.h"
#include "boundaryquad.h"
#include "vertex.h"

/******************************************************/

namespace Gascoigne
{
template<int DIM>
class CurvedShapes : public std::map<int,BoundaryFunction<DIM>* >
{
 public:

  typedef typename std::map<int,BoundaryFunction<DIM>* >::iterator iterator;

  ~CurvedShapes() {
    for (iterator p=std::map<int,BoundaryFunction<DIM>* >::begin();p!=std::map<int,BoundaryFunction<DIM>* >::end();++p)
    if (p->second) {
      std::cerr<< "not deleting shape: "<< p->second->GetName() << std::endl;
    }
  }

  const BoundaryFunction<DIM>& GetShape(int col) const { return *this->find(col)->second;}

  void AddShape(int col, BoundaryFunction<DIM>* f) {
    (*this)[col] = f;
  }
  
  void newton(int col,Vertex<DIM>& V) { GetShape(col).newton(V); }

    int Curved(int col) const { return (this->find(col)!=std::map<int,BoundaryFunction<DIM>* >::end());}

    bool empty() const {return (std::map<int,BoundaryFunction<DIM>* >::size()==0);}
};
}

/******************************************************/

#endif
