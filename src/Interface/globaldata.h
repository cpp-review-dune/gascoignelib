#ifndef  __GlobalData_h
#define  __GlobalData_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GlobalData

////
////
/////////////////////////////////////////////

#include "gascoigne.h"

using namespace std;
using namespace Gascoigne;

class GlobalData
{
public:

  GlobalNodeData      _node;
  GlobalCellData      _cell;
  GlobalParameterData _parameter;

//
////  Con(De)structor 
//

  GlobalData() {}
  ~GlobalData() {}

  void AddNodeVector(const GlobalVector* d) {
    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(nd!=NULL);
    _node.insert(nd);
  }
  void AddCellVector(const GlobalVector* d) {
    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(nd!=NULL);
    _cell.insert(nd);
  }
  void AddParameterVector(const GlobalVector* d) {
    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(nd!=NULL);
    _parameter.insert(nd);
  }
  void DeleteNodeVector(const GlobalVector* d) {
    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(nd!=NULL);
    _node.erase(nd);
  }
  void DeleteCellVector(const GlobalVector* d) {
    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(nd!=NULL);
    _cell.erase(nd);
  }
  void DeleteParameterVector(const GlobalVector* d) {
    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(nd!=NULL);
    _parameter.erase(nd);
  }

  const GlobalNodeData& GetNodeData() const {return _node;}
  const GlobalCellData& GetCellData() const {return _cell;}
  const GlobalParameterData& GetParameterData() const {return _parameter;}
};


#endif
