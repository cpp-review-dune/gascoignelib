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


class GlobalData
{
public:

  Gascoigne::GlobalNodeData      _node;
  Gascoigne::GlobalCellData      _cell;
  Gascoigne::GlobalParameterData _parameter;

//
////  Con(De)structor 
//

  GlobalData() {}
  ~GlobalData() {}

  void AddNodeVector(const Gascoigne::GlobalVector* d) {
    const Gascoigne::GlobalVector* nd = dynamic_cast<const Gascoigne::GlobalVector*>(d);
    assert(nd!=NULL);
    _node.insert(nd);
  }
  void AddCellVector(const Gascoigne::GlobalVector* d) {
    const Gascoigne::GlobalVector* nd = dynamic_cast<const Gascoigne::GlobalVector*>(d);
    assert(nd!=NULL);
    _cell.insert(nd);
  }
  void AddParameterVector(const Gascoigne::GlobalVector* d) {
    const Gascoigne::GlobalVector* nd = dynamic_cast<const Gascoigne::GlobalVector*>(d);
    assert(nd!=NULL);
    _parameter.insert(nd);
  }
  void DeleteNodeVector(const Gascoigne::GlobalVector* d) {
    const Gascoigne::GlobalVector* nd = dynamic_cast<const Gascoigne::GlobalVector*>(d);
    assert(nd!=NULL);
    _node.erase(nd);
  }
  void DeleteCellVector(const Gascoigne::GlobalVector* d) {
    const Gascoigne::GlobalVector* nd = dynamic_cast<const Gascoigne::GlobalVector*>(d);
    assert(nd!=NULL);
    _cell.erase(nd);
  }
  void DeleteParameterVector(const Gascoigne::GlobalVector* d) {
    const Gascoigne::GlobalVector* nd = dynamic_cast<const Gascoigne::GlobalVector*>(d);
    assert(nd!=NULL);
    _parameter.erase(nd);
  }

  const Gascoigne::GlobalNodeData& GetNodeData() const {return _node;}
  const Gascoigne::GlobalCellData& GetCellData() const {return _cell;}
  const Gascoigne::GlobalParameterData& GetParameterData() const {return _parameter;}
};


#endif
