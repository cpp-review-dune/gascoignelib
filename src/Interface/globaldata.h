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

namespace Gascoigne{
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

  void AddNodeVector(const std::string& name, const GlobalVector* d) {
//    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(d!=NULL);
    if(!_node.insert(std::make_pair(name,d)).second)
    {
      std::cerr << "NodeVector \"" << name << "\" already added" << std::endl;
    }
  }
  
  void AddCellVector(const std::string& name, const GlobalCellVector* d) {
//    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(d!=NULL);
    if(!_cell.insert(std::make_pair(name,d)).second)
    {
      std::cerr << "CellVector \"" << name << "\" already added" << std::endl;
    }
  }
  
  void AddParameterVector(const std::string& name, const GlobalParameterVector* d) {
//    const GlobalVector* nd = dynamic_cast<const GlobalVector*>(d);
    assert(d!=NULL);
    if(!_parameter.insert(std::make_pair(name,d)).second)
    {
      std::cerr << "ParameterVector \"" << name << "\" already added" << std::endl;
    }
  }
  
  void DeleteNodeVector(const std::string& name) {
    if(!_node.erase(name))
    {
      std::cerr << "NodeVector \"" << name << "\" cannot be deleted" << std::endl;
    }
  }
  
  void DeleteCellVector(const std::string& name) {
    if(!_cell.erase(name))
    {
      std::cerr << "NodeVector \"" << name << "\" cannot be deleted" << std::endl;
    }
  }
  
  void DeleteParameterVector(const std::string& name) {
    if(!_parameter.erase(name))
    {
      std::cerr << "NodeVector \"" << name << "\" cannot be deleted" << std::endl;
    }
  }

  const GlobalNodeData& GetNodeData() const {return _node;}
  const GlobalCellData& GetCellData() const {return _cell;}
  const GlobalParameterData& GetParameterData() const {return _parameter;}
};
}

#endif
