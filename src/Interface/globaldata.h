#ifndef  __GlobalData_h
#define  __GlobalData_h


#include "gascoigne.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments GlobalData

  ////
  ////
  /////////////////////////////////////////////

  class GlobalData
  {
    private:

    protected:

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
        assert(d!=NULL);
        if(!_node.insert(std::make_pair(name,d)).second)
        {
          std::cerr << "NodeVector \"" << name << "\" already added" << std::endl;
          abort();
        }
      }
      void AddCellVector(const std::string& name, const GlobalCellVector* d) {
        assert(d!=NULL);
        if(!_cell.insert(std::make_pair(name,d)).second)
        {
          std::cerr << "CellVector \"" << name << "\" already added" << std::endl;
          abort();
        }
      }
      void AddParameterVector(const std::string& name, const GlobalParameterVector* d) {
        assert(d!=NULL);
        if(!_parameter.insert(std::make_pair(name,d)).second)
        {
          std::cerr << "ParameterVector \"" << name << "\" already added" << std::endl;
          abort();
        }
      }
      
      void DeleteNodeVector(const std::string& name) {
        _node.erase(name);
      }
      void DeleteCellVector(const std::string& name) {
        _cell.erase(name);
      }
      void DeleteParameterVector(const std::string& name) {
        _parameter.erase(name);
      }

      const GlobalNodeData& GetNodeData() const {
        return _node;
      }
      const GlobalCellData& GetCellData() const {
        return _cell;
      }
      const GlobalParameterData& GetParameterData() const {
        return _parameter;
      }
  };
}

#endif
