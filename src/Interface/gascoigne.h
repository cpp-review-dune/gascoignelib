#ifndef  __Gascoigne_h
#define  __Gascoigne_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Gascoigne

////
////
/////////////////////////////////////////////

#include  <set>
#include  <map>
#include  <string>
#include  "compvector.h"
#include  "derivativevector.h"

namespace Gascoigne
{

  typedef CompVector<double>         GlobalVector;  
  typedef CompVector<double>         LocalVector;  

  typedef nvector<int>               IntVector;  
  typedef nvector<double>            DoubleVector;  

  typedef std::set<int>              IntSet;  

//   typedef std::map<const std::string,const NodeVector*>  GlobalNodeData;
//   typedef std::map<const std::string,const CellVector*>  GlobalCellData;
//   typedef std::map<const std::string,const ParameterVector*>  GlobalParameterData;
  typedef std::map<const std::string,const GlobalVector*>  GlobalNodeData;
  typedef std::map<const std::string,const GlobalVector*>  GlobalCellData;
  typedef std::map<const std::string,const GlobalVector*>  GlobalParameterData;

  typedef std::map<const std::string,LocalVector>       LocalData;
  typedef CompVector<double>::iterator                  VectorIterator;
  typedef std::vector<DerivativeVector>                 FemFunction;
  typedef std::map<const std::string,FemFunction>       FemData;

  typedef  DerivativeVector  TestFunction;
};


#endif
