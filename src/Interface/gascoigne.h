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
#include  "compvector.h"
#include  "derivativevector.h"

namespace Gascoigne
{

  typedef CompVector<double>         GlobalVector;  
  typedef CompVector<double>         LocalVector;  

  typedef nvector<int>               IntVector;  
  typedef nvector<double>            DoubleVector;  

  typedef std::set<int>              IntSet;  

//   typedef std::set<const NodeVector*>  GlobalNodeData;
//   typedef std::set<const CellVector*>  GlobalCellData;
//   typedef std::set<const ParameterVector*>  GlobalParameterData;
  typedef std::set<const GlobalVector*>  GlobalNodeData;
  typedef std::set<const GlobalVector*>  GlobalCellData;
  typedef std::set<const GlobalVector*>  GlobalParameterData;

  typedef std::vector<LocalVector>       LocalData;
  typedef CompVector<double>::iterator   VectorIterator;
  typedef std::vector<DerivativeVector>  FemFunction;
  typedef std::vector<FemFunction>       FemData;

  typedef  DerivativeVector  TestFunction;
};


#endif
