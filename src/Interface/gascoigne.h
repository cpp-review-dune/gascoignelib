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
#include  "nmatrix.h"
#include  "derivativevector.h"

namespace Gascoigne
{
  typedef CompVector<double>                                        GlobalVector;
  typedef CompVector<double>                                        LocalVector;
  typedef std::map<const std::string,const GlobalVector*>           GlobalNodeData;
  typedef std::map<const std::string,LocalVector>                   LocalNodeData;

  typedef CompVector<double>                                        GlobalCellVector;
  typedef CompVector<double>                                        LocalCellVector;
  typedef std::map<const std::string,const GlobalCellVector*>       GlobalCellData;
  typedef std::map<const std::string,LocalCellVector>               LocalCellData;

  typedef nvector<double>                                           GlobalParameterVector;
  typedef nvector<double>                                           LocalParameterVector;
  typedef std::map<const std::string,const GlobalParameterVector*>  GlobalParameterData;
  typedef std::map<const std::string,LocalParameterVector>          LocalParameterData;

  typedef nvector<int>                                              IntVector;
  typedef nvector<double>                                           DoubleVector;
  typedef nmatrix<double>                                           DoubleMatrix;
  typedef std::set<int>                                             IntSet;
  typedef CompVector<double>::iterator                              VectorIterator;

  typedef nmatrix<double>                                           TimePattern;

  typedef DerivativeVector                                          TestFunction;
  typedef std::vector<TestFunction>                                 FemFunction;
  typedef std::map<const std::string,FemFunction>                   FemData;

  typedef nvector<double>                                           CellFunction;
  typedef std::map<const std::string,CellFunction>                  CellData;
}


#endif
