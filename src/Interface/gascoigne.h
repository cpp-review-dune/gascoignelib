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
                                                                   
  typedef nvector<int>                                              IntVector;  
  typedef nvector<double>                                           DoubleVector;  
  typedef nmatrix<double>                                           DoubleMatrix;  
  typedef nmatrix<double>                                           TimePattern;
                                                                   
  typedef std::set<int>                                             IntSet;  

  typedef CompVector<double>::iterator                              VectorIterator;

  typedef  DerivativeVector                                         TestFunction;
  typedef std::vector<TestFunction>                                 FemFunction;
  
  typedef CompVector<double>                                        GlobalCellVector;
  typedef nvector<double>                                           GlobalParameterVector;
  typedef std::map<const std::string,const GlobalVector*>           GlobalNodeData;
  typedef std::map<const std::string,const GlobalCellVector*>       GlobalCellData;
  typedef std::map<const std::string,const GlobalParameterVector*>  GlobalParameterData;
  
  typedef CompVector<double>                                        LocalCellVector;
  typedef nvector<double>                                           LocalParameterVector;
  typedef std::map<const std::string,LocalVector>                   LocalNodeData;
  typedef std::map<const std::string,LocalCellVector>               LocalCellData;
  typedef std::map<const std::string,LocalParameterVector>          LocalParameterData;
                                                                   
  typedef std::map<const std::string,FemFunction>                   FemData;
}


#endif
