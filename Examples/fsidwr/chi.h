/*----------------------------   chi.h     ---------------------------*/
/*      $Id: chi.h,v 1.11 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __chi_H
#define __chi_H
/*----------------------------   chi.h     ---------------------------*/

#include "vertex.h"

namespace Gascoigne
{

  
  class Chi
  {
    std::string __type;
    
  public:


    void BasicInit(std::string X)
      {
	__type = X;
      }
    
    //  0: interface
    // -1: fluid
    //  1: solid

    int operator()(const Vertex3d& v) const;
    int operator()(const Vertex2d& v) const;
  };
  
  
  
}



/*----------------------------   chi.h     ---------------------------*/
/* end of #ifndef __chi_H */
#endif
/*----------------------------   chi.h     ---------------------------*/
