/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "stdloop.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include "vertex.h"
#include  "gascoignemesh2d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"


#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HASHMAP   std::tr1::unordered_map
#else
#include  <ext/hash_map>
#define HASHMAP  __gnu_cxx::hash_map
#endif

extern double BOUNDARY;


namespace Gascoigne
{


  template<int DIM>
    class Loop : public StdLoop
    {

    public:

      void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,const FunctionalContainer* FC)
      {
 //   	GetMultiLevelSolverPointer() = new MyMLS;

    	StdLoop::BasicInit(paramfile,PC,FC);
      }

      void run(const std::string& problemlabel);
    };

}


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
