#ifndef  __PatchIntegrationFormula_h
#define  __PatchIntegrationFormula_h

#include  "integrationformula.h"

namespace Gascoigne
{

/*------------------------------------------------------------*/

template<int N, class INT>
class PatchFormula2d : public IntegrationFormula2d{
public: PatchFormula2d();};

/*------------------------------------------------------------*/

template<int N, class INT>
class PatchFormula3d : public IntegrationFormula3d{
public: PatchFormula3d();};

}

#endif
