#ifndef  __IntegrationFormulaFace_h
#define  __IntegrationFormulaFace_h

#include  "integrationformula.h"

/*-----------------------------------------*/

class IntegrationFormulaEdge1 : public IntegrationFormula2d
{
protected:
 public:  
  IntegrationFormulaEdge1() : IntegrationFormula2d(4)
    {
      iw = 0.25;
      ic[0].x() = 0.5;  ic[0].y() = 0.0;
      ic[1].x() = 1.0;  ic[1].y() = 0.5;
      ic[2].x() = 0.5;  ic[2].y() = 1.0;
      ic[3].x() = 0.0;  ic[3].y() = 0.5;
    }

};


#endif
