#include  "integrationformula.h"

/*------------------------------------------------------------*/

IntegrationFormulaInterface::IntegrationFormulaInterface() 
{
}

/*------------------------------------------------------------*/

IntegrationFormulaInterface::IntegrationFormulaInterface(int n)
{
  init(n);
} 

/*------------------------------------------------------------*/

void IntegrationFormulaInterface::init(int n)
{
  in = n;
  iw.reserve(n);      
  iw.resize (n);
}

