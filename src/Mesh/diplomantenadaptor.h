#ifndef  __DiplomandenAdaptor_h
#define  __DiplomandenAdaptor_h

#include "gascoigne.h"
#include "adaptordata.h"

/*-----------------------------------------*/

namespace Gascoigne
{

/*
* Versucht den Fehler im neuen Gitter multipliziert
* mit der Zahl der Zellen zu minimieren:
*
* min:    E(neu) * N(neu)^(alpha)
*
* Dabei ist alpha = global_conv aus der Parameterdatei,
* die erreichbare Konvergenzordnung bzgl. N.
* Die Fehlerindikatoren sollen lokal mit der Konvergenz-
* ordnung beta = local_conv bzgl. h konvergieren.
*
* Der Parameter rfactor gibt zusaetzlich an, wieviele
* Zellen verfeinert werden sollen.
* rfactor = 1 bedeutet keine Nebenbedingung
* rfactor = 0.3 : Es werden mindestens die Zellen verfeinert,
*                 deren Fehlerindikatoren zusammen 70% betragen.
*/

class DiplomandenAdaptor
{
protected:

  AdaptorData&             info;
  const DoubleVector&   eta;
  int ppp;

  void analyse() const;

public:

  DiplomandenAdaptor(AdaptorData&, const DoubleVector& eta);
  void refine(IntVector& ref);
  void MalteRefine(IntVector& ref) const;
};
}

#endif
