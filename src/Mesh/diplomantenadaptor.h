#ifndef  __DiplomantenAdaptor_h
#define  __DiplomantenAdaptor_h

#include "nvector.h"
#include "adaptordata.h"

/*-----------------------------------------*/

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

/*-----------------------------------------*/

class DiplomantenAdaptor
{
protected:

  typedef nvector<double>  dvector;
  typedef nvector<int>     ivector;

  AdaptorData&             info;
  const nvector<double>&   eta;
  int ppp;

  void analyse() const;

public:

  DiplomantenAdaptor(AdaptorData&, const dvector& eta);
  void refine(ivector& ref);
  void MalteRefine(ivector& ref) const;
};

#endif
