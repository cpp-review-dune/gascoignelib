#ifndef __newton_h
#define __newton_h

#include  "nlinfo.h"

template<class Operator, class Vector>
void newton(Operator& S, Vector& u, const Vector& f, Vector& r, Vector& w, NLInfo& info);

template<class Operator, class Vector>
void newnewton(Operator& S, Vector& u, const Vector& f, Vector& r, Vector& w, NLInfo& info);

#endif
