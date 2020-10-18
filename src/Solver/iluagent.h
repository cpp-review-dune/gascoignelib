
/*----------------------------   iluagent.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __iluagent_H
#define __iluagent_H
/*----------------------------   iluagent.h     ---------------------------*/


/**
 *
 * Copyright (C) 2020 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#include  "gascoigne.h"
#include  "matrix.h"
#include  "iluinterface.h"

namespace Gascoigne
{

  /////////////////////////////////////////////
  ////
  ////@brief
  ////  stores the ilu data objects belonging to the Ilu labels
  //// 
  /////////////////////////////////////////////

  class IluAgent : public std::map<Matrix,IluInterface*>
  {
  public:

    typedef std::map<Matrix,IluInterface*>::const_iterator const_iterator;
    typedef std::map<Matrix,IluInterface*>::iterator       iterator;

    //
    ////  Con(De)structor 
    //

    IluAgent();
    ~IluAgent();

    void Register(const Matrix& mat);
    void Delete  (Matrix& mat);

    IluInterface& operator()(const Matrix& g);

    friend std::ostream& operator<<(std::ostream& os, const IluAgent& gva) {
      int i=0,n=gva.size();
      os << "IluAgent: size=" << n << ", ";
      for (auto p=gva.begin(); p!=gva.end(); p++,i++){
	os << "Ilu("<<i<<")=('"<< p->first << "',"<< p->second <<")";
	if( i <n-1 ) os << ", "; else os << ". ";
      }
      return os;
    }

  };
}



/*----------------------------   iluagent.h     ---------------------------*/
/* end of #ifndef __iluagent_H */
#endif
/*----------------------------   iluagent.h     ---------------------------*/
