/*----------------------------   functionalcontainer.h     ---------------------------*/
/*      $Id$                 */
#ifndef __functionalcontainer_H
#define __functionalcontainer_H
/*----------------------------   functionalcontainer.h     ---------------------------*/


#include <map>
#include <string> 
#include "functional.h"

namespace Gascoigne
{
  
 class  FunctionalContainer : public std::map<std::string, const Functional*>
   {
     public:

     void AddFunctional(const std::string& label, const Functional* P)
       {
	 if (find(label)!=end())
	   {
	     std::cerr << "Functional " << label << " already present!\n";
	     assert(0);
	   }
	 (*this)[label]=P;
       }
     
     void RemoveFunctional(const std::string& label)
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
	     assert(0);
	   }
	 this->erase(label);
       }
     
     const Functional* GetFunctional(const std::string& label) const
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Functional " << label << " not present!\n";
	     assert(0);
	   }
	 return find(label)->second;
       }
     
   };
 
}


/*----------------------------   functionalcontainer.h     ---------------------------*/
/* end of #ifndef __functionalcontainer_H */
#endif
/*----------------------------   functionalcontainer.h     ---------------------------*/
