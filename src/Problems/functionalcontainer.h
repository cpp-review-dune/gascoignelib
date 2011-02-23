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
             abort();
	   }
	 (*this)[label]=P;
       }
     
     void RemoveFunctional(const std::string& label)
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
             abort();
	   }
	 this->erase(label);
       }
     
     const Functional* GetFunctional(const std::string& label) const
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Functional " << label << " not present!\n";
             abort();
	   }
	 return find(label)->second;
       }

     int GetIndex(const std::string& label) const
     {
       int i = 0;
       for (FunctionalContainer::const_iterator it = this->begin(); it!=this->end();++it,++i)
	 {
	   if(it->first == label)
	     return i;
	 }
       std::cerr<<"Label not found"<<std::endl;
       abort();
     }
     
   };
 
}


/*----------------------------   functionalcontainer.h     ---------------------------*/
/* end of #ifndef __functionalcontainer_H */
#endif
/*----------------------------   functionalcontainer.h     ---------------------------*/
