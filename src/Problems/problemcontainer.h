/*----------------------------   problemcontainer.h     ---------------------------*/
/*      $Id$                 */
#ifndef __problemcontainer_H
#define __problemcontainer_H
/*----------------------------   problemcontainer.h     ---------------------------*/


#include <map>
#include <string> 
#include "problemdescriptorinterface.h"

namespace Gascoigne
{
  
 class  ProblemContainer : public std::map<std::string, const ProblemDescriptorInterface*>
   {
     public:

     void AddProblem(const std::string& label, const ProblemDescriptorInterface* P)
       {
	 if (find(label)!=end())
	   {
	     std::cerr << "Problemdescriptor " << label << " already present!\n";
	     assert(0);
	   }
	 (*this)[label]=P;
       }
     
     void RemoveProblem(const std::string& label)
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
	     assert(0);
	   }
	 this->erase(label);
       }
     
     const ProblemDescriptorInterface* GetProblem(const std::string& label) const
       {
	 if (find(label)==end())
	   {
	     std::cerr << "Problemdescriptor " << label << " not present!\n";
	     assert(0);
	   }
	 return find(label)->second;
       }
     
   };
 
}



/*----------------------------   problemcontainer.h     ---------------------------*/
/* end of #ifndef __problemcontainer_H */
#endif
/*----------------------------   problemcontainer.h     ---------------------------*/
