#ifndef  __findcheck_h
#define  __findcheck_h

/*------------------------------------------*/

#define FindMacro(where,message) \
    { CHashIt ip = ##where .find(i); \
       if(ip==##where .end()) { \
          std::string s = "##message ";\
	  cerr << "Index::" << s << " g2l()\n";\
	  cerr << "there is no " << s << " "<< i << std::endl;\
	  abort();}\
      return ip->second; }

#define CheckMacro(where) \
    { CHashIt ip = ##where .find(i);\
      if(ip==##where .end()) return -2;\
      return ip->second;}

/*------------------------------------------*/

#endif
