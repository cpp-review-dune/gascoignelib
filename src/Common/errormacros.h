#ifndef  __errormacros_h
#define  __errormacros_h

#define FILE_ERROR(FILE,NAME) \
if(!FILE.is_open()) { std::cerr << " Could not open file " << NAME << std::endl;;}

#endif
