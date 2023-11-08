#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "7282f8511f4cce76293089e1005f490667d3d34e"
#define GIT_REFSPEC "refs/heads/master"
#define GIT_LOCAL_STATUS "CLEAN"

#define SPLISHSPLASH_VERSION "2.12.4"

#ifdef DL_OUTPUT

#endif

#endif
