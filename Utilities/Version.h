#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "799b34d7a7f9cb917ca819f0824a680f6b2cc901"
#define GIT_REFSPEC "refs/heads/master"
#define GIT_LOCAL_STATUS "DIRTY"

#define SPLISHSPLASH_VERSION "2.12.4"

#ifdef DL_OUTPUT
#pragma WARNING(Local changes not committed.)
#endif

#endif
