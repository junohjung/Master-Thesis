#ifndef ZFSTYPES_H
#define ZFSTYPES_H
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This files defines the types used in ZFS
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include "zfscompiler.h"
////////////////////////////////////////////////////////////////////////////////
// Note: If you change a type, make sure to check if zfstypetraits.h needs to be
//       updated as well.
typedef double                  ZFSFloat;
typedef std::basic_string<char> ZFSString;
typedef char                    ZFSChar;
typedef int                     ZFSId;
typedef int                     ZFSInt;
#if defined(ZFS_MS_COMPILER)
typedef __int64                ZFSLong;
typedef unsigned __int64       ZFSUlong;
#else
typedef long int                ZFSLong;
typedef unsigned long int       ZFSUlong;
#endif
typedef long long int           ZFSLonglong;
typedef bool                    ZFSBool;
typedef unsigned int            ZFSUint;

typedef unsigned char           ZFSUchar;
typedef unsigned short          ZFSUshort;
typedef short                   ZFSShort;
////////////////////////////////////////////////////////////////////////////////
//! define array structures
typedef struct z {
    ZFSString name;
    ZFSId id;
    ZFSId *blocks;
    ZFSId noBlocks;
} ZFSZone;
////////////////////////////////////////////////////////////////////////////////
#endif
