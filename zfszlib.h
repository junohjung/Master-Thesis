#ifndef ZFS_ZLIB_H_
#define ZFS_ZLIB_H_

////////////////////////////////////////////////////////////////////////////////
/// \file This file wraps the zlib header file
///
/// Note: this allows disabling the use of zlib
////////////////////////////////////////////////////////////////////////////////

#ifdef WITH_ZLIB
#include "zlib.h"
#endif

#endif // ZFS_ZLIB_H_
