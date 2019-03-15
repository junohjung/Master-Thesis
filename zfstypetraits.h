#ifndef ZFSTYPETRAITS_H
#define ZFSTYPETRAITS_H

#include "zfsmpi.h"
#include "zfstypes.h"

namespace zfs
{

// Disallow usage of non-specified traits
namespace detail_ { template <class T> struct traits_error { }; }
template <class T> struct type_traits
{
  typedef typename detail_::traits_error<T>::ERROR_CANNOT_USE_NON_SPECIFIED_TRAITS error;
};

// Specialized templates for known types

template <> struct type_traits<ZFSFloat>
{
  static ZFSString name() { return "ZFSFloat"; };
  static MPI_Datatype mpiType() { return MPI_DOUBLE; };
};

template <> struct type_traits<ZFSString>
{
  static ZFSString name() { return "ZFSString"; };
  // ZFSString has no corresponding MPI data type
};

template <> struct type_traits<ZFSChar>
{
  static ZFSString name() { return "ZFSChar"; };
  static MPI_Datatype mpiType() { return MPI_CHAR; };
};

template <> struct type_traits<ZFSInt>
{
  static ZFSString name() { return "ZFSInt"; };
  static MPI_Datatype mpiType() { return MPI_INT; };
};

template <> struct type_traits<ZFSUint>
{
  static ZFSString name() { return "ZFSUint"; };
  static MPI_Datatype mpiType() { return MPI_UNSIGNED; };
};

template <> struct type_traits<ZFSLong>
{
  static ZFSString name() { return "ZFSLong"; };
  static MPI_Datatype mpiType() { return MPI_LONG; };
};

template <> struct type_traits<ZFSUlong>
{
  static ZFSString name() { return "ZFSUlong"; };
  static MPI_Datatype mpiType() { return MPI_UNSIGNED_LONG; };
};

template <> struct type_traits<ZFSBool>
{
  static ZFSString name() { return "ZFSBool"; };
  // ZFSBool has no corresponding MPI data type
};

} /* namespace zfs */


#endif /* ifndef ZFSTYPETRAITS_H */
