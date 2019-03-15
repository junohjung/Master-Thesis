#include "zfsvariables.h"

#include "zfsalloc.h"

namespace zfs { namespace fv {

////////////////////////////////////////////////////////////////////////////////
/// ConservativeVariables:

template<ZFSId nd> const ZFSId ConservativeVariables<nd>::RHO_VV[3] = { 0, 1, 2 };

template<ZFSId nd>
ConservativeVariables<nd>::ConservativeVariables(const ZFSId noSpecies)
    : m_noSpecies(noSpecies), noVariables(nd + 2 + noSpecies)  {
  zfsAlloc(RHO_Y,m_noSpecies,"zfs::fv::ConservativeVariables::RHO_Y",__CALLING_FUNCTION__ );
  for( ZFSId i = 0; i < m_noSpecies; ++i ) {
    RHO_Y[i] = RHO_Z+i;
  }
}

template<ZFSId nd>
ConservativeVariables<nd>::~ConservativeVariables() {
  zfsDeallocate(RHO_Y);
}

template struct ConservativeVariables<2>;
template struct ConservativeVariables<3>;

////////////////////////////////////////////////////////////////////////////////
/// PrimitiveVariables:

template<ZFSId nd> const ZFSId PrimitiveVariables<nd>::VV[3] = { 0, 1, 2 };

template<ZFSId nd>
PrimitiveVariables<nd>::PrimitiveVariables( const ZFSId noSpecies )
    : m_noSpecies(noSpecies), noVariables(nd + 2 + noSpecies) {
  zfsAlloc(Y,m_noSpecies,"zfs::fv::PrimitiveVariables::Y",__CALLING_FUNCTION__);
  for( ZFSId i = 0; i < m_noSpecies; ++i ) {
    Y[i] = Z + i;
  }
}

template<ZFSId nd>
PrimitiveVariables<nd>::~PrimitiveVariables() { zfsDeallocate(Y); }

template struct PrimitiveVariables<2>;
template struct PrimitiveVariables<3>;

} // namespace fv
} // namespace zfs
