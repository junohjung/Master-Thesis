#ifndef ZFS_VARIABLES_H_
#define ZFS_VARIABLES_H_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file defines the index classes for the access of within
/// ZFS's arrays of variables (e.g. Cell.m_variables, block.m_variables,...)
////////////////////////////////////////////////////////////////////////////////
#include "zfstypes.h"
////////////////////////////////////////////////////////////////////////////////

/*!\class ZFSConservativeVariables
   \date begin: 00.05.19 change 00.05.19
   \\author changed by Pascal Meysonnat to deal also with RANS

   \brief Storage of the Position of the Conservative Variables (RHO, RHO_VV, RHO_E)
   in the value vectors of the blocks ans surfaces

   If you search for a value in a value vector you only have to type in the name
   of the value you are searching (e.g. variables[RHO_E])
   \warning Be careful! Some value vectors include the conservative, others the
   primitive Variables
*/
template <ZFSInt nDim>
class ZFSConservativeVariables
{
public:
  //! Sets the position for the conservative variables
  ZFSConservativeVariables(ZFSId noSpecies, ZFSId noRANSEq=0 )
  {
    m_noSpecies = noSpecies;
    m_noRansEquations=noRANSEq;
    RHO_VV = new ZFSId[ nDim ];
    for( ZFSId i = 0;  i < nDim;  i++ ) RHO_VV[ i ] = RHO_U + i;
    RHO_Y = new ZFSId[ m_noSpecies ];
    for( ZFSId i=0; i<m_noSpecies; i++ ) RHO_Y[i] = RHO_Z+i;
    //pascalm< 
    RANS_FIRST=nDim +2; 
    if(m_noSpecies !=0){
      RANS_FIRST= nDim+2+m_noSpecies;
    }
    RANS_VAR = new ZFSId[m_noRansEquations];
    for(ZFSId i=0; i<m_noRansEquations; ++i) RANS_VAR[i]=RANS_FIRST+i;
    //>pascalm 
    noVariables = nDim + 2 + m_noSpecies + m_noRansEquations;

    rhoVVInfinity = new ZFSFloat [ nDim ];
    ransInfinity = new ZFSFloat [ nDim ];
  }
  ~ZFSConservativeVariables()
  {
    delete[] RHO_VV;
    delete[] RHO_Y;
    delete[] rhoVVInfinity;
    delete[] RANS_VAR;
    delete[] ransInfinity;
  }

  //! Position of RHO_U
  static const constexpr ZFSId   RHO_U = 0;
  //! Position of RHO_V
  static const constexpr ZFSId   RHO_V = 1;
  //! Position of RHO_W
  static const constexpr ZFSId   RHO_W = 2;
  //! Pointer for the velocities so you can use them in a loop
  ZFSId * RHO_VV;
  //! Position of RHO_E
  static const constexpr ZFSId   RHO_E = nDim;
  //! Position of RHO
  static const constexpr ZFSId   RHO = nDim + 1;
  //! Position of RHO_Z
  static const constexpr ZFSId   RHO_Z = nDim + 2;
  //! Position of RHO_C
  static const constexpr ZFSId   RHO_C = nDim + 2;
  //! Position of RHO_Yi
  ZFSId*  RHO_Y;
  //! first Position of RANS Variables
  ZFSId   RANS_FIRST;
  //! Position of RANS Variables
  ZFSId*  RANS_VAR;


  //! The Nr. of Conservative Variables (nDim + 2)
  ZFSId   noVariables;
  ZFSId   m_noSpecies;
  ZFSId   m_noRansEquations;
  ZFSFloat rhoUInfinity, rhoVInfinity, rhoWInfinity, rhoEInfinity, rhoInfinity, rhoInfinityGhost;
  ZFSFloat * ransInfinity;
  ZFSFloat * rhoVVInfinity;
};

/*!\class ZFSPrimitiveVariables
   \date begin: 00.05.19 change 00.05.19
   \brief Storage of the Position of the Primitive Variables (u, v, w, T, p) in
   the value vectors of the blocks ans surfaces

   If you search for a value in a value vector you only have to type in the name
   of the value you are searching (e.g. variables[T])
   \warning Be careful! Some value vectors include the conservative, others the
   primitive Variables
*/
template <ZFSInt nDim>
class ZFSPrimitiveVariables
{
public:
  //! Sets the position for the primitive variables
  ZFSPrimitiveVariables(ZFSId noSpecies, ZFSId noRANSEq=0 ) {
    m_noSpecies = noSpecies;
    m_noRansEquations=noRANSEq;

    VV = new ZFSId[ nDim ];
    for( ZFSId i = 0;  i < nDim;  i++ ) VV[ i ] = U + i;
    Y = new ZFSId[ m_noSpecies ];
    for( ZFSId i=0; i<m_noSpecies; i++ ) Y[i] = Z+i;
    //G = nDim+1;

    RANS_FIRST=nDim + 2; 
    if(m_noSpecies !=0){
      RANS_FIRST= nDim+2+m_noSpecies;
    }
    RANS_VAR = new ZFSId[m_noRansEquations];
    for(ZFSId i=0; i<m_noRansEquations; ++i) RANS_VAR[i]=RANS_FIRST+i;

    noGCellVariables = nDim + 1;
    noVariables = nDim + 2 + m_noSpecies + m_noRansEquations;

    VVInfinity = new ZFSFloat [ nDim ];
    ransInfinity = new ZFSFloat [ nDim ];
  }
  ~ZFSPrimitiveVariables() {
    delete[] VV;
    delete[] Y;
    delete[] VVInfinity;
    delete[] RANS_VAR;
    delete[] ransInfinity;
  }

  //! Position of U
  static const constexpr ZFSId   U = 0;
  //! Position of V
  static const constexpr ZFSId   V = 1;
  //! Position of W
  static const constexpr ZFSId   W = 2;
  //! Pointer for the velocities so you can use them in a loop
  ZFSId * VV;
  //! Position of P
  static const constexpr ZFSId   P = nDim + 1;
  //! Position of RHO (equal to P in this case)
  static const constexpr ZFSId   RHO = nDim;
  //! Position of T
  static const constexpr ZFSId  T = nDim + 1;
  //! Position of Z
  static const constexpr ZFSId  Z = nDim + 2;
  //! Position of C
  static const constexpr ZFSId  C = nDim + 2;
  //! Position of Yi
  ZFSId* Y;
  //ZFSId  G;
  //! first Position of RANS Variables
  ZFSId   RANS_FIRST;
  //! Position of RANS Variables
  ZFSId*  RANS_VAR;

  //! The Nr. of primitive variables (nDim + 2)
  ZFSId   noVariables;
  ZFSId   noGCellVariables;
  ZFSId   m_noSpecies;
  ZFSId   m_noRansEquations;
  ZFSFloat UInfinity, VInfinity, WInfinity, PInfinity, TInfinity, TInfinityGhost;
  ZFSFloat * VVInfinity;
  ZFSFloat * ransInfinity;
  ZFSFloat DthInfinity, muInfinity, DInfinity;
};

////////////////////////////////////////////////////////////////////////////////
/// Classes with constant number of space dimensions nd

namespace zfs {

template<bool, ZFSId trueValue, ZFSId falseValue> struct conditional_value { static const ZFSId value = falseValue; };
template<ZFSId trueValue, ZFSId falseValue> struct conditional_value<true,trueValue,falseValue> { static const ZFSId value = trueValue;};

namespace fv {

namespace variables {
/// Error value: produces a segmentation fault if used as index.
static const ZFSId Segfault = -1000000000;
}

/// \brief Static indices for accessing conservative variables
/// in nd spatial dimensions
template<ZFSId nd> struct ConservativeVariables {
  static constexpr ZFSId RHO_U = 0;
  static constexpr ZFSId RHO_V = 1;
  static constexpr ZFSId RHO_W = conditional_value<nd==3,2,variables::Segfault>::value;
  //constexpr for arrays is not supported by some compilers (Nov. 2015)
  static const ZFSId RHO_VV[3];
  static constexpr ZFSId RHO_E = nd;
  static constexpr ZFSId RHO = nd + 1;
  static constexpr ZFSId RHO_Z = nd + 2;
  static constexpr ZFSId RHO_C = nd + 2;

  const ZFSId m_noSpecies;
  const ZFSId noVariables;
  ZFSId*  RHO_Y;

  ConservativeVariables(const ZFSId noSpecies);
  ~ConservativeVariables();
};

/// \brief Static indices for accessing primitive variables
/// in nd spatial dimensions
template<ZFSId nd> struct PrimitiveVariables {
  static constexpr ZFSId U = 0;
  static constexpr ZFSId V = 1;
  static constexpr ZFSId W = conditional_value<nd==3,2,variables::Segfault>::value;
  //constexpr for arrays is not supported by some compilers (Nov. 2015)
  static const ZFSId VV[3];
  static constexpr ZFSId RHO = nd;
  static constexpr ZFSId P = nd + 1;
  static constexpr ZFSId T = nd + 1;
  static constexpr ZFSId Z = nd + 2;
  static constexpr ZFSId C = nd + 2;

  static constexpr ZFSId noGCellVariables = nd + 1;
  const ZFSId m_noSpecies;
  const ZFSId noVariables;

  ZFSId* Y;
  PrimitiveVariables(const ZFSId noSpecies);
  ~PrimitiveVariables();
};

} // namespace fv
} // namespace zfs
#endif
