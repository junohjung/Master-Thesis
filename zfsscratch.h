#ifndef ZFSSCRATCH_H
#define ZFSSCRATCH_H

#include <iomanip>
#include <list>
#include <ostream>
#include <sstream>
#include <cstdint>
#include <new>
#include "zfsconfig.h"
#include "zfsfunctions.h"
#include "zfsglobalvariables.h"
#include "zfsinfoout.h"
#include "zfstypes.h"

//#if defined(linux)
//#include <sys/sysinfo.h>
//const uintptr_t ALIGNMENT_BOUNDARY = sysconf(_SC_LEVEL1_DCACHE_LINESIZE); // make alignment boundary equal to L1-cache linesize
//#else
//const uintptr_t ALIGNMENT_BOUNDARY = 64; // see below
//#endif


/// \brief type traits for dealing with signed vs unsigned indices
///
/// Note: we should not do this...
/// \todo use <type_traits> instead of this hack
namespace zfs {
struct zfs_signed {};
struct zfs_unsigned {};
template<class T> struct is_unsigned {
  static const bool value = false;
  typedef zfs_signed type;
};
template<> struct is_unsigned<ZFSUint> {
  static const bool value = true;
  typedef zfs_unsigned type;
};
// used in cartesian grid constructor, comment out to see the errors:
template<> struct is_unsigned<long unsigned int> {
  static const bool value = true;
  typedef zfs_unsigned type;
};
} // namespace zfs

template < class T >
class ZFSScratchSpace;
class ZFSScratchSpaceBase;

typedef ZFSScratchSpace<ZFSString> ZFSStringScratchSpace;
typedef ZFSScratchSpace<ZFSChar> ZFSCharScratchSpace;
typedef ZFSScratchSpace<ZFSInt> ZFSIntScratchSpace;
typedef ZFSScratchSpace<ZFSId> ZFSIdScratchSpace;
typedef ZFSScratchSpace<ZFSLong> ZFSLongScratchSpace;
typedef ZFSScratchSpace<ZFSFloat> ZFSFloatScratchSpace;
typedef ZFSScratchSpace<ZFSBool> ZFSBoolScratchSpace;
typedef ZFSScratchSpace<ZFSUshort> ZFSUshortScratchSpace;

typedef ZFSScratchSpace<ZFSInt *> ZFSIntPointerScratchSpace;
typedef ZFSScratchSpace<ZFSId *> ZFSIdPointerScratchSpace;
typedef ZFSScratchSpace<ZFSLong *> ZFSLongPointerScratchSpace;
typedef ZFSScratchSpace<ZFSFloat *> ZFSFloatPointerScratchSpace;
typedef ZFSScratchSpace<ZFSBool *> ZFSBoolPointerScratchSpace;

typedef ZFSScratchSpace<ZFSFloat **> ZFSFloatPointerPointerScratchSpace;

typedef ZFSScratchSpace<MPI_Group> ZFSMPI_GroupScratchSpace;

typedef std::list<ZFSScratchSpaceBase *> ZFSScratchList;

/** \brief This class holds the complete scratch space
 *
 * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
 * \date 10.05.2011, 10.06.2011, 21.10.2011
 *
 * This class defines the overall size of the scratch space and tracks
 * the sizes and allocation of newly requested Scratch space arrays.
 *
 **/
class ZFSScratch
{

  friend class ZFSScratchSpaceBase;

  static const uintptr_t ALIGNMENT_BOUNDARY = ZFS_SCRATCH_ALIGNMENT_BOUNDARY; // see above

  public:

  /** \brief Constructor
   *
   * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
   * \date 10.05.2011, 10.06.2011, 21.10.2011
   * \note Base memory now points to an aligned memory address, Lennart Schneiders, 12.12.2012
   *
   * Allocates space of "size * sizeof(ZFSFloat) * Cells" bytes, that is "size" ZFSFloat for every cell.
   *
   * \param[in] size number of ZFSFloat for every cell.
   * \param[in] Cells number of cells.
   *
   **/
  ZFSScratch(ZFSFloat size, ZFSInt Cells);

  /** \brief Destructor
   *
   * \author Andreas Lintermann, Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   *
   **/
  ~ZFSScratch()
  {
    m_number_of_cells = 0;
    m_number_of_elements = 0;
    m_object_id = 0;
    m_usedmemsize = 0;
    m_nextfree = nullptr;
    m_scratchSpaces.clear();
    m_maxused = nullptr;
    m_report = "n/a";
    delete [] m_totalScratch;
    m_totalScratch = nullptr;
  };

  /** \brief Returns a pointer to the end of the scratch.
   *
   * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
   * \date 10.05.2011, 10.06.2011, 21.01.2011
   *
   * \return a pointer to the end of the scratch
   **/
  static char* getEndPointer() {return (char*)(size_t)(m_totalScratch + m_number_of_elements - 1);}

  /** \brief Returns the amount of available memory in scratch.
   *
   * \author Andreas Lintermann, Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   *
   * \return the amount of available memory
   **/
  //static size_t getAvailableMemory() {return (getTotalMemory() - m_usedmemsize);}//this leads to errors if the m_nextfree pointer is corrupted for whatever reason, Lennart
  static size_t getAvailableMemory() {return (&m_totalScratch[m_number_of_elements-1] - m_nextfree);}

  /** \brief Returns the amount of total available memory in scratch.
   *
   * \author Andreas Lintermann. Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   *
   * \return the amount of total available memory
   **/
  static size_t getTotalMemory() {return (getEndPointer() + 1 - m_totalScratch);}

  static ZFSString printSelfScratch();
  static ZFSString printSelf();
  static ZFSString printSelfReport();

  static size_t m_number_of_cells;
  static size_t m_number_of_elements;
  static ZFSInt m_object_id;
  static char* m_totalScratch;
  static size_t m_usedmemsize;
  static char* m_nextfree;
  static ZFSScratchList m_scratchSpaces;
  static char* m_maxused;
  static ZFSString m_report;

};



/** \brief This class is a base class for all ZFSScratchSpaces.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * This class is the base of all derived scratch space elements.
 **/
class ZFSScratchSpaceBase {

  public:

  static const uintptr_t ALIGNMENT_BOUNDARY = ZFS_SCRATCH_ALIGNMENT_BOUNDARY; // see above

  const size_t m_memsize;
  const size_t m_memsizePadded;
  const ZFSString m_calling_function;
  const ZFSString m_variable_name;
  ZFSInt m_object_id;
  bool m_nonterminal;
  bool m_destroy;

  /** \brief Constructor
   *
   * \author Andreas Lintermann, Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   * \note Each ScratchSpace now points to an aligned memory address, Lennart Schneiders, 12.12.2012
   *
   * Sets the scratch space elements properties.
   *
   * \param[num] num number of elements in array
   * \param[size] size the size of the array in bytes
   * \param[name] name the name of the calling function
   * \param[varname] varname the name of the array variable
   *
   **/
 ZFSScratchSpaceBase(size_t num, size_t size, const ZFSString& name, const ZFSString& varname) :
  m_memsize(num*size), m_memsizePadded( (num*size)+((ALIGNMENT_BOUNDARY-((num*size)%ALIGNMENT_BOUNDARY))%ALIGNMENT_BOUNDARY) ),
    m_calling_function(name), m_variable_name(varname), m_nonterminal(false), m_destroy(false)
  {
    ASSERT( ( ((uintptr_t)ZFSScratch::m_nextfree) % ALIGNMENT_BOUNDARY == 0 ), "Scratch memory is not aligned" );
  }
  virtual ZFSString printSelfReport() const = 0;
  virtual ZFSString printSelf() const = 0;
};

/** \brief This class is a ZFSScratchSpace.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * This class is a derivative of ZFSScratchSpaceBase and is responsible
 * for the management of the scratch space to be allocated. Scratch objects
 * partially model the container concept and should be usable with some,
 * but not all, STL algorithms.
 **/
template <class T>
class ZFSScratchSpace : public ZFSScratchSpaceBase
{

 public:

  typedef T               value_type;
  typedef value_type&       reference;
  typedef const value_type& const_reference;
  typedef value_type*       pointer;
  typedef const value_type* const_pointer;
  typedef pointer           iterator;
  typedef const_pointer     const_iterator;
  typedef long int        difference_type;
  typedef std::size_t     size_type;

  ZFSScratchSpace(ZFSInt num, const ZFSString& name, const ZFSString& varname); //!< 1D Scratch constructor
  ZFSScratchSpace(ZFSInt num0, ZFSInt num1, const ZFSString& name, const ZFSString& varname); //!< 2D Scratch constructor
  ZFSScratchSpace(ZFSInt num0, ZFSInt num1, ZFSInt num2, const ZFSString& name, const ZFSString& varname); //!< 3D Scratch constructor
  ~ZFSScratchSpace();

  template<class Integral> inline reference operator[](const Integral); //!< 1D Scratch access
  template<class Integral> inline reference operator()(const Integral); //!< 1D Scratch access
  template<class Integral> inline reference operator()(const Integral,const Integral); //!< 2D Scratch access
  template<class Integral> inline reference operator()(const Integral,const Integral,const Integral); //!< 3D Scratch access
  template<class Integral> inline const_reference operator[](const Integral i) const; //!< 1D Scratch const access
  template<class Integral> inline const_reference operator()(const Integral i) const; //!< 1D Scratch const access
  template<class Integral> inline const_reference operator()(const Integral,const Integral) const; //!< 2D Scratch const access
  template<class Integral> inline const_reference operator()(const Integral,const Integral,const Integral) const; //!< 3D Scratch const access
  ZFSScratchSpace<T>& operator=(ZFSScratchSpace<T>&); //!< copy entries of S to this, given there is enough space available
  template <class>
  friend std::ostream& operator<<(std::ostream& os, const ZFSScratchSpace& s);

  inline iterator  begin()        { checkForEmptyScratch(); return p; }
  inline iterator  cbegin() const { checkForEmptyScratch(); return p; }
  inline iterator  end()          { checkForEmptyScratch(); return last; }
  inline iterator  cend()   const { checkForEmptyScratch(); return last; }
  inline ZFSInt size0() const { return m_size0; }
  inline ZFSInt size1() const { return m_size1; }
  inline ZFSInt size2() const { return m_size2; }
  inline size_t getMemsize() const { return (last-p)*sizeof(T); }
  inline size_type size() const {
    ASSERT(last >= p, "Scratch cannot have negative size!");
    return last - p;
  }
  inline bool empty() const {
    ASSERT(last >= p, "Scratch cannot have negative size!");
    return (last - p == 0) ? true : false;
  }

  void fill( T val ) { std::fill( p, last, val ); } //!< fill the scratch with a given value
  ZFSString printSelf() const ;
  ZFSString printSelfReport() const;

  pointer p; ///< Deprecated: use [] instead!
  T* getPointer() const {return p;} ///< Deprecated: use begin() instead!

 private:
  // TODO: when ".p" is removed: const_iterator first, last;
  pointer last;
  // TODO: this should be std::size_t
  const ZFSInt m_size0;
  const ZFSInt m_size1;
  const ZFSInt m_size2;

  // Disable copying and heap allocation
  ZFSScratchSpace() = delete;
  ZFSScratchSpace(const ZFSScratchSpace<T>&) = delete;

  void* operator new(std::size_t)           = delete;
  void* operator new(std::size_t,void* p)   = delete;
  void* operator new[](std::size_t)         = delete;
  void* operator new[](std::size_t,void* p) = delete;
  void  operator delete(void*)              = delete;
  void  operator delete(void* p,void*)      = delete;
  void  operator delete[](void*)            = delete;
  void  operator delete[](void* p,void*)    = delete;

  /// \brief Assert the bounds of 1D acces to an array
  ///
  /// Note: the condition is tested by checkBounds,
  /// which is overloaded for signed and unsigned indices
  /// using the zfs::is_unsigned trait.
  ///
  /// \todo replace zfs::is_unsigned with <type_traits>
  ///
  /// @{
  template<class Integral> inline void testBounds(const Integral& i) const;
  template <class Integral> inline bool checkBounds(const Integral& i) const {
    typedef typename zfs::is_unsigned<Integral>::type type;
    return checkBounds(i,type());
  }
  template <class Integral>
  inline bool checkBounds(const Integral& i, zfs::zfs_unsigned) const {
    return static_cast<size_type>(i) < size();
  }
  template <class Integral>
  inline bool checkBounds(const Integral& i, zfs::zfs_signed) const {
    return i >= 0 && static_cast<size_type>(i) < size();
  }
  /// @}
  void init();
  inline void checkForEmptyScratch() const;
};


/** \brief Constructor
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * Allocates a scratch space element of a certain type.
 * If not enough memory is available provided by the ZFSScratch class, it terminates the program run.
 * Otherwise the properties of this management object are set and a pointer to this object is stored in the management list of the ZFSScratch class.
 * If a maximum of memory of the ZFSScratch is used, the state stored to be outputed at the end of the program run to enable memory analysis.
 *
 * \param[num] num number of elements in array
 * \param[name] name the name of the calling function
 * \param[varname] varname the name of the array variable
 *
 * Note: If ZFS_EXTRA_DEBUG is defined writes a SCRATCH_DOUBLE_ALLOCATION WARNING
 * at run-time if allocation of two variables with the same name is performed.
 **/
template<class T> ZFSScratchSpace<T>::ZFSScratchSpace
(ZFSInt num, const ZFSString& name, const ZFSString& varname) : ZFSScratchSpaceBase(zfsMAX(1,num), sizeof(T), name, varname),
  p(reinterpret_cast<T*>(ZFSScratch::m_nextfree)), last(p+num), m_size0(num), m_size1(1), m_size2(1) {
  this->init();
}

/** \brief Constructor for 2D scratch
 *
 * \author Lennart
 * \date 26.02.2013
 *
 * This is the 2D-extension of the 1D constructor above. Unfortunately, the order of constructor arguments
 * prohibits the use of a constructor default value for the second dimension, so most of the code here is redundant.
 *
 **/
template<class T> ZFSScratchSpace<T>::ZFSScratchSpace
(ZFSInt num0, ZFSInt num1, const ZFSString& name, const ZFSString& varname) : ZFSScratchSpaceBase(zfsMAX(1,num0*num1), sizeof(T), name, varname),
  p(reinterpret_cast<T*>(ZFSScratch::m_nextfree)), last(p+(num0*num1)), m_size0(num0), m_size1(num1), m_size2(1) {
  this->init();
}


/** \brief Constructor for 3D scratch
 *
 * \author Lennart
 * \date 26.02.2013
 *
 * This is the 3D-extension of the 1D constructor above. Unfortunately, the order of constructor arguments
 * prohibits the use of constructor default values for the second and third dimension, so most of the code here is redundant.
 *
 **/
template<class T> ZFSScratchSpace<T>::ZFSScratchSpace
(ZFSInt num0, ZFSInt num1, ZFSInt num2, const ZFSString& name, const ZFSString& varname) : ZFSScratchSpaceBase(zfsMAX(1,num0*num1*num2), sizeof(T), name, varname),
  p(reinterpret_cast<T*>(ZFSScratch::m_nextfree)), last(p+(num0*num1*num2)), m_size0(num0), m_size1(num1), m_size2(num2) {
  this->init();
}


/** \brief extended constructor functionality
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 * \brief Content was moved from the constructor body to this function since it is the same for all above constructors (Lennart)
 * \note added 'placement new' operator, as previously the use of non-primitive types resulted in undefined behavior, Lennart 2013
 *
 **/
template<class T> void ZFSScratchSpace<T>::init() {

  ZFSScratch::m_object_id++;

  m_object_id=ZFSScratch::m_object_id;
  ZFSScratch::m_scratchSpaces.push_back(this);

  /// Check if allocation would exceed the available scratch space:
  if(ZFSScratch::getAvailableMemory() < m_memsizePadded) {
    m_destroy=true;
    std::cerr << ZFSScratch::printSelf();
    zfs_log << ZFSScratch::printSelf();
    std::stringstream errorMessage, suggestionMessage;
    errorMessage << "ZFSScratch is not big enough - exceeded size in "
                 << m_calling_function << " with variable " << m_variable_name << std::endl;
    errorMessage << "           " << m_memsize/(1024.0*1024.0)
                 << " MB required, but there are only "
                 << ZFSScratch::getAvailableMemory()/(1024.0*1024.0)
                 << " MB available" << std::endl;
    errorMessage << "           Suggested to raise the property scratchSize from " <<
        ZFSFloat(ZFSScratch::m_number_of_elements)
        /(sizeof(ZFSFloat)*ZFSScratch::m_number_of_cells)
                 << " to " <<
        ZFSFloat(ZFSScratch::m_number_of_elements + m_memsize - ZFSScratch::getAvailableMemory())
        /(sizeof(ZFSFloat)*ZFSScratch::m_number_of_cells)
                 << std::endl;
    zfs_log << errorMessage.str();
    zfsTerm(1,__CALLING_FUNCTION__,errorMessage.str());
  }
  ZFSScratch::m_usedmemsize += m_memsizePadded;
  ZFSScratch::m_nextfree += m_memsizePadded;

  for ( ZFSInt i = 0; i < (signed)this->size(); i++ ) {
    ASSERT( ((char*)&p[i] >= &ZFSScratch::m_totalScratch[0]) && ((char*)&p[i] <= &ZFSScratch::m_totalScratch[ZFSScratch::m_number_of_elements-1]), "ScratchSpace '" + m_variable_name + "' out of bounds during init." );
    new (&p[i]) T(); // 'placement new' operator calling each element's default constructor (no further memory allocated)
  }

  for ( ZFSInt i = 0; i < (signed)this->size(); i++ ) {
    new (&p[i]) T(); // 'placement new' operator calling each element's default constructor
  }

  /// Check for double allocation:
  #ifdef ZFS_EXTRA_DEBUG
  for(ZFSScratchList::iterator it = ZFSScratch::m_scratchSpaces.begin();
      it != ZFSScratch::m_scratchSpaces.end(); ++it) {
    if((*it) == this) continue;
    if ((*it)->m_variable_name == this->m_variable_name) {
      zfs_log <<  "WARNING: SCRATCH_DOUBLE_ALLOCATION: double usage of variable "
           << this->m_variable_name << " allocated in "
           << (*it)->m_calling_function << " and "
           << this->m_calling_function << std::endl;
    }
  }
  #endif

  if(ZFSScratch::m_nextfree > ZFSScratch::m_maxused) {
    ZFSScratch::m_maxused = ZFSScratch::m_nextfree;
    ZFSScratch::m_report = ZFSScratch::printSelfScratch();
    ZFSScratchList::iterator iter;
    for(iter = ZFSScratch::m_scratchSpaces.begin(); iter!=ZFSScratch::m_scratchSpaces.end();
        iter++) {
      ZFSScratch::m_report += (*iter)->printSelfReport();
      ZFSScratch::m_report += "\n";
    }
    ZFSScratch::m_report += "\n\n";
  }
}



/** \brief Destructor
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * Removes this scratch space elemet from the list of all scratch elements.
 **/
template<class T> ZFSScratchSpace<T>::~ZFSScratchSpace() {
  if(!m_destroy) {
      ZFSScratchList::iterator it = ZFSScratch::m_scratchSpaces.end();
      it--;

      //! Only delete if we are at the last element (guaranteed by C++)
      if((*it)->m_object_id == m_object_id) {
        ZFSScratch::m_nextfree -= m_memsizePadded;
        ZFSScratch::m_usedmemsize -= m_memsizePadded;
        ZFSScratch::m_scratchSpaces.erase(it);
      }

      if(m_nonterminal)
        ZFSScratch::m_nextfree -= m_memsizePadded;
  }
}

/** \brief Returns a string summing up this scratch space element information.
 *
 * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
 * \date 10.05.2011, 10.06.2011, 21.10.2011
 *
 * \return a string summing up this scratch space element information
 **/
template <class T> ZFSString ZFSScratchSpace<T>::printSelf() const {
  std::stringstream id,mem,memsize;
  mem << p;

  id <<  m_object_id;
  memsize << m_memsize;

  ZFSString message = "ZFSScratchSpace ID: ";
  message += id.str();
  message += "\n-----------------------------\n";
  message += "Calling function:\t";
  message += m_calling_function;
  message += "\nVariable name:\t\t";
  message += m_variable_name;
  message += "\nMemory size:\t\t";
  message += memsize.str();
  message += "\nMemory pointer:\t\t";
  message += mem.str();
  message += "\nSize:\t\t\t";
  message += std::to_string(size());
  message += "\n";

  return message;
}

/** \brief Returns a shortened string summing up this scratch space element information.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * This function is used in the report process for tracking the occasion of maximal
 * memory usage during program execution.
 *
 * \return a shortened string summing up this scratch space element information
 **/
template <class T> ZFSString ZFSScratchSpace<T>::printSelfReport() const {
  std::stringstream id,memsize,funct,var;
  memsize << m_memsize;
  id << m_object_id;

  ZFSString message = "ID: ";
  message += id.str();
  message += "    SIZE: ";
  message += memsize.str();
  message += "    FUNCT: ";
  message += m_calling_function;
  message += "    VAR: ";
  message += m_variable_name;
  return message;
}

template<class T> template<class Integral>
inline void ZFSScratchSpace<T>::testBounds(const Integral& i) const {
  ASSERT(checkBounds(i), "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
      << m_variable_name << " allocated from " << m_calling_function
      << " with index " << i << " out of range [0," << size() << ")." << std::endl
      << "Note: if the array range = [0,0) the array is EMPTY. "
      << "You cannot dereference elements of an empty array!");
}

/// \brief Access the i-th element of the scratch object.
///
/// \return Refernence to the i-th scratch space element
///
/// \param[i] Index of Integral type.
/// \param[Integral] Type that models the Integral concept.
///
/// Note: if NDEBUG is not defined, the array bounds are checked:
/// i >= 0 and i < size.
/// If bound checking fails the program will terminate with a
/// SCRATCH_OUT_OF_BOUNDS ERROR.
///
template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::reference ZFSScratchSpace<T>::operator()(const Integral i) {
#ifndef NDEBUG
  testBounds(i);
#endif
  return p[i];
}
template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::const_reference ZFSScratchSpace<T>::operator()(const Integral i) const {
#ifndef NDEBUG
  testBounds(i);
#endif
  return p[i];
}

/// For 1D access the [] operator is also provided, the () operator should be preferred.
template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::reference ZFSScratchSpace<T>::operator[](const Integral i) {
  return operator()(i);
}
template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::const_reference ZFSScratchSpace<T>::operator[](const Integral i) const {
  return operator()(i);
}

/// \brief 2D and 3D Scratch access
///
/// Note: this would probably result in a compilation error
/// if m_sizeX
template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::reference ZFSScratchSpace<T>::operator()(const Integral i, const Integral j) {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1, "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
      << m_variable_name << " allocated from " << m_calling_function
         << " with indices " << i << "," << j << ". Dimensions are " << m_size0 << "," << m_size1 << "." );
  return p[ m_size1*i + j ];
}

template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::reference ZFSScratchSpace<T>::operator()(const Integral i, const Integral j, const Integral k) {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1 && k >= 0 && k < m_size2, "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
      << m_variable_name << " allocated from " << m_calling_function
         << " with indices " << i << "," << j << "," << k << ". Dimensions are " << m_size0 << "," << m_size1 << "," << m_size2 << "." );
  return p[ m_size1*m_size2*i + m_size2*j + k ];
}

template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::const_reference ZFSScratchSpace<T>::operator()(const Integral i, const Integral j) const {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1, "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
      << m_variable_name << " allocated from " << m_calling_function
         << " with indices " << i << "," << j << ". Dimensions are " << m_size0 << "," << m_size1 << "." );
  return p[ m_size1*i + j ];
}

template<class T> template<class Integral>
inline typename ZFSScratchSpace<T>::const_reference ZFSScratchSpace<T>::operator()(const Integral i, const Integral j, const Integral k) const {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1 && k >= 0 && k < m_size2, "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
      << m_variable_name << " allocated from " << m_calling_function
         << " with indices " << i << "," << j << "," << k << ". Dimensions are " << m_size0 << "," << m_size1 << "," << m_size2 << "." );
  return p[ m_size1*m_size2*i + m_size2*j + k ];
}

template<class T>
ZFSScratchSpace<T>& ZFSScratchSpace<T>::operator=(ZFSScratchSpace<T>& S) {
  for ( ZFSInt i = 0; i < zfsMIN(m_size0,S.m_size0); i++ ) {
    for ( ZFSInt j = 0; j < zfsMIN(m_size1,S.m_size1); j++ ) {
      for ( ZFSInt k = 0; k < zfsMIN(m_size2,S.m_size2); k++ ) {
        p[ m_size1*m_size2*i + m_size2*j + k ] = S.p[ m_size1*m_size2*i + m_size2*j + k ];
      }
    }
  }
  return *this;
}


/// \brief Print contents of scratch space object.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-13
///
/// \tparam T Element type.
/// \param[in] os Stream object to write to.
/// \param[in] s Scratch space object.
///
/// \return Stream object to allow chained writing
///
/// Note: if T is not convertible to string, the contents are not printed.
template <class T>
std::ostream& operator<<(std::ostream& os, const ZFSScratchSpace<T>& s) {
  os << s.printSelf();
  os << "Content:\n";
  for (typename ZFSScratchSpace<T>::size_type i = 0; i < s.size(); i++) {
    os << std::setw(7) << i << ": " << s[i] << "\n";
  }
  return os;
}


/// \brief Checks if scratch space is empty before returning a pointer to it.
///
/// Note: the check is only performed if NDEBUG is not defined,
/// and ZFS_EXTRA_DEBUG is defined.
///
template<class T>
inline void ZFSScratchSpace<T>::checkForEmptyScratch() const {
#ifdef ZFS_EXTRA_DEBUG
  if(empty()) {
    zfs_log << "SCRATCH WARNING: passing a pointer to the empty scratch space "
            << m_variable_name << " allocated from " << m_calling_function
            << "! Make sure this pointer will NEVER be dereferenced!" << std::endl;
  }
#endif
}

#endif // ZFSSCRATCH_H
