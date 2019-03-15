#include "zfsscratch.h"

using namespace std;

size_t ZFSScratch::m_number_of_cells;
size_t ZFSScratch::m_number_of_elements;// = 100;

ZFSScratchList ZFSScratch::m_scratchSpaces;
char* ZFSScratch::m_totalScratch;// = new ZFSInt[ZFSScratch::m_number_of_elements];;
size_t ZFSScratch::m_usedmemsize;
char* ZFSScratch::m_nextfree;
ZFSInt ZFSScratch::m_object_id=0;
char* ZFSScratch::m_maxused;
string ZFSScratch::m_report;

#if defined(ZFS_CLANG_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#endif

ZFSScratch::ZFSScratch(ZFSFloat size, ZFSInt Cells) {
      m_number_of_cells = Cells;
      m_number_of_elements = floor(size*sizeof(ZFSFloat)*m_number_of_cells) + ALIGNMENT_BOUNDARY;
      m_totalScratch = new char[m_number_of_elements];
      m_nextfree = m_totalScratch;
      m_maxused = m_totalScratch;

      m_nextfree += ((ALIGNMENT_BOUNDARY-(((uintptr_t)m_nextfree)%ALIGNMENT_BOUNDARY))%ALIGNMENT_BOUNDARY);
      ASSERT( ( ((uintptr_t)m_nextfree) % ALIGNMENT_BOUNDARY == 0 ), "Scratch memory is not aligned" );

      m_report = printSelfScratch();
}


/** \brief Returns a string summing up the scratch state information.
 *
 * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
 * \date 10.05.2011, 10.06.2011, 21.10.2011
 *
 * \return a string summing up the scratch state information
 **/
string ZFSScratch::printSelfScratch()
{
  stringstream begin,end,nextfree,number,total,totalmb,used,usedmb,available,availablemb,maxmem,maxmemmb;
  const ZFSFloat bytesInAMegaByte = 1048576.0;
  begin << (ZFSInt*)m_totalScratch;
  end << (ZFSInt*)getEndPointer();
  nextfree << (ZFSInt*)m_nextfree;
  number << m_number_of_elements;
  total << getTotalMemory();
  totalmb << (ZFSFloat)getTotalMemory()/bytesInAMegaByte;
  used << m_usedmemsize;
  usedmb << (ZFSFloat)m_usedmemsize/bytesInAMegaByte;
  available << getAvailableMemory();
  availablemb << (ZFSFloat)getAvailableMemory()/bytesInAMegaByte;
  maxmem << (m_maxused-m_totalScratch);
  maxmemmb << (ZFSFloat)(m_maxused-m_totalScratch)/bytesInAMegaByte;
  
  string message = "\n\nScratch:";
  message += "\n-----------------------------\n";
  message += "Bytes allocated:\t\t";
  message += number.str();
  message += "\nTotal memory:\t\t\t";
  message += total.str();
  message += "\t(";
  message += totalmb.str();
  message += "MB)";
  message += "\nUsed memory:\t\t\t";
  message += used.str();
  message += "\t(";
  message += usedmb.str();
  message += "MB)";
  message += "\nFree memory:\t\t\t";
  message += available.str();
  message += "\t(";
  message += availablemb.str();
  message += "MB)";
  message += "\nMax memory:\t\t\t";
  message += maxmem.str();
  message += "\t(";
  message += maxmemmb.str();
  message += "MB)";
  message += "\nScratch start pointer:\t\t";
  message += begin.str();
  message += "\nScratch end pointer:\t\t";
  message += end.str();
  message += "\nNext free pointer:\t\t";
  message += nextfree.str();
  message += "\n\n";
  
  return message;
}

/** \brief Returns a string summing up the scratch state information and all scratch space elements information.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * \return a string summing up the scratch state information and all scratch space elements information
 **/
string ZFSScratch::printSelf()
{
  
  string message = printSelfScratch();
  ZFSScratchList::iterator iter;
  
  for(iter = ZFSScratch::m_scratchSpaces.begin(); iter!=ZFSScratch::m_scratchSpaces.end();iter++)
    {
      message += (*iter)->printSelf();
      message += "\n";
    }
  
  message += "\n\n";
  return message;
}

/** \brief Returns a shortened string summing up the scratch space state information.
 *
 * \author Andreas Lintermann
 * \date 10.05.2011
 *
 * This function is used in the report process for tracking the occasion of maximal 
 * memory usage during program execution.
 *
 * \return a shortened string summing up the scratch space state information
 **/
string ZFSScratch::printSelfReport()
{
  string message = m_report;
  message += "\n\n";
  return message;
}

#if defined(ZFS_CLANG_COMPILER)
#pragma GCC diagnostic pop
#endif
