/* start of smfio.h */

/* Use #ifndef and #define to avoid double inclusion of the header file.
   If the file has been previously included, than the preprocessor will 
   skip everything up to the #endif appearing at the end of this file */
#ifndef smfio_h
#define smfio_h

/*****************************************************************************
!C-INC

!Description:

  Header file necessary for use of io.c and smf.c SDP toolkit
  integration modules.

!Input Parameters: None

!Output Parameters: None
  
!Revision History:
 $Log: smfio.h,v $
 * Revision 1.3  1996/03/08  20:21:50  fisher
 * Added memory definition for timea.c
 *
 * Revision 1.2  1995/07/14  17:57:12  fisher
 * Prologue update
 *
 * Revision 1.1  1995/07/12  18:58:26  fisher
 * Initial revision
 *
 $Id: smfio.h,v 1.3 1996/03/08 20:21:50 fisher Exp $

!Team-unique Header:
  This software has been created by the MODIS Science Data Support
  Team for the National Aeronautics and Space Administration,
  Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits

    Written by Paul S. Fisher
    Research and Data Systems Corporation
    SAIC/GSC MODIS Support Office
    7501 Forbes Blvd
    Seabrook MD 20706  

    fisher@modis-xl.gsfc.nasa.gov

!Design Notes

    This header file should be included in all relevant, integrated
    pieces of source code.  Both io.c and smf.c should be included on
    the compilation command line, or in the Makefile.

!Externals: 

    PGSd_IO_Gen_Read               (PGS_IO.h)
    PGSd_IO_Gen_Write              (PGS_IO.h)
    PGSd_IO_Gen_Append             (PGS_IO.h)
    PGSd_IO_Gen_Update             (PGS_IO.h)
    PGSd_IO_Gen_Trunc              (PGS_IO.h)
    PGSd_IO_Gen_AppendUpdate       (PGS_IO.h)
    PGSt_SMF_Status                (PGS_SMF.h)


!END
*****************************************************************************/



/* SDP Toolkit I/O header file */
#include "PGS_IO.h"

/* SDP Toolkit Status Message Facility header file */
#include "PGS_SMF.h"

/* C time functions */
#include <time.h>

#define SMF_MAX_ACT_SIZE 80              /* this is supposed to be defined
                                            as PGSd_SMF_MAX_ACT_SIZE,
                                            according to the SDP Toolkit
                                            Primer, in the SDP Toolkit;
                                            until it is, we'll define
                                            it here */

#ifndef READ
#define READ   PGSd_IO_Gen_Read            /* read "r"*/
#endif

#ifndef WRITE
#define WRITE  PGSd_IO_Gen_Write           /* write, truncate, create "w"*/
#endif

#ifndef APPEND
#define APPEND PGSd_IO_Gen_Append          /* write/append, create "a"*/
#endif

#ifndef UPDATE
#define UPDATE PGSd_IO_Gen_Update          /* read/write "r+"*/
#endif

#ifndef TRUNC  
#define TRUNC  PGSd_IO_Gen_Trunc           /* read/write, truncate, create 
					      "w+" */
#endif

#ifndef APUP   
#define APUP   PGSd_IO_Gen_AppendUpdate    /* read whole file, write/append,
					      create "a+"*/
#endif

/* Prototypes */
void modsmf(PGSt_SMF_code mnemonicstring,
	    char *infostring,
	    char *functionstring);

PGSt_SMF_status modfileopen(PGSt_PC_Logical FilePCFIndex,
			    PGSt_IO_Gen_AccessType FileAccessFlag,
			    PGSt_IO_Gen_FileHandle **FilePtr,
			    PGSt_integer VersionNumber);

PGSt_SMF_status modfileclose(PGSt_PC_Logical FilePCFIndex,
			     PGSt_IO_Gen_FileHandle *FilePtr);

PGSt_SMF_status current_time_a(char *);

/* definition for proper memory allocation of current_time_a() */

#define TIMECODEASIZE 28


#endif  /* part of #ifndef/#define statements at beginning of file */

/* End of smfio.h */
