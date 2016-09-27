/* start of smf.c */

/* Use #ifndef and #define to avoid double inclusion of the header file.
   If the file has been previously included, than the preprocessor will 
   skip everything up to the #endif appearing at the end of this file */
#ifndef smf_c 
#define smf_c 

#define AT_SCF

/* SMF module */
#include "smfio.h"

/* Global PGS seed file */
#include "PGS_MODIS_37121.h"


void modsmf(PGSt_SMF_code mnemonicstring,
	    char *infostring,
	    char *functionstring)
/****************************************************************************
!C

!Description:

  This module writes errors and messages through the SDP toolkit to
  LogStatus and causes program termination upon encountering a fatal
  error (as determined by the error mnemonic).  As input it acccepts an
  error mnemonic, string containing user information, and a string
  containing the function name and module where the error/message
  occurs.

!Input Parameters:

  PGSt_SMF_code mnemonicstring      SMF message mnemonic
  char *infostring                  String containing info to be written to 
                                    LogStatus
  char *functionstring              String containing function/module name
                                    to be written to LogStatus

!Output Parameters:
  None
  
!Revision History:
$Log: smf.c,v $
   Revision 1.2  1998/02/12  10:26     ktai
   Added AT_SCF functionality to notify user if toolkit fails
 
   Revision 1.1  1995/07/12  19:27:46  fisher
   Initial revision
  
$Id: smf.c,v 1.1 1995/07/12 19:27:46 fisher Exp $

!Team-unique Header:
  This software is modified by the MODIS Science Data Support
  Team for the National Aeronautics and Space Administration,
  Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
  Written by Paul S. Fisher

    Research and Data Systems Corporation
    SAIC/GSC MODIS Support Office
    7501 Forbes Blvd
    Seabrook MD 20706  

    301-352-2129

    fisher@modis-xl.gsfc.nasa.gov

!Design Notes:

  The determination of time is through native C function(s).

  Upon determination of a fatal error, exit code 1 is used.
  Exit code 1 is also used to indicate SDP toolkit error.  These exit functions
  may be updated at a later date if required by the processing string.

  If AT_SCF is defined, output to screen upon smf message failure.

  Designed to run with seed file format of:
    MODIS_X_MNEMONIC STRING %sStatic Message.  %s

  Example of function call:

  modsmf(MODIS_X_MNEMONIC_STRING, "user message string", "function, 
                                                     module.c string");

  -- The mnemonic string must be identical to one contained in an included 
       seed file.
  -- The user message string can contain anything that needs to be written
       to the log(s).
  -- The function, module.c string should contain "the name of the function,
       the name of the module".

					      
!Externals:
  constants:
       PGS_SMF_MAX_MSGBUF_SIZE                 (PGS_SMF.h)
       PGS_TRUE                                (PGS_SMF.h)
       PGS_FALSE                               (PGS_SMF.h)
       PGS_E_UNIX                              (PGS_SMF.h)
       PGS_E_UNDEFINED_CODE                    (PGS_SMF.h)
       PGS_SMF_E_LOGFILE                       (PGS_SMF.h)
       PGS_S_SUCCESS                           (PGS_SMF.h)
       NULL                                    <stdio.h>

  typedefs:
       PGSt_SMF_status                         (PGS_SMF.h)
       struct tm                               (smfio.h)<time.h>
       time_t                                  (smfio.h)<time.h>

  functions:
       PGS_SMF_TestFatalLevel                  (PGS_SMF_h)
       PGS_SMF_GetMsgByCode                    (PGS_SMF.h)
       PGS_SMF_SetDynamicMsg                   (PGS_SMF.h)
       exit                                    <stlib.h>
       sprintf                                 <stdio.h>
       time                                    (smfio.h)<time.h>
       localtime                               (smfio.h)<time.h>
       asctime                                 (smfio.h)<time.h>


!END
*****************************************************************************/


{
  
  struct tm *local=NULL;                         /* Time functions */
  time_t t;

  PGSt_SMF_status messageString=PGS_S_SUCCESS;  /* Holds string of message 
						   from 
						   PGS_SMF_GetMsgByCode */
	
  char message[PGS_SMF_MAX_MSGBUF_SIZE]={0};    /* holds the message string 
					       associated with
					       the error code returned 
					       by GetMsg */

  char buf[PGS_SMF_MAX_MSGBUF_SIZE]={0};        /* Holds a concatenated dynamic
					       message string */


  t=time(NULL);                             /* Calculate timestamp */
  local=(struct tm *)localtime(&t);
  
  //DEBUG: Josh
  printf("SMF 0\n");

  /* Get message from seed file based on error mnemonic */
  messageString=PGS_SMF_GetMsgByCode(mnemonicstring, message);
  if(messageString!=PGS_S_SUCCESS) 
  /* Modified for not allowing exit(99); in the code  ktai, 02/12/98 */
  {
#ifdef AT_SCF
#endif
//DEBUG: Josh
printf("SMF 1\n");
printf("UHHHHHH\n");
printf("Mnemonic String: %s \n", mnemonicstring);
printf("Uhhhh\n");
   exit(1);
  }

  /* Concatenate timestamp, and user massage w/ buffer message */
  sprintf(buf, message, asctime(local), infostring);
  
  /* Write message to LogStatus */
  messageString=PGS_SMF_SetDynamicMsg(mnemonicstring, buf, functionstring);
  //DEBUG: Josh
  printf("SMF 2\n");
  
  if(messageString!=PGS_S_SUCCESS) 
  /* Modified for not allowing exit(99); in the code   ktai, 02/12/98 */
   {
#ifdef AT_SFC
#endif
//DEBUG: Josh
printf("SMF 3\n");
   exit(1);
  }

  /* Test for fatal error */
  if(PGS_SMF_TestFatalLevel(mnemonicstring)==PGS_TRUE)
  {
  	//DEBUG: Josh
	printf("SMF 4\n");
  	exit(1);
  }  
  else
  {
  //DEBUG: Josh
  printf("SMF 5\n");
  return;

}

#endif  /* part of #ifndef/#define statements at beginning of file */
}
/* End of smf.c */
