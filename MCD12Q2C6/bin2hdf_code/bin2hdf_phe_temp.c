/***************************************
This program is used to create temporal HDF-EOS files for LC product.
The import file:
binary Land/water mask.
Previous PHE_LW_TEMP*hdf
The BEGINNINGDATE and ENDINGDATE for example:"2000-07-12","2000-12-18"
Ouput:
updated     PHE_LW_TEMP.HDF (metadata and land/water mask)

****************************************/

/*Header file that contains define 
statements and function prototypes*/

#include "bin2hdf_phe_temp.h"

char IN_FILE[1000];


int main(int argv,char *argc[])   /* changed to return int status*/
  
/******************************************************************************
!C

!Description: 

   This program is used to extract user defined element from HDF-EOS files
   and then save in BIP format.  


!Input Parameters:

    none

!Output Parameters:

    Integer exit status qualifying program termination conditon

!Revision History:


!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:

     Principal Investigator: Alan H. Strahler, Boston University
                alan@crsa.bu.edu

    Developers: Crystal Schaaf, Xiaoyang Zhang, Feng Gao, Jordan Borak
                schaaf@crsa.bu.edu, fgao@crsa.bu.edu
		(617) 353-8033

		Boston University
		Department of Geography &
		Center for Remote Sensing
		675 Commonwealth Avenue
		Boston, MA   02215  

!Design Notes:

    Program is executed with no command-line options.

Externals:

!END
*******************************************************************************/
{

  /*integer return status*/
  int ier=0;


 /*String for SMF logs*/
  char string[STR_MAX]={};		
  char BegingDate[12];
  char EndingDate[12];
   
    
  /*Declare structure of vars_and_ptrs*/
  struct vars_and_ptrs all={0};					    
    		        
	 
  /*Variables necessary for computing special SDS-level metadata*/



 if(argv!=4) {
    printf("\nUsage: %s <LW_BIN_TILE> <BEGINNINGDATE> <ENDINGDATE>   \n",argc[0]); 
     printf("\nLW_BIN_TILE--Input LW Binary file name (2400x2400)\n");
    printf("BEGINNINGDATE--The beginning date of the LC product period, format='2000-01-01'\n");
    printf("ENDINGDATE--The ending dateof the LC product period, format='2000-12-31'\n\n");
    exit(1);
  }

strcpy(IN_FILE,argc[1]);


strcpy(BegingDate, argc[2]);
strcpy(EndingDate, argc[3]);

  /*Open all required input files*/

  ier=open_files(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "open, hdf2bip.c");
    return(FAILURE);
  }

  /*Initialize the variable structure*/

  ier=init_vars(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, 
	   "init_vars, hdf2bip.c");
    return(FAILURE);
  } 

  /*Get Global ECS metadata from inputfiles */

  ier=get_ecs_met(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "get_ecs_met, hdf2bip.c");
    return(FAILURE); 
  }

 ier=process_ecs_met(&all, BegingDate);
  /***  ier=process_ecs_met(&all, &EndingDate);**/
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "process_ecs_met, hdf2bip.c");
    return(FAILURE); 
  }

  /*Read relevant SDS metadata from input files */

  ier=get_sds_met(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "get_sds_met, hdf2bip.c");
    return(FAILURE); 
  }
 
  /* Get grid metadata */

  ier = get_grid_met (&all);
  if (ier == FAILURE) {
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "get_grid_met, monthly.c");
    return(FAILURE); 
  }
	    

  /*get input array and create output file */
  
  ier=create_bip(&all);
  if (ier == FAILURE) {
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "creat_bip: hdf2bip.c");
    return(FAILURE); 
  }

 /*Attach SDS-level metadata*/

  ier=put_sds_met(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "put_sds_met, hdf2bip.c");
    return(FAILURE);
  }


 /*Attach ECS-level metadata*/

  ier=put_ecs_met(&all,BegingDate, EndingDate);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "put_ecs_met,hdf2bip.c");
    return(FAILURE);
  }	
  sprintf(string, "fixed:  ECS metadata attached");
//  modsmf(MODIS_S_SUCCESS, string, "main: hdf2bip.c");
 
	  
  /*Close all open files*/

  ier=close_files(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "close,hdf2bip.c");
    return(FAILURE);
  }  

  /*Free the allocated memory*/
  
  ier=free_memory(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "free_memory, hdf2bip.c");
    return(FAILURE);
  }  
   
  /*Normal program termination*/
      
  return(SUCCESS);
}


int open_files(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function opens all HDF input files, and creates 
    the HDF output file.  

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of this code (if the required input is missing).

Externals:

               
               
!END
*******************************************************************************/

{


  /*PGS_PC Logical numbers for all input and out files */

  PGSt_integer       VERSION = 1;
  PGSt_SMF_status returnStatus = PGS_S_SUCCESS;

  int i=0;


  char  string[STR_MAX]={};		/*string for SMF logs*/




  /***metadata in LC_temp****/

  returnStatus = PGS_PC_GetNumberOfFiles(ANC_DEM, &all->ptrs.fp.numofiles[1]);


  if (returnStatus != PGS_S_SUCCESS){	
	    
 //   modsmf(MODIS_E_FUNCTION_ERROR, "PGS_PC_GetNumberOfFiles", "hdf2bip.c");
    all->ptrs.fp.fileflag[1]=1;
    return(FAILURE);	/*if there aren't any monthly files at all*/
  }
  else{ 
	
    /* Open monthly input data files for reading */

  
   
      VERSION = 1;

      if ( PGS_PC_GetReference(ANC_DEM, &VERSION, all->ptrs.fp.names[0])
	   != PGS_S_SUCCESS ) {

        sprintf(string, "Logical Id: %d, Version: %d", 
		(int) ANC_DEM, (int) VERSION );

//	modsmf(MODIS_F_GETREF_FAIL, string, "main: hdf2bip.c");

      }
      else{


	all->ptrs.fp.dem_lw_ifp=openMODISfile(all->ptrs.fp.names[0], "r");
    
	if(all->ptrs.fp.dem_lw_ifp==NULL){
		
	  sprintf(string, "openMODISfile, %s, hdf2bip.c", 
		  all->ptrs.fp.names[0]);
	
	//  modsmf(MODIS_N_OPEN_HDF, I_ERR, string);
	 
            return(FAILURE); /*if there aren't any LC_TEMP files at all*/
	}
	    
      }
   
	       
  }



  sprintf(string, "fixed: End of open function");	
//  modsmf(MODIS_S_SUCCESS, string, "main: hdf2bip.c");

    
  return(SUCCESS);
    
}



int init_vars(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function initializes most of the variables in the 
    structure "all" or calls other functions to this end.  

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:
    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf,Xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of this code.

Externals:

               
!END
*******************************************************************************/
{

  int ier=0;					/*int ret. status*/

  char  string[STR_MAX]={};		/*string for SMF logs*/


  /*Initialize total_data_present and total_data_interp*/
    
  all->emet.total_data_present=(unsigned)0;
  all->emet.total_data_interp=(unsigned)0;
  
  /***output*****/    
 

  /*land-water values*/
  strcpy(all->smet.ln_lw, "Land_Water");
  strcpy(all->smet.unit_lw,"flags");

  all->smet.vr_lw[0]=0;
  all->smet.vr_lw[1]=7;
  all->smet.fv_lw=255;

  
strcpy(all->smet.ln_dem_lw, "Land_Water");
  strcpy(all->smet.unit_dem_lw,"flags");
  all->smet.vr_dem_lw[0]=0;
  all->smet.vr_dem_lw[1]=7;
   all->smet.fv_dem_lw=255;


  /*initialize all array name strings*/
    

  /*****LW********/

 strcpy(all->arrnm.arrnm_lw, "LW");
 /***Dynamic_LW****/	
  strcpy(all->arrnm.arrnm_dem_lw, "LW");	
 

    
  /*initialize group name*/

    
  strcpy(all->grpnm.grpnm, "\0");
    

    
  /*initialize SDS ranks*/
    
  ier=init_ranks(all);
    
  if(ier==FAILURE){
	    
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "init_ranks, hdf2bip.c");
	    	
    return(FAILURE);
  }
    
    


    	
  /*Retrieve the MODIS array dimensions for DEM and previous MOD12*/
	
  ier=get_input_ardims(all);
    
  if(ier==FAILURE){
	    
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, 
	   "get_input_ardims, hdf2bip.c");
	    
    return(FAILURE);
  }



    
  sprintf(string, "fixed: End of init_var function"); 
//  modsmf(MODIS_S_SUCCESS, string, "main: hdf2bip.c");
   
			  
    
    
    
    
        
  return(SUCCESS);
    
}




int init_ranks(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function initializes relevant SDS ranks to overestimated 
    values (DIM_MAX) in order to ensure successful retrieval of the 
    actual values via getMODISardims().  

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:
    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:

    This function always returns a successful termination status.

Externals:

               
!END
*******************************************************************************/
{


  register int counter=0;			    /*counter*/



  /**LW***/

 all->rank.rank_lw=DIM_MAX;
 /***DEM_LW***/
all->rank.rank_dem_lw=DIM_MAX;
    
  return(SUCCESS);
}





int get_input_ardims(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function retrieves relevant SDS ranks, dimensions and data types 
    via repeated calls to getMODISardims().  

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu
!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of this code.

Externals:

               
!END
*******************************************************************************/
{
    
  register int sds_counter=0;		/*SDS counter*/
    
  int ier=0;				/*int ret. status*/

  char string[STR_MAX]={};	/*string for SMF logs*/


   
  /*get MODIS array dimensions for each relevant SDS input*/

    
    
  /*Monthly data*/
	
    strcpy(all->dtype.datatype_dem_lw,"uint8");     	    

    ier=getMODISardims(all->ptrs.fp.dem_lw_ifp, 
		       all->arrnm.arrnm_dem_lw, 
		       all->grpnm.grpnm, 
		       all->dtype.datatype_dem_lw, 
		       &(all->rank.rank_dem_lw), 
		       all->dims.dimsizes_dem_lw);

    if(ier==FAILURE){
		
      sprintf(string, "getMODISardims, %s, hdf2bip.c", 
	      all->arrnm.arrnm_lw);
	    
//      modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
		       
      return(FAILURE);
		    
    }
    

	

 /*output files*/
  /*Assign attributes of the output files*/
    
  strcpy(all->dtype.datatype_lw,"uint8");
  all->rank.rank_lw=2;
 
  /*all->dims.dimsizes_lw[0]=all->dims.dimsizes_dem_lw[0];
  all->dims.dimsizes_lw[1]=all->dims.dimsizes_dem_lw[1];*/
  
  /**** modified by Bin Tan for C4-C5 convert ***/
  /*** need change back after the first run of C5   ***/
  all->dims.dimsizes_lw[0]=2400;
  all->dims.dimsizes_lw[1]=2400;
  /*** end by Bin Tan ****/
  sprintf(string, "fixed: End of get_input_ardims function");
//  modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");


 	
  return(SUCCESS);
    
}



int get_ecs_met(struct vars_and_ptrs *all)
/******************************************************************************
!C

!Description: 

    This function retrieves all relevant ECS attributes from the 
    input files.  

!Input Parameters:

    hvars	pointer to structure hdf_vars
    
!Output Parameters:

    hvars, some variable members of this structure may have changed
    
!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:

               
!END
*******************************************************************************/
{

   PGSt_integer       VERSION = 0;
   PGSt_SMF_status returnStatus = PGS_S_SUCCESS;
   char reference[PGSd_PC_FILE_PATH_MAX] = "\0";
   int i=0;
   int ier=0;					/*int return statuses*/
   int inputcntr=0;                             /*counter of inputfiles*/
   char string[STR_MAX]={};			/*SMF string*/
   char datatype[SHORT_STRING_LENGTH] = "\0";
   long nelements = SHORT_STRING_LENGTH * 1L;
   long nelements2 =  MAX_STRING_LENGTH  * 1L;
   long nel = 1L;
   long nel4 = 4L;





   /*Start retrieving metadata from input data */
	

   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "REPROCESSINGPLANNED", datatype,
	      &nelements2, &all->emet.ReprocessingPlanned);


    if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.ReprocessingPlanned);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.ReprocessingPlanned, "Null");
	}
	    


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "REPROCESSINGACTUAL", datatype,
			  &nelements, &all->emet.ReprocessingActual);




    if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.ReprocessingActual);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.ReprocessingActual, "Null");
	}
	

   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "LOCALGRANULEID", datatype,
			  &nelements2, &all->emet.LocalGranuleID);


   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.LocalGranuleID);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.LocalGranuleID, "Null");
	}
	

   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "DAYNIGHTFLAG", datatype,
			  &nelements, &all->emet.DayNightFlag);



   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.DayNightFlag);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.DayNightFlag, "Null");
	}
	


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "LOCALVERSIONID", datatype,
			  &nelements, &all->emet.LocalVersionID);


   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.LocalVersionID);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.LocalVersionID, "Null");

	}
	

   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PGEVERSION", datatype,
			  &nelements, &all->emet.PGEVersion);


   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s,hdf2bip.c", 
		    all->emet.PGEVersion);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.PGEVersion, "Null");


	}
   
   /*INPUTPOINTER --- will have to adjust this when and if 
     the snow product is used*/

  for(i=0;i<all->ptrs.fp.numofiles[1];i++) {
   
    VERSION=(PGSt_integer)(i+1);
    returnStatus = PGS_PC_GetUniversalRef(ANC_DEM,&VERSION,reference);
      if (returnStatus != PGS_S_SUCCESS) {
	 sprintf (string,
	    "Problem with GetUniversalRef with file #%d\n", 
		 ANC_DEM );

//	modsmf(MODIS_E_FUNCTION_ERROR, string, "main: hdf2bip.c");
      }

          else {
	strcpy(all->emet.InputPointer[i],reference);
      }
  }





  /**********************/
for(i=0;i<1;i++) {
   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "RANGEBEGINNINGTIME", datatype,
			  &nelements2, &all->emet.RangeBeginningTimeIn[i]);



   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.RangeBeginningTimeIn[i]);
	
//	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.RangeBeginningTimeIn[i], "Null");
	}



  strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "RANGEBEGINNINGDATE", datatype,
			  &nelements2, &all->emet.RangeBeginningDateIn[i]);



   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.RangeBeginningDateIn[i]);
	
//	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.RangeBeginningDateIn[i], "Null");
	}


  




   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "RANGEENDINGDATE", datatype,
			  &nelements2, &all->emet.RangeEndingDateIn[i]);

   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.RangeEndingDateIn[i]);
	
//	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.RangeEndingDateIn[i], "Null");

	}

   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "RANGEENDINGTIME", datatype,
			  &nelements2, &all->emet.RangeEndingTimeIn[i]);

   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.RangeEndingDateIn[i]);
	
//	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.RangeEndingTimeIn[i], "Null");
	}

  }




   strcpy (all->emet.ExclusionGRingFlag, "N"); /* to comply with ECS SPEC, ktai 4/27/1998 */
   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "EXCLUSIONGRINGFLAG.1", datatype,
			  &nelements, &all->emet.ExclusionGRingFlag);

   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.ExclusionGRingFlag);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.ExclusionGRingFlag, "Null");
	}


   strcpy (datatype, "float64");
   nel4 = 4L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "GRINGPOINTLATITUDE.1", datatype,
			  &nel4, &all->emet.GRingPointLatitude[0]);

   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %lf,hdf2bip.c", 
		    all->emet.GRingPointLatitude[0]);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    all->emet.GRingPointLatitude[0] = (float64) FILL_META_VALUE;
	    all->emet.GRingPointLatitude[1] = (float64) FILL_META_VALUE;
	    all->emet.GRingPointLatitude[2] = (float64) FILL_META_VALUE;
	    all->emet.GRingPointLatitude[3] = (float64) FILL_META_VALUE;
	}


   strcpy (datatype, "float64");
   nel4 = 4L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "GRINGPOINTLONGITUDE.1", datatype,
			  &nel4, &all->emet.GRingPointLongitude[0]);


   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %lf,hdf2bip.c", 
		    all->emet.GRingPointLongitude[0]);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    all->emet.GRingPointLongitude[0] = (float64) FILL_META_VALUE;
	    all->emet.GRingPointLongitude[1] = (float64) FILL_META_VALUE;
	    all->emet.GRingPointLongitude[2] = (float64) FILL_META_VALUE;
	    all->emet.GRingPointLongitude[3] = (float64) FILL_META_VALUE;
	}


   strcpy (datatype, "int32");
   nel4 = 4L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "GRINGPOINTSEQUENCENO.1", datatype,
			  &nel4, &all->emet.GRingPointSequenceNo[0]);
   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c", 
		    all->emet.GRingPointSequenceNo[0]);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
      all->emet.GRingPointSequenceNo[0] = (int32) FILL_META_VALUE;
      all->emet.GRingPointSequenceNo[1] = (int32) FILL_META_VALUE;
      all->emet.GRingPointSequenceNo[2] = (int32) FILL_META_VALUE;
      all->emet.GRingPointSequenceNo[3] = (int32) FILL_META_VALUE;
	}


   strcpy (datatype, "int32");
   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "AUTOMATICQUALITYFLAG.1", datatype,
			  &nelements2, &all->emet.AutomaticQualityFlag);

   if(ier==FAILURE){
	    
	    sprintf(string, "getMODISECSinfo, %s, hdf2bip.c", 
		    all->emet.AutomaticQualityFlag);
	
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);	    
	    strcpy (all->emet.AutomaticQualityFlag, "Null");

	}

   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "AUTOMATICQUALITYFLAGEXPLANATION.1", datatype,
			  &nelements2, 
			  &all->emet.AutomaticQualityFlagExplanation);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AutomaticQualityFlagExplanation);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AutomaticQualityFlagExplanation, "Null");
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			"QAPERCENTINTERPOLATEDDATA.1", datatype,
			  &nel, &all->emet.QAPercentInterpolatedData);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.QAPercentInterpolatedData);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.QAPercentInterpolatedData = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "QAPERCENTMISSINGDATA.1", datatype,
			  &nel, &all->emet.QAPercentMissingData);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.QAPercentMissingData);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.QAPercentMissingData = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "QAPERCENTOUTOFBOUNDSDATA.1", datatype,
			  &nel, &all->emet.QAPercentOutofBoundsData);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.QAPercentOutofBoundsData);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.QAPercentOutofBoundsData = (int32) FILL_META_VALUE;
   }



/*****03/12/01*Crystal suggested to set to zero*******/
all->emet.QAPercentOutofBoundsData=0;


   /*****03/09/01********/

all->emet.QAPercentCloudCover = (int32) FILL_CLOUD;

   /*****03/09/01********/



   /*QAPERCENTGOODQUALITY */
   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.1", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName1);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName1);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName1, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.1", datatype,
			  &nelements, &all->emet.ParameterValue1);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ParameterValue1);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue1, "Null");
   }


   /*QAPERCENTOTHERQUALITY */
   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.2", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName2);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName2);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName2, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.2", datatype,
			  &nelements, &all->emet.ParameterValue2);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ParameterValue2);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue2, "Null");
   }


   /*QAPERCENTNOTPRODUCEDCLOUD */
   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.3", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName3);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName3);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName3, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.3", datatype,
			  &nelements, &all->emet.ParameterValue3);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ParameterValue3);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue3, "Null");
   }


   /*QAPERCENTNOTPRODUCEDOTHER */
   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.4", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName4);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName4);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName4, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.4", datatype,
			  &nelements, &all->emet.ParameterValue4);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ParameterValue4);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue4, "Null");
   }

   /*TILEID*/
   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.5", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName5);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName5);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName5, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.5", datatype,
			  &nelements, &all->emet.ParameterValue5);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	    all->emet.ParameterValue5);  
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue5, "Null");
   }



      /*HORIZONTALTILENUMBER */

   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.6", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName6);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName6);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName6, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.6", datatype,
			  &nelements, &all->emet.ParameterValue6);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ParameterValue6);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue6, "Null");
   }



   /*VERTICALTILENUMBER */

   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "ADDITIONALATTRIBUTENAME.7", datatype,
			  &nelements, 
			  &all->emet.AdditionalAttributeName7);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AdditionalAttributeName7);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AdditionalAttributeName7, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "CoreMetadata.0",
			  "PARAMETERVALUE.7", datatype,
			  &nelements, &all->emet.ParameterValue7);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ParameterValue7);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ParameterValue7, "Null");
   }



   /*ARCHIVE METADATA */

   strcpy (datatype, "float64");
   nel = 1;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "WESTBOUNDINGCOORDINATE", datatype,
			  &nel, &all->emet.WestBoundingCoordinate);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf, hdf2bip.c",
	      all->emet.WestBoundingCoordinate);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.WestBoundingCoordinate = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "float64");
   nel = 1;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "NORTHBOUNDINGCOORDINATE", datatype,
			  &nel, &all->emet.NorthBoundingCoordinate);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf, hdf2bip.c",
	      all->emet.NorthBoundingCoordinate);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.NorthBoundingCoordinate = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "float64");
   nel = 1;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "EASTBOUNDINGCOORDINATE", datatype,
			  &nel, &all->emet.EastBoundingCoordinate);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf, hdf2bip.c",
	      all->emet.EastBoundingCoordinate);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.EastBoundingCoordinate = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "float64");
   nel = 1;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "SOUTHBOUNDINGCOORDINATE", datatype,
			  &nel, &all->emet.SouthBoundingCoordinate);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf, hdf2bip.c",
	      all->emet.SouthBoundingCoordinate);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.SouthBoundingCoordinate = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "ALGORITHMPACKAGEACCEPTANCEDATE", datatype,
			  &nelements2, 
			  &all->emet.AlgorithmPackageAcceptanceDate);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AlgorithmPackageAcceptanceDate);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AlgorithmPackageAcceptanceDate, "Null");
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "ALGORITHMPACKAGEMATURITYCODE", datatype,
			  &nelements2, 
			  &all->emet.AlgorithmPackageMaturityCode);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AlgorithmPackageMaturityCode);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AlgorithmPackageMaturityCode, "Null");
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "ALGORITHMPACKAGENAME", datatype,
			  &nelements2, &all->emet.AlgorithmPackageName);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AlgorithmPackageName);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AlgorithmPackageName, "Null");
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "ALGORITHMPACKAGEVERSION", datatype,
			  &nelements2, 
			  &all->emet.AlgorithmPackageVersion);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.AlgorithmPackageVersion);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.AlgorithmPackageVersion, "Null");
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "GEOANYABNORMAL", datatype,
			  &nelements2, &all->emet.GeoAnyAbnormal);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.GeoAnyAbnormal);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.GeoAnyAbnormal, "Null");
   }


   strcpy (datatype, "float64");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "GEOESTMAXRMSERROR", datatype,
			  &nel, &all->emet.GeoEstMaxRMSError);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf, hdf2bip.c",
	      all->emet.GeoEstMaxRMSError);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.GeoEstMaxRMSError = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "LONGNAME", datatype,
			  &nelements2, &all->emet.LongNameIn);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.LongNameIn);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.LongNameIn, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "SPSOPARAMETERS", datatype,
			  &nelements, &all->emet.SPSOParameters);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.SPSOParameters);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.SPSOParameters, "Null");
   }


   strcpy (datatype, "char *");
   nelements = SHORT_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "PROCESSINGCENTER", datatype,
			  &nelements, &all->emet.ProcessingCenter);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.ProcessingCenter);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.ProcessingCenter, "Null");
   }


   strcpy (datatype, "float64");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
		       "CHARACTERISTICBINANGULARSIZE", datatype,
			  &nel, &all->emet.CharacteristicBinAngularSize);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf, hdf2bip.c",
	      all->emet.CharacteristicBinAngularSize);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.CharacteristicBinAngularSize = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "float64");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "CHARACTERISTICBINSIZE", datatype,
			  &nel, &all->emet.CharacteristicBinSize);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %lf,hdf2bip.c",
	      all->emet.CharacteristicBinSize);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.CharacteristicBinSize = (float64) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "DATACOLUMNS", datatype,
			  &nel, &all->emet.DataColumns);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.DataColumns);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.DataColumns = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "DATAROWS", datatype,
			  &nel, &all->emet.DataRows);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.DataRows);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.DataRows = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "GLOBALGRIDCOLUMNS", datatype,
			  &nel, &all->emet.GlobalGridColumns);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      &all->emet.GlobalGridColumns);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.GlobalGridColumns = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "GLOBALGRIDROWS", datatype,
			  &nel, &all->emet.GlobalGridRows);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.GlobalGridRows);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.GlobalGridRows = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "MAXIMUMOBSERVATIONS", datatype,
			  &nel, &all->emet.MaximumObservations);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.MaximumObservations);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.MaximumObservations = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "int32");
   nel = 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "NUMBEROFGRANULES", datatype,
			  &nel, &all->emet.NumberofGranules);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %ld, hdf2bip.c",
	      all->emet.NumberofGranules);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      all->emet.NumberofGranules = (int32) FILL_META_VALUE;
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "COVERAGECALCULATIONMETHOD", datatype,
			  &nelements2, 
			  &all->emet.CoverageCalculationMethod);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.CoverageCalculationMethod);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.CoverageCalculationMethod, "Null");
   }


   strcpy (datatype, "char *");
   nelements2 = MAX_STRING_LENGTH * 1L;
   ier = getMODISECSinfo (all->ptrs.fp.dem_lw_ifp, "ArchiveMetadata.0",
			  "NADIRDATARESOLUTION", datatype,
			  &nelements2, &all->emet.NadirDataResolution);
   if(ier==FAILURE){
      sprintf(string, "getMODISECSinfo, %s, hdf2bip.c",
	      all->emet.NadirDataResolution);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      strcpy (all->emet.NadirDataResolution, "Null");
   }


       sprintf(string, "fixed: End of get_ecs function");
//       modsmf(MODIS_S_SUCCESS, string, "main: hdf2bip.c");

    
 	
    return(SUCCESS);
   
}






int process_ecs_met(struct vars_and_ptrs *all, char *FAKEENDINGDATE)
/******************************************************************************
!C

!Description: 

    This function processes all relevant ECS attributes.  

!Input Parameters:

    all	pointer to structure hdf_vars
    
!Output Parameters:

    all, some variable members of this structure may have changed
    
!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:

               
!END
*******************************************************************************/
{
  /* A double should be enough to hold a date and time of the format
     YYYYMMDDHHMMSS, which maximally will be something like
     21001231235959 */
  double time[MODERATE_ALLOC]={0};
  int i=0;
  char *str=NULL;
  char buffer[MODERATE_ALLOC];
  char buffer1[MODERATE_ALLOC];
  PGSt_integer julday = 0;
  PGSt_integer julday1 = 0;
  PGSt_integer year = 0;
  PGSt_integer month = 0;
  PGSt_integer day = 0;
  int doy = 0;
  int doy1=0;

  /* Find smallest time stamp of input observation files */
  /* parse and convert to double */
  for (i=0;i<all->ptrs.fp.numofiles[1];i++){
        str = all->emet.RangeBeginningDateIn[i];
 
    strcpy(buffer,"\0");
    while (str[0]!='\0'){
      if(str[0]!='-' && str[0]!=':') strncat(buffer,str,1);
      str++;
    }    
    str = all->emet.RangeBeginningTimeIn[i];
    while (str[0]!='\0' && str[0]!='.'){  /* cut off decimal seconds */
      if(str[0]!='-' && str[0]!=':') strncat(buffer,str,1);
      str++;
    }
    time[i] = atof(buffer);
  } /*end: loop over number of files */
  /* sort */
  all->emet.BeginTime=1.e30;
  all->emet.BeginTimeIndex=0;
  for(i=0;i<all->ptrs.fp.numofiles[1];i++){
    if(all->emet.BeginTime>time[i]){
      all->emet.BeginTime=time[i];
      all->emet.BeginTimeIndex=i;
    }
  }

  /* Store time string of earliest time for inputgranuleID metadata use */
  /* Maybe this should be the latest time???*/
  /****2-26-01****
  strcpy(buffer1,"\0");
  strncat(buffer1,&buffer[0],4);
  year = (PGSt_integer) atoi(buffer1);
  strcpy(buffer1,"\0");
  strncat(buffer1,&buffer[4],2);
  month = (PGSt_integer) atoi(buffer1);
  strcpy(buffer1,"\0");
  strncat(buffer1,&buffer[6],2);
  day = (PGSt_integer) atoi(buffer1);
  julday=PGS_TD_julday(year,month,day);
  julday1=PGS_TD_julday(year,1,1);
  doy = julday - julday1 + 1;
  doy1=(doy-1)/32;
  doy=doy1*32+1;

  if(doy<10) sprintf(all->emet.BeginnDOY,"00%1d",doy);
  else if(doy>=10 && doy<100) sprintf(all->emet.BeginnDOY,"0%2d",doy);
  else sprintf(all->emet.BeginnDOY,"%3d",doy);
  sprintf(all->emet.BeginnYear,"%d",year);
  *****2-26-01**/

  /* Find largest time stamp of input observation files */
  /* parse and convert to double */
  for (i=0;i<all->ptrs.fp.numofiles[1];i++){
    str = all->emet.RangeEndingDateIn[i];
    strcpy(buffer,"\0");
    while (str[0]!='\0'){
      if(str[0]!='-' && str[0]!=':') strncat(buffer,str,1);
      str++;
    }    
    str = all->emet.RangeEndingTimeIn[i];
    while (str[0]!='\0' && str[0]!='.'){  /* cut off decimal seconds */
      if(str[0]!='-' && str[0]!=':') strncat(buffer,str,1);
      str++;
    }
    time[i] = atof(buffer);
  } /*end: loop over number of files */

  /* sort */
  all->emet.EndTime=0.;
  all->emet.EndTimeIndex=all->ptrs.fp.numofiles[1];
  for(i=0;i<all->ptrs.fp.numofiles[1];i++){
    if(all->emet.EndTime<time[i]){
      all->emet.EndTime=time[i];
      all->emet.EndTimeIndex=i;
    }
  }


  /****2-26-01******/


 str = FAKEENDINGDATE;
    strcpy(buffer,"\0");
    while (str[0]!='\0'){
      if(str[0]!='-' && str[0]!=':') strncat(buffer,str,1);
      str++;
    }    

    /***the date is fixed in ".h" file******/

  strcpy(buffer1,"\0");
  strncat(buffer1,&buffer[0],4);
  year = (PGSt_integer) atoi(buffer1);
  strcpy(buffer1,"\0");
  strncat(buffer1,&buffer[4],2);
  month = (PGSt_integer) atoi(buffer1);
  strcpy(buffer1,"\0");
  strncat(buffer1,&buffer[6],2);
  day = (PGSt_integer) atoi(buffer1);
  julday=PGS_TD_julday(year,month,day);
  julday1=PGS_TD_julday(year,1,1);
  doy = julday - julday1 + 1;
  doy1=(doy-1)/96;
  doy=doy1*96+1;

  if(doy<10) sprintf(all->emet.BeginnDOY,"00%1d",doy);
  else if(doy>=10 && doy<100) sprintf(all->emet.BeginnDOY,"0%2d",doy);
  else sprintf(all->emet.BeginnDOY,"%3d",doy);
  sprintf(all->emet.BeginnYear,"%d",year);
  /******** *****2-26-01**/





   return (SUCCESS);
}




int get_sds_met(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function retrieves all SDS metadata from the first monthly input 
    file and from the previous MOD12 product via repeated call to 
    getMODISarinfo.

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:


    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu
!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of this code.

Externals:

               
!END
*******************************************************************************/
{
    
  int ier=0;					/*int return statuses*/

  char string[STR_MAX]={};			/*string for SMF logs*/




  /*Initialize SDS-level metadata variables*/

  ier=init_sds_met(all);
	
  if(ier==FAILURE){
	
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "init_smet, hdf2bip.c");
	    
    return(FAILURE);
  }





   

	    /*long_name*/
  /***DEM_LW***/  

	all->smet.nelementsS=(long)STR_MAX;
	ier=getMODISarinfo(all->ptrs.fp.dem_lw_ifp,
			   all->arrnm.arrnm_dem_lw, 
			   all->grpnm.grpnm,
			   all->smet.ln_name,
			   all->smet.ln_type,
			   &all->smet.nelementsS,
			   all->smet.ln_dem_lw);
	if(ier==FAILURE){
	    
	    sprintf(string, "getMODISarinfo, %s, hdf2bip.c", 
		    all->arrnm.arrnm_dem_lw);
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	    return(FAILURE);
	}
  
    

   
  return(SUCCESS);   
     
}


int init_sds_met(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function initializes all relevant SDS-level
    metadata variables.  

!Input Parameters:

    all		pointer to structure vars_and_ptrs
    
!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:


    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu
!Design Notes:

    This function always returns a successful termination status.

Externals:

               
!END
*******************************************************************************/
{


  /*array dimension labels*/
    
  strcpy(all->smet.attribute[0], "YDim");
  strcpy(all->smet.attribute[1], "XDim");
    
    
  /*long names*/
    
  strcpy(all->smet.ln_name, "long_name");		
  strcpy(all->smet.ln_type, "char *");			
    
    
    
  /*units*/
    
  strcpy(all->smet.unit_name, "units");			
  strcpy(all->smet.unit_type, "char *");		
    
  /*scale factor*/
  strcpy(all->smet.scale_factor_name, "scale_factor");
  strcpy(all->smet.scale_factor_type, "float64");

  /*add_offset*/
  strcpy(all->smet.add_offset_name, "add_offset");
  strcpy(all->smet.add_offset_type, "float64");

  /*add_offset_err*/
  strcpy(all->smet.add_offset_err_name, "add_offset_err");
  strcpy(all->smet.add_offset_err_type, "float64");


  /*calibrated_nt*/
  strcpy(all->smet.calibrated_nt_name, "calibrated_nt");
  strcpy(all->smet.calibrated_nt_type, "int32");


 
  /*valid_range*/
    
  strcpy(all->smet.vr_name, "valid_range");
  strcpy(all->smet.vr_type0, "int8");		
  strcpy(all->smet.vr_type1, "uint8");			 
  strcpy(all->smet.vr_type2, "uint16");			
  strcpy(all->smet.vr_type3, "uint32");			
  strcpy(all->smet.vr_type4, "int16");			
  strcpy(all->smet.vr_type5, "float64");
    
    
    
  /*_FillValue*/
    
  strcpy(all->smet.fv_name, "_FillValue");
  strcpy(all->smet.fv_type0, "int8");	
  strcpy(all->smet.fv_type1, "uint8");			
  strcpy(all->smet.fv_type2, "uint16");			
  strcpy(all->smet.fv_type3, "uint32");
  strcpy(all->smet.fv_type4, "int16");			
  strcpy(all->smet.fv_type5, "float64");
    			

    
  /*Special entries (specific to LC and LCC)*/
 
  /***   
  strcpy(all->smet.spec_ent1, "Percentage of pixels processed completely");
  strcpy(all->smet.spec_ent2, 
	 "Percent changed from previous product generation");
  strcpy(all->smet.spec_ent3, "Proportion of each cover type present");
    
  strcpy(all->smet.spec_type1, "uint8");			
  strcpy(all->smet.spec_type2, "char *");
  ***/ 
    
  /*number of elements*/
    
  all->smet.nelements1=1L;
  all->smet.nelements2=2L;
  all->smet.nelementsS=STR_MAX;
    
    
  /*array dim numbers (for use with getMODISarlabel())*/
    
  all->smet.dim1=1L;
  all->smet.dim2=2L;
    
    
  return(SUCCESS);
   
}



int get_grid_met(struct vars_and_ptrs *all)
/******************************************************************************
!C

!Description: 

    This function retrieves all relevant Grid attributes from the 
    input files.  

!Input Parameters:

    all	pointer to structure hdf_vars


!Output Parameters:

    all, some variable members of this structure may have changed

    
!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

    Portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.


!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Jordan Borak, Feng Gao(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:


Externals:

               
!END



*******************************************************************************/

{



							
    int ier=0;					/*int return statuses*/
    char string[STR_MAX]={};			/*SMF string*/
    char file_name[PGSd_PC_FILE_PATH_MAX] = "\0";
    PGSt_integer       VERSION = 1;
    int32 gfid=0;
    int32 ngrid=0;
    int32 bufsize=MODERATE_ALLOC;
    int32 gid=0;





  /* grid metadata will be copied from the first monthly
     input data set and assumed same for the others */
    /* The question is whether this should be grabbed from 
       MOD12_Prev instead ----- hmmmmm*/



  /* temporarily close first monthly file so it can be
     re-opened with grid tools */
  ier = closeMODISfile(&all->ptrs.fp.dem_lw_ifp);
  if (ier == FAILURE) {
      sprintf (string, "Error temporarily closing observations file\n");
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
  }



  /* open for grid i/o*/
  /* find filename */
   VERSION=(PGSt_integer) 1;


    if ( PGS_PC_GetReference(ANC_DEM, &VERSION, file_name)
	!= PGS_S_SUCCESS ) {
        sprintf(string, "Logical Id: %d, Version: %d", 
		(int) ANC_DEM, (int) VERSION );
//	modsmf(MODIS_F_GETREF_FAIL, string, "get_grid_met: quarterly.c");
    }


  gfid = GDopen(file_name,DFACC_READ);
  if(gfid==GRID_ERRCODE){
      sprintf (string,"Not successful in retrieving grid file ID/open");
	modsmf(MODIS_F_OPEN_HDF_FILE, string, "get_grid_met: quarterly.c");
   }


  /* find out about grid type */
  ngrid=GDinqgrid(file_name,all->emet.GD_gridlist,&bufsize);
  if(ngrid==GRID_ERRCODE){
      sprintf (string,"Not successful in retrieving grid name list");
	modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: monthly.c");
   }



  /* attach grid */
  gid = GDattach(gfid,all->emet.GD_gridlist);
  if(gid==GRID_ERRCODE){
      sprintf (string,"Not successful in attaching grid.");
      modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: quarterly.c");
  }



  /* get projection parameters */
  ier = GDprojinfo(gid,&all->emet.GD_projcode,
        &all->emet.GD_zonecode, &all->emet.GD_spherecode,
        all->emet.GD_projparm);
  if(ier==GRID_ERRCODE){
      sprintf (string,"Not successful in reading grid projection info.");
       modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: quarterly.c");
  }


  /* get grid origin */
  ier = GDorigininfo(gid,&all->emet.GD_origincode);
  if(ier==GRID_ERRCODE){
      sprintf (string,"Failed to read grid origin info."); 
       modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: quarterly.c");
  }



  /* get grid info */
  ier = GDgridinfo(gid,&all->emet.GD_ncols,&all->emet.GD_nrows,
	all->emet.GD_upleft,all->emet.GD_lowright);
  if(ier==GRID_ERRCODE){
      sprintf (string,"Failed to read grid info.");
        modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: quarterly.c");
  }
/**** modified by Bin Tan for C4-C5 convert ***/
  /*** need change back after the first run of C5   ***/
all->emet.GD_ncols=2400;
all->emet.GD_nrows=2400;
/***** done by bin ****/

 /* detach grid */
  ier = GDdetach(gid);
  if(ier==GRID_ERRCODE){
      sprintf (string,"Failed to detach grid.");
       modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: quarterly.c");
  }

  /* close for grid access */
  ier = GDclose(gfid);
  if(ier==GRID_ERRCODE){
      sprintf (string,"GD-file close failed.");
         modsmf(MODIS_E_FUNCTION_ERROR, string, "get_grid_met: quarterly.c");
  }



  /* reopen first reflectance file that was temporarily closed earlier */
   all->ptrs.fp.dem_lw_ifp=openMODISfile(file_name, "r");
	if(all->ptrs.fp.dem_lw_ifp==NULL){
	    sprintf(string, "openMODISfile, %s, quarterly.c", file_name);
	    modsmf(MODIS_F_OPEN_HDF_FILE, F_ERR, string);
	}
	
        sprintf(string, "fixed: End of get_grid_met function");
//	modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");


	return(SUCCESS);    
}











int create_output_arrays(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function creates the SDS's that are written to the output file.     

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon


!Revision History:


    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu
!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of this code.

Externals:

               
!END
*******************************************************************************/
{

  int ier=0;				    /*int ret. status*/

  PGSt_integer       VERSION = 1;
  char file_name[PGSd_PC_FILE_PATH_MAX];
  char string[STR_MAX]={};	    /*string for SMF logs*/
  int kk;
  int32 gfid=0;
  int32 gid=0;
  char dimlist[MAX_STRING_LENGTH];
  strcpy(all->emet.GD_gridlist,"MOD12Q1");

  
  /*Create output SDS's with HDF-EOS for output database*/

  VERSION = 1;

  if ( PGS_PC_GetReference(OUTPUT_HDF, &VERSION, file_name)
       != PGS_S_SUCCESS ) {
    sprintf(string, "Logical Id: %d, Version: %d", 
	    (int) OUTPUT_HDF, (int) VERSION );
    modsmf(MODIS_F_GETREF_FAIL, string, 
	   "create_output_arrays: quarterly.c");
    return(FAILURE);
  }



  /* open */
  gfid = GDopen(file_name,DFACC_CREATE);
  if(gfid==GRID_ERRCODE){
    sprintf (string,"Not successful in retrieving grid file ID/open");
    modsmf(MODIS_F_OPEN_HDF_FILE, string, 
	   "create_output_arrays: quarterly.c");
  }


  /* create grid */
  gid = GDcreate(gfid,all->emet.GD_gridlist,
		 all->emet.GD_ncols, all->emet.GD_nrows,
		 all->emet.GD_upleft,all->emet.GD_lowright);
  if(gfid==GRID_ERRCODE){
    sprintf (string,"Not successful in getting grid ID/ create");
    modsmf(MODIS_E_FUNCTION_ERROR, string, 
	   "create_output_arrays: quarterly.c");
  }




  /* define grid projection */
  ier = GDdefproj(gid,all->emet.GD_projcode,
		  all->emet.GD_zonecode, all->emet.GD_spherecode,
		  all->emet.GD_projparm);
  if(ier==GRID_ERRCODE){
    sprintf (string,"Not successful in defining grid projection");
    modsmf(MODIS_E_FUNCTION_ERROR, string, 
	   "create_output_arrays: quarterly.c");
  }

  /* define grid origin */
  ier = GDdeforigin(gid,all->emet.GD_origincode);
  if(ier==GRID_ERRCODE){
    sprintf (string,"Not successful in defining grid origin");

    modsmf(MODIS_E_FUNCTION_ERROR, string, 
	   "create_output_arrays: quarterly.c");
  }








/* SDS data for LW */
  strcpy(dimlist,all->smet.attribute[0]); /*Note that names YDim, XDim are
					    required*/
  strcat(dimlist,",");
  strcat(dimlist,all->smet.attribute[1]);
  ier = GDdeffield(gid,all->arrnm.arrnm_lw,
		   dimlist,DFNT_UINT8,HDFE_NOMERGE);
  if(ier==GRID_ERRCODE){
    sprintf (string,"Not successful in defining LC Type1 SDS");
    modsmf(MODIS_E_FUNCTION_ERROR, string, 
	   "create_output_arrays: quarterly.c");
  }






  /* detach grid */
  ier = GDdetach(gid);
  if(ier==GRID_ERRCODE){
    sprintf (string,"Unable to detach grid.");
    modsmf(MODIS_E_FUNCTION_ERROR, string, 
	   "create_output_arrays: quarterly.c");
  }


  /* close files for grid access */
  ier = GDclose(gfid);
  if(ier==GRID_ERRCODE){
    sprintf (string,"GD-file close LCQuarterly failed.");
    modsmf(MODIS_E_FUNCTION_ERROR, string, 
	   "create_output_arrays: quarterly.c");
  }



  /*Reopen the output file to sort out pointer to the output file*/
	

  all->ptrs.fp.ofp=openMODISfile(file_name, "a");				    
  if(all->ptrs.fp.ofp==NULL){
    sprintf(string, "openMODISfile, %s, quarterly.c", file_name);
    modsmf(MODIS_F_OPEN_HDF_FILE, F_ERR, string);
    return(FAILURE);
  }
    
    
  sprintf(string, "fixed: End of create_output_arrays function");
//  modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");


    
  return(SUCCESS);

    
}


int close_files(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function closes all HDF input files, and appends 
    relevant header information to the output file.  

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:


    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:
               
!END
*******************************************************************************/
{

 PGSt_MET_all_handles             mdHandles;
  ECSattr_names_for_all_handles   HDFattrnms;
  long                            NumHandles = 0;

  int    files=0;
  char  string[STR_MAX]={};		/*string for SMF logs*/



 

  /**DEM_LW***/
 closeMODISfile(&all->ptrs.fp.dem_lw_ifp);


completeMODISfile(&all->ptrs.fp.ofp, mdHandles, HDFattrnms, NumHandles);


  sprintf(string, "fixed: End of close function"); 
//  modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");


  return(SUCCESS);
    
}






/*create bip file */
  
int create_bip(struct vars_and_ptrs *all)

/******************************************************************************
!C

!Description: 

    read subset definition file and created bip binary file and "names" file which
    can be used in decision tree. 

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure such as landcover type 
    will be assigned to a new classification value
    

!Revision History:


!Team-unique Header:

    This software is developed for NASA under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:
               
!END
*******************************************************************************/

{

  uint16            Data_Rows=0;		 /*Lines per tile*/
  uint16            Data_Columns=0;		 /*Pixels per line*/         
  uint32            line_counter=0;	         /*loop variance for line */
  uint32            cell_counter=0;	         /*loop variance for cell */
  uint32            month_counter=0;             /*loop variance for attribute*/
  char              string[STR_MAX]={};      /*String for SMF logs*/
  int               ier=0;                       /*integer return status*/
  int               TotalFiles;
  PGSt_integer      VERSION=1;
  char              **file_date;                 /*file_date index (YearMonth)*/
  PGSt_SMF_status   returnStatus;
  PGSt_IO_Gen_AccessType   access;
  PGSt_IO_Gen_FileHandle   *F,*S,*DT;  	         /*file pointer for DT "name" file */  
  
 
  int i,j,k,tl,i1;
  int fillvi1,fillvi2; /*zhang 11/16/1999*/  
 int tempdata;
  /* uint16 text_temp;*/
  int16 text_temp;
  int16 bar_temp[7];
  /*  int8 brdf_temp[2];  */
  int16 brdf_temp[8];
/***  int16 vi_temp[2];
      uint16 lst_temp[2];****/
  int16 vi_temp;
  /* uint16 lst_temp;*/
  int16 lst_temp;
  int16 brdf_tempp[8];
  uint8 data[2400];

  int fillvi=32767; /*zhang 11/16/1999*/  

  FILE *in;

  /*initialize local Data_Columns and Data_Rows variables*/

  Data_Columns=(uint16)all->emet.GD_ncols;
  Data_Rows=(uint16)all->emet.GD_nrows;

 
  ier=get_month_order(all,file_date);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, 
	   "get_month_order, quarterly.c");
    return(FAILURE);
  } 

  /***read in binary file***/ 
  if((in=fopen(IN_FILE,"rb"))==NULL) {
    printf("can't open file %s\n",IN_FILE);
    exit(1);
  }
  
/*   printf("\nBip will save in the PCF file sequence \n");  */



  ier=allocate_arrays(all, (uint16)Data_Columns);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR,  "allocate_arrays, quarterly.c");
    return(FAILURE);
  } 

 /*Create the output MODIS arrays*/

  ier=create_output_arrays(all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR,  "create_output_arrays, quarterly.c");
    return(FAILURE);
  }


  ier=set_input_lines(all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "set_input_lines, quarterly.c");
    return(FAILURE);
  }


  /* extracting begin... loop once for each row */
 
  for(line_counter=0;line_counter<(uint32)Data_Rows;
      line_counter++){ 

    /*update appropriate starting points for input MODIS arrays
      (except weights); all dimensions, except dimension 1 always 
      start at zero*/
    
    ier=update_input_starts(all, line_counter);
    if(ier==FAILURE){
      sprintf(string, "update_starts, line %d, quarterly.c", 
	      line_counter);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      continue;
    }

    /*Retrieve the MODIS arrays for monthly 
      databases, plus DEM and previous MOD12*/   	    
  
    /***previous LW***/  
   ier=get_input_arrays(all);
    if(ier==FAILURE){
      sprintf(string, "get_input_arrays, line %d, quarterly.c", 
	      line_counter);
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      continue;
    } 

  /*read NEW LW Mask from binary file */
    fread(data,sizeof(unsigned char),Data_Columns,in);

   
    for(cell_counter=0;cell_counter<(uint32)Data_Columns;
	cell_counter++) {
     
      /**comparison is needed for two different data**/
      /********
                Shallow ocean = 0 ;
                Land (Nothing else but land) = 1;
                Ocean coastlines and lake shorelines = 2 ;
                Shallow inland water = 3 ;
                Ephemeral water = 4 ;
                Deep inland water = 5 ;
                Moderate or continental ocean = 6 ;
                Deep ocean = 7 ;
      **********/

 /*Assign fill values to cell's output values*/
	
  all->ptrs.marray.lw[0][cell_counter]=data[cell_counter];

   }

 ier=put_arrays(all);



 
    
  }
  
  sprintf(string, "fixed: End of create_bip function");
//  modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");
  
  
  return(SUCCESS);
}
  



int get_month_order(struct vars_and_ptrs *all,char **file_date)
/******************************************************************************
!C

!Description: 

    This function get right month info from PCF sequential monthly file.
    Through the relationship of months and files decision tree (*.names) 
    can get the right ordered data from monthly HDF-EOSied files.

    Monthly HDF-EOS file in *.pcf  <==> YearMonth
    DT case data in *.names        <==> YearMonth
  
    SO:  DT case data <==> Monthly HDF-EOS file
    
!Input Parameters:

    all	pointer to structure hdf_vars
    
!Output Parameters:

    file_date, the date of each file
    
!Revision History:

    Version 2: 

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:

               
!END
*******************************************************************************/
{


  return (SUCCESS);
}



int update_input_starts(struct vars_and_ptrs *all, uint32 line_counter)
     /* changed to uint32 by G. Ye, 8/8/96 */
     
/******************************************************************************
       !C

       !Description: 

       This function updates the starting locations in the input/output 
       SDS's from which to retrieve/write relevant information for a 
       given iteration of the code.

       !Input Parameters:

       all		    pointer to structure vars_and_ptrs
    
       line_counter    int value indicating which line is currently 
       being processed by the code 

       !Output Parameters:

       all, some variable members of this structure may have changed
    
       Integer return status qualifying function termination conditon

       !Revision History:

       none

       !Team-unique Header:

       This software is developed by the MODIS Science Data Support
       Team for the National Aeronautics and Space Administration,
       Goddard Space Flight Center, under contract NAS5-32373.

       !References and Credits:

       Written by Jordan S. Borak  

       Boston University
       Department of Geography &
       Center for Remote Sensing
       675 Commonwealth Avenue
       Boston, MA   02215  

       617-353-2088

       borak@crsa.bu.edu

       !Design Notes:

       This function always returns a successful termination status.

       Externals:

               
       !END
      
*******************************************************************************/

{

  
   
  /*MOD12*/
	


  all->start.start_mod12_ncol[0]=(long)line_counter;	
  all->start.start_mod12_icolstart[0]=(long)line_counter;		
  all->start.start_mod12_ncoltile[0]=(long)line_counter;	
  all->start.start_mod12_ipixstart[0]=(long)line_counter;	
    

  /****LW****/
 all->start.start_lw[0]=(long)line_counter;					
 /***DEM_LW***/
all->start.start_dem_lw[0]=(long)line_counter;	
  				

  return(SUCCESS);
}






int get_input_arrays(struct vars_and_ptrs *all)
     
/******************************************************************************
!C

!Description: 

    This function retrieves all input arrays via repeated 
    calls to getMODISarray().

!Input Parameters:

    all -- pointer to structure vars_and_ptrs
    
!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    none

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.


!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu


!Design Notes:

    Unsuccessful termination of this function results in a 
    crash of this code for this iteration.  The code proceeds 
    to the next line of data associated with this tile.

Externals:

               
!END
*******************************************************************************/
{
    
  register int month_counter=0;	    /*month counter*/
  int ier=0;			    /*int ret. status*/
  char string[STR_MAX]={};	    /*string for SMF logs*/


  /***DEM_LW***/


ier=getMODISarray(all->ptrs.fp.dem_lw_ifp, 
		      all->arrnm.arrnm_dem_lw, 
		      all->grpnm.grpnm,
		      all->start.start_dem_lw, 
		      all->dims.dimsizes_dem_lw, 
		      all->ptrs.marray.dem_lw[0]);
		      
    if(ier==FAILURE){

      sprintf(string, "getMODISarray, %s %d, quarterly.c", 
	      all->arrnm.arrnm_dem_lw);
	
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
		
      return(FAILURE);
    }
      


   	
  return(SUCCESS);    
    
}


int set_input_lines(struct vars_and_ptrs *all)
/******************************************************************************
!C

!Description: 

    This function sets the "lines" dimension variable for each SDS 
    to "1" in order to retrieve one line at a time with full 
    dimensionality.  

!Input Parameters:

    all -- pointer to structure vars_and_ptrs

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    none

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu


!Design Notes:

    This function always returns a successful termination status.

Externals:

               
!END
*******************************************************************************/
{

    register int sds_counter=0;		    /*SDS counter*/
    
  
    /*******LW****/

  all->dims.dimsizes_lw[0]=1;

  all->dims.dimsizes_dem_lw[0]=1;

    return(SUCCESS);
}




int allocate_arrays(struct vars_and_ptrs *all, uint16 Data_Columns)
     
/******************************************************************************
!C

!Description: 

    This function calls functions that allocate all dynamic memory 
    required by the code.  Each array allocated is equal in size to 
    one line of data for a given tile with full dimensionality in other 
    dimensions, if they exist.   

!Input Parameters:

    all		     pointer to structure vars_and_ptrs
    
    Data_Columns     the maximum number of valid pixels possible in one 
                     line of this tile

!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    none

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.


!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu



!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of this code.

Externals:

               
!END

*******************************************************************************/

{

  int ier=0;					    /*int ret. status*/

  register int counter=0;			    /*counter*/

  char string[STR_MAX]={};		    /*string for SMF logs*/



  /*Allocate for MODIS data file arrays; 
    on a per-line basis*/
    
    
 

  /***DEM_LW***/

ier=allocate_2d((void ***)&all->ptrs.marray.dem_lw, 1, Data_Columns, 
		  sizeof(uint8));
		    
  if(ier==FAILURE){
	    
    sprintf(string, "allocate_2d, %s quarterly.c", 
	    all->arrnm.arrnm_dem_lw);
		
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
	    
    return(FAILURE);
  }
	



  /******* allocate for LW*****/

 ier=allocate_2d((void ***)&all->ptrs.marray.lw, 1, Data_Columns, 
		  sizeof(uint8));
		    
  if(ier==FAILURE){
	    
    sprintf(string, "allocate_2d, %s quarterly.c", 
	    all->arrnm.arrnm_lw);
		
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
	    
    return(FAILURE);
  }
	


   
  return(SUCCESS);
}










int put_arrays(struct vars_and_ptrs *all)
/******************************************************************************
!C

!Description: 

    This function puts all data in the output SDS's for 
    a given tile line via repeated calls to putMODISarray().

!Input Parameters:

    all -- pointer to structure vars_and_ptrs
    
!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    none

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.


!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf,Xiaoyang ZHAng, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

  

!Design Notes:

    Unsuccessful termination of this function results in a 
    crash of the code for this iteration.  The code proceeds 
    to the next line of data associated with this tile.

Externals:

               
!END
*******************************************************************************/
{
    
    int ier=0;				/*int ret. status*/
    char string[STR_MAX]={};	/*string for SMF logs*/

   
  /*LW arrays*/
        
    ier=putMODISarray(all->ptrs.fp.ofp, 
                      all->arrnm.arrnm_lw, 
                      all->grpnm.grpnm, 
		      all->start.start_lw, 
		      all->dims.dimsizes_lw, 
		      all->ptrs.marray.lw[0]);
			      
	if(ier==FAILURE){
	    sprintf(string, "putMODISarray, %s, quarterly.c", 
		    all->arrnm.arrnm_lw);
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	    return(FAILURE);
	}
	
 	
	
    
    return(SUCCESS);
    
    
}












int put_sds_met(struct vars_and_ptrs *all)
/******************************************************************************
!C

!Description: 

    This function attaches all relevant attributes to the 
    output SDS's.  

!Input Parameters:

        all -- pointer to structure hdf_vars
    
!Output Parameters:

    hvars, some variable members of this structure may have changed

!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:

               
!END
*******************************************************************************/
{

    int ier=0;					/*int return statuses*/
    int i=0;

    char string[STR_MAX]={};			/*SMF string*/
    
    


    
    
    /*Write SDS-level metadata*/

    
    
	    /*long_name*/
	   
       
	  all->smet.nelementsS=(long)STR_MAX;
 
	  ier=putMODISarinfo(all->ptrs.fp.ofp,
			   all->arrnm.arrnm_lw, 
			   all->grpnm.grpnm,
			   all->smet.ln_name,
			   all->smet.ln_type,
			   all->smet.nelementsS,
			   all->smet.ln_lw);
	  if(ier==FAILURE){
	    
	    sprintf(string, "putMODISarinfo, %s, quarterly.c", 
		    all->arrnm.arrnm_lw);
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	    return(FAILURE);
	  }
    
	


	/*units*/


     

	  all->smet.nelementsS=(long)STR_MAX;
 
	  ier=putMODISarinfo(all->ptrs.fp.ofp,
			   all->arrnm.arrnm_lw, 
			   all->grpnm.grpnm,
			   all->smet.unit_name,
			   all->smet.unit_type,
			   all->smet.nelementsS,
			   all->smet.unit_lw);
	  if(ier==FAILURE){
	    
	    sprintf(string, "putMODISarinfo, %s, quarterly.c", 
		    all->arrnm.arrnm_lw);
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	    return(FAILURE);
	  }
    
	


	/*valid range*/


	all->smet.nelements2=2L;
        ier=putMODISarinfo(all->ptrs.fp.ofp,
			   all->arrnm.arrnm_lw, 
			   all->grpnm.grpnm,
			   all->smet.vr_name,
			   all->smet.vr_type1,
			   all->smet.nelements2,
			   all->smet.vr_lw);
	  if(ier==FAILURE){
	 
	    sprintf(string, "putMODISarinfo, %s, quarterly.c", 
		    all->arrnm.arrnm_lw);
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	    return(FAILURE);
	  }
    



	/*fill values*/

	all->smet.nelements1=1L;
        ier=putMODISarinfo(all->ptrs.fp.ofp,
			   all->arrnm.arrnm_lw, 
			   all->grpnm.grpnm,
			   all->smet.fv_name,
			   all->smet.fv_type1,
			   all->smet.nelements1,
			   &all->smet.fv_lw);
	  if(ier==FAILURE){
	 
	    sprintf(string, "putMODISarinfo, %s, quarterly.c", 
		    all->arrnm.arrnm_lw);
	    modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	    return(FAILURE);
	  }
    


  
  sprintf(string, "fixed: End of put_sds function");
//  modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");




    return(SUCCESS);
   
}









int put_ecs_met(struct vars_and_ptrs *all, char *FAKEBEGINNINGDATE, char *FAKEENDINGDATE)
/******************************************************************************
!C

!Description: 

    This function puts all ECS-level metadata in the 
    output file.  

!Input Parameters:

    all -- pointer to structure hdf_vars
    
!Output Parameters:

    all, some variable members of this structure may have changed

!Revision History:

    Version 2: 
    Original version prepared as V1 delivery 
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Xiaoyang Zhang, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:

               
!END
*******************************************************************************/
{		    
    
    int ier=0;					/*integer return status*/
    int i=0;
    int totalinputfiles=0;
    register int j=0;				/*counters*/
    char string[STR_MAX]={};		/*string for SMF logs*/

 
   char *valuestr[HUGE_ALLOC] = {};	   
   char *metastring=NULL;
   float64 *metavalue=NULL;
   int32 *metavalueint=NULL;
   int32 metavalue4int[1000] = {0};
   float64 metavalue4[1000] = {0.};
    PGSt_MET_all_handles mdHandles={};	/*PGS metadata handles*/
   char str[MAX_STRING_LENGTH] = "\0";
   char timestr[MAX_STRING_LENGTH] = "\0";
    struct tm *currtime = NULL;
   time_t t;

 
    /*Set those ECS metadata that are changed in this code*/



   /* General comment:
      in the following, some values are first converted to strings
      using sprintf; then they are copied to the output string
      using strcpy; this seems like an unnecessary 2-step procedure;
      however, if one uses sprintf to print directly to the output
      string, the metadata writes encounter segmentation faults
      in the toolkit; something's going on here. */

    /*CORE METADATA*/

   strcpy(all->emet.ReprocessingPlanned,"no further update anticipated");
   /* XL: following needs to be made active, but how? */
   strcpy(all->emet.ReprocessingActual,"processed once");
   strcpy(all->emet.LocalVersionID,"2.2.0");
   strcpy (all->emet.PGEVersion, "2.3.0");

   /***   strcpy (all->emet.RangeBeginningDate,
	  all->emet.RangeBeginningDateIn[all->emet.BeginTimeIndex]);****/
 strcpy (all->emet.RangeBeginningDate,FAKEBEGINNINGDATE);
 /**2-26-2001**/

   strcpy (all->emet.RangeBeginningTime, 
      all->emet.RangeBeginningTimeIn[all->emet.BeginTimeIndex]);
   /**   strcpy (all->emet.RangeEndingDate, 
	 all->emet.RangeEndingDateIn[all->emet.EndTimeIndex]);**/

 strcpy (all->emet.RangeEndingDate,FAKEENDINGDATE);
 /**2-26-2001**/

   strcpy (all->emet.RangeEndingTime, 
      all->emet.RangeEndingTimeIn[all->emet.EndTimeIndex]);

   strcpy (all->emet.ParameterName,"LW, NBAR, NBAR_QC, Text, Text_QC, BRDF, BRDF_QC, VI, VI_QC, SNOW, SNOW_QC, LST, LST_QC");

   strcpy (all->emet.AutomaticQualityFlag,"Passed");
   strcpy (all->emet.AutomaticQualityFlagExplanation,
      "To be set as 'passed' or 'failed' to indicate failure of PGE test.");





   /*ARCHIVE METADATA*/


   strcpy (all->emet.LongNameOut,
	   "MODIS Quarterly 1km Land Cover");


   strcpy(all->emet.SPSOParameters,"2669");    /* with ECS SPEC, Ktai, 4/27/1998 */ 
   /** strcpy(all->emet.ProcessingCenter,"EDC");**/
   strcpy(all->emet.ProcessingCenter,"MODAPS");




    /*Initialize metadata tools*/
    
    if(PGS_MET_Init(MCF_FILE, mdHandles)!=PGS_S_SUCCESS){
    
	sprintf(string, "PGS_MET_Init, quarterly.c");
	
	modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
	
	return(FAILURE);
    }


      /* write ECS metadata */

 

     /* Reprocessing*/

      /* value set in code */
      metastring =  all->emet.ReprocessingPlanned;
      ier = PGS_MET_SetAttr (mdHandles[1],"ReprocessingPlanned",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.ReprocessingPlanned);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


     /* value set in code */
      metastring =  all->emet.ReprocessingActual;
      ier = PGS_MET_SetAttr (mdHandles[1],"ReprocessingActual",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		all->emet.ReprocessingActual); 
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /*LocalGranuleID*/


    /* value set here */
     strcpy(str,ESDTNAME);
     strcat(str,".");
     strcat(str,"A");
     strcat(str,all->emet.BeginnYear); /* Granule beginning time */
     strcat(str,all->emet.BeginnDOY);
     strcat(str,".");  /* Modified to comply with ECS SPEC, Ktai 4/27/1998 */ 
     strcat(str,"h");
     strcat(str,all->emet.ParameterValue5);  /*Horizontal Tile Number*/
     strcat(str,"v");
     strcat(str,all->emet.ParameterValue6); /*Vertical Tile Number*/
     strcat(str,".");
     strcat(str,"002");
     strcat(str,".");
     if((t = time(NULL)) == -1){  /* production time */
       sprintf (string,"time() failed, affects LocalGranuleID metadata.");
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
     }
     if(!(currtime = (struct tm *)gmtime(&t))){
       /*    if(!(currtime = (struct tm *)localtime(&t))){*/
       sprintf (string,"gmtime() failed, affects LocalGranuleID metadata.");
        modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
     }
     strftime(timestr,100,"%Y%j%H%M%S",currtime);
     strcat(str,timestr);
     strcat(str,".hdf");
     strcpy(all->emet.LocalGranuleID,str);

      metastring =  all->emet.LocalGranuleID;
      ier = PGS_MET_SetAttr (mdHandles[1],"LocalGranuleID",
			     &metastring);
      if (ier != PGS_S_SUCCESS) {
	 sprintf (string,"Error completing PGS_MET_SetAttr"
		  " LocalGranuleID");
	  modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }




      metastring =  all->emet.LocalGranuleID;
      ier = PGS_MET_SetAttr (mdHandles[1],"LocalGranuleID",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.LocalGranuleID);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }






      /* value as read from input */
      metastring =  all->emet.DayNightFlag;
      ier = PGS_MET_SetAttr (mdHandles[1],"DayNightFlag",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.DayNightFlag);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


     /* value as set by code */
      metastring =  all->emet.LocalVersionID;
      ier = PGS_MET_SetAttr (mdHandles[1],"LocalVersionID",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.LocalVersionID);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as set by code */
      metastring =  all->emet.PGEVersion;
      ier = PGS_MET_SetAttr (mdHandles[1],"PGEVersion",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.PGEVersion);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }




      /* value as read from all the inputs */
      for(i=0;i<2;i++){
      totalinputfiles=totalinputfiles+all->ptrs.fp.numofiles[i];
      }
      for (j = 0; j < totalinputfiles; j++) {
	 valuestr[j] = calloc (MAX_STRING_LENGTH, sizeof (char));
	 strcpy(valuestr[j],all->emet.InputPointer[j]);
      }
      ier = PGS_MET_SetAttr (mdHandles[1],"InputPointer",
			     valuestr);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.InputPointer);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }
      for (j = 0; j < totalinputfiles; j++) {
	free(valuestr[j]);
      }



      /* value as read from input */
      metastring =  all->emet.RangeBeginningDate;
      ier = PGS_MET_SetAttr (mdHandles[1],"RangeBeginningDate",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.RangeBeginningDate);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as read from input */
      metastring =  all->emet.RangeBeginningTime;
      ier = PGS_MET_SetAttr (mdHandles[1],"RangeBeginningTime",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.RangeBeginningTime);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as computed by code from input */
      metastring =  all->emet.RangeEndingDate;
      ier = PGS_MET_SetAttr (mdHandles[1],"RangeEndingDate",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.RangeEndingDate);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as computed by code from input */
      metastring =  all->emet.RangeEndingTime;
      ier = PGS_MET_SetAttr (mdHandles[1],"RangeEndingTime",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.RangeEndingTime); 
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as read from input */
      metastring =  all->emet.ExclusionGRingFlag;
      ier = PGS_MET_SetAttr (mdHandles[1],"ExclusionGRingFlag.1",
			      &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.ExclusionGRingFlag);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as read from input */
      for (j = 0; j < 4; j++) {
	 metavalue4[j] = all->emet.GRingPointLatitude[j];
      }
      ier = PGS_MET_SetAttr (mdHandles[1],"GRingPointLatitude.1",
			     metavalue4);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c", 
		 all->emet.GRingPointLatitude[0]);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as read from input */
      for (j = 0; j < 4; j++) {
	 metavalue4[j] = all->emet.GRingPointLongitude[j];
      }
      ier = PGS_MET_SetAttr (mdHandles[1],"GRingPointLongitude.1",
			     metavalue4);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c", 
		  all->emet.GRingPointLongitude[0]);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* value as read from input */
      for (j = 0; j < 4; j++) {
	 metavalue4int[j] = 
	   /**  all->emet.GRingPointSequenceNo[j]+1;**/ /* comply with ECS SPEC, ktai 4/28/1998 */ 
	   /********keep looking******/
	   all->emet.GRingPointSequenceNo[j];  /***?????04/30/01 ***/

	 /**It should be (1,2,3,4) told by TK, so in TEMPORAL FILE, IT SHOULD BE (0,1,2,3)***/
/*	 metavalue4int[j] = all->emet.GRingPointSequenceNo[j]; */
      }
      ier = PGS_MET_SetAttr (mdHandles[1],"GRingPointSequenceNo.1",
			     metavalue4int);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c", 
		 all->emet.GRingPointSequenceNo[0]);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }



      strcpy (metastring, all->emet.ParameterName);
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterName.1",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.ParameterName); 
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As set in code */ 
      metastring =  all->emet.AutomaticQualityFlag;
      ier = PGS_MET_SetAttr (mdHandles[1],"AutomaticQualityFlag.1",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.AutomaticQualityFlag);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As set in code */ 
      metastring =  all->emet.AutomaticQualityFlagExplanation;
      ier = PGS_MET_SetAttr (mdHandles[1],"AutomaticQualityFlagExplanation.1",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.AutomaticQualityFlagExplanation);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As passed from input*/ 
      metavalueint = &all->emet.QAPercentInterpolatedData;
      ier = PGS_MET_SetAttr (mdHandles[1],"QAPercentInterpolatedData.1",
			     metavalueint);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c", 
		 all->emet.QAPercentInterpolatedData);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As passed from input*/ 
      metavalueint = &all->emet.QAPercentMissingData;
      ier = PGS_MET_SetAttr (mdHandles[1],"QAPercentMissingData.1",
			     metavalueint);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
		 all->emet.QAPercentMissingData);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As passed from input)*/ 
      metavalueint = &all->emet.QAPercentOutofBoundsData;
      ier = PGS_MET_SetAttr (mdHandles[1],"QAPercentOutofBoundsData.1",
			     metavalueint);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c", 
		 all->emet.QAPercentOutofBoundsData);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


   /******03/09/01**cloud***/
 metavalueint = &all->emet.QAPercentCloudCover;
      ier = PGS_MET_SetAttr (mdHandles[1],"QAPercentCloudCover.1",
			     metavalueint);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c", 
		 all->emet.QAPercentCloudCover);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }
 /******03/09/01**cloud***/


      /*QAPERCENTGOODQUALITY */
      metastring =  all->emet.AdditionalAttributeName1;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.1",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.AdditionalAttributeName1);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As set by nbars (or vi) */ 
      metastring =  all->emet.ParameterValue1;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.1",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.ParameterValue1);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /*QAPERCENTOTHERQUALITY */
      metastring =  all->emet.AdditionalAttributeName2;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.2",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.AdditionalAttributeName2);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As set by nbars (or vi) */ 
      metastring =  all->emet.ParameterValue2;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.2",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		 all->emet.ParameterValue2);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /*QANOTPRODUCEDCLOUD */
      metastring =  all->emet.AdditionalAttributeName3;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.3",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.AdditionalAttributeName3);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As set by nbars (or vi)*/ 
      metastring =  all->emet.ParameterValue3;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.3",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		 all->emet.ParameterValue3);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /*QANOTPRODUCEDOTHER */
      metastring =  all->emet.AdditionalAttributeName4;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.4",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		 all->emet.AdditionalAttributeName4);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* As set by nbars (or vi) */ 
      metastring =  all->emet.ParameterValue4;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.4",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		  all->emet.ParameterValue4);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


     


 /*HORIZONTALTILENUMBER */

      metastring =  all->emet.AdditionalAttributeName5;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.5",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		  all->emet.AdditionalAttributeName5);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* as read from input */
      metastring =  all->emet.ParameterValue5;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.5",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		 all->emet.ParameterValue5);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }



      /*VERTICALTILENUMBER    */

      metastring =  all->emet.AdditionalAttributeName6;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.6",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		  all->emet.AdditionalAttributeName6);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }


      /* as read from input */
      metastring =  all->emet.ParameterValue6;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.6",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		 all->emet.ParameterValue6);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }



 /*TileID*/

      metastring =  all->emet.AdditionalAttributeName7;
      ier = PGS_MET_SetAttr (mdHandles[1],"AdditionalAttributeName.7",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c", 
		  all->emet.AdditionalAttributeName7);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }

    
      /* as read from input */
      metastring =  all->emet.ParameterValue7;
      ier = PGS_MET_SetAttr (mdHandles[1],"ParameterValue.7",
			     &metastring);
      if(ier==FAILURE){
	 sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
		  all->emet.ParameterValue7);
	 modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      }




      /*ARCHIVE METADATA*/



	/*As passed from input (nbars or vi)*/
     metavalue = &all->emet.WestBoundingCoordinate;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "WestBoundingCoordinate", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
         all->emet.WestBoundingCoordinate );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }



	/*As passed from input (nbars or vi)*/
     metavalue = &all->emet.NorthBoundingCoordinate;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "NorthBoundingCoordinate", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
         all->emet.NorthBoundingCoordinate );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }



	/*As passed from input (nbars or vi)*/
     metavalue = &all->emet.EastBoundingCoordinate;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "EastBoundingCoordinate", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
         all->emet.EastBoundingCoordinate );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }



	/*As passed from input (nbars or vi)*/
     metavalue = &all->emet.SouthBoundingCoordinate;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "SouthBoundingCoordinate", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
         all->emet.SouthBoundingCoordinate );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }


	/*As passed from input*/
     metastring = (void *) all->emet.GeoAnyAbnormal;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "GeoAnyAbnormal", 
                     &metastring);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
          all->emet.GeoAnyAbnormal);
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	

	/*As passed from input*/
     metavalue = &all->emet.GeoEstMaxRMSError;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "GeoEstMaxRMSError", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
          all->emet.GeoEstMaxRMSError);
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }


	/*As set in code*/
     metastring = (void *) all->emet.LongNameOut;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "LongName", 
                     &metastring);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
         all->emet.LongNameOut );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	


	/*As set in code*/
     metastring = (void *) all->emet.SPSOParameters;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "SPSOParameters", 
                     &metastring);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
         all->emet.SPSOParameters );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	


	/*As set in code*/
     metastring = (void *) all->emet.ProcessingCenter;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "ProcessingCenter", 
                     &metastring);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
         all->emet.ProcessingCenter );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	


	/*As passed from input*/
    metavalue = &all->emet.CharacteristicBinAngularSize;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "CharacteristicBinAngularSize", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
         all->emet.CharacteristicBinAngularSize );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }

	

	/*As passed from input*/
    metavalue = &all->emet.CharacteristicBinSize;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "CharacteristicBinSize", 
                          metavalue);  
        if(ier==FAILURE){
        sprintf(string, "PGS_MET_SetAttr, %lf, quarterly.c",
          all->emet.CharacteristicBinSize);
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
        }



	/*As passed from input*/
     metavalueint = &all->emet.DataColumns;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "DataColumns", 
			  metavalueint);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
          all->emet.DataColumns);
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	

	/*As passed from input*/
     metavalueint = &all->emet.DataRows;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "DataRows", 
			  metavalueint);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
          all->emet.DataRows);
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	


	/*As passed from input*/
     metavalueint = &all->emet.GlobalGridColumns;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "GlobalGridColumns", 
			  metavalueint);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
         all->emet.GlobalGridColumns );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}
	

	/*As passed from input*/
     metavalueint = &all->emet.GlobalGridRows;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "GlobalGridRows", 
			  metavalueint);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
        all->emet.GlobalGridRows );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}



	/*As passed from input*/
     metavalueint = &all->emet.MaximumObservations;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "MaximumObservations", 
			  metavalueint);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
         all->emet.MaximumObservations );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}


	/*As passed from input*/
     metavalueint = &all->emet.NumberofGranules;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "NumberofGranules", 
			  metavalueint);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %ld, quarterly.c",
         all->emet.NumberofGranules );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}




	/*As passed from input*/
     metastring = (void *) all->emet.CoverageCalculationMethod;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "CoverageCalculationMethod", 
                     &metastring);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
          all->emet.CoverageCalculationMethod );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}


	/*As passed from input*/
     metastring = (void *) all->emet.NadirDataResolution;
     ier=PGS_MET_SetAttr( mdHandles[2], 
                     "NadirDataResolution", 
                     &metastring);  
        if(ier==FAILURE){
	sprintf(string, "PGS_MET_SetAttr, %s, quarterly.c",
         all->emet.NadirDataResolution );
       modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
	}






    /*write metadata*/


                                             /* cast by G. Ye, 8/15/96 */
	ier= PGS_MET_Write(mdHandles[1], "CoreMetadata.0", 
			   (PGSt_integer)all->ptrs.fp.ofp->sd_id);
	/*XL not sure how to switch this to a warning-- 
	  presently complains about Mandatory metadata*/
	/*
	  if (ier!=PGS_S_SUCCESS){
	
	  sprintf(string, "PGS_MET_Write, quarterly.c");
	  
	  modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
	
	  return(FAILURE);
	
	  }*/
  

	ier= PGS_MET_Write(mdHandles[2], "ArchiveMetadata.0", 
			   (PGSt_integer)all->ptrs.fp.ofp->sd_id);
	/*XL not sure how to switch this to a warning-- 
	  presently complains about Mandatory metadata*/
	/*
	  if(ier!=PGS_S_SUCCESS){
	  
	  sprintf(string, "PGS_MET_Write, quarterly.c");
	  
	  modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    
	  return(FAILURE);
	  }*/
	
      
    /*Free MET-allocated memory*/
 /*XL: Do I need this?*/
 /*    PGS_MET_Remove();*/



    
    	
    return(SUCCESS);
    
		
}




PGSt_SMF_status current_time_a(char *buffer)
/******************************************************************************
!C

!Description:
  
    This function generates the current system time in ASCII Time
    Code A format.  Function prototype is in header file smfio.h.

!Input Parameters:

    char *buffer  Character buffer allocated to TIMECODEASIZE*sizeof(char)

!Output Parameters:

    char *buffer  Character buffer containing current time in ASCII Time Code A

Returns:
  PGSt_SMF_status
    on success:    MODIS_S_SUCCESS
    on error:      MODIS_E_BAD_SYSTEM_TIME
                   MODIS_E_NULL_STRING

!Revision History:
  Revision 1.0  1996/03/08
  Paul Fisher/SDST
  Original Development

  Revision 1.1  1996/03/15
  Paul Fisher/SDST
  Efficiency modification.
  Modified to use only one call to strftime. 

!Team-unique Header:
  This software was developed by:

    MODIS Science Data Support Team for the National Aeronautics and Space 
    Administration, Goddard Space Flight Center, under contract NAS5-32373.

  Developer:
      Paul S. Fisher
      MODIS Science Data Support Team   
      Research and Data Systems Corporation
      SAIC/GSC MODIS Support Office
      7501 Forbes Blvd
      Seabrook, MD 20706  
      fisher@ltpmail.gsfc.nasa.gov     

!Design Notes:
  This module returns the current system time in ASCII Time Code A
  format.  It is designed to be used in the SDP Toolkit environment.

  Error handling is performed through the modules smf.c and the
  toolkit SMF seed file MODIS_37434.t.

  As input, this module takes a character string of length
  TIMECODEASIZE, which has been defined in smfio.h, a necessary header
  file for both this module and smf.c

  Upon error, this function will return to the calling module to
  allow exection to continue.

  If 'buffer' contains less than TIMECODEASIZE characters, module
  behavior is undefined.

!END***************************************************************************/
{

  struct tm *curtime=NULL;			    /*ptr to type tm*/
  time_t t=0;					    /*time_t return value*/

  /* check to make sure buffer has some memory allocated */
  if(!buffer){
    modsmf(MODIS_E_NULL_STRING, 
	   "buffer", 
	   "current_time_a, timea.c");
    return (MODIS_E_NULL_STRING);
  }

  /* set time */
  if((t=time(NULL))==-1){
    modsmf(MODIS_E_BAD_SYSTEM_TIME, 
	   "Bad time value returned from operating system", 
	   "current_time_a, timea.c");
    return (MODIS_E_BAD_SYSTEM_TIME);
  }

  if(!(curtime=(struct tm *)localtime(&t))){
    modsmf(MODIS_E_BAD_SYSTEM_TIME, 
	   "No time structure filled in call to 'localtime'.", 
	   "current_time_a, timea.c");
    return (MODIS_E_BAD_SYSTEM_TIME);
  }

  /* obtain CCSDS ASCII time code A */
  if((strftime(buffer, TIMECODEASIZE, "%Y-%m-%dT%H:%M:%S.000000Z", curtime))
     ==0){
    modsmf(MODIS_E_BAD_SYSTEM_TIME, "Invalid time structure "
	   "returned from operating system.", 
	   "current_time_a, timea.c");
    return (MODIS_E_BAD_SYSTEM_TIME);
  }

  return(MODIS_S_SUCCESS);
   
}















/*Dynamic memory allocation functions*/

int allocate_1d(void **i_ptr, uint16 dim1, int elsize)
/******************************************************************************
!C

!Description:

    This function allocates all memory requested as a
    contiguous vector in compliance with M-API requirements.   

!Input Parameters:

    i_ptr   pointer to a pointer to type void
    
    elsize  int size of data type for which memory is being allocated
    
    dim1    uint16 number of elements of size "elsize"
    
!Output Parameters:

    i_ptr, memory may have been allocated to this pointer
    
    Integer return status qualifying function termination conditon

!Revision History:

    This function was written by John B. Collins, Boston University, 7/95

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:

    Written by John B. Collins, Boston University
    
    Contact:   Jordan S. Borak  

		Boston University
		Department of Geography &
		Center for Remote Sensing
		675 Commonwealth Avenue
		Boston, MA   02215  

		617-353-2088

		borak@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of the code.

!Externals:

               
!END
******************************************************************************/
{
    
    void *ptr1=NULL;                            /*pointer to 1-D array*/
    ptr1=calloc((size_t)dim1, (size_t)elsize);
                                                /* cast by G. Ye, 8/8/96 */
    if(!ptr1){
        modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "allocate_1d, allocate.c");
	return(FAILURE);
    }
    *i_ptr=ptr1;
    return(SUCCESS);
}



int allocate_2d(void ***i_ptr, uint16 dim1, uint16 dim2, int elsize)
/******************************************************************************
!C

!Description:

    This function allocates all memory requested in a contiguous block 
    via a call to allocate_1d().  It then allocates a vector of pointers 
    to this memory via another call to allocate_1d() for a two dimensional 
    allocation.

!Input Parameters:

    i_ptr   pointer to a pointer to a pointer to type void
    
    elsize  int size of data type for which memory is being allocated
    
    dim1    uint16 number of elements of size "elsize" for first dimension
    
    dim2    uint16 number of elements of size "elsize" for second dimension
    
!Output Parameters:

    i_ptr, memory may have been allocated to this pointer
    
    Integer return status qualifying function termination conditon

!Revision History:

    This function was written by John B. Collins, Boston University, 7/95

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:

    Written by John B. Collins, Boston University
    
    Contact:   Jordan S. Borak  

		Boston University
		Department of Geography &
		Center for Remote Sensing
		675 Commonwealth Avenue
		Boston, MA   02215  

		617-353-2088

		borak@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of the code.

!Externals:

               
!END
*******************************************************************************/
{
    
    int ier=0;					/*error code*/
    void *ptr1=NULL;                            /*pointer to 1-D array*/
    void **ptr2=NULL;                           /*pointer to 2-D array*/
    register int counter=0;                     /*dimension counter*/

    /*allocate memory for data matrix */
    
    ier=allocate_1d((void **)&ptr1, (uint16)(dim1*dim2), elsize);
                                                /* cast by G. Ye, 8/8/96 */
	if(ier==FAILURE){
	   modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "allocate_1d, allocate.c");
	   return(FAILURE); 
	}

    /*allocate array of void ptrs for dimension 1 */
        
    ier=allocate_1d((void **)&ptr2, dim1, sizeof(void *));
	if(ier==FAILURE){
	   modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "allocate_1d, allocate.c");
	   return(FAILURE); 
	}     
    
    /*point allocated pointers at 
    appropriate data elements*/

    for(counter=0;counter<dim1;counter++){      
        ptr2[counter]=(char *)ptr1+(counter*dim2*elsize);
    }
    *i_ptr=ptr2;
    return(SUCCESS);
}


int allocate_3d(void ****i_ptr, uint16 dim1, uint16 dim2, uint16 dim3, 
                int elsize)
/******************************************************************************
!C

!Description:

    This function allocates all memory requested in a contiguous block 
    via a call to allocate_1d().  It then allocates vectors of pointers 
    to this memory via repeated calls to allocate_1d() for a three 
    dimensional allocation.

!Input Parameters:

    i_ptr   pointer to a pointer to a pointer to a pointer to type void
    
    elsize  int size of data type for which memory is being allocated
    
    dim1    uint16 number of elements of size "elsize" for first dimension
    
    dim2    uint16 number of elements of size "elsize" for second dimension
    
    dim3    uint16 number of elements of size "elsize" for third dimension
    
!Output Parameters:

    i_ptr, memory may have been allocated to this pointer
    
    Integer return status qualifying function termination conditon

!Revision History:

    This function was written by John B. Collins, Boston University, 7/95

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:

    Written by John B. Collins, Boston University
    
    Contact:   Jordan S. Borak  

		Boston University
		Department of Geography &
		Center for Remote Sensing
		675 Commonwealth Avenue
		Boston, MA   02215  

		617-353-2088

		borak@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of the code.

!Externals:

               
!END
*******************************************************************************/	
{
		
    int ier=0;					/*error code*/
    
    void *ptr1=NULL;                            /*pointer to 1-D array*/
    void **ptr2=NULL;                           /*pointer to 2-D array*/
    void ***ptr3=NULL;                          /*pointer to 3-D array*/
    register int counter=0;                     /*dimension counter*/

    
    /*allocate memory for data matrix */
    
    ier=allocate_1d((void **)&ptr1, (uint16)(dim1*dim2*dim3), elsize);
                                                /* cast by G. Ye, 8/8/96 */
	if(ier==FAILURE){
	   modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "allocate_1d, allocate.c");
	   return(FAILURE); 
	}          

    /*allocate array of void ptrs for dimension 2 */
    
    ier=allocate_1d((void **)&ptr2, (uint16)(dim1*dim2), sizeof(void *));
                                                /* cast by G. Ye, 8/8/96 */
 	if(ier==FAILURE){
	  modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "allocate_1d, allocate.c");
	  return(FAILURE); 
	}           
    
    /*allocate array of void ptrs for dimension 1 */
	
	ier=allocate_1d((void **)&ptr3, dim1, sizeof(void **));
 	if(ier==FAILURE){
	   modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "allocate_1d, allocate.c");
	   return(FAILURE); 
	}     
	    
    /*point allocated dim1 pointers at appropriate data elements*/
    
    for(counter=0;counter<dim1;counter++){      
        ptr3[counter]=ptr2+(counter*dim2);
    }
    
    /*point allocated dim2 pointers at appropriate data elements*/
    
    for(counter=0;counter<(dim1*dim2);counter++){       
        ptr2[counter]=(char *)ptr1+(counter*dim3*elsize);
    }
    *i_ptr=ptr3;
    return(SUCCESS);
}



int free_memory(struct vars_and_ptrs *all)
/******************************************************************************
!C

!Description: 

    This function frees all dynamically allocated memory 
    via repeated calls to free().

!Input Parameters:

    all -- pointer to structure vars_and_ptrs
    
!Output Parameters:

    all, some variable members of this structure may have changed
    
    Integer return status qualifying function termination conditon

!Revision History:

    none

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:

    Written by Jordan S. Borak  

		Boston University
		Department of Geography &
		Center for Remote Sensing
		675 Commonwealth Avenue
		Boston, MA   02215  

		617-353-2088

		borak@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a 
    complete crash of the code.

Externals:

               
!END
*******************************************************************************/
{

    register int counter=0;		    /*counter*/
    char  string[STR_MAX]={};		/*string for SMF logs*/

       
    /*****free LW array****/

    free(all->ptrs.marray.lw[0]);
    free(all->ptrs.marray.lw);

    /***DEM_LW***/
  free(all->ptrs.marray.dem_lw[0]);
    free(all->ptrs.marray.dem_lw);


    sprintf(string, "fixed: End of free_memory function"); 
 //   modsmf(MODIS_S_SUCCESS, string, "main: quarterly.c");
    
    return(SUCCESS);
}    



