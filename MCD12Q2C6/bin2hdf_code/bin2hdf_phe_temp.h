/******************************************************************************
!C-INC

!Description: 

    This file is the main header file for quarterly.c. It includes 
    all structure definitions and definitions of all global variables 
    used in quarterly.c in addition to prototypes for nearly all of 
    the functions called in quarterly.c.

!Input Parameters:

    none

!Output Parameters:

    none

!Revision History:

    none

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.
!References and Credits:

    Beta, V1, written by Jordan S. Borak  
    V2 updated by Crystal L.B. Schaaf, Feng Gao

                Boston University
                Department of Geography &
                Center for Remote Sensing
                675 Commonwealth Avenue
                Boston, MA   02215  

                617-353-8033

                schaaf@crsa.bu.edu

!Design Notes:

    none

!Externals:

               
!END
******************************************************************************/
#ifndef QUART_H
#define QUART_H

#include <stdio.h>		    /*Standard io header*/
#include <stdlib.h>		    /*Standard lib header*/
#include <math.h>		    /*Math header*/
#include <string.h>		    /*String header*/
#include <time.h>		    /*Time header*/
#include <ctype.h>

#include "PGS_MET.h"		    /*PGS MET functionality*/
#include "PGS_IO.h"		    /*PGS IO functionality*/
#include "PGS_MODIS_37121.h"        /*Basic mnemonic definitions*/

#include "smfio.h"		    /*Header for PGS SMF*/
#include "mapi.h"		    /*Main mapi Header file*/

#include "hdf.h"
#include "mfhdf.h"
#include "HdfEosDef.h"

/*Function return statuses (integer)*/

#define SUCCESS 0		    /*successful function execution*/
#define FAILURE 1                  /*failure in function execution,
				    changed -1 to 1 according to MODIS
				    Standards by G. Ye, 8/21/96 */

#define GRID_ERRCODE -1         /* HDF-EOS error code */



/*SMF log strings for LogStatus*/

#define I_ERR "Proceeding to next iteration" /*iteration crash*/
#define F_ERR "Exiting"		             /*full code crash*/



/*File string constants*/
/*PGS_PC Logical numbers for all input and output files*/


#define ANC_DEM	    212514 	    /*ANC_DEM file for getting DEM_LW*/

#define WEIGHTS     212515            /*weights file for neural nets*/

#define FEATURES    212516          /*selected features file for classifiers*/

#define OUTPUT	    212550          /*output file*/

#define OUTPUT_HDF	 212549          /*output HDF temporqal file*/



#define MCF_FILE    212551	    /*MCF*/

#define SITES_PERC    212561        /*Trainning sites file saved percentage*/

#define SITES_CODE    212562        /*Trainning sites file saved code*/

#define SUBSET      212563          /*defined subset elements */


  
/*ESDT Shortname*/
#define ESDTNAME "MOD12Q1_TEMP"    

						
/*Info relevant to HDF 
file formats*/



					
/*Neural net data matrix 
location constants*/
						
#define F_ROW 0					/*feature row index in array*/

#define F_COL 1					/*feature col. index in array*/
						

#define DEM_ROW -4				/*Flag to indicate an input 
						net feature is from the DEM*/
						



/*QC-related constants*/

#define NEW_AGE_QC 192				/*min. QC for new observation*/
#define FILL -1					/*fill pixel flag*/
#define OUT_OF_RANGE -2				/*out of range pixel flag*/
#define GOOD_VALUE 0				/*good data pixel flag*/
#define MISSING_VALUE 101			/*Missing value*/
#define PREVIOUS_FLAG 0				/*previous value for overall 
						QC if current is missing*/ 




/*M-API constants*/

#define DIM_MAX 10				/*maximum num. of dimensions*/
#define FILL_META_STRING Null
#define FILL_META_VALUE -999
#define FILL_CLOUD 0                /***default cloud cover*03/08/01****/ 

/*Maximum string length*/

#define STR_MAX 80
#define MAX_STRING_LENGTH 2000
#define SHORT_STRING_LENGTH 30


#define MODERATE_ALLOC  100
#define LARGE_ALLOC  1000	/* Used for allocating buffer space*/
#define HUGE_ALLOC  10000       /* Used for allocating buffer space */




/*****change beginningDATE and ENDDINGDATE every quarterly running*****/
   /****???????????????????????????????????????****/
/****   
#define FAKEBEGINNINGDATE "2000-07-12" 
#define FAKEENDINGDATE "2001-01-16"  
*****input from command line***/
/****???????????????????????????????????????****/
 /*****change beginningDATE and ENDDINGDATE every quarterly running*****/





/*Structure declarations*/


/*Structure that contains 
all file pointers*/

struct file_pointers{

    MODFILE *dem_lw_ifp;			/*DEM file pointer*/ 
    MODFILE *ofp;			/*Output file pointer*/
    int     fileflag[5];		/*file there or not*/
    int     numofiles[5];		/*Number of versions per file*/

    char names[40][PGSd_PC_FILE_PATH_MAX];   /*Monthly file labels*/
 
};


/*Structure that contains all MODLAND- 
specific global metadata variables*/

struct ms_met_vars{

    /*MODLAND-specific global metadata*/

    char attribute[11][STR_MAX];		/*attribute strings*/
    char dtype[11][STR_MAX];			/*data type strings*/
    long nelements1;				/*number of elements1*/
    long nelements2;				/*number of elements2*/
    long nelementsS;				/*number of elements, string*/
    
    float64 ctrl_merid;				/*Central Meridian*/
    float32 cbad;				/*Char. Bin Ang. Dim.*/
    float32 cbs;				/*Characteristic Bin Size*/
    char csf[STR_MAX];				/*Column Storage Format*/
    uint16 data_columns;			/*Data Columns*/
    uint16 data_rows;				/*Data Rows*/
    uint32 dulc[2];				/*Data Upper Left Coordinates*/
    uint32 ggc;					/*Global Grid Columns*/
    uint32 ggr;					/*Global Grid Rows*/
    int16 gnl;					/*Grid Nesting Level*/
    float64 sr;					/*Sphere Radius*/

};

/*Structure that contains all 
SDS-level metadata variables*/

struct sds_met_vars{

    /*array dimension labels*/

    char attribute[2][STR_MAX];			/*regular labels*/
    
    /*long names*/
    char ln[24][STR_MAX];    
    char ln_name[STR_MAX];			/*attribute name*/
    char ln_type[STR_MAX];			/*data type*/

    char ln_lw[STR_MAX];  /**LAND_WATER MASK****/	
  char ln_dem_lw[STR_MAX]; /***LAND_WATER_MASK_IN_DEM**/


    
    /*units*/

    char unit[24][STR_MAX];				/*default units value*/
    char unit_name[STR_MAX];			/*attribute name*/
    char unit_type[STR_MAX];			/*data type*/

  /******Land/Water mask****/
    char unit_lw[STR_MAX];				/*units value*/
  char unit_dem_lw[STR_MAX]; /**LW in DEM***/



   /*add_offset*/

    float64 add_offset[24];
    char add_offset_name[STR_MAX];		/*attribute name*/
    char add_offset_type[STR_MAX];		/*data type*/

/* no offset for lw and snow, KSTai 1/8/1999 */
/*    float64 add_offset_lw; */			/*add_offset value*/




    /*add_offset_err*/

    float64 add_offset_err[24];
    char add_offset_err_name[STR_MAX];		/*attribute name*/
    char add_offset_err_type[STR_MAX];		/*data type*/

/* no offset_err for lw and snow, KSTai 1/8/1999 */ 
/*    float64 add_offset_err_lw;*/		/*add_offset value*/


    /*calibrated_nt*/
    int32 calibrated_nt[24];			/*attribute name*/
    char calibrated_nt_name[STR_MAX];			/*attribute name*/
    char calibrated_nt_type[STR_MAX];		/*data type*/

/* no calibrated_nt for lw and snow, KSTai 1/8/1999 */ 
/*    int32 calibrated_nt_lw;*/ 		/*calibrated_nt value*/


    /*scale_factor*/

    float64 scale_factor[24];			/*attribute name*/
    char scale_factor_name[STR_MAX];			/*attribute name*/
    char scale_factor_type[STR_MAX];		/*data type*/

/* no scale factor for lw and snow, KSTai 1/8/1999 */
/*    float64 scale_factor_lw; */		/*scale_factor  value*/




    
    /*valid_range*/
    
    char vr_name[STR_MAX];			/*attribute name*/
    char vr_type0[STR_MAX];			/*data type0*/
    char vr_type1[STR_MAX];			/*data type1*/ 
    char vr_type2[STR_MAX];			/*data type2*/
    char vr_type3[STR_MAX];			/*data type3*/
    char vr_type4[STR_MAX];			/*data type4*/
    char vr_type5[STR_MAX];			/*data type5*/

  /***land/water mask****/
     uint8 vr_lw[2];				/*input lw range values*/
     uint8 vr_dem_lw[2]; /***LW in DEM***/



    
    int16 vr_dem[2];				/*dem range values*/
    int16 vr_ftrs[2];				/*features range values*/
    float64 vr_wts[2];				/*weights range values*/
    
    /*_FillValue*/
    

    char fv_name[STR_MAX];			/*attribute name*/
    char fv_type0[STR_MAX];			/*data type0*/
    char fv_type1[STR_MAX];			/*data type1*/
    char fv_type2[STR_MAX];			/*data type2*/
    char fv_type3[STR_MAX];			/*data type3*/
    char fv_type4[STR_MAX];			/*data type4*/
    char fv_type5[STR_MAX];			/*data type5*/


  /***land/water mask****/
      uint8 fv_lw;				/*input lw fill value*/
     uint8 fv_dem_lw;  /****LW in DEm***/


    int16 fv_dem;				/*dem fill value*/
    int16 fv_ftrs;				/*features fill value*/
    float64 fv_wts;				/*weights fill value*/
    
  
    
    long nelements1;				/*number of elements 1*/
    long nelements2;				/*number of elements 2*/
    long nelementsS;				/*number of elements, string*/
    
    
    /*dimension info for putMODISarlabel()*/
    
    long dim1;					/*array dim 1*/
    long dim2;					/*array dim 2*/
    
};

/*Structure that contains ECS metadata variables*/

struct ecs_met_vars{
    
    
    unsigned total_data_present;		/*num. of non-missing obs.*/
    unsigned total_data_interp;			/*num. of interpolated obs.*/
    unsigned total_data_possible;		/*total possible num. obs. 
				       	in output (DATA_LAYERS*pixels)*/

  /*ECS Metadata */
   float64 SizeMBECSDataGranule;
   char ReprocessingPlanned[MODERATE_ALLOC];
   char ReprocessingActual[MODERATE_ALLOC];
   char LocalGranuleID[MODERATE_ALLOC];
   char DayNightFlag[MODERATE_ALLOC];
   char LocalVersionID[MODERATE_ALLOC];
   char PGEVersion[MODERATE_ALLOC];
   char InputPointer[MODERATE_ALLOC][HUGE_ALLOC];
   char RangeBeginningDateIn[MODERATE_ALLOC][MODERATE_ALLOC];
   char RangeBeginningDate[MODERATE_ALLOC];
   char RangeBeginningTimeIn[MODERATE_ALLOC][MODERATE_ALLOC];
   char RangeBeginningTime[MODERATE_ALLOC];
   char RangeEndingDateIn[MODERATE_ALLOC][MODERATE_ALLOC];
   char RangeEndingDate[MODERATE_ALLOC];
   char RangeEndingTimeIn[MODERATE_ALLOC][MODERATE_ALLOC];
   char RangeEndingTime[MODERATE_ALLOC];
   char ExclusionGRingFlag[MODERATE_ALLOC];
   float64 GRingPointLatitude[4];
   float64 GRingPointLongitude[4];
   int32 GRingPointSequenceNo[4];
   char ParameterName[MAX_STRING_LENGTH];
   char AutomaticQualityFlag[MODERATE_ALLOC];
   char AutomaticQualityFlagExplanation[MODERATE_ALLOC];
   int32 QAPercentInterpolatedData;
   int32 QAPercentMissingData;
   int32 QAPercentOutofBoundsData;
   int32 QAPercentCloudCover;    /***03/09/01**zhang***/
   char AdditionalAttributeName1[MODERATE_ALLOC];
   char ParameterValue1[MODERATE_ALLOC];
   char AdditionalAttributeName2[MODERATE_ALLOC];
   char ParameterValue2[MODERATE_ALLOC];
   char AdditionalAttributeName3[MODERATE_ALLOC];
   char ParameterValue3[MODERATE_ALLOC];
   char AdditionalAttributeName4[MODERATE_ALLOC];
   char ParameterValue4[MODERATE_ALLOC];
   char AdditionalAttributeName5[MODERATE_ALLOC];
   char ParameterValue5[MODERATE_ALLOC];
   char AdditionalAttributeName6[MODERATE_ALLOC];
   char ParameterValue6[MODERATE_ALLOC];
   char AdditionalAttributeName7[MODERATE_ALLOC];
   char ParameterValue7[MODERATE_ALLOC];
   char AdditionalAttributeName8[MODERATE_ALLOC];
   char ParameterValue8[MODERATE_ALLOC];
   char AdditionalAttributeName9[MODERATE_ALLOC];
   char ParameterValue9[MODERATE_ALLOC];
   char AdditionalAttributeName10[MODERATE_ALLOC];
   char ParameterValue10[MODERATE_ALLOC];
   char AdditionalAttributeName11[MODERATE_ALLOC];
   char ParameterValue11[MODERATE_ALLOC];
   char AdditionalAttributeName12[MODERATE_ALLOC];
   char ParameterValue12[MODERATE_ALLOC];
   char AdditionalAttributeName13[MODERATE_ALLOC];
   char ParameterValue13[MODERATE_ALLOC];
   char AdditionalAttributeName14[MODERATE_ALLOC];
   char ParameterValue14[MODERATE_ALLOC];
   char AdditionalAttributeName15[MODERATE_ALLOC];
   char ParameterValue15[MODERATE_ALLOC];
   char AdditionalAttributeName16[MODERATE_ALLOC];
   char ParameterValue16[MODERATE_ALLOC];
   char AdditionalAttributeName17[MODERATE_ALLOC];
   char ParameterValue17[MODERATE_ALLOC];
   char AdditionalAttributeName18[MODERATE_ALLOC];
   char ParameterValue18[MODERATE_ALLOC];


   float64 WestBoundingCoordinate;
   float64 NorthBoundingCoordinate;
   float64 EastBoundingCoordinate;
   float64 SouthBoundingCoordinate;
   char GranuleBeginningDateTime[MODERATE_ALLOC];
   char GranuleEndingDateTime[MODERATE_ALLOC];
   char InstrumentMode[MODERATE_ALLOC];
   char AlgorithmPackageAcceptanceDate[MODERATE_ALLOC];
   char AlgorithmPackageMaturityCode[MODERATE_ALLOC];
   char AlgorithmPackageName[MODERATE_ALLOC];
   char AlgorithmPackageVersion[MODERATE_ALLOC];
   char GeoAnyAbnormal[MODERATE_ALLOC];
   float64 GeoEstMaxRMSError;
   char LongNameIn[MODERATE_ALLOC];
   char LongNameOut[MODERATE_ALLOC];
   char SPSOParameters[MODERATE_ALLOC];
   char SPSOParameterAlbedo[MODERATE_ALLOC];
   char SPSOParameterNadRef[MODERATE_ALLOC];
   char ProcessingCenter[MODERATE_ALLOC];
   float64 CharacteristicBinAngularSize;
   float64 CharacteristicBinSize;
   int32 DataColumns;
   int32 DataRows;
   int32 GlobalGridColumns;
   int32 GlobalGridRows;
   int32 MaximumObservations;
   int32 NumberofGranules;
   int32 TotalObservations;
   int32 TotalAdditionalObservations;
   char CoverageCalculationMethod[MODERATE_ALLOC];
   char NadirDataResolution[MODERATE_ALLOC];
   char GD_gridlist[MODERATE_ALLOC];
   int32 GD_projcode;
   int32 GD_zonecode;
   int32 GD_spherecode;
   float64 GD_projparm[LARGE_ALLOC];
   int32 GD_origincode;
   int32 GD_ncols;
   int32 GD_nrows;
   float64 GD_upleft[2];
   float64 GD_lowright[2];
   int BeginTimeIndex;
   int EndTimeIndex;
   char BeginnDOY[MODERATE_ALLOC];
   char BeginnYear[MODERATE_ALLOC];
   double BeginTime;
   double EndTime;

};

/*This structure contains all 
tile geometry variables*/

struct tile_geom{
    
    uint32 ncol_left_grid;	/*# columns left of center in the global grid*/    
    uint32 ipix_end_grid;	/*last (rightmost) ISCCP pixel in tile*/    
    uint32 ncol_left_row;	/*# valid global grid columns 
				left of center in tile row*/    
    uint32 ipix_valid_start;	/*1st (leftmost) valid pixel in a tile row*/    
    uint32 ipix_valid_end;	/*last (rightmost) valid pixel in a tile row*/    
    
};

/*Structure that contains monthly 
input data arrays*/




/*Structure that contains all 
MOD12 input data arrays*/

struct M_arrays{
  /****LW******/  
  uint8 **lw;	   /****OUTPUT LAND/WATER MASk******/
  uint8 **dem_lw;      		  /*land-water flag in DEM*/
    /*MOD12*/
    
 
   int32 *ncol_data;					/*ncol array*/
    int32 *icolstart_data;				/*icolstart array*/
    int32 *ncoltile_data;				/*ncoltile array*/
    int32 *ipixstart_data;				/*ipixstart array*/
    
    
};


/*Structure that contains all 
non-MOD12 input data arrays*/

struct arrays{


    double *scaled_dem_data;		/*scaled DEM data*/ 
    double **scaled_monthly_data;	/*array to hold scaled 
    					monthly input data*/   
    
    /*DEM*/
    
    int16 **dem_data;			/*DEM array*/
};


/*Structure that contains all pointers*/

struct all_ptrs{

    struct file_pointers fp;			/*file pointers*/
    struct M_arrays marray;			/*input MOD12 data arrays*/
    struct arrays array;			/*non-MOD12 data arrays*/

};


/*Structure that contains ISCCP 
grid coordinates for a given cell*/

struct icoords{

    uint32 x;						/*x coordinate*/
    uint32 y;						/*y coordinate*/

};


/*Structure that contains all 
M-API array name strings*/

struct arrnms{



    
  /*water mask***/

    char arrnm_lw[STR_MAX];
    
/*DEM_lw*/    
    char arrnm_dem_lw[STR_MAX];				/*DEM SDS*/

};

/*Structure that contains all 
SDS group names*/

struct grpnms{

    char grpnm[STR_MAX];				/*group name*/

};

/*This structure contains all
SDS data type strings*/

struct datatypes{

  
  /******LW*******/

  char datatype_lw[10];	
    
    
       
    /*DEM_lw*/
    
    char datatype_dem_lw[10];				/*DEM_lw*/
};


/*This structure contains 
all SDS ranks*/

struct ranks{

    
 

  /*******LW*******/    
  long rank_lw;	


    /*DEM_lw*/
    
    long rank_dem_lw;					/*DEM_lw*/
};


/*This structure contains all
SDS dimension magnitudes*/

struct dimsizes{


  /*****LW*********/
    long dimsizes_lw[DIM_MAX];	 
    
    /*DEM_lw*/
    
    long dimsizes_dem_lw[DIM_MAX];				/*DEM dims*/
};



/*This structure contains all
SDS starting points*/

struct starts{

    /*Monthly data*/
    
  /*   long start_monthly_lw[2];*/				/*lw start*/    

    /*MOD12*/
    
    long start_mod12_type1[2];			      /*cover type1 start*/    
    long start_mod12_type2[2];			      /*cover type2 start*/    
    long start_mod12_typeqc[2];			    /*cover type qc start*/ 
    long start_mod12_typeassess[2];		    /*cover type assess start*/
    
    long start_mod12_ncol[1];				/*ncol start*/    
    long start_mod12_icolstart[1];			/*icolstart start*/    
    long start_mod12_ncoltile[1];			/*ncoltile start*/    
    long start_mod12_ipixstart[1];			/*ipixstart start*/
        

  /******LW*****/
 long start_lw[2];		
    
    /*DEM_LW*/
    
    long start_dem_lw[2];					/*DEM*/
};


/*Structure that contains all M-API 
HDF variables plus all pointers*/

struct vars_and_ptrs{

    struct all_ptrs ptrs;				/*all pointers*/
    struct icoords icoord;				/*ISCCP grid coords.*/
    struct arrnms arrnm;				/*array names*/
    struct grpnms grpnm;				/*group names*/
    struct datatypes dtype;				/*datatypes*/
    struct ranks rank;					/*ranks*/
    struct dimsizes dims;				/*dimension sizes*/
    struct starts start;				/*array starts*/
    struct ms_met_vars ms_gmet;				/*MODLAND-spec. meta.*/
    struct sds_met_vars smet;				/*sds-level metadata*/
    struct ecs_met_vars emet;				/*MOD12-spec. ecs met.*/
    struct tile_geom geom;				/*tile geom. variables*/
    
};



/*Structure to hold coordinates of 
  selected features in input data
  matrix*/
  
struct matrix_coords{
    
    int16 row;	/*Row position in matrix*/
    
    int16 col;	/*Column position in matrix*/
    
};





    


/*Function prototypes*/


/*Function from which all calls 
to open files are made, changed open to open_files by G. Ye, 8/8/96*/

int open_files(struct vars_and_ptrs *all);


/*Function from which all calls 
to close files are made, changed close to close_files by G. Ye, 8/8/96*/

int close_files(struct vars_and_ptrs *all);


/*Function that initializes most 
variables in the structure "all"*/

int init_vars(struct vars_and_ptrs *all);


/*Function that initializes SDS 
ranks in the structure "all"*/

int init_ranks(struct vars_and_ptrs *all);


/*Function to allocate input 
and scaled input data arrays*/

int allocate_arrays(struct vars_and_ptrs *all, uint16 Data_Columns);


/*Function to set "lines" dimension for each SDS to 1
in order to pick up 1 line with full dimensionality 
per call to getMODISarray*/

int set_input_lines(struct vars_and_ptrs *all);


/*Function to retrieve ECS metadata*/

int get_ecs_met(struct vars_and_ptrs *all);



/*Function to retrieve ECS metadata*/


int process_ecs_met(struct vars_and_ptrs *all, char *FAKEENDINGDATE);

/*Function to retrieve array dimension info from
HDF input files; makes calls to getMODISardims*/

int get_input_ardims(struct vars_and_ptrs *all);

					
/*Function to retrieve MODIS arrays (except; 
for weights); makes calls to getMODISarray*/

int get_input_arrays(struct vars_and_ptrs *all);

					
/*Function to retrieve neural net weights*/

int get_weights(struct vars_and_ptrs *all, int *net_ok_flag);

					
/*Function to retrieve selected feature LUT's*/

int get_features(struct vars_and_ptrs *all, int *net_ok_flag);


/*Function to update starting points of MODIS 
arrays from one iteration to the next, changed to uint32 by G. Ye, 8/8/96*/

int update_input_starts(struct vars_and_ptrs *all, uint32 line_counter);

					
/*Function to free all 
allocated memory*/

int free_memory(struct vars_and_ptrs *all);				


/*Function to scale input data to interval 
[0,1] for use by neural network classifiers, changed to uint32 by G. Ye, 8/8/96*/

extern int do_scaling(struct vars_and_ptrs *all, uint32 cell_counter);


/*Function to insert fill value for all 
output fields at invalid data locations, changed to uint32 by G. Ye, 8/8/96*/

int fill_all(struct vars_and_ptrs *all, uint32 cell_counter);


/*Function to assign appropriate values 
to output SDS if a given pixel is water, changed to uint32 by G. Ye, 8/8/96*/

int do_water_init(struct vars_and_ptrs *all, uint32 cell_counter);


/*Function to check input data for water 
flags; initializes output value in function, changed to uint32 by G. Ye, 8/8/96*/

int water_check(struct vars_and_ptrs *all, uint32 cell_counter);


/*Function to set previous land cover type value, changed to uint32 by G. Ye, 8/8/96*/

int set_previous_lc(struct vars_and_ptrs *all, uint32 cell_counter, 
		uint8 *prev_type_label);

				
/*Function to calculate  
96-day output QC's, changed to uint32 by G. Ye, 8/8/96*/

int do_QC(struct vars_and_ptrs *all, uint32 cell_counter);


/*Function to insert values where
no new MOD12 is generated, changed to uint32 by G. Ye, 8/8/96*/

int if_no_gen(struct vars_and_ptrs *all, uint32 cell_counter);


/*Function to retrieve ISCCP 
grid cell coordinates, changed to uint32 by G. Ye, 8/8/96*/

int get_grid_cell_coords(struct vars_and_ptrs *all, uint32 cell_counter, 
                          uint32 line_counter);

					
/*Function to write information from 
output structure to output file*/

int do_output(struct vars_and_ptrs *all);				


/*Calculate special SDS-specific metadata values*/

int do_special_smet(struct vars_and_ptrs *all, uint32 ncol_tile_total, 
		    uint32 spec_pproc_lc_total, uint32 spec_pchange_lc_total, 
		    char *spec_prop_tmp);


/*Increment land cover class counter 
(for proportion calculations), changed to uint32 by G. Ye, 8/8/96*/

int inc_class(struct vars_and_ptrs *all, uint32 cell_counter);


/*Check for and record change in land cover type 
label (for special SDS-level metadata value), changed to uint32 by G. Ye, 8/8/96*/

int change_check(struct vars_and_ptrs *all, uint32 *spec_pchange_lc_total, 
	         uint8 prev_type_label, uint32 cell_counter);

					
/*Write MODLAND-specific global metadata*/

int put_ms_gmet(struct vars_and_ptrs *all);


/*Function to write appropriate 
arrays to output file*/

int put_arrays(struct vars_and_ptrs *all);				


/*Function to write SDS-level metadata*/

int put_sds_met(struct vars_and_ptrs *all);


/*Function to write SDS-level metadata*/

int write_labels(struct vars_and_ptrs *all);


/*Check input for proper range and fill value*/

int check_input(double val, double fv, double low_lim, double up_lim, 
		int *flag);


/*Function to create output MODIS arrays*/

int create_output_arrays(struct vars_and_ptrs *all );				

/**create output data**/
int create_bip(struct vars_and_ptrs *all);

/*Function to retrieve output SDS 
metadata from previous product*/

int get_sds_met(struct vars_and_ptrs *all);					

/*Function to retrieve all relevant grid attributes
from the input files*/

int get_grid_met(struct vars_and_ptrs *all);

/*Function to retrieve MODLAND-specific 
global metadata from previous product*/

int get_ms_gmet(struct vars_and_ptrs *all);


/*Function to initialize MODLAND-specific 
global metadata variables*/

int init_ms_gmet(struct vars_and_ptrs *all);					


/*Do ECS-level metadata*/


int put_ecs_met(struct vars_and_ptrs *all, char *FAKEBEGINNINGDATE, char *FAKEENDINGDATE);


/*Remove full path from file names retrieved from PCF*/

int remove_path(char *full_fname, char *short_fname);

/* added two ext function by Gang Ye, 8/8/96 */

extern int putMODISarlabel(MODFILE *file, char *arrayname, char *groupname,
                    long int dimension, char *label);

int init_sds_met(struct vars_and_ptrs *all);


#endif

/*get decision tree and begin classification */
int dt_classification(struct vars_and_ptrs *all);

/*get the relationship of months and PCF sequential files*/
int get_month_order(struct vars_and_ptrs *all,char *file_date[]);

/*allocate for 1-D array*/
int allocate_1d(void **i_ptr, uint16 dim1, int elsize);

/*allocate for 2-D array*/
int allocate_2d(void ***i_ptr, uint16 dim1, uint16 dim2, int elsize);

/*allocate for 3-D array*/
int allocate_3d(void ****i_ptr, uint16 dim1, uint16 dim2, uint16 dim3,int elsize);


/**get LW from BAR_QC***/

int get_lw(struct vars_and_ptrs *all, uint32 line_counter);

/***check Nbar_qc****/

int check_input(double val, double fv, double low_lim, double up_lim, int *flag);
/***get bit pattern****/
int get_bit_pattern(unsigned num, int *bitpattern);
