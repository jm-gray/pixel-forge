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



/* Shouldn't this be the same logical number as 
   my output??? Since they are both MOD12Q??? Maybe 
   just different version numbers???*/




/*ESDT Shortname*/
#define ESDTNAME "MOD12Q_PHE"    

						
/*Info relevant to HDF 
file formats*/



#define NUM_ANC_FILES 3				/*Number of ancillary files*/
											
						
#define DATA_LAYERS 10				/*number of output layers, 
						less CPGGR*/
						
#define NUM_MOD12_SDS 15			/*total number of SDS's in 
						MOD12 HDF product*/
					

#define NUM_CLASS_SDS 2				/*Number SDS's in MOD12_PREV
						used for LC classification*/
						
#define NUM_PHE_CLASSES 18			/*Total number of LC classes*/

#define NUM_UMD_CLASSES 16
#define NUM_LF_CLASSES 10  /**zhang11/29/04**/			/*Total number of  LAI/FPAR  classes*/
#define NUM_BGC_CLASSES 10			/*Total number of  BGC  classes*/

	
#define	NUMMODES 2							
/*Info about monthly input files*/		
				
  
				

/*Water flags*/

#define WATER 0					/*Water flag*/
#define NON_WATER 3				/*Non-water flag*/

/*QC-related constants*/

#define NEW_AGE_QC 192				/*min. QC for new observation*/
#define FILL -1					/*fill pixel flag*/
#define OUT_OF_RANGE -2				/*out of range pixel flag*/
#define GOOD_VALUE 0				/*good data pixel flag*/
#define MISSING_VALUE 101			/*Missing value*/
#define PREVIOUS_FLAG 0				/*previous value for overall 
						

/*M-API constants*/

#define DIM_MAX 10				/*maximum num. of dimensions*/
#define FILL_META_STRING Null
#define FILL_META_VALUE -999


/*Maximum string length*/

#define STR_MAX 80
#define MAX_STRING_LENGTH 2000
#define SHORT_STRING_LENGTH 30


#define MODERATE_ALLOC  100
#define LARGE_ALLOC  1000	/* Used for allocating buffer space*/
#define HUGE_ALLOC  10000       /* Used for allocating buffer space */


/*Structure declarations*/


/*Structure that contains 
all file pointers*/

struct file_pointers{

    MODFILE *ofp;			/*Output file pointer*/
    int     fileflag[5];		/*file there or not*/
    int     numofiles[5];		/*Number of versions per file*/
   
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

    char attribute[5][STR_MAX];			/*regular labels*/
    
    /*long names*/
    char ln[18][STR_MAX];    
    char ln_name[STR_MAX];			/*attribute name*/
    char ln_type[STR_MAX];			/*data type*/


    
    /*units*/

    char unit[18][STR_MAX];				/*default units value*/
    char unit_name[STR_MAX];			/*attribute name*/
    char unit_type[STR_MAX];			/*data type*/

 	
    /*valid_range*/
    
 uint8 vr_lct[2];	
  uint8 fv_lct;	

    char vr_name[STR_MAX];			/*attribute name*/
    char vr_type0[STR_MAX];			/*data type0*/
 

    char vr_type1[STR_MAX];			/*data type1*/ 
    char vr_type2[STR_MAX];			/*data type2*/
    char vr_type3[STR_MAX];			/*data type3*/
    char vr_type4[STR_MAX];			/*data type4*/
    char vr_type5[STR_MAX];			/*data type5*/

    uint8 vr_lw[2];    

    uint16 vr_tcv1[2];	       		/*output TCV_Detail_1 range values*/
    uint16 vr_tcv2[2];			/*output TCV_Detail_2 range values*/
    uint16 vr_tcv3[2];		        /*output TCV_Detail_3 range values*/
    uint16 vr_tcv4[2];       	        /*output TCV_Detail_4 range values*/
    uint16 vr_tcv5[2];		        /*output TCV_Detail_5 range values*/
    uint8 vr_dyqc[2];
    uint16 vr_phe1[2];		       /*output onset_greenness range values*/
    uint16 vr_phe2[2];	               /*output onset_maturity range values*/
    uint16 vr_phe3[2];		       /*output onset_senescence range values*/
    uint16 vr_phe4[2];		       /*output onset_dormancy range values*/
    uint16 vr_pkge[2];		       /*output peak_greenness range values*/
    uint16 vr_vige[2];		       /*output VI_value_onset range values*/
    uint16 vr_vima[2];		       /*output VI_Max range values*/
    uint16 vr_viar[2];		       /*output VI_Area range values*/
    uint16 vr_vtbd[2];		       /*output TBD range values*/   







  			/*input lstqc range values*/
	       /*output tempmet3 range values*/
 
    /*_FillValue*/
    

    char fv_name[STR_MAX];			/*attribute name*/
    char fv_type0[STR_MAX];			/*data type0*/

      char fv_type1[STR_MAX];			/*data type1*/
    char fv_type2[STR_MAX];			/*data type2*/
    char fv_type3[STR_MAX];			/*data type3*/
    char fv_type4[STR_MAX];			/*data type4*/
    char fv_type5[STR_MAX];			/*data type5*/

uint8 fv_lw;
    uint16 fv_tcv1;	       		/*output TCV_Detail_1 range values*/
    uint16 fv_tcv2;			/*output TCV_Detail_2 range values*/
    uint16 fv_tcv3;		        /*output TCV_Detail_3 range values*/
    uint16 fv_tcv4;       	        /*output TCV_Detail_4 range values*/
    uint16 fv_tcv5;		        /*output TCV_Detail_5 range values*/
  uint8 fv_dyqc;
    uint16 fv_phe1;		       /*output onset_greenness range values*/
    uint16 fv_phe2;	               /*output onset_maturity range values*/
    uint16 fv_phe3;		       /*output onset_senescence range values*/
    uint16 fv_phe4;		       /*output onset_dormancy range values*/
    uint16 fv_pkge;		       /*output peak_greenness range values*/
    uint16 fv_vige;		       /*output VI_value_onset range values*/
    uint16 fv_vima;		       /*output VI_Max range values*/
    uint16 fv_viar;		       /*output VI_Area range values*/
    uint16 fv_vtbd;		       /*output TBD range values*/   
 
    
   
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
  

};



/*Structure that contains all 
MOD12 input data arrays*/

struct M_arrays{
    
    /*MOD12*/
    
     uint16 **tcv1_data;	       		/*output TCV_Detail_1 range values*/
    uint16 **tcv2_data;			/*output TCV_Detail_2 range values*/
    uint16 **tcv3_data;		        /*output TCV_Detail_3 range values*/
    uint16 **tcv4_data;       	        /*output TCV_Detail_4 range values*/
    uint16 **tcv5_data;		        /*output TCV_Detail_5 range values*/
    uint8 ***dyqc_data;                   /*Time_series_assement**/
    uint16 ***phe1_data;		       /*output onset_greenness range values*/
    uint16 ***phe2_data;	               /*output onset_maturity range values*/
    uint16 ***phe3_data;		       /*output onset_senescence range values*/
    uint16 ***phe4_data;		       /*output onset_dormancy range values*/
    uint16 ***pkge_data;		       /*output peak_greenness range values*/
    uint16 ***vige_data;		       /*output VI_value_onset range values*/
    uint16 ***vima_data;		       /*output VI_Max range values*/
    uint16 ***viar_data;		       /*output VI_Area range values*/
    uint16 ***vtbd_data;		       /*output TBD range values*/   


   uint8 **type1_data;					/*type1 array*/
       
};




/*Structure that contains all pointers*/

struct all_ptrs{

    struct file_pointers fp;			/*file pointers*/
    struct M_arrays marray;			/*input MOD12 data arrays*/

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


    
    /*MOD 12*/
    
    char arrnm_mod12[NUM_MOD12_SDS][STR_MAX];		/*MOD12 SDS's*/
    


};

/*Structure that contains all 
SDS group names*/

struct grpnms{

    char grpnm[STR_MAX];				/*group name*/

};

/*This structure contains all
SDS data type strings*/

struct datatypes{
   
    /*MOD12*/
    
    char datatype_mod12[NUM_MOD12_SDS][10];		/*MOD12 data types*/
    
};


/*This structure contains 
all SDS ranks*/

struct ranks{

 
    /*MOD12*/
    
    long rank_mod12[NUM_MOD12_SDS];			/*MOD12 ranks*/

};


/*This structure contains all
SDS dimension magnitudes*/

struct dimsizes{

 
    /*MOD12*/
    
    long dimsizes_mod12[NUM_MOD12_SDS][DIM_MAX];	/*MOD12 dims*/
    

};


/*This structure contains all
SDS starting points*/

struct starts{

   /*Monthly data*/
    
  /**?????????????????????????????????????????????
   **I donot understand why I have to keep one of 'start_monthly_'**
   ????????????????????????????????????????***/

     long start_monthly_lw[2];				
 	
    	

    long start_mod12q2_tcv1[2];	       		/*output dytype1TCV_Detail_1 range values*/
    long start_mod12q2_tcv2[2];			/*output dytype1TCV_Detail_2 range values*/
    long start_mod12q2_tcv3[2];		        /*output dytype1TCV_Detail_3 range values*/
    long start_mod12q2_tcv4[2];       	        /*output dytype1TCV_Detail_4 range values*/
    long start_mod12q2_tcv5[2];		        /*output dytype1TCV_Detail_5 range values*/

    long start_mod12q2_phe1[3];		       /*output onset_greenness range values*/
    long start_mod12q2_phe2[3];	               /*output onset_maturity range values*/
    long start_mod12q2_phe3[3];		       /*output onset_senescence range values*/
    long start_mod12q2_phe4[3];		       /*output onset_dormancy range values*/
    long start_mod12q2_pkge[3];		       /*output peak_greenness range values*/
    long start_mod12q2_vige[3];		       /*output VI_value_onset range values*/
    long start_mod12q2_vima[3];		       /*output VI_Max range values*/
    long start_mod12q2_viar[3];		       /*output VI_Area range values*/
    long start_mod12q2_vtbd[3];		       /*output TBD range values*/   



    /*MOD12*/
    
    long start_mod12_type1[3];			      /*cover type1 start*/    
      
    long start_mod12_ncol[1];				/*ncol start*/    
    long start_mod12_icolstart[1];			/*icolstart start*/    
    long start_mod12_ncoltile[1];			/*ncoltile start*/    
    long start_mod12_ipixstart[1];			/*ipixstart start*/


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
  /**    struct tile_geom geom;		**/		/*tile geom. variables*/
    
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

								
/*Function to free all 
allocated memory*/

int free_memory(struct vars_and_ptrs *all);				

					
/*Function to write information from 
output structure to output file*/

int do_output(struct vars_and_ptrs *all);				



/*Function to write appropriate 
arrays to output file*/

int put_arrays(struct vars_and_ptrs *all);				


/*Function to write SDS-level metadata*/

int put_sds_met(struct vars_and_ptrs *all);



/*Function to create output MODIS arrays*/

int create_output_arrays(struct vars_and_ptrs *all );				


/*Function to retrieve output SDS 
metadata from previous product*/

int get_sds_met(struct vars_and_ptrs *all);					


/*Remove full path from file names retrieved from PCF*/

int remove_path(char *full_fname, char *short_fname);

/* added two ext function by Gang Ye, 8/8/96 */

extern int putMODISarlabel(MODFILE *file, char *arrayname, char *groupname,
                    long int dimension, char *label);

int init_sds_met(struct vars_and_ptrs *all);


#endif

/*allocate for 1-D array*/
int allocate_1d(void **i_ptr, uint16 dim1, int elsize);

/*allocate for 2-D array*/
int allocate_2d(void ***i_ptr, uint16 dim1, uint16 dim2, int elsize);

/*allocate for 3-D array*/
int allocate_3d(void ****i_ptr, uint16 dim1, uint16 dim2, uint16 dim3,int elsize);

int get_phe(struct vars_and_ptrs *all, int16 NewYear);

int set_phe_met(struct vars_and_ptrs *all,int TileRow,int TileCol);

















