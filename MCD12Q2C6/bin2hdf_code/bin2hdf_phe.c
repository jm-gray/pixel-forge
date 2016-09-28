/* bin2hdf_phe.c -- Converting Binary file to MODIS HDF file for MOD12Q2: version 1.0*/
/**Bin Tan 1/2/2008***/

/*Header file that contains define
statements and function prototypes*/


#include "bin2hdf_phe.h"


char IN_FILE[100],OUT_FILE[100];

char SDS_NAME[100]; /**the name of the input data***/


int main(int argv,char *argc[])   /* changed to return int status*/

/******************************************************************************
!C

!Description:

    Main. Classifies input data to produce a land cover
    map using various decision tree and/or neural network
    classifiers.  A given pixel is classified via a trained tree or
    net that is appropriate given that pixel's location in the grid.

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
                schaaf@crsa.bu.edu, zhang@crsa.bu.edu
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

  int TileRow,TileCol;
  int16 NewYear;
  /*String for SMF logs*/
  char string[STR_MAX];


  /*Declare structure of vars_and_ptrs*/
  struct vars_and_ptrs all={0};


  // //DEBUG: Josh
  // printf("Here 0\n");

  /*Variables necessary for computing special SDS-level metadata*/
  if(argv!=7) {
    printf("Usage: %s <phe_BIP_TILE> <phe_HDF_TILE> <phe_TYPE> <Tile_Row> <Tile_Col> <Begin Year>\n",argc[0]);
    printf("\nphe_BIP_TILE--Input Binary file name (2400x2400)\n");
    printf("phe_HDF_TILE--Output HDF file name\n");
    printf("phe_TYPE--The name of phenology parameter type. It should be one of the followings\n");
    printf("      ----Onset_Greenness_Increase------------(3d-int16)\n");
    printf("      ----Onset_Greenness_Maximum-------------(3d-int16)\n");
    printf("      ----Onset_Greenness_Decrease-----------(3d-int16)\n");
    printf("      ----Onset_Greenness_Minimum-------------(3d-int16)\n");
    printf("      ----NBAR_EVI_Onset_Greenness_Maximum-------------(3d-int16)\n");
    printf("      ----NBAR_EVI_Onset_Greenness_Minimum---------(3d-int16)\n");
    printf("      ----NBAR_EVI_Area--------------------(3d-uint16)\n");
    printf("      ----Time_Series_Assessment-----(3d-uint8, 1-madatory, 2-assessment)\n");
    printf("Tile_Row--The number of row in the tile\n");
    printf("Tile_Col--The number of column in the tile\n");
    printf("BeginYear--The beginning year of calculated phenology\n");
    exit(1);
  }

  // //DEBUG: Josh
  // printf("Here 1\n");

  strcpy(IN_FILE,argc[1]);
  strcpy(OUT_FILE,argc[2]);
  strcpy(SDS_NAME, argc[3]);
  TileRow=atoi(argc[4]);
  TileCol=atoi(argc[5]);
  NewYear=atoi(argc[6]);

  // //DEBUG: Josh
  // printf("Here 1.0\n");

  /*Open all required input files*/


  /*Initialize the variable structure*/

  ier=init_vars(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR,
	   "init_vars,bin2hdf_phe.c");
    return(FAILURE);
  }
  // //DEBUG: Josh
  // printf("Here 1.1\n");

  /*Get Global ECS metadata from inputfiles */




  /*Read relevant SDS metadata from input files */


  ier=get_sds_met(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "get_sds_met, bin2hdf_phe.c");
    return(FAILURE);
  }
  // //DEBUG: Josh
  // printf("Here 2\n");


 /* reset some metadata */
  /*ATTENTION: If you can get right metadata from your input file you don't
    need to call this function. I can't get other information except binary phe
    image. So I need it to calculate some coordinate info. -- Feng */


 ier=set_phe_met(&all,TileRow,TileCol);
  if (ier == FAILURE) {
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "set_phe_met: bin2hdf_phe.c");
    return(FAILURE);
  }
  //   //DEBUG: Josh
  // printf("Here 3\n");





  /* get binary land cover data */
  ier=get_phe(&all, NewYear);
  if (ier == FAILURE) {
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "get_phe: bin2hdf_phe.c");
    return(FAILURE);
  }
  // //DEBUG: Josh
  // printf("Here 4\n");


  /*Attach SDS-level metadata*/

  ier=put_sds_met(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "put_sds_met, bin2hdf_phe.c");
    return(FAILURE);
  }
  // //DEBUG: Josh
  // printf("Here 5\n");


  /*Attach ECS-level metadata*/



  /*Close all open files*/

  ier=close_files(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "close, bin2hdf_phe.c");
    return(FAILURE);
  }
  // //DEBUG: Josh
  // printf("Here 6\n");

  /*Free the allocated memory*/

  ier=free_memory(&all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "free_memory, bin2hdf_phe.c");
    return(FAILURE);
  }
  // //DEBUG: Josh
  // printf("Here 7\n");

  /*Normal program termination*/

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

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:

    Unsuccessful termination of this function results in a
    complete crash of this code.

Externals:


!END
*******************************************************************************/
{

  int ier=0;					/*int ret. status*/
  int sds_counter=0;
  char  string[STR_MAX];		/*string for SMF logs*/


  //  //DEBUG: Josh
  //  printf("Init 0\n");

  /*Initialize total_data_present and total_data_interp*/

  all->emet.total_data_present=(unsigned)0;
  all->emet.total_data_interp=(unsigned)0;


  /*initialize all output longnames*/
    strcpy(all->smet.ln[0], "TCV_Detail_1");
    strcpy(all->smet.ln[1], "TCV_Detail_2");
    strcpy(all->smet.ln[2], "TCV_Detail_3");
    strcpy(all->smet.ln[3], "TCV_Detail_4");
    strcpy(all->smet.ln[4], "TCV_Detail_5");
    strcpy(all->smet.ln[5], "Time_Series_Assessment");
    strcpy(all->smet.ln[6], "Onset_Greenness_Increase");
    strcpy(all->smet.ln[7], "Onset_Greenness_Maximum");
    strcpy(all->smet.ln[8], "Onset_Greenness_Decrease");
    strcpy(all->smet.ln[9], "Onset_Greenness_Minimum");
    strcpy(all->smet.ln[10], "Peak_Greenness");
    strcpy(all->smet.ln[11], "NBAR_EVI_Onset_Greenness_Minimum");
    strcpy(all->smet.ln[12], "NBAR_EVI_Onset_Greenness_Maximum");
    strcpy(all->smet.ln[13], "NBAR_EVI_Area");
    strcpy(all->smet.ln[14], "VI_TBD");
  // # layer_names <- c("num_cycles", "fill_code", "evi_area_cycle1", "evi_amp_cycle1", "evi_min_cycle1", "frac_filled_gup_cycle1", "frac_filled_gdown_cycle1", "length_gup_cycle1", "length_gdown_cycle1", "ogi_cycle1", "midgup_cycle1", "mat_cycle1", "peak_cycle1", "sen_cycle1", "midgdown_cycle1", "dor_cycle1", "ogi_qual_cycle1", "midgup_qual_cycle1", "mat_qual_cycle1", "peak_qual_cycle1", "sen_qual_cycle1", "midgdown_qual_cycle1", "dor_qual_cycle1", "evi_area_cycle2", "evi_amp_cycle2", "evi_min_cycle2", "frac_filled_gup_cycle2", "frac_filled_gdown_cycle2", "length_gup_cycle2", "length_gdown_cycle2", "ogi_cycle2", "midgup_cycle2", "mat_cycle2", "peak_cycle2", "sen_cycle2", "midgdown_cycle2", "dor_cycle2", "ogi_qual_cycle2", "midgup_qual_cycle2", "mat_qual_cycle2", "peak_qual_cycle2", "sen_qual_cycle2", "midgdown_qual_cycle2", "dor_qual_cycle2")

  // /*initialize all output longnames*/
  //     strcpy(all->smet.ln[0], "TCV_Detail_1");
  //     strcpy(all->smet.ln[1], "TCV_Detail_2");
  //     strcpy(all->smet.ln[2], "TCV_Detail_3");
  //     strcpy(all->smet.ln[3], "TCV_Detail_4");
  //     strcpy(all->smet.ln[4], "TCV_Detail_5");
  //     strcpy(all->smet.ln[5], "Time_Series_Assessment");
  //     strcpy(all->smet.ln[6], "Greenup_Twenty_Percent");
  //     strcpy(all->smet.ln[7], "Greenup_Fifty_Percent");
  //     strcpy(all->smet.ln[8], "Greenup_Eighty_Percent");
  //     strcpy(all->smet.ln[9], "Peak_Greenness");
  //     strcpy(all->smet.ln[10], "Greendown_Eighty_Percent");
  //     strcpy(all->smet.ln[11], "Greendown_Fifty_Percent");
  //     strcpy(all->smet.ln[12], "Greendown_Twenty_Percent");
  //     // strcpy(all->smet.ln[10], "Peak_Greenness");
  //     strcpy(all->smet.ln[11], "NBAR_EVI_Onset_Greenness_Minimum");
  //     strcpy(all->smet.ln[12], "NBAR_EVI_Onset_Greenness_Maximum");
  //     strcpy(all->smet.ln[13], "NBAR_EVI_Area");
  //     strcpy(all->smet.ln[14], "VI_TBD");


  //   //DEBUG: Josh
  //  printf("Init 1\n");




    /*initialize all output units*/

  strcpy(all->smet.unit[0], "vector");
    strcpy(all->smet.unit[1], "vector");
    strcpy(all->smet.unit[2], "vector");
    strcpy(all->smet.unit[3], "vector");
    strcpy(all->smet.unit[4], "vector");
    strcpy(all->smet.unit[5], "concatenated flags");
    strcpy(all->smet.unit[6], "date");
    strcpy(all->smet.unit[7], "date");
    strcpy(all->smet.unit[8], "date");
    strcpy(all->smet.unit[9], "date");
    strcpy(all->smet.unit[10], "date");

    strcpy(all->smet.unit[11], "VI value");
    strcpy(all->smet.unit[12], "VI value");
    strcpy(all->smet.unit[13], "VI area");
    strcpy(all->smet.unit[14], "flags for now");
    //DEBUG: Josh
   printf("Init 2\n");

    all->smet.vr_tcv1[0]=0;
    all->smet.vr_tcv1[1]=32766;
    all->smet.vr_tcv2[0]=0;
    all->smet.vr_tcv2[1]=32766;
    all->smet.vr_tcv3[0]=0;
    all->smet.vr_tcv3[1]=32766;
    all->smet.vr_tcv4[0]=0;
    all->smet.vr_tcv4[1]=32766;
    all->smet.vr_tcv5[0]=0;
    all->smet.vr_tcv5[1]=32766;
    all->smet.vr_dyqc[0]=0;
    all->smet.vr_dyqc[1]=254;

    all->smet.vr_phe1[0]=0;
    all->smet.vr_phe1[1]=32766;   /*-356----356?**/
    all->smet.vr_phe2[0]=0;
    all->smet.vr_phe2[1]=32766;   /*-356----356?**/
    all->smet.vr_phe3[0]=0;
    all->smet.vr_phe3[1]=32766;   /*-356----356?**/
    all->smet.vr_phe4[0]=0;
    all->smet.vr_phe4[1]=32766;   /*-356----356?**/
    all->smet.vr_pkge[0]=0;
    all->smet.vr_pkge[1]=32766;   /*-356----356?**/

    all->smet.vr_vige[0]=0;
    all->smet.vr_vige[1]=32766;
    all->smet.vr_vima[0]=0;
    all->smet.vr_vima[1]=32766;
     all->smet.vr_viar[0]=0;
    all->smet.vr_viar[1]=32766;

    all->smet.vr_vtbd[0]=0;
    all->smet.vr_vtbd[1]=32766;   /*TBD**/
    //DEBUG: Josh
   printf("Init 3\n");


    /*initialize all output fill values*/

    all->smet.fv_tcv1=32767;
    all->smet.fv_tcv2=32767;
    all->smet.fv_tcv3=32767;
    all->smet.fv_tcv4=32767;
    all->smet.fv_tcv5=32767;

    all->smet.fv_dyqc=255;
    all->smet.fv_phe1=32767;
    all->smet.fv_phe2=32767;
    all->smet.fv_phe3=32767;
    all->smet.fv_phe4=32767;

    all->smet.fv_pkge=32767;
    all->smet.fv_vige=32767;
    all->smet.fv_vima=32767;

    all->smet.fv_viar=32767;
    all->smet.fv_vtbd=32767;


    all->smet.vr_lct[0]=0;
    all->smet.vr_lct[1]=254;
    all->smet.fv_lct=255;
    //DEBUG: Josh
   printf("Init 4\n");


    /* output scalefactor, calibrated_nt, add_offset, add_offset_err*/
    /*Since all of our units are only flags at present, none of these
      are provided*/




    /*initialize igbp class names*/

      strcpy(all->arrnm.arrnm_mod12[0], "TCV_Detail_1");
  strcpy(all->arrnm.arrnm_mod12[1], "TCV_Detail_2");
  strcpy(all->arrnm.arrnm_mod12[2], "TCV_Detail_3");
  strcpy(all->arrnm.arrnm_mod12[3], "TCV_Detail_4");
  strcpy(all->arrnm.arrnm_mod12[4], "TCV_Detail_5");
  strcpy(all->arrnm.arrnm_mod12[5], "Time_Series_Assessment");
  strcpy(all->arrnm.arrnm_mod12[6], "Onset_Greenness_Increase");
  strcpy(all->arrnm.arrnm_mod12[7], "Onset_Greenness_Maximum");
  strcpy(all->arrnm.arrnm_mod12[8], "Onset_Greenness_Decrease");
  strcpy(all->arrnm.arrnm_mod12[9], "Onset_Greenness_Minimum");
  strcpy(all->arrnm.arrnm_mod12[10], "Peak_Greenness");
  strcpy(all->arrnm.arrnm_mod12[11], "NBAR_EVI_Onset_Greenness_Minimum");
  strcpy(all->arrnm.arrnm_mod12[12], "NBAR_EVI_Onset_Greenness_Maximum");
  strcpy(all->arrnm.arrnm_mod12[13], "NBAR_EVI_Area");
  strcpy(all->arrnm.arrnm_mod12[14], "VI_TBD");
  //DEBUG: Josh
   printf("Init 5\n");




  /*initialize group name*/


  strcpy(all->grpnm.grpnm, "\0");
  //DEBUG: Josh
   printf("Init 6\n");



  for(sds_counter=0;sds_counter<NUM_MOD12_SDS;sds_counter++){

    /* this following bit was used to getaround not having
       a MOD12Prev the first time ---- once you've got a MOD12Prev then
       this doesn't matter (the return(FAILURES) have to be commented
       out in the MOD12Prev open and the section on Previous MOD12Q
       below needs to be commented out*/


    strcpy(all->dtype.datatype_mod12[sds_counter], "uint8");
    all->rank.rank_mod12[sds_counter]=2;
    all->dims.dimsizes_mod12[sds_counter][0]=2400;
    all->dims.dimsizes_mod12[sds_counter][1]=2400;


  }
  //DEBUG: Josh
   printf("Init 7\n");

    strcpy(all->dtype.datatype_mod12[0], "uint16");
    all->rank.rank_mod12[0]=2;
    strcpy(all->dtype.datatype_mod12[1], "uint16");
    all->rank.rank_mod12[1]=2;
    strcpy(all->dtype.datatype_mod12[2], "uint16");
    all->rank.rank_mod12[2]=2;
    strcpy(all->dtype.datatype_mod12[3], "uint16");
    all->rank.rank_mod12[3]=2;
    strcpy(all->dtype.datatype_mod12[4], "uint16");
    all->rank.rank_mod12[4]=2;

    strcpy(all->dtype.datatype_mod12[5], "uint8");
    all->rank.rank_mod12[5]=3;
     all->dims.dimsizes_mod12[5][2]=NUMMODES;

     //DEBUG: Josh
   printf("Init 8\n");

    strcpy(all->dtype.datatype_mod12[6], "uint16");
    all->rank.rank_mod12[6]=3;
    all->dims.dimsizes_mod12[6][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[7], "uint16");
    all->rank.rank_mod12[7]=3;
   all->dims.dimsizes_mod12[7][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[8], "uint16");
    all->rank.rank_mod12[8]=3;
   all->dims.dimsizes_mod12[8][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[9], "uint16");
    all->rank.rank_mod12[9]=3;
   all->dims.dimsizes_mod12[9][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[10], "uint16");
    all->rank.rank_mod12[10]=3;
   all->dims.dimsizes_mod12[10][2]=NUMMODES;
   //DEBUG: Josh
   printf("Init 9\n");

    strcpy(all->dtype.datatype_mod12[11], "uint16");
    all->rank.rank_mod12[11]=3;
   all->dims.dimsizes_mod12[11][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[12], "uin16");
    all->rank.rank_mod12[12]=3;
   all->dims.dimsizes_mod12[12][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[13], "uint16");
    all->rank.rank_mod12[13]=3;
   all->dims.dimsizes_mod12[13][2]=NUMMODES;
    strcpy(all->dtype.datatype_mod12[14], "uint16");
    all->rank.rank_mod12[14]=3;
   all->dims.dimsizes_mod12[14][2]=NUMMODES;

   //DEBUG: Josh
   printf("Init 10\n");




  sprintf(string, "fixed: End of init_var function");
//  modsmf(MODIS_S_SUCCESS, string, "main: bin2hdf_phe.c");
  //DEBUG: Josh
//   printf("Init 11\n");



  return(SUCCESS);

}







int get_sds_met(struct vars_and_ptrs *all)

{

  int ier=0;					/*int return statuses*/

  char string[STR_MAX];			/*string for SMF logs*/




  /*Initialize SDS-level metadata variables*/

  ier=init_sds_met(all);

  if(ier==FAILURE){

    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, "init_smet, bin2hdf_phe.c");

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
  strcpy(all->smet.attribute[2], "Num_QC_Words"); /*zhang/11/24/1999**/
  strcpy(all->smet.attribute[3], "Num_Modes");

  /*long names*/

  strcpy(all->smet.ln_name, "long_name");
  strcpy(all->smet.ln_type, "char *");



  /*units*/

  strcpy(all->smet.unit_name, "units");
  strcpy(all->smet.unit_type, "char *");


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





  /*number of elements*/

  all->smet.nelements1=1L;
  all->smet.nelements2=2L;
  all->smet.nelementsS=STR_MAX;


  /*array dim numbers (for use with getMODISarlabel())*/

  all->smet.dim1=1L;
  all->smet.dim2=2L;




  return(SUCCESS);

}






/*set some metadata for defined tile */

int set_phe_met(struct vars_and_ptrs *all,int TileRow,int TileCol)
{

  double xUpperLeftGrid = -20015109.354;
  double yUpperLeftGrid = 10007554.677;
  double xSize = 926.62543305/2;
  double ySize = 926.62543305/2;


  /* sepecify the corners of tile */
  all->emet.GD_upleft[0]=xUpperLeftGrid+TileCol*xSize*2400;
  all->emet.GD_upleft[1]=yUpperLeftGrid-TileRow*ySize*2400;
  all->emet.GD_lowright[0]=xUpperLeftGrid+(TileCol+1)*xSize*2400;
  all->emet.GD_lowright[1]=yUpperLeftGrid-(TileRow+1)*ySize*2400;


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
  int i=0;
  PGSt_integer       VERSION = 1;
  char file_name[PGSd_PC_FILE_PATH_MAX];
  char string[STR_MAX];	    /*string for SMF logs*/

  int32 gfid=0;
  int32 gid=0;
  char dimlist[MAX_STRING_LENGTH];
  strcpy(all->emet.GD_gridlist,"MOD12Q2_PHE");


  /*Create output SDS's with HDF-EOS for output database*/



  strcpy(file_name,OUT_FILE);

  /* open */
  gfid = GDopen(file_name,DFACC_CREATE);
  if(gfid==GRID_ERRCODE){
    sprintf (string,"Not successful in retrieving grid file ID/open");
    modsmf(MODIS_F_OPEN_HDF_FILE, string,
	   "create_output_arrays: bin2hdf_phe.c");
  }


  /* create grid */
  gid = GDcreate(gfid,all->emet.GD_gridlist,
		 all->emet.GD_ncols, all->emet.GD_nrows,
		 all->emet.GD_upleft,all->emet.GD_lowright);
  if(gfid==GRID_ERRCODE){
    sprintf (string,"Not successful in getting grid ID/ create");
    modsmf(MODIS_E_FUNCTION_ERROR, string,
	   "create_output_arrays: bin2hdf_phe.c");
  }



	all->emet.GD_projcode=99;
	all->emet.GD_zonecode=-1;
	all->emet.GD_spherecode=-1;
	all->emet.GD_projparm[0]=6371007.181000;
	all->emet.GD_projparm[1]=0;
	all->emet.GD_projparm[2]=0;
	all->emet.GD_projparm[3]=0;
	all->emet.GD_projparm[4]=0;
	all->emet.GD_projparm[5]=0;
	all->emet.GD_projparm[6]=0;
	all->emet.GD_projparm[7]=0;
	all->emet.GD_projparm[8]=86400;
	all->emet.GD_projparm[9]=0;
	all->emet.GD_projparm[10]=1;
	all->emet.GD_projparm[11]=0;
	all->emet.GD_projparm[12]=0;



  /* define grid projection */


  ier = GDdefproj(gid,all->emet.GD_projcode,
		  all->emet.GD_zonecode, all->emet.GD_spherecode,
		  all->emet.GD_projparm);
  if(ier==GRID_ERRCODE){
    sprintf (string,"Not successful in defining grid projection");
    modsmf(MODIS_E_FUNCTION_ERROR, string,
	   "create_output_arrays: bin2hdf_phe.c");
  }


  /* define grid origin */


  all->emet.GD_origincode=0;
  ier = GDdeforigin(gid,all->emet.GD_origincode);
  if(ier==GRID_ERRCODE){
    sprintf (string,"Not successful in defining grid origin");

    modsmf(MODIS_E_FUNCTION_ERROR, string,
	   "create_output_arrays: bin2hdf_phe.c");
  }





/*PHE  SDS*/
	for(i=0;i<5;i++)
     	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
  		{


			/* TCV1 SDS */
  			strcpy(dimlist,all->smet.attribute[0]); /*Note that names YDim, XDim are
					    required*/
  			strcat(dimlist,",");
 			strcat(dimlist,all->smet.attribute[1]);
  			ier = GDdeffield(gid,all->arrnm.arrnm_mod12[i],
			                  dimlist,DFNT_UINT16,HDFE_NOMERGE);
  			if(ier==GRID_ERRCODE)
			{
    				sprintf (string,"Not successful in defining TCV1 SDS");
    				modsmf(MODIS_E_FUNCTION_ERROR, string,"create_output_arrays: quarterly2.c");
  			}
 		}
      }

	/* Time_Series_Assessment SDS *2layers***/
	if(strcmp(SDS_NAME,all->smet.ln[5])==0)
	{

		ier=GDdefdim(gid,all->smet.attribute[3],all->dims.dimsizes_mod12[5][2]);
		if(ier==GRID_ERRCODE){
			sprintf (string,"Not successful in defining nbar_bands grid dimension");
			modsmf(MODIS_E_FUNCTION_ERROR, string,
				"create_output_arrays: phe_bin2hdf.c");
		}
	/* time series*/
		strcpy(dimlist,all->smet.attribute[0]); /*Note that names YDim, XDim are
		      required*/
		strcat(dimlist,",");
		strcat(dimlist,all->smet.attribute[1]);
		strcat(dimlist,",");
		strcat(dimlist,all->smet.attribute[3]);

		ier = GDdeffield(gid,all->arrnm.arrnm_mod12[5],
		dimlist,DFNT_UINT8,HDFE_NOMERGE);
		if(ier==GRID_ERRCODE){
			sprintf (string,"Not successful in defining Dynamic QC SDS");
			modsmf(MODIS_E_FUNCTION_ERROR, string,
				"create_output_arrays: quarterly2.c");
		}
	}


	for(i=6;i<11;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{

			/* define grid dimensions for those SDSs that have more than
			2 dimensions *** change*****/
			ier=GDdefdim(gid,all->smet.attribute[3],all->dims.dimsizes_mod12[i][2]);
			if(ier==GRID_ERRCODE){
				sprintf (string,"Not successful in defining nbar_bands grid dimension");
				modsmf(MODIS_E_FUNCTION_ERROR, string,
					"create_output_arrays: phe_bin2hdf.c");
			}



			/* Onset_Greenness */
			strcpy(dimlist,all->smet.attribute[0]); /*Note that names YDim, XDim are
				      required*/
			strcat(dimlist,",");
			strcat(dimlist,all->smet.attribute[1]);
			strcat(dimlist,",");
			strcat(dimlist,all->smet.attribute[3]);

			ier = GDdeffield(gid,all->arrnm.arrnm_mod12[i],
			dimlist,DFNT_INT16,HDFE_NOMERGE);
			if(ier==GRID_ERRCODE){
				sprintf (string,"Not successful in defining PHE1 SDS");
				modsmf(MODIS_E_FUNCTION_ERROR, string,
					"create_output_arrays: quarterly2.c");
			}
		}
	}



	for(i=11;i<13;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{


			/* define grid dimensions for those SDSs that have more than
			2 dimensions *** change*****/
			ier=GDdefdim(gid,all->smet.attribute[3],all->dims.dimsizes_mod12[i][2]);
			if(ier==GRID_ERRCODE){
				sprintf (string,"Not successful in defining nbar_bands grid dimension");
				modsmf(MODIS_E_FUNCTION_ERROR, string,
					"create_output_arrays: phe_bin2hdf.c");
			}


			/*VI_Greenness */
			strcpy(dimlist,all->smet.attribute[0]); /*Note that names YDim, XDim are
				      required*/
			strcat(dimlist,",");
			strcat(dimlist,all->smet.attribute[1]);
			strcat(dimlist,",");
			strcat(dimlist,all->smet.attribute[3]);

			ier = GDdeffield(gid,all->arrnm.arrnm_mod12[i],
			dimlist,DFNT_INT16,HDFE_NOMERGE);

			if(ier==GRID_ERRCODE){
				sprintf (string,"Not successful in defining VIGE SDS");
				modsmf(MODIS_E_FUNCTION_ERROR, string,
			"create_output_arrays: quarterly2.c");
			}
		}
	}



	for(i=13;i<15;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{
		/*VI_Area*/

		/* define grid dimensions for those SDSs that have more than
		2 dimensions *** change*****/
			ier=GDdefdim(gid,all->smet.attribute[3],all->dims.dimsizes_mod12[i][2]);
			if(ier==GRID_ERRCODE){
				sprintf (string,"Not successful in defining nbar_bands grid dimension");
				modsmf(MODIS_E_FUNCTION_ERROR, string,
					"create_output_arrays: phe_bin2hdf.c");
			}



			strcpy(dimlist,all->smet.attribute[0]); /*Note that names YDim, XDim are
				      required*/
			strcat(dimlist,",");
			strcat(dimlist,all->smet.attribute[1]);
			strcat(dimlist,",");
			strcat(dimlist,all->smet.attribute[3]);

			ier = GDdeffield(gid,all->arrnm.arrnm_mod12[i],
			dimlist,DFNT_UINT16,HDFE_NOMERGE);
			if(ier==GRID_ERRCODE){
				sprintf (string,"Not successful in defining VIAR SDS");
				modsmf(MODIS_E_FUNCTION_ERROR, string,
					"create_output_arrays: quarterly2.c");
			}

		}
	}


	/* detach grid */
	ier = GDdetach(gid);
	if(ier==GRID_ERRCODE){
		sprintf (string,"Unable to detach grid.");
		modsmf(MODIS_E_FUNCTION_ERROR, string,
			"create_output_arrays: lc_bin2hdf.c");
	}


	/* close files for grid access */
	ier = GDclose(gfid);
	if(ier==GRID_ERRCODE){
		sprintf (string,"GD-file close LCQuarterly failed.");
		modsmf(MODIS_E_FUNCTION_ERROR, string,
			"create_output_arrays: lc_bin2hdf.c");
	}



	/*Reopen the output file to sort out pointer to the output file*/


	all->ptrs.fp.ofp=openMODISfile(file_name, "a");
	if(all->ptrs.fp.ofp==NULL){
		sprintf(string, "openMODISfile, %s, lc_bin2hdf.c", file_name);
		modsmf(MODIS_F_OPEN_HDF_FILE, F_ERR, string);
		return(FAILURE);
	}


	sprintf(string, "fixed: End of create_output_arrays function");
//	modsmf(MODIS_S_SUCCESS, string, "main: lc_bin2hdf.c");



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
  char  string[STR_MAX];		/*string for SMF logs*/



  /* Close the output file and
     insert proper header info*/


  completeMODISfile(&all->ptrs.fp.ofp, mdHandles, HDFattrnms, NumHandles);



  sprintf(string, "fixed: End of close function");
//  modsmf(MODIS_S_SUCCESS, string, "main: bin2hdf_phe.c");


  return(SUCCESS);

}




int allocate_arrays(struct vars_and_ptrs *all, uint16 Data_Columns)

/******************************************************************************
!C

!Description:




!Design Notes:

    Unsuccessful termination of this function results in a
    complete crash of this code.

Externals:


!END

*******************************************************************************/

{

  int ier=0;					    /*int ret. status*/

  register int counter=0;			    /*counter*/

  char string[STR_MAX];		    /*string for SMF logs*/



  /*Allocate for MODIS data file arrays;
    on a per-line basis*/

  /***MOD12Q2**/

  ier=allocate_2d((void ***)&all->ptrs.marray.tcv1_data, 1, Data_Columns,
		  sizeof(uint16));
  if(ier==FAILURE){
    sprintf(string, "allocate_2d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[0]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
  }

 ier=allocate_2d((void ***)&all->ptrs.marray.tcv2_data, 1, Data_Columns,
		  sizeof(uint16));
  if(ier==FAILURE){
    sprintf(string, "allocate_2d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[1]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
  }


 ier=allocate_2d((void ***)&all->ptrs.marray.tcv3_data, 1, Data_Columns,
		  sizeof(uint16));
  if(ier==FAILURE){
    sprintf(string, "allocate_2d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[2]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
  }

 ier=allocate_2d((void ***)&all->ptrs.marray.tcv4_data, 1, Data_Columns,
		  sizeof(uint16));
  if(ier==FAILURE){
    sprintf(string, "allocate_2d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[3]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
  }

 ier=allocate_2d((void ***)&all->ptrs.marray.tcv5_data, 1, Data_Columns,
		  sizeof(uint16));
  if(ier==FAILURE){
    sprintf(string, "allocate_2d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[4]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
  }






  ier=allocate_3d((void ****)&all->ptrs.marray.dyqc_data,1,
		  Data_Columns,NUMMODES, sizeof(uint8));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[5]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
  }
  /**dyqc is based on the Time_Series_Assessment format**/



  /*phenology data is stored in a BIP format**/
    ier=allocate_3d((void ****)&all->ptrs.marray.phe1_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[6]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }


    ier=allocate_3d((void ****)&all->ptrs.marray.phe2_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[7]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }

    ier=allocate_3d((void ****)&all->ptrs.marray.phe3_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[8]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }

    ier=allocate_3d((void ****)&all->ptrs.marray.phe4_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[9]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }



    ier=allocate_3d((void ****)&all->ptrs.marray.pkge_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[10]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }

    ier=allocate_3d((void ****)&all->ptrs.marray.vige_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[11]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }

    ier=allocate_3d((void ****)&all->ptrs.marray.vima_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[12]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }


    ier=allocate_3d((void ****)&all->ptrs.marray.viar_data,
		    1, Data_Columns, NUMMODES, sizeof(int16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[13]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }

    ier=allocate_3d((void ****)&all->ptrs.marray.vtbd_data,
		    1, Data_Columns, NUMMODES, sizeof(uint16));
  if(ier==FAILURE){
    sprintf(string, "allocate_3d, %s, quarterly2.c",
	    all->arrnm.arrnm_mod12[14]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
    return(FAILURE);
    }
  /*Allocate for MODIS data file arrays;
    on a per-line basis*/



  /*MOD12*/


  ier=allocate_2d((void ***)&all->ptrs.marray.type1_data, 1, Data_Columns,
		  sizeof(uint8));
  if(ier==FAILURE){
    sprintf(string, "allocate_2d, %s, bin2hdf_phe.c",
	    all->arrnm.arrnm_mod12[0]);
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR, string);
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

    Developers: Crystal Schaaf, Feng Gao, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:


!END
*******************************************************************************/
{

    int ier=0;					/*int return statuses*/
    int i=0;


    char string[STR_MAX];			/*SMF string*/


    /*Write SDS-level metadata*/


	for(i=0;i<15;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{

			/*long_name*/

			/**  all->smet.nelementsS=(long)strleSTR_MAX;**/
			all->smet.nelementsS=strlen(all->smet.ln[i]);
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.ln_name,
			all->smet.ln_type,
			all->smet.nelementsS,
			all->smet.ln[i]);
			if(ier==FAILURE){
				sprintf(string, "putMODISarinfo, %s, phe_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}


			/*units*/


			/**  all->smet.nelementsS=(long)STR_MAX;**/
			all->smet.nelementsS=strlen(all->smet.unit[i]);
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.unit_name,
			all->smet.unit_type,
			all->smet.nelementsS,
			all->smet.unit[i]);
			if(ier==FAILURE){
				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}


		}
	}
	for(i=0;i<5;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{

		/*valid range*/


			all->smet.nelements2=2L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.vr_name,
			all->smet.vr_type2,
			all->smet.nelements2,
			all->smet.vr_tcv1);
			if(ier==FAILURE){
				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}



			/*fill values*/

			all->smet.nelements1=1L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.fv_name,
			all->smet.fv_type2,
			all->smet.nelements1,
			&all->smet.fv_tcv1);

			if(ier==FAILURE){
				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}

		}
	}



	if(strcmp(SDS_NAME,all->smet.ln[5])==0)
	{
		/*valid range*/
		all->smet.nelements2=2L;
		ier=putMODISarinfo(all->ptrs.fp.ofp,
		all->arrnm.arrnm_mod12[5],
		all->grpnm.grpnm,
		all->smet.vr_name,
		all->smet.vr_type1,
		all->smet.nelements2,
		all->smet.vr_lct);
		if(ier==FAILURE){

			sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
				all->arrnm.arrnm_mod12[i]);
			modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
			return(FAILURE);
		}



		/*fill values*/

		all->smet.nelements1=1L;
		ier=putMODISarinfo(all->ptrs.fp.ofp,
		all->arrnm.arrnm_mod12[5],
		all->grpnm.grpnm,
		all->smet.fv_name,
		all->smet.fv_type1,
		all->smet.nelements1,
		&all->smet.fv_lct);
		if(ier==FAILURE){

			sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
				all->arrnm.arrnm_mod12[i]);
			modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
			return(FAILURE);
		}
	}




	for(i=6;i<11;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{

		/*valid range*/


			all->smet.nelements2=2L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.vr_name,
			all->smet.vr_type4,
			all->smet.nelements2,
			all->smet.vr_phe1);

			if(ier==FAILURE){

				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}



			/*fill values*/

			all->smet.nelements1=1L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.fv_name,
			all->smet.fv_type4,
			all->smet.nelements1,
			&all->smet.fv_phe1);

			if(ier==FAILURE){

				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}

		}
	}





	for(i=11;i<13;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{

		/*valid range*/


			all->smet.nelements2=2L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.vr_name,
			all->smet.vr_type4,
			all->smet.nelements2,
			all->smet.vr_vige);

			if(ier==FAILURE){

				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}



			/*fill values*/

			all->smet.nelements1=1L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.fv_name,
			all->smet.fv_type4,
			all->smet.nelements1,
			&all->smet.fv_vige);
			if(ier==FAILURE){

				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}

		}
	}



	for(i=13;i<15;i++)
	{
		if(strcmp(SDS_NAME,all->smet.ln[i])==0)
		{

		/*valid range*/


			all->smet.nelements2=2L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.vr_name,
			all->smet.vr_type2,
			all->smet.nelements2,
			all->smet.vr_viar);
			if(ier==FAILURE){

				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}



			/*fill values*/

			all->smet.nelements1=1L;
			ier=putMODISarinfo(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->smet.fv_name,
			all->smet.fv_type2,
			all->smet.nelements1,
			&all->smet.fv_viar);

			if(ier==FAILURE){

				sprintf(string, "putMODISarinfo, %s, lc_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}

		}
	}



	sprintf(string, "fixed: End of put_sds function");
//	modsmf(MODIS_S_SUCCESS, string, "main: phe_bin2hdf.c");

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





int remove_path(char *full_fname, char *short_fname)
/******************************************************************************
!C

!Description:

    This function truncates the full path from a PCF-retrieved
    physical file name.

!Input Parameters:

    full_name	-- full string
    short_name	-- shortened string

!Output Parameters:

    short_name should have changed


!Revision History:

    Version 2:
    Original version prepared as V1 delivery
    Jordan Borak Boston University

!Team-unique Header:

    This software is developed for NASA under contract NAS5-31369.

!References and Credits:
    Principal Investigators: Alan H. Strahler, Boston University
                         alan@crsa.bu.edu

    Developers: Crystal Schaaf, Jordan Borak(Boston University)
                         schaaf@crsa.bu.edu

!Design Notes:



!Externals:


!END
*******************************************************************************/
{
    register int i=0, j=0;				/*counters*/


    for(i=0,j=0;;i++,j++){



	if(full_fname[i]=='/'){		/*If sub-string is a directory name...*/

	    memset(short_fname, '\0', PGSd_PC_FILE_NAME_MAX);

	    j=-1;			/*reset short_fname counter
					before incrementation in loop*/
	    continue;
	}

	if(full_fname[i]=='\0') break;	/*If it's the end of the string...*/

	short_fname[j]=full_fname[i];	/*Assign full filename character
					to short file name*/
    }

    return(SUCCESS);

}




int do_output(struct vars_and_ptrs *all)
/******************************************************************************

!END
*******************************************************************************/
{

    int ier=0;				/*int ret. status*/

    char string[STR_MAX];	/*string for SMF logs*/


    /*Write arrays to output file*/

    ier=put_arrays(all);
    if(ier==FAILURE){
      sprintf(string, "put_arrays, bin2hdf_phe.c");
      modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
      return(FAILURE);
    }
    return(SUCCESS);
}





int put_arrays(struct vars_and_ptrs *all)
/******************************************************************************
!C


!END
*******************************************************************************/
{

    int ier=0;				/*int ret. status*/
    char string[STR_MAX];	/*string for SMF logs*/
   int i=0;
    /*land cover classification arrays*/


	for(i=0;i<5;i++)
	{
		if(strcmp(SDS_NAME, all->smet.ln[i])==0)
		{
			ier=putMODISarray(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->start.start_mod12q2_tcv1,
			all->dims.dimsizes_mod12[i],
			all->ptrs.marray.tcv1_data[0]);

			if(ier==FAILURE){
				sprintf(string, "putMODISarray, %s, phe_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}


		}
	}



	for(i=6;i<11;i++)
	{
		if(strcmp(SDS_NAME, all->smet.ln[i])==0)
		{

		/** printf("333 phe=%d\n",all->ptrs.marray.phe1_data[0][3][10]); **/

			ier=putMODISarray(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->start.start_mod12q2_phe1,
			all->dims.dimsizes_mod12[i],
			all->ptrs.marray.phe1_data[0][0]);

			if(ier==FAILURE){
				sprintf(string, "putMODISarray, %s, phe_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}

		}
	}





	for(i=13;i<15;i++)
	{

		if(strcmp(SDS_NAME, all->smet.ln[i])==0)
		{

			ier=putMODISarray(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->start.start_mod12q2_viar,
			all->dims.dimsizes_mod12[i],
			all->ptrs.marray.viar_data[0][0]);

			if(ier==FAILURE){
				sprintf(string, "putMODISarray, %s, phe_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}

		}
	}


	for(i=11;i<13;i++)
	{

		if(strcmp(SDS_NAME, all->smet.ln[i])==0)
		{

			ier=putMODISarray(all->ptrs.fp.ofp,
			all->arrnm.arrnm_mod12[i],
			all->grpnm.grpnm,
			all->start.start_mod12q2_vige,
			all->dims.dimsizes_mod12[i],
			all->ptrs.marray.vige_data[0][0]);

			if(ier==FAILURE){
				sprintf(string, "putMODISarray, %s, phe_bin2hdf.c",
					all->arrnm.arrnm_mod12[i]);
				modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
				return(FAILURE);
			}
		}
	}





	if(strcmp(SDS_NAME, all->smet.ln[5])==0)
	{

		ier=putMODISarray(all->ptrs.fp.ofp,
		all->arrnm.arrnm_mod12[5],
		all->grpnm.grpnm,
		all->start.start_mod12_type1,
		all->dims.dimsizes_mod12[5],
		all->ptrs.marray.dyqc_data[0][0]);

		if(ier==FAILURE){
			sprintf(string, "putMODISarray, %s, lc_bin2hdf.c",
				all->arrnm.arrnm_mod12[5]);
			modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
			return(FAILURE);
		}


	}


	return(SUCCESS);


}


/* get classification results from binary file and set all information to fit tile 0412 */
/* Purpose: convert binary classification result to HDF-EOS file which looks like a real output
   from MOD12
   */



int get_phe(struct vars_and_ptrs *all, int16 newyear)
{
  uint16            Data_Rows=2400;		 /*Lines per tile*/
  uint16            Data_Columns=2400;		 /*Pixels per line*/
  uint32            line_counter=0;	         /*loop variance for line */
  char              string[STR_MAX];      /*String for SMF logs*/
  int               ier=0;                       /*integer return status*/
  int i,j,imod;

  uint8 data[2400];
  uint16 datau16[2400];
  uint8 datau8[2400];
  int16 data16[4800];
  uint8 datau83[4800];
  uint16 dataui16[4800];


  PGSt_integer julday = 0;
  PGSt_integer julday1 = 2000;  /**year 2000**/
  PGSt_integer year = 0;
  PGSt_integer month = 0;
  PGSt_integer day = 0;
  uint16 doy = 0;
  int doy1=0;
  /**calculating the days from january 1, year 2000****/
  FILE *in;

  if((in=fopen(IN_FILE,"rb"))==NULL)
  {
    printf("can't open file %s\n",IN_FILE);
    exit(1);
  }

  /*Allocate the arrays to hold the input
    and output data, one line at a time*/

  ier=allocate_arrays(all, (uint16)Data_Columns);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR,
	   "allocate_arrays, bin2hdf_phe.c");
    return(FAILURE);
  }


  /*Create the output MODIS arrays*/
  all->emet.GD_ncols=2400;
  all->emet.GD_nrows=2400;

  for(i=0;i<NUM_MOD12_SDS;i++){
    all->dims.dimsizes_mod12[i][0]=1;
    all->dims.dimsizes_mod12[i][1]=2400;
  }

 for(i=5;i<15;i++){

   all->dims.dimsizes_mod12[i][2]=NUMMODES;
  }

  ier=create_output_arrays(all);
  if(ier==FAILURE){
    modsmf(MODIS_E_FUNCTION_ERROR, F_ERR,
	   "create_output_arrays, bin2hdf_phe.c");
    return(FAILURE);
  }


	for(line_counter=0;line_counter<(uint32)Data_Rows;line_counter++)
	{

		all->start.start_mod12_type1[0]=(long)line_counter;
		all->start.start_mod12q2_tcv1[0]=(long)line_counter;
		all->start.start_mod12q2_phe1[0]=(long)line_counter;
		all->start.start_mod12q2_vige[0]=(long)line_counter;
		all->start.start_mod12q2_viar[0]=(long)line_counter;

		/*read PHE types from binary file */

		if(strcmp(SDS_NAME,all->smet.ln[0])==0||strcmp(SDS_NAME,all->smet.ln[1])==0||strcmp(SDS_NAME,all->smet.ln[2])==0||strcmp(SDS_NAME,all->smet.ln[3])==0||strcmp(SDS_NAME,all->smet.ln[4])==0)
		{
			fread(datau16,sizeof(uint16),Data_Columns,in);
			/****two dims*uint16*/
			for(i=0;i<Data_Columns;i++)
			{

				all->ptrs.marray.tcv1_data[0][i]=datau16[i];

			}
		}


		if(strcmp(SDS_NAME,all->smet.ln[5])==0)   /***2 layers of time series**/
		{

			fread(datau83,sizeof(uint8),2*Data_Columns,in);

			j=0;
			for(i=0;i<Data_Columns;i++)
			{
				for(imod=0;imod<2;imod++)
				{
					all->ptrs.marray.dyqc_data[0][i][imod]=datau83[j];
					j=j+1;
				}
			}
		}


		if(strcmp(SDS_NAME,all->smet.ln[6])==0||strcmp(SDS_NAME,all->smet.ln[7])==0||strcmp(SDS_NAME,all->smet.ln[8])==0||strcmp(SDS_NAME,all->smet.ln[9])==0||strcmp(SDS_NAME,all->smet.ln[10])==0)

		{
		/***calculate days**/

			year=(PGSt_integer)newyear;
			julday=PGS_TD_julday(2000,1,1);   /**first day of year 2000**/
			julday1=PGS_TD_julday(year,1,1); /**first day of concerned year**/
			doy=(uint16)(julday1-julday);    /**days first day of year 2000**/

			fread(data16,sizeof(int16),2*Data_Columns,in);
			/****three dims*int16*/
			j=0;
			for(i=0;i<Data_Columns;i++)
			{

				for(imod=0;imod<2;imod++)
				{

					if((data16[j]==all->smet.fv_phe1)||(data16[j]==0))  /**fill value*or water*/
					{
						all->ptrs.marray.phe1_data[0][i][imod]=data16[j];
					}
					else
					{

						/***if input data for phenology calculation range three years (2000, 2001, 2003),
						the calculated DOY is from -1 to -368 for 2000;
						the calculated DOY is from 1 to  368 for 2001---current year the NewYear;
						the calculated DOY is from 369 to 730   for 2002;
						*****xyz***2005**/

						if((data16[j]<-365)&&(data16[j]>-369))
							data16[j]=-365;
						/*   if((data16[j]<0)&&(data16[j]>-366))
						{
						all->ptrs.marray.phe1_data[0][i][imod]=doy-365+abs(data16[j]);
						}
						*/
						if((data16[j]<369)&&(data16[j]>365))
							data16[j]=365;/*
						if((data16[j]>0)&&(data16[j]<366))
							all->ptrs.marray.phe1_data[0][i][imod]=data16[j]+doy;
						*/
						/*  if(data16[j]>368)   /**368-365=3**/
						/* all->ptrs.marray.phe1_data[0][i][imod]=data16[j]+doy-3;*/
						if(data16[j]>368)   /**368-365=3**/
							data16[j]=data16[j]-3;
						all->ptrs.marray.phe1_data[0][i][imod]=data16[j]+doy-3;

					}

					j=j+1;
				}

			}
		}


		if(strcmp(SDS_NAME,all->smet.ln[11])==0||strcmp(SDS_NAME,all->smet.ln[12])==0)
		{

			/**fread(datau83,1,2*Data_Columns,in);**/
			/****three dims*uint8*/
			fread(dataui16,sizeof(uint16),2*Data_Columns,in); /****three dims*int16*/
			j=0;
			for(i=0;i<Data_Columns;i++)
			{

				for(imod=0;imod<2;imod++)
				{

					all->ptrs.marray.vige_data[0][i][imod]=dataui16[j];

					j=j+1;
				}

			}
		}

		if(strcmp(SDS_NAME,all->smet.ln[13])==0||strcmp(SDS_NAME,all->smet.ln[14])==0)

		{
			fread(dataui16,sizeof(uint16),2*Data_Columns,in);
			/****three dims*int16*/
			j=0;
			for(i=0;i<Data_Columns;i++)
			{
				for(imod=0;imod<2;imod++)
				{

					all->ptrs.marray.viar_data[0][i][imod]=dataui16[j];
					j=j+1;
				}

			}
		}


		ier=do_output(all);
		if(ier==FAILURE){
			sprintf(string, "do_output, phe_bin2hdf.c");
			modsmf(MODIS_E_FUNCTION_ERROR, I_ERR, string);
			continue;
		}

	}

	sprintf(string, "fixed: End of get_phe  function");
//	modsmf(MODIS_S_SUCCESS, string, "main: phe_bin2hdf.c");

	fclose(in);
	return(SUCCESS);

}




/*Dynamic memory allocation functions*/

int allocate_1d(void **i_ptr, uint16 dim1, int elsize)
/******************************************************************************
!C

!Description:



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



!Design Notes:

    Unsuccessful termination of this function results in a
    complete crash of the code.

Externals:


!END
*******************************************************************************/
{

    register int counter=0;		    /*counter*/
    char  string[STR_MAX];		/*string for SMF logs*/



        free(all->ptrs.marray.tcv1_data[0]);
	free(all->ptrs.marray.tcv1_data);
        free(all->ptrs.marray.tcv2_data[0]);
	free(all->ptrs.marray.tcv2_data);
        free(all->ptrs.marray.tcv3_data[0]);
	free(all->ptrs.marray.tcv3_data);
        free(all->ptrs.marray.tcv4_data[0]);
	free(all->ptrs.marray.tcv4_data);
        free(all->ptrs.marray.tcv5_data[0]);
	free(all->ptrs.marray.tcv5_data);

        free(all->ptrs.marray.phe1_data[0][0]);
        free(all->ptrs.marray.phe1_data[0]);
        free(all->ptrs.marray.phe1_data);
        free(all->ptrs.marray.phe2_data[0][0]);
        free(all->ptrs.marray.phe2_data[0]);
        free(all->ptrs.marray.phe2_data);
        free(all->ptrs.marray.phe3_data[0][0]);
        free(all->ptrs.marray.phe3_data[0]);
        free(all->ptrs.marray.phe3_data);
        free(all->ptrs.marray.phe4_data[0][0]);
        free(all->ptrs.marray.phe4_data[0]);
        free(all->ptrs.marray.phe4_data);
        free(all->ptrs.marray.pkge_data[0][0]);
        free(all->ptrs.marray.pkge_data[0]);
        free(all->ptrs.marray.pkge_data);
        free(all->ptrs.marray.vige_data[0][0]);
        free(all->ptrs.marray.vige_data[0]);
        free(all->ptrs.marray.vige_data);
        free(all->ptrs.marray.vima_data[0][0]);
        free(all->ptrs.marray.vima_data[0]);
        free(all->ptrs.marray.vima_data);
        free(all->ptrs.marray.viar_data[0][0]);
        free(all->ptrs.marray.viar_data[0]);
        free(all->ptrs.marray.viar_data);
        free(all->ptrs.marray.vtbd_data[0][0]);
        free(all->ptrs.marray.vtbd_data[0]);
        free(all->ptrs.marray.vtbd_data);


    /*free the MOD12 arrays*/

    free(all->ptrs.marray.type1_data[0]);
    free(all->ptrs.marray.type1_data);

    sprintf(string, "fixed: End of free_memory function");
 //   modsmf(MODIS_S_SUCCESS, string, "main: lc_bin2hdf.c");

    return(SUCCESS);
}
