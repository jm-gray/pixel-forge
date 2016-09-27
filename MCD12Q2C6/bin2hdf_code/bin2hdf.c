/********************************************************************************/
/*This program does the following step:                               		*/
/*   1. Generate Land/Water layer HDF file, which contains important metadata	*/
/*	info									*/
/*   2. Generate SDS layers for each landcover classification scheme		*/
/*										*/
/*input:    1. Previous Land/Water layer HDF files, which are used to extract	*/
/*		metadata info (e.g. projection parameters etc.)			*/
/*	    2. binary Landcover classification results				*/
/*										*/
/*output:   1. HDF format landcover files. Each SDS will be stored in a single 	*/
/*		HDF file.							*/
/*										*/
/* Author: Bin Tan (tanbin@bu.edu)					        */
/* Last updated: Mar 15, 2007                                         		*/
/********************************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/param.h>
#include <sys/types.h>
#include <ctype.h>

#define TOTAL_INPUT_AMOUNT 9  /***the amount of effective input lines in file bin2hdf_paths.txt *****/

int get_line(FILE *file,char* input_string);
int get_line_2(FILE *file,char* input_string);
int get_random_name(char* name);
void generate_HDF(char* type, char* processYear, char* tile_id,char* pre_hdf_directory,char* pre_hdf_name_stem,char* output_directory,char* output_hdf_name_stem,char* binary_directory, char* binary_name_stem,char* lw_full_name);
void prepare_pcf(char* tile_id,char* binary_full_name,char* pre_hdf_directory, char* pre_hdf_name_stem);
int main(int argc, char *argv[])
{       
	int read_count,status;
	FILE *f_paths,*f_tile_list;
	char pre_hdf_directory[300],pre_hdf_name_stem[100], output_directory[300], output_hdf_name_stem[100];
	char binary_directory[300], binary_name_stem[100], lw_directory[300], lw_name_stem[100];
	char input_string[500],command_line[500],tile_id[10];
	char pre_hdf_full_name[400],lw_full_name[1000];
	char process_year[20];
	
	/****open the input info files (paths.txt) and read all path and file name info one by one******/
	if ((f_paths=fopen("bin2hdf_paths.txt","r"))==NULL)
	{
		printf("Cannot open the input file --- bin2hdf_paths.txt! Abort!\n");
		return 0;
	}
	read_count=0;
	while(get_line(f_paths,input_string))
	{
		read_count++;
		if (read_count==1)
			sprintf(pre_hdf_directory,"%s",input_string);
		else if (read_count==2)
			sprintf(pre_hdf_name_stem,"%s",input_string);		
		else if (read_count==3)
			sprintf(output_directory,"%s",input_string);
		else if (read_count==4)
			sprintf(output_hdf_name_stem,"%s",input_string);		
		else if (read_count==5)
			sprintf(binary_directory,"%s",input_string);
		else if (read_count==6)
			sprintf(binary_name_stem,"%s",input_string);		
		else if (read_count==7)
			sprintf(lw_directory,"%s",input_string);
		else if (read_count==8)
			sprintf(lw_name_stem,"%s",input_string);
		else if (read_count==9)
			sprintf(process_year,"%s",input_string);
			
	}
	fclose(f_paths);
		
		
	/******If the input info file (paths.txt) is not in the expecting format, quit the program ******/
	
	if (read_count!=TOTAL_INPUT_AMOUNT)
	{
		printf("The number of inputs %d or the format of the input file is not correct!\n",read_count);
		return 0;	
	}
	status=setenv("PGS_PC_INFO_FILE","phe_bin2hdf.pcf",1);
	if (status!=0)
	{
		printf("Failed to set up the environment variable -- PGS_PC_INFO_FILE\n");
		return 0;
	}
	
	/*****finish reading paths.txt **************/
	/**** prepare output direcotries ****/
	sprintf(command_line,"mkdir %s/Greenup",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/Maturity",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/Senescence",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/Dormancy",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/EVI_Minimum",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/EVI_Maximum",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/EVI_Area",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/Assessment",output_directory);
	system(command_line);
	sprintf(command_line,"mkdir %s/TEMP",output_directory);
	system(command_line);
	
	if ((f_tile_list=fopen("gltiles.txt","r"))==NULL)
	{
		printf("Cannot open the input file --- gltiles.txt! Abort!\n");
		return 0;
	}
	
	/**** begin to generate the LC_LW_TEMP.* one tile by one tile ******/
	while(fscanf(f_tile_list,"%s",tile_id)!=EOF)
	{
		printf("Processing tile %s\n",tile_id);
		sprintf(lw_full_name,"%s/%s.%s.05.bin",lw_directory,lw_name_stem,tile_id);
		printf("\t..Layer: TEMP\n");
		generate_HDF("TEMP",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: Greenup\n");
		generate_HDF("Greenup",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: Maturity\n");
		generate_HDF("Maturity",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: Senescence\n");
		generate_HDF("Senescence",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: Dormancy\n");
		generate_HDF("Dormancy",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: EVI_Minimum\n");
		generate_HDF("EVI_Minimum",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: EVI_Maximum\n");
		generate_HDF("EVI_Maximum",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: EVI_Area\n");
		generate_HDF("EVI_Area",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
		printf("\t..Layer: Assessment\n");
		generate_HDF("Assessment",process_year,tile_id,pre_hdf_directory, pre_hdf_name_stem,output_directory,output_hdf_name_stem,binary_directory,binary_name_stem,lw_full_name);
	}
	fclose(f_tile_list);
	/****** clean up all log and temporary files *****/
	/*** These files could be kept when debugging ****/
	sprintf(command_line,"rm -rf Log*");
	system(command_line);
	sprintf(command_line,"rm -rf PHE_LW_TEMP.hdf.met");
	system(command_line);
	sprintf(command_line,"rm -rf ShmMem");
	system(command_line);
	sprintf(command_line,"rm -rf GetAttr.temp");
	system(command_line);
	/*** finish clean up ****/
	return 1;
}


/****** read a valid line from the ASCII file *********/
int get_line(FILE *file,char* string_input)
{
	int indicator=1;
	char char_input;
	int string_count=0;
	int blank_count=0;
	int initial_sign=0;
	while ((char_input=fgetc(file))!=EOF)
	{
	
		if (initial_sign==0)  /***the first letter of a line****/
		{
			if (char_input=='#')
				indicator=0;
			else
				indicator=1;
			initial_sign=1;
		}
		if (char_input=='\n')  /****the end of a line *****/ 
		{
			if (indicator==1 && string_count>0)  /**** there is valid line ******/
			{
				if (blank_count==string_count)  /***** one line contains only space, ignore it *****/
				{	
					string_count=0;
					blank_count=0;
					initial_sign=0;
					continue;
				}
				string_input[string_count]='\0';
				return 1;
			}
			initial_sign=0;
			continue;
		}
		
		if (indicator==1)
		{
			if (char_input==' ')
				blank_count++;
			string_input[string_count]=char_input;
			string_count++;
		}
	}
	return 0;
}
/****This function is used to generate a random file name for temporary file******/
/****This file must be deleted after using *********/
int get_random_name(char* name)
{
	struct tm *ptm;
	time_t rawtime;
	time (&rawtime);
	ptm = gmtime(&rawtime);
	sprintf(name,"%02d-%02d-%02d-%02d-%02d.temp",ptm->tm_mon,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
	return 1;
}

/**** This function converts binary file to hdf format *********/
/**** This is the primary processing function in this program *****/
void generate_HDF(char* type,char* processYear, char* tile_id,char* pre_hdf_directory,char* pre_hdf_name_stem,char* output_directory,char* output_hdf_name_stem,char* binary_directory, char* binary_name_stem,char* lw_full_name)
{
	char binary_full_name[1000],command[1000],output_full_name[400], phe_layer_name[100];
	int tile_row, tile_column;
	char begin_date[20], end_date[20], outputNameStem[20];
	
	sprintf(begin_date,"%s-01-01",processYear);
	sprintf(end_date,"%s-12-31",processYear);
	
	sscanf(tile_id,"h%02dv%02d",&tile_column,&tile_row);
	
	
	if (strcmp(type,"Greenup")==0)
	{
		sprintf(binary_full_name,"%s/Greenup/%s_Greenup_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"Onset_Greenness_Increase");
		sprintf(outputNameStem,"PheGre");
	}
	else if (strcmp(type,"Maturity")==0)
	{
	 	sprintf(binary_full_name,"%s/Maturity/%s_Maturity_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"Onset_Greenness_Maximum");
		sprintf(outputNameStem,"PheMat");
	}
	else if (strcmp(type,"Senescence")==0)
	{
		sprintf(binary_full_name,"%s/Senescence/%s_Senescence_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"Onset_Greenness_Decrease");
		sprintf(outputNameStem,"PheSen");
	}
	else if (strcmp(type,"Dormancy")==0)
	{
	 	sprintf(binary_full_name,"%s/Dormancy/%s_Dormancy_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"Onset_Greenness_Minimum");
		sprintf(outputNameStem,"PheDor");
	}
	else if (strcmp(type,"EVI_Minimum")==0)
	{
		sprintf(binary_full_name,"%s/EVI_Minimum/%s_EVI_Minimum_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"NBAR_EVI_Onset_Greenness_Minimum");
		sprintf(outputNameStem,"EVIMin");
	}
	else if (strcmp(type,"EVI_Maximum")==0)
	{
		sprintf(binary_full_name,"%s/EVI_Maximum/%s_EVI_Maximum_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"NBAR_EVI_Onset_Greenness_Maximum");
		sprintf(outputNameStem,"EVIMax");
	}
	else if (strcmp(type,"EVI_Area")==0)
	{
		sprintf(binary_full_name,"%s/EVI_Area/%s_EVI_Area_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"NBAR_EVI_Area");
		sprintf(outputNameStem,"EVIare");
	}
	else if (strcmp(type,"Assessment")==0)
	{
		sprintf(binary_full_name,"%s/Fill_View/%s_Fill_View_%s.bin",binary_directory,binary_name_stem,tile_id);
		sprintf(phe_layer_name,"Time_Series_Assessment");
		sprintf(outputNameStem,"PheAss");
	}
	else if (strcmp(type,"TEMP")==0)
	{
		sprintf(binary_full_name,"%s",lw_full_name);
	}
	else
	{
		printf("Unknown Land cover layer. Exit!\n");
		exit(0);
	}
	if (strcmp(type,"TEMP")==0)
	{
		prepare_pcf(tile_id,binary_full_name,pre_hdf_directory,pre_hdf_name_stem);
		sprintf(command,"./bin2hdf_phe_temp.exe %s %s %s",binary_full_name, begin_date,end_date);
		system(command);
		sprintf(command,"mv PHE_LW_TEMP.hdf %s/%s/PHE_LW_TEMP.%s.%s.hdf", output_directory,type, output_hdf_name_stem, tile_id);
		system(command);		
	}
	else
	{
		
		
		sprintf(output_full_name,"%s/%s/%s.%s.%s.hdf",output_directory, type,outputNameStem,output_hdf_name_stem, tile_id);
		sprintf(command,"./bin2hdf_phe.exe %s %s %s %d %d %s",binary_full_name,output_full_name,phe_layer_name,tile_row,tile_column, processYear);
		// DEBUG: Josh
		printf("Executing bin2hdf_phe.exe\n");
		
		system(command);
		
		//DEBUG: Josh
		printf("Done with bin2hdf.exe\n");			
	}
		
}

/****** this function prepare *.pcf file which is requried to generate LW hdf file ******/
void prepare_pcf(char* tile_id,char* lw_full_name, char* pre_hdf_directory, char* pre_hdf_name_stem)
{
	FILE *f_pcf, *f_pcf_templateA, *f_pcf_templateB;
	char pre_hdf_name[300], input_a_line[1000];
	
	sprintf(pre_hdf_name,"PHE_LW_TEMP.%s.%s.hdf",pre_hdf_name_stem, tile_id);
	
	if ((f_pcf=fopen("phe_bin2hdf.pcf","w"))==NULL)
	{
		printf("Cannot open pcf file to write. \n");
		exit(0);
	}
	if ((f_pcf_templateA=fopen("phe_bin2hdf.pcf.templateA","r"))==NULL)
	{
		printf("Cannot open pcf template A file to read. \n");
		exit(0);
	}
	if ((f_pcf_templateB=fopen("phe_bin2hdf.pcf.templateB","r"))==NULL)
	{
		printf("Cannot open pcf template B file to read. \n");
		exit(0);
	}
	
	while(get_line_2(f_pcf_templateA, input_a_line))
		fprintf(f_pcf,"%s\n",input_a_line);
		
	fprintf(f_pcf,"212514|%s|%s/TEMP||UR_%s|%s|1\n",pre_hdf_name, pre_hdf_directory, pre_hdf_name, pre_hdf_name);
	
	while(get_line_2(f_pcf_templateB, input_a_line))
		fprintf(f_pcf,"%s\n",input_a_line);
	
	fclose(f_pcf_templateB);
	fclose(f_pcf_templateA);
	fclose(f_pcf);
}
/****** read a valid line from the ASCII file *********/
int get_line_2(FILE *file,char* string_input)
{
	char char_input;
	int string_count=0;
	while ((char_input=fgetc(file))!=EOF)
	{
	
		if (char_input=='\n')  /****the end of a line *****/ 
		{
			string_input[string_count]='\0';
			return 1;
		}
		
		string_input[string_count]=char_input;
		string_count++;
	}
	return 0;
}
