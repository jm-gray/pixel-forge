Function: 
	  1. Convert binary Landcover files to HDF format
	  2. Put the necessary metadata info into HDF files, especially the hdf files in TEMP directory.
      
input:    1. IGBP, BVN, EVD, PVA, and AG remap result
          2. Land-water mask
	  

output:   1. HDF files in corresponding directories
          
Usage: 
         The souce code including the followiing files:
	 bin2hdf.c
	 bin2hdf_lc.c
	 bin2hdf_lc.h
	 bin2hdf_lc_temp.c
	 bin2hdf_lc_temp.h
	 bin2hdf_paths.txt
	 FILL_TILE
	 gltiles.txt
	 lc_bin2hdf.pcf.templateA
	 lc_bin2hdf.pcf.templateB
	 make_all_files
	 makefile_bin2hdf
	 makefile_bin2hdf_lc
	 makefile_bin2hdf_lc_temp
	 MOD12Q.mcf
	 smf.c
 	 smfio.h

Steps to run:
	 1. Put the following line in your .cshrc file (if already did in last
	 step, ignore this)
	 
	 if ($OS == 'Linux') then
		set path = ($path /data/modwork20/tanbin/f90_9.1.039/bin)		
		if ($HOST == 'modis3.bu.edu') then

			#GCTP lib
			setenv LIBGCTP /data/modwork20/tanbin/GCTPC/gctpc/source
			setenv SRCGCTP /data/modwork20/tanbin/GCTPC/gctpc/source

		else  
       			#setup HDF4 path
			setenv HDFBIN /data/modwork20/tanbin/hdf-tools/hdf4.2/bin
			setenv HDFINC /data/modwork20/tanbin/hdf-tools/hdf4.2/include
			setenv HDFLIB /data/modwork20/tanbin/hdf-tools/hdf4.2/lib

			#setup HDF5 path
			setenv HDF5BIN /data/modwork20/tanbin/hdf-tools/hdf5-1.6.2/bin
			setenv HDF5INC /data/modwork20/tanbin/hdf-tools/hdf5-1.6.2/include
			setenv HDF5LIB /data/modwork20/tanbin/hdf-tools/hdf5-1.6.2/lib

			#source HDFEOS
				source /data/modwork20/tanbin/hdf-tools/hdfeos/bin/linux/hdfeos_env.csh  
			#source HDFEOS5
			source /data/modwork20/tanbin/hdf-tools/hdfeos5/bin/linux/hdfeos_env.csh  

			#source TOOLKIT
			set my_path = ($path)
			source /data/modwork20/tanbin/hdf-tools/TOOLKIT/bin/linux/pgs-dev-env.csh
			set path = ($my_path $pgs_path)	
			setenv PGS_PC_INFO_FILE $HOME/PCF.relB0

			#GCTP lib
			setenv LIBGCTP /data/modwork20/tanbin/GCTPC/gctpc/source
			setenv SRCGCTP /data/modwork20/tanbin/GCTPC/gctpc/source
			#MAPI lib
			setenv API_INC /data/modwork20/tanbin/MAPI/mapi2.3.4/h
			setenv API_LIB /data/modwork20/tanbin/MAPI/mapi2.3.4/lib

			setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH":"$HDF5LIB":"$HDFLIB":"$HDFEOS_LIB":"$HDFEOS5_LIB":"$API_LIB":"$PGSLIB
		endif
	endif

	 2. Login to a 32-bits Linux server, e.g. modis1, modis2 etc. Exclude MODIS3 because the HDF library does not support this 64-bits server now.
	 3. Run this command "./make_all_files " to compile all C programs
	 (after compiling the program, you can use 64-bits Linux server, such as MODIS3).
	 4. Edit bin2hdf_paths.txt to change all parameters.
	 5. Edit gltiles.txt to include the tiles need to be processed
	 6. Run the command "./bin2hdf.exe" 
	
Author: Bin Tan
Last updated: 2007-3-15
