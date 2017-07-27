/*


    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <sys/time.h>
#include "reproject.h"
#include "io.h"

int main(int argc, char ** argv) {
	hid_t output_file = H5Fcreate("misr_modis_test_repro.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char* file_path = "/projects/TDataFus/kent/temp/40-orbit-file/Jun15.2/TERRA_BF_L1B_O69365_F000_V000.h5";
	hid_t file;
	if(0 > (file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	char* m_250_list[2] = {"1", "2"};
		char* m_500_list[5] = {"3", "4", "5", "6", "7"};
	char* km_1_ref_list[15] = {"8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "26"};
	char* kme_1_list[16] = {"20", "21", "22", "23", "24", "25", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"};
	
	int nCellMODIS;
	int band_index;
	int file_size;
	int h;
	for(h = 0; h < 15; h++){
		char* d_name = get_modis_filename("_1KM", km_1_ref_list[h], &band_index);
		double * MODIS_rad = get_modis_rad_by_band(file, "_1KM", d_name, &band_index, &file_size);
		printf("MODIS_rad: %f\n", MODIS_rad[0]);
		printf("MODIS_rad: %f\n", MODIS_rad[2748620]);
	}
	
	herr_t ret = af_close(file);


	return 0;
}
