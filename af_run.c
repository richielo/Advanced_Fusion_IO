/*
    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <sys/time.h>
#include "reproject.h"
#include "io.h"

#define MAX_MODIS_BANDS 38

int main(int argc, char ** argv) {
	//Input variables
	char* file_path;
	char* output_file;
	char* project_instrument;
	char* base_instrument; 
	char* misr_args[3][50];
	char* modis_args[1][50];
	char* modis_bands[MAX_MODIS_BANDS][50];

	if(argc < 2){
		printf("Usage: ./af_run input_parameters.txt\n");
		return -1;
	}
	
	//Open input parameters file
	FILE* file = fopen(argv[1], "r");
	if(file == NULL){
		printf("Input parameters file does not exist\n");
		return -1;
	}
	//Read in input parameters
	char line[256];
	int line_count = 0;
	while(fgets(line, sizeof(line), file)){
		char* arg = strchr(line,'=');
		arg++;
		char* value = arg;
		printf("Found value: %s\n", value);
	}
	
	//Close input parameters file
	fclose(file);
	return 0;

}
