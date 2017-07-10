#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define FALSE   0

herr_t file_info(hid_t loc_id, const char *name, void* opdata);
hid_t af_open(char* file_path);
herr_t af_close(hid_t file);
double* af_read(hid_t file, char* dataset_name);
void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size);
int af_misr_handler(char* file_path, char* camera_angle, char* resolution, char* radiance);

//Hard-coded information to test MISR IO operations
char* dataset_name;

int af_misr_handler(char* file_path, char* camera_angle, char* resolution, char* radiance){
	//Path to dataset proccessing 
	char* instrument = "MISR";
	char* d_fields = "Data_Fields";
	char* location;
	if(strcmp(resolution, "H") == 0){
		location = "HRGeolocation";
	}
	else{
		location = "Geolocation";
	}
	char* lat = "GeoLatitude";
	char* longitude = "GeoLongitude";
	const char* arr[] = {instrument, camera_angle, d_fields, radiance};
	const char* arr2[] = {instrument, location, lat};
	const char* arr3[] = {instrument, location, longitude};
	char* rad_dataset_name;
	concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(camera_angle) + strlen(d_fields) + strlen(radiance) + 4, 4);

	//Open file
	hid_t file = af_open(file_path);
	/*Dimensions - 180 blocks, 512 x 2048 ordered in 1D Array*/
	//Retrieve radiance dataset and dataspace
	double* data = af_read(file, rad_dataset_name);
	printf("rad_data: %f\n", data[0]);
	
	//Retrieve geolocation dataset and dataspace
	char* lat_dataset_name;
	concat_by_sep(&lat_dataset_name, arr2, "/", strlen(instrument) + strlen(location) + strlen(lat) + 4, 3);
	double* lat_data = af_read(file, lat_dataset_name);
	char* long_dataset_name;
	concat_by_sep(&long_dataset_name, arr3, "/", strlen(instrument) + strlen(location) + strlen(longitude) + 4, 3);
	double* long_data = af_read(file, long_dataset_name);
    /* Close file */
    free(data);
    herr_t ret = af_close(file);
	
	return ret;
}

int af_modis_handler(char* file_path, char* resolution){
	//Get all granule file names
}

double* af_read(hid_t file, char* dataset_name){
	hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset);
	const int ndims = H5Sget_simple_extent_ndims(dataspace);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	double * data = calloc ( dims[0] * dims[1] * dims[2] , sizeof(double) );
	hid_t memspace = H5Screate_simple(ndims,dims,NULL);
	herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, memspace, H5P_DEFAULT, data);
	printf("read status: %d\n", status);
	return data;
}

hid_t af_open(char* file_path){
	hid_t f = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
	return f;
}

herr_t af_close(hid_t file){
	herr_t ret = H5Fclose(file);
	return ret;
}

int main (int argc, char *argv[]){
	//Preset filename here for easy testing 
	char* file_path = "/projects/TDataFus/kent/temp/40-orbit-file/Jun15.2/TERRA_BF_L1B_O69365_F000_V000.h5";
	argv[1] = file_path;
	if(argc < 3){
		printf("Usage: %s filename instrument_name extra_arguments\n", argv[0]);
		return 0;
	}
	else if(strcmp(argv[2], "MISR") == 0){
		if(argc < 6){
			printf("MISR Usage: %s filename MISR camera_angle resolution(H/L) radiance\n");
			return 0; 
		}
		else{
			//MISR input requirements fulfilled 
			af_misr_handler(argv[1], argv[3], argv[4], argv[5]);
			return 0;
		}
	}

	return 0;
}

//Helper Functions
//String helper
void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size){
	int i;
	*source = calloc(length+1, sizeof(char));
	for(i = 0; i < arr_size; i++){
		printf("%d %s\n", i, w[i]);
		if(i == 0){
			strncpy(*source, sep, strlen(sep));
			strncat(*source, w[i], strlen(w[i]));
		}
		else{
			strncat(*source, sep, strlen(sep));
			strncat(*source, w[i], strlen(w[i]));
		}
	}
}
