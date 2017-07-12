#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#define FALSE   0

hid_t af_open(char* file_path);
herr_t af_close(hid_t file);
double* af_read(hid_t file, char* dataset_name);
hsize_t* af_read_size(hid_t file, char* dataset_name);

int af_misr_handler(char* file_path, char* camera_angle, char* resolution, char* radiance);
int af_modis_handler(char* file_path, char* resolution, char* d_name);

int alpha_compare(const void *a, const void *b);
void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size);
double dim_sum(hsize_t* dims, int arr_len);
double float_to_double(float f);

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
	if(file < 0){
		printf("File not found\n");
		return -1;
	}
	/*Dimensions - 180 blocks, 512 x 2048 ordered in 1D Array*/
	//Retrieve radiance dataset and dataspace
	double* data = af_read(file, rad_dataset_name);
	if(data == NULL){
		return -1;
	}
	printf("rad_data: %f\n", data[0]);
	
	//Retrieve geolocation dataset and dataspace
	char* lat_dataset_name;
	concat_by_sep(&lat_dataset_name, arr2, "/", strlen(instrument) + strlen(location) + strlen(lat) + 4, 3);
	double* lat_data = af_read(file, lat_dataset_name);
	if(lat_data == NULL){
		return -1;
	}
	printf("lat_data: %f\n", lat_data[0]);
	
	char* long_dataset_name;
	concat_by_sep(&long_dataset_name, arr3, "/", strlen(instrument) + strlen(location) + strlen(longitude) + 4, 3);
	double* long_data = af_read(file, long_dataset_name);
	if(long_data == NULL){
		return -1;
	}
	printf("long_data: %f\n", long_data[0]);
    /* Close file */
    herr_t ret = af_close(file);
	
	return 0;
}

int af_modis_handler(char* file_path, char* resolution, char* d_name){
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	char* lat = "Latitude";
	char* longitude = "Longitude";
	
	//Open file
	hid_t file = af_open(file_path);
	if(file < 0){
		printf("File not found\n");
		return -1;
	}
	//Get all granule file names
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return -1;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][20];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(20*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 20);
		strcpy(&names[i], name);
		free(name);
	}
	
	int h;
	double* data;
	double* lat_data;
	double* long_data;
	double curr_size;
	double curr_lat_size;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, resolution, d_fields, d_name};
		const char* lat_arr[] = {instrument, name, resolution, location, lat};
		const char* long_arr[] = {instrument, name, resolution, location, longitude};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 5);
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(lat), 5);
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(longitude), 5);
		printf("granule_name: %s\n", name);
		if(read_first < 0){
			data = af_read(file, dataset_name);
			if(data == NULL){
				continue;
			}
			curr_size = dim_sum(af_read_size(file, dataset_name), 3); 
			data = af_read(file, dataset_name);
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			lat_data = af_read(file, lat_dataset_name);
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_data = af_read(file, dataset_name);
			if(adding_data == NULL){
				continue;
			}
			double new_d_size = dim_sum(af_read_size(file, dataset_name), 3);
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			
			//Reallocating arrays of data
			data = realloc(data, sizeof(double)*(curr_size + new_d_size));
			memcpy(&data[(int)curr_size], adding_data, sizeof(double)*new_d_size);
			curr_size += new_d_size;

			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;
			
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;
			
			free(adding_data);
			free(adding_lat);
			free(adding_long);
		}
	}
	
	/*Print statements to verify data's existence	
	printf("test data: %f\n", data[0]);
	printf("test_data (next_page): %f\n", data[2748620]);
	printf("test data (next_granule): %f\n", data[41229300]);
	printf("test data (last_granule): %f\n", data[453725399]);
	printf("test_lat_data (next_granule): %f\n", lat_data[0]);
	printf("test_lat_data (next_granule): %f\n", lat_data[2748620]);
	printf("test_lat_data (next_granule): %f\n", lat_data[5510780]);
	*/
	
	
	herr_t ret = af_close(file);
	
	return 0;
}

hsize_t* af_read_size(hid_t file, char* dataset_name){
	hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
	if(dataset < 0){
		printf("Dataset open error\n");
		return NULL; 
	}
	hid_t dataspace = H5Dget_space(dataset);
	if(dataspace < 0){
		printf("Dataspace open error\n");
		return NULL;	
	}
	const int ndims = H5Sget_simple_extent_ndims(dataspace);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	return dims;
}

double* af_read(hid_t file, char* dataset_name){
	hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
	if(dataset < 0){
		printf("Dataset open error\n");
		return NULL; 
	}
	hid_t dataspace = H5Dget_space(dataset);
	if(dataspace < 0){
		printf("Dataspace open error\n");
		return NULL;	
	}
	
	const int ndims = H5Sget_simple_extent_ndims(dataspace);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	hid_t memspace = H5Screate_simple(ndims,dims,NULL);
	hid_t dtype = H5Dget_type(dataset);
	hid_t ndtype = H5Tget_native_type(dtype, H5T_DIR_DESCEND);
	float * data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(ndtype) );
	double* converted_data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(double) );
	
	herr_t status = H5Dread(dataset, ndtype, memspace, memspace, H5P_DEFAULT, data);
	int i;
	for(i=0;i < dim_sum(dims, sizeof(dims)/sizeof(hsize_t)); i++){
		converted_data[i] = float_to_double(data[i]);
	}
	free(data);	
	close(dtype);
	close(ndtype);
	printf("read status: %d\n", status);
	return converted_data;
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
			printf("MISR Usage: %s filename MISR camera_angle resolution(H/L) radiance\n", argv[0]);
			return 0; 
		}
		else{
			//MISR input requirements fulfilled 
			af_misr_handler(argv[1], argv[3], argv[4], argv[5]);
			return 0;
		}
	}
	else if(strcmp(argv[2], "MODIS") == 0){
		if(argc < 4){
			printf("MODIS Usage: %s filename MODIS resolution(1KM/500m/250m)\n", argv[0]);
		}
		else{
			char* resolution = argv[3];
			char* d_name = "";
			if(strcmp(resolution, "1KM") == 0){
				resolution = "_1KM";
				d_name = "EV_1KM_RefSB";
			}
			else if(strcmp(resolution, "250M") == 0){
				resolution = "_250m";
				d_name = "EV_250_RefSB";
			}
			else if(strcmp(resolution, "500M") == 0){
				resolution = "_500m";
				d_name = "EV_500_RefSB";
			}
			else{
				printf("Wrong resolution, choose from 1KM, 500M or 250M\n");
			}
			af_modis_handler(argv[1], resolution, d_name);
		}
	}

	return 0;
}

//Helper Functions
//String helper
void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size){
	int i;
	*source = calloc(length+20, sizeof(char));
	for(i = 0; i < arr_size; i++){
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
//Summing up dimensions
double dim_sum(hsize_t* dims, int arr_len){
	double sum = 0.0;
	int i;
	for(i = 0; i < arr_len; i++){
		if(i == 0){
			sum = (double)dims[i];
		}
		else{
			sum *= (double)dims[i];
		}
	}
	return sum;
}
//Turning float to double
double float_to_double(float f){
	char buf[50];
	sprintf(buf, "%.7g", f);
	return atof(buf);
}
