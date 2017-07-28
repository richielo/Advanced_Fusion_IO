/*


    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu


*/

#include "io.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#define FALSE   0

//Band constants for MODIS
char* m_250_list[2] = {"1", "2"};
char* m_500_list[5] = {"3", "4", "5", "6", "7"};
char* km_1_ref_list[15] = {"8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "26"};
char* kme_1_list[16] = {"20", "21", "22", "23", "24", "25", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"};


double* get_misr_rad(hid_t file, char* camera_angle, char* resolution, char* radiance, int* size){
	//Path to dataset proccessing 
	int down_sampling = 0;
	char* instrument = "MISR";
	char* d_fields = "Data_Fields";
	const char* arr[] = {instrument, camera_angle, d_fields, radiance};
	
	//Dataset name parsing
	char* rad_dataset_name;
	concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(camera_angle) + strlen(d_fields) + strlen(radiance) + 4, 4);
		//Check for correct specification
	if(strcmp(camera_angle, "AN") != 0 && strcmp(radiance, "Red_Radiance") != 0 && strcmp(resolution, "H") == 0){
		printf("Error: Your specification does not support high resolution.\n");
		return NULL;
	}
	else if((strcmp(camera_angle, "AN") == 0 || strcmp(radiance, "Red_Radiance") == 0) && strcmp(resolution, "L") == 0){
		//Downsampling has to be done
		down_sampling = 1;
	}
	
	printf("Reading MISR\n");
	/*Dimensions - 180 blocks, 512 x 2048 ordered in 1D Array*/
	//Retrieve radiance dataset and dataspace
	double* data = af_read(file, rad_dataset_name);
	*size = dim_sum(af_read_size(file, rad_dataset_name), 3);
	
	if(data == NULL){
		return NULL;
	}
	printf("Reading successful\n");
	//Variable containing down sampled data
	double* down_data;
	if(down_sampling == 1){
		printf("Undergoing downsampling\n");
		hsize_t* dims = af_read_size(file, rad_dataset_name);
		*size = dims[0] * (dims[1]/4) * (dims[2]/4);
		down_data = malloc(dims[0] * (dims[1]/4) * (dims[2]/4) * sizeof(double));
		int i, j, k;
		for(i = 0; i < dims[0]; i++){
			for(j = 0; j < dims[1]; j = j + 4){
				for(k = 0; k < dims[2]; k = k + 4){
					//Retrieving 4x4 window for averaging
					//Formula for converting i, j and k to index in data array
					//int index = i*dims[1]*dims[2] + j*dims[2] + k;
					int a,b;
					int max_x = j + 4;
					int max_z = k + 4; 
					int* index_array = malloc(16*sizeof(int));
					int index_iter = 0;
					for(a = j; a < max_x; a++){
						for(b = k; b < max_z; b++){
							index_array[index_iter] = i*dims[1]*dims[2] + a*dims[2] + b;
							index_iter += 1;
						}
					}
					double* window = malloc(16*sizeof(double));
					int c;
					for(c = 0; c < 16; c++){
						window[c] = data[index_array[c]];
					}
					//Window Retrieved, get average and assign to new data grid
					double average = misr_averaging(window);
					int new_index = i*dims[1]/4*dims[2]/4 + (j/4)*dims[2]/4 + k/4;
					down_data[new_index] = average;
					free(index_array);
					free(window);
				}
			}
		}
		free(data);
		printf("Downsampling done\n");
	}
	
	if(down_sampling == 1){
		printf("rad_data: %f\n", down_data[0]);
		return down_data;
	}
	else{
		printf("rad_data: %f\n", data[0]);
		return data;	
	}
	
}

double* get_misr_lat(hid_t file, char* resolution, int* size){
	//Path to dataset proccessing 
	char* instrument = "MISR";
	char* location;
	if(strcmp(resolution, "H") == 0){
		location = "HRGeolocation";
	}
	else{
		location = "Geolocation";
	}
	char* lat = "GeoLatitude";
	const char* arr2[] = {instrument, location, lat};
	
	//Dataset names parsing
	char* lat_dataset_name;
	concat_by_sep(&lat_dataset_name, arr2, "/", strlen(instrument) + strlen(location) + strlen(lat) + 4, 3);
	
	printf("Retrieveing latitude data for MISR\n");
	//Retrieve latitude dataset and dataspace
	double* lat_data = af_read(file, lat_dataset_name);
	*size = dim_sum(af_read_size(file, lat_dataset_name), 3);
	if(lat_data == NULL){
		return NULL;
	}
	printf("lat_data: %f\n", lat_data[0]);
	return lat_data;
}

double* get_misr_long(hid_t file, char* resolution, int* size){
	//Path to dataset proccessing 
	char* instrument = "MISR";
	char* location;
	if(strcmp(resolution, "H") == 0){
		location = "HRGeolocation";
	}
	else{
		location = "Geolocation";
	}
	char* longitude = "GeoLongitude";
	const char* arr3[] = {instrument, location, longitude};
	
	//Dataset names parsing
	char* long_dataset_name;
	concat_by_sep(&long_dataset_name, arr3, "/", strlen(instrument) + strlen(location) + strlen(longitude) + 4, 3);
	
	printf("Retrieveing longitude data for MISR\n");
	//Retrieve longitude dataset and dataspace
	double* long_data = af_read(file, long_dataset_name);
	*size = dim_sum(af_read_size(file, long_dataset_name), 3);
	if(long_data == NULL){
		return NULL;
	}
	printf("long_data: %f\n", long_data[0]);
	return long_data;
}

//geo - 0:not geolocation attributes, 1:lat, 2:long
void* get_misr_attr(hid_t file, char* camera_angle, char* resolution, char* radiance, char* attr_name, int geo, void* attr_pt){
	//Path variables
	char* instrument = "MISR";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";

	//Dataset name parsing
	char* rad_dataset_name;
	if(geo == 0){
		const char* arr[4] = {instrument, camera_angle, d_fields, radiance};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(camera_angle) + strlen(d_fields) + strlen(radiance) + 4, 4);
	}
	else if(geo == 1){
		char* lat = "GeoLatitude";
		const char* arr[3] = {instrument, location, lat};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(location) + strlen(lat) + 4, 3);
	}
	else if(geo == 2){
		char* longitude = "GeoLongitude";
		const char* arr[3] = {instrument, location, longitude};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(location) + strlen(longitude) + 4, 3);
	}
	else{
		printf("Wrong geo number");
		return NULL;
	}
	
	//Get attribute
	hid_t attr = H5Aopen_by_name(file, rad_dataset_name, attr_name, H5P_DEFAULT, H5P_DEFAULT);
	if(attr < 0){
		printf("Attribute %s does not exists\n", attr_name);
	}
	hid_t attr_type = H5Aget_type(attr);
	if(strcmp(attr_name, "Units") == 0 || strcmp(attr_name, "units") == 0){
		attr_pt = malloc(sizeof(char) * 50);
		H5Aread(attr, attr_type, attr_pt);
	}
	else if(strcmp(attr_name, "_FillValue") == 0){
		attr_pt = malloc(sizeof(float));
		H5Aread(attr, attr_type, attr_pt);
	}
	H5Aclose(attr);
	
	return attr_pt;
}

double* get_modis_rad(hid_t file, char* resolution, char* bands[], int band_size, int* size){
	printf("Reading MODIS rad\n");
	
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	
	//Get dataset names from bands
	printf("Retreving dataset names\n");
	char* dnames[band_size][50];
	int band_indices[band_size];
	int j;
	for(j = 0; j < band_size; j++){
		char* dname = get_modis_filename(resolution, bands[j], &band_indices[j]);
		if(dname == NULL){
			printf("Band %s is not supported for %s resolution\n", bands[j], resolution);
			return NULL;
		}
		printf("dname: %s\n", dname);
		strcpy(&dnames[j], dname);
	}
	
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int m;
	int total_size = 0;
	for(m = 0; m < band_size; m++){
		for(k = 0; k < num_groups; k++){
			char* name = names[k];
			const char* d_arr[] = {instrument, name, resolution, d_fields, dnames[m]};
			char* dataset_name;
			concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(dnames[m]), 5);
			hsize_t* curr_dim = af_read_size(file, dataset_name);
			if(curr_dim == NULL){
				continue;
			}
			total_size += curr_dim[1]*curr_dim[2];
		}
	}
	
	double* result_data = calloc(total_size, sizeof(double));
	int start_point = 0;
	
	//Start reading data
	int n;
	for(n = 0; n < band_size; n++){
		int file_size;
		double * MODIS_rad = get_modis_rad_by_band(file, resolution, &dnames[n], &band_indices[n], &file_size);
		memcpy(&result_data[start_point], MODIS_rad, file_size);
		start_point += file_size;
	}
	
	if(total_size == start_point){
		printf("Final size validated\n");
	}
	
	*size = total_size;
	
	return result_data;
}

double* get_modis_rad_by_band(hid_t file, char* resolution, char* d_name, int* band_index, int* size){
	printf("Reading MODIS rad by band\n");
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int total_size = 0;
	for(k = 0; k < num_groups; k++){
		char* name = names[k];
		const char* d_arr[] = {instrument, name, resolution, d_fields, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 5);
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		if(curr_dim == NULL){
			continue;
		}
		total_size += curr_dim[1]*curr_dim[2];
	}
	
	//Allocate data size
	double* result_data = calloc(total_size, sizeof(double));
	
	//Retreving data
	int h;
	int curr_size = 0;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		double* data;
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, resolution, d_fields, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 5);
		printf("granule_name: %s\n", name);
		data = af_read(file, dataset_name);
		if(data == NULL){
			continue;
		} 
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		int band_length = curr_dim[1] * curr_dim[2];
		int read_offset = (*band_index)*band_length;
		memcpy(&result_data[curr_size], &data[read_offset], band_length*sizeof(double));
		curr_size += band_length;
		free(data);
	}
	*size = curr_size;
	
	assert(curr_size == total_size);
	printf("Size validated\n");
	
	return result_data;
}

double* get_modis_lat(hid_t file, char* resolution, char* d_name, int* size){
	printf("Reading MODIS lat\n");
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	char* lat = "Latitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, resolution, d_fields, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 4);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* lat_arr[] = {instrument, name, resolution, location, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(lat), 5);
		
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;
			
			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	
	if(lat_data != NULL){
	printf("test_lat_data: %f\n", lat_data[0]);
	printf("test_lat_data: %f\n", lat_data[2748620]);
	printf("test_lat_data: %f\n", lat_data[5510780]);
	}
	
	return lat_data;
}

double* get_modis_long(hid_t file, char* resolution, char* d_name, int* size){
	printf("Reading MODIS long\n");
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	
	int h;
	int valid_granule_count = 0;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, resolution, d_fields, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 4);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		valid_granule_count += 1;
		const char* long_arr[] = {instrument, name, resolution, location, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(longitude), 5);
		
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	
	if(long_data != NULL){
		printf("test_long_data: %f\n", long_data[0]);
		printf("test_long_data: %f\n", long_data[1]);
		printf("test_long_data: %f\n", long_data[1353]);
		printf("test_long_data: %f\n", long_data[1354]);
		printf("test_long_data: %f\n", long_data[2748620]);
		printf("test_long_data: %f\n", long_data[5510780]);
	}
	
	return long_data;
}

//geo - 0:not geolocation attributes, 1:lat, 2:long
double* get_modis_attr(hid_t file, char* resolution, char* d_name, char* attr_name, int geo, void* attr_pt){
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	
	//Get one group name, assuming all attributes across granules are the same
	printf("Retrieving granule group name\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* rad_dataset_name;
	char* name = malloc(50*sizeof(char));
	int h;
	for(h = 0; h < num_groups; h++){
		H5Gget_objname_by_idx(group, (hsize_t)h, name, 50);
		const char* arr[] = {instrument, name, resolution, d_fields, d_name};
		//Dataset name parsing
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name) + 5, 5);
		memmove(&rad_dataset_name[0], &rad_dataset_name[1], strlen(rad_dataset_name));
		htri_t status = H5Lexists(group, rad_dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		else{
			break;
		}
		
	}
	
	if(geo == 1){
		char* lat = "Latitude";
		const char* arr[] = {instrument, name, resolution, location, lat};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(lat) + 5, 5);
	}
	else if(geo == 2){
		char* longitude = "Longitude";
		const char* arr[] = {instrument, name, resolution, location, longitude};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(longitude) + 5, 5);
	}
	
	//Get attribute 	
	hid_t attr = H5Aopen_by_name(file, rad_dataset_name, attr_name, H5P_DEFAULT, H5P_DEFAULT);
	if(attr < 0){
		printf("Attribute %s does not exists\n", attr_name);
	}
	hid_t attr_type = H5Aget_type(attr);
	if(strcmp(attr_name, "units") == 0){
		attr_pt = malloc(sizeof(char) * 50);
		H5Aread(attr, attr_type, attr_pt);
	}
	else if(strcmp(attr_name, "_FillValue") == 0){
		attr_pt = malloc(sizeof(float));
		H5Aread(attr, attr_type, attr_pt);
	}
	else if(strcmp(attr_name, "valid_min") == 0){
		attr_pt = malloc(sizeof(float));
		H5Aread(attr, attr_type, attr_pt);
	}
	H5Aclose(attr);
	
	return attr_pt;
	
}

char* get_modis_filename(char* resolution, char* band, int* band_index){
	if(strcmp(resolution, "_1KM") == 0){
		int i;
		for(i = 0; i < 15; i++){
			if(strcmp(band, km_1_ref_list[i]) == 0){
				*band_index = i;
				return "EV_1KM_RefSB";
			}
		}
		int j;
		for(j = 0; j < 16; j++){
			if(strcmp(band, kme_1_list[j]) == 0){
				*band_index = j;
				return "EV_1KM_Emissive";
			}
		}
		int k;
		for(k = 0;k < 2; k++){
			if(strcmp(band, m_250_list[k]) == 0){
				*band_index = k;
				return "EV_250_Aggr1km_RefSB";
			}
		}
		int a;
		for(a = 0; a < 5; a++){
			if(strcmp(band, m_500_list[a]) == 0){
				*band_index = a;
				return "EV_500_Aggr1km_RefSB";
			}
		}
		return NULL; 
	}
	else if(strcmp(resolution, "_250m") == 0){
		int i;
		for(i = 0; i < 2; i++){
			if(strcmp(band, m_250_list[i]) == 0){
				*band_index = i;
				return "EV_250_RefSB";
			}
		}
		return NULL;
	}
	else if(strcmp(resolution, "_500m") == 0){
		int i;
		for(i = 0; i < 5; i++){
			if(strcmp(band, m_500_list[i]) == 0){
				*band_index = i;
				return "EV_500_RefSB";
			}
		}
		int j;
		for(j = 0; j < 2; j++){
			if(strcmp(band, m_250_list[j]) == 0){
				*band_index = j;
				return "EV_250_Aggr500_RefSB";
			}
		}
		return NULL;
	}
}

double* get_ceres_rad(hid_t file, char* camera, char* d_name, int* size){
	printf("Reading CERES radiance\n");
	//Path variables
	char* instrument = "CERES";
	char* rad = "Radiances";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* data;
	double curr_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, camera, rad, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(camera) + strlen(rad) + strlen(d_name) + 5, 5);
		printf("granule_name: %s\n", name);
		if(read_first < 0){
			data = af_read(file, dataset_name);
			if(data == NULL){
				continue;
			}
			curr_size = dim_sum(af_read_size(file, dataset_name), 1); 
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_data = af_read(file, dataset_name);
			if(adding_data == NULL){
				continue;
			}
			double new_d_size = dim_sum(af_read_size(file, dataset_name), 1);
			//Reallocating arrays of data
			data = realloc(data, sizeof(double)*(curr_size + new_d_size));
			memcpy(&data[(int)curr_size], adding_data, sizeof(double)*new_d_size);
			curr_size += new_d_size;
			
			free(adding_data);
		}
	}
	*size = curr_size;
	
	//Print statements to verify data's existence	
	if(data != NULL){
		printf("test data: %f\n", data[0]);
		printf("test_data: %f\n", data[1]);
		printf("test data: %f\n", data[2]);
	}	
	return data;
}

double* get_ceres_lat(hid_t file, char* camera, char* d_name, int* size){
	printf("Reading CERES lat\n");
	//Path variables
	char* instrument = "CERES";
	char* rad = "Radiances";
	char* tp = "Time_and_Position";
	char* lat = "Latitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, camera, rad, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(camera) + strlen(rad) + strlen(d_name), 4);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* lat_arr[] = {instrument, name, camera, tp, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(camera) + strlen(tp) + strlen(lat) + 5, 5);
		
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 1);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 1);
			
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;
			
			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	//Print statements to verify data's existence
	if(lat_data != NULL){
		printf("test_lat_data: %f\n", lat_data[0]);
		printf("test_lat_data: %f\n", lat_data[1]);
		printf("test_lat_data: %f\n", lat_data[2]);
	}
	
	return lat_data;
}

double* get_ceres_long(hid_t file, char* camera, char* d_name, int* size){
	printf("Reading CERES long\n");
	//Path variables
	char* instrument = "CERES";
	char* rad = "Radiances";
	char* tp = "Time_and_Position";
	char* longitude = "Longitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, camera, rad, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(camera) + strlen(rad) + strlen(d_name), 4);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* long_arr[] = {instrument, name, camera, tp, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(camera) + strlen(tp) + strlen(longitude) + 5, 5);
		
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 1);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 1);
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	
	if(long_data != NULL){
	printf("test_long_data: %f\n", long_data[0]);
	printf("test_long_data: %f\n", long_data[1]);
	printf("test_long_data: %f\n", long_data[2]);
	}
	
	return long_data;
}

double* get_mop_rad(hid_t file, int* size){
	printf("Reading MOPITT radiance\n");
	//Path variables
	char* instrument = "MOPITT";
	char* d_field = "Data_Fields";
	char* rad = "MOPITTRadiances";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* data;
	double curr_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, d_field, rad};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(d_field) + strlen(rad) + 4, 4);
		printf("granule_name: %s\n", name);
		if(read_first < 0){
			data = af_read(file, dataset_name);
			if(data == NULL){
				continue;
			}
			curr_size = dim_sum(af_read_size(file, dataset_name), 5); 
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_data = af_read(file, dataset_name);
			if(adding_data == NULL){
				continue;
			}
			double new_d_size = dim_sum(af_read_size(file, dataset_name), 1);
			//Reallocating arrays of data
			data = realloc(data, sizeof(double)*(curr_size + new_d_size));
			memcpy(&data[(int)curr_size], adding_data, sizeof(double)*new_d_size);
			curr_size += new_d_size;
			
			free(adding_data);
		}
	}
	*size = curr_size;
	
	//Print statements to verify data's existence	
	if(data != NULL){
		printf("test data: %f\n", data[0]);
		printf("test_data: %f\n", data[1]);
		printf("test data: %f\n", data[2]);
	}	
	return data;
}

double* get_mop_lat(hid_t file, int* size){
	printf("Reading MOPITT lat\n");
	//Path variables
	char* instrument = "MOPITT";
	char* d_field = "Data_Fields";
	char* rad = "MOPITTRadiances";
	char* location = "Geolocation";
	char* lat = "Latitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, d_field, rad};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(d_field) + strlen(rad), 3);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		const char* lat_arr[] = {instrument, name, location, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(location) + strlen(lat) + 4, 4);
		
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 3);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 3);
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;

			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	//Print statements to verify data's existence
	if(lat_data != NULL){
		printf("test_lat_data: %f\n", lat_data[0]);
		printf("test_lat_data: %f\n", lat_data[1]);
		printf("test_lat_data: %f\n", lat_data[2]);
	}
	
	return lat_data;
}

double* get_mop_long(hid_t file, int* size){
	printf("Reading MOPITT longitude\n");
	//Path variables
	char* instrument = "MOPITT";
	char* d_field = "Data_Fields";
	char* rad = "MOPITTRadiances";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	
	int h;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, d_field, rad};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(d_field) +strlen(rad) + 3, 3);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* long_arr[] = {instrument, name, location, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(location) + strlen(longitude) + 4, 4);
		
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 3);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 3);
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	//Print statements to verify data's existence
	if(long_data != NULL){
		printf("test_long_data: %f\n", long_data[0]);
		printf("test_long_data: %f\n", long_data[1]);
		printf("test_long_data: %f\n", long_data[2]);
	}
	return long_data;
}

double* get_ast_rad(hid_t file, char* subsystem, char* d_name, int*size){
	printf("Reading ASTER radiance\n");
	//Path variables
	char* instrument = "ASTER";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* data;
	double curr_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, subsystem, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(d_name) + 5, 4);
		printf("granule_name: %s\n", name);
		if(read_first < 0){
			data = af_read(file, dataset_name);
			if(data == NULL){
				continue;
			}
			curr_size = dim_sum(af_read_size(file, dataset_name), 2); 
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimention
			double* adding_data = af_read(file, dataset_name);
			if(adding_data == NULL){
				continue;
			}
			double new_d_size = dim_sum(af_read_size(file, dataset_name), 2);
			//Reallocating arrays of data
			data = realloc(data, sizeof(double)*(curr_size + new_d_size));
			memcpy(&data[(int)curr_size], adding_data, sizeof(double)*new_d_size);
			curr_size += new_d_size;
			
			free(adding_data);
		}
	}
	*size = curr_size;
	
	//Print statements to verify data's existence
	if(data != NULL){
		printf("test data: %f\n", data[0]);
		printf("test_data: %f\n", data[1]);
		printf("test data: %f\n", data[2]);
	}	
	
	return data;
}

double* get_ast_lat(hid_t file, char* subsystem, char* d_name, int*size){
	printf("Reading ASTER lat\n");
	//Path variables
	char* instrument = "ASTER";
	char* location = "Geolocation";
	char* lat = "Latitude";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, subsystem, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(subsystem) + strlen(d_name) + 4, 3);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		const char* lat_arr[] = {instrument, name, subsystem, location, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(location) + strlen(lat) + 5, 5);
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;

			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	//Print statements to verify data's existence
	if(lat_data != NULL){	
		printf("test_lat_data: %f\n", lat_data[0]);
		printf("test_lat_data: %f\n", lat_data[1]);
		printf("test_lat_data: %f\n", lat_data[2]);
	}
	return lat_data;
}

double* get_ast_long(hid_t file, char* subsystem, char* d_name, int* size){
	printf("Reading ASTER long\n");
	//Path variables
	char* instrument = "ASTER";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(&names[i], name);
		free(name);
	}
	int h;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, subsystem, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(subsystem) + strlen(d_name) + 4, 3);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		const char* long_arr[] = {instrument, name, subsystem, location, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(location) + strlen(longitude) + 5, 5);
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and it's dimention
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	//Print statements to verify data's existence
	if(long_data != NULL){
		printf("test_long_data: %f\n", long_data[0]);
		printf("test_long_data: %f\n", long_data[1]);
		printf("test_long_data: %f\n", long_data[2]);
	}
	
	return long_data;
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
	hsize_t* dims = malloc(sizeof(hsize_t) * ndims);
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	H5Dclose(dataset);	
	H5Sclose(dataspace);
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
	if(strstr(dataset_name, "ASTER") != NULL && strstr(dataset_name, "Geolocation") != NULL){
		//Special case for ASTER geolocation because they are 64bit floating point numbers
		double* data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(double) );
		herr_t status = H5Dread(dataset, ndtype, memspace, memspace, H5P_DEFAULT, data);
		H5Dclose(dataset);	
		H5Sclose(dataspace);
		H5Tclose(dtype);
		H5Tclose(ndtype);
		if(status < 0){
			printf("read error: %d\n", status);
		}
		return data;
	}
	else{
		float* data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(ndtype) );
		double* converted_data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(double) );
		herr_t status = H5Dread(dataset, ndtype, memspace, memspace, H5P_DEFAULT, data);
		int i;
		for(i=0;i < dim_sum(dims, sizeof(dims)/sizeof(hsize_t)); i++){
			converted_data[i] = (double) data[i];
		}
		free(data);
		H5Dclose(dataset);	
		H5Sclose(dataspace);
		H5Tclose(dtype);
		H5Tclose(ndtype);
		if(status < 0){
			printf("read error: %d\n", status);
		}
		return converted_data;
	}
}

int af_write_misr_on_modis(hid_t output_file, double* misr_out, double* modis, int modis_size, int misr_size){
	//Create datafield group
	hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	//Write MODIS first
	hsize_t modis_dim[3];
	modis_dim[0] = 15;
	modis_dim[1] = 1354;
	modis_dim[2] = (modis_size)/15/1354;
	hid_t modis_dataspace = H5Screate_simple(3, modis_dim, NULL);
	hid_t	modis_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t  modis_status = H5Tset_order(modis_datatype, H5T_ORDER_LE);  
    hid_t modis_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", modis_datatype, modis_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    modis_status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, modis);
    H5Sclose(modis_dataspace);
	H5Tclose(modis_datatype);
	H5Dclose(modis_dataset);
    if(modis_status < 0){
    	printf("MODIS write error\n");
    	return -1;
	}
    
    //Write 
    hsize_t misr_dim[2];
	misr_dim[0] = (misr_size) / 1354;
	misr_dim[1] = 1354;
	hid_t misr_dataspace = H5Screate_simple(2, misr_dim, NULL);
	hid_t misr_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t misr_status = H5Tset_order(misr_datatype, H5T_ORDER_LE);  
    hid_t misr_dataset = H5Dcreate2(output_file, "/Data_Fields/misr_out", misr_datatype, misr_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    misr_status = H5Dwrite(misr_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, misr_out);
    H5Sclose(misr_dataspace);
	H5Tclose(misr_datatype);
	H5Dclose(misr_dataset);
	if(misr_status < 0){
		printf("MISR write error\n");
		return -1;
	}
	
	return 1;
	
}

int af_write_mm_geo(hid_t output_file, int geo_flag, double* geo_data, int geo_size){
	//Check if geolocation group exists
	herr_t status = H5Gget_objinfo(output_file, "/Geolocation", 0, NULL);
	if(status != 0){
		hid_t group_id = H5Gcreate2(output_file, "/Geolocation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	char* d_name;
	if(geo_flag == 0){
		d_name = "/Geolocation/Latitude";
	}
	else if(geo_flag == 1){
		d_name = "/Geolocation/Longitude";
	}
	hsize_t     geo_dim[2];
	geo_dim[0] = (geo_size) / 1354;
	geo_dim[1] = 1354;
	hid_t geo_dataspace = H5Screate_simple(2, geo_dim, NULL);
	hid_t geo_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t geo_status = H5Tset_order(geo_datatype, H5T_ORDER_LE);  
    hid_t geo_dataset = H5Dcreate2(output_file, d_name, geo_datatype, geo_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(geo_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, geo_data);
    H5Sclose(geo_dataspace);
	H5Tclose(geo_datatype);
	H5Dclose(geo_dataset);
	
	if(geo_status < 0){
		printf("Geo Data write error\n");
		return -1;
	}
	return 1;
}

hid_t af_open(char* file_path){
	hid_t f = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
	return f;
}

herr_t af_close(hid_t file){
	herr_t ret = H5Fclose(file);
	return ret;
}

/*int main (int argc, char *argv[]){
	//Preset filename here for easy testing 
	char* file_path = "/projects/TDataFus/kent/temp/40-orbit-file/Jun15.2/TERRA_BF_L1B_O69365_F000_V000.h5";
	hid_t file;
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
			//Open file
			file = af_open(file_path);
			if(file < 0){
				printf("File not found\n");
				return -1;
			}
						
			int d_size = 0;
			int* data_pt = &d_size;
			double* data = get_misr_rad(file, argv[3], argv[4], argv[5], data_pt);
			printf("Data size: %d\n", *data_pt);
			
			int lat_size = 0;
			int* lat_pt = &lat_size;
			double* lat_data = get_misr_lat(file, argv[4], lat_pt);
			printf("Lat size: %d\n", *lat_pt);
			
			int long_size = 0;
			int* long_pt = &long_size;
			double* long_data = get_misr_long(file, argv[4], long_pt);
			printf("Long size: %d\n", *long_pt);
			
			if(data != NULL && lat_data != NULL, long_data != NULL){
				printf("MISR Data retrieval successful\n");
			}
			else{
				printf("MISR Data retrieval failed\n");
			}
			
			printf("Retrieving MISR attributes\n");
			void* unit_attr;
			unit_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "Units", 0, unit_attr);
			printf("Unit attr: %s\n", (char*)unit_attr);
			void* fill_attr;
			fill_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "_FillValue", 0, fill_attr);
			printf("FillValue: %f\n", *(float*) fill_attr);
			void* lat_attr;
			lat_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "units", 1, lat_attr);
			printf("lat_units: %s\n", (char*)lat_attr);
			void* long_attr;
			long_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "units", 2, long_attr);
			printf("long_units: %s\n", (char*)long_attr);
			
			herr_t ret = af_close(file);
		}
	}
	else if(strcmp(argv[2], "MODIS") == 0){
		if(argc < 4){
			printf("MODIS Usage: %s filename MODIS resolution(1KM/500m/250m)\n", argv[0]);
		}
		else{
			file = af_open(file_path);
			if(file < 0){
				printf("File not found\n");
				return -1;
			}
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
			
			int data_size = 0;
			int* data_pt = &data_size;
			double* data = get_modis_rad(file, resolution, d_name, data_pt);
			printf("Data size: %d\n", *data_pt);
			
			int lat_size = 0;
			int* lat_pt = &lat_size;
			double* lat_data = get_modis_lat(file, resolution, d_name, lat_pt);
			printf("Lat size: %d\n", *lat_pt);
			
			int long_size = 0;
			int* long_pt = &long_size;
			double* long_data = get_modis_long(file, resolution, d_name, long_pt);
			printf("Long size: %d\n", *long_pt);

			
			if(data != NULL && lat_data != NULL, long_data != NULL){
				printf("MODIS retrieval successful\n");
			}
			else{
				printf("MODIS retrieval failed\n");
			}
			
		
			printf("Retrieving MODIS attributes\n");
			void* unit_attr;
			unit_attr = get_modis_attr(file, resolution, d_name, "units", 0, unit_attr);
			void* fill_attr;
			fill_attr = get_modis_attr(file, resolution, d_name, "_FillValue", 0, fill_attr);
			void* min_attr;
			min_attr = get_modis_attr(file, resolution, d_name, "valid_min", 0, min_attr);
			void* lat_attr;
			lat_attr = get_modis_attr(file, resolution, d_name, "units", 1, lat_attr);
			void* long_attr;
			long_attr = get_modis_attr(file, resolution, d_name, "units", 2, long_attr);
			printf("Unit attr: %s\n", (char*)unit_attr);
			printf("FillValue: %f\n", *(float*) fill_attr);
			printf("valid_min: %f\n", *(float*) min_attr);
			printf("lat unit: %s\n", (char*) lat_attr);
			printf("long unit: %s\n", (char*) long_attr); 
			
			herr_t ret = af_close(file);
		}
	}
	else if(strcmp(argv[2], "CERES") == 0){
		if(argc < 6){
			printf("CERES Usage: %s filename CERES camera radiance(LW/SW/WN/TOT) filtered/unfiltered(F/U)");
		}
		else{
				file = af_open(file_path);
				if(file < 0){
					printf("File not found\n");
					return -1;
				}
				char* d_name = calloc(30, sizeof(char));
				if(strcmp(argv[5], "F") == 0){
					char* f = "_Filtered";
					char* r = "_Radiance";
					strcpy(d_name, argv[4]);
					strncat(d_name, f, strlen(f));
					strncat(d_name, r, strlen(r));
				}
				else{
					char* r = "_Radiance";
					strcpy(d_name, argv[4]);
					strncat(d_name, r, strlen(r));
				}
				
				int data_size = 0;
				int* data_pt = &data_size;
				double* data = get_ceres_rad(file, argv[3], d_name, data_pt);
				printf("Data size: %d\n", *data_pt);
				
				int lat_size = 0;
				int* lat_pt = &lat_size;
				double* lat_data = get_ceres_lat(file, argv[3], d_name, lat_pt);
				printf("Lat size: %d\n", *lat_pt);

				int long_size = 0;
				int* long_pt = &long_size;
				double* long_data = get_ceres_long(file, argv[3], d_name, long_pt);
				printf("Long size: %d\n", *long_pt);
				
				herr_t ret = af_close(file);
		}
	}
	else if(strcmp(argv[2], "MOPITT") == 0){
		file = af_open(file_path);
		if(file < 0){
			printf("File not found\n");
			return -1;
		}
		int data_size = 0;
		int* data_pt = &data_size;
		double* data = get_mop_rad(file, data_pt);
		printf("Data size: %d\n", *data_pt);
		
		int lat_size = 0;
		int* lat_pt = &lat_size;
		double* lat_data = get_mop_lat(file, lat_pt);
		printf("Lat size: %d\n", *lat_pt);

		int long_size = 0;
		int* long_pt = &long_size;
		double* long_data = get_mop_long(file, long_pt);
		printf("Long size: %d\n", *long_pt);
		
		herr_t ret = af_close(file);
	}
	else if(strcmp(argv[2], "ASTER") == 0){
		if(argc < 5){
			printf("ASTER Usage: %s filename ASTER subsystem(TIR/VNIR/SWIR) dataset_name\n");
		}
		else{
				file = af_open(file_path);
				if(file < 0){
					printf("File not found\n");
					return -1;
				}
				
				int data_size = 0;
				int* data_pt = &data_size;
				double* data = get_ast_rad(file, argv[3], argv[4], data_pt);
				if(data != NULL){
					printf("Data size: %d\n", *data_pt);	
				}
				int lat_size = 0;
				int* lat_pt = &lat_size;
				double* lat_data = get_ast_lat(file, argv[3], argv[4], lat_pt);
				printf("Lat size: %d\n", *lat_pt);
				
				int long_size = 0;
				int* long_pt = &long_size;
				printf("going into ast_long\n");
				double* long_data = get_ast_long(file, argv[3], argv[4], long_pt);
				printf("Long size: %d\n", *long_pt);	
				
				herr_t ret = af_close(file);	
		}
	}
	else{
		printf("Invalid instrument\n");
	}
	
	return 0;
}*/

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

double misr_averaging(double window[16]){
	double sum = 0.0;
	double count = 0.0;
	int i;
	for(i = 0; i < 16; i++){
		if(window[i] < 0){
			return -999.0;
		}
		else{
			sum += window[i];
			count += 1;
		}
	}
	return sum/count;
}
