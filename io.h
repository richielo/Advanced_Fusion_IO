/*


    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu


*/
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#define FALSE   0

//HDF5 API operations wrapper
hid_t af_open(char* file_path);
herr_t af_close(hid_t file);
double* af_read(hid_t file, char* dataset_name);
double* af_read_hyperslab(hid_t file, char*dataset_name, int x_offset, int y_offset, int z_offset);
hsize_t* af_read_size(hid_t file, char* dataset_name);
int af_write_misr_on_modis(hid_t output_file, double* misr_out, double* modis, int modis_size, int misr_size);
int af_write_mm_geo(hid_t output_file, int geo_flag, double* geo_data, int geo_size);

//Instrument data retrieval functions
double* get_misr_rad(hid_t file, char* camera_angle, char* resolution, char* radiance, int* size);
double* get_misr_lat(hid_t file, char* resolution, int* size);
double* get_misr_long(hid_t file, char* resolution, int* size);
void* get_misr_attr(hid_t file, char* camera_angle, char* resolution, char* radiance, char* attr_name, int geo, void* attr_pt);
double* get_modis_rad(hid_t file, char* resolution, char* bands[], int band_size, int* size);
double* get_modis_rad_by_band(hid_t file, char* resolution, char* d_name, int* band_index, int* size);
double* get_modis_lat(hid_t file, char* resolution, char* d_name, int* size);
double* get_modis_long(hid_t file, char* resolution, char* d_name, int* size);
double* get_modis_attr(hid_t file, char* resolution, char* d_name, char* attr_name, int geo, void* attr_pt);
char* get_modis_filename(char* resolution, char* band, int* band_index);
double* get_ceres_rad(hid_t file, char* camera, char* d_name, int* size);
double* get_ceres_lat(hid_t file, char* camera, char* d_name, int* size);
double* get_ceres_long(hid_t file, char* camera, char* d_name, int* size);
double* get_mop_rad(hid_t file, int* size);
double* get_mop_lat(hid_t file, int*size);
double* get_mop_long(hid_t file, int* size);
double* get_ast_rad(hid_t file, char* subsystem, char* d_name, int*size);
double* get_ast_lat(hid_t file, char* subsystem, char* d_name, int*size);
double* get_ast_long(hid_t file, char* subsystem, char* d_name, int*size);

//Helper functions
void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size);
double dim_sum(hsize_t* dims, int arr_len);
double float_to_double(float f);
double misr_averaging(double window[16]);
