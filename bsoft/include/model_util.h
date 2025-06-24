/**
@file	model_util.h
@brief	Library routines used for model processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250617
**/

#include "Bmodel.h"

// Function prototypes
//long		models_process(Bmodel* model, long (modfunc)(Bmodel*));
//long		models_process(Bmodel* model, long i, long (modfunc)(Bmodel*, long));
//long		models_process(Bmodel* model, double d, long (modfunc)(Bmodel*, double));
//long		models_process(Bmodel* model, string str, long (modfunc)(Bmodel*, string str));
long		models_component_count(Bmodel* model);
long		models_list(Bmodel* model);
long		models_list_comp(Bmodel* model);
Bmodel*		models_copy(Bmodel* model);
Bmodel*		models_copy(Bmodel* model, string order);
void		models_copy(Bmodel** model, string order);
long		model_merge(Bmodel* model);
long		model_number_ids(Bmodel* model);
long		model_rename(Bmodel* model, char first_name);
long		models_rename_components(Bmodel* model);
long		models_reorder_circular(Bmodel** model, string first);
int			model_associate(Bmodel* model, string associate_type, string associate_file);
int			model_associate_mass(Bmodel* model, string associate_type, double mass);
int			models_set_comptype_filenames(Bmodel* model, string filename);
long		models_set_component_radius(Bmodel* model, double comprad);
int			models_set_map_filenames(Bmodel* model, string mapfile);
int			models_set_type(Bmodel* model, string set_type);
int			model_change_type(Bmodel* model, string change_type);
double		model_mass(Bmodel* model);
long		model_mass_all(Bmodel* model);
long		models_mass(Bmodel* model);
Vector3<double>	models_center_of_coordinates(Bmodel* model);
Vector3<double>	model_geometric_median(Bmodel* model);
double		model_gyration_radius(Bmodel* model);
double		model_effective_thickness(Bmodel* model);
Vector3<double> 	model_principal_axes(Bmodel* model, Vector3<double>* eigenvec);
Vector3<double> 	model_principal_axes(Bmodel* model, Matrix& eigenvec);
long 		model_principal_axes(Bmodel* model);
long		model_radial_distribution(Bmodel* model, double interval);
long		model_average_components(Bmodel* model, int number);
Bcomponent**	component_get_array(Bmodel* model, long& ncomp);
Vector3<double>	component_plane(vector<Bcomponent*>& comparray, double& offset);
vector<Vector3<double>>	models_calculate_bounds(Bmodel* model);
vector<vector<Bcomponent*>>	models_component_grid(Bmodel* model, Vector3<long>& size,
				Vector3<double>& origin, Vector3<double>& sampling);
double		models_volume(Bmodel* model);
double		models_density(Bmodel* model);
vector<Bcomponent*>	models_get_component_array(Bmodel* model);
vector<vector<Bcomponent*>>	model_split_into_slices(Bmodel* model, double bottom, double top, double thickness);
vector<Bmodel*>	model_split_into_slice_models(Bmodel* model, double bottom, double top, double thickness);
int			model_insert(Bmodel* model, Bmodel* modinsert, double distance) ;

