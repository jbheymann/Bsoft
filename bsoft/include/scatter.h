/**
@file	scatter.h
@brief	Functions for calculating electron scattering profiles
@author Bernard Heymann
@date	Created: 20190521
@date	Modified: 20230507
**/

//#include "rwmodel.h"
#include "rwmodel_param.h"
//#include "Bmaterial.h"
#include "ctf.h"
#include "string_util.h"
//#include "utilities.h"

// Function prototypes 
vector<string>	all_elements(map<string,Bcomptype>& types);
map<string, vector<double>>	calculate_elastic_scattering_curves(
				map<string,Bcomptype>& types, double ds, double scut);
map<string, vector<double>>	calculate_potential_curves(map<string,Bcomptype>& types,
				double dr, double rmax);
double		elastic_cross_section_integrated(Bcomptype& ct, CTFparam& ctf);
double		elastic_cross_section(map<string,Bcomptype>& types, double volt);
double		elastic_cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf);
double		elastic_excluded_cross_section(map<string,Bcomptype>& types, CTFparam& ctf);
double		inelastic_cross_section_langmore(long Z, double volt);
double		inelastic_cross_section_langmore(long Z, CTFparam& ctf);
double		inelastic_excluded_cross_section_langmore(long Z, CTFparam& ctf);
double		inelastic_cross_section_langmore(map<string,Bcomptype>& types, double volt);
double		inelastic_cross_section_langmore(map<string,Bcomptype>& types, CTFparam& ctf);
double		inelastic_excluded_cross_section_langmore(map<string,Bcomptype>& types, CTFparam& ctf);
double		cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf);
double		cross_section_half_maximal_frequency(Bcomptype& ct);
double		cross_section_half_maximal_frequency(map<string,Bcomptype>& at);
int			material_cross_section_half_maximal_frequencies(Bmaterial& material);
double		material_elastic_cross_section(Bmaterial& material, CTFparam& ctf);
double		material_cross_section(Bmaterial& material, CTFparam& ctf);
double		material_show_cross_section(Bmaterial& material, CTFparam& ctf);
double		material_excluded_cross_section(Bmaterial& material, CTFparam& ctf);
//double		effective_mean_free_path(map<string,Bcomptype>& types, CTFparam& ctf);
double		effective_mean_free_path(Bmaterial& material, CTFparam& ctf);
double		effective_mean_free_path(vector<Bmaterial>& material, vector<double> fractions, CTFparam& ctf);
double		signal_intensity(Bmaterial& material, double thickness, CTFparam& ctf);
vector<double>	signal_intensity(Bmaterial& material, double thick_step, double thick_max, CTFparam& ctf);
int			write_scattering_curves(string paramfile, string outfile, string selection, double resolution);
long		write_scattering_curves(string outfile, map<string,Bcomptype>& types, double resolution);
long		write_scattering_curves(string outfile, map<string, vector<double>>& scat, double ds);
int			write_potential_curves(string paramfile, string outfile, string selection, double radius);
int			write_potential_curves(string outfile, map<string,Bcomptype>& types, double radius);
