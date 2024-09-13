/**
@file	mg_xtal.h
@brief	Header file for functions to process crystallographic data.
@author 	Bernard Heymann and Samuel Payne
@date	Created: 20061110
@date	Modified: 20240306
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
int			mg_unitcell_vectors(Bmicrograph* mg);
long		mg_generate_reflections(Bmicrograph* mg, Vector3<double> real_size, double resolution);
int			img_mask_reflections(Bimage* p, Bstrucfac* sflist, double radius);
vector<Vector3<double>>	img_find_reflections(Bimage* p, double ref_res, double threshold, long kedge);
vector<Vector3<double>>	symmetrize_reflections(vector<Vector3<double>> loc, long sym, double symdist);
Vector3<double>	analyze_reflections(vector<Vector3<double>>& loc, double s);

