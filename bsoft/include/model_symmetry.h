/**
@file	model_symmetry.h
@brief	Library routines used for model symmetry operations
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250618
**/

#include "Bmodel.h"
#include "symmetry.h"
#include "View2.h"
#include "Transform.h"
#include "UnitCell.h"

// Function prototypes
long		model_find_asymmetric_unit(Bmodel* model, Bsymmetry& sym);
long		models_find_asymmetric_unit(Bmodel* model, Bsymmetry& sym);
//long 		model_apply_point_group(Bmodel* model, Bsymmetry& sym,
//					Vector3<double> origin, View2<double> ref_view, int flags=0);
long		models_apply_point_group(Bmodel* model, Bsymmetry& sym,
					Vector3<double> origin, View2<double> ref_view, int flags=0);
long 		models_symmetrize(Bmodel* model, Bsymmetry& sym);
Transform 	model_find_standard_view(Bmodel* model, Bsymmetry& sym, View2<double> ref_view);
int 		model_orient_to_standard_view(Bmodel* model, Bsymmetry& sym, View2<double> ref_view);
vector<Vector3<double>>	model_symmetry_axes(Bmodel* model);
double		model_symmetry_RMSD(Bmodel* model, Bsymmetry& sym, View2<double> ref_view);
double		model_symmetry_B(Bmodel* model, Bsymmetry& sym, View2<double> ref_view);
int 		models_generate_lattice(Bmodel* model, UnitCell uc, Vector3<long> lattice);

