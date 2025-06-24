/**
@file	model_mol.h
@brief	Library routines for processing molecular models
@author 	Bernard Heymann
@date	Created: 20220215
@date	Modified: 20250623
**/

#include "Bmodel.h"
#include "Bsequence.h"
#include "rwresprop.h"
#include "json.h"

// Function prototypes
string		component_element(Bcomponent* comp, map<string,Bcomptype>& atompar);
JSvalue		model_elements_json(Bmodel* model, map<string,Bcomptype>& atompar);
Bmaterial	model_elements(Bmodel* model, map<string,Bcomptype>& atompar);
Bmaterial	material_from_model(Bmodel* model, string& atompropfile);
Bmaterial	material_from_model(Bmodel* model, string& atompropfile, double density, DensityUnit units);
vector<Bsequence>	models_sequence(Bmodel* model);
long		model_select_residues(Bmodel* model, string res_select);
long		models_select_residues(Bmodel* model, string res_select);
long		model_select_corresponding_CA(Bmodel* model, Bmodel* refmod);
long		models_select_corresponding_CA(Bmodel* model, Bmodel* refmod);
long		model_select_aligned_CA(Bmodel* model, Bmodel* refmod, double gapopen,
				double gapextend, Bresidue_matrix& simat);
long		models_select_aligned_CA(Bmodel* model, Bmodel* refmod, double gapopen,
				double gapextend, Bresidue_matrix& simat);
double		models_fit_CA(Bmodel* model, Bmodel* refmod, int type);
int			models_compare_orientations(Bmodel* model);

