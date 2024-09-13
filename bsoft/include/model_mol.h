/**
@file	model_mol.h
@brief	Library routines for processing molecular models
@author 	Bernard Heymann
@date	Created: 20220215
@date	Modified: 20230821
**/

#include "Bmodel.h"
#include "json.h"

// Function prototypes
string		component_element(Bcomponent* comp, map<string,Bcomptype>& atompar);
JSvalue		model_elements_json(Bmodel* model, map<string,Bcomptype>& atompar);
Bmaterial	model_elements(Bmodel* model, map<string,Bcomptype>& atompar);
Bmaterial	material_from_model(Bmodel* model, string& atompropfile);
Bmaterial	material_from_model(Bmodel* model, string& atompropfile, double density, DensityUnit units);
