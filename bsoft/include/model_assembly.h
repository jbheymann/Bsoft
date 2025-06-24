/**
@file	model_assembly.h
@brief	Library routines used for model processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250318
**/

#include "rwmodel.h"
#include "rwmolecule.h"

// Function prototypes
//Bmolgroup*	model_assemble(Bmodel* model, string paramfile, int separate);
Bmodel*		model_assemble(Bmodel* model, string paramfile);
//Bmodel*		model_generate_com(Bmolgroup* molgroup);
Bmodel*		model_generate_assembly(vector<string> file_list, string paramfile);
//long		model_find_views(Bmodel* model, Bstring& reffile, Bstring& paramfile);

