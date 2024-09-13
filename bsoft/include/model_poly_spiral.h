/**
@file	model_poly_spiral.h
@brief	Functions to generate polyhedra using the spiral algorithm.
@author 	Bernard Heymann
@date	Created: 20071127
@date	Modified: 20210124
**/

#include "Bmodel.h"

// Function prototypes 
Bmodel*		model_poly_spiral(string seq, int valence, int requirements);
int			model_polyhedron_check(Bmodel* model, int valence);
Bmodel*		model_poly_gen_sequence(string seq, int valence, int enantiomorph, int requirements, int nm);
Bmodel*		model_poly_gen_sequence(string seq, int valence, int enantiomorph, int requirements, int nm, vector<double>& table);
Bmodel*		model_poly_gen_permutations(int vertices, int valence, int enantiomorph);
Bmodel*		model_poly_gen_cone(int tip, int body, int base, int valence, int enantiomorph, int requirements);
Bmodel*		model_poly_gen_lozenge(int ttop, int tbody, int valence, int enantiomorph, int requirements);
Bmodel*		model_poly_gen_coffin(int ttop, int tbody, int tbase, int valence, int enantiomorph, int requirements);
Bmodel*		model_poly_gen_coffin_loose(int ttop, int tbody, int tbase, int valence, int enantiomorph, int requirements);
Bmodel*		model_poly_gen_coffin_jiggle(int ttop, int tbody, int tbase, int valence, int enantiomorph, int requirements);
Bmodel*		model_poly_gen_3part(string stip, string sbody, string sbase, int valence, int enantiomorph, int requirements);
Bmodel*		model_poly_gen_move_pentagons(string seq, int valence, int enantiomorph, int requirements);

