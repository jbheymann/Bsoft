/**
@file	model_mol.cpp
@brief	Library routines for processing molecular models
@author 	Bernard Heymann
@date	Created: 20220215
@date	Modified: 20250623
**/

#include "model_mol.h"
#include "model_transform.h"
#include "model_compare.h"
#include "seq_align.h"
#include "seq_util.h"
#include "rwmodel_param.h"
#include "Transform.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Determines the element from the type identifier.
@param 	*comp			component.
@param	atompar			atom type parameters.
@return string			element, empty if not found.
**/
string		component_element(Bcomponent* comp, map<string,Bcomptype>& atompar)
{
	string		cel = comp->element();

	if ( atompar.find(cel) == atompar.end() ) {
		cerr << "Warning: Element " << cel << " not found!" << endl;
		cel.clear();
	}

	return cel;
}

/**
@brief 	Calculates the elemental composition.
@param 	*model			model.
@param	atompar			atom type parameters.
@return JSvalue			composition.
**/
JSvalue		model_elements_json(Bmodel* model, map<string,Bcomptype>& atompar)
{
	JSvalue			el(JSobject);
	
	if ( atompar.size() < 1 ) {
		cerr << "Error in model_elements_json: No parameters!" << endl;
		bexit(-1);
	}

	long			nf(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;
	string			cel;
	
	if ( verbose )
		cout << "Adding atom type parameters based on the element" << endl;

    for ( mp = model; mp; mp = mp->next ) {
//    	comp_type_list_kill(mp->type);
 //   	mp->type = NULL;
		mp->clear_types();
		for( comp = mp->comp; comp; comp = comp->next ) {
			cel = component_element(comp, atompar);
			if ( cel.length() ) {
				if ( el.exists(cel) )
					el[cel] += 1;
				else
					el[cel] = 1;
				if ( atompar.find(cel) != atompar.end() ) {
//					cout << "Before conversion:" << endl;
//					atompar[cel].show();
					ct = model->add_type(&atompar[cel]);
					comp->type(ct);
//					cout << "After conversion:" << endl;
//					ct->show();
//					cout << cel << tab << "Z = " << comp->type()->index() << endl;
				} else nf++;
			} else nf++;
		}
	}

	if ( verbose )
		cout << el << endl;
		
	if ( nf )
		cerr << "Component types not found:      " << nf << endl << endl;
	
	return el;
}

Bmaterial	model_elements(Bmodel* model, map<string,Bcomptype>& atompar)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_elements: atompar size = " << atompar.size() << endl;
	
	Bmaterial				material(model->identifier());
	map<string,Bcomptype>&	comp = material.composition();

	long			nf(0);
	Bmodel*			mp;
	Bcomponent*		c;
	Bcomptype*		ct;
	string			cel;

    for ( mp = model; mp; mp = mp->next ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG model_elements: model " << mp->identifier() << tab << mp->component_count_selected() << endl;
//    	comp_type_list_kill(mp->type);
//    	mp->type = NULL;
		mp->clear_types();
		for( c = mp->comp; c; c = c->next ) {
//			cout << c->identifier() << tab << c->description().size() << endl;
			cel = component_element(c, atompar);
//			cout << cel << endl;
			if ( cel.length() ) {
				if ( atompar.find(cel) != atompar.end() ) {
					ct = mp->add_type(&atompar[cel]);
					c->type(ct);
					if ( comp.find(cel) == comp.end() )
						comp[cel] = atompar[cel];
				} else nf++;
				if ( comp.find(cel) == comp.end() )
					comp[cel] = Bcomptype(cel);
				comp[cel].component_count_increment();
			}
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_elements: types = " << comp.size() << endl;

	if ( material.type_count() < 1 ) {
		cerr << "Warning: Material composition is empty!" << endl;
//		bexit(-1);
	}

	if ( nf )
		cerr << "Component types not found:       " << nf << endl << endl;
	
	return material;
}

/**
@brief 	Calculates the elemental composition from a model.
@param 	*model			model.
@param	&atompropfile	file with elemental properties.
@return Bmaterial			material with composition.
**/
Bmaterial	material_from_model(Bmodel* model, string& atompropfile)
{
	map<string,Bcomptype>	atompar = read_atom_properties(atompropfile);
	
	Bmaterial				material = model_elements(model, atompar);

//	if ( verbose )
//		material.show();
		
	return material;
}

/**
@brief 	Calculates the elemental composition from a model.
@param 	*model			model.
@param	&atompropfile	file with elemental properties.
@param	density			density.
@param	units			density units.
@return Bmaterial			material with composition.
**/
Bmaterial	material_from_model(Bmodel* model, string& atompropfile, double density, DensityUnit units)
{
	if ( density < 0.01 ) {
		density = RHO;
		units = DA_A3;
	}

	Bmaterial				material = material_from_model(model, atompropfile);

	material.density(density, units);
	density = material.unit(DA_A3);

	if ( verbose )
		material.show();
	
	if ( material.type_count() < 1 ) {
		cerr << "Error: No component types found!" << endl << endl;
		bexit(-1);
	}
	
	return material;
}

Bsequence	model_sequence(Bmodel* model)
{
	Bsequence	seq(model->identifier());
	seq.description(model->description());
	
	char		g('-');
	string		s, at, res;
	long		i(0), ds, resnum(0);
	Bcomponent*	comp;
	
	map<string, char>	rc = get_res_code_3_1();
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		ds = comp->description().size();
		if ( ds > 1 ) at = comp->description()[1];
		if ( at == "CA" ) {
			if ( ds > 2 ) res = comp->description()[2];
			if ( ds > 4 ) resnum = to_integer(comp->description()[4]);
			if ( i ) {		// Previous resnum
				for ( ++i; i<resnum; ++i )
					s.push_back(g);
			} else {		// Set to first resnum
				for ( i=1; i<resnum; ++i )
					s.push_back(g);
			}
			if ( rc.find(res) != rc.end() )
				s.push_back(rc[res]);
			else
				s.push_back(g);
		}
	}
	
	seq.sequence(s);
	
	if ( verbose )
		cout << model->identifier() << tab << s.length() << tab << s << endl;
	
	if ( verbose & VERB_FULL ) {
		i = 0;
		for ( auto c: s )
			cout << ++i << tab << c << endl;
	}
	
	return seq;
}

/**
@brief 	Retrieves the sequence from from a molecular model.
@param 	*model			model.
@return vector<string>	set of sequence strings.
**/
vector<Bsequence>	models_sequence(Bmodel* model)
{
	vector<Bsequence>	seqs;
	Bmodel*				mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		seqs.push_back(model_sequence(mp));
	
	return seqs;
}

/**
@brief 	Selects residues in a molecule.
@param 	model			model structure.
@param 	res_select		string specifying selection.
@return long				number of selected components.

	Only the first model is used for selection.
	
**/
long		model_select_residues(Bmodel* model, string res_select)
{
	long		nres(0), resnum;
	Bcomponent*	comp;
	
	if ( verbose )
		cout << "Selecting residues:             " << res_select << endl;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		resnum = to_integer(comp->description()[4]);
		if ( nres < resnum ) nres = resnum;
	}
	
	if ( verbose )
		cout << "Number of residues:             " << nres << endl;

	vector<int>	num = select_numbers(res_select, nres);

	for ( nres=0, comp = model->comp; comp; comp = comp->next ) {
		resnum = to_integer(comp->description()[4]);
		if ( num[resnum] ) comp->select(1);
		else comp->select(0);
		nres += comp->select();
	}

	if ( verbose )
		cout << "Number of residues selected:    " << nres << endl << endl;
	
	return model->component_count_selected();
}

/**
@brief 	Selects residues in a set of molecules.
@param 	model			model structure.
@param 	res_select		string specifying selection.
@return long				number of selected components.
	
**/
long		models_select_residues(Bmodel* model, string res_select)
{
	long				nsel(0);
	Bmodel*				mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		nsel += model_select_residues(mp, res_select);

	return nsel;
}

/**
@brief 	Selects the CA atoms that correspond in the  model and reference model.
@param 	model			model structure.
@param 	refmod			reference model.
@return long				number of components selected.

	Only the first model and first reference model is used for selection.
	The reference model is assumed to be the same protein.
	
**/
long		model_select_corresponding_CA(Bmodel* model, Bmodel* refmod)
{
	Bcomponent*		comp = model->comp;
	Bcomponent*		refcomp = refmod->comp;
	long			rnum(0), refnum(0);
	string			at, res, s1, s2;
	
	comp->deselect_all();
	refcomp->deselect_all();

	map<string, char>	rcode = get_res_code_3_1();
	
	// Find the first corresponding CA atoms 
	for ( comp = model->comp; comp ; comp = comp->next ) {
		at = comp->description()[1];
		if ( at == "CA" ) {
			rnum = to_integer(comp->description()[4]);
			res = comp->description()[2];
			for ( ; refcomp && refnum<rnum; refcomp = refcomp->next ) {
				at = refcomp->description()[1];
				if ( at == "CA" ) {
					refnum = to_integer(refcomp->description()[4]);
					if ( refnum == rnum ) {
						comp->select(1);
						refcomp->select(1);
						s1.push_back(rcode[res]);
						res = refcomp->description()[2];
						s2.push_back(rcode[res]);
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Selecting corresponding CA atoms:" << endl;
		cout << "Model:                          " << model->identifier() << endl;
		cout << "Components selected:            " << model->component_count_selected() << endl;
		cout << "Reference:                      " << refmod->identifier() << endl;
		cout << "Reference components selected:  " << refmod->component_count_selected() << endl;
		if ( verbose & VERB_FULL )
			cout << s1 << endl << s2 << endl;
	}
		
	return model->component_count_selected();
}

/**
@brief 	Fits multiple models to a reference model based on the CA atoms for corresponding residues.
@param 	model				list of models.
@param 	refmod				reference model.
@return long				number of components selected.

	Each selected model is processed relative to the reference model.
	
**/
long		models_select_corresponding_CA(Bmodel* model, Bmodel* refmod)
{
	long		nsel(0);
	Bmodel*		mp;
	Bmodel*		mpr;

	if ( model->count() == refmod->count() ) {
		for ( mp = model, mpr = refmod; mp; mp = mp->next, mpr = mpr->next )
			if ( mp->select() )
				nsel += model_select_corresponding_CA(mp, mpr);
	} else {
		for ( mp = model; mp; mp = mp->next )
			if ( mp->select() )
				nsel += model_select_corresponding_CA(mp, refmod);
	}
	
	return nsel;
}

/**
@brief 	Selects the CA atoms that align in the model and reference model.
@param 	model			model structure.
@param 	refmod			reference model.
@param 	gapopen			gap opening penalty.
@param 	gapextend		gap extension penalty.
@param 	&simat			residue similarity matrix.
@return long				number of selected components.

	Only the first model and first reference model is used for selection.
	The two models are first aligned.
	
**/
long		model_select_aligned_CA(Bmodel* model, Bmodel* refmod, double gapopen,
				double gapextend, Bresidue_matrix& simat)
{
	long			i(0), j(0), ia(0), ja(0);
	Bcomponent*		comp = model->comp;
	Bcomponent*		refcomp = refmod->comp;
	string			at;
	
	string			seq1 = model_sequence(model).sequence();
	string			seq2 = model_sequence(refmod).sequence();
	
	pair<string,string>	seq_aln = seq_pair_align(seq1, seq2, gapopen, gapextend, simat);
	
	string			seq1a = seq_aln.first;
	string			seq2a = seq_aln.second;
	
	comp->deselect_all();
	refcomp->deselect_all();
	
	// Find the first corresponding CA atoms 
	for ( comp = model->comp; comp ; comp = comp->next ) {
		at = comp->description()[1];
		if ( at == "CA" ) {
			while ( seq1a[ia] != seq1[i] ) ia++;
			for ( ; refcomp && ja<ia; refcomp = refcomp->next ) {
				at = refcomp->description()[1];
				if ( at == "CA" ) {
					while ( seq2a[ja] != seq2[j] ) ia++;
					if ( seq1a[ia] == seq2a[ja] ) {
						comp->select(1);
						refcomp->select(1);
						cout << comp->description()[4] << tab << refcomp->description()[4] << endl;
					}
					j++;
				}
			}
			i++;
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Selecting aligned CA atoms:" << endl;
		cout << "Components selected:            " << model->component_count_selected() << endl;
		cout << "Reference components selected:  " << refmod->component_count_selected() << endl;
	}
		
	return model->component_count_selected();
}

/**
@brief 	Fits multiple models to a reference model based on the CA atoms for aligned residues.
@param 	model			list of models.
@param 	refmod			reference model.
@param 	gapopen			gap opening penalty.
@param 	gapextend		gap extension penalty.
@param 	&simat			residue similarity matrix.
@return long				number of selected components.

	Each selected model is processed relative to the reference model.
	
**/
long		models_select_aligned_CA(Bmodel* model, Bmodel* refmod, double gapopen,
				double gapextend, Bresidue_matrix& simat)
{
	long		nsel(0);
	Bmodel*		mp;
	Bmodel*		mpr;

	if ( model->count() == refmod->count() ) {
		for ( mp = model, mpr = refmod; mp; mp = mp->next, mpr = mpr->next )
			if ( mp->select() )
				nsel += model_select_aligned_CA(mp, mpr, gapopen, gapextend, simat);
	} else {
		for ( mp = model; mp; mp = mp->next )
			if ( mp->select() )
				nsel += model_select_aligned_CA(mp, refmod, gapopen, gapextend, simat);
	}
	
	return nsel;
}

/**
@brief 	Fits multiple models to a reference model based on the CA atoms for corresponding residues.
@param 	model			list of models.
@param 	refmod			reference model.
@param 	type			type of fitting and transformation.
@return double			RMSD.

	The default is to fit each model the reference model and transform it.
	Type:
		0	fit first model to the first reference, transform all.
		1	fit each model to the first reference, transform each.
		2	fit corresponding models and references, transform each.
	
**/
double		models_fit_CA(Bmodel* model, Bmodel* refmod, int type)
{
	Bmodel*		mp;
	Bmodel*		mpr;
	
	if ( verbose )
		cout << "Fitting models to " << refmod->identifier() << endl;

	Transform		t;
	double			rmsd(0);

	if ( type == 0 ) {
		model_select_corresponding_CA(model, refmod);
		t = model_find_transform(model, refmod);
		if ( verbose ) {
			cout << "Transform for " << model->identifier() << ": " << endl;
			cout << "Shift:                          " << t.trans << " angstrom" << endl;
			cout << "Axis:                           " << t.axis << endl;
			cout << "Angle:                          " << t.angle*180.0/M_PI << " degrees" << endl;
		}
		models_rotate(model, t);
		rmsd = model_compare(model, refmod);
	} else {
		mpr = refmod;							// First reference
		for ( mp = model; mp && mpr; mp = mp->next ) if ( mp->select() ) {
			model_select_corresponding_CA(mp, mpr);
			t = model_find_transform(mp, mpr);
			if ( verbose ) {
				cout << "Transform for " << mp->identifier() << ": " << endl;
				cout << "Shift:                          " << t.trans << " angstrom" << endl;
				cout << "Axis:                           " << t.axis << endl;
				cout << "Angle:                          " << t.angle*180.0/M_PI << " degrees" << endl;
			}
			model_rotate(mp, t);
			rmsd = model_compare(mp, mpr);
			if ( type > 1 ) mpr = mpr->next;	// Corresponding reference
		}
	}
	

	return rmsd;
}

/**
@brief 	Finds the relative orientations between all models in a linked list.
@param 	model			list of models.
@return double			RMSD.

	All the models should be of the same molecule type.
	
**/
int			models_compare_orientations(Bmodel* model)
{
	long			nsel(0);
	double			R;
	Transform		t;
	Bmodel*			mp;
	Bmodel*			mp2;
	Bmodel*			mpc;
	
	cout << "Model1\tModel2\tRMSD" << endl;
	for ( mp = model; mp->next; mp = mp->next ) {
		for ( mp2 = mp->next; mp2; mp2 = mp2->next ) {
			nsel = model_select_corresponding_CA(mp2, mp);
			t = model_find_transform(mp2, mp);
			mpc = mp2->copy();
			model_rotate(mpc, t);
			R = model_compare_corresponding(mp, mpc);
			delete mpc;
			cout << mp->identifier() << tab << mp2->identifier() << tab << nsel << tab
				<< t.axis << tab << t.angle*180.0/M_PI << tab << R << endl;
		}
	}
	
	return 0;
}

