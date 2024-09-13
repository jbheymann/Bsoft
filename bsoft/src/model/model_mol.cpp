/**
@file	model_mol.cpp
@brief	Library routines for processing molecular models
@author 	Bernard Heymann
@date	Created: 20220215
@date	Modified: 20230821
**/

#include "model_mol.h"
#include "rwmodel_param.h"
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
//	cout << cel << endl;

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

	long			nf(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;
	string			cel;
	
	if ( verbose )
		cout << "Adding atom type parameters based on the element" << endl;

    for ( mp = model; mp; mp = mp->next ) {
    	comp_type_list_kill(mp->type);
    	mp->type = NULL;
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
    	comp_type_list_kill(mp->type);
    	mp->type = NULL;
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
