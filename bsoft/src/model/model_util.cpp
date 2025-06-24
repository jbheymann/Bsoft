/**
@file	model_util.cpp
@brief	Library routines used for model processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250617
**/

#include "model_util.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_links.h"
#include "symmetry.h"
#include "matrix_linear.h"
#include "matrix_util.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	modfunc			function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, long (modfunc)(Bmodel*))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next, ++n )
		(modfunc)(mp);
	
	if ( verbose & VERB_PROCESS )
		cout << "Models processed:               " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	i				an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, long i, long (modfunc)(Bmodel*, long))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next, ++n )
		(modfunc)(mp, i);
	
	if ( verbose & VERB_PROCESS )
		cout << "Models processed:               " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	d				an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, double d, long (modfunc)(Bmodel*, double))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next, ++n )
		(modfunc)(mp, d);
	
	if ( verbose & VERB_PROCESS )
		cout << "Models processed:               " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	str				an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long				aggregate number returned by function.

**/
long		models_process(Bmodel* model, string str, long (modfunc)(Bmodel*, string str))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next, ++n )
		(modfunc)(mp, str);
	
	if ( verbose & VERB_PROCESS )
		cout << "Models processed:               " << n << endl << endl;
	
	return n;
}

/**
@brief 	Counts all the components.
@param 	*model		model parameters.
@return long			number of components.
**/
long		models_component_count(Bmodel* model)
{
	long		ncomp(0);
	Bmodel*		mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += mp->component_count();
	
	return ncomp;
}

/**
@brief 	Lists models in table form.
@param 	*model		model parameters.
@return long			number of models.
**/
long		models_list(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			nmod(0), ncomp, nct(0);
	Bmodel*			mp;
	string			type_id, sym, fmap;
	Vector3<double>	coc;

	cout << "Model\tType\tNcomp\tPG\tHand\tFOM\tSelect\tMap\tNumber\tCenter" << endl;
	for ( mp = model; mp; mp = mp->next, nmod++ ) {
		ncomp = mp->component_count();
		nct += ncomp;
		if ( mp->model_type().length() ) type_id = mp->model_type();
		else type_id = "?";
		if ( mp->symmetry().length() ) sym = mp->symmetry();
		else sym = "?";
		if ( mp->mapfile().length() ) fmap = mp->mapfile();
		else fmap = "?";
		coc = mp->center_of_coordinates();
		cout << mp->identifier() << tab << type_id << tab << ncomp << tab << sym << 
			tab << mp->handedness() << tab << mp->FOM() << tab << mp->select() << 
			tab << fmap << tab << mp->image_number() << 
			tab << setprecision(2) << coc << endl;
	}
	cout << "Total\t\t" << nct << endl << endl;
	
	return nmod;
}

/**
@brief 	Lists models with component counts in table form.
@param 	*model			model parameters.
@return long			number of models.
**/
long		models_list_comp(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			i, nmod(0), nct(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;
	string			type_id, sym, fmap;
	vector<string>	comptypelist;
	
	for ( mp = model; mp; mp = mp->next, nmod++ ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) {
			if ( find(comptypelist.begin(), comptypelist.end(), ct->identifier()) == comptypelist.end() ) {
				comptypelist.push_back(ct->identifier());
				nct++;
			}
		}
	}
    
	vector<int>		n(nct, 0);
	vector<int>		ntot(nct, 0);

	cout << "Model";
	for ( auto t: comptypelist ) cout << tab << t;
	cout << endl;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( i=0; i<nct; i++ ) n[i] = 0;
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			i = 0;
			for ( auto t = comptypelist.begin(); t != comptypelist.end() && *t != comp->type()->identifier(); ++t ) ++i;
			if ( i < nct ) n[i]++;
		}
		cout << mp->identifier();
		for ( i=0; i<nct; ++i ) {
			ntot[i] += n[i];
			cout << tab << n[i];
		}
		cout << endl;
		nmod++;
	}

	cout << "Total";
	for ( i=0; i<nct; ++i )
		cout << tab << ntot[i];
	cout << endl << endl;

	return nmod;
}

/**
@brief 	Copies a linked list of models.
@param 	*model		linked list of models.
@return Bmodel*		new linked list of models.
**/
Bmodel*		models_copy(Bmodel* model)
{
	if ( verbose & VERB_FULL )
		cout << "Copying models" << endl;

	vector<Bmodel*>		marr = model->array();
	long				n(marr.size());

#ifdef HAVE_GCD
	__block	vector<Bmodel*>	mnuarr(n);
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t i){
		mnuarr[i] = marr[i]->copy();
	});
#else
	vector<Bmodel*>	mnuarr(n);
#pragma omp parallel for
	for ( long i=0; i<n; ++i )
		mnuarr[i] = marr[i]->copy();
#endif

	Bmodel*		modnu = mnuarr[0];
	Bmodel*		mp = modnu;
		
	for ( long i=1; i<n; ++i, mp = mp->next )
		mp->next = mnuarr[i];
	
	return modnu;
}

/**
@brief 	Copies selected models in a given order.
@param 	*model		linked list of models.
@param 	order		list of identifiers to in order.
@return Bmodel*		new linked list of models.

	Only the models in the order string are copied.
 
**/
Bmodel*		models_copy(Bmodel* model, string order)
{
	Bmodel*			model_nu = NULL;
	Bmodel*			mp;
	Bmodel*			mp_nu = NULL;

	if ( verbose )
		cout << "Copying models:      " << order << endl;

	vector<string>  vs = split(order,',');
	
	for ( auto s: vs ) {
		for ( mp = model; mp && mp->identifier() != s; mp = mp->next ) ;
		if ( mp ) {
			cout << mp->identifier() << endl;
			if ( mp_nu ) mp_nu->add(mp->copy());
			else model_nu = mp_nu = mp->copy();
		}
	}

	if ( model_nu->count() != vs.size() )
		cerr << "Warning in models_copy: Not all models have been copied!" << endl;
	
	return model_nu;
}

/**
@brief 	Copies selected models in a given order.
@param 	**model		pointer to linked list of models.
@param 	order		list of identifiers to in order.

	Only the models in the order string are copied.
 
**/
void		models_copy(Bmodel** model, string order) {
	Bmodel* 	m = models_copy(*model, order);
	delete *model;
	*model = m;
}

/**
@brief 	Merges components from all models into one.
@param 	*model		model parameters.
@return long		number of components.
**/
long		model_merge(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			nct, ncomp, nlink;
	Bmodel*			mp;
	Bcomptype*		ct = model->type;
	Bcomponent*		comp = model->comp;
	Blink*			link = model->link;

	if ( comp ) for ( ; comp->next; comp = comp->next ) ;
	if ( link ) for ( ; link->next; link = link->next ) ;
	
	for ( mp = model->next; mp; mp = mp->next ) {
		if ( ct ) ct->next = mp->type;
		else model->type = ct = mp->type;
		mp->type = NULL;
		if ( ct ) for ( ; ct->next; ct = ct->next ) ;
		if ( comp ) comp->next = mp->comp;
		else model->comp = comp = mp->comp;
		mp->comp = NULL;
		if ( comp ) for ( ; comp->next; comp = comp->next ) ;
		if ( link ) link->next = mp->link;
		else model->link = link = mp->link;
		mp->link = NULL;
		if ( link ) for ( ; link->next; link = link->next ) ;
	}
	
	delete model->next;
	model->next = NULL;
	
	for ( nct=0, ct = model->type; ct; ct = ct->next ) nct++;

	// Renumber all the components
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next )
		comp->identifier() = to_string(++ncomp);

	for ( nlink=0, link = model->link; link; link = link->next ) nlink++;

	if ( verbose & VERB_PROCESS ) {
		cout << "Merged component types:         " << nct << endl;
		cout << "Merged components:              " << ncomp << endl;
		cout << "Merged links:                   " << nlink << endl;
		cout << endl;
	}
		
	return ncomp;
}

/**
@brief 	Add a number to each model id.
@param 	*model		model parameters.
@return long			number of models.

	The intention is to give unique id's to models.

**/
long		model_number_ids(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			n(0);
	Bmodel*			mp;
	string			id;
	
	for ( mp = model; mp; mp = mp->next ) {
		id = mp->identifier() + "_" + to_string(++n);
		mp->identifier() = id;
	}
	
	return n;
}

/**
@brief 	Rename models with alphabetical letters.
@param 	*model		model parameters.
@param 	first_name	letter of first model.
@return long			number of models.

**/
long		model_rename(Bmodel* model, char first_name)
{
	if ( first_name == 0 ) return 0;
	
	long			nmod(0);
	char			letter(first_name);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next, nmod++ ) {
		mp->identifier() = letter++;
		if ( letter > 'Z' ) letter = 'A';
	}
	
	return nmod;
}

/**
@brief 	Rename components.
@param 	*model		model parameters.
@return long			number of components.

	The number of links to a component determines its new name.
	Only the first model is processed.

**/
long		models_rename_components(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	long			i, n(0);
	Bmodel*			mp;
	Bcomponent*		comp;

	string			ctstr[10];
	ctstr[0] = "NIL";
	ctstr[1] = "MON";
	ctstr[2] = "DI";
	ctstr[3] = "TRI";
	ctstr[4] = "TET";
	ctstr[5] = "PEN";
	ctstr[6] = "HEX";
	ctstr[7] = "HEP";
	ctstr[8] = "OCT";
	ctstr[9] = "NON";
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = model->comp; comp; comp = comp->next, n++ ) {
			for ( i=0; comp->link[i]; i++ ) ;
			if ( i > 9 ) i = 9;
//			comp->type = model_add_type_by_id(model, ctstr[i]);
			comp->type(model->add_type(ctstr[i]));
		}
	}
	
	return n;
}

/**
@brief     Reorders models starting from a specified model and looping back.
@param     **model        model parameters, replaced starting from new first model.
@param     first        identifier of new first model.
@return long            number of models.

**/
long        models_reorder_circular(Bmodel** model, string first)
{
    Bmodel*            mod_nu = *model;
    Bmodel*            mp;
    long            nmod(mod_nu->count());

    for ( mp = mod_nu; mp->next && mp->next->identifier() != first; mp = mp->next ) ;
    
    if ( mp->next ) {
        mod_nu = mp->next;
        mp->next = NULL;
        if ( verbose )
            cout << "Reordering models, starting from " << mod_nu->identifier() << endl;
        for ( mp = mod_nu; mp->next; mp = mp->next ) ;
        mp->next = *model;
        *model = mod_nu;
    }
    
    if ( mod_nu->count() != nmod )
        cerr << "Error in models_reorder_circular: The number of models changed!" << endl;
    
    return nmod;
}


/**
@brief 	Associates a model file with a component type.
@param 	*model			the model.
@param 	associate_type	component type.
@param 	associate_file	component file name.
@return int				number of types associated.

	Model files can be coordinates or maps.

**/
int			model_associate(Bmodel* model, string associate_type, string associate_file)
{
	if ( !model ) return 0;
	
	int				n(0);
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Associating component " << associate_type << " with file " << associate_file << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) if ( ct->select() ) {
			if ( ct->identifier() == associate_type ) {
				ct->file_name(associate_file);
				n++;
			}
		}
	}
	
	return  n;
}

/**
@brief 	Associates a mass with a component type.
@param 	*model			model list.
@param 	associate_type	component type.
@param 	mass			component type mass.
@return int				number of types associated.
**/
int			model_associate_mass(Bmodel* model, string associate_type, double mass)
{
	if ( !model ) return 0;
	
	int				n(0);
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Associating component " << associate_type << " with mass = " << mass << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) if ( ct->select() ) {
			if ( ct->identifier() == associate_type ) {
				ct->mass(mass);
				n++;
			}
		}
	}
	
	return  n;
}

/**
@brief 	Sets the filenames of all selected component types to the given string.
@param 	*model		model parameters.
@param 	filename	component file name.
@return int			number of component types set.

	The image numbers are sequentially set as well.

**/
int			models_set_comptype_filenames(Bmodel* model, string filename)
{
	if ( !model ) return 0;
	
	int				n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		n += mp->set_type_filenames(filename);
	
	return  n;
}


/**
@brief 	Set the display radius for all components to a specific value.
@param 	*model		model parameters.
@param 	comprad		component display radius.
@return long			number of components selected.
**/
long		models_set_component_radius(Bmodel* model, double comprad)
{
	if ( !model ) return 0;

	int				n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		n += mp->set_component_radius(comprad);

	return n;
}

/**
@brief 	Sets all the map file names of selected models.
@param 	*model		model parameters.
@param 	mapfile	map file name.
@return int			0.
**/
int			models_set_map_filenames(Bmodel* model, string mapfile)
{
	Bmodel*		mp;
	
	for ( mp = model; mp; mp = mp->next )
		if ( mp->select()) mp->mapfile(mapfile);
	
	return 0;
}

/**
@brief 	Reset the component types.
@param 	*model		model.
@param 	set_type	component type.
@return int			number of models.

	Sets all the component types to the given string.

**/
int			models_set_type(Bmodel* model, string set_type)
{
	if ( !model ) return 0;
	
	int				n;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
//		comp_type_list_kill(mp->type);
//		mp->type = NULL;
		mp->clear_types();
		ct = mp->add_type(set_type);
		for ( comp = mp->comp; comp; comp = comp->next )
			comp->type(ct);
	}
	
	return  n;
}

/**
@brief 	Change a component type name.
@param 	*model			model.
@param 	change_type	component type.
@return int				number of models.

	Sets all the component types to the given string.

**/
int			model_change_type(Bmodel* model, string change_type)
{
	if ( !model ) return 0;
	
	int				n;
	Bmodel*			mp;
	Bcomptype*		ct;
	
	vector<string>	t = split(change_type, ',');

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		for ( ct = mp->type; ct; ct = ct->next )
			if ( ct->identifier() == t[0] ) ct->identifier(t[1]);
	}
	
	return  n;
}

/**
@brief 	Calculates the mass of a model from component masses.
@param 	*model		model parameters.
@return double		model mass.

	The component type masses must be provided.
	Only the first model in the list is processed.

**/
double		model_mass(Bmodel* model)
{
	double			mass(0);

	if ( !model ) return mass;
	if ( !model->select() ) return mass;
	
	Bcomptype*		ct = NULL;
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next ) {
		ct = comp->type();
		if ( ct ) mass += ct->mass();
	}
	
	model->mass(mass);
	
	if ( verbose & VERB_FULL )
		cout << "Model: " << model->identifier() << "  Mass = " << mass << endl;
	
	return mass;
}

/**
@brief 	Calculates the masses of all the models in the list.
@param 	*model		linked list of model parameters.
@return long		number of selected models.

	The component type masses must be provided.

**/
long		model_mass_all(Bmodel* model)
{
	long			nmod(0);
	Bmodel*			mp;
	
	if ( verbose )
		cout << "Model\tMass" << endl;
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		if ( verbose )
			cout << mp->identifier() << tab << mp->mass() << endl;
		nmod++;
	}
	if ( verbose )
		cout << endl;
	
	return nmod;
}

/**
@brief 	Calculates the masses of all the models in the list.
@param 	*model		linked list of model parameters.
@return double		mass in Dalton.

	The component type masses must be provided.

**/
long		models_mass(Bmodel* model)
{
	double			mass(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		if ( mp->select() )
			mass += mp->mass();
	
	return mass;
}

/**
@brief 	Calculates the center-of-coordinates of a list of models.
@param 	*model			model parameters.
@return Vector3<double>	center-of-mass.

**/
Vector3<double>	models_center_of_coordinates(Bmodel* model)
{
	Vector3<double>	com;

	if ( !model ) return com;
	if ( !model->select() ) return com;
	
	long		ncomp(0);
	Bmodel* 	mp;
	Bcomponent*	comp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			com += comp->location();
			ncomp++;
		}
	}
	
	com /= ncomp;
	
	return com;
}

/**
@brief 	Calculates the sum of distances of model components to a location.
@param 	*model			model parameters.
@param 	loc				reference location.
@return double			sum of distances.

	Only the first model in the list is processed.

**/
double			model_distance_sum(Bmodel* model, Vector3<double> loc)
{
	long			n(0);
	double			ds(0);
	Bcomponent*		comp;
	
	for ( comp = model->comp; comp; comp = comp->next, ++n )
		ds += loc.distance(comp->location());

	if ( n ) ds /= n;
	
	return ds;
}

/**
@brief 	Calculates the geometric median estimate of model components.
@param 	*model			model parameters.
@param 	pgm				previous geometric median estimate.
@return Vector3<double>	new geomtric median estimate.

	Only the first model in the list is processed.

**/
Vector3<double>	model_geometric_median_estimate(Bmodel* model, Vector3<double> pgm)
{
	double			d, ds(0);
	Vector3<double>	gm;
	Bcomponent*		comp;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		d = 1.0/pgm.distance(comp->location());
		if ( d > 0.001 ) {
			ds += d;
			gm += comp->location() * d;
		}
	}
	
	if ( ds ) gm /= ds;

	return gm;
}

/**
@brief 	Calculates the geometric median of model components.
@param 	*model			model parameters.
@return Vector3<double>	geometric median.

	Only the first model in the list is processed.
	Based on Weiszfeld’s method.

**/
Vector3<double>	model_geometric_median(Bmodel* model)
{
	long			i, iter(1000);
	double			tol(0.01), dd, ddd(1), ds;
	Vector3<double>	gm, pgm;
	
	pgm = model->center_of_coordinates();
	ddd = dd = pgm.length();
	ds = model_distance_sum(model, gm);
	if ( verbose & VERB_PROCESS ) {
		cout << "Iter\tGMx\tGMy\tGMz\t∆d\tDistSum" << endl;
		cout << 0 << tab << pgm << tab << ddd << tab << ds << endl;
	}
	
	for ( i=0; i<iter && ddd > tol; ++i ) {
		ddd = dd;
		gm = model_geometric_median_estimate(model, pgm);
		dd = gm.distance(pgm);
		ddd -= dd;
		pgm = gm;
		ds = model_distance_sum(model, gm);
		if ( verbose & VERB_PROCESS )
			cout << i+1 << tab << gm << tab << ddd << tab << ds << endl;
	}
	
	return gm;
}

/**
@brief 	Calculates the radius of gyration for a model.
@param 	*model		model parameters.
@return double		radius of gyration.

	Only the first model in the list is processed.

**/
double		model_gyration_radius(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	int				n;
	double			d, R;
	Bcomponent*		comp;
	
	Vector3<double>	com = model->center_of_coordinates();
	
	for ( n = 0, R = 0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		d = (comp->location() - com).length();
		R += d*d;
		n++;
	}
	
	if ( n ) R = sqrt(R/n);
	
	return R;
}

/**
@brief 	Calculates the effective thickness for a model.
@param 	*model		model parameters.
@return double		effective thickness.

	The variance in z is related to the effective thickness:
		sigma^2 = (1/12)*thickness^2

**/
double		model_effective_thickness(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	int				n(0);
	double			za(0), zv(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			za += comp->location()[2];
			zv += comp->location()[2]*comp->location()[2];
			n++;
		}
	}
	
	if ( n ) {
		za /= n;
		zv = zv/n - za*za;
	}
	
	return sqrt(12*zv);
}


/**
@brief	Calculates the principal axes of a model.
@param 	*model			model structure.
@param 	*eigenvec		eigen vectors (can be NULL).
@return Vector3<double>	3-valued vector of principal axes.

	Only the first model in the list is processed.

**/
Vector3<double> 	model_principal_axes(Bmodel* model, Vector3<double>* eigenvec)
{
	Vector3<double>	eigenval;

	if ( !model  ) return eigenval;
	if ( !model->select() ) return eigenval;
	
	double			summass(0);
	Vector3<double>	loc, vec, vec2, vecx;
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		loc = comp->location();
		vec += loc;					// Sums
		vec2 += loc*loc;		// Square sums
		vecx[0] += loc[0]*loc[1];	// Cross-term sums
		vecx[1] += loc[0]*loc[2];
		vecx[2] += loc[1]*loc[2];
		summass += 1;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Model:                          " << model->identifier() << endl;
		cout << "Number of components:           " << summass << endl;
	}
	
	if ( summass < 1 ) {
		cout << "Error: No components found!" << endl << endl;
		return vec;
	}
	
	vec /= summass;
	vec2 /= summass;
	vecx /= summass;
	
	return principal_axes(vec, vec2, vecx, eigenvec);
}

Vector3<double> 	model_principal_axes(Bmodel* model, Matrix& eigenvec)
{
	Vector3<double>			pax;

	if ( !model  ) return pax;
	if ( !model->select() ) return pax;

	vector<Vector3<double>>	coor;
	Vector3<double>			loc;
	Bcomponent*				comp;

	for ( comp = model->comp; comp; comp = comp->next )
		if ( comp->select() ) {
			loc = comp->location();
			coor.push_back(loc);
		}
	
	if ( verbose & VERB_PROCESS )
		cout << "Model:                          " << model->identifier() << endl;
	
	if ( coor.size() < 1 ) {
		cout << "Error: No components found for model " <<
			model->identifier() << "!" << endl << endl;
		return pax;
	}
	
	return principal_axes(coor, eigenvec);
}

long		model_principal_axes(Bmodel* model)
{
	if ( !model  ) return -1;

	Vector3<double>		pax;
	Bmodel*				mp;
	Matrix				eigenvec(3,3);
	
	if ( verbose )
		cout << "Principal axes:" << endl << "Model\tMajor\tMiddle\tMinor" << endl;
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		pax = model_principal_axes(mp, eigenvec);
		if ( verbose )
			cout << mp->identifier() << tab << pax[0] << tab << pax[1] << tab << pax[2] << endl;
	}
	
	if ( verbose )
		cout << endl;
	
	return 0;
}


/**
@brief	Calculates the radial distribution function of a model.
@param 	*model		model structure.
@param 	interval	interval between bins.
@return long		0.

	Only the first model in the list is processed.

**/
long		model_radial_distribution(Bmodel* model, double interval)
{
	if ( interval <= 0 ) interval = 1;	
	if ( !model  ) return 0;
	if ( !model->select() ) return 0;
	
	Vector3<double>	min(1e10,1e10,1e10), max(-1e10,-1e10,-1e10);
	Bcomponent*		comp1, *comp2;

	for ( comp1 = model->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		min = min.min(comp1->location());
		max = max.max(comp1->location());
	}
	
	int				r, rmax = (int) (max.distance(min)/interval + 1);
	double			d;
	int*			rdf = new int[rmax];
	for ( r=0; r<rmax; r++ ) rdf[r] = 0;
	
	for ( comp1 = model->comp; comp1->next; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( comp2 = comp1->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			d = comp1->location().distance(comp2->location());
			r = (int) (d/interval + 0.5);
			if ( r < rmax ) rdf[r]++;
		}
	}
	
//	if ( verbose ) {
		cout << "Radius\tCount" << endl;
		for ( r=0; r<rmax; r++ )
			cout << r*interval << tab << rdf[r] << endl;
		cout << endl;
//	}
	
	delete[] rdf;
	
	return 0;
}

/**
@brief     Averages sequential components.
@param 	*model		model structure to be modified.
@param 	number		number of components to average.
@return long		number of remaining components.

	Only the first component in each set with modified coordinates is kept.

**/
long		model_average_components(Bmodel* model, int number)
{
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp_avg;

	if ( verbose & VERB_PROCESS )
		cout << "Averaging every " << number << " of selected components" << endl;
		
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		comp_avg = NULL;
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( !comp_avg ) {
				comp_avg = comp;
				comp_avg->select(1);
			} else {
				comp->select(0);
				comp_avg->shift(comp->location());
				comp_avg->select_increment();
				if ( comp_avg->select() == number ) {
					comp_avg->scale(1.0/comp_avg->select());
					comp_avg = NULL;
				}
			}
		} else comp->select(1);
		if ( comp_avg ) comp_avg->location(comp_avg->location() / comp_avg->select());
	}
	
	return models_delete_non_selected(&model);
}

/**
@brief	Generates an array of pointers to component structures.
@param 	*model			model structure.
@param 	&n				pointer to number of comps found.
@return Bcomponent**	array of pointers to components.
**/
Bcomponent**	component_get_array(Bmodel* model, long& n)
{
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) n++;

	Bcomponent**	comparray = new Bcomponent*[n];

	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() )
			comparray[n++] = comp;

	return comparray;
}

/**
@brief 	Calculates a plane through an array of components.
@param 	**comparray		array of components.
@param 	&offset			offset from plane.
@return Vector3<double>	plane normal.

	A plane is fit through the polygon vertices and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.
	The polygon planarity is defined as the root-mean-square-deviation from 
	the fitted plane.

**/
Vector3<double>	component_plane(vector<Bcomponent*>& comparray, double& offset)
{
    long     			i, j, n;
	vector<double>		b(3);
	Matrix				a(3,3);
	Vector3<double>		loc, normal, center;

	for ( i=0; i<3; i++ ) b[i] = 0;
	
	center = 0;
	for ( n=0; n<comparray.size() && comparray[n]; n++ ) {
		loc = comparray[n]->location();
		center += loc;
		for ( i=0; i<3; i++ ) {
			for ( j=0; j<=i; j++ )
				a[3][i] += loc[i]*loc[j];
			b[i] += loc[i];
		}
	}
	center /= n;
	
	for ( i=1; i<3; i++ )
		for ( j=0; j<i; j++ )
			a[3][j] = a[3][i];
	
	offset = 0;
	if ( a[0][0] == 0 ) {
		normal = Vector3<double>(1, 0, 0);
	} else if ( a[1][1] == 0 ) {
		normal = Vector3<double>(0, 1, 0);
	} else if ( a[2][2] == 0 ) {
		normal = Vector3<double>(0, 0, 1);
	} else {
		a.LU_decomposition(b);
		normal = Vector3<double>(b[0], b[1], b[2]);
		normal.normalize();
	}

	offset = center.scalar(normal);

	return normal;
}

/**
@brief	Calculates the bounds of a list of models.
@param 	*model	 	model list.
@return vector<Vector3<double>>	minimum and maximum bounds.

**/
vector<Vector3<double>>	models_calculate_bounds(Bmodel* model)
{
	vector<Vector3<double>> bounds;
	Vector3<double>			mn(1e30,1e30,1e30), mx(-1e30,-1e30,-1e30);
	Bmodel*					mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		mp->calculate_bounds();
		mn = mn.min(mp->minimum());
		mx = mx.max(mp->maximum());
	}
	
	bounds.push_back(mn);
	bounds.push_back(mx);
	return bounds;
}

/**
@brief	Generates lists of atoms based on a grid.
@param 	*model	 	model list.
@param 	size		size of grid.
@param 	origin		origin of grid.
@param 	sampling	spacing in each dimension.
@return vector<vector<Bcomponent*>>	array of component arrays.

	The goal is to fit all the components within the grid boundaries.
	Components located outside the grid will be added to the edges.

**/
vector<vector<Bcomponent*>>	models_component_grid(Bmodel* model, Vector3<long>& size,
			Vector3<double>& origin, Vector3<double>& sampling)
{
	if ( sampling.volume() < 1 ) {
		cerr << "Error in models_component_grid: sampling must be specified!" << endl;
		bexit(-1);
	}

	long			i;
	Vector3<double>	mn(1e30,1e30,1e30), mx(-1e30,-1e30,-1e30), s1, loc;
	
	Bmodel*			mp;
	Bcomponent*		comp;
	
	if ( size.volume() < 1 ) {
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			mp->calculate_bounds();
			mn = mn.min(mp->minimum());
			mx = mx.max(mp->maximum());
		}
		s1 = mx - mn;
		origin = -mn;
		size = s1/sampling;
	}

	long			gridvol = (long) size.volume();
	Vector3<long>	size1 = size - 1;
	vector<vector<Bcomponent*>>	grid(gridvol);

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			loc = (comp->location() + origin)/sampling;
			loc = loc.max(0);
			loc = loc.min(size1);
			i = (loc[2]*size[1] + loc[1])*size[0] + loc[0];
			grid[i].push_back(comp);
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG models_component_grid: Done!" << endl;
	
	return grid;
}

/**
@brief	Calculates an estimate of the volume of a model.
@param 	*model 				model.
@return double				volume in angstrom^3.
**/
double		models_volume(Bmodel* model)
{
	double			vol(0);
	Vector3<long> 	size;
	Vector3<double>	origin, sampling(1,1,1);
	
	vector<vector<Bcomponent*>>	grid = models_component_grid(model, size, origin, sampling);
	
	for ( auto c: grid )
		vol += c.size();
	
	return vol;
}

/**
@brief	Calculates an estimate of the density of a model.
@param 	*model 				model.
@return double				density in Dalton/angstrom^3.
**/
double		models_density(Bmodel* model)
{
	return	models_mass(model)/models_volume(model);
}


/**
@brief	Generates an array of pointers to model components.
@param 	*model 				model.
@return vector<Bcomponent*>	array of pointers to components.
**/
vector<Bcomponent*>	models_get_component_array(Bmodel* model)
{
	Bmodel*				mp;
	Bcomponent*			comp;
	vector<Bcomponent*>	carr;

    for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		for( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() )
			carr.push_back(comp);

	return carr;
}

/**
@brief	Splits the components of models into slices.
@param 	*model 				model.
@param 	bottom 				minimum coordinate in z.
@param 	top 					maximum coordinate in z.
@param 	thickness 			slice thickness.
@return vector<vector<Bcomponent*>>	sets of arrays of pointers to components.
**/
vector<vector<Bcomponent*>>	model_split_into_slices(Bmodel* model, double bottom, double top, double thickness)
{
	long			i, n((top - bottom)/thickness);
	vector<vector<Bcomponent*>>	comp_slice(n);
	
	if ( verbose ) {
		cout << "Splitting model into slices:" << endl;
		cout << "Bottom:                         " << bottom << " A" << endl;
		cout << "Top:                            " << top << " A" << endl;
		cout << "Slice thickness:                " << thickness << " A" << endl;
		cout << "Number of slices:               " << n << endl;
	}
	
	Bmodel*			mp;
	Bcomponent*		comp;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			i = (comp->location()[2] - bottom)/thickness;
			if ( i>=0 && i<n )
				comp_slice[i].push_back(comp);
		}
	}
	
	if ( verbose ) {
		long		ncomp(0);
		cout << "Slice\tComponents" << endl;
		for ( i=0; i<n; ++i ) {
			ncomp += comp_slice[i].size();
			cout << i+1 << tab << comp_slice[i].size() << endl;
		}
		cout << "Total number of components:     " << ncomp << endl << endl;
	}
	
	return comp_slice;
}

/**
@brief	Splits the components of models into slices and write to new models.
@param 	*model 				model.
@param 	bottom 				minimum coordinate in z.
@param 	top 					maximum coordinate in z.
@param 	thickness 			slice thickness.
@return vector<Bmodel*>	sets of arrays of pointers to components.
**/
vector<Bmodel*>	model_split_into_slice_models(Bmodel* model, double bottom, double top, double thickness)
{
	long			i(0), n((top - bottom)/thickness);
	vector<Bmodel*>	model_slice(n);
	
	if ( verbose ) {
		cout << "Splitting model into slices:" << endl;
		cout << "Bottom:                         " << bottom << " A" << endl;
		cout << "Top:                            " << top << " A" << endl;
		cout << "Slice thickness:                " << thickness << " A" << endl;
		cout << "Number of slices:               " << n << endl;
	}
	
	for ( auto& m: model_slice ) m = new Bmodel(++i);
	
	Bmodel*			mp;
	Bcomponent*		comp;
	vector<Bcomponent*>	comp_slice(n, NULL);

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			i = (comp->location()[2] - bottom)/thickness;
			if ( i>=0 && i<n ) {
				if ( comp_slice[i] ) comp_slice[i] = comp_slice[i]->add(comp);
				else comp_slice[i] = model_slice[i]->add_component(comp);
				comp_slice[i]->description() = comp->description();
//				cout << comp_slice[i]->description().size() << endl;
			}
		}
	}
	
	if ( verbose ) {
		long		ncomp(0);
		cout << "Slice\tComponents" << endl;
		for ( auto& m: model_slice ) {
			ncomp += m->component_count();
			cout << m->identifier() << tab << m->component_count() << endl;
		}
		cout << "Total number of components:     " << ncomp << endl << endl;
	}
	
	return model_slice;
}


/**
@brief	Inserts one model into another.
@param 	*model		 	model to be modified.
@param 	*modinsert 		model to insert. (deallocated)
@param 	distance		cutoff distance to remove atoms.
@return int				0.

	Components overlapping in the receiving molecule group are deleted.
	The footprint of the models being inserted is calculated on a grid
	and all components within this footprint is tested for deletion.
	Note: The model list is transferred from the insertion group to 
		the main group and the insertion group is deallocated.

**/
int			model_insert(Bmodel* model, Bmodel* modinsert, double distance) 
{
	if ( distance <= 0 ) distance = 2;  // Default
	if ( distance < 1 ) distance = 1;   // Limits on cutoff distance in angstrom
	if ( distance > 5 ) distance = 5;
	
	long			i, delete_comp, ncompdel(0), nmoddel(0);
	long			ii, x, y, z, xx, yy, zz, ix, iy, iz;
	Vector3<double>	sampling(distance, distance, distance);
	Bmodel			*m, *pm;
	Bcomponent		*c, *pc;
	Vector3<double>	box = model->maximum() - model->minimum();
	Vector3<long>	gridsize((long) (box[0]/sampling[0] + 0.001), 
		(long) (box[1]/sampling[1] + 0.001), (long) (box[2]/sampling[2] + 0.001));
	gridsize = gridsize.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/gridsize[i] + 0.001;
//	long	gridvol = (long) gridsize.volume();
	
	if ( verbose )
		cout << "Inserting a model and deleting overlapping components" << endl;
	
	if ( verbose & VERB_PROCESS )
		cout << "Distance cutoff:                " << distance << " A" << endl;

	Vector3<double>	gridori;
	vector<vector<Bcomponent*>>	grid = models_component_grid(modinsert, gridsize, gridori, sampling);
	
	// Find the atoms under the footprint to be deleted
	for ( m = pm = model; m; ) {
		for( c = pc = m->comp; c; ) {
			delete_comp = 0;
			x = (long) ((c->location()[0] - model->minimum()[0])/sampling[0]);
			y = (long) ((c->location()[1] - model->minimum()[1])/sampling[1]);
			z = (long) ((c->location()[2] - model->minimum()[2])/sampling[2]);
			if ( x >=0 && x < gridsize[0] && y >= 0 && y < gridsize[1] && z >= 0 && z < gridsize[2] ) {
				i = (z*gridsize[1] + y)*gridsize[0] + x;
				for ( zz=z-1; zz<=z+1; zz++ ) {
					iz = zz;
					if ( iz < 0 ) iz += gridsize[2];
					if ( iz >= gridsize[2]) iz -= gridsize[2];
					for ( yy=y-1; yy<=y+1; yy++ ) {
						iy = yy;
						if ( iy < 0 ) iy += gridsize[1];
						if ( iy >= gridsize[1] ) iy -= gridsize[1];
						for ( xx=x-1; xx<=x+1; xx++ ) {
							ix = xx;
							if ( ix < 0 ) ix += gridsize[0];
							if ( ix >= gridsize[0] ) ix -= gridsize[0];
							ii = (iz*gridsize[1] + iy)*gridsize[0] + ix;
							for ( auto& c2: grid[ii] ) {
								if ( c->location().distance(c2->location()) < distance )
									delete_comp = 1;
							}
						}
					}
				}
			}
			if ( delete_comp ) {
				if ( verbose & VERB_FULL )
					cout << "Removing a component" << endl;
				if ( c == m->comp ) {
					m->comp = pc = c->next;
					delete c;
					c = m->comp;
				} else {
					pc->next = c->next;
					delete c;
					c = pc->next;
				}
				ncompdel++;
			} else {
				pc = c;
				if ( c ) c = c->next;
			}
		}
		if ( !m->comp ) {
			if ( verbose & VERB_FULL )
				cout << "Removing a model" << endl;
			if ( model == m ) {
				model = pm = m->next;
				delete m;
				m = model;
			} else {
				pm->next = m->next;
				delete m;
				m = pm->next;
			}
			nmoddel++;
		} else {
			pm = m;
			if ( m ) m = m->next;
		}
	}

	grid.clear();

	if ( verbose & VERB_PROCESS )
		cout << "Models and components deleted: " << nmoddel << " " << ncompdel << endl;
	
	// Add the new molecules
	if ( model ) {
		for ( m = model; m->next; m = m->next ) ;
		m->next = modinsert;
	} else  model = modinsert;
	modinsert = NULL;
	
	return 0;
}

