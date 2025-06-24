/**
@file	model_monte.cpp
@brief	Functions for a monte carlo metroplis algorithm to energy minimize model positions.
@author 	Bernard Heymann
@date	Created: 20041230
@date	Modified: 20250620
**/

#include "model_mechanics.h"
#include "model_monte.h"
#include "model_transform.h"
#include "model_util.h"
#include "Bimage.h"
#include "Matrix.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			display_energy_headers(Bmodparam& md)
{
	cout << "#";
	if ( md.Klink ) cout << "\tElink";
	if ( md.Kangle ) cout << "\tEangle";
	if ( md.Kdistance ) cout << "\tEdistance";
	if ( md.Kelec ) cout << "\tEelec";
	if ( md.Ksep ) cout << "\tEsep";
	if ( md.Kpoint ) cout << "\tEpoint";
	if ( md.Kmap ) cout << "\tEmap";
	cout << "\tEpot\tdEpot" << endl;
	
	return 0;
}

int			display_energy(long iter, Bmodparam& md, double dE)
{
	cout << iter;
	if ( md.Klink ) cout << tab << md.Elink;
	if ( md.Kangle ) cout << tab << md.Eangle;
	if ( md.Kdistance ) cout << tab << md.Edistance;
	if ( md.Kelec ) cout << tab << md.Eelec;
	if ( md.Ksep ) cout << tab << md.Esep;
	if ( md.Kpoint ) cout << tab << md.Epoint;
	if ( md.Kmap ) cout << tab << md.Emap;
	cout << tab << md.Epot << tab << dE << endl;
	
	return 0;
}

/**
@brief 	Monte Carlo fit of a model to a map.  
@param 	*model			model list.
@param 	&md				model dynamics structure.
@param 	*map			map.
@param 	beta			equivalent of 1/kT.
@param 	max_angle		maximum allowed angular step size.
@param 	max_shift		maximum allowed shift.
@param 	max_iter		maximum number of iterations.
@fn 	(Efunc)(Bmodel*, Bimage*, Bmd*)	energy function.
@fn 		(Tfunc)(Bmodel*, double, double)		transformation function.
@return	Bmodel*			transformed coordinates.

	Using metropolis criterion:
		if ( ∆E > 0 & exp(-ß∆E) > random[0,1] ) accept change.

**/
Bmodel*	monte_carlo_metropolis(Bmodel* model, Bmodparam& md, Bimage* map, 
				double beta, double max_angle, double max_shift, long max_iter, 
				double (Efunc)(Bmodel*, Bimage*, Bmodparam&),
				int (Tfunc)(Bmodel*, double, double, int))
{
	if ( !map ) {
		cerr << "Error: No map specified!" << endl;
		bexit(-1);
	}

	int				done(0), accept, select_min(0);
	long			i, isel(0), cycle(0), naccept(0);
	double			r, E, pE, bestE, Z(0), Ze, Zi, Emax;
	double			irm = 1.0/get_rand_max();
	Vector3<double>	gloc = models_center_of_coordinates(model);
	Bmodel*			modcopy = models_copy(model);
	Bmodel*			mp = NULL;
	Bmodel*			mpold = NULL;
	
	vector<double>	Elist(max_iter+1,0);

	long			nmod(0);
	for ( mp = modcopy; mp; mp = mp->next ) nmod++;
	vector<Vector3<double>>	iloc(nmod);
	for ( i=0, mp = model; mp; mp = mp->next, ++i )
		iloc[i] = mp->center_of_coordinates();

	if ( verbose & VERB_RESULT ) {
		cout << "Monte carlo metropolis energy minimization:" << endl;
		cout << "Number of models:               " << nmod << endl;
		cout << "Map file:                       " << map->file_name() << endl;
		cout << "Coordinate minima and maxima:   " << model->minimum() << " " << model->maximum() << endl;
		cout << "Bounding box size:              " << md.max - md.min << endl;
		cout << "Bounding box center:            " << (md.max + md.min)*0.5 << endl;
		cout << "Klink:                          " << md.Klink << endl;
		cout << "Kangle:                         " << md.Kangle << endl;
		cout << "Kelectrostatic:                 " << md.Kelec << endl;
		cout << "Kdistance:                      " << md.Kdistance << endl;
		cout << "Kmap:                           " << md.Kmap << endl;
		cout << "Kseparation:                    " << md.Ksep << endl;
		cout << "Kpoint:                         " << md.Kpoint << " (" << md.pointdecay << ")" << endl;
		cout << "location:                       " << md.point << endl;
		cout << "beta constant:                  " << beta << endl;
		cout << "Cutoff distance:                " << md.cutoff << endl;
		cout << "Separation distance:            " << md.sepdist << endl;
		cout << "Maximum shift:                  " << max_shift << endl;
		cout << "Maximum angular increment:      " << max_angle*180.0/M_PI << endl;
		cout << "Steps along link:               " << md.linksteps << endl;
//		cout << "Number of links:                " << nlink << endl;
		if ( md.rigid < 1 ) cout << "Treating the whole ensemble as a rigid body" << endl;
		else if ( md.rigid == 1 ) cout << "Treating each model as a rigid body" << endl;
		else cout << "Modifying each component's coordinates" << endl;
		cout << endl;
		model_calculate_deviations(model, md);
	}
	
	pE = (Efunc)(modcopy, map, md);

	bestE = E = pE;
	Elist[0] = pE;
	
	if ( verbose & VERB_RESULT ) {
		display_energy_headers(md);
		display_energy(0, md, E-pE);
	}
	while ( !done && ( cycle < max_iter ) ) {
		cycle++;
		mpold = models_copy(modcopy);
		if ( md.rigid == 1 ) {		// Rigid models
			if ( select_min ) {
				for ( i=0, Emax=-1e30, mp=modcopy; mp; mp=mp->next, i++ ) {
					if ( Emax < mp->FOM() ) {
						Emax = mp->FOM();
						isel = i;
					}
				}
			} else {
				isel = (long) ((nmod-0.001)*irm*random());
			}
			for ( i=0, mp=modcopy; mp; mp=mp->next, i++ ) {
				if ( i == isel ) mp->select(1);
				else mp->select(0);
			}
		} else if ( md.rigid == 2 ) {
			if ( select_min ) {
				for ( i=0, Emax=-1e30, mp=modcopy; mp; mp=mp->next, i++ ) {
					if ( Emax < mp->FOM() ) {
						Emax = mp->FOM();
						isel = i;
					}
				}
			} else {
				isel = (long) ((nmod-0.001)*irm*random());
			}
			for ( i=0, mp=modcopy; mp; mp=mp->next, i++ ) {
				if ( i == isel ) mp->select(1);
				else mp->select(0);
			}
		}
		(Tfunc)(modcopy, max_angle, max_shift, md.rigid);
		E = (Efunc)(modcopy, map, md);
		Elist[cycle] = E;
//		display_energy(cycle, md, E-pE);
		if ( bestE > E ) {
			bestE = E;
			delete model;
			model = models_copy(modcopy);
			if ( verbose & VERB_RESULT )
				display_energy(cycle, md, E-pE);
		}
		accept = 1;
		if ( E > 1e30 ) accept = 0;		// High energy indicates the model might be outside the box
		if ( E > pE ) {
			r = random()*irm;
			if ( exp(beta*(pE-E)) < r ) accept = 0;
		}
		if ( accept ) {
			pE = E;
			naccept++;
			delete mpold;
		} else {
			delete modcopy;
			modcopy = mpold;
		}
//		model_test_if_within_box(modcopy, md.min, md.max);
	}
	
	delete modcopy;

	E = (Efunc)(model, map, md);

	for ( mp = model; mp; mp = mp->next, ++i )
		mp->select(1);
	
	for ( i=0, Z=Ze=0; i<=cycle; i++ ) {
		Zi = exp(-beta*(Elist[i] - bestE));
		Z += Zi;
		if ( Elist[i] - bestE < 0.001 ) Ze += Zi;
	}
	if ( Z <= 0 ) Z = 1;
	
//	md.Epot = bestE;
		
	if ( verbose & VERB_RESULT ) {
		gloc = models_center_of_coordinates(model) - gloc;
		cout << "Total shift:                        " << gloc << endl;
		for ( i=0, mp = model; mp; mp = mp->next, ++i ) {
			iloc[i] = mp->center_of_coordinates() - iloc[i];
			E = model_map_energy(mp, map, md.Kmap);
			cout << tab << mp->identifier() << tab << iloc[i] << tab << E << endl;
		}
		cout << "Fraction of perturbations accepted: " << naccept*1.0/cycle << endl;
		cout << "Final energy:                       " << bestE << endl;
		cout << "Partition function:                 " << Z/cycle << endl;
		cout << "Probability of final orientation:   " << Ze/Z << endl << endl;
		model_calculate_deviations(model, md);
	}
	
	return model;
}

/**
@brief 	Sets the figure-of-merit for all models to zero.  
@param 	*model		list of models.
@return long			0.
**/
long		model_zero_fom(Bmodel* model)
{
	Bmodel*		mp;
	
	for ( mp = model; mp; mp = mp->next ) mp->FOM(0);
	
	return 0;
}

/**
@brief 	Tests if a model overlaps with a defined box.  
@param 	*model		list of models.
@param 	min			start of box.
@param 	max			end of box.
@return long			number of components selected.
**/
long		model_test_if_within_box(Bmodel* model, Vector3<double> min, Vector3<double> max)
{
	Bmodel*			mp;	
	long			ncomp(0), nsel(0);
	
	if ( verbose ) {
		cout << "Testing if model coordinates fall within the box:" << endl;
		cout << "Box minimum:                    " << min << endl;
		cout << "Box maximum:                    " << max << endl;
		cout << "Model\tMinimum\tMaximum" << endl;
	}
	
    for ( mp = model; mp; mp = mp->next ) {
    	if ( verbose )
    		cout << mp->identifier() << tab << mp->minimum() << tab << mp->maximum() << endl;
    	ncomp += mp->component_count();
    	nsel += mp->select_within_bounds(min, max);
	}
	
	if ( verbose )
		cout << "Fraction within bounding box:   " << nsel*1.0/ncomp << endl;
	
	return nsel;
}

/**
@brief 	Calculates the potential energy for rigid body fitting.
@param 	*model		model structure.
@param 	*map			density map.
@param 	*md				model dynamics structure.
@return double			potential energy.

	The energy is the sum of the overlap, map, and point force energies.

**/
double		monte_rigid_body_fit_energy(Bmodel* model, Bimage* map, Bmodparam& md)
{
	model_zero_fom(model);
	
	md.Epot = 0;
	
	if ( map && md.Kmap ) {
		vector<Bmodel*>		marr = model->array();
#ifdef HAVE_GCD
		dispatch_apply(marr.size(), dispatch_get_global_queue(0, 0), ^(size_t i){
			model_map_energy(marr[i], map, md.Kmap);
		});
#else
#pragma omp parallel for
		for ( long i=0; i<marr.size(); ++i )
			model_map_energy(marr[i], map, md.Kmap);
#endif
		md.Emap = 0;
		for ( long i=0; i<marr.size(); ++i )
			md.Emap += marr[i]->FOM();
		md.Epot += md.Emap;
	}

	if ( md.Ksep )
		md.Epot += models_component_overlap(model, md);

	if ( md.Kdistance )
		md.Epot += models_clashes(model, md);

	if ( md.Kpoint ) {
		md.Epoint = md.Kpoint*md.point.distance(model->center_of_coordinates());
		md.Epot += md.Epoint;
	}
	
	if ( md.Klink ) {
		md.Elink = model_link_energy(model, md.Klink, md.wrap, md.box());
		md.Epot += md.Elink;
	}
	
	return md.Epot;
}

/**
@brief 	Calculates the potential energy for fitting components to a map.
@param 	*model		model structure.
@param 	*map			density map.
@param 	*md				model dynamics structure.
@return double			potential energy.

	The energy is the sum of the link, angle, and map energies.

**/
double		monte_component_fit_energy(Bmodel* model, Bimage* map, Bmodparam& md)
{
	model_zero_forces(model);
	
	Vector3<double>		box(1e10,1e10,1e10);
	md.Elink = model_link_energy(model, md.Klink, md.wrap, box);
	
//	md.Eangle = md_angular_forces(model, md.Kangle, md.wrap);
	
	model_electrostatic_energy(model, md);
	model_distance_energy(model, md);
	
	md.Epoint = 0;
	if ( md.Kpoint )
		md.Epoint = md.Kpoint*md.point.distance(model->center_of_coordinates());
	
	md.Emap = 0;
	if ( map ) md.Emap = model_map_energy(model, map, md.Kmap);
	
	md.Epot = md.Elink + md.Eangle + md.Edistance + md.Eelec + md.Epoint + md.Emap;
	
	return md.Epot;
}

/**
@brief 	Calculates the potential energy for fitting links to a map.
@param 	*model		model structure.
@param 	*map			density map.
@param 	*md				model dynamics structure.
@return double			potential energy.

	The energy is the sum of the link, angle, and map energies.

**/
double		monte_link_fit_energy(Bmodel* model, Bimage* map, Bmodparam& md)
{
	model_zero_forces(model);
	
	Vector3<double>		box(1e10,1e10,1e10);
	md.Elink = model_link_energy(model, md.Klink, md.wrap, box);
	
//	md.Eangle = md_angular_forces(model, md.Kangle, md.wrap);
	
	model_electrostatic_energy(model, md);
	model_distance_energy(model, md);
	
	md.Epoint = 0;
	if ( md.Kpoint )
		md.Epoint = md.Kpoint*md.point.distance(model->center_of_coordinates());
	
	md.Emap = 0;
	if ( map ) md.Emap = model_link_map_energy(model, map, md.Kmap, md.linksteps);
	
	md.Epot = md.Elink + md.Eangle + md.Edistance + md.Eelec + md.Epoint + md.Emap;
	
	return md.Epot;
}


/**
@brief 	Calculates an energy term based on component overlap.
@param 	*model		model structure.
@param 	*md			model dynamics structure.
@return double		total overlap energy.

	The energy is defined as linear decay to the reference separation distance
	and zero beyond:
		Esep = Ksep * (1 - d/dsep)  for  d < dsep, zero otherwise

**/
double		models_component_overlap(Bmodel* model, Bmodparam& md)
{
	if ( md.rigid < 1 || md.Ksep <= 0 ) {
		md.Esep = 0;
		return 0;
	}
	
	long			n(0), i;
	Vector3<double>	min = md.min;
	Vector3<double>	max = md.max;
	Vector3<double>	box = max - min;
	Vector3<double>	sampling(md.sepdist, md.sepdist, md.sepdist);
	Vector3<long>	size = box/sampling + 0.001;
	size = size.max(1);
	sampling = box/size + 0.001;
	long			boxsize = (long) size.volume();
	Vector3<long>	coor;

	Bmodel*			mp;
	Bcomponent*  	comp;

	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	p->sampling(sampling);
	p->origin(-min/sampling);
	float*			data = (float *) p->data_pointer();
	p->next = new Bimage(Float, TSimple, p->size(), p->images());
	float*			dataone = (float *) p->next->data_pointer();
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			coor = p->image->image_coordinates(comp->location());
			if ( p->within_boundaries(coor) ) {
				i = p->index(coor, 0);
				dataone[i] = 1;
			}
			if ( md.rigid > 1 ) {		// Individual components
				for ( i=0; i<boxsize; i++ ) {
					data[i] += dataone[i];
					dataone[i] = 0;
				}
			}
		}
		if ( md.rigid == 1 ) {		// Rigid models
			for ( i=0; i<boxsize; i++ ) {
				data[i] += dataone[i];
				dataone[i] = 0;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG ) {
		Bstring			filename("overlap.map");
		write_img(filename, p, 0);
	}
	
	md.Esep = 0;
	
	for ( i=n=0; i<boxsize; i++ ) {
		if ( data[i] ) {
			n++;
			if ( data[i] > 1 ) md.Esep += data[i] - 1;
		}
	}
	
	md.Esep *= md.Ksep * md.sepdist/n;
	
	delete p;
	
	return md.Esep;
}
	
/**
@brief 	Calculates the distance potential based on clashes between molecules.
@param 	*model			linked list of models.
@param 	*md				model dynamics structure.
@return double			distnace energy.

**/
double		models_clashes(Bmodel* model, Bmodparam& md)
{
	double			d;
	Vector3<double>	box(md.box());
	Bmodel*			mp, *mp2;
	Bcomponent*		comp, *comp2;
	
	for ( mp = model; mp->next; mp = mp->next ) {
		for ( mp2 = mp->next; mp2; mp2 = mp2->next ) {
			for ( comp = mp->comp; comp; comp = comp->next ) {
				for ( comp2 = mp2->comp; comp2; comp2 = comp2->next ) {
					d = comp->location().distance(comp2->location());
					if ( d < md.sepdist )
						md.Edistance += component_distance_potential(comp, comp2, md.Kdistance, md.distancetype, md.wrap, box);
				}
			}
		}
	}
	
	return md.Edistance;
}

/**
@brief 	Randomly transforms models.
@param 	*model			model structure.
@param 	max_angle		maximum rotation angle.
@param 	shift_std		gaussian length for shift vector.
@param	rigid			which part for rigid transformation: 0=all, 1=each model.
@return int				0.

	The transformation is calculated as a random angular rotation and a
	random shift. The shift is sampled from a random vector with a
	gaussian length distribution.
	
	Only one model is rotated if rigid is one.

**/
int			model_rigid_body_transform(Bmodel* model, double max_angle, double shift_std, int rigid)
{
	double			irm = 1.0/get_rand_max();
	Transform		t;

	t.axis = vector3_random(-1, 1);
	t.axis.normalize();	
	t.angle = max_angle*(random()*irm*2 - 1);	
	t.trans = vector3_random_gaussian(0, shift_std);
		
	if ( rigid < 1 ) {
		t.origin = models_center_of_coordinates(model);
		models_rotate(model, t);
	} else if ( rigid == 1 ) {
		for ( ; model; model=model->next ) if ( model->select() ) {
			t.origin = model->center_of_coordinates();
			model_rotate(model, t);
		}
	}
	
	return 0;
}

/**
@brief 	Move components random distances down the energy gradient.
@param 	*model			model structure.
@param 	max_angle		(not used).
@param 	max_shift		maximum shift for each component.
@param	rigid			which part for rigid transformation (here 2=none).
@return double			0.

	The distance of movement is limited to the maximum shift.

**/
int			model_move_components_down_energy(Bmodel* model, double max_angle, double max_shift, int rigid)
{
	Bmodel*			mp;
	Bcomponent*  	comp;
	
	double			irm = 1.0/get_rand_max();
	Vector3<double>	shift;
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = model->comp; comp; comp = comp->next ) {
			shift = comp->force() * (random()*irm);
			shift = shift.min(max_shift);
			comp->location() += shift;
		}
	}
	
	return 0;
}

