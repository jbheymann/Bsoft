/**
@file	model_water.cpp
@brief	Generating and managing water
@author 	Bernard Heymann
@date	Created: 20001014
@date	Modified: 20230704
**/

#include "Bmodel.h"
//#include "mol_bonds.h"
#include "model_util.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates one water molecule at a given location.
@param 	&watername		molecule name.
@param 	Ocoord			oxygen coordinates.
@return Bmodel*			new water molecule.
**/
Bmodel*		model_generate_one_water(string& watername, Vector3<double> Ocoord)
{
	double			rand_max = 1.0*get_rand_max();	
	double			OHbond(0.9572);
	double			HOHangle(109.47*M_PI/180.0);
	string			restype("HOH");
	Vector3<double>	H1vec, H2vec;
	Matrix3			mat(1);
	
	Bmodel*			water = new Bmodel(watername);
	Bcomptype*		cto = water->add_type("O");
	Bcomptype*		cth = water->add_type("H");
	cto->index(7);
	cth->index(0);
	cto->mass(16);
	cth->mass(1);
	
	Bcomponent*		O = water->add_component(1);
	O->type(cto);
	O->description("O");
	O->add_description("O");
	O->add_description(restype);
	O->charge(-0.834);
	O->location() = Ocoord;
	
	Bcomponent*		H1 = water->add_component(2);
	H1->type(cth);
	H1->description("H");
	H1->add_description("HO1");
	H1->add_description(restype);
	H1->charge(0.417);
	H1vec[0] = random()/rand_max - 0.5;
	H1vec[1] = random()/rand_max - 0.5;
	H1vec[2] = random()/rand_max - 0.5;
	H1vec.normalize();
	H1->location() = (H1vec * OHbond) + Ocoord;
	water->comp->find_and_add_links(O->identifier(), H1->identifier());
	
	Bcomponent*		H2 = water->add_component(3);
	H2->type(cth);
	H2->description("H");
	H2->add_description("HO2");
	H2->add_description(restype);
	H2->charge(0.417);
	H2vec[0] = random()/rand_max - 0.5;
	H2vec[1] = random()/rand_max - 0.5;
	H2vec[2] = random()/rand_max - 0.5;
	H2vec.normalize();

	H2vec = H1vec.cross(H2vec);
	mat = Matrix3(H2vec, HOHangle);
	H2vec = mat * H1vec;
	H2->location() = (H2vec * OHbond) + Ocoord;
	water->comp->find_and_add_links(O->identifier(), H2->identifier());

	return water;
}

/**
@brief 	Generates a block of water based on a regular lattice.
@param 	size			size of block.
@param 	type			type of lattice, 2=rectangular, 3=tetrahedral.
@return Bmodel*	new list of water molecules.


	The number of water molecules generated is calculated as:
		n = volume * 0.03346.

**/
Bmodel*	model_generate_regular_water(Vector3<double> size, int type)
{
	Bmodel*			waters = NULL;

	if ( size.volume() < 1 ) {
		cerr << "Error: A box size must be specified!" << endl;
		return waters;
	}
		
	Vector3<double>	Ocoord, start;
	Vector3<double>	interval(3.1, 3.1, 3.1);			// Interval between water molecules in each direction
	if ( type == 3 ) {
		interval[0] = 3.48;
		interval[1] = interval[0]*sqrt(3.0)/2.0;
		interval[2] = interval[0]*sqrt(2.0/3.0);
	}
	
	size = interval*(size/interval) + 0.001;

	double			volume(size.volume());
	long			n(volume*0.03346);
	char			c('A');			// Character for molecule ID
	string			watername("A");
	double			x, y, z;
	Bmodel*			w1 = NULL;

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating regular water" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Volume:                         " << size << " = " << volume << " A3" << endl;
		cout << "Interval between molecules:     " << interval << " A" << endl;
		cout << "Estimated number of waters:     " << n << endl << endl;
	}
	
	n = 0;
	for ( z=interval[2]/2; z<size[2]; z+=interval[2] ) {
		if ( type == 3 ) start[1] = 1 - start[1];
		start[0] = start[1];
		for ( y=interval[1]/2+start[1]*interval[0]/(2*sqrt(3.0)); y<size[1]; y+=interval[1] ) {
			if ( type == 3 ) start[0] = 1 - start[0];
			for ( x=(2-start[0])*interval[0]/2; x<size[0]; x+=interval[0] ) {
				watername[0] = c;
				Ocoord = Vector3<double>(x, y, z);
				if ( w1 ) w1 = w1->next = model_generate_one_water(watername, Ocoord);
				else waters = w1 = model_generate_one_water(watername, Ocoord);
				n++;
				c++;
				if ( c > 90 ) c = 65;
			}
		}
	}

	vector<Vector3<double>>		bounds = models_calculate_bounds(waters);
	size = bounds[1] - bounds[0];
	volume = size.volume();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of waters:               " << n << endl;
		cout << "Volume:                         " << size << " = " << volume << " A3" << endl;
		cout << "Density:                        " << n*1.0/volume << 
			" molecules/A3 (" << n/(volume*0.03346) << ")" << endl << endl;
	}
	
	return waters;
}

/**
@brief 	Generates a block of water with random placement.
@param 	size			size of block.
@return Bmodel*	new list of water molecules.

	The number of water molecules generated is calculated as:
		n = volume * 0.03346.

**/
Bmodel*		model_generate_random_water(Vector3<double> size)
{
	Bmodel*	waters = NULL;
	
	if ( size.volume() < 1 ) {
		cerr << "Error: A box size must be specified!" << endl;
		return waters;
	}
	
	random_seed();
	
	int				i;
	char			c('A');		// Character for molecule ID
	string			watername("A");
	double			rand_max = 1.0*get_rand_max();
	
	double			volume(size.volume());
	long			n(volume*0.03346);
	Vector3<double>	Ocoord;
	Bmodel*			w1 = NULL;
	Bcomponent*		comp;
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating random water" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Volume:                         " << size << " = " << volume << " A3" << endl;
		cout << "Number of waters:               " << n << endl << endl;
	}
	
	for ( i=0; i<n; i++ ) {
		watername[0] = c;
		Ocoord[0] = size[0]*random()/rand_max;
		Ocoord[1] = size[1]*random()/rand_max;
		Ocoord[2] = size[2]*random()/rand_max;
		if ( w1 ) w1 = w1->next = model_generate_one_water(watername, Ocoord);
		else waters = w1 = model_generate_one_water(watername, Ocoord);
		for ( comp = w1->comp; comp; comp = comp->next )
			comp->location() = vector3_set_PBC(comp->location(), size);
		c++;
		if ( c > 90 ) c = 65;
	}
	
	vector<Vector3<double>>		bounds = models_calculate_bounds(waters);
	size = bounds[1] - bounds[0];
	volume = size.volume();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of waters:               " << n << endl;
		cout << "Volume:                         " << size << " = " << volume << " A3" << endl;
		cout << "Density:                        " << n*1.0/volume << 
			" molecules/A3 (" << n/(volume*0.03346) << ")" << endl << endl;
	}
	
	return waters;
}


/*
@brief 	Generates a bond angle list for a block of waters.
@param 	*molgroup	molecule group.
@return Bangle*				new bond angle list.
*/
/*Bangle*		water_angle_list(Bmolgroup* molgroup)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Bangle*			angle = NULL;
	double			OHangle = 109.47*M_PI/180.0;
	int				nangle = 0;

	molgroup->angle = NULL;
	
    for ( nangle = 0, mol = molgroup->mol; mol; mol = mol->next ) {
		res = mol->res;
		atom = res->atom;
		angle = angle_add(&angle, atom->next, atom, atom->next->next, OHangle, 1);
		if ( !molgroup->angle ) molgroup->angle = angle;
		nangle ++;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of angles generated:     " << nangle << endl;
	
	return molgroup->angle;
}
*/

/**
@brief 	Calculates a radial distribution function for water molecules.
@param 	*waters			list of water molecules.
@param 	interval		interval between bins.
@param 	cutoff			distance cutoff.
@return int				0.
**/
int			model_calc_water_rdf(Bmodel* waters, double interval, double cutoff)
{
	vector<Vector3<double>>		bounds = models_calculate_bounds(waters);
	Vector3<double>	box(bounds[1]-bounds[0]);

	int				i, ii, j, t, x, y, z, xx, yy, zz, ix, iy, iz;
	double			dist;
	Vector3<double>	d, sam(cutoff, cutoff, cutoff), ori;
	Vector3<long>	size((long) (box[0]/sam[0] + 0.001),
		(long) (box[1]/sam[1] + 0.001), (long) (box[2]/sam[2] + 0.001));
	size = size.max(1);
	for ( i=0; i<3; i++ ) sam[i] = box[i]/size[i] + 0.001;

	vector<vector<Bcomponent*>>	grid = model_component_grid(waters, size, ori, sam);

	double			mult(1.0/interval);
	long			n(mult*cutoff);
	vector<long>	rdf(n*3, 0);
	Bcomptype*		ct1;
	Bcomptype*		ct2;
	
	for ( i=z=0; z<size[2]; z++ ) {
		for ( y=0; y<size[1]; y++ ) {
			for ( x=0; x<size[0]; x++, i++ ) {
				for ( auto& c1: grid[i] ) {
					ct1 = c1->type();
					for ( zz=z-1; zz<=z+1; zz++ ) {
						iz = zz;
						if ( iz < 0 ) iz += size[2];
						if ( iz >= size[2]) iz -= size[2];
						for ( yy=y-1; yy<=y+1; yy++ ) {
							iy = yy;
							if ( iy < 0 ) iy += size[1];
							if ( iy >= size[1] ) iy -= size[1];
							for ( xx=x-1; xx<=x+1; xx++ ) {
								ix = xx;
								if ( ix < 0 ) ix += size[0];
								if ( ix >= size[0] ) ix -= size[0];
								ii = (iz*size[1] + iy)*size[0] + ix;
								for ( auto& c2: grid[ii] )
										if ( c2 != c1 ) {
									ct2 = c2->type();
									d = vector3_difference_PBC(c2->location(), c1->location(), box);
									dist = d.length();
									if ( dist < cutoff ) {
										j = (int) (mult*dist + 0.5);
										if ( j < n ) {
											if ( ct1->identifier()[0] == 'H' && ct2->identifier()[0] == 'H' ) t = 0;
											else if ( ct1->identifier()[0] == 'H' || ct2->identifier()[0] == 'H' ) t = 1;
											else t = 2;
											rdf[3*j+t] += 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	grid.clear();

	long		nwater = model_component_count(waters);
			
	cout << "Calculating the radial distribution function:" << endl;
	
	double		vol_per_water = 180/6.022;	// Volume occupied by water = 29.9 A^3
//	double		vol_within_cutoff = (4.0/3.0)*M_PI*cutoff*cutoff*cutoff;
//	double		water_within_cutoff = vol_within_cutoff/vol_per_water;
	long		ndist[3] = {0,0,0};
	double		rdftot[3] = {0,0,0};
	double		norm[3], vol_surface;
	double		vol_surf_fac = (4.0/3.0)*M_PI*interval*interval*interval;
	cout << endl << "r(A)\tH-H\tRDF(HH)\tH-O\tRDF(HO)\tO-O\tRDF(OO)" << endl;
	for ( i=1; i<n; i++ ) {
		cout << i*interval;
		vol_surface = vol_surf_fac*(3*i*i + 0.25);
		norm[2] = vol_per_water/(vol_surface*2*nwater);		// Normalizer for RDF(OO)
		norm[0] = norm[1] = 0.25*norm[2];					// Volume per H is half and 2 times more H's
		for ( t=0; t<3; t++ ) {
			cout << tab << rdf[3*i+t] << tab << rdf[3*i+t]*norm[t];
			ndist[t] += rdf[3*i+t];
			rdftot[t] += rdf[3*i+t]*norm[t];
		}
		cout << endl;
	}
	cout << "Number of distances included = " << ndist[0] << " " << ndist[1] << " " << ndist[2] << endl;
	cout << "RDF sum = " << rdftot[0] << " " << rdftot[1] << " " << rdftot[2] << endl;
	cout << "RDF avg = " << rdftot[0]/n << " " << rdftot[1]/n << " " << rdftot[2]/n << endl << endl;
	
	return 0;
}

