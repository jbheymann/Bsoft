/**
@file	rwmodel_vega.cpp
@brief	Library routines to read and write Vega model parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 20221115
**/

#include "rwmodel.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Vega model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_vega(vector<string> file_list)
{
	bool			readflag(0);
	int				n, m(0), l1, l2, l3;
	float			x, y, z;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	string			s, id;
	string			comptype("VER");
	ifstream		fmod;

	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		fmod.open(filename.c_str());
		if ( fmod.fail() ) return  NULL;
		readflag = 0;
		m = 0;
		while ( !fmod.eof() ) {
			getline(fmod, s);
			if ( readflag ) {
//				sscanf(aline, "%d %f %f %f %d %d %d", &n, &x, &y, &z, &l1, &l2, &l3);
    			istringstream	ss(s);
    			ss >> n >> x >> y >> z >> l1 >> l2 >> l3;
				if ( n > 0 ) {
					id = to_string(n);
					if ( comp ) comp = comp->add(n);
					else model->comp = comp = new Bcomponent(n);
					comp->location(Vector3<float>(x, y, z));
					comp->flag[0] = l1;	// The flags carry the indices of the linked components
					comp->flag[1] = l2;
					comp->flag[2] = l3;
					comp->select(n);
					comp->type(mp->add_type(comptype));
					if ( m < n ) m = n;
				} else readflag = 0;
			}
			if ( s.find("writegraph3d planar") != string::npos ) {
				link = NULL;
				comp = NULL;
				if ( mp ) mp = mp->add(base(filename));
				else model = mp = new Bmodel(base(filename));
				mp->model_type(mp->identifier());
				readflag = 1;
			}
		}
		fmod.close();
		m++;
		vector<Bcomponent*>		car(m);
		for ( comp = mp->comp; comp; comp = comp->next ) car[comp->select()] = comp;
		for ( comp = mp->comp; comp; comp = comp->next ) {
			for ( n=0; n<comp->link.size() && comp->flag[n]>0; n++ ) {
				comp->link[n] = car[comp->flag[n]];
				if ( comp->select()< comp->flag[n] ) {
					link = link_add(&link, comp, comp->link[n], 0, 1);
					if ( !mp->link ) mp->link = link;
				}
				comp->flag[n] = 1;	// flag=1 indicates a link
			}
			comp->select(1);
		}
	}

	return model;
}

/**
@brief 	Writes Vega model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int			models written.
**/
int			write_model_vega(string& filename, Bmodel* model, int splt)
{
	int				n;
	Bmodel*			mp = NULL;
	Bcomponent*		comp;
	string			onename;

	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
//			onename = filename.pre_rev('.') + string(n+1, format) + filename.post_rev('.');
			onename = insert(filename, n+1, splt);
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << ">>writegraph3d planar <<" << endl;
		for ( comp = mp->comp; comp; comp = comp->next )
			fmod << stol(comp->identifier()) << " " <<
				comp->location() << " " <<
				stol(comp->link[0]->identifier()) << " " <<
				stol(comp->link[1]->identifier()) << " " <<
				stol(comp->link[2]->identifier()) << endl;
		fmod << "0" << endl;
		fmod.close();
	}
	
	return 0;
}
