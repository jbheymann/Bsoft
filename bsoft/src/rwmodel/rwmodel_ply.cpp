/**
@file	rwmodel_ply.cpp
@brief	Library routines to read and write PLY model parameters
@author 	Bernard Heymann
@date	Created: 20250417
@date	Modified: 20250417
**/

#include "rwmodel.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
	The properties are encoded in a map with the property name first and its type second
*/
class PlyElement {				// description of an element
public:
	string				name;	// element name
	int					num;	// number of elements in this object
	int					size;	// size of element (bytes) or -1 if variable
	int					nprops;	// number of properties for this element
	int					list;	// if properties are stored as a list
	map<string,string>	props;	// list of properties in the file
	void		clear() {
		name = "";
		num = 0;
		size = 0;
		nprops = 0;
		list = 0;
		props.clear();
	}
	void		show () {
		cout << "element " << name << tab << num << endl;
		for ( auto pr: props ) {
			if ( list )			
				cout << "property list uchar " << pr.second << tab << pr.first << endl;
			else
				cout << "property " << pr.second << tab << pr.first << endl;
		}
	}
} ;


/**
@brief 	Reads PLY model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_ply(vector<string> file_list)
{
	int					i, j, n(0), m(0);
	Bmodel*				model = NULL;
	Bmodel*				mp = NULL;
	Bcomponent*			comp = NULL;
	Bcomponent*			comp2 = NULL;
	Blink*				link = NULL;
	Bpolygon*			poly = NULL;
	string				s, id, v1, v2;
	string				comptype("VER");
	ifstream			fmod;
	PlyElement			el;
	vector<PlyElement>	elements;
	Vector3<double>		loc;
	RGBA<double>		rgba(255,255,255,255);

	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		fmod.open(filename.c_str());
		if ( fmod.fail() ) return  NULL;
		getline(fmod, s);
		if ( s.find("ply") != 0 ) {
			cerr << "Error: The file is not PLY format!" << endl;
			return NULL;
		}
		getline(fmod, s);
		if ( s.find("ascii") == string::npos ) {
			cerr << "Error: The file is not ASCII format!" << endl;
			return NULL;
		}
		
		// Header
		while ( !fmod.eof() ) {
			getline(fmod, s);
			vector<string> vs = split(s);
			if ( vs[0] == "comment" ) {
			} else if ( vs[0] == "element" ) {
				if ( el.name.size() > 1 ) {
					elements.push_back(el);
					el.clear();
				}
				el.name = vs[1];
				el.num = to_integer(vs[2]);
				el.list = 0;
				el.nprops = 0;
				el.size = 0;	
			} else if ( vs[0] == "property" ) {
				if ( vs[1] == "list" ) {
					el.props[vs[4]] = vs[3];
					el.list = 1;
				} else {
					el.props[vs[2]] = vs[1];
				}
				el.nprops++;
			} else if ( vs[0] == "end_header" ) {
				if ( el.name.size() > 1 ) {
					elements.push_back(el);
					el.clear();
				}
				break;
			}
		}
		
//		for ( auto el: elements ) el.show();

		// Content transferred to model
		if ( mp ) mp = mp->add(base(filename));
		else model = mp = new Bmodel(base(filename));
		mp->model_type(mp->identifier());
		comp = NULL;
		link = NULL;
		poly = NULL;
		m = n = 0;
		while ( !fmod.eof() && m < elements.size() ) {
			getline(fmod, s);
			vector<string> vs = split(s);
			if ( elements[m].name == "vertex" ) {
				i = 0;
				for ( auto p: elements[m].props ) {
					if ( p.first == "x" ) loc[0] = to_real(vs[i]);
					else if ( p.first == "y" ) loc[1] = to_real(vs[i]);
					else if ( p.first == "z" ) loc[2] = to_real(vs[i]);
					else if ( p.first == "red" ) rgba[0] = to_real(vs[i]);
					else if ( p.first == "green" ) rgba[1] = to_real(vs[i]);
					else if ( p.first == "blue" ) rgba[2] = to_real(vs[i]);
					i++;
				}
				id = to_string(n);
				if ( comp ) comp = comp->add(n);
				else mp->comp = comp = new Bcomponent(n);
				comp->location(loc);
				comp->color(rgba/255);
				comp->FOM(1);
				comp->select(1);
				comp->type(mp->add_type(comptype));
				n++;
			} else if ( elements[m].name == "edge" ) {
				i = 0;
				for ( auto p: elements[m].props ) {
					if ( p.first == "vertex1" ) v1 = vs[i];
					else if ( p.first == "vertex2" ) v2 = vs[i];
					else if ( p.first == "red" ) rgba[0] = to_real(vs[i]);
					else if ( p.first == "green" ) rgba[1] = to_real(vs[i]);
					else if ( p.first == "blue" ) rgba[2] = to_real(vs[i]);
					i++;
				}
				comp = mp->comp->find(v1);
				comp2 = mp->comp->find(v2);
				if ( link ) link = link->add(comp, comp2, 0, 1);
				else mp->link = link = new Blink(comp, comp2, 0, 1);
				link->color(rgba/255);
				link->FOM(1);
				link->select(1);
				n++;
			} else if ( elements[m].name == "face" ) {
				i = 0;
				if ( elements[m].list ) {
					if ( poly ) poly = poly->add();
					else mp->poly = poly = new Bpolygon();
					j = to_integer(vs[0]);
					for ( i=1; i<=j; ++i ) {
						comp = mp->comp->find(vs[i]);
						poly->comp.push_back(comp);
					}
					poly->closed(1);
					poly->FOM(1);
					poly->select(1);
				}
				n++;
			}
//			cout << m << tab << elements[m].num << tab << n << tab << s << endl;
			if ( n >= elements[m].num ) {
				m++;
				n = 0;
			}
		}
	}

	return model;
}

/**
@brief 	Writes PLY model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int			models written.
**/
int			write_model_ply(string& filename, Bmodel* model, int splt)
{
	long			n, ncomp(0), nlink(0), npoly(0);
	Bmodel*			mp = NULL;
	Bcomponent*		comp;
	Blink*			link;
	Bpolygon*		poly;
	string			onename;
	RGBA<short>		rgba(255,255,255,255);

	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
			onename = insert(filename, n+1, splt);
		else
			onename = filename;
		ncomp = mp->component_count();
		nlink = mp->link_count();
		npoly = mp->polygon_count();
		vector<string> 	com = split(model->comment(), '\n');
		
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		
		fmod << "ply" << endl;
		fmod << "format ascii 1.0" << endl;
		
		for ( auto c: com )
			fmod << "comment " << c << endl;
		
		fmod << "element vertex " << ncomp << endl;
		fmod << "property float x" << endl;
		fmod << "property float y" << endl;
		fmod << "property float z" << endl;
		fmod << "property uchar red" << endl;
		fmod << "property uchar green" << endl;
		fmod << "property uchar blue" << endl;
		
		if ( nlink ) {
			fmod << "element edge " << nlink << endl;
			fmod << "property int vertex1" << endl;
			fmod << "property int vertex2" << endl;
			fmod << "property uchar red" << endl;
			fmod << "property uchar green" << endl;
			fmod << "property uchar blue" << endl;
		}
		
		if ( npoly ) {
			fmod << "element face " << npoly << endl;
			fmod << "property list uchar int vertex_index" << endl;
		}
		
		fmod << "end_header" << endl;
		
		for ( comp = mp->comp; comp; comp = comp->next ) {
			rgba = comp->color()*255;
			fmod << comp->location()[0] << " "
				<< comp->location()[1] << " "
				<< comp->location()[2] << " "
				<< rgba[0] << " "
				<< rgba[1] << " "
				<< rgba[2] << endl;
		}
		
		for ( link = mp->link; link; link = link->next ) {
			rgba = comp->color()*255;
			fmod << link->comp[0]->identifier() << " "
				<< link->comp[1]->identifier() << " "
				<< rgba[0] << " "
				<< rgba[1] << " "
				<< rgba[2] << endl;
		}
		
		for ( poly = mp->poly; poly; poly = poly->next ) {
			fmod << poly->comp.size();
			for ( auto& c: poly->comp )
				fmod << " " << c->identifier();
			fmod << endl;
		}
		
		fmod.close();
	}
	
	return 0;
}
