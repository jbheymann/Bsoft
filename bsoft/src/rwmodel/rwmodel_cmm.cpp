/**
@file	rwmodel_cmm.cpp
@brief	Library routines to read and write Chimera marker model parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 20230621
**/

#include "rwmodel.h"
//#include "file_util.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

vector<pair<string,string>>	tag_value_parse(string s)
{
//	cout << s << endl;
	vector<pair<string,string>>	tv;
	
	size_t		i = s.find("<") + 1;
	s = s.substr(i, s.find(">") - i);
	
	vector<string>	slist = split(s);
	string			tag, value;
	
	for ( auto sp: slist ) {
		i = sp.find("=");
//		cout << "sp = " << sp << endl;
		if ( i != string::npos ) {
			tag = sp.substr(0, i);
			if ( sp[++i] == '\"' ) i++;
			value = sp.substr(i, sp.rfind("\"") - i);
		} else {
			tag = sp;
			value = "";
		}
		tv.push_back(make_pair(tag, value));
//		cout << "tv: " << tag << tab << value << endl;
	}
	
	return tv;
}


/**
@brief 	Reads Chimera marker model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_chimera(vector<string> file_list)
{
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp1 = NULL;
	Bcomponent*		comp2 = NULL;
	Blink*			link = NULL;
	string			s, id, path;
	ifstream		fmod;
	RGBA<float>		rgba(1,1,1,1);	// Default white
	vector<pair<string,string>>	tvlist;

	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		fmod.open(filename.c_str());
		if ( fmod.fail() ) return NULL;
		path = filename.substr(filename.rfind("/")-1);
		while ( !fmod.eof() ) {
			getline(fmod, s);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG read_cmm: " << s << endl;
			tvlist = tag_value_parse(s);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG read_cmm: tag=" << tvlist[0].first << endl;
			if ( tvlist[0].first == "marker_set" ) {
				link = NULL;
				comp = NULL;
				if ( model ) mp = mp->add("1");
				else model = mp = new Bmodel("1");
				for ( auto tv: tvlist ) if ( tv.second.length() ) {
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG read_cmm: tag=" << tv.first << " value=" << tv.second << endl;
					if ( tv.first == "name" ) mp->identifier(tv.second);
					if ( tv.first == "type" ) mp->model_type(tv.second);
					if ( tv.first == "hand" ) mp->handedness(to_integer(tv.second));
					if ( tv.first == "symmetry" ) mp->symmetry(tv.second);
					if ( tv.first == "file" ) mp->mapfile(tv.second);
					if ( tv.first == "img_num" ) mp->image_number(to_integer(tv.second));
					if ( tv.first == "fom" ) mp->FOM(to_real(tv.second));
					if ( tv.first == "select" ) mp->select(to_integer(tv.second));
				}
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG read_cmm: model id = " << mp->identifier() << endl;
			} else if ( tvlist[0].first == "/marker_set" ) {
				break;
			} else if ( tvlist[0].first == "comment" ) {
				while ( getline(fmod, s) && s.find("</comment>") == string::npos )
					mp->comment(mp->comment() + s + "\n");
			} else if ( tvlist[0].first == "type" ) {
				for ( auto tv: tvlist )
					if ( tv.first == "id" ) id = tv.second;
				ct = mp->add_type(id);
				if ( !mp->type ) mp->type = ct;
				for ( auto tv: tvlist ) {
					if ( tv.first == "file" ) ct->file_name(tv.second);
					if ( tv.first == "num" ) ct->image_number(to_integer(tv.second));
					if ( tv.first == "mass" ) ct->mass(to_real(tv.second));
					if ( tv.first == "fom" ) ct->FOM(to_real(tv.second));
					if ( tv.first == "select" ) ct->select(to_integer(tv.second));
				}
			} else if ( tvlist[0].first == "marker" ) {
				for ( auto tv: tvlist )
					if ( tv.first == "id" ) id = tv.second;
				if ( comp ) comp = comp->add(id);
				else comp = mp->add_component(id);
				for ( auto tv: tvlist ) {
					if ( tv.first == "type" )
						comp->type(mp->add_type(tv.second));
					if ( tv.first == "x" ) comp->location()[0] = to_real(tv.second);
					if ( tv.first == "y" ) comp->location()[1] = to_real(tv.second);
					if ( tv.first == "z" ) comp->location()[2] = to_real(tv.second);
					if ( tv.first == "vx" ) comp->view()[0] = to_real(tv.second);
					if ( tv.first == "vy" ) comp->view()[1] = to_real(tv.second);
					if ( tv.first == "vz" ) comp->view()[2] = to_real(tv.second);
					if ( tv.first == "va" ) comp->view()[3] = to_real(tv.second)*M_PI/180.0;
					if ( tv.first == "radius" ) comp->radius(to_real(tv.second));
					if ( tv.first == "r" ) comp->color()[0] = to_real(tv.second);
					if ( tv.first == "g" ) comp->color()[1] = to_real(tv.second);
					if ( tv.first == "b" ) comp->color()[2] = to_real(tv.second);
					if ( tv.first == "a" ) comp->color()[3] = to_real(tv.second);
					if ( tv.first == "density" ) comp->density(to_real(tv.second));
					if ( tv.first == "fom" ) comp->FOM(to_real(tv.second));
					if ( tv.first == "select" ) comp->select(to_integer(tv.second));
				}
			} else if ( tvlist[0].first == "link" ) {
				for ( auto tv: tvlist ) {
					if ( tv.first == "id1" )
						for ( comp = mp->comp; comp; comp = comp->next )
							if ( comp->identifier() == tv.second ) comp1 = comp;
					if ( tv.first == "id2" )
						for ( comp = mp->comp; comp; comp = comp->next )
							if ( comp->identifier() == tv.second ) comp2 = comp;
				}
				link = link_add(&link, comp1, comp2, 0, 1);
				if ( !mp->link ) mp->link = link;
				for ( auto tv: tvlist ) {
					if ( tv.first == "radius" ) link->radius(to_real(tv.second));
					if ( tv.first == "r" ) rgba[0] = to_real(tv.second);
					if ( tv.first == "g" ) rgba[1] = to_real(tv.second);
					if ( tv.first == "b" ) rgba[2] = to_real(tv.second);
					if ( tv.first == "a" ) rgba[3] = to_real(tv.second);
					if ( tv.first == "fom" ) link->FOM(to_real(tv.second));
					if ( tv.first == "select" ) link->select(to_integer(tv.second));
				}
				link->color(rgba);
			}
		}
		fmod.close();
	}

	return model;
}

/**
@brief 	Writes Chimera marker model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int			models written.
**/
int			write_model_chimera(string& filename, Bmodel* model, int splt)
{
	int				n;
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	string			onename;
//	char			format[32];

//	snprintf(format, 32, "_%%0%dd.", split);

	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
//			onename = filename.pre_rev('.') + string(n+1, format) + filename.post_rev('.');
			onename = insert(filename, n+1, splt);
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << "<?xml version=\"1.0\"?>" << endl;
		fmod << "<marker_set name=\"" << mp->identifier() << "\" type=\"" << mp->model_type()
			<< "\" hand=\"" << mp->handedness() << "\" symmetry=\"" << mp->symmetry()
			<< "\" file=\"" << mp->mapfile() << "\" img_num=\"" << mp->image_number()
			<< "\" fom=\"" << mp->FOM()<< "\" select=\"" << mp->select()<< "\">" << endl;
		fmod << "<comment>" << endl << mp->comment() << endl << "</comment>" << endl;
		for ( ct = mp->type; ct; ct = ct->next ) {
			fmod << "<type id=\"" << ct->identifier() << "\" file=\"" << ct->file_name()
				<< "\" num=\"" << ct->image_number() << "\" mass=\"" << ct->mass()
				<< "\" fom=\"" << ct->FOM()<< "\" select=\"" << ct->select()<< "\"/>" << endl;
		}
		for ( comp = mp->comp; comp; comp = comp->next ) {
			fmod << "<marker id=\"" << comp->identifier() << "\" type=\"" << comp->type()->identifier()
				<< "\" x=\"" << comp->location()[0] << "\" y=\"" << comp->location()[1]
				<< "\" z=\"" << comp->location()[2] << "\" vx=\"" << comp->view()[0]
				<< "\" vy=\"" << comp->view()[1] << "\" vz=\"" << comp->view()[2]
				<< "\" va=\"" << comp->view().angle()*180.0/M_PI
				<< "\" r=\"" << comp->color()[0] << "\" g=\"" << comp->color()[1]
				<< "\" b=\"" << comp->color()[2] << "\" a=\"" << comp->color()[3]
				<< "\" radius=\"" << comp->radius() << "\" density=\"" << comp->density()
				<< "\" fom=\"" << comp->FOM()<< "\" select=\"" << comp->select()<< "\"/>" << endl;
		}
		for ( link = mp->link; link; link = link->next ) {
			fmod << "<link id1=\"" << link->comp[0]->identifier() << "\" id2=\"" << link->comp[1]->identifier()
				<< "\" r=\"" << link->color()[0] << "\" g=\"" << link->color()[1]
				<< "\" b=\"" << link->color()[2] << "\" a=\"" << link->color()[3]
				<< "\" radius=\"" << link->radius() << "\" fom=\"" << link->FOM()
				<< "\" select=\"" << link->select()<< "\"/>" << endl;
		}
		fmod << "</marker_set>" << endl;
		fmod.close();
	}
	
	return  n;
}


