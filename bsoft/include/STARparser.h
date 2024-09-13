/**
@file	STARparser.h
@author	Bernard Heymann
@date	20170302 - 20240402
**/

#include <iomanip>

#include "json.h"
#include "string_util.h"
#include "utilities.h"

using namespace std;

#ifndef _STAR_


class STARparser {
private:
	string		s, f;
	void		get_clean_line(ifstream& fstar) {
		if ( fstar.eof() ) return;
		getline(fstar, s);
		s.erase(0, s.find_first_not_of(" \t"));

	}
	string		extract_tag() {
		size_t		i(s.find_first_of(" \t"));
		string		tag(s.substr(1,i-1));
		s = s.substr(i+1);
		s.erase(0, s.find_first_not_of(" \t"));
		return tag;
	}
	string		extract_loop_tag(ifstream& fstar) {
		get_clean_line(fstar);
		size_t		i(s.find_first_of("._", 1));
		return s.substr(1,i-1);
	}
	JSvalue		multiline(ifstream& fstar) {
		string		ml;
		while ( !fstar.eof() ) {
			get_clean_line(fstar);
			if ( s[0] == ';' ) break;
			ml += s;
		}
		return JSvalue(ml);
	}
	JSvalue		value(string vs) {
		long			l;
		double			d;
		if ( isdigit(vs[0]) || vs[0] == '.' || vs[0] == '+' || vs[0] == '-' ) {
			if ( s.find_first_of('.') != string::npos ) {
				try { d = stod(vs); }
				catch (...) { fail("No conversion to a real number."); }
				return JSvalue(d);
			} else {
				try { l = stol(vs); }
				catch (...) { fail("No conversion to an integer."); }
				return JSvalue(l);
			}
		}
		return JSvalue(vs);
	}
	void		fail(string msg) {
		cerr << "STAR parser error: " << msg << endl;
		exit(-1);
	}
public:
	STARparser() {};
	STARparser(string filename) { f = filename; }
	JSvalue		parse() {
		if ( f.length() < 1 ) fail("No filename specified!");
		return parse(f);
	}
	JSvalue		parse(string filename) {
		f = filename;
		cout << "Reading " << f << endl;
		ifstream	fstar(f);
		if ( fstar.fail() ) fail("File " + f + " not opened");
	
		JSvalue					root(JSobject);
		JSvalue					comment(JSarray);
		JSvalue					block(JSarray);

		while ( !fstar.eof() ) {
			get_clean_line(fstar);
			if ( s.find("data_") == 0 ) break;
			if ( s.length() ) comment.push_back(JSvalue(s));
		}

		root["comment"] = comment;
		
		while ( !fstar.eof() ) {
			if ( s.find("data_") == 0 )
				block.push_back(read_block(fstar));
		}
		
//		cout << s << endl;
		cout << "Number of blocks: " << block.size() << endl;
		
		root["data"] = block;
			
		get_clean_line(fstar);

		fstar.close();
	
		cout << "done " << f << endl;
	
		return root;
	}
	JSvalue		read_block(ifstream& fstar) {
		string					tag;
		JSvalue					item(JSobject);

//		cout << s << endl;
		
		get_clean_line(fstar);
	
		while ( !fstar.eof() && s.find("data_") != 0 ) {
			if ( s[0] == '_' ) {
				tag = extract_tag();
				if ( s.length() ) item[tag] = value(s);
			} else if ( s[0] == ';' ) {
				item[tag] = multiline(fstar);
			} else if ( s.find("loop_") == 0 ) {
				tag = extract_loop_tag(fstar);
				item[tag] = read_loop(fstar);
			}
			get_clean_line(fstar);
		}
//		cout << s << endl;

		return item;
	}
	JSvalue		read_loop(ifstream& fstar) {
		JSvalue					loop(JSarray);
		vector<string>			tags, sv;
		JSvalue					item(JSobject);

		while ( !fstar.eof() && s[0] == '_' ) {
			tags.push_back(extract_tag());
			get_clean_line(fstar);
		}

		while ( !fstar.eof() ) {
			sv = split(s);
			if ( sv.size() < tags.size() ) break;
			for ( int i=0; i<tags.size(); ++i )
				item[tags[i]] = value(sv[i]);
			loop.push_back(item);
			item.clear();
			get_clean_line(fstar);
		}
		
		return loop;
	}
};



/**
@brief 	Reads paramaters and data into a JSON data base from a STAR file.
@param	&fstar		reference to output stream.
@param	&loop		reference to loop.

	Every data block is read separately and comments are preserved as far 
	as possible.

	The split flag allows the user to output data blocks in separate files
	in stead of one big file.
	The line length field allows the user to output long lines without
	wrapping it around.
	The comments are ignored but output to the a new file - this can be used
	to document the history of the file.
	The STAR database is a hierarchy consisting of blocks, each with a set
	of items.

**/

int			write_loop(ofstream& fstar, vector<JSvalue>& loop)
{
	int					err(0);
//	map<string, JSvalue>	item(loop[0].object());
	JSvalue&			item(loop[0]);
	
	fstar << endl << "loop_" << endl;
	
	for ( auto it = item.object_begin(); it != item.object_end(); ++it )
		fstar << "_" << it->first << endl;
	
	for ( auto il = loop.begin(); il != loop.end(); ++il ) {
		for ( auto it = il->object_begin(); it != il->object_end(); ++it )
			fstar << it->second << " ";
		fstar << endl;
	}
	
	fstar << endl;
	
	return err;
}

int			write_block(ofstream& fstar, JSvalue& block)
{
	int					err(0), line_length(80);
	
	fstar << endl << "data_" << endl << endl;

	for ( auto it = block.object_begin(); it != block.object_end(); ++it ) 
			if ( it->second.type() != JSarray ) {
		fstar << left << "_" << setw(38) << it->first;
		if ( it->second.type() != JSstring ) {
			fstar << " " << it->second << endl;
		} else if ( it->second.size() < line_length ) {
			fstar << " " << it->second << endl;
		} else {
			fstar << endl << ";" << endl;
			for ( size_t i = 0; i<it->second.value().length(); i += line_length )
				fstar << it->second.value().substr(i, line_length) << endl;
			fstar << ";" << endl;
		}
	}

	for ( auto it = block.object_begin(); it != block.object_end(); ++it ) 
		if ( it->second.type() == JSarray ) err += write_loop(fstar, it->second.array());

	return err;
}


int			write_block(string filename, JSvalue& block)
{
	int				err(0);
	
	cout << "writing " << filename << endl;

	string			comment = "# Written by Bsoft2\n";
	
	ofstream		fstar(filename.c_str());
	if ( fstar.fail() ) {
		cerr << "Error: Not able to write " << filename << endl;
		return -1;
	}

	fstar << comment << endl;

	write_block(fstar, block);
	
	fstar.close();
	
	return err;
}

/**
@brief 	Writes a JSON data base to a STAR format file.
@param	filename	file name.
@param	root		JSON database
@return int			error code (<0 means failure).

**/
int			writeSTAR(string filename, JSvalue& root)
{
	int				err(0);
	
	cout << "writing " << filename << endl;
	
	JSvalue&			comment(root["comment"]);
	JSvalue&			block(root["data"]);

	ofstream		fstar(filename.c_str());
	if ( fstar.fail() ) {
		cerr << "Error: Not able to write " << filename << endl;
		return -1;
	}

	for ( auto it = comment.begin(); it != comment.end(); ++it )
		fstar << it->value() << endl;
	fstar << "# " << filename << " written by Bsoft" << endl;

	for ( auto ib = block.begin(); ib != block.end(); ++ib )
		err += write_block(fstar, *ib);
	
	fstar.close();
	
	return err;
}

#define _STAR_
#endif


