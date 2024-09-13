/**
@file	rwmodel_bild.h
@brief	Header file for reading and writing Chimera BILD files
@author 	Bernard Heymann
@date	Created: 20140706
@date	Modified: 20221115
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_bild(vector<string> file_list);
int			write_model_bild(string filename, Bmodel* model, int split);
int			model_to_bild_orientations(string filename, Bmodel* model, int vec_type, int color_type);
int			model_to_bild_view_sphere(string filename, Bmodel* model, int color_type);
int			model_to_bild_force_vectors(string filename, Bmodel* model, int color_type);
int			model_to_bild_view_polygons(string filename, Bmodel* model, int order, int color_type);
int			model_to_bild_polygons(string filename, Bmodel* model, int color_type);
int			model_to_bild_neighbor_planes(string filename, Bmodel* model, int color_type);


