/**
@file	model_select.h
@brief	Header file for reading and writing atomic model files
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250623
**/

#include "Bmodel.h"
#include "Bimage.h"

/* Function prototypes */
long		models_selection_stats(Bmodel* model);
long		models_count_selected(Bmodel* model);
long 		models_show_selection(Bmodel* model);
long		models_show_component_types(Bmodel* model);
long 		models_select(Bmodel* model, Bstring& selstr);
long 		models_select(Bmodel* model, string& selstr);
long		models_select(Bmodel* model, long number);
long		models_select_all(Bmodel* model);
long		models_select_unknowns(Bmodel* model);
long		models_unset_selection(Bmodel* model);
long		models_invert_selection(Bmodel* model, int type);
long		models_select_sets(Bmodel* model, int size, int flag);
long		models_select_within_bounds(Bmodel* model, Vector3<double>& start, Vector3<double>& end);
long		models_select_number_of_components(Bmodel* model, long ncomp_min, long ncomp_max);
int			model_select_random(Bmodel* model, long number);
long		models_select_closed(Bmodel* model, int closure_rule, int val_order);
long		models_select_fullerene(Bmodel* model);
long		models_select_non_fullerene(Bmodel* model);
long		models_select_valence(Bmodel* model, int valence);
long		models_select_polygons(Bmodel* model, int order);
long 		models_select_first(Bmodel* model, int first);
long 		models_select_within_shell(Bmodel* model, Vector3<double> center, double minrad, double maxrad);
long 		models_select_in_mask(Bmodel* model, Bimage* pmask);
long		models_select_slices(Bmodel* model, double bottom, double top, double thickness);
long 		models_delete(Bmodel** model);
long 		models_delete_comp_type(Bmodel* model, string comptype);
long 		models_delete_non_selected(Bmodel** model);
long		model_type_from_selection(Bmodel* model, Bstring* comp_type, string filename);
long		model_fom_deselect(Bmodel* model, double fom_cutoff);
long		model_fom_max_fraction_deselect(Bmodel* model, double fom_fraction);
long		model_fom_histogram(Bmodel* model, double fom_step);
long		model_fom_ranking(Bmodel* model, int nrank);
long		models_delete_overlapped_components(Bmodel** model, double distance);
long		model_average_overlapped_components(Bmodel* model, double distance);
long		model_find_overlap(Bmodel* model, string reffile, double distance);
long		models_prune(Bmodel* model, int prune_type, double distance);


