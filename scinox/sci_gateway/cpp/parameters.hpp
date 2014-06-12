#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

enum  type_check {CHECK_NONE, CHECK_SIZE, CHECK_MIN, CHECK_MAX, CHECK_BOTH, CHECK_VALUES};

int find_label_partial(char ** LabelList, int nblabels, char const * const LabelToFind);
int find_label(char ** LabelList, int nblabels, char const * const LabelToFind);

int init_parameters(int Pos, int ** param_addr);
int check_parameters(int * param_addr);
int has_label_partial(int * param_addr, const char * Label);
int get_int_parameter(int * param_addr, const char * Label, int * value, int * found, int default_value, type_check check, ...);
int get_double_parameter(int * param_addr, const char * Label, double * value, int * found, double default_value, type_check check, ...);
int get_vec_int_parameter(int * param_addr, const char * Label, int * value, int * found, int default_value, int * size, type_check check, ...);
int get_vec_double_parameter(int * param_addr, const char * Label, double * value, int * found, double default_value, int * size, type_check check, ...);
int get_string_parameter(int * param_addr, const char * Label, char ** value, int * found, const char * default_value, type_check check, ...);
#endif
