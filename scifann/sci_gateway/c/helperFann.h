/*
 * Scilab interface for the FANN library
 * Author: Dirk Gorissen <dirk.gorissen@ua.ac.be>
 * Licence: GPL version 2 or later
 */

#ifndef HELPERFANN_INCLUDE
#define HELPERFANN_INCLUDE


#include <fann.h>
#include <fann_train.h>
#include <fann_cascade.h>
#include <fann_data.h>

struct fann_train_data *read_from_array(const double *din,
					const double *dout,
					const unsigned int num_data,
					const unsigned int num_input,
					const unsigned int num_output);

struct fann* createNetwork(const unsigned int numLayers,
			   const unsigned int* layers,
			   const float connectionRate);

struct fann* trainNetwork(struct fann *ann,
			  struct fann_train_data *data,
			  const float desiredError,
			  const unsigned int maxEpochs);

void evaluateNetwork(struct fann *ann, 
		     const double *input, 
		     double* output,
		     const unsigned int numData);

int detect_fannlist(int StackPos);
int detect_fanntraindatalist(int StackPos);
int detect_fannerrorlist(int StackPos);

int createScilabFannStructFromCFannStruct(struct fann* ann, unsigned int StackPos);
int createScilabFannTrainDataStructFromCFannTrainDataStruct(struct fann_train_data* ann_data, unsigned int StackPos);
int createScilabFannErrorStructFromCFannErrorStruct(struct fann_error* ann, FILE * log_file, unsigned int StackPos);

struct fann * createCFannStructFromScilabFannStruct(unsigned int StackPos, int * res);
struct fann_train_data * createCFannTrainDataStructFromScilabFannTrainDataStruct(unsigned int StackPos, int * res);
struct fann_error * createCFannErrorStructFromScilabFannErrorStruct(unsigned int StackPos, int * res);
FILE * createFILEFromScilabFannErrorStruct(unsigned int StackPos, int * res);
#endif
