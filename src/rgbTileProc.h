/* rgbTileProc.h
*  implementation of TILE COMPRESSION
*/

//typedef.h
#ifndef _TYPE_DEFS
#define _TYPE_DEFS

typedef unsigned char   	uint8;
typedef unsigned short		uint16;
typedef unsigned long		uint32;
typedef signed char     	int8;
typedef signed short		int16;
typedef signed long		int32;
typedef float			float32;
typedef double			float64;
// typedef enum{0, 1}	bool;

#endif

//parameter.h
#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LEN	100

struct parameters{
	//char input_file[MAX_FILENAME_LEN];	// input filename
	//char output_file[MAX_FILENAME_LEN];	// output filename
	bool decoding_flag;			// encoding/decoding flag
	bool verbose;				// verbose flag
	uint8 NEAR;				// difference bound for near lossless coding
	uint16 MAXVAL;				// max image sample value
	uint16 T1, T2, T3;			// thresholds for local gradients
	bool specified_T;
	uint16 RESET;				// threshold value at which A,B, and N are halved
	uint8 ILV;				// interleave     
} typedef parameters;

parameters coding_parameters(bool decoding_flag);

#endif

//codingvars.h
#ifndef __CODINGVARS_H__
#define __CODINGVARS_H__

#include <math.h>

#define		CONTEXTS	365

struct codingvars{

	bool RunModeProcessing;		// Regular/Run mode processing flag
	uint16 RANGE;			// range of prediction error representation
	uint8 bpp;			// number of bits needed to represent MAXVAL
	uint8 qbpp;			// number of bits needed to represent a mapped error value
	uint8 LIMIT;			// max length in bits of Golomb codewords in regular mode

	uint16 comp, row, col;

	int32 N[CONTEXTS + 2];		// occurrences counters for each context
	uint32 A[CONTEXTS + 2];		// prediction error magnitudes accumulator
	int32 B[CONTEXTS];		// bias values
	int32 C[CONTEXTS];		// prediction correction values
	uint8 RUNindex	;		// index for run mode order
	uint8 RUNindex_val;	
	uint16 RUNval;			// repetitive reconstructed sample value in a run
	uint16 RUNcnt;			// repetetive sample count for run mode
	uint32 TEMP;			// auxiliary variable for Golomb variable calculation in run interruption coding
	uint8 map;			// auxiliary variable for error mapping at run interruption
	uint8 RItype;			// index for run interruption coding
	
	uint8* J;			// order of run-length codes
	int32 Nn[2];			// counters for negative prediction error for run interruption

	uint16 Ra, Rb, Rc, Rd, prevRa, Ix;	// pixels used in the causal template
	int32 Rx;				// reconstructed value of the current sample
	int32 Px;				// predicted value for the current sample
	int32 Errval;				// prediction error
	uint32 MErrval;				// Errval mapped to a non-negative integer
	uint32 EMErrval;			// Errval mapped to non-negative integers in run interruption mode 

	uint16 Q;				// context
	int8 SIGN;				// sign of the current context

	uint8 k;				// Golomb coding variable

	int32 MAX_C;				// maximum allowed value for C[0..364]
	int32 MIN_C;				// minimum allowed value for C[0..364]

	uint16 BASIC_T1;			// basic default values for gradient quantization thresholds
	uint16 BASIC_T2;
	uint16 BASIC_T3;

	uint16 FACTOR;

} typedef codingvars;

void init_codingvars(codingvars* vars, parameters params);
void update_codingvars(codingvars* vars, parameters params);

#endif

//create_image.h
#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct image_data {
    uint8*** image;		// pointer to image data
    uint16 		width;		// max supported image size is 65535x65535 pixels (~4.3*10^9 pixels)
    uint16 		height;
    uint16 		maxval;		// max precision is 16bpp for each component
    uint8 		n_comp;		// number of components (max 255)
} typedef image_data;

image_data* create_image();
image_data* allocate_image(image_data* im_data);
void print_image(image_data* im_data);
bool compare_image(image_data* im_data, image_data* im_data_out);

#endif

//bitstream.h
#ifndef __BITSTREAM_H__
#define __BITSTREAM_H__

#include <stdio.h>
#include <stdlib.h>

struct bitstream{
	uint8* pTile;
	uint8 byte_bits;
	uint16 tot_bytes;
	uint16 byte_count;
	uint8 buffer;
	uint32 tot_bits;
} typedef bitstream;

// void init_bitstream(uint8* outTile, char mode, int pTile_size);
void init_bitstream_w(uint8* outTile);
void init_bitstream_r(const uint8* outTile,uint8* tempTile, int pTile_size);
void print_bpp(image_data* im_data);

void append_bit(uint8 bit);
void append_bits(uint32 value, uint8 n_bits);

void print_bitstream();

uint8 read_bit();
uint16 read_bits(int num);

void deal_EOF();


#endif

//golomb.h
#ifndef __GOLOMB_H__
#define __GOLOMB_H__

#include <stdio.h>
#include <stdlib.h>

void limited_length_Golomb_encode(uint32 MErrval, uint8 k, uint8 LIMIT, uint8 qbpp);
uint32 limited_length_Golomb_decode(uint8 k, uint8 LIMIT, uint8 qbpp);

#endif

//predictivecoding.h
#ifndef __PREDICTIVECODING_H__
#define __PREDICTIVECODING_H__

#include <math.h>

void context_determination(codingvars* vars, parameters params, image_data* im_data);
void predict_sample_value(codingvars* vars, parameters params);
void encode_prediction_error(codingvars* vars, parameters params, image_data* im_data);
void decode_prediction_error(codingvars* vars, parameters params, image_data* im_data);
void encode_run(codingvars* vars, parameters params, image_data* im_data);
void decode_run(codingvars* vars, parameters params, image_data* im_data);

#endif

#ifndef _RGBTILEPROC_H_
#define _RGBTILEPROC_H_

void tileSetSize(int nTileWidth, int nTileHeight);

bool judge_lossless(const unsigned char* pClrBlk,const unsigned char* pTile, int nTileSize);

/* compress ARGB data to tile
*  param:
*    pClrBlk      -- IN, pixel's ARGB data
*    pTile        -- OUT, tile data
*    pTileSize    -- OUT, tile's bytes
*  return:
*    0  -- succeed
*   -1  -- failed
*/
int argb2tile(const unsigned char *pClrBlk, unsigned char *pTile, int *pTileSize);

/* decompress tile data to ARGB
*  param:
*    pTile        -- IN, tile data
*    pTileSize    -- IN, tile's bytes
*    pClrBlk      -- OUT, pixel's ARGB data
*  return:
*    0  -- succeed
*   -1  -- failed
*/
int tile2argb(const unsigned char* pTile, int nTileSize, unsigned char* pClrBlk);

#endif
