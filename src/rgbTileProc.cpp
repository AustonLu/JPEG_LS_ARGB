/* rgbTileProc.cpp
*  implementation of TILE COMPRESSION
*/
#include <stdio.h>
#include <math.h>
#include<stdio.h>
#include <memory.h>
#include <assert.h>
#include "rgbTileProc.h"

static int g_nTileWidth = 0;
static int g_nTileHeight = 0;
static int lossless_count=0;
static int nonlossless_count=0;

void tileSetSize(int nTileWidth, int nTileHeight)
{
	g_nTileWidth = nTileWidth;
	g_nTileHeight = nTileHeight;
}

/*Main Public Calculation*/

static inline int32 CLAMP(int32 a, int32 b, int32 c) { return (a > c || a < b) ? b : a; }
static inline uint32 min(uint32 a,uint32 b){ return (a > b) ? b : a; }
static inline uint32 max(uint32 a,uint32 b){ return (a > b) ? a : b; }

/*parameter.c*/
parameters coding_parameters(bool decoding_flag)
{
	parameters params;
	params.verbose = false;
	params.NEAR = 0;
	params.RESET = 64;
	params.specified_T = false;
	params.decoding_flag = decoding_flag;
	params.T1 = 3;
	params.T2 = 7;
	params.T3 = 21;


	return params;
}

/*codingvars.c*/

static uint8 J_values[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,5,5,6,6,7,7,8,9,10,11,12,13,14,15};
void init_codingvars(codingvars* vars, parameters params)
{
	/* A.2 Initializations and conventions */

	/* A.2.1 Initializations */

	uint16 i;
	uint32 A_init_value;

	(*vars).RUNindex = 0;		
	(*vars).J = J_values;
	(*vars).prevRa = 0;
	(*vars).MAX_C=127;			
	(*vars).MIN_C=-128;			
	(*vars).BASIC_T1 = 3;			
	(*vars).BASIC_T2 = 7;
	(*vars).BASIC_T3 = 21;

	(*vars).RANGE = floor( (float64)(params.MAXVAL + 2*params.NEAR) / (2*params.NEAR + 1) ) + 1;
	(*vars).qbpp = ceil( log2((*vars).RANGE) );
	(*vars).bpp = max( 2, log2(params.MAXVAL + 1) );
	(*vars).LIMIT = 2*( (*vars).bpp + max(8,(*vars).bpp) );

	A_init_value = max( 2, floor( (float64)((*vars).RANGE + 32)/64 ) );

	for(i=0; i<CONTEXTS; i++)
	{
		(*vars).A[i] = A_init_value;
		(*vars).N[i] = 1;
		(*vars).B[i] = 0;
		(*vars).C[i] = 0;
	}
	(*vars).A[CONTEXTS] = A_init_value;
	(*vars).A[CONTEXTS+1] = A_init_value;
	(*vars).N[CONTEXTS] = 1;
	(*vars).N[CONTEXTS+1] = 1;
	(*vars).Nn[0] = 0;
	(*vars).Nn[1] = 0;
}
void update_codingvars(codingvars* vars, parameters params)
{
	/* A.6.1 Update */

	(*vars).B[(*vars).Q] += (*vars).Errval*(2*params.NEAR + 1);
	(*vars).A[(*vars).Q] += abs((*vars).Errval);
	if((*vars).N[(*vars).Q] == params.RESET)
	{
		(*vars).A[(*vars).Q] = (*vars).A[(*vars).Q]>>1;
		if((*vars).B[(*vars).Q]>=0)
			(*vars).B[(*vars).Q] = (*vars).B[(*vars).Q]>>1;
		else
			(*vars).B[(*vars).Q] = -((1-(*vars).B[(*vars).Q])>>1);
		(*vars).N[(*vars).Q] = (*vars).N[(*vars).Q]>>1;
	}
	(*vars).N[(*vars).Q] += 1;

	/* A.6.2 Bias computation */

	if((*vars).B[(*vars).Q] <= -(*vars).N[(*vars).Q])
	{
		(*vars).B[(*vars).Q] += (*vars).N[(*vars).Q];
		if((*vars).C[(*vars).Q] > (*vars).MIN_C)
			(*vars).C[(*vars).Q]--;
		if((*vars).B[(*vars).Q] <= -(*vars).N[(*vars).Q])
			(*vars).B[(*vars).Q] = -(*vars).N[(*vars).Q] + 1;
	}
	else if((*vars).B[(*vars).Q] > 0)
	{
		(*vars).B[(*vars).Q] -= (*vars).N[(*vars).Q];
		if((*vars).C[(*vars).Q] < (*vars).MAX_C)
			(*vars).C[(*vars).Q]++;
		if((*vars).B[(*vars).Q] > 0)
			(*vars).B[(*vars).Q] = 0;
	}
}

/*create_image.c*/
image_data* create_image(const unsigned char* pClrBlk_in) {

	image_data* im_data;
	im_data = (image_data*)malloc(sizeof(image_data));
	im_data->n_comp = 4;
	int input_size = sizeof(pClrBlk_in)/ im_data->n_comp;
	// im_data->width = (int)sqrt(input_size);
	// im_data->height = (int)sqrt(input_size);
	im_data->width = 8;
	im_data->height = 8;
	im_data->maxval = 255;

	//convert into 3dimagge
	// allocate image space
	//uint8 pClrBlk_in_3D[im_data->n_comp][im_data->width][im_data->height];
	uint8*** pClrBlk_in_3D;
	pClrBlk_in_3D = (uint8***)malloc(im_data->n_comp * sizeof(uint8**));
	if (pClrBlk_in_3D == NULL)
	{
		fprintf(stderr, "Error in memory allocation\n");
		exit(EXIT_FAILURE);
	}

	for (int c = 0; c < im_data->n_comp; c++)
	{
		pClrBlk_in_3D[c] = (uint8**)malloc(im_data->height * sizeof(uint8*));
		if (pClrBlk_in_3D[c] == NULL)
		{
			fprintf(stderr, "Error in memory allocation\n");
			exit(EXIT_FAILURE);
		}
		for (int row = 0; row < im_data->height; row++)
		{
			pClrBlk_in_3D[c][row] = (uint8*)malloc(im_data->width * sizeof(uint8));
			if (pClrBlk_in_3D[c][row] == NULL)
			{
				fprintf(stderr, "Error in memory allocation\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	for (int i = 0; i < im_data->height; i++)
	{
		for (int j = 0; j < im_data->width; j++)
		{
			pClrBlk_in_3D[0][i][j] = pClrBlk_in[(i * im_data->width + j)*4];
			pClrBlk_in_3D[1][i][j] = pClrBlk_in[(i * im_data->width + j)*4 + 1];
			pClrBlk_in_3D[2][i][j] = pClrBlk_in[(i * im_data->width + j)*4 + 2];
			pClrBlk_in_3D[3][i][j] = pClrBlk_in[(i * im_data->width + j)*4 + 3];
		}
	}
	im_data->image = pClrBlk_in_3D;
	return im_data;
}
image_data* allocate_image()
{
	//create a image_data with the same parameter with the input image(use for decode)
	uint8 c;
	uint16 row;
	image_data* im_data_out;
	im_data_out = (image_data*)malloc(sizeof(image_data));

	im_data_out->n_comp = 4;
	im_data_out->height = 8;
	im_data_out->width =8;
	im_data_out->maxval = 255;

	uint8*** pClrBlk_3D;
	pClrBlk_3D = (uint8***)malloc(im_data_out->n_comp * sizeof(uint8**));
	if (pClrBlk_3D == NULL)
	{
		fprintf(stderr, "Error in memory allocation\n");
		exit(EXIT_FAILURE);
	}

	for (int c = 0; c < im_data_out->n_comp; c++)
	{
		pClrBlk_3D[c] = (uint8**)malloc(im_data_out->height * sizeof(uint8*));
		if (pClrBlk_3D[c] == NULL)
		{
			fprintf(stderr, "Error in memory allocation\n");
			exit(EXIT_FAILURE);
		}
		for (int row = 0; row < im_data_out->height; row++)
		{
			pClrBlk_3D[c][row] = (uint8*)malloc(im_data_out->width * sizeof(uint8));
			if (pClrBlk_3D[c][row] == NULL)
			{
				fprintf(stderr, "Error in memory allocation\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	im_data_out->image = pClrBlk_3D;
	return im_data_out;

}
void print_image(image_data* im_data) {
	// im_data->n_comp
	printf("{");
	{
		// printf("[");
		for (int row = 0; row < im_data->height; row++) {
			// printf(" ");
			for (int col = 0; col < im_data->width; col++)
			{
				for (int comp = 0; comp < im_data->n_comp; comp++) 
					printf("%d,", im_data->image[comp][row][col]);
			}
			// printf(" \n");
		}
		// printf("]\n");
	}
	printf("}\n");
	// printf("{");
	// for (int comp = 0; comp < im_data->n_comp; comp++) {
	// 	// printf("[");
	// 	for (int row = 0; row < im_data->height; row++) {
	// 		// printf(" ");
	// 		for (int col = 0; col < im_data->width; col++)
	// 			printf("%d,", im_data->image[comp][row][col]);
	// 		// printf(" \n");
	// 	}
	// 	// printf("]\n");
	// }
	// printf("}\n");
		
			
}
bool compare_image(image_data* im_data, image_data* im_data_out) {
	for (int comp = 0; comp < im_data->n_comp; comp++) {
		for (int row = 0; row < im_data->height; row++) {
			for (int col = 0; col < im_data->width; col++) {
				if (im_data->image[comp][row][col] != im_data_out->image[comp][row][col])
					return false;
			}
		}
	}
	return true;
}


/*bitstream.c*/
bitstream bs;
void init_bitstream_w(uint8* outTile){
	bs.pTile = outTile;
	bs.tot_bits = 0;
	bs.byte_bits = 0;
	bs.tot_bytes = 0;
	bs.byte_count=0;
	bs.buffer=0;
}
void init_bitstream_r(const uint8* outTile,uint8* tempTile, int pTile_size){
	bs.pTile = tempTile;
	bs.tot_bits = 0;
	bs.byte_bits = 0;
	bs.byte_count = 0;
	bs.buffer=0;
	bs.tot_bytes = pTile_size;
	for(int i=0;i<pTile_size;i++){
		bs.pTile[i]=outTile[i];
		// printf("%d,",outTile[i]);
	}
	// printf("\n");
}
void print_bpp(image_data* im_data)
{
	float32 bpp = 	(float32)bs.tot_bits / (im_data->height * im_data->width * im_data->n_comp);
	fprintf(stdout, "\n%ld bits,%.2f bpp\n", bs.tot_bits,bpp);
}
void append_bit(uint8 bit)
{
	bs.tot_bits++;
	if(bs.byte_bits==8)
	{
		bs.pTile[bs.tot_bytes] = bs.buffer;
		bs.buffer = 0;
		bs.byte_bits = 0;
		//printf("!");
		bs.tot_bytes++;
		//printf("\n");
	}

	bs.byte_bits++;
	bs.buffer = bs.buffer * 2 + bit;
	// printf("%d",bit);

}
void append_bits(uint32 value, uint8 bits)
{
	while(bits > 0)
	{
		// msb first
		append_bit((value>>(bits-1))&0x1);
		bits--;
	}
}
uint8 read_bit()   
{
	uint8 bit;

	if (bs.byte_bits == 0)
	{
		bs.buffer = bs.pTile[bs.byte_count];
		bs.byte_bits = 8;
		//printf("!");
		bs.byte_count++;
		//printf("\n");
	}

	// msb first
	bit = bs.buffer / (int)pow(2,bs.byte_bits-1);
	bs.buffer= bs.buffer % (int)pow(2, bs.byte_bits - 1);
	bs.byte_bits--;
	// printf("%d", bit);

	return bit;
}
uint16 read_bits(int num)
{
	uint16 bitsresult=0;
	for (int i = 0; i < num; i++) {
		bitsresult *= 2;
		bitsresult += read_bit();
	}
	
	return bitsresult;

}
void print_bitstream() {
	printf("-----------------------\n");
	for (int byte = 0; byte < bs.tot_bytes+1; byte++) {
		printf("%d,", bs.pTile[byte]);
	}
	printf("\nTotle Bytes:%d\n",bs.tot_bytes);
}
void deal_EOF() {
	if (bs.byte_bits != 0) {
		bs.buffer = bs.buffer * pow(2, 8 - bs.byte_bits);
		bs.pTile[bs.tot_bytes] = bs.buffer;
	}
}


/*golomb.c*/
void limited_length_Golomb_encode(uint32 MErrval, uint8 k, uint8 LIMIT, uint8 qbpp)
{
	// A.5.3 Mapped-error encoding

	// bits generated are packed into 8-bit bytes
	// first output bit is on the msb

	uint32 hMErrval;
	uint8 lim = LIMIT-qbpp-1;

	hMErrval = MErrval>>k;

	//printf("%d,",k);

	if(hMErrval < lim)
	{
		// unary part		
		while(hMErrval>0)
		{
			append_bit(0);
			hMErrval--;
		}
		append_bit(1);
		//binary part
		while(k>0)
		{
			append_bit((MErrval>>(k-1))&0x1);
			k--;		
		}
	}
	else
	{
		// unary part		
		while(lim>0)
		{
			append_bit(0);
			lim--;
		}
		append_bit(1);
		//binary part
		while(qbpp>0)
		{
			append_bit(((MErrval-1)>>(qbpp-1))&0x1);
			qbpp--;		
		}
	}
	//printf("\n");

	return;
}
uint32 limited_length_Golomb_decode(uint8 k, uint8 LIMIT, uint8 qbpp)
{
	
	uint32 MErrval = 0;
	uint8 lim = LIMIT-qbpp-1;
	//printf("%d,", k);

	while(read_bit()==0)
		MErrval++;


	if(MErrval<lim)
	{
		while(k>0)
		{
			MErrval = MErrval * 2 + read_bit();
			/*MErrval = (MErrval<<1)|read_bit();*/
			k--;
		}
	}
	else
	{
		MErrval = 0;
		while(qbpp>0)
		{
			MErrval = MErrval * 2 + read_bit();
			/*MErrval = (MErrval<<1)|read_bit();*/
			qbpp--;
		}
		MErrval += 1;
	}
	//printf("\n");

	return MErrval;

}


/*predictivecoding.c*/
void context_determination(codingvars* vars, parameters params, image_data* im_data)
{

	int32 D1, D2, D3;			// local gradients
	int8 Q1, Q2, Q3;			// region numbers to quantize local gradients

	// causal template construction
	// c b d
	// a x

	(*vars).Rb = ((*vars).row==0) ? 0 : im_data->image[(*vars).comp][(*vars).row-1][(*vars).col];
	(*vars).Rd = ((*vars).col==im_data->width-1 || (*vars).row==0) ? (*vars).Rb : im_data->image[(*vars).comp][(*vars).row-1][(*vars).col+1];
	if((*vars).col>0)
		(*vars).Rc = ((*vars).row==0) ? 0 : im_data->image[(*vars).comp][(*vars).row-1][(*vars).col-1];
	else
		(*vars).Rc = (*vars).prevRa;
	if ((*vars).col == 0) {
		(*vars).Ra = (*vars).Rb;
		(*vars).prevRa = (*vars).Ra;
	}
	else
		(*vars).Ra = im_data->image[(*vars).comp][(*vars).row][(*vars).col - 1];
	


	(*vars).Ix = im_data->image[(*vars).comp][(*vars).row][(*vars).col];

	/* A.3.1 Local gradient computation */

	D1 = (*vars).Rd - (*vars).Rb;
	D2 = (*vars).Rb - (*vars).Rc;
	D3 = (*vars).Rc - (*vars).Ra;

	/* A.3.2 Mode selection */

	bool flag;
	

	if((abs(D1) <= params.NEAR)&&(abs(D2) <= params.NEAR)&&(abs(D3) <= params.NEAR))
		flag = true;
	else
		flag = false;

	if ((*vars).comp == 0 && (*vars).row == 0 && (*vars).col == 0) {
		flag = false;
	}

	(*vars).RunModeProcessing = flag;

	if(!(*vars).RunModeProcessing)	// regular mode
	{
		/* A.3.3 Local gradients quantization */

		if      (D1 <= -params.T3)   	Q1 =-4;
		else if (D1 <= -params.T2)   	Q1 =-3;
		else if (D1 <= -params.T1)   	Q1 =-2;
		else if (D1 <  -params.NEAR) 	Q1 =-1;
		else if (D1 <=  params.NEAR) 	Q1 = 0;
		else if (D1 <   params.T1)   	Q1 = 1;
		else if (D1 <   params.T2)   	Q1 = 2;
		else if (D1 <   params.T3)   	Q1 = 3;
		else                  		Q1 = 4;

		if      (D2 <= -params.T3)   	Q2 =-4;
		else if (D2 <= -params.T2)   	Q2 =-3;
		else if (D2 <= -params.T1)   	Q2 =-2;
		else if (D2 <  -params.NEAR) 	Q2 =-1;
		else if (D2 <=  params.NEAR) 	Q2 = 0;
		else if (D2 <   params.T1)   	Q2 = 1;
		else if (D2 <   params.T2)   	Q2 = 2;
		else if (D2 <   params.T3)   	Q2 = 3;
		else                  		Q2 = 4;

		if      (D3 <= -params.T3)   	Q3=-4;
		else if (D3 <= -params.T2)   	Q3=-3;
		else if (D3 <= -params.T1)   	Q3=-2;
		else if (D3 <  -params.NEAR) 	Q3=-1;
		else if (D3 <=  params.NEAR) 	Q3= 0;
		else if (D3 <   params.T1)   	Q3= 1;
		else if (D3 <   params.T2)   	Q3= 2;
		else if (D3 <   params.T3)   	Q3= 3;
		else 				Q3= 4;

		/* A.3.4 Quantized gradient merging */

		if( (Q1<0) || (Q1==0 && Q2<0) || (Q1==0 && Q2==0 && Q3<0) )
		{
			(*vars).SIGN = -1;
			Q1 = -Q1;
			Q2 = -Q2;
			Q3 = -Q3;
		}
		else
			(*vars).SIGN = 1;

		// one-to-one mapping of the vector (Q1,Q2,Q3) to the integer Q
		if (Q1 == 0)
		{
			if (Q2 == 0)
				(*vars).Q=360+Q3;
			else
				(*vars).Q=324+(Q2-1)*9+(Q3+4);
		}
		else
			(*vars).Q=(Q1-1)*81+(Q2+4)*9+(Q3+4);
	}
}
void predict_sample_value(codingvars* vars, parameters params)
{
	/* A.4.1 Edge-detecting predictor */

	if((*vars).Rc>=max((*vars).Ra,(*vars).Rb))
		(*vars).Px = min((*vars).Ra,(*vars).Rb);
	else
	{
		if((*vars).Rc<=min((*vars).Ra,(*vars).Rb))
			(*vars).Px = max((*vars).Ra,(*vars).Rb);
		else
			(*vars).Px = (*vars).Ra + (*vars).Rb - (*vars).Rc;
	}

	/* A.4.2 Prediction correction */

	(*vars).Px += ((*vars).SIGN == 1)? (*vars).C[(*vars).Q] : -(*vars).C[(*vars).Q];

	if((*vars).Px>params.MAXVAL)
		(*vars).Px = params.MAXVAL;
	else if((*vars).Px<0)
		(*vars).Px = 0;
}
void encode_prediction_error(codingvars* vars, parameters params, image_data* im_data)
{
	/* A.4.2 Computation of prediction error */

	(*vars).Errval = (*vars).Ix - (*vars).Px;
	if((*vars).SIGN==-1)
		(*vars).Errval = -(*vars).Errval;

	/* A.4.4 Error quantization for near-lossless coding, and reconstructed value computation */
	(*vars).Rx = (*vars).Ix;
	im_data->image[(*vars).comp][(*vars).row][(*vars).col] = (*vars).Rx;

	// modulo reduction of the error
	if((*vars).Errval<0)
		(*vars).Errval = (*vars).Errval + (*vars).RANGE;
	if((*vars).Errval>=(((*vars).RANGE + 1)/2))
		(*vars).Errval = (*vars).Errval - (*vars).RANGE;

	/* A.5 Prediction error encoding */

	/* A.5.1 Golomb coding variable computation */

	for((*vars).k=0;((*vars).N[(*vars).Q]<<(*vars).k)<(*vars).A[(*vars).Q];(*vars).k++);

	/* A.5.2 Error mapping */

	if((params.NEAR==0)&&((*vars).k==0)&&(2*(*vars).B[(*vars).Q]<=-(*vars).N[(*vars).Q]))
		if((*vars).Errval>=0)
			(*vars).MErrval = 2*(*vars).Errval + 1;
		else
			(*vars).MErrval = -2*((*vars).Errval + 1);
	else
		if((*vars).Errval>=0)
			(*vars).MErrval = 2*(*vars).Errval;
		else
			(*vars).MErrval = -2*(*vars).Errval -1;

	/* A.5.3 Mapped-error encoding */
	limited_length_Golomb_encode((*vars).MErrval, (*vars).k, (*vars).LIMIT, (*vars).qbpp);
	// printf("\tREGU\n");
}
void decode_prediction_error(codingvars* vars, parameters params, image_data* im_data)
{
	/* A.5.1 Golomb coding variable computation */

	int32 ErrvalAfterMR;

	for((*vars).k=0;((*vars).N[(*vars).Q]<<(*vars).k)<(*vars).A[(*vars).Q];(*vars).k++);

	/* Mapped-error decoding */

	(*vars).MErrval = limited_length_Golomb_decode((*vars).k, (*vars).LIMIT, (*vars).qbpp);

	/* Inverse Error mapping */

	if((params.NEAR==0)&&((*vars).k==0)&&(2*(*vars).B[(*vars).Q]<=-(*vars).N[(*vars).Q]))
		if((*vars).MErrval%2==0)
			(*vars).Errval = -((int32)(*vars).MErrval / 2) - 1;
		else
			(*vars).Errval = ((int32)(*vars).MErrval - 1) / 2;
	else
		if((*vars).MErrval%2==0)
			(*vars).Errval = (int32)(*vars).MErrval / 2;
		else
			(*vars).Errval = -((int32)(*vars).MErrval + 1) / 2;

	ErrvalAfterMR = (*vars).Errval;	

	(*vars).Errval = (*vars).Errval * (int32)(2*params.NEAR + 1);

	if((*vars).SIGN==-1)
		(*vars).Errval = -(*vars).Errval;

	
	(*vars).Rx = ((*vars).Errval + (*vars).Px) ;
	
	im_data->image[(*vars).comp][(*vars).row][(*vars).col] = (*vars).Rx;
	(*vars).Errval = ErrvalAfterMR;

}
void encode_run(codingvars* vars, parameters params, image_data* im_data)
{
	/* A.7.1 Run scanning and run-length coding */
			
	(*vars).RUNval = (*vars).Ra;
	(*vars).RUNcnt = 0;
	while(abs((*vars).Ix - (*vars).RUNval) <= params.NEAR)
	{
		(*vars).RUNcnt += 1;
		(*vars).Rx = (*vars).RUNval;
		if((*vars).col == (im_data->width-1))
			break;
		else
		{
			(*vars).col++;
			(*vars).Ix = im_data->image[(*vars).comp][(*vars).row][(*vars).col];
		}
	}

	/* A.7.1.2 Run-length coding */

	while((*vars).RUNcnt >= (1<<(*vars).J[(*vars).RUNindex]))
	{
		append_bit(1);
		(*vars).RUNcnt -= (1<<(*vars).J[(*vars).RUNindex]);
		if((*vars).RUNindex<31)
			(*vars).RUNindex += 1;
	}

	if(abs((*vars).Ix - (*vars).RUNval) > params.NEAR)
	{
		append_bit(0);
		append_bits((*vars).RUNcnt,(*vars).J[(*vars).RUNindex]);
		if((*vars).RUNindex > 0)
			(*vars).RUNindex -= 1;

		/* A.7.2 Run interruption sample encoding */

	// index computation
		(*vars).Rb = ((*vars).row == 0) ? 0 : im_data->image[(*vars).comp][(*vars).row - 1][(*vars).col];
		if ((*vars).col == 0) {
			(*vars).Ra = (*vars).Rb;
			(*vars).prevRa = (*vars).Ra;
		}
		else
			(*vars).Ra = im_data->image[(*vars).comp][(*vars).row][(*vars).col - 1];

		if (abs((*vars).Ra - (*vars).Rb) <= params.NEAR)
			(*vars).RItype = 1;
		else
			(*vars).RItype = 0;

		// prediction error for a run interruption sample
		if ((*vars).RItype == 1)
			(*vars).Px = (*vars).Ra;
		else
			(*vars).Px = (*vars).Rb;
		(*vars).Errval = (*vars).Ix - (*vars).Px;

		// error computation for a run interruption sample
		if (((*vars).RItype == 0) && ((*vars).Ra > (*vars).Rb))
		{
			(*vars).Errval = -(*vars).Errval;
			(*vars).SIGN = -1;
		}
		else
			(*vars).SIGN = 1;

		if (params.NEAR > 0)
		{
			// error quantization
			if ((*vars).Errval > 0)
				(*vars).Errval = ((*vars).Errval + params.NEAR) / (2 * params.NEAR + 1);
			else
				(*vars).Errval = -(params.NEAR - (*vars).Errval) / (2 * params.NEAR + 1);

			// reconstructed value computation
			(*vars).Rx = (*vars).Px + (*vars).SIGN * (*vars).Errval * (2 * params.NEAR + 1);
			if ((*vars).Rx < 0)
				(*vars).Rx = 0;
			else if ((*vars).Rx > params.MAXVAL)
				(*vars).Rx = params.MAXVAL;
		}
		else
			(*vars).Rx = (*vars).Ix;

		// modulo reduction of the error
		if ((*vars).Errval < 0)
			(*vars).Errval = (*vars).Errval + (*vars).RANGE;
		if ((*vars).Errval >= (((*vars).RANGE + 1) / 2))
			(*vars).Errval = (*vars).Errval - (*vars).RANGE;

		// computation of the auxiliary variable TEMP
		if ((*vars).RItype == 0)
			(*vars).TEMP = (*vars).A[365];
		else
			(*vars).TEMP = (*vars).A[366] + ((*vars).N[366] >> 1);

		// Golomb coding variable computation
		(*vars).Q = (*vars).RItype + 365;
		for ((*vars).k = 0; ((*vars).N[(*vars).Q] << (*vars).k) < (*vars).TEMP; (*vars).k++);

		// computation of map for Errval mapping
		if (((*vars).k == 0) && ((*vars).Errval > 0) && (2 * (*vars).Nn[(*vars).Q - 365] < (*vars).N[(*vars).Q]))
			(*vars).map = 1;
		else if (((*vars).Errval < 0) && (2 * (*vars).Nn[(*vars).Q - 365] >= (*vars).N[(*vars).Q]))
			(*vars).map = 1;
		else if (((*vars).Errval < 0) && ((*vars).k != 0))
			(*vars).map = 1;
		else
			(*vars).map = 0;

		// Errval mapping for run interruption sample
		(*vars).EMErrval = 2 * abs((*vars).Errval) - (*vars).RItype - (*vars).map;

		// limited length Golomb encoding
		limited_length_Golomb_encode((*vars).EMErrval, (*vars).k, (*vars).LIMIT - (*vars).J[(*vars).RUNindex] - 1, (*vars).qbpp);

		// update of variables for run interruption sample
		if ((*vars).Errval < 0)
			(*vars).Nn[(*vars).Q - 365] = (*vars).Nn[(*vars).Q - 365] + 1;
		(*vars).A[(*vars).Q] = (*vars).A[(*vars).Q] + (((*vars).EMErrval + 1 + (*vars).RItype) >> 1);
		if ((*vars).N[(*vars).Q] == params.RESET)
		{
			(*vars).A[(*vars).Q] = (*vars).A[(*vars).Q] >> 1;
			(*vars).N[(*vars).Q] = (*vars).N[(*vars).Q] >> 1;
			(*vars).Nn[(*vars).Q - 365] = (*vars).Nn[(*vars).Q - 365] >> 1;
		}
		// printf(" R0 ,%d,%d,%d,%d,", (*vars).comp, (*vars).row, (*vars).col, im_data->image[(*vars).comp][(*vars).row][(*vars).col]);
		// printf("\tRUN\n");

	}
	else if ((*vars).RUNcnt > 0) {
		append_bit(1);
		// printf(" R1 ,%d,%d,%d,%d,", (*vars).comp, (*vars).row, (*vars).col, im_data->image[(*vars).comp][(*vars).row][(*vars).col]);
		// printf("\tRUN\n");
	}
}
void decode_run(codingvars* vars, parameters params, image_data* im_data)
{
		if (read_bit() == 1)
		{
			/* F.1 18 a) If R = '1' */

			(*vars).RUNcnt = pow(2, (*vars).J[((*vars).RUNindex)]);
			while ((*vars).RUNcnt > 0)
			{
				(*vars).RUNcnt -= 1;
				im_data->image[(*vars).comp][(*vars).row][(*vars).col] = (*vars).Ra;
				if ((*vars).col == (im_data->width - 1)) {
					if ((*vars).RUNcnt == 0 && (*vars).RUNindex < 31) {
						(*vars).RUNindex += 1;
					}
					return;
				}
					
				else
					(*vars).col += 1;    
			}
			if ((*vars).RUNcnt == 0 && (*vars).RUNindex < 31) {
				(*vars).RUNindex += 1;
			}
			decode_run(vars, params, im_data);
		}
		else
		{
			/* F.1 18 b) If R = '0' */
			(*vars).RUNcnt = read_bits((*vars).J[((*vars).RUNindex)]);
			while ((*vars).RUNcnt > 0)
			{
				(*vars).RUNcnt -= 1;
				im_data->image[(*vars).comp][(*vars).row][(*vars).col] = (*vars).Ra;
				(*vars).col++;
			}
			if ((*vars).RUNindex > 0) {
				(*vars).RUNindex -= 1;
			}

			/* The reverse process of A.7.2 */
			(*vars).Rb = ((*vars).row == 0) ? 0 : im_data->image[(*vars).comp][(*vars).row - 1][(*vars).col];
			if ((*vars).col == 0) {
				(*vars).Ra = (*vars).Rb;
				(*vars).prevRa = (*vars).Ra;
			}
			else
				(*vars).Ra = im_data->image[(*vars).comp][(*vars).row][(*vars).col - 1];

			// index computation
			if (abs((*vars).Ra - (*vars).Rb) <= params.NEAR)
				(*vars).RItype = 1;
			else
				(*vars).RItype = 0;

			// prediction value for a run interruption sample
			if ((*vars).RItype == 1)
				(*vars).Px = (*vars).Ra;
			else
				(*vars).Px = (*vars).Rb;

			// Computation of the auxiliary variable TEMP
			if ((*vars).RItype == 0)
				(*vars).TEMP = (*vars).A[365];
			else
				(*vars).TEMP = (*vars).A[366] + ((*vars).N[366] >> 1);

			// Golomb coding variable computation
			(*vars).Q = (*vars).RItype + 365;
			for ((*vars).k = 0; ((*vars).N[(*vars).Q] << (*vars).k) < (*vars).TEMP; (*vars).k++);

			// Decode EMErrval through Golomb coding
			(*vars).EMErrval = limited_length_Golomb_decode((*vars).k, (*vars).LIMIT - (*vars).J[(*vars).RUNindex] - 1, (*vars).qbpp);

			bool caseA = (*vars).k == 0 && (2 * (*vars).Nn[(*vars).Q - 365] < (*vars).N[(*vars).Q]);
			(*vars).map = ((*vars).EMErrval + (*vars).RItype) % 2;
			if (((*vars).map && caseA) || (!(*vars).map && !caseA))
				(*vars).Errval = ((*vars).EMErrval + (*vars).RItype + (*vars).map) / 2;
			else
				(*vars).Errval = -((*vars).EMErrval + (*vars).RItype + (*vars).map) / 2;


			if (((*vars).RItype == 0) && ((*vars).Ra > (*vars).Rb))
			{
				(*vars).Errval = -(*vars).Errval;
				(*vars).SIGN = -1;
			}
			else
				(*vars).SIGN = 1;

			//Compute the raw pixel
			(*vars).Ix = (*vars).Px + (*vars).Errval;
			im_data->image[(*vars).comp][(*vars).row][(*vars).col] = (*vars).Ix;

			// update of variables for run interruption sample
			if ((*vars).Errval < 0)
				(*vars).Nn[(*vars).Q - 365] = (*vars).Nn[(*vars).Q - 365] + 1;
			(*vars).A[(*vars).Q] = (*vars).A[(*vars).Q] + (((*vars).EMErrval + 1 + (*vars).RItype) >> 1);
			if ((*vars).N[(*vars).Q] == params.RESET)
			{
				(*vars).A[(*vars).Q] = (*vars).A[(*vars).Q] >> 1;
				(*vars).N[(*vars).Q] = (*vars).N[(*vars).Q] >> 1;
				(*vars).Nn[(*vars).Q - 365] = (*vars).Nn[(*vars).Q - 365] >> 1;
			}

	}
}






/* compress ARGB data to tile
*  param:
*    pClrBlk      -- IN, pixel's ARGB data
*    pTile        -- OUT, tile data
*    pTileSize    -- OUT, tile's bytes
*  return:
*    0  -- succeed
*   -1  -- failed
*/

bool judge_lossless(const unsigned char* pClrBlk,const unsigned char* pTile, int nTileSize,image_data* im_raw){

	unsigned char* tempTilea = (uint8*)malloc(nTileSize * sizeof(uint8));
	bool decoing_flag = true;
	init_bitstream_r(pTile, tempTilea,nTileSize);
	// print_bitstream();
	
	parameters params = coding_parameters(decoing_flag);
	codingvars vars;
	image_data* im_data_out = NULL;
	im_data_out=allocate_image();
	
	params.MAXVAL = im_data_out->maxval;

	/* A.2 Initializations and conventions */
	init_codingvars(&vars, params);

	// setting parameters
	if (params.specified_T == false)
	{
		/* C.2.4.1.1 Default threshold values */
		if (params.MAXVAL >= 128)
		{
			vars.FACTOR = floor((float64)(min(params.MAXVAL, 4095) + 128) / 256);
			params.T1 = CLAMP(vars.FACTOR * (vars.BASIC_T1 - 2) + 2 + 3 * params.NEAR, params.NEAR + 1, params.MAXVAL);
			params.T2 = CLAMP(vars.FACTOR * (vars.BASIC_T2 - 3) + 3 + 5 * params.NEAR, params.T1, params.MAXVAL);
			params.T3 = CLAMP(vars.FACTOR * (vars.BASIC_T3 - 4) + 4 + 7 * params.NEAR, params.T2, params.MAXVAL);
		}
		else
		{
			vars.FACTOR = floor(256.0 / (params.MAXVAL + 1));
			params.T1 = CLAMP(max(2, floor((float64)vars.BASIC_T1 / vars.FACTOR) + 3 * params.NEAR), params.NEAR + 1, params.MAXVAL);
			params.T2 = CLAMP(max(2, floor((float64)vars.BASIC_T2 / vars.FACTOR) + 5 * params.NEAR), params.T1, params.MAXVAL);
			params.T3 = CLAMP(max(2, floor((float64)vars.BASIC_T3 / vars.FACTOR) + 7 * params.NEAR), params.T2, params.MAXVAL);
		}
	}


	for (vars.comp = 0; vars.comp < im_data_out->n_comp; vars.comp++)
		for (vars.row = 0; vars.row < im_data_out->height; vars.row++)
			for (vars.col = 0; vars.col < im_data_out->width; vars.col++)
			{
				/* A.3 Context determination */
				
				context_determination(&vars, params, im_data_out);

				if (!vars.RunModeProcessing)
				{
					// regular mode

					/* A.4 Prediction*/
					predict_sample_value(&vars, params);

						// decoding
					decode_prediction_error(&vars, params, im_data_out);

					/* A.6 Update variables */
					update_codingvars(&vars, params);

				}
				else
				{
					// run mode
						// decoding
					decode_run(&vars, params, im_data_out);
				}
			}
	// printf("%d,%d\n",im_data_out->height,im_data_out->width); 
	// for (int i = 0; i < im_data_out->height; i++)
	// {
	// 	for (int j = 0; j < im_data_out->width ; j++)
	// 	{
	// 		if(pClrBlk[(i * im_data_out->width + j)*4]!=im_data_out->image[0][i][j]) return false;
	// 		if(pClrBlk[(i * im_data_out->width + j)*4+1]!=im_data_out->image[1][i][j]) return false;
	// 		if(pClrBlk[(i * im_data_out->width + j)*4+2]!=im_data_out->image[2][i][j]) return false;
	// 		if(pClrBlk[(i * im_data_out->width + j)*4+3]!=im_data_out->image[3][i][j]) return false;
	// 		// printf("%d,",im_data_out->image[0][i][j]);
	// 	}
	// 	// printf("\n");
	// }
	// return true;
	if(compare_image(im_data_out,im_raw)){
		return true;
	}
	else{
		print_image(im_data_out);
		print_image(im_raw);
		return false;
	}
}

int argb2tile(const unsigned char* pClrBlk, unsigned char* pTile, int* pTileSize)
{
	bool decoing_flag = false;
	int MAX_BYTE = 1024;   //表示输出最大的字节数
	// uint8*** pClrBlk_3d;

	parameters params;
	codingvars vars;
	image_data* im_data = NULL;

	// parsing command line parameters and/or JPEG-LS header
	params = coding_parameters(decoing_flag);
	// encoding process

	// loading image data
	im_data = create_image(pClrBlk);

	// bitstream initialization
	unsigned char* tempTile = (uint8*)malloc(MAX_BYTE * sizeof(uint8));
	init_bitstream_w(tempTile);
	params.MAXVAL = im_data->maxval;
	// print_image(im_data);

	/* A.2 Initializations and conventions */   
	init_codingvars(&vars, params);

	// setting parameters
	if (params.specified_T == false)
	{
		/* C.2.4.1.1 Default threshold values */
		if (params.MAXVAL >= 128)
		{
			vars.FACTOR = floor((float64)(min(params.MAXVAL, 4095) + 128) / 256);
			params.T1 = CLAMP(vars.FACTOR * (vars.BASIC_T1 - 2) + 2 + 3 * params.NEAR, params.NEAR + 1, params.MAXVAL);
			params.T2 = CLAMP(vars.FACTOR * (vars.BASIC_T2 - 3) + 3 + 5 * params.NEAR, params.T1, params.MAXVAL);
			params.T3 = CLAMP(vars.FACTOR * (vars.BASIC_T3 - 4) + 4 + 7 * params.NEAR, params.T2, params.MAXVAL);
		}
		else
		{
			vars.FACTOR = floor(256.0 / (params.MAXVAL + 1));
			params.T1 = CLAMP(max(2, floor((float64)vars.BASIC_T1 / vars.FACTOR) + 3 * params.NEAR), params.NEAR + 1, params.MAXVAL);
			params.T2 = CLAMP(max(2, floor((float64)vars.BASIC_T2 / vars.FACTOR) + 5 * params.NEAR), params.T1, params.MAXVAL);
			params.T3 = CLAMP(max(2, floor((float64)vars.BASIC_T3 / vars.FACTOR) + 7 * params.NEAR), params.T2, params.MAXVAL);
		}
	}

	// tileSetSize(im_data->width,im_data->height);

	for (vars.comp = 0; vars.comp < im_data->n_comp; vars.comp++)
		for (vars.row = 0; vars.row < im_data->height; vars.row++)
			for (vars.col = 0; vars.col < im_data->width; vars.col++)
			{
				/* A.3 Context determination */

				context_determination(&vars, params, im_data);

				if (!vars.RunModeProcessing)
				{
					// regular mode

					/* A.4 Prediction*/
					predict_sample_value(&vars, params);
					encode_prediction_error(&vars, params, im_data);

					/* A.6 Update variables */
					update_codingvars(&vars, params);

				}
				else
				{	
						encode_run(&vars, params, im_data);
				}
			}

	deal_EOF();
	// print_bitstream();
	*pTileSize=bs.tot_bytes+1;      //Our calculation is from 0, so it needs add 1
	for(int i=0;i<*pTileSize;i++){
		pTile[i]=bs.pTile[i];
		// printf("%d,",pTile[i]);
	}
	// printf("\n");
	// printf("Bytes:%d\n",*pTileSize);
	

	if(judge_lossless(pClrBlk,tempTile,*pTileSize,im_data)){
		lossless_count++;
		printf("\nLossless!");
	}
	else {
		printf("NOT Lossless!");
		nonlossless_count++;
	}
	printf("Lossless=%d,NOT Lossless=%d\n",lossless_count,nonlossless_count);
	return 0;
}

/* decompress tile data to ARGB
*  param:
*    pTile        -- IN, tile data
*    pTileSize    -- IN, tile's bytes
*    pClrBlk      -- OUT, pixel's ARGB data
*  return:
*    0  -- succeed
*   -1  -- failed
*/

int tile2argb(const unsigned char* pTile, int nTileSize, unsigned char* pClrBlk)
{
	// printf("%d\n",nTileSize);
	unsigned char* tempTile = (uint8*)malloc(nTileSize * sizeof(uint8));
	bool decoing_flag = true;
	init_bitstream_r(pTile, tempTile,nTileSize);
	// print_bitstream();
	
	parameters params = coding_parameters(decoing_flag);
	codingvars vars;
	image_data* im_data_out = NULL;
	im_data_out=allocate_image();
	
	params.MAXVAL = im_data_out->maxval;

	/* A.2 Initializations and conventions */
	init_codingvars(&vars, params);

	// setting parameters
	if (params.specified_T == false)
	{
		/* C.2.4.1.1 Default threshold values */
		if (params.MAXVAL >= 128)
		{
			vars.FACTOR = floor((float64)(min(params.MAXVAL, 4095) + 128) / 256);
			params.T1 = CLAMP(vars.FACTOR * (vars.BASIC_T1 - 2) + 2 + 3 * params.NEAR, params.NEAR + 1, params.MAXVAL);
			params.T2 = CLAMP(vars.FACTOR * (vars.BASIC_T2 - 3) + 3 + 5 * params.NEAR, params.T1, params.MAXVAL);
			params.T3 = CLAMP(vars.FACTOR * (vars.BASIC_T3 - 4) + 4 + 7 * params.NEAR, params.T2, params.MAXVAL);
		}
		else
		{
			vars.FACTOR = floor(256.0 / (params.MAXVAL + 1));
			params.T1 = CLAMP(max(2, floor((float64)vars.BASIC_T1 / vars.FACTOR) + 3 * params.NEAR), params.NEAR + 1, params.MAXVAL);
			params.T2 = CLAMP(max(2, floor((float64)vars.BASIC_T2 / vars.FACTOR) + 5 * params.NEAR), params.T1, params.MAXVAL);
			params.T3 = CLAMP(max(2, floor((float64)vars.BASIC_T3 / vars.FACTOR) + 7 * params.NEAR), params.T2, params.MAXVAL);
		}
	}


	for (vars.comp = 0; vars.comp < im_data_out->n_comp; vars.comp++)
		for (vars.row = 0; vars.row < im_data_out->height; vars.row++)
			for (vars.col = 0; vars.col < im_data_out->width; vars.col++)
			{
				/* A.3 Context determination */
				
				context_determination(&vars, params, im_data_out);

				if (!vars.RunModeProcessing)
				{
					// regular mode

					/* A.4 Prediction*/
					predict_sample_value(&vars, params);

						// decoding
					decode_prediction_error(&vars, params, im_data_out);

					/* A.6 Update variables */
					update_codingvars(&vars, params);

				}
				else
				{
					// run mode
						// decoding
					decode_run(&vars, params, im_data_out);
				}
			}
	// printf("%d,%d\n",im_data_out->height,im_data_out->width); 
	for (int i = 0; i < im_data_out->height; i++)
	{
		for (int j = 0; j < im_data_out->width ; j++)
		{
			pClrBlk[(i * im_data_out->width + j)*4]=im_data_out->image[0][i][j] ;
			pClrBlk[(i * im_data_out->width + j)*4+1]=im_data_out->image[1][i][j] ;
			pClrBlk[(i * im_data_out->width + j)*4+2]=im_data_out->image[2][i][j] ;
			pClrBlk[(i * im_data_out->width + j)*4+3]=im_data_out->image[3][i][j] ;
			// printf("%d,",im_data_out->image[0][i][j]);
		}
		// printf("\n");
	}
	return 0;
}

