#include "stdio.h"
#include "fir.h"

void fir(const float input[], float output[]){
#pragma HLS INTERFACE m_axi port=input bundle=gmem0 offset=slave depth=1024
#pragma HLS INTERFACE m_axi port=output bundle=gmem1 offset=slave depth=1024
#pragma HLS INTERFACE s_axilite port=input bundle=control
#pragma HLS INTERFACE s_axilite port=output bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

	static float shift_reg[NUM_TAPS];
	#pragma HLS ARRAY_PARTITION variable=shift_reg complete dim=0

	for(int j = 0; j < SIZE; j ++ ) {
		#pragma HLS pipeline II=1

		float acc = 0;
		for (int i = NUM_TAPS - 1; i > 0; i--) {
			shift_reg[i] = shift_reg[i - 1];
			acc += shift_reg[i] * taps[i];
		}

		acc += input[j] * taps[0];
		shift_reg[0] = input[j];
		output[j] = acc;
	}
}

// void fir_fixed(const float input[], float output[]){
// #pragma HLS INTERFACE m_axi port=input bundle=gmem0 offset=slave depth=1024
// #pragma HLS INTERFACE m_axi port=output bundle=gmem1 offset=slave depth=1024
// #pragma HLS INTERFACE s_axilite port=input bundle=control
// #pragma HLS INTERFACE s_axilite port=output bundle=control
// #pragma HLS INTERFACE s_axilite port=return bundle=control

// 	static ap_fixed<W,I> shift_reg_fixed[NUM_TAPS];
// 	#pragma HLS ARRAY_PARTITION variable=shift_reg_fixed complete dim=0

// 	for(int j = 0; j < SIZE; j ++ ) {
// 		#pragma HLS pipeline II=1

// 		ap_fixed<W,I> acc = 0;
// 		for (int i = NUM_TAPS - 1; i >= 0; i--) {
// 			shift_reg_fixed[i] = shift_reg_fixed[i - 1];
// 			acc += shift_reg_fixed[i] * taps_fixed[i];
// 		}

// 		acc += (ap_fixed<W,I>)input[j] * taps_fixed[0];
// 		shift_reg_fixed[0] = (ap_fixed<W,I>)input[j];
// 		output[j] = (float)acc;
// 	}
// }
