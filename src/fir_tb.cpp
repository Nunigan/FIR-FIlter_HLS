#include "stdio.h"
#include "fir.h"
#include <math.h>       /* fabs */

void fir_fixed(const float input[], float output[]);
void fir(const float input[], float output[]);

int main() {
	
	float out[SIZE];
	float out_fixed[SIZE];

	fir(x, out);

	for (int i = 0; i < SIZE; i++) {
		if (fabs(out[i]-filtered_x[i]) > 0.000001)
			std::cout << "float fir wrong " << i << ": " << out[i] << " ,"<< filtered_x[i] << std::endl;
	}

	fir_fixed(x, out_fixed);

		for (int i = 0; i < SIZE; i++) {
		if (fabs(out[i]-out_fixed[i]) > 0.000001)
			std::cout << "fixed fir wrong " << i << ": " << out[i] << " ,"<< filtered_x[i] << std::endl;
	}
}



