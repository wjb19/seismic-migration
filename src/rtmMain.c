//driver for rtm 
//wjb 03/11

#include <stdio.h>
#include <math.h>
#include <rtmHeader.h>
#include <time.h>
#include <omp.h>


int main()
{
	struct rtmGrid2d input;
	input.numHorizSample=240;
	input.numDepthSample=2400;
	input.numTimeSample=2400;
	input.deltaHorizSample=10.0f;
	input.deltaDepthSample=10.0f;
	input.deltaTimeSample=0.002f;


	float** vel = (float**) malloc (sizeof(float*) * input.numDepthSample);
	float** mod = (float**) malloc (sizeof(float*) * input.numDepthSample);
	float* dat = (float*) malloc (sizeof(float) * input.numTimeSample * input.numHorizSample);

	for (int i=0; i<input.numDepthSample; i++){
		vel[i] = (float*) malloc (sizeof(float) * input.numHorizSample);
		mod[i] = (float*) malloc (sizeof(float) * input.numHorizSample);
	}
	

	for (int i=0; i<input.numDepthSample; i++){
		for (int j=0; j<input.numHorizSample; j++){
			vel[i][j] = 1500.0f;
		}
	}
	

	for (int i=0; i<input.numTimeSample * input.numHorizSample; i++)
		dat[i] = 1e-32f;

	dat[input.numTimeSample * input.numHorizSample / 2] = 1.0f;
	dat[input.numTimeSample * input.numHorizSample / 4] = 1.0f;
	dat[3*input.numTimeSample * input.numHorizSample / 4] = 1.0f;

	//wavelet
	
	//fm
	//forwardTime2d(mod,dat,vel,&input);

	//rm
	
	reverseTime2d(dat, mod, vel, &input);

	free(vel);
	free(mod);
	free(dat);

	for (int i=0; i<input.numDepthSample; i++)
		for (int j=0; j<input.numHorizSample; j++)
			printf("%f\n",mod[i][j]);

	return 0;

}
