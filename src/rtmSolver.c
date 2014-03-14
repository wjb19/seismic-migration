//test kernel for rtm2d
//take some trace data and migrate back in time
//WJB 03/11

#include<stdio.h>
#include "rtmHeader.h"
#include<math.h>
#include<fenv.h>

void reverseTime2d(float * data,
	float ** model,
	float ** velocity,
	struct rtmGrid2d* inputVariables){


	int horizSamples 	= inputVariables->numHorizSample;
	int depthSamples 	= inputVariables->numDepthSample;
	int timeSamples 	= inputVariables->numTimeSample;
	float deltaHoriz 	= inputVariables->deltaHorizSample;
	float deltaDepth 	= inputVariables->deltaDepthSample;
	float deltaTime 	= inputVariables->deltaTimeSample;
	
	double ***finiteDifferences = (double***) malloc(sizeof(double**) * depthSamples);
	
	for (int i=0; i<depthSamples; i++)
		finiteDifferences[i] = (double**) malloc(sizeof(double*) * horizSamples);
	
	for (int i=0; i<depthSamples; i++)
		for (int j=0; j<horizSamples; j++)
			finiteDifferences[i][j] = (double*) malloc(sizeof(double) * 3);

	double **waveCoefficients = (double**) malloc(sizeof(double*) * depthSamples);	

	for (int i=0; i<depthSamples; i++)
		waveCoefficients[i] = (double*) malloc(sizeof(double) * horizSamples);

	//not quite zero to avoid fpe
	for (int i=0; i<depthSamples; i++){
		for (int j=0; j<horizSamples; j++){
			for (int k=0; k<3; k++){
				finiteDifferences[i][j][k]=1e-32f;
			}
		}
	}

	//absorbing boundary
	
	double* boundary = (double*) malloc(sizeof(double) * 20);
	for (int i=1; i<21; i++){
		double kern = (0.015f*(20.0f-(double) i));
		double kernal = - kern * kern;
		boundary[i-1] = pow(exp( kernal ),10.0f);
	}

	//init
	for (int i=0; i<horizSamples; i++){
		finiteDifferences[0][i][0] = data[i * timeSamples + timeSamples-1];
		finiteDifferences[0][i][1] = data[i * timeSamples + timeSamples-2];
		finiteDifferences[0][i][2] = data[i * timeSamples + timeSamples-3];
	}


	//coefficients
	for (int i=0; i<depthSamples; i++){
		for (int j=0; j<horizSamples; j++){
			waveCoefficients[i][j]= (double) (velocity[i][j] * deltaTime / deltaDepth);
			waveCoefficients[i][j]*= (double) (velocity[i][j] * deltaTime / deltaDepth);
		}
	}
	
	int cz=1;
	int bz=0;

	//int feenableexcept();
        //feenableexcept(FE_OVERFLOW);
	//feenableexcept(FE_UNDERFLOW);
	
	//race problems
	//#pragma omp parallel for shared(finiteDifferences,waveCoefficients,model,data)

	for (int ii=timeSamples-1; ii>0; --ii){

		cz+=1;
		bz=min(cz,depthSamples-1);
	
		#pragma omp parallel for shared(finiteDifferences)
	
		//apply absorbing boundary on left/right sides
		for (int j=0; j<bz; ++j){
			for (int k=0; k<20; ++k){
				finiteDifferences[j][k][0] *= boundary[19-k];
				finiteDifferences[j][k][1] *= boundary[19-k];
				finiteDifferences[j][horizSamples -20 + k][0] *= boundary[k];
				finiteDifferences[j][horizSamples -20 + k][1] *= boundary[k];
			}
		}

	
		
		//apply absorbing boundary at depth nz
		if (bz >= depthSamples-19){
		#pragma omp parallel for shared(finiteDifferences)
			for (int j=depthSamples-19; j<bz; ++j){
				for (int k=0; k<horizSamples; ++k){
                                	finiteDifferences[j][k][0] *= boundary[depthSamples-j];
                                	finiteDifferences[j][k][1] *= boundary[depthSamples-j];
                        	}
			}
                }

		//computing grid depth (extend in z to solve)

		int ez= (bz==depthSamples-1) ? depthSamples-2 : bz;
	
		//time extrapolation 

		#pragma omp parallel for shared(finiteDifferences)
		for (int i=0; i<bz; i++){
			for (int j=1; j<horizSamples-2; j++){
				finiteDifferences[i][j][2] -= finiteDifferences[i][j][0];
			}
		}
		
		float b;

		#pragma omp parallel for shared(finiteDifferences)
		for (int i=1; i<ez; i++){
			for (int j=1; j<horizSamples-2; j++){
				b = 2.0f - 4.0f*waveCoefficients[i][j];
				finiteDifferences[i][j][1] += (b*finiteDifferences[i][j][0] +
					waveCoefficients[i][j+1]*finiteDifferences[i][j+1][0] +
					waveCoefficients[i][j-1]*finiteDifferences[i][j-1][0] +
					waveCoefficients[i+1][j]*finiteDifferences[i+1][j][0] +
					waveCoefficients[i-1][j]*finiteDifferences[i-1][j][0]);
			
				//printf("%i %i %i %f %f\n",ii,i,j,b,finiteDifferences[i][j][1]);
			}
		}

		#pragma omp parallel for shared(finiteDifferences)
		for (int j=1; j<horizSamples-2; j++){
			b = 2.0f - 4.0f*waveCoefficients[0][j];
			finiteDifferences[0][j][1] += (b*finiteDifferences[0][j][0] +
				waveCoefficients[0][j+1]*finiteDifferences[0][j+1][0] +
				waveCoefficients[0][j-1]*finiteDifferences[0][j-1][0] +
				waveCoefficients[1][j]*finiteDifferences[1][j][0]);

		}

		//extrapolation at depthSamples
		if (bz==depthSamples-1){

		#pragma omp parallel for shared(finiteDifferences)
			for (int j=1; j<horizSamples-2; j++){
				b = 2.0f - 4.0f*waveCoefficients[bz][j];
				finiteDifferences[bz][j][1] += (b*finiteDifferences[bz][j][0] +
					waveCoefficients[bz][j+1]*finiteDifferences[bz][j+1][0] +
					waveCoefficients[bz][j-1]*finiteDifferences[bz][j-1][0] +
					waveCoefficients[bz-1][j]*finiteDifferences[bz-1][j][0]);
			}
		

		b = 2.0f - 4.0f*waveCoefficients[bz][0];
		finiteDifferences[bz][0][1] += (b*finiteDifferences[bz][0][0] +
			waveCoefficients[bz][1]*finiteDifferences[bz][1][0] + 
			waveCoefficients[bz-1][0]*finiteDifferences[bz-1][0][0]);
			
		}

	
		#pragma omp parallel for shared(finiteDifferences)
		for (int i=1; i<ez; i++){
			b = 2.0f - 4.0f*waveCoefficients[i][0];
			finiteDifferences[i][0][1] += (b*finiteDifferences[i][0][0] +
				waveCoefficients[i][1]*finiteDifferences[i][1][0] +
				waveCoefficients[i+1][0]*finiteDifferences[i+1][0][0] +
				waveCoefficients[i-1][0]*finiteDifferences[i-1][0][0]);
		}

		int ss = horizSamples-1;

		
		#pragma omp parallel for shared(finiteDifferences)
		for (int i=1; i<ez; i++){
			b = 2.0f - 4.0f*waveCoefficients[i][ss];
			finiteDifferences[i][ss][1] += (b*finiteDifferences[i][ss][0] +
				waveCoefficients[i][ss-1]*finiteDifferences[i][ss-1][0] +
				waveCoefficients[i+1][ss]*finiteDifferences[i+1][ss][0] +
				waveCoefficients[i-1][ss]*finiteDifferences[i-1][ss][0]);
		}
		
		//corners
		
		b = 2.0f - 4.0f*waveCoefficients[0][0];
		finiteDifferences[0][0][1] += (b*finiteDifferences[0][0][0] +
			waveCoefficients[0][1]*finiteDifferences[0][1][0] +
			waveCoefficients[1][0]*finiteDifferences[1][0][0]);
		

		b = 2.0f - 4.0f*waveCoefficients[0][ss];
		finiteDifferences[0][ss][1] = (b*finiteDifferences[0][ss][0] +
			waveCoefficients[0][ss-1]*finiteDifferences[0][ss-1][0] +
			waveCoefficients[1][ss]*finiteDifferences[1][ss][0]);


		//init for next

		#pragma omp parallel for shared(finiteDifferences)
		for (int i=0; i<depthSamples; i++){
			for (int j=0; j<horizSamples; j++){
				finiteDifferences[i][j][0] = finiteDifferences[i][j][1];
				finiteDifferences[i][j][1] = finiteDifferences[i][j][2];
			}
		}

		//surface boundary wavefield

		if (ii>1){

	
		#pragma omp parallel for shared(finiteDifferences)
			for (int i=1; i<depthSamples; i++){
				for (int j=0; j<horizSamples; j++){
					finiteDifferences[i][j][2] = 0.0f;
				}
			}
			
		#pragma omp parallel for shared(finiteDifferences)
			for (int j=0; j<horizSamples; j++){
				finiteDifferences[0][j][2] = data[j * timeSamples + ii - 1];
			}
		}
	}

		#pragma omp parallel for shared(finiteDifferences)
	for (int i=0; i<depthSamples; i++){
		for (int j=0; j<horizSamples; j++){
			model[i][j]=finiteDifferences[i][j][0];
		}
	}

	free(finiteDifferences);
	free(waveCoefficients);
}

