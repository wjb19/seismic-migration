// CPU code for psktm
//
// loop in time/pseudo-depth, trace records
//
// wjb 11/09, 12/09, 06/10


#include <stdio.h>
#include <math.h>
#include <omp.h>


void ktmMigrationCPU(struct imageGrid* imageX, 
	struct imageGrid* imageY, 
	struct imageGrid* imageZ, 
	struct jobParams* config, 
	float* midX, 
	float* midY, 
	float* offX, 
	float* offY, 
	float* traces, 
	float* slowness, 
	float* image){
	
	int timeIndex;
	int Li		=imageX->imagePnts;
	int Lj		=imageY->imagePnts;
	int Ll		=imageZ->imagePnts;
	float sRate	=config->traceSampleRate;	

	//malloc for image axes
	
	int imXYZ = Li*Lj*Ll;
	float* imageXX = (float*) malloc(Li*sizeof(float));	
	float* imageYY = (float*) malloc(Lj*sizeof(float));
	float* imageZZ = (float*) malloc(Ll*sizeof(float));
	

	float travelTime, weight, temp, tempA, tempB, tempC, tempD; 	
	initThreeAxes(imageX, imageY, imageZ, imageXX, imageYY, imageZZ);

	float* tauS = (float*) malloc(sizeof(float)*Ll);
	float* slownessS = (float*) malloc(sizeof(float)*Ll);

	//precompute squares
	for (int i=0; i<Ll; i++){
		tauS[i] = imageZZ[i] * imageZZ[i];
		slownessS[i] = 1.0f / (slowness[i]*slowness[i]);
	}

	//clear image mem
	for(int i=0; i<imXYZ; i++)
		image[i]=0.0f;


	#pragma omp parallel for shared(image, traces,Li,Lj,Ll,config,imageXX,imageYY,midX,midY,offX,offY,tauS,slownessS,sRate) private(temp,tempA,tempB,tempC,tempD,timeIndex)
	
	//loop over trace records
	for (int k=0; k<config->traceNo; k++){
		
			//loop over imageX
	
			for(int i=0; i<Li; i++){

				tempC = ( midX[k] - imageXX[i]-offX[k]) * (midX[k]- imageXX[i]-offX[k]);
				tempD = ( midX[k] - imageXX[i]+offX[k]) * (midX[k]- imageXX[i]+offX[k]);

				//loop over imageY
				for(int j=0; j<Lj; j++){
				
					tempA = tempC + ( midY[k] - imageYY[j]-offY[k]) * (midY[k]- imageYY[j]-offY[k]);
					tempB = tempD + ( midY[k] - imageYY[j]+offY[k]) * (midY[k]- imageYY[j]+offY[k]);
			
        
					//loop over imageZ                		
					for (int l=0; l<Ll; l++){

						temp = sqrtf(tauS[l] + tempA * slownessS[l]); 
						temp += sqrtf(tauS[l] + tempB * slownessS[l]); 
						timeIndex = (int) (temp / sRate);

                        			if ((timeIndex < config->tracePts) && (timeIndex > 0)){

                                			image[i*Lj*Ll + j*Ll + l] += 
								traces[timeIndex + k * config->tracePts] * temp *sqrtf(tauS[l] / temp);
						}
			
					} //imageZ
				} //imageY
			} //imageX
		}//input traces

	free(imageXX);
	free(imageYY);	
	free(imageZZ);
	free(tauS);
	free(slownessS);
	return;
	
	
} //end function 

