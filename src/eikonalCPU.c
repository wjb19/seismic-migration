// host code for CUDA eikonal solver
#include <math.h>
#include <stdio.h>
#include <omp.h>        		

void eikUpdateCPU(float* ttHeap,
	float* ttKeys,
	float* radius,
	float* ttime,
	int sizeTheta,
	int sizePhi){



	minHeapKeys(ttHeap, ttKeys, sizeTheta*sizePhi);

	int index = ttKeys[0];
	printf("index : %i\n",index);
	radius[index] += 1.0f;
	ttime[index] = ttHeap[0];

}


void eikSolverCPU(float *ttInput,
   	float *ttOutput,
   	float *ttHeap,
	float *keys,
        float *radius,
        float *slowness,
        int sizeR,
        float deltaR,
        float startR,
        int sizeTheta,
        float deltaTheta,
        float startTheta,
        int sizePhi,
        float deltaPhi,
        float startPhi){

	float theta = startTheta;
        float delRSqInv         = 1 / (deltaR * deltaR);
	#pragma omp parallel for
	for (int i=0; i<sizeTheta; i++){
        	theta     += deltaTheta;
		
		for (int j=0; j<sizePhi; j++){
			
			float secondOrder=0;
			float firstOrder=0;
			float zerothOrder=0;

			float slowSq    	= slowness[j+i*sizePhi];
        		float rad    		= radius[j+i*sizePhi] * deltaR;
        		float delTheRadInv      = 1 / (rad * deltaTheta);
        		float delThePhiRadInv   = 1 / (rad * deltaPhi * sin(theta));

			//printf("%e %e %e %e %e\n",theta,slowSq,rad,delTheRadInv,delThePhiRadInv);

        		//solve
        		secondOrder    +=      delRSqInv;
        		secondOrder    +=      delTheRadInv;
        		secondOrder    +=      delThePhiRadInv;
        
			firstOrder      += (j < sizePhi) ? - 2.0f * delRSqInv * ttInput[j+1+i*sizePhi] :
                -2.0f * delRSqInv * ttInput[i*sizePhi];

        		firstOrder      += (j < sizePhi) ? - 2.0f * delTheRadInv * ttInput[j+1+i*sizePhi] :
                -2.0f * delTheRadInv * ttInput[i*sizePhi];

        		firstOrder      += (j < sizePhi) ? - 2.0f * delThePhiRadInv * ttInput[j+1+i*sizePhi] :
                -2.0f * delThePhiRadInv * ttInput[i*sizePhi];

        		zerothOrder     += (i > 0) ? delRSqInv * ttInput[j+(i-1)*sizePhi] :
                delRSqInv * ttInput[j+(sizeTheta-1)*sizePhi];

        		zerothOrder     += (j > 0) ? delTheRadInv * ttInput[j-1+i*sizePhi] :
                delTheRadInv * ttInput[(i+1)*sizePhi-1];

        		zerothOrder     += (j > 0) ? delThePhiRadInv * ttInput[j-1+i*sizePhi] :
                delThePhiRadInv * ttInput[(i+1)*sizePhi-1];

        		zerothOrder     += - slowSq;

        		ttHeap[j+i*sizePhi]     = - firstOrder + sqrt(firstOrder * firstOrder - 4 * secondOrder * zerothOrder) / ( 2.0f * secondOrder);
        		ttOutput[j+i*sizePhi]     = - firstOrder + sqrt(firstOrder * firstOrder - 4 * secondOrder * zerothOrder) / ( 2.0f * secondOrder);
			keys[j+i*sizePhi] = j+ i*sizePhi;
			
		}
	}
}

void kdmEikonalSolverCPU(struct procGrid* waveRadius, 
	struct procGrid* waveTheta, 
	struct procGrid* wavePhi,  
	float* slowness, 
	float* ttime,
	int solveSteps){
	
	// grab and copy scalars from structs

	int WAVE_R_SIZE 	= waveRadius->procPnts;
	float WAVE_R_STEP 	= waveRadius->procStep;
	float WAVE_R_START 	= waveRadius->procStart;

	int WAVE_THETA_SIZE 	= waveTheta->procPnts;
	float WAVE_THETA_STEP 	= waveTheta->procStep;
	float WAVE_THETA_START 	= waveTheta->procStart;
	
	int WAVE_PHI_SIZE 	= wavePhi->procPnts;
	float WAVE_PHI_STEP 	= wavePhi->procStep;
	float WAVE_PHI_START 	= wavePhi->procStart;

	int memSizeAliveTable 	= WAVE_PHI_SIZE*WAVE_THETA_SIZE*sizeof(float);
	int memSizeTTime 	= WAVE_PHI_SIZE*WAVE_THETA_SIZE*WAVE_R_SIZE*sizeof(float);
	int memSizeSlowness 	= memSizeTTime;
	int memSizeRadius	= memSizeAliveTable;


        float* slownessS        = (float*) malloc(memSizeSlowness);
        float* aliveTableInput 	= (float*) malloc(memSizeAliveTable);
        float* aliveTableOutput = (float*) malloc(memSizeAliveTable);
        float* aliveTableHeap 	= (float*) malloc(memSizeAliveTable);
        float* aliveTableKeys 	= (float*) malloc(memSizeAliveTable);
        float* radius           = (float*) malloc(memSizeRadius);



        //init
        for (int i=0; i<WAVE_R_SIZE*WAVE_THETA_SIZE*WAVE_PHI_SIZE; i++){
                slownessS[i]    = slowness[i] * slowness[i];
                ttime[i] = 0.0f;
		
        }

        //init with analytic solution close to origin
        for (int i=0; i<WAVE_THETA_SIZE*WAVE_PHI_SIZE; i++){
                aliveTableInput[i]   = slowness[i] * 0.00001f;       //zeroth radius
                radius[i] = 1.0f;
        }


        for (int i=0; i<solveSteps/2; i++){

                //solve
                eikSolverCPU(aliveTableInput,
                        aliveTableOutput,
			aliveTableHeap,
			aliveTableKeys,
                        radius,
                        slownessS,
                        WAVE_R_SIZE,
                        WAVE_R_STEP,
                        WAVE_R_START,
                        WAVE_THETA_SIZE,
                        WAVE_THETA_STEP,
                        WAVE_THETA_START,
                        WAVE_PHI_SIZE,
                        WAVE_PHI_STEP,
                        WAVE_PHI_START);


		eikUpdateCPU(aliveTableHeap,
			aliveTableKeys,
			radius,
			ttime,
			WAVE_THETA_SIZE,
			WAVE_PHI_SIZE);

         	eikSolverCPU(aliveTableOutput,
                        aliveTableInput,
			aliveTableHeap,
			aliveTableKeys,
                        radius,
                        slownessS,
                        WAVE_R_SIZE,
                        WAVE_R_STEP,
                        WAVE_R_START,
                        WAVE_THETA_SIZE,
                        WAVE_THETA_STEP,
                        WAVE_THETA_START,
                        WAVE_PHI_SIZE,
                        WAVE_PHI_STEP,
                        WAVE_PHI_START);


		eikUpdateCPU(aliveTableHeap,
			aliveTableKeys,
			radius,
			ttime,
			WAVE_THETA_SIZE,
			WAVE_PHI_SIZE);


	}

	

        free(slownessS);
        free(aliveTableInput);
        free(aliveTableOutput);
        free(aliveTableHeap);
        free(aliveTableKeys);
        free(radius);




}
