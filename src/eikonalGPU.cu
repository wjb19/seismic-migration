// host code for CUDA eikonal solver
//wjb 03/11

#include <cuda.h>
#include <stdio.h>



extern "C"{
void kdmEikonalSolverGPU(struct procGrid* waveRadius, 
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

	int memSizeAliveTable	= sizeof(float) * WAVE_THETA_SIZE*WAVE_PHI_SIZE;
	int memSizeSlowness	= sizeof(float) * WAVE_R_SIZE*WAVE_THETA_SIZE*WAVE_PHI_SIZE;
	int memSizeRadius	= memSizeAliveTable;
	int memSizeTTable	= memSizeSlowness;
	int memSizeReduction2k	= sizeof(float) * 2 * 1024;
	//ship off variables to cmem

	//cudaMemcpyToSymbol("deltaTheta",&WAVE_THETA_STEP,sizeof(float));
	//cudaMemcpyToSymbol("deltaPhi",&WAVE_PHI_STEP,sizeof(float));
	//cudaMemcpyToSymbol("deltaR",&WAVE_R_STEP,sizeof(float));
	//cudaMemcpyToSymbol("sizeR",&WAVE_R_SIZE,sizeof(int));
	//cudaMemcpyToSymbol("sizeTheta",&WAVE_THETA_SIZE,sizeof(int));
	//cudaMemcpyToSymbol("sizePhi",&WAVE_PHI_SIZE,sizeof(int));
	//cudaMemcpyToSymbol("startR",&WAVE_R_START,sizeof(float));
	//cudaMemcpyToSymbol("startTheta",&WAVE_THETA_START,sizeof(float));
	//cudaMemcpyToSymbol("startPhi",&WAVE_PHI_START,sizeof(float));
	
	float* slownessS	= (float*) malloc(memSizeSlowness);
	float* aliveTable	= (float*) malloc(memSizeAliveTable);
	float* radius		= (float*) malloc(memSizeRadius);
	float* reduction	= (float*) malloc(4 * sizeof(float));

		

	//init
	for (int i=0; i<WAVE_R_SIZE*WAVE_THETA_SIZE*WAVE_PHI_SIZE; i++){
		slownessS[i]	= slowness[i] * slowness[i];
		ttime[i] = 0.0f;
	}
	
	//init with analytic solution close to origin
	for (int i=0; i<WAVE_THETA_SIZE*WAVE_PHI_SIZE; i++){
		aliveTable[i]	= slowness[i] * 0.00001f;	//zeroth radius
		radius[i] = 0.0f;
	}


        float *devAliveTableInput, *devAliveTableOutput, *devSlowness, *devRadius, *devTTable;
	float *dev2kReductionKeys, *dev2kReductionValues, *dev4ReductionKeys, *dev4ReductionValues;
	

        cudaMalloc((void**) &devAliveTableInput,     memSizeAliveTable);
        cudaMalloc((void**) &devAliveTableOutput,     memSizeAliveTable);
        cudaMalloc((void**) &devSlowness,     	memSizeSlowness);
        cudaMalloc((void**) &devRadius, 	memSizeRadius);
	cudaMalloc((void**) &devTTable,		memSizeTTable);

	cudaMalloc((void**) &dev2kReductionKeys, memSizeReduction2k);
	cudaMalloc((void**) &dev2kReductionValues, memSizeReduction2k);
	cudaMalloc((void**) &dev4ReductionKeys, sizeof(float) * 4);
	cudaMalloc((void**) &dev4ReductionValues, sizeof(float) * 4);
	
	//init reduction
	for (int i=0; i<4; i++) reduction[i]=0.0f;

	// Load device memory
        cudaMemcpy(devAliveTableInput, 	aliveTable, 	memSizeAliveTable,	cudaMemcpyHostToDevice);
	cudaMemcpy(devAliveTableOutput, aliveTable, 	memSizeAliveTable,	cudaMemcpyHostToDevice);       	
	cudaMemcpy(devSlowness, 	slownessS, 	memSizeSlowness,	cudaMemcpyHostToDevice);
       	cudaMemcpy(devRadius,		radius,   	memSizeRadius,		cudaMemcpyHostToDevice);
        cudaMemcpy(devTTable, 		ttime, 		memSizeTTable,		cudaMemcpyHostToDevice);
	cudaMemcpy(dev4ReductionKeys, 	reduction, 	sizeof(float)*4,	cudaMemcpyHostToDevice);

	dim3 solverThreads,solverBlocks;
	solverThreads.x=EIK_SOLVER_THREADS_X; 		
	solverThreads.y=EIK_SOLVER_THREADS_Y;
	solverBlocks.x=EIK_SOLVER_BLOCK_X; 		
	solverBlocks.y=EIK_SOLVER_BLOCK_Y;

	dim3 reductionThreads;
	reductionThreads.x=512;
	reductionThreads.y=1;

	dim3 reduction2kBlocks;
	reduction2kBlocks.x=2*1024;
	reduction2kBlocks.y=1;

	dim3 reduction4Blocks;
	reduction4Blocks.x=4;
	reduction4Blocks.y=4;

	for (int i=0; i<solveSteps/2; i++){
		
		//solve
		eikSolverCUDA<<<solverBlocks, solverThreads>>>(devAliveTableInput,
			devAliveTableOutput,
			devRadius,
			devSlowness,
			WAVE_R_SIZE,
			WAVE_R_STEP,
			WAVE_R_START,
			WAVE_THETA_SIZE,
			WAVE_THETA_STEP,
			WAVE_THETA_START,
			WAVE_PHI_SIZE,
			WAVE_PHI_STEP,
			WAVE_PHI_START);

		//advance wavefront: find global min	
		//spread across many MP
		reductionGlobalMin<<<reduction2kBlocks, reductionThreads>>>(devAliveTableOutput, 
			dev2kReductionKeys, 
			dev2kReductionValues,
			WAVE_THETA_SIZE,
			WAVE_PHI_SIZE);
	
		//fits on single MP: find local min & increment radius index in mesh
		//& write out a ttime
		advanceWave<<<reduction4Blocks, reductionThreads>>>(dev2kReductionKeys,
			dev2kReductionValues,
			dev4ReductionKeys,
			dev4ReductionValues,
			devRadius,
			devTTable,
			WAVE_THETA_SIZE,
			WAVE_PHI_SIZE);

		//solve, having swapped arrays
		eikSolverCUDA<<<solverBlocks, solverThreads>>>(devAliveTableOutput,
			devAliveTableInput,
			devRadius,
			devSlowness,
			WAVE_R_SIZE,
			WAVE_R_STEP,
			WAVE_R_START,
			WAVE_THETA_SIZE,
			WAVE_THETA_STEP,
			WAVE_THETA_START,
			WAVE_PHI_SIZE,
			WAVE_PHI_STEP,
			WAVE_PHI_START);

		//~ global min
		reductionGlobalMin<<<reduction2kBlocks, reductionThreads>>>(devAliveTableInput, 
			dev2kReductionKeys, 
			dev2kReductionValues,
			WAVE_THETA_SIZE,
			WAVE_PHI_SIZE);
	
		//local min & increment radius index in mesh
		//& write out a ttime
		advanceWave<<<reduction4Blocks, reductionThreads>>>(dev2kReductionKeys,
			dev2kReductionValues,
			dev4ReductionKeys,
			dev4ReductionValues,
			devRadius,
			devTTable,
			WAVE_THETA_SIZE,
			WAVE_PHI_SIZE);

		
		//printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	}


	//transfer back to host
	cudaMemcpy(ttime,
		devTTable,
		memSizeTTable,           
		cudaMemcpyDeviceToHost);


	//cleanup

	cudaFree(devAliveTableInput);
	cudaFree(devAliveTableOutput);
	cudaFree(devTTable);
	cudaFree(devSlowness);
	cudaFree(devRadius);
	cudaFree(dev2kReductionKeys);
	cudaFree(dev2kReductionValues);
	cudaFree(dev4ReductionKeys);
	cudaFree(dev4ReductionValues);


	free(radius);
	free(slownessS);
	free(aliveTable);

} //end host function 
} //end extern
