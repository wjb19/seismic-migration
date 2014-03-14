// CUDA for psktm
// grid dims:
//
// 64k blocks (max)
// (128,1,1) threads / block
//
// x thread -> trace x/y
// block -> one image x/y
// loop in image time/pseudo-depth
// 
//
// wjb 11/09, 12/09, 06/10


#include <cuda.h>
#include <stdio.h>


//data variables/vectors

__device__ 	__constant__ float devSlownessS[2048];
__device__ 	__constant__ float devTauS[2048];

__device__ 	__constant__ float devImageXStart;
__device__ 	__constant__ float devImageYStart;

__device__ 	__constant__ float devImageXStep;
__device__	__constant__ int devImageXSize;
__device__	__constant__ float devImageYStep;
__device__ 	__constant__ int devImageYSize;

__device__  	__constant__ int devTracePts;
__device__  	__constant__ int devTraceNo;
__device__ 	__constant__ float devTraceSampleRate;
__device__	__constant__ int devImageDepths;

texture<float,1,cudaReadModeElementType> tracesTexture;

//calculate min value in shared array
//unused 06/10

__device__ int min(int *array){

	int min=devTracePts;
	for (int i=0; i<THREADS_X; i++) {if (array[i]<min)  {min=array[i];}}

	return min;
}


//calculate max value in shared array
//unused 06/10
__device__ int max(int *array){
	
	int max=0.0; 
	for (int i=0; i<THREADS_X; i++) {if ((array[i]>max) & (array[i]<devTracePts)) {max=array[i];}}


     	return max;
}

// calculate travelTime on device

__device__ float initTravelTime(float spaceS, float spaceR, int depth){
	
	float travelTime 	= 0.0f;
	travelTime += sqrtf(devTauS[depth]  + spaceS* devSlownessS[depth]);
	travelTime += sqrtf(devTauS[depth]  + spaceR* devSlownessS[depth]);
	
	return travelTime;
}
// calculate weight on device

__device__ float initWeight(float travelTime, int depth)
{

	float weight	=0.0f;
	weight	= travelTime * sqrtf(devTauS[depth] / travelTime);
	
	return weight;

}

//CUDA kernel

__global__ void ktmCUDA(float* devMidX, 
	float* devMidY, 
	float* devOffX, 
	float* devOffY, 
	float* devTraces, 
	float* devImage){

	int offI, offJ, timeIndex;

	const int traceBlocks 		= devTraceNo / THREADS_X;
	const int traceCol		= THREADS_X * devTracePts;
	const int imageBlocks		= devImageDepths / THREADS_Y;
	const int thread_x		= threadIdx.x;
	const int thread_y		= threadIdx.y;
	const int block_x		= blockIdx.x;


	float travelTime, weight, traceBlock;

	__shared__ int flag;
	__shared__ float imageBlock[THREADS_X][THREADS_Y];

	float imageXBlock = devImageXStart +  ((float)(block_x / devImageYSize)* devImageXStep); 
	float imageYBlock = devImageYStart + ((float)(block_x % devImageYSize) * devImageYStep); 
	float tempA, tempB, spaceS, spaceR;

	//loop over input trace space
	for (int i=0; i<traceBlocks; i++){

		offI = thread_x + THREADS_X*i; 
	
		tempA = devMidX[offI] - imageXBlock;
		tempB = devMidY[offI] - imageYBlock;

		spaceS = (tempA - devOffX[offI])*(tempA - devOffX[offI]);
		spaceS += (tempB - devOffY[offI])*(tempB - devOffY[offI]);
       	 
		spaceR = (tempA + devOffX[offI])*(tempA + devOffX[offI]);
		spaceR +=(tempB + devOffY[offI])*(tempB + devOffY[offI]);

		//loop over image time/depth
		for (int j=0; j<imageBlocks; j++){

			offJ = thread_y + THREADS_Y*j;
			
			//calculate ttime
			travelTime=initTravelTime(spaceS,spaceR,offJ);
			timeIndex = (int) (travelTime / devTraceSampleRate); 
			

			//if (__any((timeIndex < devTracePts) & (timeIndex > 0))){
			//	flag = 1;
			//}

			flag=1;
			__syncthreads();

			if (flag!=0){
				imageBlock[thread_x][thread_y]=0.0f;

				if ((timeIndex < devTracePts) & (timeIndex > 0)) {
					weight     = initWeight(travelTime, offJ);
					traceBlock = devTraces[timeIndex + thread_x*devTracePts + i*traceCol];
					imageBlock[thread_x][thread_y] += weight * traceBlock;
				}

				//horizontal x/y reduction

				__syncthreads();

			#if THREADS_X > 511
				if (thread_x<256)
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+256][thread_y];

				__syncthreads();
			#endif

			#if THREADS_X > 255	
				if (thread_x<128)
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+128][thread_y];

				__syncthreads();

			#endif

			#if THREADS_X > 127
				if (thread_x<64)
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+64][thread_y];

				__syncthreads();

			#endif

			#if THREADS_X > 63
				if (thread_x<32){
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+32][thread_y];
			
			#elif THREADS_X==32
				if (thread_x<16){
			#endif
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+16][thread_y];
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+8][thread_y];
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+4][thread_y];
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+2][thread_y];
					imageBlock[thread_x][thread_y] += imageBlock[thread_x+1][thread_y];
				}

				__syncthreads();
		
				//write image point to global

				if ((thread_x==0) & imageBlock[0]>0){
		
					devImage[block_x * devImageDepths + offJ] += imageBlock[0][thread_y];
				}


			} flag=0;
		} //end for image time/depth
	} //end for trace x/y space 
} //end kernel



//extern "C"{
__host__ void ktmMigrationGPU(struct imageGrid* imageX, 
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
	
	// grab and copy scalars from structs
	//float temp;

	int IMAGE_X_SIZE = imageX->imagePnts;
	float IMAGE_X_STEP = imageX->imageStep;
	float IMAGE_X_START = imageX->imageStart;

	int IMAGE_Y_SIZE = imageY->imagePnts;
	float IMAGE_Y_STEP = imageY->imageStep;
	float IMAGE_Y_START = imageY->imageStart;
	

	int IMAGE_Z_SIZE = imageZ->imagePnts;
	float IMAGE_Z_STEP = imageZ->imageStep;
	float IMAGE_Z_START = imageZ->imageStart;


	 //clear image mem
	for(int i=0; i<IMAGE_X_SIZE*IMAGE_Y_SIZE*IMAGE_Z_SIZE; i++)
		image[i]=0.0f;
					        

	int TRACE_PTS 	= config->tracePts;
	int TRACE_NO	= config->traceNo;
	float SAMPLE_RATE = config->traceSampleRate;

	cudaMemcpyToSymbol("devTracePts",&TRACE_PTS,sizeof(int));
	cudaMemcpyToSymbol("devTraceNo",&TRACE_NO,sizeof(int));
	cudaMemcpyToSymbol("devTraceSampleRate",&SAMPLE_RATE,sizeof(float));
	cudaMemcpyToSymbol("devImageDepths",&IMAGE_Z_SIZE,sizeof(int));
	cudaMemcpyToSymbol("devImageXStart",&IMAGE_X_START,sizeof(float));
	cudaMemcpyToSymbol("devImageXStep",&IMAGE_X_STEP,sizeof(float));
	cudaMemcpyToSymbol("devImageYSize",&IMAGE_Y_SIZE,sizeof(int));
	cudaMemcpyToSymbol("devImageYStart",&IMAGE_Y_START,sizeof(float));
	cudaMemcpyToSymbol("devImageYStep",&IMAGE_Y_STEP,sizeof(float));
	//vectors

	float* tauS 		= (float*) malloc(IMAGE_Z_SIZE*sizeof(float));
	float* slownessS	= (float*) malloc(IMAGE_Z_SIZE*sizeof(float));
	
	for (int i=0; i<IMAGE_Z_SIZE; i++){
		tauS[i] 	= (IMAGE_Z_START + (float) i * IMAGE_Z_STEP)*(IMAGE_Z_START + (float) i * IMAGE_Z_STEP);
		slownessS[i]	= 1.0f / (slowness[i] * slowness[i]);
	}
	

	cudaMemcpyToSymbol(devSlownessS,slownessS,sizeof(float)*IMAGE_Z_SIZE);
	cudaMemcpyToSymbol(devTauS,tauS,sizeof(float)*IMAGE_Z_SIZE);

        float *devImage, *devTraces, *devMidX, *devMidY, *devOffX, *devOffY;
	int memSizeImage	= sizeof(float) * IMAGE_X_SIZE * IMAGE_Y_SIZE * IMAGE_Z_SIZE;
	int memSizeTraces	= sizeof(float) * TRACE_PTS * TRACE_NO;
	int memSizeMidX		= sizeof(float) * TRACE_NO;
	int memSizeMidY		= memSizeMidX;
	int memSizeOffX		= memSizeMidX;
	int memSizeOffY		= memSizeMidX;
		

        cudaMalloc((void**) &devImage,          	memSizeImage);
        cudaMalloc((void**) &devTraces,     	memSizeTraces);
        cudaMalloc((void**) &devMidX, 		memSizeMidX);
        cudaMalloc((void**) &devMidY,          	memSizeMidY);
	cudaMalloc((void**) &devOffX, 		memSizeOffX);
	cudaMalloc((void**) &devOffY, 		memSizeOffY);


	// Load device memory
        cudaMemcpy(devImage, 	image, 	memSizeImage,	cudaMemcpyHostToDevice);
       	cudaMemcpy(devTraces, 	traces, memSizeTraces,	cudaMemcpyHostToDevice);
       	cudaMemcpy(devMidX,	midX,   memSizeMidX,	cudaMemcpyHostToDevice);
        cudaMemcpy(devMidY, 	midY, 	memSizeMidY,	cudaMemcpyHostToDevice);
        cudaMemcpy(devOffX, 	offX, 	memSizeOffX,	cudaMemcpyHostToDevice);
        cudaMemcpy(devOffY, 	offY, 	memSizeOffY,	cudaMemcpyHostToDevice);


	if ((IMAGE_X_SIZE * IMAGE_Y_SIZE) > BLOCK_LIMIT){
		printf("ERROR: desired image x*y size exceeds block limit of %i\n",BLOCK_LIMIT);
		
		return;
	}


	dim3 threads,blocks;
	threads.x=THREADS_X; 				threads.y=THREADS_Y;
	blocks.x=IMAGE_X_SIZE*IMAGE_Y_SIZE; 		blocks.y=1;


	ktmCUDA<<<blocks, threads>>>(devMidX, devMidY, devOffX, devOffY, devTraces, devImage);

	printf("%s\n", cudaGetErrorString(cudaGetLastError()));


	//transfer back to host
	cudaMemcpy(image,      devImage,          memSizeImage,           cudaMemcpyDeviceToHost);


	//cleanup

	cudaFree(devImage);
	cudaFree(devTraces);
	cudaFree(devMidX);
	cudaFree(devMidY);
	cudaFree(devOffX);
	cudaFree(devOffY);
	cudaThreadExit();

	free(tauS);
	free(slownessS);

} //end host function 
//} //end extern
