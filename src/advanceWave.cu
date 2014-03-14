//a reduction on a single MP with a little bit 'o ugly
//at the end to advance the wavefront
//WJB 03/11

#include <cuda.h>
#include <stdio.h>


//__constant__ int sizeTheta;
//__constant__ int sizePhi;
//


__global__ void advanceWave(float *devReductionKeysIn, 
		float *devReductionValuesIn,
		float *devReductionKeysOut,
		float *devReductionValuesOut,
		float *devRadius,
		float *devTTable,
		int sizeTheta,
		int sizePhi){	
	


	const int thread_x = blockIdx.x*blockDim.x + threadIdx.x;

	__shared__ float valBlock[LOC_RED_THREADS];
	__shared__ float keyBlock[LOC_RED_THREADS];

	//__shared__ int flag;

	valBlock[threadIdx.x] = devReductionValuesIn[thread_x];
	keyBlock[threadIdx.x] = devReductionKeysIn[thread_x];

	#if LOC_RED_THREADS > 511
		if (threadIdx.x<256){
		
			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+256]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+256];
				keyBlock[threadIdx.x] = thread_x+256;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}

		}

		__syncthreads();
	#endif

	#if LOC_RED_THREADS > 255	
		if (threadIdx.x<128){

			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+128]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+128];
				keyBlock[threadIdx.x] = thread_x+128;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}
		}

		__syncthreads();
	#endif

	#if LOC_RED_THREADS > 127
		if (threadIdx.x<64){

			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+64]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+64];
				keyBlock[threadIdx.x] = thread_x+64;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}
		}

		__syncthreads();
	#endif

	#if LOC_RED_THREADS > 63
		if (threadIdx.x<32){
			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+32]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+32];
				keyBlock[threadIdx.x] = thread_x+32;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}
		}
	#elif LOC_RED_THREADS==32
		if (threadIdx.x<16){
	#endif
			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+16]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+16];
				keyBlock[threadIdx.x] = thread_x+16;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}

			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+8]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+8];
				keyBlock[threadIdx.x] = thread_x+8;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}

			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+4]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+4];
				keyBlock[threadIdx.x] = thread_x+4;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}
			
			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+2]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+2];
				keyBlock[threadIdx.x] = thread_x+2;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}

			if (valBlock[threadIdx.x] > valBlock[threadIdx.x+1]){
				valBlock[threadIdx.x] = valBlock[threadIdx.x+1];
				keyBlock[threadIdx.x] = thread_x+1;			} else{

				keyBlock[threadIdx.x] = thread_x;
			}
			
			

	#if LOC_RED_THREADS==32
		}
	#endif

		__syncthreads();
		
		//write key/value pair to global

		if ((threadIdx.x==0)){
		
			devReductionKeysOut[blockIdx.x]=keyBlock[threadIdx.x];
			devReductionValuesOut[blockIdx.x]=valBlock[threadIdx.x];

		//if we're the lowest thread in the last block, increment wavefront radius at key position

			if ((devReductionKeysOut[0]==0.0f) || (devReductionKeysOut[1]==0.0f) || (devReductionKeysOut[2]==0.0f)
				|| (devReductionKeysOut[3]==0.0f)){
				return;
			} else {

				int key = sizeTheta*sizePhi -1;
				int value = devReductionValuesOut[0];
				
				if (devReductionValuesOut[1] < value ){
					key=devReductionKeysOut[1];
					value=devReductionValuesOut[1];
				}

				if (devReductionValuesOut[2] < value ){
					key=devReductionKeysOut[2];
					value=devReductionValuesOut[2];
				}

				if (devReductionValuesOut[3] < value ){
					key=devReductionKeysOut[3];
					value=devReductionValuesOut[3];
				}

				float oldRadius = devRadius[(key / sizeTheta) * sizePhi + key % sizePhi];
				devRadius[(key / sizeTheta) * sizePhi + key % sizePhi] = oldRadius + 1.0f;
				devTTable[(key / sizeTheta) * sizePhi + key % sizePhi + (int) oldRadius* sizePhi*sizeTheta]= value;

				devReductionKeysOut[0]=0.0f;
				devReductionKeysOut[1]=0.0f;
				devReductionKeysOut[2]=0.0f;
				devReductionKeysOut[3]=0.0f;
				
			}
				
		}



} //end kernel
