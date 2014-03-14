#include <cuda.h>
#include <stdio.h>



__global__ void eikSolverCUDA(float *devAliveTableInput,
		float *devAliveTableOutput,
                float *devRadius,
		float *devSlowness,
		int sizeR,
		float deltaR,
		float startR,
		int sizeTheta,
		float deltaTheta,
		float startTheta,
		int sizePhi,
		float deltaPhi,
		float startPhi){

        const int ix		= blockIdx.x*blockDim.x + threadIdx.x;
	const int iy		= blockIdx.y*blockDim.y + threadIdx.y;
	const int itile		= ix + iy*sizePhi;	
	const int iwest		= sizePhi + iy*sizePhi - 1;
	const int isouth	= ix;
	const int ieast		= iy*sizePhi;
	const int inorth	= ix + (sizeTheta-1)*sizePhi;


	//vel and grid data 
        float slowSq 	= devSlowness[itile]; 
	float radius	= devRadius[itile] * deltaR;
	float theta	= iy * deltaTheta;

	//traveltime data
        __shared__ float travelTime[EIK_SOLVER_THREADS_X][EIK_SOLVER_THREADS_Y];
        __shared__ float travelTimeWest[EIK_SOLVER_THREADS_Y];
        __shared__ float travelTimeSouth[EIK_SOLVER_THREADS_X];
        __shared__ float travelTimeEast[EIK_SOLVER_THREADS_Y];
        __shared__ float travelTimeNorth[EIK_SOLVER_THREADS_X];

        //load main tile
	travelTime[threadIdx.x][threadIdx.y] = devAliveTableInput[itile];

	//boundaries
	if (threadIdx.x==0) travelTimeWest[threadIdx.y] = devAliveTableInput[iwest];  
	if (threadIdx.y==blockDim.y-1) travelTimeSouth[threadIdx.x] = devAliveTableInput[isouth];  
	if (threadIdx.x==blockDim.x-1) travelTimeEast[threadIdx.y] = devAliveTableInput[ieast];  
	if (threadIdx.y==0) travelTimeNorth[threadIdx.x] = devAliveTableInput[inorth];  

        //coefficients
        float secondOrder 	= 0;
        float firstOrder 	= 0;
        float zerothOrder 	= 0;

	float delRSqInv 	= 1 / (deltaR * deltaR);
	float delTheRadInv	= 1 / (radius * deltaTheta);
	float delThePhiRadInv	= 1 / (radius * deltaPhi * sin(theta));

	//solve
	secondOrder 	+=	delRSqInv;
	secondOrder 	+=	delTheRadInv;
	secondOrder	+=	delThePhiRadInv;
	
	firstOrder 	+= (threadIdx.x < blockDim.x-1) ? - 2.0f * delRSqInv * travelTime[threadIdx.x+1][threadIdx.y] :
		-2.0f * delRSqInv * travelTimeEast[threadIdx.y];
 
	firstOrder 	+= (threadIdx.y < blockDim.y-1) ? - 2.0f * delTheRadInv * travelTime[threadIdx.x][threadIdx.y+1] :
		-2.0f * delTheRadInv * travelTimeSouth[threadIdx.x];

	firstOrder 	+= (threadIdx.y < blockDim.y-1) ? - 2.0f * delThePhiRadInv * travelTime[threadIdx.x][threadIdx.y+1] :
		-2.0f * delThePhiRadInv * travelTimeSouth[threadIdx.x];	

	zerothOrder	+= (threadIdx.x > 0) ? delRSqInv * travelTime[threadIdx.x-1][threadIdx.y] :
		delRSqInv * travelTimeWest[threadIdx.y];

	zerothOrder 	+= (threadIdx.y > 0) ? delTheRadInv * travelTime[threadIdx.x][threadIdx.y-1] : 
		delTheRadInv * travelTimeNorth[threadIdx.x];

	zerothOrder	+= (threadIdx.y > 0) ? delThePhiRadInv * travelTime[threadIdx.x][threadIdx.y-1] :
		delThePhiRadInv * travelTimeNorth[threadIdx.x];

	zerothOrder	+= - slowSq;

	float ttime 	= - firstOrder + sqrt(firstOrder * firstOrder - 4 * secondOrder * zerothOrder) / ( 2.0f * secondOrder);

	//write out to flat array
	devAliveTableOutput[itile]=ttime;


}
