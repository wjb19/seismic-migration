#include <stdio.h>
// to convert IBM floats to IEEE

void ktmIBM32toIEEE(unsigned int* input, float* output, int n){


	//ibm32 format; 31bit = sign, 30:24bits = exponent, 23:0bits = fraction
	//exponent is base 16, bias 64
	//bit31 is leftmost

	//temp vars

	
	float 		 sign	=0.0f;
	unsigned int fraction 	=0;
	unsigned int exponent 	=0;
 
	for (int i=0; i<n; i++){

		//4 byte; sign and exponent			
			
		sign = (-2.0f * (float) ((input[i] & SIGN_MASK) > 0)) + 1.0f;

		exponent = (input[i] & EXPONENT_MASK);
		

		//3,2,1 byte; fraction

		fraction = (input[i] & FRACTION_MASK);
		// bit shift 24 places & remove bias in exp 
		float exp  = (float) ((exponent >> 24) & 127u) - 64.0f;

		output[i]=sign * ((float) fraction / (float) FRACTION_MASK) * powf(16.0f,exp);

		sign 		=0.0f;
		exponent	=0;
		fraction	=0;
	}
	return;
}

// file i/o functions

int fileSizeFourBytes(FILE *fin){

	int z=-1; 
 	float x;

  	while(!feof(fin)){
    		fread(&x, sizeof(unsigned int), 1, fin); 
      		z++;
	}

	fseek(fin, 0, SEEK_SET);
  	return z;						    
}



unsigned int* binReadUints(FILE *fin, int n){
	

	unsigned int * x_ptr = (unsigned int*) malloc(n * sizeof(unsigned int));
 	fread(x_ptr, sizeof(unsigned int), n, fin);
	fseek(fin, 0, SEEK_SET);

 	return x_ptr;
}

float * binReadFloats(FILE *fin, int n){

	float * x_ptr = (float*) malloc(n*sizeof(float));
	fread(x_ptr, sizeof(float), n, fin);
	fseek(fin, 0, SEEK_SET);

	return x_ptr;

}

void fileWrite(FILE *fout, float *ptr, int t){
  
	fwrite(ptr, sizeof(float), t, fout);
  	return;
}



void errorMessage(char *str){
    printf("Can't open %s.\n",str);
    return;
}


void loadInputData(struct fileInputPtrs* input){

	input->flag=0;

	FILE *fptr1;  
  	char filename1[]="../data/inputData.bin";

  	if ((fptr1 = fopen(filename1, "r+b")) == NULL){
    		errorMessage(filename1);
		input->flag=1;	
  	} else{
		input->inputData = fptr1;
	}


 	FILE *fptr2;
  	char filename2[]="../data/inputSrcX.bin";

  	if ((fptr2 = fopen(filename2, "r+b")) == NULL){
    		errorMessage(filename2);
 		input->flag=1; 
  	} else{
		input->inputSrcX= fptr2;
	}


  	FILE *fptr3;
  	char filename3[]="../data/inputSrcY.bin";

  	if ((fptr3 = fopen(filename3, "r+b")) == NULL){
    		errorMessage(filename3);
  		input->flag=1;

	} else {
		input->inputSrcY=fptr3;
	}

 	FILE *fptr4;
  	char filename4[]="../data/inputRecX.bin";

  	if ((fptr4 = fopen(filename4, "r+b")) == NULL){
    		errorMessage(filename4);
  		input->flag=1;
	} else {
		input->inputRecX=fptr4;
	}

 	FILE *fptr5;
  	char filename5[]="../data/inputRecY.bin";

  	if ((fptr5 = fopen(filename5, "r+b")) == NULL){
    		errorMessage(filename5);
		input->flag=1;  
	} else {
		input->inputRecY=fptr5;
	}



}

void loadRTM2dConfig(){

	//read sdtdin and setup a struct



}

void loadKTMConfig(){
	//ditto


}


void loadKDMConfig(){

	//ditto


}
