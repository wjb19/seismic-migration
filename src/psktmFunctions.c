void initTwoAxes(struct imageGrid* imageY,
        struct imageGrid* imageZ,
        float* imageYY,
        float* imageZZ){



        for(int i=0; i< imageY->imagePnts; i++)
                imageYY[i] = imageY->imageStart + (float) i * imageY->imageStep;

        for(int i=0; i< imageZ->imagePnts; i++)
                imageZZ[i] = imageZ->imageStart + (float) i * imageZ->imageStep;


        return;

}

void initThreeAxes(struct imageGrid* imageX,
        struct imageGrid* imageY,
        struct imageGrid* imageZ,
        float* imageXX,
        float* imageYY,
        float* imageZZ){


        for (int i=0; i< imageX->imagePnts; i++)
                imageXX[i] = imageX->imageStart + (float) i * imageX->imageStep;

        for(int i=0; i< imageY->imagePnts; i++)
                imageYY[i] = imageY->imageStart + (float) i * imageY->imageStep;

        for(int i=0; i< imageZ->imagePnts; i++)
                imageZZ[i] = imageZ->imageStart + (float) i * imageZ->imageStep;


        return;

}

