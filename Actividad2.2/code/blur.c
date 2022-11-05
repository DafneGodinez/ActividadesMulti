#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define NUM_THREADS 10
 
void blurring(char outputImage[], int sizeM);
void blurring_inv(char outputImageName[], int sizeM);
 
int main(){
 #pragma omp parallel
   {
       #pragma omp sections
       {
       #pragma omp section
           blurring("blurred1_sample.bmp", 3);
       #pragma omp section
           blurring("blurred2_sample.bmp", 5);
       #pragma omp section
           blurring("blurred3_sample.bmp", 7);
       #pragma omp section
           blurring("blurred4_sample.bmp", 9);
       #pragma omp section
           blurring("blurred5_sample.bmp", 11);
       #pragma omp section
           blurring("blurred6_sample.bmp", 13);
       #pragma omp section
           blurring("blurred7_sample.bmp", 15);
       #pragma omp section
           blurring("blurred8_sample.bmp", 17);
       #pragma omp section
           blurring("blurred9_sample.bmp", 19);
       #pragma omp section
           blurring("blurredI0_sample.bmp", 21);
       #pragma omp section
           blurring_inv("blurred_inv1_sample.bmp", 21);  
       #pragma omp section
           blurring_inv("blurred_inv2_sample.bmp", 19);
       #pragma omp section
           blurring_inv("blurred_inv3_sample.bmp", 17);
       #pragma omp section
           blurring_inv("blurred_inv4_sample.bmp", 15);
       #pragma omp section
           blurring_inv("blurred_inv5_sample.bmp", 13);
       #pragma omp section
           blurring_inv("blurred_inv6_sample.bmp", 11);
       #pragma omp section
           blurring_inv("blurred_inv7_sample.bmp", 9);
       #pragma omp section
           blurring_inv("blurred_inv8_sample.bmp", 7); 
       #pragma omp section
           blurring_inv("blurred_inv9_sample.bmp", 5);
       #pragma omp section
           blurring_inv("blurred_inv10_sample.bmp", 3);    
       }
   }
 return 0;
}
 
void blurring(char outputImageName[], int sizeM){
 FILE *image, *outputImage, *lecturas;
 image = fopen("sample.bmp","rb");                 //Imagen original a transformar
 outputImage = fopen(outputImageName,"wb");      //Imagen transformada
 long ancho;
 long alto;
 unsigned char r, g, b;    //Pixel
 unsigned char* ptr;
 
 unsigned char xx[54];
 long cuenta = 0, anchoR = 0, altoR = 0, anchoM = 0, altoM = 0;
 long sum;
 int iR, jR;
 for(int i = 0; i < 54; i++) {
   xx[i] = fgetc(image);
   fputc(xx[i], outputImage);      //Copia cabecera a nueva imagen
 }
 ancho = (long)xx[20]*65536 + (long)xx[19]*256 + (long)xx[18];
 alto = (long)xx[24]*65536 + (long)xx[23]*256 + (long)xx[22];
 printf("largo img %li\n",alto);
 printf("ancho img %li\n",ancho);
 
 ptr = (unsigned char*)malloc(alto*ancho*3* sizeof(unsigned char));
 omp_set_num_threads(NUM_THREADS);
  
 unsigned char foto[alto][ancho], fotoB[alto][ancho];
 unsigned char pixel;
    for(int i = 0; i < alto; i++){
     for(int j = 0; j < ancho; j++){
       b = fgetc(image);
       g = fgetc(image);
       r = fgetc(image);
 
       pixel = 0.21*r + 0.72*g + 0.07*b;
       foto[i][j] = pixel;
       fotoB[i][j] = pixel;
     }
   }
 
   anchoR = ancho/sizeM;
   altoR = alto/sizeM;
   anchoM = ancho%sizeM;
   altoM = alto%sizeM;
 
   int inicioX,inicioY,cierreX,cierreY,ladoX,ladoY;
   //casi toda la imagen se representa aquí
   for(int i = 0; i < alto; i++){
     //iR=i*sizeM;
     for(int j = 0; j < ancho; j++){
       //jR=j*sizeM;
       if(i < sizeM){
         inicioX = 0;
         cierreX = i+sizeM;
         ladoX = i+sizeM;
       } else if(i >= alto-sizeM){
         inicioX = i-sizeM;
         cierreX = alto;
         ladoX = alto-i+sizeM;
       }else{
         inicioX = i-sizeM;
         cierreX = i+sizeM;
         ladoX = sizeM*2+1;
       }
  
       if(j < sizeM){
         inicioY = 0;
         cierreY = j+sizeM;
         ladoY = j+sizeM;
       }else if(j >= ancho-sizeM){
         inicioY = j-sizeM;
         cierreY = ancho;
         ladoY = ancho-j+sizeM;
       }else{
         inicioY = j-sizeM;
         cierreY = j+sizeM;
         ladoY = sizeM*2+1;
       }
       sum = 0;
       for(int x = inicioX; x < cierreX; x++){
         for(int y = inicioY; y < cierreY; y++){
           sum += foto[x][y];
         }
       }
       sum = sum/(ladoX*ladoY);
       fotoB[i][j] = sum;
     }
   }
 
 //asignacion a imagen
 cuenta = 0;
 for(int i = 0; i < alto; i++){
   for(int j = 0; j < ancho; j++){
     ptr[cuenta] = fotoB[i][j]; //b
     ptr[cuenta+1] = fotoB[i][j]; //g
     ptr[cuenta+2] = fotoB[i][j]; //r
     cuenta++;
   }
 }      
 
 //Grises
 const double startTime = omp_get_wtime();
 #pragma omp parallel
 {
   #pragma omp for schedule(dynamic)
   for (int i = 0; i < alto*ancho; ++i) {
     fputc(ptr[i], outputImage);
     fputc(ptr[i+1], outputImage);
     fputc(ptr[i+2], outputImage);
   }
 }
 
 const double endTime = omp_get_wtime();
 free(ptr);
 fclose(image);
 fclose(outputImage);
 printf("Tiempo= %f\n", endTime-startTime);
}
 
void blurring_inv(char outputImageName[], int sizeM){
 FILE *image, *outputImage, *lecturas;
 image = fopen("sample.bmp","rb");                 //Imagen original a transformar
 outputImage = fopen(outputImageName,"wb");      //Imagen transformada
 long ancho;
 long alto;
 unsigned char r, g, b;    //Pixel
 unsigned char* ptr;
 
 unsigned char xx[54];
 long cuenta = 0, anchoR = 0, altoR = 0, anchoM = 0, altoM = 0;
 long sum;
 int iR, jR;
 for(int i = 0; i < 54; i++) {
   xx[i] = fgetc(image);
   fputc(xx[i], outputImage);      //Copia cabecera a nueva imagen
 }
 ancho = (long)xx[20]*65536 + (long)xx[19]*256 + (long)xx[18];
 alto = (long)xx[24]*65536 + (long)xx[23]*256 + (long)xx[22];
 printf("largo img %li\n",alto);
 printf("ancho img %li\n",ancho);
 
 ptr = (unsigned char*)malloc(alto*ancho*3* sizeof(unsigned char));
 omp_set_num_threads(NUM_THREADS);
  
 unsigned char foto[alto][ancho], fotoB[alto][ancho];
 unsigned char pixel;
    for(int i = 0; i < alto; i++){
     for(int j = 0; j < ancho; j++){
       b = fgetc(image);
       g = fgetc(image);
       r = fgetc(image);
 
       pixel = 0.21*r + 0.72*g + 0.07*b;
       foto[i][j] = pixel;
       fotoB[i][j] = pixel;
     }
   }
 
   anchoR = ancho/sizeM;
   altoR = alto/sizeM;
   anchoM = ancho%sizeM;
   altoM = alto%sizeM;
 
   int inicioX,inicioY,cierreX,cierreY,ladoX,ladoY;
   //casi toda la imagen se representa aquí
   for(int i = alto; i > 0; i--){
     //iR=i*sizeM;
     for(int j = ancho; j > 0; j--){
       //jR=j*sizeM;
       if(i < sizeM){
         inicioX = alto;
         cierreX = 0;
         ladoX = i+sizeM;
       } else if(i >= alto-sizeM){
         inicioX = i+sizeM;
         cierreX = 0;
         ladoX = alto-i+sizeM;
       }else{
         inicioX = i-sizeM;
         cierreX = i+sizeM;
         ladoX = sizeM*2+1;
       }
  
       if(j < sizeM){
         inicioY = ancho;
         cierreY = j+sizeM;
         ladoY = j+sizeM;
       }else if(j >= ancho-sizeM){
         inicioY = j-sizeM;
         cierreY = ancho;
         ladoY = ancho-j+sizeM;
       }else{
         inicioY = j-sizeM;
         cierreY = j+sizeM;
         ladoY = sizeM*2+1;
       }
       sum = 0;
       for(int x = inicioX; x < cierreX; x++){
         for(int y = inicioY; y < cierreY; y++){
           sum += foto[x][y];
         }
       }
       sum = sum/(ladoX*ladoY);
       fotoB[i][j] = sum;
     }
   }
 
 //asignacion a imagen
 cuenta = 0;
 
 /*Matriz inversa(Matriz m) {
 
 m->ptr[cuenta] = m->fotoB[i][j];
 m->ptr[cuenta+1] = m->fotoB[i][j];
 m->ptr[cuenta+2] = m->fotoB[i][j];
 cuenta ++
 
 
 return m;
 */
 
 for(int i = 0; i<alto ; i++){
   for(int j = 0; j<ancho; j++){
     ptr[cuenta] = fotoB[i][ancho-j]; //b
     ptr[cuenta+1] = fotoB[i][ancho-j]; //g
     ptr[cuenta+2] = fotoB[i][ancho-j]; //r
     cuenta++;
   }
 }      
              
 //Grises
 const double startTime = omp_get_wtime();
 #pragma omp parallel
 {
   #pragma omp for schedule(dynamic)
   for (int i = 0; i < alto*ancho; ++i) {
     fputc(ptr[i], outputImage);
     fputc(ptr[i+1], outputImage);
     fputc(ptr[i+2], outputImage);
   }
 }
 
 const double endTime = omp_get_wtime();
 free(ptr);
 fclose(image);
 fclose(outputImage);
 printf("Tiempo para imagen invertida = %f\n", endTime-startTime);
}
