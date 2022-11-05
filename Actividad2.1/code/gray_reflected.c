#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NUM_THREADS 20 
void HorizontalRotation();

void VerticalRotation();

void GrayScale();

int main()
{

  double start;
  double end;
  start = omp_get_wtime();

  omp_set_num_threads(NUM_THREADS);
  HorizontalRotation();
  VerticalRotation();
  GrayScale();


  end = omp_get_wtime();
  printf("Tardó %f seg\n", end - start);
  return 0;
}

void VerticalRotation()
{
  FILE *image, *outputImage, *lecturas;
  int row, col;
  image = fopen("toy.bmp", "rb");        // Imagen original a transformar
  outputImage = fopen("toy_vrotated.bmp", "wb"); // Imagen transformada
  long ancho;
  int counter = 0;
  long alto;
  unsigned char r, g, b, Pixel, Pixel1, Pixel2, Pixel3; // Pixel
  // unsigned char* ptr;
  unsigned char xx[54];

  for (int i = 0; i < 54; i++)
  {
    xx[i] = fgetc(image);
    fputc(xx[i], outputImage); // Copia cabecera a nueva imagen
  }

  ancho = (long)xx[20] * 65536 + (long)xx[19] * 256 + (long)xx[18];
  alto = (long)xx[24] * 65536 + (long)xx[23] * 256 + (long)xx[22];
  printf("largo img %li\n", alto);
  printf("ancho img %li\n", ancho);

  // Guarda la imagen en una matriz paa hacer la inversión
  unsigned char mat[ancho * 3][alto];


#pragma omp for schedule(static)
  for (row = 0; row < alto; row++)
  {
    for (col = 0; col < ancho * 3; col = col + 3)
    {
      b = fgetc(image);
      g = fgetc(image);
      r = fgetc(image);

      unsigned char pixel = 0.21 * r + 0.72 * g + 0.07 * b;

      mat[col][row] = pixel;
      mat[col + 1][row] = pixel;
      mat[col + 2][row] = pixel;
    }
  }
#pragma omp for schedule(static) 
  for (row = alto; row > 0; row--)
  {
    for (col = 0; col < ancho * 3; col++)
    {
      Pixel = mat[col][row];
      fputc(Pixel, outputImage);
      counter++;
    }
  }

  printf("%d\n", counter);

  fclose(image);
  fclose(outputImage);
}

void HorizontalRotation()
{
  FILE *image, *outputImage, *lecturas;
  int row, col;
  image = fopen("toy.bmp", "rb");        // Imagen original a transformar
  outputImage = fopen("toy_hrotated.bmp", "wb"); // Imagen transformada
  long ancho;
  int counter = 0;
  long alto;
  unsigned char r, g, b, Pixel, Pixel1, Pixel2, Pixel3; // Pixel
  // unsigned char* ptr;
  unsigned char xx[54];
  for (int i = 0; i < 54; i++)
  {
    xx[i] = fgetc(image);
    fputc(xx[i], outputImage); // Copia cabecera a nueva imagen
  }
  ancho = (long)xx[20] * 65536 + (long)xx[19] * 256 + (long)xx[18];
  alto = (long)xx[24] * 65536 + (long)xx[23] * 256 + (long)xx[22];
  printf("largo img %li\n", alto);
  printf("ancho img %li\n", ancho);

  //unsigned char threshold = 100;

  // Guardar imagen en una matriz
  unsigned char mat[ancho * 3][alto];

#pragma omp for schedule(static) 
  for (row = 0; row < alto; row++)
  {
    for (col = 0; col < ancho * 3; col = col + 3)
    {
      b = fgetc(image);
      g = fgetc(image);
      r = fgetc(image);

      unsigned char pixel = 0.21 * r + 0.72 * g + 0.07 * b;

      mat[col][row] = pixel;
      mat[col + 1][row] = pixel;
      mat[col + 2][row] = pixel;
    }
  }

#pragma omp for schedule(static) 
  for (row = 0; row < alto; row++)
  {
    for (col = ancho * 3; col > 0; col--)
    {
      Pixel = mat[col][row];
      fputc(Pixel, outputImage);
      counter++;
    }
  }

  printf("%d\n", counter);

  fclose(image);
  fclose(outputImage);
}

void GrayScale(){
    // Obtiene el número de threads
    omp_set_num_threads(NUM_THREADS);
    // Declare pointers for image & new image
    FILE *image, *outputImage, *lecturas;
  int row, col;
  image = fopen("toy.bmp", "rb");        // Imagen original a transformar
  outputImage = fopen("toy_gray.bmp", "wb"); // Imagen transformada
  long ancho;
  int counter = 0;
  long alto;
    int nthreads;

    unsigned char r, g, b;               //Pixel

    // Read first 54 header's line and move img memory counter
    unsigned char xx[54];
    for(int i=0; i<54; i++){
        xx[i] = fgetc(image);
        fputc(xx[i], outputImage);   //Copia cabecera a nueva imagen
    }

    // Calculate the width & height of the original image    
    ancho = (long)xx[20]*65536 + (long)xx[19]*256 + (long)xx[18];
    alto = (long)xx[24]*65536 + (long)xx[23]*256 + (long)xx[22];
    long n = ancho * alto * 3;
    // printf("n %li\n", n);

    // Measure initial time (startTime)
	const double startTime = omp_get_wtime();

    // Convert original img to gray scale
    nthreads = omp_get_num_threads();
    #pragma omp parallel for      
        for (int i = 0; i < n; i++){
        //while(!feof(image)){                                        //Grises
            b = fgetc(image);
            g = fgetc(image);
            r = fgetc(image);
            
            // Create new pixel in gray scale
            unsigned char pixel = 0.21*r + 0.72*g + 0.07*b;
            fputc(pixel, outputImage);  //b
            fputc(pixel, outputImage);  //g
            fputc(pixel, outputImage);  //r
        }

    // Measure final time (endTime)
	const double endTime = omp_get_wtime();

	// Print how much time did it take
	printf("Tiempo imagen gris = %f\n", endTime-startTime);

    // Close original image
    fclose(image);
    // Close new image
    fclose(outputImage);
}