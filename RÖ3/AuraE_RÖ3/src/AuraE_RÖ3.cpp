/*
  Smooth test program.
  Compile with c++ -O3 smooth.cpp -o smooth
*/

#include <cmath>
#include <cstdio>
#include <numeric>
#include <cassert>
#include <time.h>

/* Macros to access elements in arrays as a three-dimensional matrix */
#define GET_AT(x,y,z,ptr) *(ptr+(z)*inIncZ+(y)*inIncY+(x)*inIncX)
#define SET_AT_OUT(x,y,z,ptr,val) *(ptr+(z)*outIncZ+(y)*outIncY+(x)*outIncX)=val

/* Print out the matrices if debug is set to 1
   Use a small matrix for debugging, otherwise you gel LOTS of output
*/
#define debug 0

/* Define function to compute the square of x */
template <typename T>
static inline T mysqr(T x) {
  return x*x;
}

float* makeKernel(double radius,int*ksize) {
  int size = (int)radius*2+1;
  float *kernel = new float[size*size];
  double v;
  for (int y=0; y<size; y++) {
    for(int x=0; x<size; x++) {
      v = (double) std::exp(-0.5*((mysqr((x-radius)/((double)radius*2))+
				   mysqr((y-radius)/((double)radius*2)))/mysqr(0.2)));
      if(v<0.0005) v = 0.0;
      kernel[y*size+x]=v;
      /* If a kernel of all ones is used, the convolution has no effect.
	 Input and output matrices will be identical after the operation.
	 This kernel can be used for testing the implementation
      */
      /* kernel[y*size + x] = 1.0f;  */
    }
  }
  *ksize = size;
  return kernel;
}

void smooth(unsigned char* inPtr, unsigned char*outPtr, int, int ext[6], float*kernel,
	    double scale, int size, int inIncX, int inIncY, int inIncZ, int outIncX,
	    int outIncY, int outIncZ) {

  int uc,vc;
  uc = size / 2;
  vc = size / 2;

  /* Check that the kernel has an odd size */
  if ((size&1)!=1) {
    printf("\n\n\n******* Error, convolution kernel size not odd *******\n\n\n");
  }

  int xmin,xmax,ymin,ymax;
  int z = ext[4];
  xmin=ext[0];
  xmax=ext[1];
  ymin=ext[2];
  ymax=ext[3];

  double sum;
  int height = ymax, width = xmax;
  int i=0;
  unsigned char val;

  for(int y=ymin; y<=ymax; y++) {
    for(int x=xmin; x<=xmax; x++) {
      sum = 0.0;
      i = 0;

      for(int v = -vc; v<=vc; v++) {
	int ny=y+v;
	for(int u = -uc; u<=uc; u++) {
	  int nx=x+u;
	  if (nx<=0) nx = 0;
	  else if (nx>width) nx = width;
	  if (ny<=0) ny = 0;
	  else if (ny>height) ny = height;
	  val = GET_AT(nx,ny,z,inPtr);
	  sum += val*kernel[i++];
	}
      }
      SET_AT_OUT(x,y,z,outPtr,(unsigned char)(scale*sum));
    }
  }

}

int main(){
  /* Size of the images (use an even value) */
  const int width  = 2000;  /* Usa a small matrix if debug is set to 1 */
  const int height = 2000;
  const int iterations = 10;
  time_t start, stop;
  int i;

  unsigned char *input = new unsigned char[width*height];
  unsigned char *output = new unsigned char[width*height];

  // init pixels to some values
  //std::fill(input,input + width*height, 10);
  for (i=0; i<width*height; i+=2) {
    input[i] = (unsigned char)10;
    input[i+1] = (unsigned char)8;
  }

  /* Print the input matrix */
  if (debug) {
    printf("Input matrix:\n");
    int x, y;
    for (y=0; y<width; y++) {
      for(x=0; x<height; x++) {
	printf("%i ",input[y*width+x]);
      }
      printf("\n");
    }
  }

  int uext[6];
  uext[0] = 0;
  uext[1] = width-1;
  uext[2] = 0;
  uext[3] = height-1;
  uext[4] = 0;
  uext[5] = 0;

  double scale;
  int  size, inIncX,inIncY,inIncZ,outIncX,outIncY,outIncZ;
  int kernelRadius = 3;

  float *kernel = makeKernel(kernelRadius, &size);

  if (debug) {
    printf("Kernel:\n");
    int x, y;
    for (y=0; y<size; y++) {
      for(x=0; x<size; x++) {
	printf("%6.4f ",kernel[y*size+x]);
      }
      printf("\n");
    }
  }

  scale = 1.0 / std::accumulate(kernel,kernel +mysqr(size),0.0);

  printf("Smooth program\n");
  printf("Input matrix size is %d by %d\n", width, height);
  printf("Kernel size is %d\n", size);
  printf("Scale is %3.2f\n", scale);

  start = clock();  /* Start measuring time */

  // Run the smooth operation 10 times to get reliable timing
  for (int i=0; i<iterations; i++) {
    smooth(input, output, 0 ,uext, kernel, scale, size,
	   inIncX = 1, inIncY = height, inIncZ = 0, outIncX = 1,
	   outIncY = height, outIncZ = 0);
  }
  stop = clock();    /* Stop timer */
  printf("\nClock time for smooth operation %6.1f seconds\n\n",
	 ((double) (stop-start))/ CLOCKS_PER_SEC);

  /* Print the result */
  if (debug) {
    printf("Output matrix:\n");
    int x, y;
    for (y=0; y<width; y++) {
      for(x=0; x<height; x++) {
	printf("%i ",output[y*size+x]);
      }
      printf("\n");
    }
  }

  /*
    Assert that input and output images are same
    Can be used together with the kernel of all 1:s for testing
  */
  /* assert(std::equal(input,input+width*height,output)); */

  //std::ofstream out("pixels.dat");
  //out.write((const char*) output, sizeof(unsigned char) * width * height);
  //out.flush();

  delete[] kernel;
  delete[] input;
  delete[] output;

  return 0;
}
