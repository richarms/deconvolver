#include <iostream>
#include <fstream>
#include <cstdlib>
#include <random>
#include <fftw3.h>
#include <complex>
#include <math.h>
#include <string.h>
#include "types.h"

using namespace std;

int height = 100;
int width = height;
int sub_height = round(sqrt(height));
int sub_width = sub_height;
double fftw_scale = 1.0/sub_height;
double threshold = 0.000001;

vis dft(unsigned int col, unsigned int row){
  fl angle = -2*M_PI*((col%sub_width)*(row%sub_height)+(col/sub_width)*(row/sub_height))*fftw_scale;
  return std::polar(1., angle);
}

vis dft2(unsigned int col, unsigned int row){
  row += row/(sub_height/2+1)*(sub_height/2-1);
  fl angle = -2*M_PI*((col%sub_width)*(row%sub_height)+(col/sub_width)*(row/sub_height))*fftw_scale;
  return std::polar(1., angle);
}

/*vis dft2(unsigned int col, unsigned int row){
  row += row/(sub_height/2+1)*(sub_height/2-1);
  row = height - row;
  fl angle = 2*M_PI*((col%sub_width)*(row%sub_height)+(col/sub_width)*(row/sub_height)+(row%sub_height != 0 ? (col/sub_width) : 0))*fftw_scale;
  return std::polar(1., angle);
}*/

void fft2(vis * img, vis * uv_data){
  fftw_plan forward_fft = fftw_plan_dft_2d(sub_width, sub_height, reinterpret_cast<fftw_complex*>(img), reinterpret_cast<fftw_complex*>(uv_data), FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(forward_fft);
  fftw_destroy_plan(forward_fft);
}

int omega(vis a){
  double o = 2*M_PI*fftw_scale;
  int i;
  if (abs(a) < threshold)
    return -1;
  for (i=0; i<sub_height; ++i)
    if (abs(a - std::polar(1., o*i)) < threshold)
      break;
  return i;
}

void print_matrix(pix * a){
  for (int j=0; j<height; ++j){
    for (int i=0; i<width; ++i)
      printf("%f ", a[j*width + i]);
    printf("\n");
  }
}

void print_matrix(vis * a){
  for (int j=0; j<height; ++j){
    for (int i=0; i<width; ++i)
      printf("%s%i ", omega(a[j*width + i]) >=0 ? " ":"", omega(a[j*width + i]));
    printf("\n");
  }
}

void transpose(vis * a){
  for (int i=0; i<width; ++i)
    for (int j=0; j<i; ++j){
      vis tmp = a[j*width + i];
      a[j*width + i] = a[i*width + j];
      a[i*width + j] = tmp;
    }
}

int main(int argc, char ** argv) {
  using namespace std;
  
  // create identity matrix
  vis * eye = new vis[height*width];
  memset(eye, 0, height*width*sizeof(*eye));
  for (int i=0; i<height; ++i){
    eye[i*width + i] = vis(1,0);
  }

  vis * cmp = new vis[height*width];
  memset(cmp, 0, height*width*sizeof(*cmp));

  // fft
  for (int i=0; i<width; ++i) {
    fft2(&eye[i*height], &cmp[i*height]);
  }

  // transpose cmp
  transpose(cmp);

  // subtract dft components
  for (int i=0; i<width; ++i){
    for (int j=0; j<height; ++j){
      cmp[j*width + i] -= dft(i, j);
    }
  }

  /*for (int i=0; i<height*width; ++i)
    if (omega(cmp[i]) != 100)
      printf("ERROR");*/

  print_matrix(cmp);
}