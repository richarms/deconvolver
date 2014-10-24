#include <iostream>
#include <fstream>
#include <cstdlib>
#include <random>
#include <fftw3.h>
#include "image.h"

using namespace std;

int main(int argc, char ** argv) {
  using namespace std;
  if (argc < 4){
    printf("Needs args\n");
    exit(2);
  }
  string mode(argv[1]);
  if (mode == "create"){
    if (argc < 5)
      printf("create needs 3 args (dest, width, height)\n");
    string dest(argv[2]);
    int width = atoi(argv[3]);
    int height = atoi(argv[4]);
    printf("Creating blank %ix%i image \n", width, height);
    image img(height, width);
    img.save(dest);
    img.save_pgm(dest+".pgm");
  } 

  else if (mode == "add"){
    if (argc < 7)
      printf("add needs 5 args (source, dest, x, y, intensity)\n");
    string source(argv[2]);
    string dest(argv[3]);
    int x = atoi(argv[4]);
    int y = atoi(argv[5]);
    float intensity = atof(argv[6]);
    printf("Adding source at (%i, %i) with intensity %f Jy\n", x, y, intensity);
    image img(source);
    img.setpix(x, y, vis(intensity, 0));
    img.save(dest);
    img.save_pgm(dest+".pgm");
  }

  else if (mode == "add_circ"){
    if (argc < 8)
      printf("add needs 5 args (source, dest, x, y, r, intensity)\n");
    string source(argv[2]);
    string dest(argv[3]);
    int x = atoi(argv[4]);
    int y = atoi(argv[5]);
    int r = atoi(argv[6]);
    float intensity = atof(argv[7]);
    printf("Adding circular extended source at (%i, %i) with radius %i and intensity %f Jy\n", x, y, r, intensity);
    image img(source);
    for (int i=0; i<img.getWidth(); ++i)
      for (int j=0; j<img.getHeight(); ++j)
        if ((i-x)*(i-x) + (j-y)*(j-y) < r*r)
          img.setpix(i, j, vis(intensity, 0));
    img.save(dest);
    img.save_pgm(dest+".pgm");
  }

  else if (mode == "transform"){
    if (argc < 4)
      printf("add needs 2 args (source, dest)\n");
    string source(argv[2]);
    string dest(argv[3]);
    printf("Transforming image to Fourier plane\n");
    image img(source);  
    pix * img_data = new pix[img.getHeight()*img.getWidth()];
    for (int i=0; i<img.getHeight()*img.getWidth(); ++i)
      img_data[i] = real(img[i]);
    image dest_img(img.getHeight(), img.getWidth()/2+1);
    fftw_plan forward_fft = fftw_plan_dft_r2c_2d(img.getWidth(), img.getHeight(), (img_data), reinterpret_cast<fftw_complex*>(dest_img.data()), FFTW_MEASURE);
    fftw_execute(forward_fft);
    fftw_destroy_plan(forward_fft);
    for (int i=0; i<img.getHeight()*(img.getWidth()/2+1); ++i)
      dest_img[i] = dest_img[i]/double(img.getHeight());
    dest_img.save(dest);
    dest_img.save_pgm(dest+".pgm");
  }

  else if (mode == "noise"){
    if (argc < 5)
      printf("add needs 3 args (source, dest, SNR)\n");
    string source(argv[2]);
    string dest(argv[3]);
    float sigma = atof(argv[4]);
    printf("Adding noise with SNR of %f", sigma);
    image img(source);
    if (sigma != 0){
      //determine variance
      vis mean = 0;
      for (int i=0; i<img.getHeight()*(img.getWidth()/2+1); ++i)
        mean += img[i];
      mean /= img.getHeight()*(img.getWidth()/2+1);
      vis var = 0;
      for (int i=0; i<img.getHeight()*(img.getWidth()/2+1); ++i)
        var += vis(pow(real(mean - img[i]), 2), pow(imag(mean - img[i]), 2));
      var /= img.getHeight()*(img.getWidth()/2+1);
      float stddev_r = sqrt(real(var)/sigma);
      float stddev_i = sqrt(imag(var)/sigma);
      printf(" (OR %f, OI %f, VR %f, VI %f)", stddev_r, stddev_i, real(var), imag(var));
      //add noise
      default_random_engine de(time(0));
      normal_distribution<double> ndr(0,stddev_r);
      normal_distribution<double> ndi(0,stddev_i);
      for (int i=0; i<img.getHeight()*img.getWidth(); ++i)
        img[i] = img[i] + vis(ndr(de), ndi(de));
    }
    printf("\n");
    img.save(dest);
    img.save_pgm(dest+".pgm");
  }

  else if (mode == "mask"){
    if (argc < 5)
      printf("add needs 3 args (source, dest, mask)\n");
    string source(argv[2]);
    string dest(argv[3]);
    string mask(argv[4]);
    printf("Applying uv mask\n");
    image img(source);
    image uv(mask);
    image target(img.getHeight(), img.getWidth());
    for (int i=0; i<img.getWidth(); ++i)
      for (int j=0; j<img.getHeight(); ++j)
        if (abs(uv[j*img.getWidth() + i]) != 0)
          target[j*(img.getWidth()) + i] = img[j*(img.getWidth()) + i];
    target.save(dest);
    target.save_pgm(dest+".pgm");
  }
  
  else if (mode == "dirty"){
    if (argc < 4)
      printf("add needs 2 args (source, dest)\n");
    string source(argv[2]);
    string dest(argv[3]);
    printf("Constructing Dirty Image\n");
    image img(source);  
    pix * img_data = new pix[img.getHeight()*(img.getWidth()-1)*2];
    image dest_img(img.getHeight(), (img.getWidth()-1)*2);
    fftw_plan forward_fft = fftw_plan_dft_c2r_2d(dest_img.getWidth(), dest_img.getHeight(), reinterpret_cast<fftw_complex*>(img.data()), (img_data), FFTW_MEASURE);
    fftw_execute(forward_fft);
    fftw_destroy_plan(forward_fft);
    for (int i=0; i<dest_img.getHeight()*dest_img.getWidth(); ++i)
      dest_img[i] = img_data[i]/double(img.getHeight());
    dest_img.save(dest);
    dest_img.save_pgm(dest+".pgm");
  }

  else if (mode == "diff"){
    if (argc < 4)
      printf("add needs 2 args (lhs, rhs)\n");
    string lhs_path(argv[2]);
    string rhs_path(argv[3]);
    image lhs(lhs_path);
    image rhs(rhs_path);
    for (int i=0; i<lhs.getHeight(); ++i){
      for (int j=0; j<lhs.getWidth(); ++j)
        printf("(%f, %f) ", abs(lhs[i*lhs.getWidth() + j]), abs(rhs[i*lhs.getWidth() + j]));
      printf("\n");
    }
  }

  else {
    printf("Unknown option: %s\n", mode.c_str());
  }
}