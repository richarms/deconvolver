#include <cstdio>
#include <string>
#include "imager/image.h"
#include "deconvolver.h"

using namespace std;

int main(int argc, char ** argv) {
  if (argc < 3)
    printf("Needs an input file and an output file\n");
  string in_path(argv[1]);
  string out_path(argv[2]);
  printf("Cleaning %s\n", in_path.c_str());
  deconvolver dc = deconvolver(in_path);
  dc.deconvolve(8);
  printf("Total time: %fs, FFT time: %fs, MatrixMult time: %f, LeastSquares time: %f\n", dc.time(), dc.fft_time(), dc.matrixMult_time(), dc.leastSquares_time());
  image out(dc.height, dc.width, dc.img);
  if (argc >= 3) {
    string sky_path(argv[3]);
    image sky(sky_path);
    printf("L2 Error of image with %s: %f\n", sky_path.c_str(), out.l2Difference(sky));
  }
  out.save_pgm(out_path+".pgm");
}