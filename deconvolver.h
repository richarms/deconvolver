#ifndef DECONVOLVER_H_
#define DECONVOLVER_H_

#include "types.h"
#include "timer.h"
#include <fftw3.h>
#include <vector>
#include <string>

class deconvolver{
public:
  int width, uv_width, height;
  pix * img;

  deconvolver(std::string path);

  ~deconvolver(void);

  pix * deconvolve(int max_iterations = 0);

  double time();

  double fft_time();

  double matrixMult_time();

  double leastSquares_time();
private:
  vis * uv_data;
  std::vector<int> uv_mask;
  fl error;
  fl scale;
  fl fftw_scale;
  fl uv_scale;
  timer total_timer, fft_timer, matrix_timer, ls_timer;
  fftw_plan forward_fft, reverse_fft;

  fl l1(vis * visibilities);

  fl l2(vis * visibilities);

  fl l1(pix * img);

  fl l2(pix * img);

  int max(pix * img);

  int max(vis * img);

  void diff(vis * lhs, vis * rhs, vis * result);

  void mask(vis * uv_data);

  void fft2(pix * img, vis * uv_data);

  void fft2(vis * uv_data, pix * img);

  void fft2b(pix * img, vis * uv_data);

  void ifft2(vis * uv_data, pix * img);

  void dftt2(vis * uv_data, pix * result);

  vis dft(unsigned int col, unsigned int row);

  vis transposeDft(unsigned int col, unsigned int row);

  vis dftProduct(unsigned int col1, unsigned int row1, unsigned int col2, unsigned int row2);

  fl realDftProduct(unsigned int col1, unsigned int row1, unsigned int col2, unsigned int row2);

  void transpose(vis * uv_data);

  void fft_convolve(pix * img);

  void convolve_cycle(pix * img);
};

#endif