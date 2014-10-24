#include "deconvolver.h"
#include <math.h>
#include <string.h>
#include <typeinfo>
#include "imager/image.h"

extern "C" void dgesv_(int * N, int * NRHS, double * A, int * LDA, int * IPIV, double * B, int * LBD, int * INFO);
extern "C" void dposv_(char * uplo, int * N, int * NRHS, double * A, int * LDA, double * B, int * LBD, int * INFO);
extern "C" void dsgesv_(int * N, int * NRHS, double * A, int * LDA, int * IPIV, double * B, int * LBD, double * X, int * LDX, double * WORK, float * SWORK, int * ITER, int * INFO);

deconvolver::deconvolver(std::string path){
  image in(path);
  uv_width = in.getWidth();
  width = uv_width*2-2;
  height = in.getHeight();
  uv_data = (vis *) fftw_malloc(height*uv_width*sizeof(vis));
  img = (pix *) fftw_malloc(height*width*sizeof(pix));
  forward_fft = fftw_plan_dft_r2c_2d(width, height, img, reinterpret_cast<fftw_complex*>(uv_data), FFTW_ESTIMATE);
  reverse_fft = fftw_plan_dft_c2r_2d(width, height, reinterpret_cast<fftw_complex*>(uv_data), img, FFTW_ESTIMATE);
  memcpy(uv_data, in.data(), height*uv_width*sizeof(*in.data()));
  uv_mask = std::vector<int>();
  for (int i=0; i<height*uv_width; ++i)
    if (abs(uv_data[i]) != 0)
      uv_mask.push_back(i);
  memset(img, 0, height*width*sizeof(*img));
  error = 0.0;
  scale = fl(1)/fl(height * width);
  fftw_scale = fl(1)/sqrt(height*width);
  uv_scale = fl(1)/sqrt(uv_mask.size());
}

deconvolver::~deconvolver(void){
  fftw_free(uv_data);
  fftw_free(img);
}

pix * deconvolver::deconvolve(int max_iterations){
  //start timer
  total_timer.start();
  std::vector<int> atoms;
  std::vector<fl> weights;
  // initialize image
  memset(img, 0, height*width*sizeof(*img));
  // set initial residue equal to uv_data
  vis * residue = (vis *) fftw_malloc(height*uv_width*sizeof(vis));
  memcpy(residue, uv_data, height*uv_width*sizeof(*uv_data));
  // Get dirty image 2d fft of uv_data
  vis * uv_temp = (vis *) fftw_malloc(height*uv_width*sizeof(vis));
  memcpy(uv_temp, uv_data, height*uv_width*sizeof(*uv_data));
  pix * dirty = (pix *) fftw_malloc(height*width*sizeof(pix));
  ifft2(uv_temp, dirty);
  fftw_free(uv_temp);

  // Get the inverse fft of the uv mask (for the convolution theorem)
  vis * full_mask = (vis *) fftw_malloc(height*uv_width*sizeof(vis));
  memset(full_mask, 0, height*uv_width*sizeof(*full_mask));
  for (int i : uv_mask)
    full_mask[i] = 1;
  pix * fft_mask = (pix *) fftw_malloc(height*width*sizeof(pix));
  ifft2(full_mask, fft_mask);
  for (int i=0; i<height*width; ++i)
    fft_mask[i] *= scale;
  fft_convolve(fft_mask);
  //convolve_cycle(fft_mask_cpx);
  fftw_free(full_mask);
  
  // set initial projection to the dirty image
  pix * projection = (pix *) fftw_malloc(height*width*sizeof(pix));
  memcpy(projection, dirty, height*width*sizeof(*projection));
  for (int i=0; i<width*height; ++i)
    dirty[i] *= fftw_scale;

  // get initatial iteration time
  double time = total_timer.time();

  printf("Completed initatialisation (%fs), residual error (L2): %f\n", total_timer.time(), l2(residue));

  // begin main loop
  for (int i=0; !((l2(residue) < uv_mask.size()*error) || (max_iterations != 0 && i >= max_iterations)); ++i){
    // --- project residue into image space ---
    ifft2(residue, projection);
    //fft_convolve(projection);

    // --- select best unmatched atom ---
    // deselect any selected rows
    for (unsigned int row=0; row < atoms.size(); ++row)
      projection[atoms[row]] = 0;
    // select max
    atoms.push_back(max(projection));
    weights.push_back(real(projection[atoms.back()])*fftw_scale);

    // --- renew weights ---
    // find matrix for LHS of the normal equation (A^T A)
    matrix_timer.start();
    fl * left_norm = new fl[atoms.size()*atoms.size()];
    memset(left_norm, 0, atoms.size()*atoms.size()*sizeof(*left_norm));
    for (unsigned int row=0; row < atoms.size(); ++row){
      for (unsigned int col=0; col < atoms.size(); ++col){
        left_norm[row * atoms.size() + col] = fft_mask[((height - atoms[row]/width + atoms[col]/width)%height)*width + (width - atoms[row]%width + atoms[col]%width)%width];
        //printf("%f ", left_norm[row * atoms.size() + col]);
      }
      //printf("\n");
    }
    matrix_timer.stop();

    for (unsigned int row=0; row < atoms.size(); ++row)
      weights[row] = real(dirty[atoms[row]]);

    
    int dim = atoms.size();
    int nrhs = 1;
    int info;
    int * ipiv = new int[dim];
    ls_timer.start();
    //for general matrix
    dgesv_(&dim, &nrhs, left_norm, &dim, ipiv, &weights[0], &dim, &info);

    // for positive definite matrix
    //char uplo = 'U';
    //dposv_(&uplo, &dim, &nrhs, left_norm, &dim, &weights[0], &dim, &info);

    // for iterative solving using single precision arithmetic
    /*double * x = new double[dim*nrhs];
    double * work = new double[dim*nrhs];
    float * swork = new float[dim*(dim+nrhs)];
    int iter;
    dsgesv_(&dim, &nrhs, left_norm, &dim, ipiv, &weights[0], &dim, x, &dim, work, swork, &iter, &info);*/
    ls_timer.stop();

    fftw_free(ipiv);
    fftw_free(left_norm);

    // --- update residue ---
    // construct image
    memset(img, 0, height*width*sizeof(*img));
    for (unsigned int j=0; j<atoms.size(); ++j){
//printf("%i: %f\t\t", atoms[j], weights[j]);
      img[atoms[j]] = weights[j];
    }
//printf("\n");

    // calculate residue
    fft2(img, residue);
    mask(residue);
    diff(uv_data, residue, residue);
    printf("Completed iteration %i (%fs), residual error (L2): %f\n", i, total_timer.time()-time, l2(residue));
    time = total_timer.time();
  }
  fftw_free(projection);
  fftw_free(dirty);
  fftw_free(residue);
  //stop timer
  total_timer.stop();
  return img;
} 

fl deconvolver::l1(vis * visibilities){
  float sum = 0;
  for (unsigned int i=0; i<uv_mask.size(); ++i)
    sum += std::abs(visibilities[uv_mask[i]]);
  return sum;
}

fl deconvolver::l2(vis * visibilities){
  float sum = 0;
  for (unsigned int i=0; i<uv_mask.size(); ++i)
    sum += pow(abs(visibilities[uv_mask[i]]), 2);
  return sum;
}

fl deconvolver::l1(pix * img){
  float sum = 0;
  for (int i=0; i<width*height; ++i)
    sum += abs(img[i]);
  return sum;
}

fl deconvolver::l2(pix * img){
  float sum = 0;
  for (int i=0; i<width*height; ++i)
    sum += pow(abs(img[i]), 2);
  return sum;
}

int deconvolver::max(pix * img){
  int ix = 0;
  float max = abs(img[0]);
  for (int i=1; i<width*height; ++i)
    if (abs(img[i]) > max){
      ix = i;
      max = abs(img[i]);
    }
  return ix;
}

int deconvolver::max(vis * img){
  int ix = 0;
  float max = abs(real(img[0]));
  for (int i=1; i<uv_width*height; ++i)
    if (abs(real(img[i])) > max){
      ix = i;
      max = abs(real(img[i]));
    }
  return ix;
}

void deconvolver::diff(vis * lhs, vis * rhs, vis * result){
  for (int i=0; i<uv_width*height; ++i)
    result[i] = lhs[i] - rhs[i]*fftw_scale;
}

void deconvolver::mask(vis * uv_data){
  auto it = uv_mask.begin();
  for (int i=0; i<height*uv_width; ++i){
    if (i - *it == 0)
      ++it;
    else
      uv_data[i] = 0;
  }
}

void deconvolver::fft2(pix * img, vis * uv_data){
  fft_timer.start();
  fftw_execute_dft_r2c(forward_fft, img, reinterpret_cast<fftw_complex*>(uv_data));
  fft_timer.stop();
}

void deconvolver::fft2(vis * uv_data, pix * img){
  fft_timer.start();
  fftw_execute_dft_c2r(reverse_fft, reinterpret_cast<fftw_complex*>(uv_data), img);
  fft_convolve(img);
  fft_timer.stop();
}

void deconvolver::fft2b(pix * img, vis * uv_data){
  fft_timer.start();
  memset(uv_data, 0, height*width*sizeof(*uv_data));
  for (unsigned int j=0; j<uv_mask.size(); ++j)
    for (int inc=0; inc<height*width; ++inc)
      uv_data[uv_mask[j]] += img[inc]/fftw_scale * dft(inc, uv_mask[j]);
  fft_timer.stop();
}

void deconvolver::ifft2(vis * uv_data, pix * img){
  fft_timer.start();
  fftw_execute_dft_c2r(reverse_fft, reinterpret_cast<fftw_complex*>(uv_data), img);
  fft_timer.stop();
}

void deconvolver::dftt2(vis * uv_data, pix * result){
  fft_timer.start();
  memset(result, 0, height*width*sizeof(*result));
  for (int inc=0; inc<height*width; ++inc)
    for (unsigned int j=0; j<uv_mask.size(); ++j)
      result[inc] += real(uv_data[uv_mask[j]] * dft(inc, uv_mask[j]));
  fft_timer.stop();
}

vis deconvolver::dft(unsigned int col, unsigned int row){
  fl angle = -2*M_PI*((col%width)*(row%height)+(col/width)*(row/height))*fftw_scale;
  return std::polar(1., angle)*fftw_scale;
}

vis deconvolver::transposeDft(unsigned int col, unsigned int row){
  return dft(row, col);
}

vis deconvolver::dftProduct(unsigned int col1, unsigned int row1, unsigned int col2, unsigned int row2){
  fl angle = 2*M_PI*((col1%width)*(row1%height)+(col2%width)*(row2%height))*fftw_scale;
  return std::polar(1., angle)*scale;
}

fl deconvolver::realDftProduct(unsigned int col1, unsigned int row1, unsigned int col2, unsigned int row2){
  return std::cos(2*M_PI*((col1%width)*(row1%height)+(col1/width)*(row1/height)+(col2%width)*(row2%height)+(col2/width)*(row2/height))*fftw_scale)*scale;
}

void deconvolver::transpose(vis * uv_data){
  for (int i=0; i<width; ++i)
    for (int j=0; j<i; ++j){
      vis tmp = uv_data[j*width + i];
      uv_data[j*width + i] = uv_data[i*width + j];
      uv_data[i*width + j] = tmp;
    }
}

void deconvolver::fft_convolve(pix * uv_data){
  for (int i=1; i<width; ++i){
    for (int j=1; j<height; ++j){
      if (j*width + i < height*(width+1)/2){
        pix tmp = uv_data[j*width + i];
        uv_data[j*width + i] = uv_data[height*(width+1) - (j*width + i)];
        uv_data[height*(width+1) - (j*width + i)] = tmp;
      }
    }
  }
}

void deconvolver::convolve_cycle(pix * img){
  pix * tmp = (pix *) fftw_malloc(height*width*sizeof(pix));
  memcpy(tmp, img, width*height*sizeof(*tmp));
  for (int i=0; i<width; ++i){
    for (int j=0; j<height; ++j){
      img[j*width + i] = tmp[((j+height/2)%height)*width + (i+width/2)%width];
    }
  }
  fftw_free(tmp);
}

double deconvolver::time(){
  return total_timer.time();
}

double deconvolver::fft_time(){
  return fft_timer.time();
}

double deconvolver::matrixMult_time(){
  return matrix_timer.time();
}

double deconvolver::leastSquares_time(){
  return ls_timer.time();
}