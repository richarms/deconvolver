#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "transform.h"

#define INV_SQRT_2 0.707106781f
#define CENTER 128
#define min(a, b) (a < b ? a : b)
#define min3(a, b, c) (a < b ? (a < c ? a : c) : (b < c ? b : c))

pix * transform::haarTransform(pix * img) {
  pix * intermediate = new pix[width * height];
  pix * output = new pix[width * height];
  for (int i = 0; i < width*height; ++i)
    intermediate[i] = img[i];
  for (unsigned int iteration=0; iteration < min3(log2(width), log2(height), times); ++iteration){
    int area_height = height / (1 << iteration);
    int area_width = width / (1 << iteration);
    std::cout << area_height << std::endl;
    for (int i=0; i < area_height; ++i)
      for (int j=0; j < area_width/2; ++j) {
        output[i*width + j] = (intermediate[i*width + 2*j] + intermediate[i*width + 2*j + 1])/2;
        output[i*width + area_width/2 + j] = CENTER + (intermediate[i*width + 2*j] - intermediate[i*width + 2*j + 1])/2;
      }
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height/2; ++j) {
        intermediate[j*width + i] = (output[2*j*width + i] + output[(2*j+1)*width + i])/2;
        intermediate[(j + area_height/2)*width + i] = CENTER + (output[2*j*width + i] - output[(2*j+1)*width + i])/2;
      }
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height; ++j)
      output[j*width + i] = intermediate[j*width + i];
  }
  delete [] intermediate;
  return output;
}

pix * transform::haarReconstruct(pix * img) {
  pix * intermediate = new pix[width * height];
  pix * output = new pix[width * height];
  for (int i = 0; i < width*height; ++i)
    output[i] = img[i];
  for (int iteration=min3(log2(width), log2(height), times)-1; iteration >= 0; --iteration){
    int area_height = height / (1 << iteration);
    int area_width = width / (1 << iteration);
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height/2; ++j) {
        intermediate[j*2*width + i] = output[j*width + i] + (output[(j+area_height/2)*width + i] - CENTER);
        intermediate[(j*2+1)*width + i] = output[j*width + i] - (output[(j+area_height/2)*width + i] - CENTER);
      }
    for (int i=0; i < area_height; ++i)
      for (int j=0; j < area_width/2; ++j) {
        output[i*width + j*2] = intermediate[i*width + j] + (intermediate[i*width + area_width/2 + j] - CENTER);
        output[i*width + j*2 + 1] = intermediate[i*width + j] - (intermediate[i*width + area_width/2 + j] - CENTER);
      }
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height; ++j)
        intermediate[j*width + i] = output[j*width + i];
  }
  delete [] intermediate;
  return output;
}

pix * transform::liftingHaarTransform(pix * img) {
  pix * intermediate = new pix[width * height];
  pix * output = new pix[width * height];
  for (int i = 0; i < width*height; ++i)
    intermediate[i] = img[i];
  for (unsigned int iteration=0; iteration < min3(log2(width), log2(height), times); ++iteration){
    int area_height = height / (1 << iteration);
    int area_width = width / (1 << iteration);
    std::cout << area_height << std::endl;
    for (int i=0; i < area_height; ++i)
      for (int j=0; j < area_width/2; ++j) {
        output[i*width + area_width/2 + j] = intermediate[i*width + 2*j + 1] - intermediate[i*width + 2*j];
        output[i*width + j] = intermediate[i*width + 2*j] + output[i*width + area_width/2 + j]/2;
      }
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height/2; ++j) {
        intermediate[(j + area_height/2)*width + i] = output[(2*j+1)*width + i] - output[2*j*width + i];
        intermediate[j*width + i] = output[2*j*width + i] + intermediate[(j + area_height/2)*width + i]/2;
      }
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height; ++j)
      output[j*width + i] = intermediate[j*width + i];
  }
  delete [] intermediate;
  return output;
}

pix * transform::liftingHaarReconstruct(pix * img) {
  pix * intermediate = new pix[width * height];
  pix * output = new pix[width * height];
  for (int i = 0; i < width*height; ++i)
    output[i] = img[i];
  for (int iteration=min3(log2(width), log2(height), times)-1; iteration >= 0; --iteration){
    int area_height = height / (1 << iteration);
    int area_width = width / (1 << iteration);
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height/2; ++j) {
        intermediate[j*2*width + i] = output[j*width + i] - output[(j+area_height/2)*width + i]/2;
        intermediate[(j*2+1)*width + i] = output[(j+area_height/2)*width + i] + intermediate[j*2*width + i];
      }
    for (int i=0; i < area_height; ++i)
      for (int j=0; j < area_width/2; ++j) {
        output[i*width + j*2] = intermediate[i*width + j] - intermediate[i*width + area_width/2 + j]/2;
        output[i*width + j*2 + 1] = intermediate[i*width + area_width/2 + j] + output[i*width + j*2];
      }
    for (int i=0; i < area_width; ++i)
      for (int j=0; j < area_height; ++j)
        intermediate[j*width + i] = output[j*width + i];
  }
  delete [] intermediate;
  return output;
}