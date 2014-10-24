#ifndef TRANSFORM_H_
#define TRANSFORM_H_

#include <limits>
#include "types.h"

class transform {
  public:
    transform(int h, int w, int t = std::numeric_limits<int>::max()) : height(h), width(w), times(t) { };
    pix * haarTransform(pix * img);
    pix * haarReconstruct(pix * img);
    pix * liftingHaarTransform(pix * img);
    pix * liftingHaarReconstruct(pix * img);
  private:
    int height, width, times;
};

#endif