#include <iostream>
#include <fstream>
#include "pgm.h"
#include "transform.h"

using namespace std;

int main(void) {
  cout << "READ FILE:" << endl;
  pgm image = pgm("data/lena.pgm");
  cout << "TRANSFORMING IMAGE:" << endl;
  transform t = transform(image.getHeight(), image.getWidth());
  pgm transformed_image = pgm(image.getHeight(), image.getWidth(), t.liftingHaarTransform(image.data()));
  cout << "RECONSTRUCTING IMAGE:" << endl;
  pgm output_image = pgm(image.getHeight(), image.getWidth(), t.liftingHaarReconstruct(transformed_image.data()));
  cout << "ERROR: " << output_image.totalDifference(image) << endl;
  cout << "WRITE FILE:" << endl;
  transformed_image.save("transform.pgm");
  output_image.save("output.pgm");
}