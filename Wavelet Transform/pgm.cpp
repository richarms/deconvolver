#include "pgm.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

pgm::pgm(string path) {
  img = new pix[1];
  ifstream image_file(path);
  image_file >> *this;
  image_file.close();
}

void pgm::save(string path) {
  ofstream image_file(path);
  image_file << *this;
  image_file.close();
}

void pgm::reallocate(int h, int w) {
  delete [] img;
  height = h;
  width = w;
  img = new pix[height*width];
  for (int i=0; i<height*width; ++i)
    img[i] = 0;
}

double pgm::totalDifference(pgm other) const {
  if (height != other.getHeight() || width != other.getWidth())
    return -1;
  double sum = 0;
  for (int i=0; i < width*height; ++i)
    sum += abs(img[i] - other[i]);
  return sum;
}

ostream & operator<<(ostream & os, const pgm & p){
  //overloads output stream opreator "<<"
  os << "P5" << endl <<
      p.width << " " << p.height << endl <<
      255 << endl;
  for (int i=0; i<p.height*p.width; ++i)
    os.put(p.img[i]);
  os.flush();
  return os;
}

istream & operator>>(istream & is, pgm & p){
  //overloads input stream opreator ">>"
  string ignore;
  int width, height;
  is >> ws >> ignore >> ws;
  if (is.peek() == '#')
    getline(is, ignore);
  is >> ws >> width >> ws >> height >> ws;
  getline(is, ignore);
  p.reallocate(height, width);
  for (int i=0; i<height*width; ++i)
    p.img[i] = is.get();
  return is;
}
