#include "image.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

image::image(string path) {
  img = new vis[1];
  ifstream image_file(path);
  image_file >> *this;
  image_file.close();
}

void image::save(string path) {
  ofstream image_file(path);
  image_file << *this;
  image_file.close();
}

void image::save_pgm(string path) {
  ofstream os(path);
  double min = abs(img[0]);
  double max = abs(img[0]);
  for (int i=1; i<height*width; ++i){
    max = (max > abs(img[i]) ? max : abs(img[i]));
    min = (min < abs(img[i]) ? min : abs(img[i]));
  }
  os << "P5" << endl <<
      width << " " << height << endl <<
      255 << endl;
  for (int i=0; i<height*width; ++i)
    os.put((unsigned char)round(((abs(img[i]))-min)/(max-min)*255));
  os.flush();
  os.close();
}

void image::reallocate(int h, int w) {
  delete [] img;
  height = h;
  width = w;
  img = new vis[height*width];
  for (int i=0; i<height*width; ++i)
    img[i] = 0;
}

double image::totalDifference(image other) const {
  if (height != other.getHeight() || width != other.getWidth())
    return -1;
  double sum = 0;
  for (int i=0; i < width*height; ++i)
    sum += abs(img[i] - other[i]);
  return sum;
}

double image::l2Difference(image other) const {
  if (height != other.getHeight() || width != other.getWidth())
    return -1;
  double sum = 0;
  for (int i=0; i < width*height; ++i)
    sum += pow(abs(img[i] - other[i]), 2);
  return sqrt(sum);
}

ostream & operator<<(ostream & os, const image & p){
  //overloads output stream opreator "<<"
  os << "Complex_image_data" << endl <<
      p.width << " " << p.height << endl;
  os.write((char *)(&p.img[0]), p.height*p.width*sizeof(p.img[0])/sizeof(char));
  os.flush();
  return os;
}

istream & operator>>(istream & is, image & p){
  //overloads input stream opreator ">>"
  string ignore;
  int width, height;
  is >> ws >> ignore >> ws;
  if (is.peek() == '#')
    getline(is, ignore);
  is >> ws >> width >> ws >> height >> ws;
  p.reallocate(height, width);
  is.read((char *)(&p.img[0]), height*width*sizeof(p.img[0])/sizeof(char));
  return is;
}
