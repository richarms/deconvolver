#ifndef PGM_H_
#define PGM_H_

#include <iostream>
#include <string>
#include "types.h"

using namespace std;

//main class containing pgm functions
class pgm {
  public:
    //default constructor
    pgm(int h=10, int w=10) : height(h), width(w), img(new pix[height*width]) {
      for (int i=0; i<height*width; ++i)
        img[i] = 0;
    }

    //raw data constructor
    pgm(int h, int w, pix * i) : height(h), width(w), img(i) { }

    //file constructor
    pgm(std::string path);
  
    //destructor
    ~pgm() { delete [] img; }

    //copy constructor
    pgm(const pgm & cpfrm) : height(cpfrm.getHeight()), width(cpfrm.getWidth()), img(new pix[height*width]){
      for (int i=0; i<height*width; ++i)
        img[i] = cpfrm[i];
    }

    //copy assignment constructor
    pgm & operator=(const pgm & cpfrm){
      if (this != &cpfrm){
        height = cpfrm.getHeight();
        width = cpfrm.getWidth();
        delete [] img;
        img = new pix[height*width];
        for (int i=0; i<height*width; ++i)
          img[i] = cpfrm[i];
      }
      return *this;
    }

    void save(std::string path);

    //subscript operator - returns pixel at position in array
    pix & operator[](int pos) const { return img[pos]; }

    //returns pixel at (x, y)
    pix getpix (int x, int y) const { return img[x + y*width]; }

    //sets pixel at (x, y)
    void setpix(int x, int y, pix c) { img[x + y*width] = c; }

    //returns the height
    int getHeight() const { return height; }

    //returns the width
    int getWidth() const { return width; }
    
    //returns the total pixel difference with another pgm image
    double totalDifference(pgm other) const;

    //gets the raw imgae data
    pix * data() const { return img; }

    //friend operators extend iostream "<<" and ">>" operators
    friend ostream & operator<<(ostream & os, const pgm & p);
    friend istream & operator>>(istream & is, pgm & p);

  private:
    void reallocate(int h, int w);

    int height, width;
    pix * img;
};

ostream & operator<<(ostream & os, const pgm & p);
istream & operator>>(istream & is, pgm & p);

#endif /* PGM_H_ */
