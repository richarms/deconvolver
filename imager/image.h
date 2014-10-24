#ifndef PGM_H_
#define PGM_H_

#include <iostream>
#include <string>
#include "types.h"

using namespace std;

//main class containing image functions
class image {
  public:
    //default constructor
    image(int h=10, int w=10) : height(h), width(w), img(new vis[height*width]) {
      for (int i=0; i<height*width; ++i)
        img[i] = 0;
    }

    //raw data constructor
    image(int h, int w, vis * source) : height(h), width(w) {
      img = new vis[height*width];
      for (int i=0; i<height*width; ++i)
        img[i] = source[i];
    }

    //real raw data constructor
    image(int h, int w, pix * reals) : height(h), width(w), img(new vis[height*width]){
      for (int i=0; i<height*width; ++i)
        img[i] = vis(reals[i], 0);
    }

    //file constructor
    image(std::string path);
  
    //destructor
    ~image() { delete [] img; }

    //copy constructor
    image(const image & cpfrm) : height(cpfrm.getHeight()), width(cpfrm.getWidth()), img(new vis[height*width]){
      for (int i=0; i<height*width; ++i)
        img[i] = cpfrm[i];
    }

    //copy assignment constructor
    image & operator=(const image & cpfrm){
      if (this != &cpfrm){
        height = cpfrm.getHeight();
        width = cpfrm.getWidth();
        delete [] img;
        img = new vis[height*width];
        for (int i=0; i<height*width; ++i)
          img[i] = cpfrm[i];
      }
      return *this;
    }

    void save(std::string path);

    void save_pgm(std::string path);

    //subscript operator - returns pixel at position in array
    vis & operator[](int pos) const { return img[pos]; }

    //returns pixel at (x, y)
    vis getpix (int x, int y) const { return img[x + y*width]; }

    //sets pixel at (x, y)
    void setpix(int x, int y, vis c) { img[x + y*width] = c; }

    //returns the height
    int getHeight() const { return height; }

    //returns the width
    int getWidth() const { return width; }
    
    //returns the total pixel difference with another image image
    double totalDifference(image other) const;

    //returns the L2 difference between images
    double l2Difference(image other) const;

    //gets the raw imgae data
    vis * data() const { return img; }

    //prevents clearing of data
    void no_free() { img = NULL; }

    //friend operators extend iostream "<<" and ">>" operators
    friend ostream & operator<<(ostream & os, const image & p);
    friend istream & operator>>(istream & is, image & p);

  private:
    void reallocate(int h, int w);

    int height, width;
    vis * img;
};

ostream & operator<<(ostream & os, const image & p);
istream & operator>>(istream & is, image & p);

#endif /* PGM_H_ */
