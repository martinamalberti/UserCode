// c++ -o unit01 `root-config --glibs --cflags` -lm hFactory.cc hChain.cc unit01.cpp

#include "hFactory.h"
#include "hFunctions.h"

int main ()
{
  hFactory factory ("test.root", true) ;
  ///---- histogram with 200 bins, in the range [-10,10], total number of steps to fill = 3
  factory.add_h1 ("test", "jet tags eta prod", 200, -10, 10, 3) ;
  ///---- fill the step 0 of the histogram "test" with the value 1.
  factory.Fill ("test", 0, 1.) ;
  ///---- fill the step 1 of the histogram "test" with the value 2.5
  factory.Fill ("test", 1, 2.5) ;
  int colors[3] = {92,98,97} ;
  factory.applyToAll (setH1Colors (std::vector<int> (colors, colors + 3))) ;
}
