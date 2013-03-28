// c++ -o unit02 `root-config --glibs --cflags` -lm h2Factory.cc hChain.cc h2Chain.cc unit02.cpp

#include "h2Factory.h"
#include "hFunctions.h"

int main ()
{
  h2Factory factory ("test2.root", true) ;
  ///---- histogram with 6 bins, in the range [-3,3] (x), 4 bins, in the range [-2,2] (y),  total number of steps to fill = 2
  factory.add_h2 ("test", "jet tags eta prod", 6, -3., 3., 4, -2., 2., 2) ;
  ///---- fill the step 0 of the histogram "test" with the pair (1.,1.) = (x,y)
  factory.Fill ("test", 0, 1., 1.) ;
  ///---- fill the step 1 of the histogram "test" with the pair (2.,1.) = (x,y)
  factory.Fill ("test", 1, 2., 1.) ;
  int colors[2] = {9,98} ;
  factory.applyToAll (setH2Colors (std::vector<int> (colors, colors + 2))) ;
}