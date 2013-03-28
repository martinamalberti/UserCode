///==== include ====


#include <cmath>

float BSweight(float diff){
  
  //--- parameters for beam spot reweighting (Matt):https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/993.html
  
  // --- parameters for 4.8 cm target BS
  float newBSmean1  = 9.9391e-02;
  float newBSmean2  = 1.8902e-01;
  float newBSnorm1  = 5.3210e+00;
  float newBSnorm2  = 4.1813e+01;
  float newBSsigma1 = 9.7530e-01;
  float newBSsigma2 = 7.0811e+00;
  
  float oldBSmean1  = 7.2055e-02;
  float oldBSmean2  = 4.9986e-01;
  float oldBSnorm1  = 3.5411e+00;
  float oldBSnorm2  = 4.0258e+01;
  float oldBSsigma1 = 7.9678e-01;
  float oldBSsigma2 = 8.5356e+00;
  // use just 1 gaussian for Z->mumu as there are no conversions
  //  newBSnorm1=0.;
  //  oldBSnorm1=0.;

  // --- parameters for 5.0 cm target BS
  //   float newBSmean1  = 9.549e-02;
  //   float newBSmean2  = 2.334e-01;
  //   float newBSnorm1  = 5.067e+00;
  //   float newBSnorm2  = 4.159e+01;
  //   float newBSsigma1 = 9.498e-01;
  //   float newBSsigma2 = 7.289e+00;
  
  //   float oldBSmean1  = 7.2055e-02;
  //   float oldBSmean2  = 4.9986e-01;
  //   float oldBSnorm1  = 3.5411e+00;
  //   float oldBSnorm2  = 4.0258e+01;
  //   float oldBSsigma1 = 7.9678e-01;
  //   float oldBSsigma2 = 8.5356e+00;

  // --- parameters for beam spot reweighting derived directly on Z->mumu data
  //   float newBSmean1  = 0.;
  //   float newBSmean2  = -0.013;
  //   float newBSnorm1  = 0.;
  //   float newBSnorm2  = 0.005875;
  //   float newBSsigma1 = 1;
  //   float newBSsigma2 = 6.951;
  
  //   float oldBSmean1  = 0;
  //   float oldBSmean2  = 0.02159;
  //   float oldBSnorm1  = 0.;
  //   float oldBSnorm2  = 0.004739;
  //   float oldBSsigma1 = 1;
  //   float oldBSsigma2 = 8.577;

  float newBSgaus1,  newBSgaus2 , oldBSgaus1 , oldBSgaus2, w; 

  w = 1.;
  
  // -- re-weight only "wrong vertices"
  if (fabs(diff)>0.1){
    newBSgaus1 = newBSnorm1*exp(-0.5*pow((diff-newBSmean1)/newBSsigma1,2));
    newBSgaus2 = newBSnorm2*exp(-0.5*pow((diff-newBSmean2)/newBSsigma2,2));
    oldBSgaus1 = oldBSnorm1*exp(-0.5*pow((diff-oldBSmean1)/oldBSsigma1,2));
    oldBSgaus2 = oldBSnorm2*exp(-0.5*pow((diff-oldBSmean2)/oldBSsigma2,2));
    w =  1.1235 *(newBSgaus1+newBSgaus2)/(oldBSgaus1+oldBSgaus2); // 2 gaussians, 4.8 cm
    //bsweight =  1.13242 *(newBSgaus1+newBSgaus2)/(oldBSgaus1+oldBSgaus2); // if using 1 gaussian
    //bsweight =  1.10332 *(newBSgaus1+newBSgaus2)/(oldBSgaus1+oldBSgaus2); // 2 gaussians , 5 cm
  }

  return (w);

}
