
///==== include ====
#include "ntpleUtils.h"


// simple electron ID WP95

float eleId95 (float pt, float tkIso, float emIso, float hadIso, float combIso, float isEB, float sigmaIetaIeta,
	     float DetaIn,float DphiIn, float HOverE, int mishits   )
{
  int myid = 0;
  
  if(isEB){
    
    myid = ( (pt            > 20   )  &&
	     (combIso/pt    < 0.15 )  && 
	     (sigmaIetaIeta < 0.010)  && 
	     (fabs(DphiIn)  < 0.8  )  && 
	     (fabs(DetaIn)  < 0.007)  &&  
	     (HOverE        < 0.15 )  && 
	     (mishits       < 2    )
	   ) ;
  }

  
  if(!isEB){
    
    myid = ( (pt            > 20   )  &&
	     (combIso/pt    < 0.10 )  && 
	     (sigmaIetaIeta < 0.030)  && 
	     (fabs(DphiIn)  < 0.7  )  && 
	     (fabs(DetaIn)  < 0.010)  &&  
	     (HOverE        < 0.07 )  && 
	     (mishits       < 2    )
	   ) ;
  }


  return (float(myid));

}

