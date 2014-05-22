//-*-mode:c++; mode:font-lock;-*-

#include "grating.h"

Grating::Grating()
{
  fAngleBlazeRad = kBlazeAngleDeg * (M_PI/180.); //in radians
  fIntegral.resize(kMaxNumIntervals);
  fNumSlits=kNumSlits;
  fDeltaNM=kDeltaNM;
}
// ************************************************************************

Grating::Grating(double BlazeAngleDeg, int NumSlits, double DeltaNM)
{
  fAngleBlazeRad =  BlazeAngleDeg * (M_PI/180.); //in radians
  fNumSlits=NumSlits;
  fDeltaNM=DeltaNM;
  fIntegral.resize(kMaxNumIntervals);

}
// ************************************************************************

Grating::~Grating()
{
  //Nothing yet
}
// ************************************************************************

double Grating::calculateFunction(double angleInRad, double angleOutRad,
                  double wavelengthNM)
//*************************************************
//    calculate the value of the intensity function
//**************************************************
{
  double v_prime;
  double v;
  v_prime = M_PI*fDeltaNM * (sin(angleOutRad) + sin(angleInRad)) / wavelengthNM;
  v = M_PI*fDeltaNM*cos(fAngleBlazeRad) * (sin(angleOutRad - fAngleBlazeRad) + sin(angleInRad - fAngleBlazeRad) ) / wavelengthNM;
  return pow(sin(fNumSlits * v_prime) / (fNumSlits * sin(v_prime)),2.0) * pow(sin(v)/v,2.0);
}
// ***************************************************************************

void Grating::makeTable(double angleInRad,  double wavelengthNM)
//*********************************************************
//    make a table of integration based on the angleOutRad. 
//    so far I don't know if there's a more efficient way
//    to do the integration. For the numerical method of 
//    integration, I think the bottleneck is actually from 
//    some of the operations like sin or square. Some other
//    numerical method may have less kStepRad for each wavelengthNM
//    but may result in more inefficient operations. Anyway,
//    I'll try other methods later.
//*********************************************************
{
  double angleOutRad;
  //double center = 2 * fAngleBlazeRad - angleInRad;
  //  double center = - angleInRad;
  int i;
    
  fIntegral.clear();
  fIntegral.resize(kMaxNumIntervals);
  
  //angleOutRad = center - (M_PI/2.0);
  angleOutRad = -(M_PI/2.0);
  for(i = 0; i < kMaxNumIntervals; i++){
    angleOutRad += kStepRad;
    if(i == 0){
      fIntegral.at(i) = kStepRad * 
                 calculateFunction(angleInRad,angleOutRad,wavelengthNM);
    }
    else{
      fIntegral.at(i) =fIntegral.at(i-1) + kStepRad * 
                 calculateFunction(angleInRad,angleOutRad,wavelengthNM);
    }
  }
  return;
}
// *************************************************************************

int Grating::binarySearch(double goal)
// **************************************************************************
//  return the index of the interval which correspond to the upper_limit of 
//  the integration.
// **************************************************************************
{ 
  int l = 0;
  int r = kMaxNumIntervals - 1;
  int mid;
    
  while(l <= r){
    mid = l + ((r - l) >> 1);
    if(goal > fIntegral.at(mid) ){
      l = mid + 1;
    }
    else {
      r = mid - 1;
    }
  }
  return mid;
}
// *************************************************************************



double Grating::calculateAngle(double angleInRad,double wavelengthNM)
// *********************************************************
// calculate the angleOutRad of each photon after grating
// *********************************************************
{
  double angleOutRad;
  double upper_limit;
  //double center = 2 * fAngleBlazeRad - angleInRad;
  //double center = - angleInRad;
    
  makeTable(angleInRad,wavelengthNM);
  upper_limit = rand()/(RAND_MAX + 1.0) * fIntegral.at(kMaxNumIntervals-1);
  
  angleOutRad = -(M_PI/2.0)+kStepRad * binarySearch(upper_limit);
  return angleOutRad;
}
// ********************************************************************




void Grating::diffract(double vxIn, double vyIn, double vzIn, 
               double vxGratingNormal, double vyGratingNormal,
               double vzGratingNormal, double& vxOut,
               double& vyOut, double& vzOut, double wavelengthNM)
// *********************************************************
// Main method: Calculates the direction out of a photon after it interacts
// with the grating.
// **********************************************************
// vxIn,vyIn,vzIn is input direction of the photon
// vxOut,vyOut,vzOut is output direction of the photon.
{

  setGratingNormal(vxGratingNormal, vyGratingNormal, vzGratingNormal);


  // **************************************************************************
  // To use the grating, we need at first transform the phothons from the lab
  // frame to the optic frame, and then call THIS function. After that,
  // we need to transform the coordinates back to the lab frame.
  // (call transform_inverse)
  // **************************************************************************
    
  // **********************************************************
  //     o_hat = a * i_hat + b * n_hat;
  //     vnDotvi = n_hat * i_hat = cos(angleInRad)
  //     I'm not quite sure about these geometries because so far 
  //   I didn't test it in the whole LSST simulation. Please help
  //   me double check that.
  // ******************************************************
  // GHS: We need negative of vin vector to get angleInRad
  // Put a test in here later to see if photon is hitting bottom of grating
  // ***************************************************************

  double vnDotvi = -(vxIn*fVxGratingNormal + vyIn*fVyGratingNormal + 
            vzIn*fVzGratingNormal); //This should be just -vzIn if N=0,0,1
  //double vnDotvi = (vxIn*fVxGratingNormal + vyIn*fVyGratingNormal + 
  //            vzIn*fVzGratingNormal);//todo
  double angleInRad = acos(vnDotvi);    
    
  // *******************************************************************
  // Change sign on angleout for our calculation. Put in test later.
  //If it comes back positive there was a problem
  // *******************************************************************

  double angleOutRad = calculateAngle(angleInRad, wavelengthNM);
  
  // *******************************************************************
  // Now generate the outgoing vector
  // *******************************************************************
  double a=0;
  double b=0;
  //if(angleOutRad>0){        //Negative mode (m<0 in grating eq.)
  //  a = -sin(-angleOutRad)/sin(angleInRad);
  //  b = +cos(angleOutRad)+a*cos(angleInRad);
  // }
  //else{
    a = sin(-angleOutRad) / sin(angleInRad);
    b = cos(-angleOutRad) + a*cos(angleInRad);
    // }
  
  vxOut = a*vxIn + b*fVxGratingNormal;
  vyOut = a*vyIn + b*fVyGratingNormal;
  vzOut = a*vzIn + b*fVzGratingNormal;

  // *************************
  // Debug prints
  // *************************
  //double mode=(sin(angleInRad)+sin(angleOutRad))*(fDeltaNM/wavelengthNM);
  //int m;
  //if (mode<0){
  //  m=mode-.5;
  //}
  //else{
  //   m=mode+.5;
  //}
  //std::cout<<m<<" "<<wavelengthNM<<" "<<angleInRad<<" "<<angleOutRad<<" "<<fAngleBlazeRad<<std::endl;
  return;
}

