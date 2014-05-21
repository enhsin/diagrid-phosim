///
/// @package phosim
/// @file instrument.cpp
/// @brief instrument  application.
///
/// @brief Created by
/// @author Nathan Todd (Purdue)
///
/// @brief Modified by
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
/// @author Glenn Sembroski (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "instrument.h"

int getModeFromTag(std::string tag)
// Find the index for this tag
{ int num=-1;
    //The allowed tags.
    if      (tag=="phi")   num=1;
    else if (tag=="psi")   num=2;
    else if (tag=="theta") num=3;
    else if (tag=="xdis")  num=4;
    else if (tag=="ydis")  num=5;
    else if (tag=="zdis")  num=6;
    else if (tag=="z0")    num=7;
    else if (tag=="z1")    num=8;
    else if (tag=="z2")    num=9;
    else if (tag=="z3")    num=10;
    else if (tag=="z4")    num=11;
    else if (tag=="z5")    num=12;
    else if (tag=="z6")    num=13;
    else if (tag=="z7")    num=14;
    else if (tag=="z8")    num=15;
    else if (tag=="z9")    num=16;
    else if (tag=="z10")   num=17;
    else if (tag=="z11")   num=18;
    else if (tag=="z12")   num=19;
    else if (tag=="z13")   num=20;
    else if (tag=="z14")   num=21;
    else if (tag=="z15")   num=22;
    else if (tag=="z16")   num=23;
    else if (tag=="z17")   num=24;
    else if (tag=="z18")   num=25;
    else if (tag=="z19")   num=26;
    else if (tag=="z20")   num=27;
    else if (tag=="z21")   num=28;
    // Note that we need to have gTotalNumTags equal to the
    // highest value num can take here. (Set in instrument.h)
    // If you add more tags, be sure to increase gTotalNumTags
    return num;
}

void updateControl(PhosimParser& control,PhosimParser& pars)
//updates control PhosimParser map with values from input command file in pars
// May add or may replace existing entries.
{
    // First we need to search through all the pars keys for the control type
    // They will be the only keys with a space " " in them. See PhosimParser
    // readCommandStream
    int numKeys=pars.getNumberKeys();
    std::string actuatorValue;
    std::string actuatorKey;
    bool haveActuator=false;

    for (int i=0;i<numKeys;i++){
        std::string key;
        bool goodKey=pars.getKey(i,key);
        if(goodKey){
            std::string value=pars[key];
            // Search the key for a space
            std::string::size_type idx=key.find(" ");
            if(idx!= std::string::npos){
                // Check for the special case of there being an "actuator" command
                // in the pars file. If there is re enter it in the pars map with the
                // proper key instead of the concatinated key. Obviously don't add it
                // to the control pars map.
                // Could also be a body command Those are handled in overrideBofy function
                // Otherwise assume its  a control.txt kind of command, in
                // which case add it to the control.txt PhosimPar control

                std::string cmdKey;
                std::istringstream iss(key);
                iss>>cmdKey;
                if (cmdKey == "actuator"){
                    haveActuator=true;
                    std::string line;
                    iss>>line;
                    actuatorValue=line + " " + value;
                    actuatorKey=key;
                }
                // See if this is not a body command it must be a control command
                else if(cmdKey != "body"){
                    // Add (or replace) this control command.
                    control.removeKey(key);
                    control.set(key,value);
                }
            }
        }
        else{
            std::cout<<"#Error: Unexpected(never should happen) error in updateControl!"
                     <<std::endl;
            return;  //opps should never be an error here.
        }
    }
    //We delay setting the new actuator key in pars so that we don't change 
    // the map while in the middle of acessing it in the above for loop.
    if(haveActuator){
        pars.removeKey(actuatorKey);
        actuatorKey="actuator";
        pars.removeKey(actuatorKey);
        pars.set(actuatorKey,actuatorValue);
    }

    return;
}

void createBody(std::map< int, std::vector<double>*>&  bodyMap,
                PhosimParser& control,
                PhosimParser& pars,
                std::vector < std::vector < double> >& actuatorMatrix,
                std::vector < double >& actuatorDistance,
                std::vector < double >& actuatorError)

// Fill the body vectors using the commands in the control PhosimPars map and
// the actuator arrays.
{
    // Use a map here for the body data to allow for any number of different devices in
    // any order appearing in the control.txt file. We can add them as they show up
    // Note we are using pointers to vectors on the heap for the data in the body map.

    //Init the random number generator
    int seed=(int)pars["obsseed"];
    if (seed==-1){
        RngSetSeedFromTime();
    }
    else{
        RngSetSeed32(seed);
    }
    RngUnwind(10000);
    RngSetSeed32_reseed(1000);

    // Get all the values we will need to determine the body values
    double scaleOne=1.0;
    if( pars.has_key("scalefactorone")){
        scaleOne=pars["scalefactorone"];
    }

    double scaleTwo=1;
    if( pars.has_key("scalefactortwo")){
        scaleTwo=pars["scalefactortwo"];
    }

    double seeing=0.0;
    if( pars.has_key("constrainseeing")){
        seeing=pars["constrainseeing"];
    }

    double domeseeing=0.1;
    if( pars.has_key("domeseeing")){
        domeseeing=pars["domeseeing"];
    }

    double zenithRad=0.0;
    if( pars.has_key("zenith")){
        double zenithDeg=pars["zenith"];
        zenithRad = zenithDeg*kPi/180.;
    }

    double temperature0=20.0;
    double temperatureC=20.0;
    if( pars.has_key("temperature")){
        temperatureC=pars["temperature"];
    }

    double temperatureVarC=0.0;
    if( pars.has_key("tempvar")){
        temperatureVarC=pars["tempvar"];
    }

    double pressure0=520.0;
    double pressuremmHg=520.0;
    if( pars.has_key("pressure")){
        pressuremmHg=pars["pressure"];
    }

    double pressureVarmmHg=0.0;
    if( pars.has_key("pressvar")){
        pressureVarmmHg=pars["pressvar"];
    }
    //elevation
    double altitudeDeg=90.0;
    if( pars.has_key("altitude")){  //Which (zenith or alititude) has precedence?
        altitudeDeg=pars["altitude"];
        zenithRad = (90.-altitudeDeg)*kPi/180.;
    }

    double altitudeVarDeg=0.0;
    if( pars.has_key("altvar")){
        altitudeVarDeg=pars["altvar"];
    }

    // Iterate through the control map
    int numKeys=control.getNumberKeys();
    for (int keyNum=0;keyNum<numKeys;keyNum++){
        std::string key;
        bool gotKey=control.getKey(keyNum,key);
        if(!gotKey){
            std::cout<<"#Error in createBody"<<std::endl;
            return;
        }

        // Parse the key
        std::string deviceStr;
        std::string tag;
        std::istringstream iss(key);
        iss>>deviceStr;
        iss>>tag;
        int mode=getModeFromTag(tag);
        if(mode==-1){
            std::cout<<"Illegal mode in control.txt or command file: "<<key
                     <<std::endl;
            continue;
        }
        // Ok, we can now fill in an optics vector though we will throw most of it
        // away.
        std::vector< double > optics;
        optics.resize(gDevCol,0);

        optics.at(0)=mode;

        std::string value=control[key];
        std::istringstream isss(value);
        double bodyValue;

        for(int i=1;i<gDevCol;i++){
            isss>>bodyValue;
            optics.at(i)=bodyValue;
        }

        //Above got us the "mandatory" stuff, now get any possible extra stuff
        // Note that the next values are first the device id for this line
        // followed by the device id's that are "connected" to this device
        // (ex 2 sides of a lens)
        while(isss>>bodyValue){
            optics.push_back(bodyValue);
        }

        // Now we get into the filling of the body array

        double bodyQuantity=0.0;

        // Fill the body as needed
        // I don't know what mostof this stuff is
        if ((int)optics.at(1)==1){
            bodyQuantity=optics.at(2)*random_gaussian()+optics.at(3);
        }
        else if((int)optics.at(1)==2){
            bodyQuantity=(optics.at(3)-optics.at(2))*RngDouble()+optics.at(2);
        }

        if ((int)optics.at(4)==1){
            bodyQuantity*=scaleOne*sqrt(seeing*seeing+domeseeing*domeseeing)/
                0.67*pow(1/cos(zenithRad),3./5.);
        }
        else if ((int)optics.at(4)==2) {
            bodyQuantity*=scaleTwo;
        }
        //fabriaction error
        bodyQuantity+=optics.at(5)*RngDouble_reseed();

        bodyQuantity+=                                        //thermal
            optics.at(6)*(temperatureC-temperature0+temperatureVarC*random_gaussian());

        bodyQuantity+=                                        //pressure
            optics.at(7)*(pressuremmHg-pressure0+pressureVarmmHg*random_gaussian());

        bodyQuantity+=optics.at(8)*(altitudeDeg+altitudeVarDeg*random_gaussian());

        bodyQuantity+=optics.at(9)*random_gaussian();         //hidden variable

        // Now apply the actuator stuff
        int numActuator=actuatorError.size();
        for (int i=0;i<numActuator;i++){
            bodyQuantity+=
                actuatorMatrix.at(keyNum).at(i)*(actuatorDistance.at(i)+
                                                 actuatorError.at(i)*random_gaussian());
        }


        int modeIndex=(int)(optics.at(0)-1);  // for things like xdis z12 phi etc.

        int devIndex =(int)optics.at(gDevCol); // for things like M1 M13 camera etc

        // See if we already have this device in our pBodyMap, if not make
        // a new one
        std::vector < double >*  pBodyDevice;
        std::map< int, std::vector< double >* >::iterator bodyPos =
            bodyMap.find(devIndex);
        if(bodyPos==bodyMap.end()){
            // Need to make a new body vector on the heap for this device and get a
            // pointer to it
            // I'm worried about the hardwiring of 28 here.
            pBodyDevice = new std::vector< double >;
            pBodyDevice->resize(gTotalNumTags,-5.0);
            bodyMap[devIndex]=pBodyDevice;
        }
        else{
            pBodyDevice=bodyPos->second;
        }

        pBodyDevice->at(modeIndex)=bodyQuantity;

        //Once body array is done see if we this is a zernike type command
        // and this device is "attached" to others (so that they have the same
        // body commands. See above for comments on how this works.
        int numOptics=optics.size();
        if(numOptics-1>=gDevCol+1){
            for(int k=gDevCol+1;k<numOptics;k++) {
                int zDevIndex =(int)optics.at(k);


                std::map< int, std::vector< double >* >::iterator zBodyPos =
                    bodyMap.find(zDevIndex);
                if(zBodyPos==bodyMap.end()){
                    pBodyDevice = new std::vector< double >;
                    pBodyDevice->resize(gTotalNumTags,-5.0);
                    bodyMap[zDevIndex]=pBodyDevice;
                }
                else{
                    pBodyDevice=zBodyPos->second;
                }
                pBodyDevice->at(modeIndex)=bodyQuantity;
            }
        }

    }
    return;
}

int main(void) {

    std::cout<<"-------------------------------------------------------------"
        "-----------------------------"<<std::endl;
    std::cout<<"Instrument Configuration"<<std::endl;
    std::cout<<"-------------------------------------------------------------"
        "-----------------------------"<<std::endl;
    PhosimParser pars;

    // Set some default values.
    pars.set("instrdir", "../data/lsst");
    pars.set("outputdir", ".");
    pars.set("actuator","\0");
    pars.set("vistime","15.0");

    // Read parameters from stdin.
    int parsNumberExtra=0;
    pars.readCommandStream(std::cin,parsNumberExtra);

    // Read optics_0.txt file to get map of surface names to surface ID
    instrumentFiles.makeSurfaceMap(pars);

    //Create and write the "tracking" file.
    instrumentFiles.makeTrackingFile(pars);

    PhosimParser controlPars;
    instrumentFiles.readControlFile(pars,controlPars);

    // Update the control map with any control commands that were found in the
    // input command file.
    // These values will over ride any existing or add to the body commands in
    // control
    if(parsNumberExtra>0){
        updateControl(controlPars,pars);
    }

    // Readf in the actuator file
    std::vector < std::vector < double> > actuatorMatrix;
    std::vector < double >  actuatorDistance;
    std::vector < double >  actuatorError;

    instrumentFiles.readActuatorFile(pars, controlPars, actuatorMatrix,
                                     actuatorDistance, actuatorError);
    // Define the body matrix map.
    // Remember to always refer to bodyMap by reference in subroutine calls for
    // efficency/performance
    std::map< int, std::vector<double>*>  bodyMap;

    // Now go through our control PhosimParser map and fill the body array
    createBody(bodyMap,controlPars, pars, actuatorMatrix, actuatorDistance,actuatorError);

    instrumentFiles.writeBodyFile(bodyMap, pars);

    //readout, hot pixel
    instrumentFiles.readoutPars(pars);

    //body, chipangle, izernike, qevariation
    instrumentFiles.focalPlanePars(pars);

    return 0;
}
