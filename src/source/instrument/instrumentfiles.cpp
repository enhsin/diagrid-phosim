///
/// @package phosim
/// @file instrumentfiles.cpp
/// @brief instrument  applications.
///
/// @brief Created by
/// @author Glenn Sembroski (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "instrumentfiles.h"

InstrumentFiles::InstrumentFiles() {
    //Nothing to do
}

InstrumentFiles::~InstrumentFiles() {
    //nothing to do
}

void InstrumentFiles::makeTrackingFile(std::string trackingFileName, double vistime, double jittertime) {

    std::ofstream outputDataFile(trackingFileName.c_str());

    double starttime = -0.5*vistime;

    double tempf1 = 0;
    double tempf2 = 0;
    double tempf3 = 0;
    double tempf4 = 0;
    double tempf5 = 0;
    double tempf6 = 0;
    double tempf7 = 0;
    double tempf8 = 0;
    double tempf9 = 0;
    long prevj = 0;
    long currj = 0;

    for (long j = 0; j < 10000; j++) {
        double time = starttime + vistime/10000.0*j;  //Step through the time
        if (j == currj) {
            tempf1 = RngDouble()*jittertime;
            if (tempf1 < vistime/10000.0){
                tempf1 = vistime/10000.0;
            }
            prevj = currj;
            currj = static_cast<long>(currj + (tempf1/vistime)*10000.0);
            if (currj == prevj) currj += 1;
            tempf5 = tempf2;
            tempf6 = tempf3;
            tempf7 = tempf4;
            tempf8 = tempf9;
            tempf2 = random_gaussian()*sqrt((jittertime/2.0)/vistime) + tempf5;
            tempf3 = random_gaussian()*sqrt((jittertime/2.0)/vistime) + tempf6;
            tempf4 = random_gaussian()*sqrt((jittertime/2.0)/vistime) + tempf7;
            tempf9 = random_gaussian()*sqrt((jittertime/2.0)/vistime) + tempf8;
        }
        double jitterele = tempf2*(j - prevj)/(currj - prevj) + tempf5*(currj - j)/(currj - prevj);
        double jitterazi = tempf3*(j - prevj)/(currj - prevj) + tempf6*(currj - j)/(currj - prevj);
        double jitterrot = tempf4*(j - prevj)/(currj - prevj) + tempf7*(currj - j)/(currj - prevj);
        double jitterwind = tempf9*(j - prevj)/(currj - prevj) + tempf8*(currj - j)/(currj - prevj);

        outputDataFile << std::fixed << std::setprecision(6) << time << " " << jitterele
                       << " " << jitterazi << " "<< jitterrot << " " << jitterwind << std::endl;
    }
    outputDataFile.close();
}


int InstrumentFiles::readControlFile(PhosimParser& pars, PhosimParser& controlPars) {
// *************************************************************************
// This method reads in the control.txt file to get the optics configuration
// Later this will be used to fill the body and zernike matices.
// Results are left in the pars map. Return numnber of entries in that map
// ************************************************************************

    // *********************************************************************
    // Open the control.txt files. Read it in first
    // We are going to do this the easy way and use rhwe PhosimParser class to
    // decode it. Note we had to add a comment ignoror to the PhosimParser class
    // *********************************************************************
    std::string directory = pars["instrdir"];
    std::string controlFileName = directory + "/control.txt";
    std::ifstream control(controlFileName.c_str());

    int numExtra = 0;
    controlPars.readCommandStream(control,numExtra);
    return numExtra;
}
// *************************************************************************


void InstrumentFiles::readActuatorFile(PhosimParser& pars,
                                       PhosimParser& controlPars,
                                       std::vector < std::vector < double> >& actuatorMatrix,
                                       std::vector < double >& actuatorDistance,
                                       std::vector < double >& actuatorError)

// *************************************************************************
// This method reads in the actuator.txt file to fill in the acutator arrays
// ************************************************************************
{
    std::string directory=pars["instrdir"];
    std::string actuatorFileName=directory +"/actuator.txt";
    std::ifstream actuator(actuatorFileName.c_str());
    // *******************************************************************
    // Use the PhosimParsar class method to read in this file. We have defined
    //  a readActuatorStream there.
    // ********************************************************************
    // First get the number of lines in the control.txt file. zernike, body and
    // camera lines combined that were found in the control.txt file
    // ********************************************************************
    int numLines = controlPars.getNumberKeys();

    // ************************************************
    // Now read in the actuator file. This file is suposed to have
    // numLines +2 values per line. That is , on value for each line in the
    // control.txt file plus an acutator distance and an error.
    // We will check that this is so. If not at least this many we will put out
    // an error message. However, we will take as duistance and error the last
    // 2 values in the line no matter how many values are listed. The missing
    // values in the acutator matrix will be set to 0.
    // ************************************************
    // Find number of lines in this file (so we can dimension our vectors)
    // I have no idea what the number of lines means. Semmes to be 20 for LSST
    // Save the 'good' lines in a string vector
    // ***********************************************
    std::vector< std::string > actuatorLines;
    std::string line;
    std::vector < int > numberOfTokens;
    while (std::getline(actuator, line, '\n')) {
        // ********************************************************************
        // seach the line for a # sign It deniotes all following chars are a
        // comment. Delete the # and all followin chars from the string.
        // *********************************************************************
        std::string::size_type idx=line.find("#");
        if(idx != std::string::npos){           //ok there is a comment in the line
            if(idx != 0){
                line = line.substr(0, idx - 1);
            }
            else{
                continue;               //Skip the line, its all comment
            }
        }
        // *********************************************************
        //This isn't real efficent but its clear
        // at least. Duplicates effort, see below, but what the hey!
        // *********************************************************
        std::istringstream iss(line);

        std::vector<std::string> tokens;
        std::string value;
        while(iss >> value){     //This uses white space (see above) to
            tokens.push_back(value);  //seperate values
        }
        int numTokens = tokens.size();
        if (numTokens == 0){         //blank line
            continue;
        }
        actuatorLines.push_back(line);
        numberOfTokens.push_back(numTokens);
    }
    int numActuator = actuatorLines.size(); //No idea what determines how many
    //lines we should have.

    // ********************************************************************
    // Now we can dimension our vectors
    // ********************************************************************

    actuatorDistance.resize(numActuator,0.0);  //Init to 0.0
    actuatorError.resize(numActuator,0.0);
    actuatorMatrix.resize(numLines,actuatorError); //actuatorError is a convient
    //vector to init the matrix
    //with 0.0's
    // *************************************************************
    // Now fill.  As far as I can tell the actuator file we have should have
    // at least numLines+2 entries per line. If it doesn't, fill with what we
    // have but always assume if we have too few that the the last 2 values are
    // for actuator distance and error. If we have too many than use the 2 after
    // number of lines in control.txt
    // **************************************************************

    for(int i=0;i<numActuator;i++){
        std::istringstream iss(actuatorLines.at(i)); //Pick up the line
        double value;
        // *********************************************************************
        // Hom many values do we have on this line
        // *********************************************************************
        int numMax = numLines;
        if(numberOfTokens.at(i)<numLines+2){
            numMax = numberOfTokens.at(i)-2;
        }
        for(int jj = 0; jj < numMax; jj++){
            iss >> value;
            actuatorMatrix.at(jj).at(i) = value;
        }
        iss >> value;
        actuatorDistance.at(i)=value;

        iss >> value;
        actuatorError.at(i)=value;
    }

    // Make sure we have actuatorstr

    if(pars.has_key("actuator")){
        std::string act=pars["actuator"];
        if(act.size()>0){
            int controlOpt=0;
            if(pars.has_key("control")){
                pars.get("control",controlOpt);
            }
            if(controlOpt==0){  //We have a command override and re fill all error
                // values.
                std::istringstream iss(act);
                double value=0.0;
                int ii=0;
                while(iss>>value){
                    actuatorError.at(ii)=value;
                    ii++;
                }
            }
        }
    }
    return;
}
// ***********************************************************************

void InstrumentFiles::writeBodyFile(std::map< int, std::vector< double>* >& bodyMap, std::string opticsFileName, std::string obsID) {
// Write the body and zernike stuff to a file

    std::ofstream ofs(opticsFileName.c_str());

    ofs << "trackingfile tracking_" << obsID << ".pars" << std::endl;
    std::map< int, std::vector< double >* >::iterator bodyPos;
    std::vector < double >*  pBodyDevice;//pointer to vector we will find on heap

    // Iterate through the devices
    for (bodyPos = bodyMap.begin(); bodyPos != bodyMap.end(); ++bodyPos){
        int devIndex = bodyPos->first;
        pBodyDevice = bodyPos->second;
        int numBodyValues=pBodyDevice->size();
        for (int i = 0; i < numBodyValues; i++){
            if (pBodyDevice->at(i) != -5){
                if (i < 6)
                    ofs << "body " << devIndex << " " << i << " "
                        << std::scientific << std::setprecision(8) << pBodyDevice->at(i)
                        << std::endl;
                else if (i>=6)
                    ofs << "izernike " << devIndex << " " << i-6 << " "
                        << std::scientific << std::setprecision(8) << pBodyDevice->at(i)
                        << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

bool InstrumentFiles::getDeviceIndex(std::string deviceStr, int& deviceIndex) {

// Look for deviceStr in our Surface map and if found return index.
// If not found false

    fSurfaceMapPos = fSurfaceMap.find(deviceStr);
    if (fSurfaceMapPos == fSurfaceMap.end()){
        return false;
    } else {
        deviceIndex = fSurfaceMapPos->second;
        return true;
    }
}




void InstrumentFiles::readoutPars(std::string focalPlaneLayoutFileName, std::string segmentationFileName, std::string readoutString, int camConfig) {

// Read in the focalplanelayout.txt file and the segmentation.txt file and
// from those make up and write out the readout_*.pars

    RngSetSeed32_reseed(1000);

    // ****************************************************************
    //Setup and read in the focalplanelayout.txt file
    // ****************************************************************
    std::ifstream fpLayout(focalPlaneLayoutFileName.c_str());
    PhosimParser fpLPars;
    fpLPars.readStream(fpLayout);

    // *****************************************************************
    //Setup and read in the segemtation.txt file
    // *****************************************************************
    std::ifstream segmentation(segmentationFileName.c_str());
    PhosimParser segmentationPars;
    segmentationPars.readSegmentation(segmentation);


    // iterate through sets of chips
    bool useGroup0=false;
    bool useGroup1=false;
    bool useGroup2=false;

    if(camConfig & 1){
        useGroup0=true;
    }
    if(camConfig & 2){
        useGroup1=true;
    }
    if(camConfig & 4){
        useGroup2=true;
    }


    // *****************************************************************
    // Go through the fpLayout map
    // *****************************************************************
    int numKeys=fpLPars.getNumberKeys();
    for (int keyNum=0;keyNum<numKeys;keyNum++){
        std::string chipID;
        bool gotKey=fpLPars.getKey(keyNum,chipID);
        if(!gotKey){
            std::cout << "Error#1 in readoutPars" << std::endl;
            return;
        }
        // **********************************************
        // Look for the Group designation
        // **********************************************

        std::string line=fpLPars[chipID];
        std::string::size_type idx=line.find("Group");
        bool writeOut=false;
        if(idx!=std::string::npos){
            std::string groupID=line.substr(idx+5,1);//gets character at end of Group
            if( (groupID=="0" && useGroup0) || (groupID=="1" && useGroup1) ||
                (groupID=="2" && useGroup2)){
                writeOut=true;
            }
        }

        // *********************************************************************
        // If we are writing out this chip find the key in the segmentation file
        // (its so easy now!)
        // This key will have as data the number of amplifies for this chip that
        // follow in the segmentation.txt file. They have keys like R00_S21_C06
        // *********************************************************************
        if (writeOut){
            bool gotKey = segmentationPars.has_key(chipID);//Better have this key in the
            if (!gotKey){                            //segmentation.txt file!
                std::cout << "Error#2 in readoutPars. Cannot find key:" << chipID
                          << " in segmentation.txt  file" << std::endl;
                return;
            }
            // ********************************************************************
            // We have the key. Note that this is a chip ID of the form R00_S21.
            // there is no amplifier designation.
            // Get the number of amplifiers for this chip. We make a seperate
            // output file for each chip
            // ********************************************************************
            std::string outputChipFileName= readoutString + chipID + ".pars";
            std::ofstream outChipStream(outputChipFileName.c_str());

            std::vector<std::string> amplifiers;
            segmentationPars.getNameList(chipID, amplifiers);
            for (size_t j(0); j<amplifiers.size()-1; j++) {
                std::istringstream iss(amplifiers[j+1]);
                std::string ampName;
                iss>>ampName;
                std::vector< double > tokens;
                double value;
                while(iss>>value){
                    tokens.push_back(value);
                }
                int serialread   = tokens.at(4);
                int parallelread = tokens.at(5);

                double mean1 = tokens.at(6);
                double var1  = tokens.at(7);
                double gain=mean1*(1+var1*random_gaussian_reseed()/100);

                double mean2 = tokens.at(8);
                double var2  = tokens.at(9);
                double bias=mean2*(1+var2*random_gaussian_reseed()/100);

                double mean3 = tokens.at(10);
                double var3  = tokens.at(11);
                double readnoise=mean3*(1+var3*random_gaussian_reseed()/100);

                double mean4 = tokens.at(12);
                double var4  = tokens.at(13);
                double darkcurrent=mean4*(1+var4*random_gaussian_reseed()/100);

                int parallelPrescan   = tokens.at(14);
                int serialOverscan  = tokens.at(15);
                int serialPrescan   = tokens.at(16);
                int parallelOverscan  = tokens.at(17);
                double hotpixel  = tokens.at(18);
                double hotcolumn = tokens.at(19);

                // Write it all out
                outChipStream << "serialread    " << j << " " << serialread << std::endl;
                outChipStream << "parallelread  " << j << " " << parallelread << std::endl;
                outChipStream << "gain          " << j << " " << gain << std::endl;
                outChipStream << "bias          " << j << " " << bias << std::endl;
                outChipStream << "readnoise     " << j << " " << readnoise << std::endl;
                outChipStream << "darkcurrent   " << j << " " << darkcurrent << std::endl;
                outChipStream << "parallelprescan   " << j << " " << parallelPrescan << std::endl;
                outChipStream << "serialoverscan    " << j << " " << serialOverscan << std::endl;
                outChipStream << "serialprescan     " << j << " " << serialPrescan << std::endl;
                outChipStream << "paralleloverscan  " << j << " " << parallelOverscan << std::endl;
                outChipStream << "hotpixelrate  " << j << " " << hotpixel << std::endl;
                outChipStream << "hotcolumnrate " << j << " " << hotcolumn << std::endl;
            }
            outChipStream.close();
        } //end of  writeout if
    }   //end of chip loop
    return;
}
void InstrumentFiles::makeSurfaceMap(std::string opticsFile) {

// Go through optics_0.txt file thats in the specified instrument directory
// and create a map of surface names to surface index. Ignore "none" in
// both key and surface type
// Also get last device and last surface numbers

    readText opticsPars(opticsFile);

    std::string surfaceName;
    std::string surfaceType;
    fLastSurface = -1;
    int dev = 0;
    for (size_t t(0); t < opticsPars.getSize(); t++) {
        std::istringstream iss(opticsPars[t]);
        iss >> surfaceName;
        iss >> surfaceType;

        if (surfaceName != "none" && surfaceType != "none"){
            fSurfaceMapPos = fSurfaceMap.find(surfaceName);
            if (fSurfaceMapPos != fSurfaceMap.end()){
                fSurfaceMap[surfaceName] = dev;
                dev++;
            }
        }

        if (surfaceType != "none"){
            fLastSurface++;
        }
    }
    return;

}




void InstrumentFiles::focalPlanePars(std::string focalPlaneLayoutFileName, std::string outChipString, int camConfig)
//  ***************************************************************************
// This function makes the chip_999999_R00_S22_C1.pars type files that have
// all the body and zenike values for the chip. That info comes form the
// focalplanelayout file.  We do watch the camConfig value to decide which
// chip files to mek (cut on Group setting).
// ****************************************************************************
{

    // ************************************************************************
    // Now go through the focalplanelayout.txt file making a chip file for
    // any chip in an acceptable Group
    // ***********************************************************************
    //Setup and read in the focalplanelayout.txt file into a PhosiumParser map.
    // ****************************************************************
    std::ifstream fpLayout(focalPlaneLayoutFileName.c_str());
    PhosimParser fpLPars;
    fpLPars.readStream(fpLayout);
    // **********************************************************************
    // Next we need to iterate thorough the lines (a line per chip)  from the
    // focalplanelayout.txt we just parsed and find those chips that match our
    // camConfig specified Group.
    // **********************************************************************
    // We use the bit settings in camConfig to specify the allowed groups
    // Ie. camConfig= 7 meas groups 0 and 1 and 2.
    // decode camcomfig
    // **********************************************************************
    bool useGroup0=false;
    bool useGroup1=false;
    bool useGroup2=false;

    if (camConfig & 1){
        useGroup0=true;
    }
    if (camConfig & 2){
        useGroup1=true;
    }
    if (camConfig & 4){
        useGroup2=true;
    }

    // *****************************************************************
    // Go through the fpLayout map
    // *****************************************************************
    int numKeys = fpLPars.getNumberKeys();
    for (int keyNum = 0; keyNum < numKeys; keyNum++){
        std::string chipID;
        bool gotKey = fpLPars.getKey(keyNum,chipID);
        if (!gotKey){
            std::cout << "Error#1 in focalPlanePars" << std::endl;
            return;
        }
        // **********************************************
        // Look for the Group designation
        // **********************************************

        std::string line = fpLPars[chipID];
        std::string::size_type idx = line.find("Group");
        bool writeOut = false;
        if (idx != std::string::npos) {
            std::string groupID = line.substr(idx+5,1);
            //gets character at end of Group
            if( (groupID == "0" && useGroup0) || (groupID == "1" && useGroup1) ||
                (groupID == "2" && useGroup2)) {
                writeOut=true;
            }
        }

        if (writeOut) {
            std::string outChipFileName = outChipString + chipID + ".pars";
            std::ofstream outChipFile(outChipFileName.c_str());

            // ********************************************************************
            // Write out the body commands for this chip
            // ********************************************************************
            line = line.substr(idx + 6);
            std::istringstream iss(line);

            //Body values first
            for (int i = 0; i < 3; i++) {
                double bodyValue;
                iss >> bodyValue;
                outChipFile << "body " << (fLastSurface+1)  << " " << i << " "
                            << std::fixed << std::setprecision(7) << bodyValue*M_PI/180.0 << std::endl;
            }
            for (int i = 0; i < 3; i++) {
                double bodyValue;
                iss >> bodyValue;
                outChipFile << "body " << (fLastSurface+1)  << " " <<  i+3  << " "
                            << std::fixed << std::setprecision(7) << bodyValue << std::endl;
            }

            std::string pertType;
            iss >> pertType;

            //zernike values next
            if (pertType == "zern") {
                for (int i = 0; i < NZERN; i++) {
                    double zernikeValue;
                    iss >> zernikeValue;
                    outChipFile << "izernike " << fLastSurface << " " << i << " "
                                << std::scientific << std::setprecision(6) << zernikeValue/1000.
                                << std::endl;
                }
            } else if (pertType == "chebyshev") {
                for (int i = 0; i < NCHEB; i++) {
                    double chebyshevValue;
                    iss >> chebyshevValue;
                    outChipFile << "ichebyshev " << fLastSurface << " " << i << " "
                                << std::scientific << std::setprecision(6) << chebyshevValue/1000.
                                << std::endl;
                }
            } else {
                std::cout << "Error: Unknown perturbation type " << pertType << std::endl;
                return;
            }

            // QE variation last
            double QEVar;
            iss >> QEVar;
            outChipFile << "qevariation " << std::fixed << std::setprecision(6) << QEVar
                        << std::endl;
        } //Group test
    } //chip loop
    return;
}
