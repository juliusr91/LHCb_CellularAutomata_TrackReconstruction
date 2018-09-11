//
// Created by julius on 3/5/18.
//


#include "../include/computation.h"
#include "../include/Tools.h"


//calculate the compatability of two hits in the dx/dz dimension
bool areCompatibleInX(float &x1, float &z1, float &x2,float &z2){
    float hitDistance = std::abs(z1 - z2);
    float dxmax = MAX_SLOPE_X * hitDistance;
    return (std::abs(x1 - x2) < dxmax);
}

//calculate the compatability of two hits in the dy/dz dimension
bool areCompatibleInY(float &y1, float &z1,float &y2, float &z2){
    float hitDistance = std::abs(z1 - z2);
    float dxmax = MAX_SLOPE_Y * hitDistance;
    return (std::abs(y1 - y2) < dxmax);
}

//calculate the compatability of two hits in the dx/dz and dy/dz dimension
bool areCompatible(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2){
    bool result = areCompatibleInX(x1,z1,x2,z2) && areCompatibleInY(y1,z1,y2,z2);
    return result;
}

//making all the doublets between the different hits
void makeDoublets (const EventInfo & eventInfo, std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS> & ListDoublets){
    //function makes all possible doublets.
    int numHitsCurrentModule = 0;
    int numHitsNeighbourModule = 0;
    int neighbourHitIndex = 0;
    int currentHitIndex = 0;
    int doubletCounter = 0;

    // -2 in the loop due to structure of detector and z location of sensors.
    // loops over all modules
    for (int module = 0; module < NUMBER_OF_SENSORS - NEIGHBOUR_INDEX; module++) {
        numHitsCurrentModule = eventInfo.module_hitNums[module];
        numHitsNeighbourModule = eventInfo.module_hitNums[module + NEIGHBOUR_INDEX];

        //loop over all hits in the right neighbours
        doubletCounter = 0;
        for (int neighbourHit = 0; neighbourHit < numHitsNeighbourModule; neighbourHit++){
            neighbourHitIndex = eventInfo.module_hitStarts[module + NEIGHBOUR_INDEX] + neighbourHit;

            //loop over all hits in the current module
            for (int hit = 0; hit < numHitsCurrentModule; hit++) {
                currentHitIndex = eventInfo.module_hitStarts[module] + hit;

                //checks compatability of the two hits and makes them into a doublet
                if (std::abs(eventInfo.phi[currentHitIndex] - eventInfo.phi[neighbourHitIndex]) < MAX_RAD_DIFF) {
                    if (areCompatible(eventInfo.hit_Xs[currentHitIndex], eventInfo.hit_Ys[currentHitIndex],
                                      eventInfo.hit_Zs[currentHitIndex],
                                      eventInfo.hit_Xs[neighbourHitIndex], eventInfo.hit_Ys[neighbourHitIndex],
                                      eventInfo.hit_Zs[neighbourHitIndex])) {
//                        ListDoublets[module][doubletCounter] = Doublet(currentHitIndex, neighbourHitIndex, eventInfo.phi[currentHitIndex], eventInfo.phi[neighbourHitIndex]);
                        ListDoublets[module].emplace_back(Doublet(currentHitIndex, neighbourHitIndex, eventInfo.phi[currentHitIndex], eventInfo.phi[neighbourHitIndex]));

                        doubletCounter++;
                    }
                }
            }
        }
    }
}

//calculate the Chi2 between three hits.
//currently used for testing how straight the line between three pointws is
//TODO: replace this by using the 2*Area^2 calculation
float calculateChi2BasedHitID(const EventInfo & eventInfo, const int & hit0, const int & hit1, const int & hit2){
    float x0 = eventInfo.hit_Xs[hit0];
    float y0 = eventInfo.hit_Ys[hit0];
    float z0 = eventInfo.hit_Zs[hit0];

    float x1 = eventInfo.hit_Xs[hit1];
    float y1 = eventInfo.hit_Ys[hit1];
    float z1 = eventInfo.hit_Zs[hit1];

    float x2 = eventInfo.hit_Xs[hit2];
    float y2 = eventInfo.hit_Ys[hit2];
    float z2 = eventInfo.hit_Zs[hit2];

    float td= 1.0f/(z1 -z0);
    float txn = x1 - x0;
    float tyn = y1 - y0;
    float tx = txn * td;
    float ty = tyn * td;

    float dz= z2 - z0;
    float x_prediction = x0 + tx * dz;
    float dx = std::abs(x_prediction - x2);

    float y_prediction = y0 + ty * dz;
    float dy = std::abs(y_prediction - y2);

    float scatterNum = (dx * dx) + (dy * dy);
    float scatterDenom = 1.0f / (z2 - z1);
    float scatter = scatterNum * scatterDenom * scatterDenom;

    return (scatter);
}

void findNeighbours (const EventInfo & eventInfo,  std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS> & ListDoublets){
    //finds all the suitable neighbours of a given doublet
    //first three loops to loop overall sensors, hits and their doublets
    std::vector<Doublet>::iterator index, index2;
    int doubletNumLeft, numDoubletsNeighbour;
    ulong numDoublets;

    for (int sensor = NEIGHBOUR_INDEX; sensor < NUMBER_LAYERS_DOUBLETS; sensor++){
        numDoublets =ListDoublets[sensor].size();



        //loop over hits in sensor
        for (int doubletNum = 0; doubletNum < numDoublets; doubletNum++){

            //find where the right neighbouring hits start based on the hit index
            index = std::lower_bound(ListDoublets[sensor - NEIGHBOUR_INDEX].begin(),ListDoublets[sensor - NEIGHBOUR_INDEX].end(),ListDoublets[sensor][doubletNum].hit1,search_comparator());
            //stop if none are found
            if (index == ListDoublets[sensor - NEIGHBOUR_INDEX].end()){
                break;
            }

            doubletNumLeft = int(index - ListDoublets[sensor - NEIGHBOUR_INDEX].begin());
            //decide on how far to search; either to the end of the vector or MAXNUM_NEIGHBOURS_POINT further
            numDoubletsNeighbour = std::min(int(ListDoublets[sensor - NEIGHBOUR_INDEX].size()), doubletNumLeft + MAXNUM_NEIGHBOURS_POINT);

            //loop over selected neighbours
            for (doubletNumLeft; doubletNumLeft < numDoubletsNeighbour; doubletNumLeft++) {
                float phi_diff = std::abs(ListDoublets[sensor][doubletNum].phi2 - ListDoublets[sensor -NEIGHBOUR_INDEX][doubletNumLeft].phi1);
                if (phi_diff < MAX_RAD_DIFF) {
                    Point<float> a{eventInfo.hit_Xs[ListDoublets[sensor - NEIGHBOUR_INDEX][doubletNumLeft].hit1],
                                   eventInfo.hit_Ys[ListDoublets[sensor - NEIGHBOUR_INDEX][doubletNumLeft].hit1],
                                   eventInfo.hit_Zs[ListDoublets[sensor - NEIGHBOUR_INDEX][doubletNumLeft].hit1]};
                    Point<float> b{eventInfo.hit_Xs[ListDoublets[sensor - NEIGHBOUR_INDEX][doubletNumLeft].hit2],
                                   eventInfo.hit_Ys[ListDoublets[sensor - NEIGHBOUR_INDEX][doubletNumLeft].hit2],
                                   eventInfo.hit_Zs[ListDoublets[sensor - NEIGHBOUR_INDEX][doubletNumLeft].hit2]};


                    Point<float> c{eventInfo.hit_Xs[ListDoublets[sensor][doubletNum].hit2],
                                   eventInfo.hit_Ys[ListDoublets[sensor][doubletNum].hit2],
                                   eventInfo.hit_Zs[ListDoublets[sensor][doubletNum].hit2]};

                    if (two_area_sq_v(a, b, c) < 10) {

                        Neighbour newNeighbour{(sensor - NEIGHBOUR_INDEX), doubletNumLeft};
                        ListDoublets[sensor][doubletNum].makeNewNeighbour(newNeighbour);
                    }
                }
            }
        }
    }
}

//function of the Cellular Automaton evolution
int checkNeighbours(Doublet & currentDoublet, std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS>  & ListDoublets){
    int neighbourModuleNum = 0;
    int neighbourDoubletNum = 0;

    //loops over own neighbours
    for (int numNeighbour = 0; numNeighbour < MAXNUM_NEIGHBOURS_DOUBLET; numNeighbour++){
        neighbourModuleNum = currentDoublet.neighbours[numNeighbour].module;
        neighbourDoubletNum = currentDoublet.neighbours[numNeighbour].doubletNum;
        if (neighbourModuleNum != -1) {
            if (currentDoublet.status == ListDoublets[neighbourModuleNum][neighbourDoubletNum].status) {
                currentDoublet.newStatus++;
                return 1;
            }
        } else{
            return 0;
        }
    }
    return 0;
}

//starts the evolution
void runCA (std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS>  & ListDoublets){
    //runs the actual CA and updates the status
    int numChanges = 100000;

    //continues running until no more updates
    while (numChanges!=0){
        numChanges = 0;
        for (int sensor = NEIGHBOUR_INDEX; sensor < NUMBER_LAYERS_DOUBLETS; sensor++) {
            ulong numDoublets =ListDoublets[sensor].size();
            for (int doubletNum = 0; doubletNum < numDoublets; doubletNum++) {
                    numChanges += checkNeighbours(ListDoublets[sensor][doubletNum], ListDoublets);

            }
        }

        //updates the status variable from the new status
        for (int sensor = NEIGHBOUR_INDEX; sensor < NUMBER_LAYERS_DOUBLETS; sensor++) {
            ulong numDoublets =ListDoublets[sensor].size();
            for (int doubletNum = 0; doubletNum < numDoublets; doubletNum++) {
                    ListDoublets[sensor][doubletNum].status = ListDoublets[sensor][doubletNum].newStatus;
            }
        }
    }
}

//fourth step of extracting the most feasible tracks
//recursively calls itself to extract more and more tracks
//extracts one track per starting doublet depending on the chi2
void extractMostFeasibleTrack(Track & startTrack, std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS> & ListDoublets, Doublet & doublet, EventInfo & eventInfo){
    float chi2 = 0;
    float oldChi2 = 10000;
    ulong bestNeighbour;
    Doublet bestDoublet;

    int lengthNeighbours = 0;

    for (ulong neighbourNum = 0; neighbourNum < MAXNUM_NEIGHBOURS_DOUBLET; neighbourNum++){
        if (doublet.neighbours[neighbourNum].module!=-1) {
            Doublet neighbourDoublet = ListDoublets[doublet.neighbours[neighbourNum].module][doublet.neighbours[neighbourNum].doubletNum];
            if (neighbourDoublet.status == (doublet.status - 1)) {
                startTrack.length++;
                int index1 = neighbourDoublet.hit1; //index

                chi2 = calculateChi2BasedHitID(eventInfo, index1, startTrack.listHits.end()[-1],
                                               startTrack.listHits.end()[-2]);
                if (chi2 < oldChi2) {
                    oldChi2 = chi2;
                    bestNeighbour = index1;
                    bestDoublet = neighbourDoublet;
                    lengthNeighbours++;
                }
            }
        }else{
            break;
        }
    }

    if (lengthNeighbours>0){
        startTrack.listHits.emplace_back(bestNeighbour);
        extractMostFeasibleTrack(startTrack, ListDoublets, bestDoublet, eventInfo);
    }
}

//makes the start of a track and then extracts the straightest track.
void extractTracks(std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS>  & ListDoublets, EventInfo & eventInfo, std::vector<Track> & AllCollectedTracks){
    for (ulong sensor = (NUMBER_LAYERS_DOUBLETS - 1); sensor >= NEIGHBOUR_INDEX ; sensor--) {
        ulong numDoublets =ListDoublets[sensor].size();
        for (int doubletNum = 0; doubletNum < numDoublets; doubletNum++) {
            if (ListDoublets[sensor][doubletNum].status > 1) {
                //make tracks
                Track startTrack{ListDoublets[sensor][doubletNum].hit2, ListDoublets[sensor][doubletNum].hit1};
                startTrack.length = 2;

                extractMostFeasibleTrack(startTrack, ListDoublets, ListDoublets[sensor][doubletNum], eventInfo);
                AllCollectedTracks.emplace_back(startTrack);

            }

        }
    }
}

//remove tracks that are shorter than MIN_LENGTH_TRACK
void removeShortTracks (std::vector<Track> & AllCollectedTracks){
    int numTracks = AllCollectedTracks.size();
    for (int trackNum = 0; trackNum < numTracks; trackNum++){
        if (AllCollectedTracks[trackNum].length < MIN_LENGTH_TRACK){ //min size of track
            AllCollectedTracks.erase(AllCollectedTracks.begin() + trackNum);
            trackNum --;
        }
    }
}

//reduces the number of clone and ghost tracks (false positive)
void removeClonesAndGhosts (std::vector<Track> & AllCollectedTracks, std::array<Track, MAXNUM_TRACKS> & longTracks){
    //sorts on length and then Chi2
    std::sort(AllCollectedTracks.begin(), AllCollectedTracks.end(), sortingOnLengthAndChi2());

    int trackIndex = 0;
    std::set<int> usedHits;
    int numTracks = AllCollectedTracks.size();

    //only take into account MAXNUM_TRACKS number of tracks
    numTracks = std::min(numTracks, MAXNUM_TRACKS);

    for (int trackNum = 0; trackNum < numTracks; trackNum++){
        int counter = 0;

        int numHitsInTrack = AllCollectedTracks[trackNum].length;
        for (int hitNum = 0; hitNum < numHitsInTrack; hitNum++) {
            bool contains_Hit =
                    usedHits.find(AllCollectedTracks[trackNum].listHits[hitNum]) != usedHits.end();
            if (contains_Hit) {
                counter++;
            }
        }

        //determines how many of the hits are already used in another track
        float percentageOfUsedHits = float(counter)/float(AllCollectedTracks[trackNum].length);
        if (percentageOfUsedHits< MAX_PERCENTAGE_OF_USED_HITS){
            for (int hitNum = 0; hitNum <  numHitsInTrack; hitNum++) {
                usedHits.emplace(AllCollectedTracks[trackNum].listHits[hitNum]);
            }
            longTracks[trackIndex] = AllCollectedTracks[trackNum];
            trackIndex++;
        }
    }
}

//function to write all the tracks to a binary file.
//same format as needed for the evaluation by Daniel
//my version can be found here:
//https://github.com/juliusr91/velopix_tracking

void writeToBinary (std::array<Track, MAXNUM_TRACKS> longTracks, int eventNumber, EventInfo & eventInfo){

    int32_t numberTracks = 0;

    for (int i = 0; i < MAXNUM_TRACKS; i++){
        if (longTracks[i].length != 0){
            numberTracks++;
        }
    }

    std::ofstream outfile (std::string("../results3/") + std::to_string(eventNumber) + std::string(".bin"), std::ios::binary);
    outfile.write((char*) &numberTracks, sizeof(int32_t));

    for (int trackNum = 0; trackNum < numberTracks; trackNum++){

        const ulong numberOfHits = longTracks[trackNum].listHits.size();
        outfile.write((char*) &numberOfHits, sizeof(uint32_t));

        for (int hitNum = 0; hitNum < numberOfHits; hitNum++){
            const uint32_t hitID = eventInfo.hit_IDs[longTracks[trackNum].listHits[hitNum]];
            outfile.write((char*) &hitID, sizeof(uint32_t));
        }
    }
    outfile.close();
}

//main function that runs all the different steps
void runCellularAutomaton (std::vector<std::vector<uint8_t>> & input){
    int numberOfEvents = input.size();

    int eventNumber;
    for (eventNumber = 0; eventNumber < numberOfEvents; eventNumber++) {

        auto eventInfo{EventInfo(input[eventNumber])};

        //1. Make Doublets
        std::array<std::vector<Doublet>, NUMBER_LAYERS_DOUBLETS> ListDoublets;
        makeDoublets(eventInfo, ListDoublets);


        //2. Find Neighbours
        findNeighbours(eventInfo, ListDoublets);

//        3. Run Cellular Automaton
        runCA(ListDoublets);

//        4. Reconstruct all tracks
        std::vector<Track> AllCollectedTracks;
        extractTracks(ListDoublets, eventInfo, AllCollectedTracks);


//        //5. remove short tracks
        removeShortTracks(AllCollectedTracks);


//        std::cout << "im here "<< eventNumber <<std::endl;
//        6. remove ghost and clone tracks
        std::array<Track, MAXNUM_TRACKS> longTracks;
        removeClonesAndGhosts(AllCollectedTracks, longTracks);


//        //do binary
        if (PRINT_BINARY){
            writeToBinary(longTracks, eventNumber, eventInfo);
        }
    }
}
