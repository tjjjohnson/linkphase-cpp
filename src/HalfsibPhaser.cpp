#include "AnimalInfo.cpp"

enum ToPhase {
    no = 0,
    yes = 1,
    alreadyPhased = 2
};


class MarkerPhasingInfo {
public:
    ToPhase toPhase = ToPhase::yes;
    int haplotypeCounts[3] = { 0 };
    int phase1 = 0;
    int phase2 = 0;

    // bool phased = false; use toPhase instead
};

class HalfsibPhaser {
public:
    AnimalInfo *parent;
    Parameters parameters;

    int homozygousCount = 0;
    int heterozygousCount = 0;

    vector<MarkerPhasingInfo> markerPhasingInfoVec;

    HalfsibPhaser(Parameters someParameters, AnimalInfo *aParent) {
        parent = aParent;
        parameters = someParameters;

        for(int i=0; i<parent->gen1.size(); i++) {
            MarkerPhasingInfo markerPhasingInfo;
            markerPhasingInfoVec.push_back(markerPhasingInfo);
        }
    }

    void identifyMarkersToPhase() {
        // was called listtophase
        for(int markerIndex=0; markerIndex < parent->gen1.size(); markerIndex++) {
            MarkerPhasingInfo *markerInfo = &markerPhasingInfoVec[markerIndex];

            // TODO: if male and outside PAR region set tophase=0 and mark as phased and set phase1 and phase2

            for(int offspring=0; offspring < parent->offspring.size(); offspring++) {
                AnimalInfo *offspringInfo = parent->offspring[offspring];
                if( offspringInfo->hap[parent->sex][markerIndex] != 0) markerPhasingInfoVec[markerIndex].haplotypeCounts[0]++;
                if( offspringInfo->hap[parent->sex][markerIndex] == 1) markerPhasingInfoVec[markerIndex].haplotypeCounts[1]++;
                if( offspringInfo->hap[parent->sex][markerIndex] == 2) markerPhasingInfoVec[markerIndex].haplotypeCounts[2]++;
            }
            if(parent->gen1[markerIndex] == parent->gen2[markerIndex]) {
                markerInfo->toPhase = ToPhase::no;
                if(parent->gen1[markerIndex] != missingGenotype) homozygousCount++;
            } else {
                heterozygousCount++;
                if( parent->hap[0][markerIndex] != missingGenotype) {
                    markerInfo->toPhase = ToPhase::alreadyPhased;
                    markerInfo->phase1 = parent->hap[0][markerIndex];
                    markerInfo->phase2 = parent->hap[1][markerIndex];
                    continue;
                }
            }
            if(parent->id==18) {
                //cerr << "Animal 18 offspring count = " << parent->offspring.size() << endl;
                cerr << "Haplotype counts\t" << parent->id << "\t" << markerIndex << "\t" << markerPhasingInfoVec[markerIndex].haplotypeCounts[0]
                     << "\t" << markerPhasingInfoVec[markerIndex].haplotypeCounts[1] << "\t" << markerPhasingInfoVec[markerIndex].haplotypeCounts[2] << endl;
            }
        }

    }

    void run() {
        if(parent->offspring.size() > 0) {
            //allocate nrec(noffspring)
            //nrec = 0
            //animalInfo.se
            //allocate(hap2(noffspring, nmarq)
            vector<int> offSpringHaplotypes;
            identifyMarkersToPhase();

            if (parameters.halfsibPhasing && (parent->offspring.size() > 2 || parent->prePhased())) {
                // initialize hmm with homozygous markers and offspring

                for (int offspring = 0; offspring < parent->offspring.size(); offspring++) {
                    // phase1 and phase2 used in multiple function pass around or create class?
                    AnimalInfo *offspringInfo = parent->offspring[offspring];
                    if (offspring > parameters.nTemplates - 1)
                        break;
                    for (int markerIndex = 0; markerIndex < parent->gen1.size(); markerIndex++) {
                        if (offspringInfo->prephaseInfo[markerIndex] == Prephase::none
                            || offspringInfo->prephaseInfo[markerIndex] > Prephase::both_informative)
                            continue; //use only mendelian marker
                        if (parent->sex == Sex::male && parent->prephaseInfo[markerIndex] == Prephase::sire_hom)
                            continue;
                        if (parent->sex == Sex::female && parent->prephaseInfo[markerIndex] == Prephase::dam_hom)
                            continue;

                        if (parent->isHet(markerIndex) and
                            offspringInfo->hap[parent->sex][markerIndex] != missingGenotype) {
                            // set phase1 and phase2 which are used in multiple functions

                        }
                    }

                    //cerr << "offspring hap = " << animalInfo->offspring[offspring]->hap[animalInfo->sex][5] << endl;
                    //offSpringHaplotypes.push_back()
                }
            }
        }
    }
};