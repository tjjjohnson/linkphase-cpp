//
// Created by thjoh0 on 7/1/22.
//

#pragma once

#include <algorithm>
#include <string>
#include <string.h>
#include <vector>
#include <iostream>
#include <map>
using namespace std;

const int missingGenotype = 0;

struct Parameters {
    string pedFile;
    string genotypeFile;
    string markerFile;
    bool halfsibPhasing; // if(paramline=='yes')phase_option=2
    bool hmmPhasing;
    bool checkPrephasing; //if(paramline=='yes') map_option=1
    // if(paramline=='yes' .and. phase_option==2)phase_option=3
    // if(paramline=='yes' .and. phase_option==1)stop 'Error in parameter file, HMM_PHASING must be preceded by HALFSIB_PHASING'
    bool sexChrom; //autosome=0
    bool sexMap;
    bool columns;
    int nTemplates;
    int nIterations;
    double gerr;

};

enum Prephase {
    none = 0,
    hom = 1,
    sire_hom = 2,
    dam_hom = 3,
    both_informative = 4
};

enum Sex {
    unknown = -1,
    male = 0,
    female = 1
};

class AnimalInfo {
public:
    int id;
    AnimalInfo *sire;
    AnimalInfo *dam;
    vector<AnimalInfo *> offspring;
    bool haplotyped;
    bool genotyped;
    Sex sex;

    vector <int> hap[2];
    vector <int> gen1;
    vector <int> gen2;
    vector<Prephase> prephaseInfo;

    AnimalInfo() {
        sire = NULL;
        dam = NULL;
        sex = Sex::unknown;
    }

    int sireGen1(int markerIndex) {
        if (sire != NULL)
            return sire->gen1[markerIndex];
        else
            return missingGenotype;
    }

    int sireGen2(int markerIndex) {
        if (sire != NULL)
            return sire->gen2[markerIndex];
        else
            return missingGenotype;
    }

    int damGen1(int markerIndex) {
        if (dam != NULL)
            return dam->gen1[markerIndex];
        else
            return missingGenotype;
    }

    int damGen2(int markerIndex) {
        if (dam != NULL)
            return dam->gen2[markerIndex];
        else
            return missingGenotype;
    }

    bool isHet(int markerIndex) {
        if(gen1[markerIndex] == missingGenotype || gen2[markerIndex] == missingGenotype)
            return false;
        return gen1[markerIndex] != gen2[markerIndex];
    }

    bool prePhased() {
        if (dam != NULL && dam->genotyped)
            return true;
        else if(sire != NULL && sire->genotyped)
            return true;
        else
            return false;
    }

    void phaseMendelian() {
        if (!genotyped)
            return;
        //cerr << "here";
        if (sire == NULL && dam == NULL)
            return;
        for (int markerIndex = 0; markerIndex < gen1.size(); markerIndex++) {
            if (gen1[markerIndex] == missingGenotype
                || hap[0][markerIndex] != missingGenotype)
                continue;

            int sireAllele1 = sireGen1(markerIndex);
            int sireAllele2 = sireGen2(markerIndex);
            int damAllele1 = damGen1(markerIndex);
            int damAllele2 = damGen2(markerIndex);
            int allele1 = gen1[markerIndex];
            int allele2 = gen2[markerIndex];

            if (sireAllele1 != missingGenotype && damAllele1 != missingGenotype) {
                if (sireAllele1 == sireAllele2 && damAllele1 == damAllele2) {
                    // both parents homozygous
                    cerr << "Both parents homozygous marker=" << markerIndex + 1<< "\t" << id << "\t" << sire->id << "\t" << dam->id << endl;
                    if (allele1 == sireAllele1 && allele2 == damAllele1) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                        haplotyped = true;
                        prephaseInfo[markerIndex] = both_informative;
                        //cerr << "MendelPhased - Parents hom 12" << endl;
                    } else if (allele2 == sireAllele1 && allele1 == damAllele1) {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                        haplotyped = true;
                        prephaseInfo[markerIndex] = both_informative;
                        //cerr << "MendelPhased - Parents hom 21" << endl;
                    }
                    continue; // incompatibility - don't prephase
                } else if (sireAllele1 == sireAllele2 && damAllele1 != damAllele2) {
                    // sire homozygous
                    if (allele1 == sireAllele1) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                        haplotyped = true;
                        prephaseInfo[markerIndex] = Prephase::sire_hom;
                        //cerr << "MendelPhased - Sire hom 12" << endl;
                    } else if (allele2 == sireAllele1) {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                        haplotyped = true;
                        prephaseInfo[markerIndex] = Prephase::sire_hom;
                        //cerr << "MendelPhased - Sire hom 21" << endl;
                    }
                    continue; // incompatibility - but offspring should be prephased as homozygous (in case of SNPs)
                } else if (damAllele1 == damAllele2) {
                    // dam homozygous
                    if (allele2 == damAllele1) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                        haplotyped = true;
                        prephaseInfo[markerIndex] = Prephase::dam_hom;
                        //cerr << "MendelPhased - Dam hom 12 marker=" << markerIndex << "\t" << id << "\t" << sire->id << "\t"<< dam->id << endl;
                    } else if (allele1 == damAllele1) {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                        haplotyped = true;
                        prephaseInfo[markerIndex] = Prephase::dam_hom;
                        //cerr << "MendelPhased - Dam hom 21 marker=" << markerIndex << "\t" << id << "\t" << sire->id << "\t" << dam->id << endl;
                    }
                    continue;
                }

                // parents both heterozygous
                if (sireAllele1==damAllele1 && sireAllele2==damAllele2) continue;
                if (sireAllele1==damAllele2 && sireAllele2==damAllele1) continue;
                // will we ever get past here?
                cerr << "Multiallelic handling isn't implemented in the c++ version of linkphase" << endl;
                exit(1);


            } else if (sireAllele1 != missingGenotype && damAllele1 == missingGenotype) {
                if (sireAllele1 == sireAllele2) {
                    prephaseInfo[markerIndex] = Prephase::sire_hom;
                    if (allele1 == sireAllele1) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                    } else {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                    }
                    haplotyped = true;
                    continue;
                } else {
                    // sire heterozygous
                    if(sireAllele1 == allele1 && sireAllele2 == allele2) continue;
                    if(sireAllele1 == allele2 && sireAllele2 == allele1) continue;

                    if(allele1 == sireAllele1 || allele2 == sireAllele2) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                        haplotyped = true;
                        continue;
                    } else if(allele2 == sireAllele1 || allele2 == sireAllele2) {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                        haplotyped = true;
                        continue;
                    }
                }
            } else if (sireAllele1 == missingGenotype && damAllele1 != missingGenotype) {
                // only dam is genotyped
                if (damAllele1 == damAllele2) {
                    prephaseInfo[markerIndex] = Prephase::dam_hom;
                    if (allele1 == damAllele1) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                    } else {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                    }
                    haplotyped = true;
                    continue;
                } else {
                    // sire heterozygous
                    if(damAllele1 == allele1 && damAllele2 == allele2) continue;
                    if(damAllele1 == allele2 && damAllele2 == allele1) continue;

                    if(allele1 == damAllele1 || damAllele2 == sireAllele2) {
                        hap[0][markerIndex] = allele1;
                        hap[1][markerIndex] = allele2;
                        haplotyped = true;
                        continue;
                    } else if(allele2 == damAllele1 || damAllele2 == sireAllele2) {
                        hap[0][markerIndex] = allele2;
                        hap[1][markerIndex] = allele1;
                        haplotyped = true;
                        continue;
                    }
                }
            }
        }


    }

    string toString() {
        string sireString = (sire != NULL) ? to_string(sire->id) : "None";
        string damString = (dam != NULL) ? to_string(dam->id) : "None";
        return to_string(id) + "<-(" + sireString + "," + damString + ")";
    }

    //allocate(hap(nani,2,nmarq),typ(nani,2*nmarq),genotyped(nani))
    //allocate(sire(nani),dam(nani),haplotyped(nani),hapin(nmarq),
};