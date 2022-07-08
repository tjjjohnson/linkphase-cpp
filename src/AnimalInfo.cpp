//
// Created by thjoh0 on 7/1/22.
//

#include <algorithm>
#include <string>
#include <string.h>
#include <vector>
#include <iostream>
#include <map>
using namespace std;

const int missingGenotype = 0;

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
                if (allele1==sireAllele1 && allele2==damAllele1) {
                    hap[1][markerIndex]=allele1;
                    hap[2][markerIndex]=allele2;
                    haplotyped = true;
                    continue;
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