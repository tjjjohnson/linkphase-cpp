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


enum Genotype {
    missing, ref, alt
};
enum Prephase {
    none = 0,
    hom = 1,
    sire_hom = 2,
    dam_hom = 3,
    both_informative = 4
};

class AnimalInfo {
public:
    int id;
    AnimalInfo *sire;
    AnimalInfo *dam;
    bool haplotyped;
    bool genotyped;

    vector <Genotype> hap1;
    vector <Genotype> hap2;
    vector <Genotype> gen1;
    vector <Genotype> gen2;
    vector<int> prephaseInfo;

    AnimalInfo() {
        sire = NULL;
        dam = NULL;
    }

    Genotype sireGen1(int markerIndex) {
        if (sire != NULL)
            return sire->gen1[markerIndex];
        else
            return Genotype::missing;
    }

    Genotype sireGen2(int markerIndex) {
        if (sire != NULL)
            return sire->gen2[markerIndex];
        else
            return Genotype::missing;
    }

    Genotype damGen1(int markerIndex) {
        if (dam != NULL)
            return dam->gen1[markerIndex];
        else
            return Genotype::missing;
    }

    Genotype damGen2(int markerIndex) {
        if (dam != NULL)
            return dam->gen2[markerIndex];
        else
            return Genotype::missing;
    }

    void phaseMendelian() {
        if (!genotyped)
            return;
        //cerr << "here";
        if (sire == NULL && dam == NULL)
            return;
        for (int markerIndex = 0; markerIndex < gen1.size(); markerIndex++) {
            if (gen1[markerIndex] == Genotype::missing
                || hap1[markerIndex] != Genotype::missing)
                return;

            Genotype sireAllele1 = sireGen1(markerIndex);
            Genotype sireAllele2 = sireGen2(markerIndex);
            Genotype damAllele1 = damGen1(markerIndex);
            Genotype damAllele2 = damGen2(markerIndex);
            Genotype allele1 = gen1[markerIndex];
            Genotype allele2 = gen2[markerIndex];

            if (sireAllele1 != Genotype::missing && damAllele1 != Genotype::missing) {
                if (sireAllele1 == sireAllele2 && damAllele1 == damAllele2) {
                    // poth parents homozygous
                    phaseBothParentsHomozygous(markerIndex, sireAllele1, damAllele1, allele1, allele2);
                    continue; // incompatibility - but offspring should be prephased as homozygous (in case of SNPs)
                }
            } else if (sireAllele1 == sireAllele2 && damAllele1 != damAllele2) {
                // sire homozygous
                if (allele1 == sireAllele1) {
                    hap1[markerIndex] = allele1;
                    hap2[markerIndex] = allele2;
                    haplotyped = true;
                    prephaseInfo[markerIndex] = Prephase::sire_hom;
                } else if (allele2 == sireAllele1) {
                    hap1[markerIndex] = allele2;
                    hap2[markerIndex] = allele1;
                    haplotyped = true;
                    prephaseInfo[markerIndex] = Prephase::sire_hom;
                }
            } else if (damAllele1 == damAllele2) {
                // dam homozygous
                if (allele2 == damAllele1) {
                    hap1[markerIndex] = allele1;
                    hap2[markerIndex] = allele2;
                    haplotyped = true;
                    prephaseInfo[markerIndex] = Prephase::dam_hom;
                } else if (allele1 == damAllele1) {
                    hap1[markerIndex] = allele2;
                    hap2[markerIndex] = allele1;
                    haplotyped = true;
                    prephaseInfo[markerIndex] = Prephase::dam_hom;
                }
            }
        }


    }

    void  phaseBothParentsHomozygous(int markerIndex, const Genotype &sireAllele1, const Genotype &damAllele1,
                               Genotype &allele1,
                               Genotype &allele2) {
        if (allele1 == sireAllele1 && allele2 == damAllele1) {
            hap1[markerIndex] = allele1;
            hap2[markerIndex] = allele2;
            haplotyped = true;
            prephaseInfo[markerIndex] = both_informative;
            cerr << "MendelPhased - Parents hom 12" << endl;
        } else if (allele2 == sireAllele1 && allele1 == damAllele1) {
            hap1[markerIndex] = allele2;
            hap2[markerIndex] = allele1;
            haplotyped = true;
            prephaseInfo[markerIndex] = both_informative;
            cerr << "MendelPhased - Parents hom 21" << endl;
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