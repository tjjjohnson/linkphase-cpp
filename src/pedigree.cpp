#include <map>
#include <vector> 
using namespace std;

struct PEDIGREE {
    int animal_key;
    int sire_key;
    int dam_key;
    int animal_code;
    int sire_code;
    int dam_code;
    bool processed;
};

void pedigreeTree(std::map<int, PEDIGREE> &pedMap, int indvID, vector<PEDIGREE> &orderedPedigree) {
    if(indvID == 0) {
        return;
    }
    PEDIGREE pedigree = pedMap[indvID];

    if(pedigree.processed) {
        
        return;
    }
    
    pedigree.processed = true;
    pedMap[indvID] = pedigree;
    pedigreeTree(pedMap, pedigree.sire_key, orderedPedigree);
    pedigreeTree(pedMap, pedigree.dam_key, orderedPedigree);
    
    pedigree.animal_code = orderedPedigree.size() + 1;
    pedigree.sire_code = pedMap[pedigree.sire_key].animal_code;
    pedigree.dam_code = pedMap[pedigree.dam_key].animal_code;
    pedigree.processed = true;
    pedMap[indvID] = pedigree;
    orderedPedigree.push_back(pedigree);
    //printf("%d\t%d\t%d\n", pedigree.animal_key, pedigree.sire_key, pedigree.dam_key);
}

vector<PEDIGREE> read_and_sort_pedigree(const char *pedFileName)
{
    int     ak,sire_ak,dam_ak;
    FILE    *pedfile;
    std::map<int, PEDIGREE> pedMap;

    if((pedfile=fopen(pedFileName,"r"))==NULL){
        fprintf(stderr,"Cannot fopen file %s\n", pedFileName);
        exit(0);
    }
    fprintf(stderr, "Reading pedigree from %s\n", pedFileName);
    while (fscanf(pedfile,"%d%d%d",&ak,&sire_ak,&dam_ak)!=EOF){
        PEDIGREE pedigree = {};
        pedigree.animal_key=ak;
        pedigree.sire_key=sire_ak;
        pedigree.dam_key=dam_ak;
        pedigree.processed = false;
        pedMap[pedigree.animal_key] = pedigree;
        while (fgetc(pedfile)!='\n');
    }
    fclose(pedfile);

    vector<PEDIGREE> orderedPedigree;
    orderedPedigree.reserve(pedMap.size());
    fprintf(stderr,"Pedmap size = %ld\n", pedMap.size());
    pedigreeTree(pedMap, 23934418, orderedPedigree);

    for ( const auto &myPair : pedMap ) {
        pedigreeTree(pedMap, myPair.first, orderedPedigree);
    }
    return orderedPedigree;
}