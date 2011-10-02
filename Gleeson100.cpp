//standard header libraries - allow us to call various basic prewritten functions (not all are used in this program, but it's not harmful to leave them in)
#include<iostream>//for reading & writing from the screen
#include<fstream>//for reading & writing to files
#include<math.h>//for performing mathematical operations such as exponents and square roots
#include<sstream>//for converting integers to strings (sets of characters)
#include<vector>//vectors are lists of things - very useful for many things
#include<map>//maps are like vectors, except the indices can be anything, not just sequential numbers (see codon table for more)
#include<stdio.h>//input formatting
#include<stdlib.h>//random useful functions
#include<limits.h>//defines maximum and minimum integers, etc. (useful for sentinels)
#include<assert.h>//if something should be true, but you want to make sure it's true, assert is useful
#include<algorithm>//useful for predefined functions such as list sorting
#include<time.h>//useful for timing things and for random number generators
#include<string.h>//useful for string operations

using namespace std;//ignore for now - just allows one to use easy shortcuts for functions

//string identifiers of all patients
vector<string> patientList;

//base call data structure
//used as part of the SNP data structure to represent genotype/depth/quality triplets
struct baseCall
{
    int genotype;
    int depth;
    double quality;
    string AD;
    string PL;
};

class dbSNP
{
    public:
        int id;
        char refNT;
        char mutNT;
        friend bool operator< (const dbSNP & s1, const dbSNP & s2);
};

bool operator< (const dbSNP & s1, const dbSNP & s2)
{
    return s1.id < s2.id;
}

//SNP data structure
//used for pretty much all the operations involving the Broad SNP data
struct SNP
{
    string chr;
    int genomicLoc, codonPos, codonStart, protPos, totAA, spliceDist;
    char refNT, mutNT, refAA, mutAA;
    vector<baseCall> calls;
    string info, dbSNP, geneName, transcriptID, polyPhen, codon, fullName;
    bool strand, is_splice, is_utr;
    string shares, unsure, full_shares;
    double GERP, phastCons, SIFTProb;
    string SIFT, OMIM, Domains, GTCounts, AlleleCounts;
    string HGMD_SNP, HGMD_Gene;
};

struct indel
{
    string chr;
    int start;
    string ref, alt;
    string ac, dp;
    bool is_unique;
    string shares, het, full_shares;
    bool is_utr;
    bool is_splice;
    string gene;
    string transcript;
    string Domains;
    int ProtStart, ProtEnd;
};

//exon data structure
//used as part of gene data structure (see below) to help determine whether or not exons are part of low coverage regions
struct exon
{
    int start;
	int end;
};

struct domain : public exon
{
    string name;
};

//gene data structure
//used in low coverage region analyses to characterize RefSeq info
struct gene
{
	string name;
	string chr;
	int start;
	int end;
	bool strand;
	int cds_start;
	int cds_end;
    vector<exon> exonList;
};

//see below for description of functions
void parseForHomMap(const vector<SNP> & SNPList);
void parseForSift(const vector<SNP> & SNPList);
void parseForPlink(const vector<SNP> & SNPList);
void parseForGVS(const vector<SNP> & SNPList);
void outputNonSyn(vector<SNP> & SNPList);
void inputData(vector<SNP> & SNPList);
void filterCoverages(vector<SNP> & SNPSamples);
void filterPlink(vector<SNP> & SNPList, const int & sampleID, const char & Plink_level);
void filterPlink(vector<vector<indel> > & indelList);
void filterHomMap(const vector<SNP> & SNPList, vector<vector<SNP> > & SNPSamples);
void outputSNP(ofstream & fout, const vector<SNP> & SNPList, const int & i, const int & j);
void filterGERP(vector<SNP> & SNPList);
string classifyAA(const char & AA);
char ThreeToOne(const string & AA);
void inputPolyphen(vector<SNP> & SNPList);
void inputSIFT(vector<vector<SNP> > & SNPSamples, const string & ID);
void inputSIFT(vector<SNP> & SNPList, const string & ID);
void plinkPipeline(vector<SNP> & SNPList);
void parsePipeline(vector<SNP> & SNPList);
void homMapPipeline(vector<SNP> & SNPList);
void indelPipeline();
void analyzeLowCovData();
void generateLowCovData();
void inputGeneList(vector<gene> & geneList);
void inputIndels(vector<vector<indel> > & indelList);
void filterIndels(vector<vector<indel> > & indelList, const double & freq);
void filterUnique(vector<vector<indel> > & indelList, const double & freq);
void outputIndels(ofstream & fout, vector<vector<indel> > & indelList);
void dbSNPFilter(vector<SNP> & SNPList);
char flipGenotypes(const char & c);
void randomTrial(const vector<SNP> & SNPList, const string & header);
void inputDomains(vector<SNP> & SNPList);
void inputDomains(vector<vector<indel> > & indelList);
void insertAlleleFreqs(vector<SNP> & SNPList, const vector<dbSNP> & rsIDs);
void filterGenic(vector<vector<indel> > & indelList, const vector<gene> & geneList);
void inputHGMD(vector<SNP> & SNPList);
int chr2int(string chr);

//main function
//first thing the program looks for - everything runs from the first line here to the last line here (function defined by everything in {})
int main()
{
    //indelPipeline();
    //generateLowCovData();
    //analyzeLowCovData();
    vector<SNP> SNPList;//list of SNPs
    cout << "Parse or filter pipeline (p/f)?\n";
    char c;
    cin >> c;
    inputData(SNPList);//input list of SNPs
    //randomTrial(SNPList, "WMD");
    if(c == 'p' || c == 'P')
        parsePipeline(SNPList);//parse list of SNPs to formats palatable to various online programs
    else if(c == 'f' || c == 'F')
        plinkPipeline(SNPList);//use Plink as part of filters
    else cout << "Invalid option\n";
    return 0;
}

//sandboxing
void myTest(){
    cout << "Hello" << endl;
}

//function parsePipeline()
//Converts SNP lists into files that are viable inputs for filters downstream
//Input: list of SNPs (SNPList)
//Output: files for homozygosity mapper, plink, and polyphen
void parsePipeline(vector<SNP> & SNPList)
{
    myTest();
    //parseForPlink(SNPList);
    //parseForGVS(SNPList);
    //parseForSift(SNPList);
}

//function indelPipeline()
//Processes lists of indels to find potentially relevant ones
//Input: None (file containing indel list)
//Output: file containing SIFT-parsed indels that passed filters
void indelPipeline()
{
    cout << "Indel pipeline\n";
    vector<vector<indel> > indelList;
    inputIndels(indelList);
    //filterIndels(indelList, 0.75);
    //filterUnique(indelList, 0.75);
    filterUnique(indelList, 0);
    vector<gene> geneList;
    inputGeneList(geneList);
    filterGenic(indelList, geneList);
    //filterPlink(indelList);
    inputDomains(indelList);
    ofstream fout("./Output/IndelList_PlinkRegs.txt");
    fout << "SIFT Input\tUnique\tHeterozygous\tDescriptor\tGene\tProtein Loc\tDomain";
    for(int i = 0; i < (int)patientList.size(); i++)
        fout << "\t" << patientList[i];
    fout << "\n";
    outputIndels(fout, indelList);
}

//function plinkPipeline()
//Filters list of SNPs based on GERP conservation, coverage constraints, and Plink runs of homozygosity
//Sift parser here too, but may be movable
//Then outputs resulting list of SNPs
//Input: list of SNPs (SNPList)
//Output: file with filtered list of SNPs
void plinkPipeline(vector<SNP> & SNPList)
{
    inputSIFT(SNPList, "7038");
    filterGERP(SNPList);
    dbSNPFilter(SNPList);
    inputDomains(SNPList);
    inputHGMD(SNPList);
    //inputPolyphen(SNPList);
    filterCoverages(SNPList);
    //filterPlink(SNPSamples);
}

//function randomTrial()
//Finds all instances in which *ONLY* patients between positions min & max are homozygous mutant
//Not in the pipeline, but just a random aside that interested me
//Input: list of SNPs (SNPList)
//Output: SNP(s) that fit this criteria to screen
void randomTrial(const vector<SNP> & SNPList, const string & header)
{
    vector<int> max_id;
    int max_count = 0;
    int count = 0;
    for(int i = 0; i < (int)SNPList.size(); i++)
    {//for all SNPs
        int j;
        count = 0;
        for(j = 0; j < (int)SNPList[i].calls.size(); j++)
        {//for all individuals
            if(((int)patientList[j].find(header) != -1 && SNPList[i].calls[j].genotype == 2 && SNPList[i].calls[j].depth >= 8))
            {
                count++;
            }
            else if((int)patientList[j].find(header) == -1 && SNPList[i].calls[j].genotype == 2 && SNPList[i].calls[j].depth >= 8)
            {
                count = 0;
                break;
            }
        }
        if(count > max_count)
        {
            max_id.clear();
            max_count = count;
            max_id.push_back(i);
        }
        else if(count == max_count)
            max_id.push_back(i);
    }
    cout << max_count << "\n";
    for(int i = 0; i < (int)max_id.size(); i++)
    {
        cout << SNPList[max_id[i]].chr << "\t" << SNPList[max_id[i]].genomicLoc << "\n";
    }
}

char flipGenotypes(const char & c)
{
    if(c == 'T')
        return 'A';
    if(c == 'G')
        return 'C';
    if(c == 'C')
        return 'G';
    if(c == 'A')
        return 'T';
    return '-';
}

void inputHGMD(vector<SNP> & SNPList)
{
    cout << "Inputting HGMD Annotations\n";
    ifstream fin("./NCBIFiles/HGMD19_dbSNP.txt");
    string line;
    map<string, string> dbSNP_info;
    getline(fin, line);
    getline(fin, line);
    while(!fin.eof())
    {
        dbSNP_info[line.substr(0, line.find('\t'))] = line.substr(line.find('\t') + 1);
        getline(fin, line);
    }
    fin.close();

    fin.open("./NCBIFiles/HGMD19_Genes.txt");
    map<string, string> gene_info;
    getline(fin, line);
    getline(fin, line);
    while(!fin.eof())
    {
        if(gene_info.find(line.substr(0, line.find('\t'))) == gene_info.end())
            gene_info[line.substr(0, line.find('\t'))] = line.substr(line.find('\t') + 1);
        else gene_info[line.substr(0, line.find('\t'))] += ("; " + line.substr(line.find('\t') + 1));
        getline(fin, line);
    }
    fin.close();

    fin.open("./NCBIFiles/HGMD19_SNPs.txt");
    getline(fin, line);
    getline(fin, line);
    string chr;
    int loc;
    chr = line.substr(0, line.find(':'));
    line = line.substr(line.find(':') + 1);
    loc = atoi(line.substr(0, line.find('\t')).c_str());
    line = line.substr(line.find('\t') + 1);
    int i = 0;

    while((i < (int)SNPList.size()) && !fin.eof())
    {
        if((chr2int(chr) > chr2int(SNPList[i].chr)) || (chr == SNPList[i].chr && loc > SNPList[i].genomicLoc))
        {
            if(dbSNP_info.find(SNPList[i].dbSNP) != dbSNP_info.end())
            {
                SNPList[i].HGMD_SNP = dbSNP_info[SNPList[i].dbSNP];
            }
            if(gene_info.find(SNPList[i].geneName) != gene_info.end())
            {
                SNPList[i].HGMD_Gene = gene_info[SNPList[i].geneName];
            }
            i++;
            continue;
        }
        if(chr == SNPList[i].chr && loc == SNPList[i].genomicLoc)
        {
            if(gene_info.find(SNPList[i].geneName) != gene_info.end())
            {
                SNPList[i].HGMD_Gene = gene_info[SNPList[i].geneName];
            }
            SNPList[i].HGMD_SNP = line;
            i++;
            continue;
        }
        getline(fin, line);
        chr = line.substr(0, line.find(':'));
        line = line.substr(line.find(':') + 1);
        loc = atoi(line.substr(0, line.find('\t')).c_str());
        line = line.substr(line.find('\t') + 1);
    }
    while(i < (int)SNPList.size())
    {
        if(dbSNP_info.find(SNPList[i].dbSNP) != dbSNP_info.end())
        {
            SNPList[i].HGMD_SNP = dbSNP_info[SNPList[i].dbSNP];
        }
        if(gene_info.find(SNPList[i].geneName) != gene_info.end())
        {
            SNPList[i].HGMD_Gene = gene_info[SNPList[i].geneName];
        }
        i++;
    }
    fin.close();
}

int chr2int(string chr)
{
    if(chr == "chrX")
        return 23;
    else if(chr == "chrY")
        return 24;
    else return atoi(chr.substr(3).c_str());
}

void dbSNPFilter(vector<SNP> & SNPList)
{
    cout << "Filter dbSNP genotype frequency of 0 (y/n)?\n";
    char choice;
    cin >> choice;
    vector<dbSNP> rsIDs;//list of snp ids to go through
    vector<string> GTFreqs;
    dbSNP dbs;
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(SNPList[i].dbSNP.size() > 2)
        {
            dbs.id = atoi((SNPList[i].dbSNP.substr(2)).c_str());
            dbs.mutNT = SNPList[i].mutNT;
            dbs.refNT = SNPList[i].refNT;
            rsIDs.push_back(dbs);
        }
    }
    GTFreqs.resize(rsIDs.size(), "");
    sort(rsIDs.begin(), rsIDs.end());
    vector<bool> isGood(rsIDs.size(), true);
    map<string, string> codes;
    ifstream fin_code("./NCBIFiles/UniGty.bcp");
    string line;
    getline(fin_code, line);
    string ID;
    while(!fin_code.eof())
    {
        ID = line.substr(0, line.find('\t'));
        line = line.substr(line.find('\t') + 1);
        codes[ID] = line.substr(0, line.find('\t'));
        getline(fin_code, line);
    }

    ifstream fin("./NCBIFiles/SNPGtyFreq.bcp");
    getline(fin, line);
    int count = 0, id;
    vector<string> genotypes;
    string s, sflip;
    bool is_flipped;
    while(!fin.eof() && count < (int)rsIDs.size())
    {
        id = atoi(line.substr(0, line.find('\t')).c_str());//rsID
        if(id > rsIDs[count].id)
        {//if the rsID is too big, compare list of codes from when they were equal
            is_flipped = true;
            if(genotypes.size() > 0)
            {
                s = rsIDs[count].refNT;
                s += "/";
                s += rsIDs[count].refNT;
                for(int i = 0; is_flipped && i < (int)genotypes.size(); i++)
                {
                    if((int)genotypes[i].find(rsIDs[count].refNT) != -1)
                    {//if we see the reference nucleotide in the set of genotypes
                        is_flipped = false;
                    }
                }
                s = rsIDs[count].mutNT;
                s += "/";
                s += rsIDs[count].mutNT;
                if(!is_flipped)
                {
                    for(int i = 0; i < (int)genotypes.size(); i++)
                    {
                        if(genotypes[i] == s)
                        {
                            isGood[count] = false;
                        }
                    }
                }
                else
                {
                    sflip = flipGenotypes(rsIDs[count].mutNT);
                    sflip += "/";
                    sflip += sflip[0];
                    for(int i = 0; i < (int)genotypes.size(); i++)
                    {
                        if(genotypes[i] == s || genotypes[i] == sflip)
                        {
                            isGood[count] = false;
                        }
                    }
                }
                genotypes.resize(0);
            }
            count++;
            continue;
        }
        if(id < rsIDs[count].id)
        {//if the rsID is too small, go on to the next line
            getline(fin, line);
            continue;
        }
        if(id == rsIDs[count].id)
        {//if the rsID fits, add to stack
            line = line.substr(line.find('\t')+1);
            s = line.substr(0, line.find('\t'));
            if(codes.find(s) != codes.end())
            {
                genotypes.push_back(codes[s]);
                line = line.substr(line.find('\t') + 1);
                if(GTFreqs[count] != "")
                    GTFreqs[count] += ",";
                GTFreqs[count] += codes[s] + ":" + line.substr(0, line.find('\t'));
            }
            getline(fin, line);
        }
    }
    vector<SNP> newSNPList;
    bool is_inserted;
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(SNPList[i].dbSNP.size() < 2)
        {
            newSNPList.push_back(SNPList[i]);
        }
        else
        {
            is_inserted = false;
            id = atoi((SNPList[i].dbSNP.substr(2)).c_str());
            for(int j = 0; (!is_inserted) && j < (int)rsIDs.size(); j++)
            {
                if(rsIDs[j].id == id)
                {
                    is_inserted = true;
                    if(isGood[j] || choice == 'n' || choice == 'N')
                    {
                        SNPList[i].GTCounts = GTFreqs[j];
                        newSNPList.push_back(SNPList[i]);
                    }
                }
            }
        }
    }
    SNPList = newSNPList;
    cout << SNPList.size() << "\n";
    insertAlleleFreqs(SNPList, rsIDs);
    cout << SNPList.size() << " remaining SNPs\n";
}

void insertAlleleFreqs(vector<SNP> & SNPList, const vector<dbSNP> & rsIDs)
{
    vector<SNP> newSNPList;
    cout << "Inserting Allele Frequencies\n";
    cout << "Enter allele frequency threshold:\n";
    double threshold;
    cin >> threshold;
    vector<string> AlleleFreqs;
    AlleleFreqs.resize(rsIDs.size(), "");
    map<string, string> codes;
    ifstream fin_code("./NCBIFiles/Allele.bcp");
    string line;
    getline(fin_code, line);
    string ID;
    while(!fin_code.eof())
    {
        ID = line.substr(0, line.find('\t'));
        line = line.substr(line.find('\t') + 1);
        codes[ID] = line.substr(0, line.find('\t'));
        getline(fin_code, line);
    }
    ifstream fin("./NCBIFiles/SNPAlleleFreq.bcp");
    getline(fin, line);
    int count = 0, id;
    string s;
    while(!fin.eof() && count < (int)rsIDs.size())
    {
        id = atoi(line.substr(0, line.find('\t')).c_str());//rsID
        if(id > rsIDs[count].id)
        {//if the rsID is too big, compare list of codes from when they were equal
            count++;
            continue;
        }
        if(id < rsIDs[count].id)
        {//if the rsID is too small, go on to the next line
            getline(fin, line);
            continue;
        }
        if(id == rsIDs[count].id)
        {//if the rsID fits, add to stack
            line = line.substr(line.find('\t')+1);
            s = line.substr(0, line.find('\t'));
            if(codes.find(s) != codes.end())
            {
                line = line.substr(line.find('\t') + 1);
                if(AlleleFreqs[count] != "")
                    AlleleFreqs[count] += ",";
                AlleleFreqs[count] += codes[s] + ":" + line.substr(0, line.find('\t'));
            }
            getline(fin, line);
        }
    }
    bool is_inserted;
    double sum;
    map<char, double> alleleCounts;
    char c;
    string alleles;
    double mutCount;
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(SNPList[i].dbSNP.size() < 2)
        {
            newSNPList.push_back(SNPList[i]);
            continue;
        }
        is_inserted = false;
        id = atoi((SNPList[i].dbSNP.substr(2)).c_str());
        for(int j = 0; (!is_inserted) && j < (int)rsIDs.size(); j++)
        {
            if(rsIDs[j].id == id)
            {
                is_inserted = true;
                SNPList[i].AlleleCounts = AlleleFreqs[j];
                alleles = AlleleFreqs[j];
                alleleCounts.clear();
                if(alleles.size() == 0)
                {
                    newSNPList.push_back(SNPList[i]);
                }
                else if((int)alleles.find(SNPList[i].refNT) != -1 && (int)alleles.find(SNPList[i].mutNT) == -1)
                {
                    newSNPList.push_back(SNPList[i]);
                }
                else if((int)alleles.find(flipGenotypes(SNPList[i].refNT)) != -1 && (int)alleles.find(flipGenotypes(SNPList[i].mutNT)) == -1)
                {
                    newSNPList.push_back(SNPList[i]);
                }
                else
                {
                    sum = 0;
                    while((int)alleles.find(',') != -1)
                    {
                        c = alleles[0];
                        alleles = alleles.substr(2);
                        alleleCounts[c] = atof(alleles.substr(0, alleles.find(',')).c_str());
                        sum += alleleCounts[c];
                        alleles = alleles.substr(alleles.find(',') + 1);
                    }
                    c = alleles[0];
                    alleles = alleles.substr(2);
                    alleleCounts[c] = atof(alleles.c_str());
                    sum += alleleCounts[c];
                    if(alleleCounts.find(SNPList[i].mutNT) != alleleCounts.end())
                    {
                        mutCount = alleleCounts[SNPList[i].mutNT];
                        if(SNPList[i].refNT == flipGenotypes(SNPList[i].mutNT) && alleleCounts.find(SNPList[i].refNT) != alleleCounts.end())
                        {
                            mutCount = min(mutCount, alleleCounts[SNPList[i].refNT]);
                        }
                        if(mutCount / sum <= threshold)
                        {
                            newSNPList.push_back(SNPList[i]);
                        }
                    }
                    else if(alleleCounts.find(flipGenotypes(SNPList[i].mutNT)) != alleleCounts.end())
                    {
                        mutCount = alleleCounts[flipGenotypes(SNPList[i].mutNT)];
                        if(SNPList[i].refNT == flipGenotypes(SNPList[i].mutNT) && alleleCounts.find(SNPList[i].refNT) != alleleCounts.end())
                        {
                            mutCount = min(mutCount, alleleCounts[SNPList[i].refNT]);
                        }
                        else if(SNPList[i].refNT == flipGenotypes(SNPList[i].mutNT))
                        {
                            mutCount = 0;
                        }
                        if(mutCount / sum <= threshold)
                        {
                            newSNPList.push_back(SNPList[i]);
                        }
                    }
                }
            }
        }
    }
    SNPList = newSNPList;
    cout << SNPList.size() << "\n";
}

void inputDomains(vector<SNP> & SNPList)
{
    cout << "Inputting Domains\n";
    map<string, vector<domain> > domainList;
    ifstream fin("./NCBIFiles/human.protein.gpff");
    string line;
    string transcriptID;
    domain d;
    while(!fin.eof())
    {
        while((int)line.find("REFSEQ: accession") == -1)
        {
            getline(fin, line);
        }
        transcriptID = line.substr(line.find("accession") + 10);
        if((int)transcriptID.find(".") != -1)
            transcriptID = transcriptID.substr(0, transcriptID.find("."));
        while(!fin.eof() && (int)line.find("LOCUS") == -1)
        {
            while(!fin.eof() && (int)line.find("LOCUS") == -1 && (int)line.find("FEATURES") == -1)
            {
                getline(fin, line);
            }
            while(!fin.eof() && (int)line.find("LOCUS") == -1)
            {
                while(!fin.eof() && (int)line.find("LOCUS") == -1 && ((int)line.find("Region") == -1 || line.find("Region") != line.find_first_not_of(" ")))
                    getline(fin, line);
                if(fin.eof() || (int)line.find("LOCUS") != -1)
                    continue;
                line = line.substr(line.find("Region")+6);
                line = line.substr(line.find_first_of("0123456789"));
                if((int)line.find_first_not_of("0123456789") != -1)
                {
                    d.start = atoi(line.substr(0, line.find_first_not_of("0123456789")).c_str());
                    line = line.substr(line.find_first_not_of("0123456789"));
                    line = line.substr(line.find_first_of("0123456789"));
                    d.end = atoi(line.c_str());
                }
                else
                {
                    d.start = atoi(line.c_str());
                    d.end = atoi(line.c_str());
                }
                getline(fin, line);
                if((int)line.find("region_name")!=-1)
                {
                    line = line.substr(line.find("\"")+1);
                    d.name = line.substr(0, line.find("\""));
                }
                getline(fin, line);
                domainList[transcriptID].push_back(d);
            }
        }
    }
    cout << domainList.size() << "\n";
    vector<domain> dList;
    string s, e;
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(domainList.find(SNPList[i].transcriptID) != domainList.end())
        {
            dList = domainList[SNPList[i].transcriptID];
            for(int j = 0; j < (int)dList.size(); j++)
            {
                if(SNPList[i].protPos >= dList[j].start && SNPList[i].protPos <= dList[j].end)
                {
                    stringstream ss1, ss2;
                    ss1 << dList[j].start;
                    ss1 >> s;
                    ss2 << dList[j].end;
                    ss2 >> e;
                    SNPList[i].Domains += "," + dList[j].name + "(" + s + "-" + e + ")";
                }
            }
        }
        if(SNPList[i].Domains.size() > 0)
        {
            SNPList[i].Domains = SNPList[i].Domains.substr(1);
        }
    }
}

void inputDomains(vector<vector<indel> > & indelList)
{
    cout << "Inputting Domains\n";
    map<string, vector<domain> > domainList;
    ifstream fin("./NCBIFiles/human.protein.gpff");
    string line;
    string transcriptID;
    domain d;
    while(!fin.eof())
    {
        while((int)line.find("REFSEQ: accession") == -1)
        {
            getline(fin, line);
        }
        transcriptID = line.substr(line.find("accession") + 10);
        if((int)transcriptID.find(".") != -1)
            transcriptID = transcriptID.substr(0, transcriptID.find("."));
        while(!fin.eof() && (int)line.find("LOCUS") == -1)
        {
            while(!fin.eof() && (int)line.find("LOCUS") == -1 && (int)line.find("FEATURES") == -1)
            {
                getline(fin, line);
            }
            while(!fin.eof() && (int)line.find("LOCUS") == -1)
            {
                while(!fin.eof() && (int)line.find("LOCUS") == -1 && ((int)line.find("Region") == -1 || line.find("Region") != line.find_first_not_of(" ")))
                    getline(fin, line);
                if(fin.eof() || (int)line.find("LOCUS") != -1)
                    continue;
                line = line.substr(line.find("Region")+6);
                line = line.substr(line.find_first_of("0123456789"));
                if((int)line.find_first_not_of("0123456789") != -1)
                {
                    d.start = atoi(line.substr(0, line.find_first_not_of("0123456789")).c_str());
                    line = line.substr(line.find_first_not_of("0123456789"));
                    line = line.substr(line.find_first_of("0123456789"));
                    d.end = atoi(line.c_str());
                }
                else
                {
                    d.start = atoi(line.c_str());
                    d.end = atoi(line.c_str());
                }
                getline(fin, line);
                if((int)line.find("region_name")!=-1)
                {
                    line = line.substr(line.find("\"")+1);
                    d.name = line.substr(0, line.find("\""));
                }
                getline(fin, line);
                domainList[transcriptID].push_back(d);
            }
        }
    }
    cout << domainList.size() << "\n";
    vector<domain> dList;
    string s, e;
    for(int i = 0; i < (int)indelList.size(); i++)
    {
        for(int j = 0; j < (int)indelList[i].size(); j++)
        {
            if(domainList.find(indelList[i][j].transcript) != domainList.end())
            {
                dList = domainList[indelList[i][j].transcript];
                for(int k = 0; k < (int)dList.size(); k++)
                {
                    if((indelList[i][j].ProtStart >= dList[k].start && indelList[i][j].ProtStart <= dList[k].end)
                       || (indelList[i][j].ProtEnd >= dList[k].start && indelList[i][j].ProtEnd <= dList[k].end))//assumption: indel length << exon length
                    {
                        stringstream ss1, ss2;
                        ss1 << dList[k].start;
                        ss1 >> s;
                        ss2 << dList[k].end;
                        ss2 >> e;
                        indelList[i][j].Domains += "," + dList[k].name + "(" + s + "-" + e + ")";
                    }
                }
            }
            if(indelList[i][j].Domains.size() > 0)
            {
                indelList[i][j].Domains = indelList[i][j].Domains.substr(1);
            }
        }
    }
}

void filterGenic(vector<vector<indel> > & indelList, const vector<gene> & geneList)
{
    vector<vector<indel> > newIndelList(indelList.size());
    int k, kprime;
    bool reached_chr, is_inserted;
    int start, end;
    int cds_start, cds_end;
    int begin, close;
    cout << "Filtering Genic Regions\n";
    for(int i = 0; i < (int)indelList.size(); i++)
    {
        k = 0;
        for(int j = 0; j < (int)indelList[i].size(); j++)
        {
            start = indelList[i][j].start+1;
            end = indelList[i][j].ref.size() > indelList[i][j].alt.size() ? indelList[i][j].start+indelList[i][j].ref.size() : indelList[i][j].start+1;
            reached_chr = (k < (int)geneList.size()) ? (geneList[k].chr == indelList[i][j].chr) : true;
            is_inserted = false;
            while(k < (int)geneList.size() && (!reached_chr || (geneList[k].chr == indelList[i][j].chr && geneList[k].end < start)))
            {
                if(!reached_chr && geneList[k].chr == indelList[i][j].chr)
                {
                    reached_chr = true;
                }
                k++;
            }
            kprime = k;
            cds_start = 0;
            cds_end = 0;
            while(!is_inserted && kprime < (int)geneList.size() && geneList[kprime].start <= start)
            {
                if(geneList[kprime].strand)
                {
                    for(int l = 0; !is_inserted && l < (int)geneList[kprime].exonList.size(); l++)
                    {
                        if(geneList[kprime].exonList[l].start-5 <= end && geneList[kprime].exonList[l].end+5 >= start && (l < (int)geneList[kprime].exonList.size() || geneList[kprime].exonList[l].end >= start))
                        {
                            //cout << "Indel on chromosome " << indelList[i][j].chr << ", positions " << start << "-" << end << " overlaps with exon " << l << " on gene " << geneList[kprime].name << " (positions " << geneList[kprime].exonList[l].start << "-" << geneList[kprime].exonList[l].end << ")\n";
                            if(geneList[kprime].exonList[l].start > end || geneList[kprime].exonList[l].end < start)
                                indelList[i][j].is_splice = true;
                            else if(end < geneList[kprime].cds_start || start > geneList[kprime].cds_end)
                                indelList[i][j].is_utr = true;
                            else
                            {
                                cds_end = cds_start;
                                begin = max(geneList[kprime].exonList[l].start, geneList[kprime].cds_start);
                                close = min(geneList[kprime].exonList[l].end, start);
                                cds_start += max(0, close - begin);
                                indelList[i][j].ProtStart = cds_start/3+1;
                                begin = max(geneList[kprime].exonList[l].start, geneList[kprime].cds_start);
                                close = min(geneList[kprime].exonList[l].end, end);
                                cds_end += max(0, close - begin);
                                indelList[i][j].ProtEnd = cds_end/3+1;
                            }
                            indelList[i][j].gene = geneList[kprime].name.substr(0, geneList[kprime].name.find('('));
                            indelList[i][j].transcript = geneList[kprime].name.substr(geneList[kprime].name.find('(') + 1);
                            indelList[i][j].transcript = indelList[i][j].transcript.substr(0, indelList[i][j].transcript.size() - 1);//closing )
                            newIndelList[i].push_back(indelList[i][j]);
                            is_inserted = true;
                        }
                        else
                        {
                            begin = max(geneList[kprime].exonList[l].start, geneList[kprime].cds_start);
                            close = min(geneList[kprime].exonList[l].end, geneList[kprime].cds_end);
                            cds_start += max(0, close - begin);
                        }
                    }
                }
                else
                {
                    for(int l = geneList[kprime].exonList.size() - 1; !is_inserted && l >= 0; l--)
                    {
                        if(geneList[kprime].exonList[l].start-5 <= end && geneList[kprime].exonList[l].end+5 >= start && (l > 0 || geneList[kprime].exonList[l].start <= start))
                        {
                            //cout << "Indel on chromosome " << indelList[i][j].chr << ", positions " << start << "-" << end << " overlaps with exon " << l << " on gene " << geneList[kprime].name << " (positions " << geneList[kprime].exonList[l].start << "-" << geneList[kprime].exonList[l].end << ")\n";
                            if(geneList[kprime].exonList[l].start > end || geneList[kprime].exonList[l].end < start)
                                indelList[i][j].is_splice = true;
                            else if(end < geneList[kprime].cds_start || start > geneList[kprime].cds_end)
                                indelList[i][j].is_utr = true;
                            else
                            {
                                cds_end = cds_start;
                                begin = max(geneList[kprime].exonList[l].start, end);
                                close = min(geneList[kprime].exonList[l].end, geneList[kprime].cds_end);
                                cds_start += max(0, close - begin);
                                indelList[i][j].ProtStart = cds_start/3+1;
                                begin = max(geneList[kprime].exonList[l].start, start);
                                close = min(geneList[kprime].exonList[l].end, geneList[kprime].cds_end);
                                cds_end += max(0, close - begin);
                                indelList[i][j].ProtEnd = cds_end/3+1;
                            }
                            indelList[i][j].gene = geneList[kprime].name.substr(0, geneList[kprime].name.find('('));
                            indelList[i][j].transcript = geneList[kprime].name.substr(geneList[kprime].name.find('(') + 1);
                            indelList[i][j].transcript = indelList[i][j].transcript.substr(0, indelList[i][j].transcript.size() - 1);//closing )
                            newIndelList[i].push_back(indelList[i][j]);
                            is_inserted = true;
                        }
                        else
                        {
                            begin = max(geneList[kprime].exonList[l].start, geneList[kprime].cds_start);
                            close = min(geneList[kprime].exonList[l].end, geneList[kprime].cds_end);
                            cds_start += max(0, close - begin);
                        }
                    }
                }
                kprime++;
            }
        }
    }
    indelList = newIndelList;
}

//function filterUnique()
//Filters out indels present in multiple patients with different conditions
//Input: List of lists of indels indelList (first dimension = patients, 2nd dimension = indels)
//Output: Reduced list of indels
void filterUnique(vector<vector<indel> > & indelList, const double & freq)
{
    cout << "Filtering indel uniqueness\n";
    vector<vector<indel> > newIndelList(indelList.size());
    bool passChr;
    for(int i = 0; i < (int)indelList.size(); i++)
    {
        cout << i << "\n";
        for(int j = 0; j < (int)indelList[i].size(); j++)
        {//for each indel
            for(int k = indelList.size() - 1; indelList[i][j].is_unique && k > i; k--)
            {
                passChr = false;
                for(int l = 0; !passChr && l < (int)indelList[k].size(); l++)
                {//for each indel after it
                    if(indelList[i][j].chr == indelList[k][l].chr && indelList[i][j].start < indelList[k][l].start)
                        passChr = true;
                    if(indelList[i][j].chr == indelList[k][l].chr && indelList[i][j].start == indelList[k][l].start && indelList[i][j].ref == indelList[k][l].ref
                       && indelList[i][j].alt == indelList[k][l].alt)
                    {//if it's identical, toss out
                        indelList[k][l].is_unique = false;
                        if((atof(indelList[k][l].ac.c_str())/atof(indelList[k][l].dp.c_str())) >= freq)
                        {
                            indelList[i][j].shares = "," + patientList[k] + indelList[i][j].shares;
                        }
                        else
                        {
                            indelList[i][j].het = "," + patientList[k] + indelList[i][j].het;
                        }
                        indelList[i][j].full_shares = indelList[i][j].full_shares.substr(0, k) + indelList[k][l].ac + "/" + indelList[k][l].dp + indelList[i][j].full_shares.substr(k);
                    }
                }
            }
            if(indelList[i][j].is_unique)
            {
                if((atof(indelList[i][j].ac.c_str())/atof(indelList[i][j].dp.c_str())) >= freq)
                {
                    indelList[i][j].shares = patientList[i] + indelList[i][j].shares;
                    if(indelList[i][j].het != "")
                        indelList[i][j].het = indelList[i][j].het.substr(1);
                }
                else
                {
                    indelList[i][j].het = patientList[i] + indelList[i][j].het;
                    if(indelList[i][j].shares != "")
                        indelList[i][j].shares = indelList[i][j].shares.substr(1);
                }
                indelList[i][j].full_shares = indelList[i][j].full_shares.substr(0, i) + indelList[i][j].ac + "/" + indelList[i][j].dp + indelList[i][j].full_shares.substr(i);
                if(indelList[i][j].shares != "")
                {
                    newIndelList[i].push_back(indelList[i][j]);
                }
            }
        }
    }
    indelList = newIndelList;
}

//function outputIndels()
//Outputs final list of indels
//Input: Output stream fout, List of lists of indels indelList (first dimension = patients, 2nd dimension = indels)
//Output: Outputted list of indels in SIFT-legible format
void outputIndels(ofstream & fout, vector<vector<indel> > & indelList)
{
    for(int i = 0; i < (int)indelList.size(); i++)
    {
        for(int j = 0; j < (int)indelList[i].size(); j++)
        {
            fout << indelList[i][j].chr.substr(3) << ",";
            fout << indelList[i][j].start << ",";
            if(indelList[i][j].ref.size() <= indelList[i][j].alt.size())//insertion
            {
                assert(indelList[i][j].ref.size() == 1);
                fout << indelList[i][j].start << ",1,";
                fout << indelList[i][j].alt.substr(1) << ",";
            }
            else
            {
                assert(indelList[i][j].alt.size() == 1);
                fout << indelList[i][j].start + indelList[i][j].ref.size() - 1 << ",1,/,";//deletion
            }
            fout << "\t" << indelList[i][j].shares << "\t" << indelList[i][j].het << "\t";
            if(indelList[i][j].is_utr)
                fout << "utr\t";
            else if(indelList[i][j].is_splice)
                fout << "splice\t";
            else fout << "CDS\t";
            fout << indelList[i][j].gene << "(" << indelList[i][j].transcript << ")\t" << indelList[i][j].ProtStart << "-" << indelList[i][j].ProtEnd << "\t" << indelList[i][j].Domains << "\t";
            fout << indelList[i][j].full_shares << "\n";
        }
    }
}

//function filterIndels()
//Filters out non-homozygous indels
//Input: List of lists of indels indelList (first dimension = patients, 2nd dimension = indels), frequency at which we can classify indel as homozygous (~75%?)
//Output: Reduced list of lists of indels
void filterIndels(vector<vector<indel> > & indelList, const double & freq)
{
    cout << "Filtering indel frequencies\n";
    vector<vector<indel> > newIndelList(indelList.size());
    for(int i = 0; i < (int)indelList.size(); i++)
    {
        for(int j = 0; j < (int)indelList[i].size(); j++)
        {
            if((atof(indelList[i][j].ac.c_str())/atof(indelList[i][j].dp.c_str())) >= freq)
                newIndelList[i].push_back(indelList[i][j]);
        }
    }
    indelList = newIndelList;
}

//function inputIndels()
//Inputs indels into indelList
//Input: Empty list of lists of indels indelList (first dimension = patients, 2nd dimension = indels)
//Output: Filled indelList
void inputIndels(vector<vector<indel> > & indelList)
{
    string line = "";
    ifstream fin("./Broad.vcf");
    while((int)line.find("#CHROM") == -1)
    {
        getline(fin, line);
    }
    for(int i = 0; i < 9; i++)
    {
        line = line.substr(line.find('\t') + 1);
    }
    string full_shares = "";
    while((int)line.find('\t') != -1)
    {
        full_shares += "\t";
        patientList.push_back(line.substr(0, line.find('\t')));
        line = line.substr(line.find('\t') + 1);
    }
    cout << patientList.size() << " patients\n";
    indelList.resize(patientList.size());
    indel l;
    fin.close();
    for(int i = 0; i < (int)patientList.size(); i++)
    {//for all patients
        cout << "Inputting " << patientList[i] << " indels\n";
        fin.open(("./IntermediateFiles/Ciliopathies_Whole_Exome_Gleeson_" + patientList[i] + "/" + patientList[i] + ".indels.vcf").c_str());
        while((int)line.find("#CHROM") == -1)
        {
            getline(fin, line);
        }
        getline(fin, line);
        while(!fin.eof())
        {
            l.chr = "chr" + line.substr(0, line.find('\t'));
            line = line.substr(line.find('\t') + 1);//chrom
            l.start = atoi(line.substr(0, line.find('\t')).c_str());
            line = line.substr(line.find('\t') + 1);//pos
            line = line.substr(line.find('\t') + 1);//ID
            l.ref = line.substr(0, line.find('\t'));
            line = line.substr(line.find('\t') + 1);//ref
            l.alt = line.substr(0, line.find('\t'));
            l.is_unique = true;
            line = line.substr(line.find("AC=") + 3);
            l.ac = line.substr(0, line.find(','));
            line = line.substr(line.find("DP=") + 3);
            l.dp = line.substr(0, line.find(';'));
            l.is_splice = false;
            l.is_utr = false;
            l.full_shares = full_shares;
            l.ProtStart = 0;
            l.ProtEnd = 0;
            indelList[i].push_back(l);
            getline(fin, line);
        }
        fin.close();
    }
}

//function inputGeneList()
//Helper function of analyzeLowCovData
//Takes RefSeq's gene list & puts them in a data structure (geneList)
//Input: File containing list of RefSeq genes
//Output: list of genes (geneList)
void inputGeneList(vector<gene> & geneList)
{
    cout << "Inputting Genes\n";
    ifstream fin_genes("./NCBIFiles/RefSeqGeneList_hg19.txt");
    string line_gene, sub_gene, exon_starts, exon_ends;
    getline(fin_genes, line_gene);
    getline(fin_genes, line_gene);
    int exon_count;
    while(!fin_genes.eof())
    {//input gene name, chromosome, etc. from line
        gene g;
        sub_gene = line_gene.substr(line_gene.find('\t') + 1);
        g.name = "(" + sub_gene.substr(0, sub_gene.find('\t')) + ")";
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        g.chr = sub_gene.substr(0, sub_gene.find('\t'));
        g.strand = (sub_gene[sub_gene.find('\t') + 1] == '+');
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 3);
        g.start = atoi(sub_gene.substr(0, sub_gene.find('\t')).c_str());
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        g.end = atoi(sub_gene.substr(0, sub_gene.find('\t')).c_str());
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        g.cds_start = atoi(sub_gene.substr(0, sub_gene.find('\t')).c_str());
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        g.cds_end = atoi(sub_gene.substr(0, sub_gene.find('\t')).c_str());
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        exon_count = atoi(sub_gene.substr(0, sub_gene.find('\t')).c_str());
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        exon_starts = sub_gene.substr(0, sub_gene.find('\t'));
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        exon_ends = sub_gene.substr(0, sub_gene.find('\t'));
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        sub_gene = sub_gene.substr(sub_gene.find('\t') + 1);
        g.name = sub_gene.substr(0, sub_gene.find('\t')) + g.name;
        for(int i = 0; i < exon_count; i++)
        {//convert exons from comma separated strings to list
            exon e;
            e.start = atoi(exon_starts.substr(0, exon_starts.find(',')).c_str());
            exon_starts = exon_starts.substr(exon_starts.find(',') + 1);
            e.end = atoi(exon_ends.substr(0, exon_ends.find(',')).c_str());
            exon_ends = exon_ends.substr(exon_ends.find(',') + 1);
            g.exonList.push_back(e);
        }
        geneList.push_back(g);//add to the gene list
        getline(fin_genes, line_gene);
    }
}

//function inputSIFT()
//Inputs SIFT annotations into SNPSamples
//Input: List of SNPs SNPList, ID representing where the file is for SIFT
//Output: SNPList with better annotation
void inputSIFT(vector<SNP> & SNPList, const string & ID)
{
    cout << "Inputting SIFT\n";
    ifstream fin(("./" + ID + "_predictions.tsv").c_str());
    string line;
    getline(fin, line);
    int loc;
    string chr;
    getline(fin, line);
    int j = 0;
    string name;
    string dbsnp;
    while(!fin.eof() && j < (int)SNPList.size())
    {//while there are more SNPs
        chr = "chr" + line.substr(0, line.find(','));//input
        line = line.substr(line.find(',') + 1);
        loc = atoi(line.substr(0, line.find(',')).c_str());
        while(chr != SNPList[j].chr || loc != SNPList[j].genomicLoc)
        {
            j++;
        }
        if(j < (int)SNPList.size() && chr == SNPList[j].chr && loc == SNPList[j].genomicLoc)
        {//if match found, modify SIFT
            line = line.substr(line.find('\t') + 1);//coordinates
            line = line.substr(line.find('\t') + 1);//codons
            line = line.substr(line.find('\t') + 1);//transID
            line = line.substr(line.find('\t') + 1);//protID
            line = line.substr(line.find('\t') + 1);//subs
            line = line.substr(line.find('\t') + 1);//region
            dbsnp = line.substr(0, line.find('\t'));
            if((SNPList[j].dbSNP == "." && (int)dbsnp.find("rs") != -1))
            {
                dbsnp = dbsnp.substr(dbsnp.find("rs"));
                SNPList[j].dbSNP = dbsnp.substr(0, dbsnp.find_first_of("A:TCGEn\t"));
            }
            line = line.substr(line.find('\t') + 1);//dbSNP
            line = line.substr(line.find('\t') + 1);//type
            SNPList[j].SIFT = line.substr(0, line.find('\t'));
            line = line.substr(line.find('\t') + 1);//SIFT Prediction
            SNPList[j].SIFTProb = atof(line.substr(0, line.find('\t')).c_str());
            line = line.substr(line.find('\t') + 1);//Score
            line = line.substr(line.find('\t') + 1);//Median
            line = line.substr(line.find('\t') + 1);//#Seqs
            line = line.substr(line.find('\t') + 1);//Gene ID
            name = line.substr(0, line.find('\t'));
            line = line.substr(line.find('\t') + 1);//Gene Name
            if(line[0] != '\t' && name == SNPList[j].geneName)
            {
                SNPList[j].fullName = line.substr(0, line.find_first_of("[(\t") - 1);
                line = line.substr(0, line.find_last_of('\t'));//Comment
                line = line.substr(0, line.find_last_of('\t'));//CEU Freqs
                line = line.substr(0, line.find_last_of('\t'));//Avg Freqs
                SNPList[j].OMIM = line.substr(line.find_last_of('\t') + 1);
            }
            else
            {
                line = line.substr(0, line.find_last_of('\t'));//Comment
                line = line.substr(0, line.find_last_of('\t'));//CEU Freqs
                line = line.substr(0, line.find_last_of('\t'));//Avg Freqs
                if(line[0] != '\t')
                    SNPList[j].OMIM = line.substr(line.find_last_of('\t') + 1) + "(" + name + "?)";
            }
            j++;
            while(j > 0 && j < (int)SNPList.size() && SNPList[j].chr == SNPList[j-1].chr && SNPList[j].genomicLoc == SNPList[j-1].genomicLoc)
            {
                SNPList[j].dbSNP = SNPList[j-1].dbSNP;
                SNPList[j].SIFT = SNPList[j-1].SIFT;
                SNPList[j].fullName = SNPList[j-1].fullName;
                if(name == SNPList[j].geneName)
                    SNPList[j].OMIM = SNPList[j-1].OMIM.substr(0, SNPList[j-1].OMIM.find('('));
                else
                    SNPList[j].OMIM = SNPList[j-1].OMIM;
                j++;
            }
        }
        getline(fin, line);
    }
}


//function inputData()
//Inputs list of SNPs from file
//Input: file containing SNPs, file containing codons
//Output: List of SNPs with corrected codon annotations
void inputData(vector<SNP> & SNPList)
{
    map<string, char> codonTable;//maps 3 letter codons to 1 letter code
    ifstream fin_codon("./codontable.txt");
    cout << "Inputting Codon Table\n";
	string aa;
	getline(fin_codon, aa);
	while(!fin_codon.eof())
	{//while there are more amino acids to go
		string code = aa.substr(0, aa.find('\t'));
		if(code[0] != '*' && (code[0] < 'A' || code[0] > 'Z'))
            code[0] = 'I';//isoleucine is messed up in the file - fixed manually
		aa = aa.substr(aa.find('\t') + 1);
		string codon;
		while((int)aa.find(',') != -1)
		{//while there are more codons leading to the same amino acid
			codon = aa.substr(0, aa.find(','));
			codonTable[codon] = code[0];//add association to table
			aa = aa.substr(aa.find(',') + 2);
		}
		codonTable[aa] = code[0];
		getline(fin_codon, aa);
	}

    ifstream fin("./Broadsnps.vcf");
    cout << "Inputting SNP Data\n";
    string line;
    getline(fin, line);
    patientList.clear();
    while((int)line.find("CHROM") == -1)
    {//skip all the header files
        getline(fin, line);
    }
    vector<string> headers;
    int i = 0;
    while((int)line.find('\t') != -1)
    {//keep the column labels
        i++;
        if(i > 9)
        {
            patientList.push_back(line.substr(0, line.find('\t')));
        }
        headers.push_back(line.substr(0, line.find('\t')));
        line = line.substr(line.find('\t') + 1);
    }
    cout << patientList.size() << " patients\n";
    headers.push_back(line);
    getline(fin, line);
    string sub, sub2;
    baseCall bc;
    int num_records;
    string num;
    SNP s_prev;
    int tot = 0;
    while(!fin.eof())
    {//while there is a SNP, input all the SNP info
        SNP s;
        s.chr = "chr" + line.substr(0, line.find('\t'));
        sub = line.substr(line.find('\t') + 1);
        s.genomicLoc = atoi((sub.substr(0, sub.find('\t'))).c_str());
        sub = sub.substr(sub.find('\t') + 1);
        s.dbSNP = sub.substr(0, sub.find('\t'));
        sub = sub.substr(sub.find('\t') + 1);
        s.refNT = sub[0];
        s.mutNT = sub[2];
        sub = sub.substr(4);
        sub = sub.substr(sub.find('\t') + 1);//qual
        sub = sub.substr(sub.find('\t') + 1);//filter
        s.info = sub.substr(0, sub.find('\t'));
        s.is_splice=false;
        s.GERP = 0;
        s.polyPhen = "N/A";
        s.codon = "";
        s.shares = "";
        s.unsure = "";
        s.protPos = 0;
        s.SIFT = "Unknown";
        s.OMIM = "";
        sub2 = sub.substr(sub.find("GT:AD:DP:GQ:PL") + 15);
        s.full_shares = sub2;
        while((int)sub2.find('\t') != -1)
        {//input all the genotpe info
            string zygosityPair = sub2.substr(0,3);
            if(zygosityPair == "./.")
            {
                bc.genotype = 0;
                bc.depth = 0;
                bc.quality = 0;
                bc.AD = "";
                bc.PL = "";
                s.calls.push_back(bc);
            }
            else if(zygosityPair == "0/0")
            {
                bc.genotype = 1;
                sub2 = sub2.substr(4);
                bc.AD = sub2.substr(0, sub2.find(':'));
                sub2 = sub2.substr(sub2.find(':')+1);
                bc.depth = atoi(sub2.substr(0, sub2.find(':')).c_str());
                sub2 = sub2.substr(sub2.find(':')+1);
                bc.quality = atof(sub2.substr(0, sub2.find(':')).c_str());
                sub2 = sub2.substr(sub2.find(':') + 1);
                bc.PL = sub2.substr(0, sub2.find('\t'));
                s.calls.push_back(bc);
            }
            else if(zygosityPair == "1/1")
            {
                bc.genotype = 2;
                sub2 = sub2.substr(4);
                bc.AD = sub2.substr(0, sub2.find(':'));
                sub2 = sub2.substr(sub2.find(':')+1);
                bc.depth = atoi(sub2.substr(0, sub2.find(':')).c_str());
                sub2 = sub2.substr(sub2.find(':')+1);
                bc.quality = atof(sub2.substr(0, sub2.find(':')).c_str());
                sub2 = sub2.substr(sub2.find(':') + 1);
                bc.PL = sub2.substr(0, sub2.find('\t'));
                s.calls.push_back(bc);
            }
            else if(zygosityPair == "0/1" || zygosityPair == "1/0")
            {
                bc.genotype = 3;
                sub2 = sub2.substr(4);
                bc.AD = sub2.substr(0, sub2.find(':'));
                sub2 = sub2.substr(sub2.find(':')+1);
                bc.depth = atoi(sub2.substr(0, sub2.find(':')).c_str());
                sub2 = sub2.substr(sub2.find(':')+1);
                bc.quality = atof(sub2.substr(0, sub2.find(':')).c_str());
                sub2 = sub2.substr(sub2.find(':') + 1);
                bc.PL = sub2.substr(0, sub2.find('\t'));
                s.calls.push_back(bc);
            }
            else
            {
                cout << "Problem: zygosityPair: " << zygosityPair << endl;
                //cout << sub2 << "\n";//shouldn't be reached
            }
            sub2 = sub2.substr(sub2.find('\t') + 1);
        }
        if(sub2.find("./.") == 0)
        {//take care of the last genotype
            bc.genotype = 0;
            bc.depth = 0;
            bc.quality = 0;
            bc.AD = "";
            bc.PL = "";
            s.calls.push_back(bc);
        }
        else if(sub2.find("0/0") == 0)
        {
            bc.genotype = 1;
            sub2 = sub2.substr(4);
            bc.AD = sub2.substr(0, sub2.find(':'));
            sub2 = sub2.substr(sub2.find(':')+1);
            bc.depth = atoi(sub2.substr(0, sub2.find(':')).c_str());
            sub2 = sub2.substr(sub2.find(':')+1);
            bc.quality = atof(sub2.substr(0, sub2.find(':')).c_str());
            sub2 = sub2.substr(sub2.find(':') + 1);
            bc.PL = sub2.substr(0, sub2.find('\t'));
            s.calls.push_back(bc);
        }
        else if(sub2.find("1/1") == 0)
        {
            bc.genotype = 2;
            sub2 = sub2.substr(4);
            bc.AD = sub2.substr(0, sub2.find(':'));
            sub2 = sub2.substr(sub2.find(':')+1);
            bc.depth = atoi(sub2.substr(0, sub2.find(':')).c_str());
            sub2 = sub2.substr(sub2.find(':')+1);
            bc.quality = atof(sub2.substr(0, sub2.find(':')).c_str());
            sub2 = sub2.substr(sub2.find(':') + 1);
            bc.PL = sub2.substr(0, sub2.find('\t'));
            s.calls.push_back(bc);
        }
        else if(sub2.find("0/1") == 0 || sub2.find("1/0") == 0)
        {
            bc.genotype = 3;
            sub2 = sub2.substr(4);
            bc.AD = sub2.substr(0, sub2.find(':'));
            sub2 = sub2.substr(sub2.find(':')+1);
            bc.depth = atoi(sub2.substr(0, sub2.find(':')).c_str());
            sub2 = sub2.substr(sub2.find(':')+1);
            bc.quality = atof(sub2.substr(0, sub2.find(':')).c_str());
            sub2 = sub2.substr(sub2.find(':') + 1);
            bc.PL = sub2.substr(0, sub2.find('\t'));
            s.calls.push_back(bc);
        }
        else
        {
            cout << sub2.substr(0, 3) << "\n";
        }
       
	//assert(s.calls.size() == patientList.size());
        if((int)sub.find("refseq.numMatchingRecords") == -1)
        {
            tot++;
            if((int)sub.find("refseq.name=") == -1)
                s.transcriptID = "";
            else
            {
                sub2 = sub.substr(sub.find("refseq.name="));
                sub2 = sub2.substr(sub2.find('=')+1);
                s.transcriptID = sub2.substr(0, sub2.find(';'));
            }
            if((int)sub.find("refseq.name2=") == -1)
            {
                s.geneName = "";
            }
            else
            {
                sub2 = sub.substr(sub.find("refseq.name2="));
                sub2 = sub2.substr(sub2.find('=')+1);
                s.geneName = sub2.substr(0, sub2.find(';'));
            }
            if((int)sub.find("refseq.transcriptStrand=")==-1)
            {
                s.strand = true;
            }
            else
            {
                sub2 = sub.substr(sub.find("refseq.transcriptStrand="));
                sub2 = sub2.substr(sub2.find('=')+1);
                s.strand = (sub2[0] == '+') ? true : false;
            }
            if((int)sub.find("refseq.codingCoordStr") == -1 || (int)sub.find("refseq.positionType=CDS")==-1)
            {
                s.codonPos = 0;
                s.codonStart = 0;
                if((int)sub.find("haplotypeReference") != -1 && (int)sub.find("haplotypeAlternate") != -1)
                {
                    sub2 = sub.substr(sub.find("haplotypeReference"));
                    sub2 = sub2.substr(sub2.find('=')+1);
                    sub2 = sub.substr(sub.find("haplotypeAlternate"));
                    sub2 = sub2.substr(sub2.find('=') + 1);
                }
            }
            else
            {
                sub2 = sub.substr(sub.find("refseq.codingCoordStr"));
                sub2 = sub2.substr(sub2.find("=c.") + 3);
                s.codonPos = atoi(sub2.substr(0, sub2.find_first_not_of("0123456789")).c_str());
                sub2 = sub2.substr(sub2.find_first_not_of("0123456789"));
                s.codonStart = s.codonPos - (s.codonPos % 3);
            }
            if((int)sub.find("refseq.proteinCoordStr") == -1 || (int)sub.find("refseq.positionType=CDS") == -1)
            {
                s.refAA = 'X';
                s.mutAA = 'X';
                s.protPos = 0;
            }
            else
            {
                sub2 = sub.substr(sub.find("refseq.proteinCoordStr"));
                sub2 = sub2.substr(sub2.find("=p.") + 3);
                s.refAA = sub2[0];
                sub2 = sub2.substr(1);
                s.protPos = atoi(sub2.substr(0, sub2.find_first_not_of("0123456789")).c_str());
                sub2 = sub2.substr(sub2.find_first_not_of("0123456789"));
                s.mutAA = sub2[0];
            }
            if((int)sub.find("refseq.referenceCodon") != -1)
            {
                sub2 = sub.substr(sub.find("refseq.referenceCodon"));
                sub2 = sub2.substr(sub2.find("=") + 1);
                s.codon = sub2.substr(0,3);
            }
            s.is_splice = false;
            if((int)sub2.find("refseq.spliceDist") != -1)
            {
                sub2 = sub.substr(sub.find("refseq.spliceDist"));
                sub2 = sub2.substr(sub2.find('=') + 1);
                s.spliceDist = atoi(sub2.substr(0, sub2.find(';')).c_str());
                if(s.spliceDist <= 5 && s.spliceDist >= -5)
                {
                    s.is_splice = true;
                }
            }
            s.is_utr = false;
            if((int)sub.find("refseq.positionType=utr") != -1)
            {
                s.is_utr = true;
            }
            SNPList.push_back(s);
        }
        else
        {
            sub2 = sub.substr(sub.find("refseq.numMatchingRecords=") + 26);
            num_records = atoi(sub2.substr(0, sub2.find(';')).c_str());
            tot += num_records;
            for(int j = 1; j <= num_records; j++)
            {
                stringstream ss;
                ss << j;
                ss >> num;
                if((int)sub.find("refseq.name_" + num + "=") == -1)
                {
                    s.transcriptID = "";
                }
                else
                {
                    sub2 = sub.substr(sub.find("refseq.name_" + num + "="));
                    sub2 = sub2.substr(sub2.find('=')+1);
                    s.transcriptID = sub2.substr(0, sub2.find(';'));
                }
                if((int)sub.find("refseq.name2_" + num + "=") == -1)
                {
                    s.geneName = "";
                }
                else
                {
                    sub2 = sub.substr(sub.find("refseq.name2_" + num + "="));
                    sub2 = sub2.substr(sub2.find('=')+1);
                    s.geneName = sub2.substr(0, sub2.find(';'));
                }
                if((int)sub.find("refseq.transcriptStrand_" + num + "=") == -1)
                {
                    s.strand = true;
                }
                else
                {
                    sub2 = sub.substr(sub.find("refseq.transcriptStrand_" + num + "="));
                    sub2 = sub2.substr(sub2.find('=')+1);
                    s.strand = (sub2[0] == '+') ? true : false;
                }
                if((int)sub.find("refseq.codingCoordStr_" + num + "=") == -1 || (int)sub.find("refseq.positionType_" + num + "=CDS") == -1)
                {
                    s.codonPos = 0;
                    s.codonStart = 0;
                    if((int)sub.find("haplotypeReference_"+num+"=") != -1 && (int)sub.find("haplotypeAlternate_"+num+"=") != -1)
                    {
                        sub2 = sub.substr(sub.find("haplotypeReference_" + num+"="));
                        sub2 = sub2.substr(sub2.find('=')+1);
                        sub2 = sub.substr(sub.find("haplotypeAlternate_" + num+"="));
                        sub2 = sub2.substr(sub2.find('=') + 1);
                    }

                }
                else
                {
                    sub2 = sub.substr(sub.find("refseq.codingCoordStr_" + num + "="));
                    sub2 = sub2.substr(sub2.find("=c.") + 3);
                    s.codonPos = atoi(sub2.substr(0, sub2.find_first_not_of("0123456789")).c_str());
                    sub2 = sub2.substr(sub2.find_first_not_of("0123456789"));
                    s.codonStart = s.codonPos - (s.codonPos % 3);
                }
                if((int)sub.find("refseq.proteinCoordStr_" + num + "=") == -1 || (int)sub.find("refseq.positionType_" + num + "=CDS") == -1)
                {
                    s.refAA = 'X';
                    s.mutAA = 'X';
                    s.protPos = 0;
                }
                else
                {
                    sub2 = sub.substr(sub.find("refseq.proteinCoordStr_" + num + "="));
                    sub2 = sub2.substr(sub2.find("=p.") + 3);
                    s.refAA = sub2[0];
                    sub2 = sub2.substr(1);
                    s.protPos = atoi(sub2.substr(0, sub2.find_first_not_of("0123456789")).c_str());
                    sub2 = sub2.substr(sub2.find_first_not_of("0123456789"));
                    s.mutAA = sub2[0];
                }
                if((int)sub.find("refseq.referenceCodon_" + num + "=") != -1)
                {
                    sub2 = sub.substr(sub.find("refseq.referenceCodon_" + num + "="));
                    sub2 = sub2.substr(sub2.find("=") + 1);
                    s.codon = sub2.substr(0,3);
                }
                s.is_splice = false;
                if((int)sub2.find("refseq.spliceDist_"+num+"=") != -1)
                {
                    sub2 = sub.substr(sub.find("refseq.spliceDist_" + num + "="));
                    sub2 = sub2.substr(sub2.find('=') + 1);
                    s.spliceDist = atoi(sub2.substr(0, sub2.find(';')).c_str());
                    if(s.spliceDist <= 5 && s.spliceDist >= -5)
                    {
                        s.is_splice = true;
                    }
                }
                s.is_utr = false;
                if((int)sub.find("refseq.positionType_" + num + "=utr") != -1)
                {
                    s.is_utr = true;
                }
                SNPList.push_back(s);
            }
        }
        getline(fin, line);
    }
    cout << "Inputted " << SNPList.size() << " SNPs\n";
    assert((int)SNPList.size() == tot);
}

//function parseForSift()
//Converts SNP List into a SIFT-compatible input
//Input: list containing all SNPs
//Output: List of SNPs that SIFT's genomic coordinate reader can interpret
void parseForSift(const vector<SNP> & SNPList)
{
    ofstream fout("./Output/Broad_WES_Data_SIFT_Input.vcf");
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(i == 0 || SNPList[i].genomicLoc != SNPList[i-1].genomicLoc)
        {
            fout << SNPList[i].chr.substr(3) << "," << SNPList[i].genomicLoc << ",1," << SNPList[i].refNT << "/" << SNPList[i].mutNT << "\n";
        }
    }
}

//function parseForPolyPhen()
//Converts SNP List into a PolyPhen-compatible input
//Input: list containing all SNPs
//Output: List of SNPs that PolyPhen can interpret
void parseForGVS(const vector<SNP> & SNPList)
{
    ofstream fout("./Output/Broad_WES_Data_GVS.txt");
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(i == 0 || SNPList[i].genomicLoc != SNPList[i-1].genomicLoc)
        {
            fout << SNPList[i].chr << "\t" << SNPList[i].genomicLoc << "\t" << SNPList[i].refNT << "\t" << SNPList[i].mutNT << "\n";
        }
    }
}

//function parseForPlink()
//Converts SNP List into a Plink-compatible input
//Input: list containing all SNPs
//Output: List of SNPs that Plink can interpret
void parseForPlink(const vector<SNP> & SNPList)
{
    ofstream fout_ped("./Output/Broad_WES_Data_Plink_Input.ped");
    ofstream fout_map("./Output/Broad_WES_Data_Plink_Input.map");
    string line;
    cout << "Parsing For Plink\n";
    for(int i = 0; i < (int)SNPList[0].calls.size(); i++)
    {//for all SNPs
        fout_ped << "F" << (i+1) << "\t" << (i+1) << "\t0\t0\t0\t0";
        for(int j = 0; j < (int)SNPList.size(); j++)
        {//create ped file
            if(j == 0 || SNPList[j].genomicLoc != SNPList[j-1].genomicLoc)
            {
                if(SNPList[j].calls[i].genotype == 1 && SNPList[j].calls[i].depth >= 8)
                    fout_ped << " " << SNPList[j].refNT << " " << SNPList[j].refNT;
                else if(SNPList[j].calls[i].genotype == 2 && SNPList[j].calls[i].depth >= 8)
                    fout_ped << " " << SNPList[j].mutNT << " " << SNPList[j].mutNT;
                else if(SNPList[j].calls[i].genotype == 3 && SNPList[j].calls[i].depth >= 15)
                    fout_ped << " " << SNPList[j].refNT << " " << SNPList[j].mutNT;
                else fout_ped << " 0 0";
                if(i == 0)
                {//create map file (using just the first patient)
                    fout_map << SNPList[j].chr.substr(3) << "\t";
                    if(SNPList[j].dbSNP == ".")
                        fout_map << "SNP" << j << "\t";
                    else fout_map << SNPList[j].dbSNP << "\t";
                    fout_map << "\t0\t" << SNPList[j].genomicLoc << "\n";
                }
            }
        }
        fout_ped << "\n";
    }
}


//function filterCoverages()
//Filters list of list of SNPs based on coverages per sample
//Input: list of list of SNPs, 1st dim = patients, 2nd dim = SNPs
//Output: List of list of SNPs in similar format as input w/ SNPs not passing coverage filter tossed out
void filterCoverages(vector<SNP> & SNPSamples)
{
    char c, p, m;
    cout << "Count low coverage SNPs? (y/n)\n";
    cin >> c;
    cout << "What level Plink? (n,1,3,5)\n";
    cin >> p;
    cout << "Enter mode of inheritance (m = homozyg, t = heterozyg):\n";
    cin >> m;
    for(int i = 0; i < (int)patientList.size(); i++)
    {//for each patient
        vector<SNP> newSNPList(0);
        for(int j = 0; j < (int)SNPSamples.size(); j++)
        {//for each SNP
            for(int k = 0; i == 0 && k < (int)SNPSamples[j].calls.size(); k++)
            {//check other genotypes
                if(((SNPSamples[j].calls[k].genotype == 2 && SNPSamples[j].calls[k].depth >= 8) ||
                    ((m == 't' || m == 'T') && SNPSamples[j].calls[k].genotype == 3 &&
                     SNPSamples[j].calls[k].depth >= 15)))
                {//if the other genotype is homozygous mutant, add to "shares" list
                    if(SNPSamples[j].shares == "")
                        SNPSamples[j].shares = patientList[k];
                    else SNPSamples[j].shares += "," + patientList[k];
                }
                else if((SNPSamples[j].calls[k].depth < 8 || (SNPSamples[j].calls[k].genotype == 3 && SNPSamples[j].calls[k].depth < 15)))
                {//else if the other genotype has low coverage, add to "unsure" list
                    if(SNPSamples[j].unsure == "")
                        SNPSamples[j].unsure = patientList[k];
                    else SNPSamples[j].unsure += "," + patientList[k];
                }
            }
            if(m == 'm' || m == 'M')
            {
                if(((SNPSamples[j].calls[i].genotype == 2 || SNPSamples[j].calls[i].genotype == 0 || (SNPSamples[j].calls[i].genotype == 1 && SNPSamples[j].calls[i].depth < 8) ||
                    (SNPSamples[j].calls[i].genotype == 3 && SNPSamples[j].calls[i].depth < 15)) && (c == 'y' || c == 'Y')) ||
                    ((c == 'n' || c == 'N') && SNPSamples[j].calls[i].genotype == 2 && SNPSamples[j].calls[i].depth >= 8))//if we're pretty sure the SNP is homozygous mutant
                {
                    newSNPList.push_back(SNPSamples[j]);
                }
            }
            else if(m == 'T' || m == 't')
            {
                if(((SNPSamples[j].calls[i].genotype == 2 || SNPSamples[j].calls[i].genotype == 0 || (SNPSamples[j].calls[i].genotype == 1 && SNPSamples[j].calls[i].depth < 8) ||
                    SNPSamples[j].calls[i].genotype == 3) && (c == 'y' || c == 'Y')) ||
                    ((c == 'n' || c == 'N') && ((SNPSamples[j].calls[i].genotype == 2 && SNPSamples[j].calls[i].depth >= 8)||
                     (SNPSamples[j].calls[i].genotype == 3 && SNPSamples[j].calls[i].depth >= 15))))//if we're pretty sure the SNP is homozygous mutant
                {
                    newSNPList.push_back(SNPSamples[j]);
                }
            }
        }
        cout << "Patient " << patientList[i] << " SNP size:" << newSNPList.size() << "\n";
        filterPlink(newSNPList, i, p);
        cout << "Patient " << patientList[i] << " SNP size:" << newSNPList.size() << "\n";
        ofstream fout(("./Output/" + patientList[i] + ".vcf").c_str());//open file/output header
        fout << "Chr\tPos\tNew_Loc\tID\tStrand\tRefNT\tAltNT\tUnique\tUnique w/ Low Cov.\tGT/Family Info\tCodon Change\tAA Change\tSide Chain Polarity/Charge\tAA Position\t"
        << "Gene Symbol\tGene Name\tTranscript ID\tSIFT\tSIFT P-val\tGERP\tPhastCons\tPolyPhen\tOMIM\tHGMD_SNP\tHGMD_Gene\tDomains\tdbSNP Genotype Counts\tAllele Counts";
        for(int j = 0; j < (int)patientList.size(); j++)
        {
            fout << "\t" << patientList[j];
        }
        fout << "\n";
        for(int j = 0; j < (int)newSNPList.size(); j++)
        {//for all SNPs, output it
            outputSNP(fout, newSNPList, j, i);
        }
    }
}

//function outputSNP()
//Outputs given SNP in proper format
//Input: List containing all SNPs, output file ID, SNP location in this list, patient ID
//Output: SNP corresponding to location outputted to the output file ID in proper format
void outputSNP(ofstream & fout, const vector<SNP> & SNPList, const int & i, const int & j)
{
    string s;
    fout << SNPList[i].chr << "\t" << SNPList[i].genomicLoc << "\t";
    if(i == 0 || SNPList[i].genomicLoc != SNPList[i-1].genomicLoc)
        fout << "Yes\t";
    else fout << "No\t";
    fout << SNPList[i].dbSNP << "\t" << (SNPList[i].strand ? '+' : '-') << "\t" << SNPList[i].refNT << "\t" << SNPList[i].mutNT << "\t";
    s = SNPList[i].shares;
    if(s == "" || s == patientList[j])
        fout << "Yes\t";
    else if((int)s.find(patientList[j]) == -1)
        fout << s << "\t";
    else if((int)s.find(patientList[j]) == 0)//first one
        fout << s.substr(patientList[j].size() + 1) << "\t";
    else
    {
        fout << s.substr(0, s.find(patientList[j]) - 1);
        if((int)s.find(patientList[j] + ",") != -1)
        {//not the last one
            fout << s.substr(s.find(patientList[j])+patientList[j].size());
        }
        fout << "\t";
    }
    s = SNPList[i].unsure;
    if(s == "" || s == patientList[j])
        fout << "Yes\t";
    else if((int)s.find(patientList[j]) == -1)
        fout << s << "\t";
    else if((int)s.find(patientList[j]) == 0)
        fout << s.substr(patientList[j].size() + 1) << "\t";
    else
    {
        fout << s.substr(0, s.find(patientList[j]) - 1);
        if((int)s.find(patientList[j] + ",") != -1)
        {//not the last one
            fout << s.substr(s.find(patientList[j])+patientList[j].size());
        }
        fout << "\t";
    }
    if(SNPList[i].calls[j].genotype == 0)//output only patient j's genotype genotypes
        fout << "./.\t";
    else if(SNPList[i].calls[j].genotype == 1)
        fout << "0/0:" << SNPList[i].calls[j].depth << ":" << SNPList[i].calls[j].quality << "\t";
    else if(SNPList[i].calls[j].genotype == 2)
        fout << "1/1:" << SNPList[i].calls[j].depth << ":" << SNPList[i].calls[j].quality << "\t";
    else if(SNPList[i].calls[j].genotype == 3)
        fout << "0/1:" << SNPList[i].calls[j].depth << ":" << SNPList[i].calls[j].quality << "\t";
    if(SNPList[i].is_splice && SNPList[i].refAA == SNPList[i].mutAA)
    {//if it's a splice variant, don't bother wrt codon changes, etc.
        fout << "splice\tsplice\tN/A\tN/A\t";
    }
    else if(SNPList[i].is_utr)
    {//if it's a splice variant, don't bother wrt codon changes, etc.
        fout << "utr\tutr\tN/A\tN/A\t";
    }
    else
    {//else output cDNA & protein info
        string codon = SNPList[i].codon;
        if(codon == "")
        {//if the codon is not defined, output N/A (the case if Broad didn't call it as coding, but GVS did)
            fout << "N/A\t";
        }
        else
        {//find the mutated BP
            char c = 'N';
            if(SNPList[i].strand)
                c = SNPList[i].mutNT;
            else
            {
                if(SNPList[i].mutNT == 'G')
                    c = 'C';
                else if(SNPList[i].mutNT == 'C')
                    c = 'G';
                else if(SNPList[i].mutNT == 'A')
                    c = 'T';
                else if(SNPList[i].mutNT == 'T')
                    c = 'A';
            }
            codon = codon.substr(0, SNPList[i].codonPos - SNPList[i].codonStart) + c + codon.substr(SNPList[i].codonPos - SNPList[i].codonStart + 1);
            if(SNPList[i].is_splice)
                fout << "splice/";
            fout << "c.(" << SNPList[i].protPos * 3 - 2 << "-" << SNPList[i].protPos * 3 << ")" << SNPList[i].codon << ">" << codon << "\t";
        }
        if(SNPList[i].is_splice)
            fout << "splice/";
        fout << "p." << SNPList[i].refAA << SNPList[i].protPos << SNPList[i].mutAA << "\t";
        fout << classifyAA(SNPList[i].refAA) << ">" << classifyAA(SNPList[i].mutAA) << "\t";
        fout << SNPList[i].protPos << "/" << SNPList[i].totAA << "\t";
    }
    fout << SNPList[i].geneName << "\t" << SNPList[i].fullName << "\t" << SNPList[i].transcriptID << "\t" << SNPList[i].SIFT << "\t" << SNPList[i].SIFTProb << "\t"
     << SNPList[i].GERP << "\t" << SNPList[i].phastCons << "\t" << SNPList[i].polyPhen << "\t" << SNPList[i].OMIM << "\t" << SNPList[i].HGMD_SNP << "\t"
     << SNPList[i].HGMD_Gene << "\t" << SNPList[i].Domains << "\t" << SNPList[i].GTCounts << "\t" << SNPList[i].AlleleCounts << "\t" << SNPList[i].full_shares << "\n";
}

//function filterPlink()
//Filters list of list of SNPs based on if it's in a Plink runs of homozygosity region
//Input: List of list of SNPs, 1st dim = patients, 2nd dim = SNPs
//Output: Same list as above, with SNPs not in regions of homozygosity filtered out
void filterPlink(vector<SNP> & SNPList, const int & sampleID, const char & Plink_level)
{
    ifstream fin;
    if(Plink_level == '1')
        fin.open("./plink.hom");
    else if(Plink_level == '3')
        fin.open("./plink_3mb.hom");
    else if(Plink_level == '5')
        fin.open("./plink_5mb.hom");
    else return;
    cout << "Filtering based on Plink homozygosity regions\n";
    string line;
    getline(fin, line);
    getline(fin, line);
    vector<SNP> newSNPList(0);
    int i = 0, start, end;
    string chr;
    while(i <= sampleID && !fin.eof())
    {//parse chromosome, starting pos, ending pos from each line
        line = line.substr(line.find('F') + 1);
        i = atoi(line.substr(0, line.find(' ')).c_str()) - 1;
        if(i == sampleID)
        {
            line = line.substr(line.find('.') + 4);
            line = line.substr(line.find_first_of("0123456789"));
            chr = "chr" + line.substr(0, line.find(' '));
            if(chr == "chr23")
                chr = "chrX";
            else if(chr == "chr24")
                chr = "chrY";
            line = line.substr(line.find_first_of("rS") + 3);
            line = line.substr(line.find_first_of("rS") + 2);
            line = line.substr(line.find(' '));
            line = line.substr(line.find_first_of("0123456789"));
            start = atoi(line.substr(0, line.find(' ')).c_str());
            line = line.substr(line.find(' '));
            line = line.substr(line.find_first_of("0123456789"));
            end = atoi(line.substr(0, line.find(' ')).c_str());
            for(int j = 0; j < (int)SNPList.size(); j++)
            {
                if(SNPList[j].chr == chr && SNPList[j].genomicLoc >= start && SNPList[j].genomicLoc <= end)
                {//if the SNP is in a region of homozygosity, put it back on the list
                    newSNPList.push_back(SNPList[j]);
                }
            }
        }
        getline(fin, line);
    }
    SNPList = newSNPList;
}

void filterPlink(vector<vector<indel> > & indelList)
{
    ifstream fin("./plink.hom");
    cout << "Filtering based on Plink homozygosity regions\n";
    string line;
    getline(fin, line);
    getline(fin, line);
    vector<vector<indel> > newIndelList(indelList.size());
    int i, start, end;
    string chr;
    while(!fin.eof())
    {//parse chromosome, starting pos, ending pos from each line
        line = line.substr(line.find('F') + 1);
        i = atoi(line.substr(0, line.find(' ')).c_str()) - 1;
        line = line.substr(line.find('.') + 4);
        line = line.substr(line.find_first_of("0123456789"));
        chr = "chr" + line.substr(0, line.find(' '));
        if(chr == "chr23")
            chr = "chrX";
        else if(chr == "chr24")
            chr = "chrY";
        line = line.substr(line.find_first_of("rS") + 3);
        line = line.substr(line.find_first_of("rS") + 2);
        line = line.substr(line.find(' '));
        line = line.substr(line.find_first_of("0123456789"));
        start = atoi(line.substr(0, line.find(' ')).c_str());
        line = line.substr(line.find(' '));
        line = line.substr(line.find_first_of("0123456789"));
        end = atoi(line.substr(0, line.find(' ')).c_str());
        for(int j = 0; j < (int)indelList[i].size(); j++)
        {
            if(indelList[i][j].chr == chr && ((indelList[i][j].start >= start && indelList[i][j].start <= end) || ((int)(indelList[i][j].start + indelList[i][j].ref.size()) >= start && (int)(indelList[i][j].start + indelList[i][j].ref.size()) <= end)))
            {//if the SNP is in a region of homozygosity, put it back on the list
                newIndelList[i].push_back(indelList[i][j]);
            }
        }
        getline(fin, line);
    }
    indelList = newIndelList;
}


//function filterGERP()
//Filters list of list of SNPs based on GERP conservation score
//Input: List of list of SNPs, 1st dim = patients, 2nd dim = SNPs
//Output: Same list as above, with SNPs with negative GERP scores filtered out
void filterGERP(vector<SNP> & SNPList)
{
    ifstream fin("./SeattleSeqAnnotation131.Broad_WES_Data_Polyphen.194278187585.txt");
    string line;
    getline(fin, line);
    cout << "Filtering GERP\n";
    vector<SNP> newSNPList;
    int distToSplice;
    string prev_loc;
    string transcript;
    int iprime;
    char ncSNPs;
    double GERP;
    cout << "Include splice/UTR SNPs (y/n)?\n";
    cin >> ncSNPs;
    cout << "Enter GERP threshold (for ns coding SNPs):\n";
    cin >> GERP;
    for(int i = 0; i < (int)SNPList.size(); i++)
    {
        if(i > 0 && SNPList[i].genomicLoc == SNPList[i-1].genomicLoc)
        {
            SNPList[i].dbSNP = SNPList[i-1].dbSNP;
            SNPList[i].phastCons = SNPList[i-1].phastCons;
            SNPList[i].polyPhen = SNPList[i-1].polyPhen;
            SNPList[i].GERP = SNPList[i-1].GERP;
            if(newSNPList.size() > 0 && SNPList[i].genomicLoc == newSNPList[newSNPList.size() - 1].genomicLoc && ((SNPList[i].mutAA != SNPList[i].refAA) || SNPList[i].is_splice || SNPList[i].mutAA == '*' || SNPList[i].is_utr))//GERP + Nonsyn
                newSNPList.push_back(SNPList[i]);
            continue;
        }
        getline(fin, line);
        line = line.substr(line.find('\t') + 1);//inDBSNP
        line = line.substr(line.find('\t') + 1);//chromosome
        if(line.substr(0, line.find('\t')) == prev_loc)
        {//if the location doesn't change, update total amino acids and go on
            i--;
            if(i < (int)SNPList.size() - 1 && SNPList[i+1].genomicLoc == SNPList[i].genomicLoc)
            {
                iprime = i+1;
                line = line.substr(line.find('\t') + 1);//position
                assert(line[0] == SNPList[iprime].refNT);//make sure we're referring to the same SNP
                line = line.substr(line.find('\t') + 1);//refBase
                assert(line[0] == SNPList[iprime].mutNT);
                line = line.substr(line.find('\t') + 1);//sampleGT
                line = line.substr(line.find('\t') + 1);//sampleAlleles
                line = line.substr(line.find('\t') + 1);//dbSNP alleles
                transcript = line.substr(0, line.find('\t'));
                while(iprime < (int)SNPList.size() && SNPList[iprime].genomicLoc == SNPList[i].genomicLoc && SNPList[iprime].transcriptID != transcript)
                {
                    iprime++;
                }
                if(iprime < (int)SNPList.size() && SNPList[iprime].transcriptID == transcript)
                {
                    SNPList[iprime].totAA = line.substr(line.find_last_of("\t")).size() - 1;//1 = * at the end
                }
            }
            else
            {
                iprime = newSNPList.size() - 1;
                line = line.substr(line.find('\t') + 1);//position
                line = line.substr(line.find('\t') + 1);//refBase
                line = line.substr(line.find('\t') + 1);//sampleGT
                line = line.substr(line.find('\t') + 1);//sampleAlleles
                line = line.substr(line.find('\t') + 1);//dbSNP alleles
                transcript = line.substr(0, line.find('\t'));
                while(iprime >= 0 && newSNPList[iprime].genomicLoc == SNPList[i].genomicLoc)
                {
                    if(newSNPList[iprime].transcriptID == transcript)
                    {
                        newSNPList[iprime].totAA = line.substr(line.find_last_of("\t")).size() - 1;//1 = tab, 1 = * at the end
                    }
                    iprime--;
                }
            }
            continue;
        }
        assert(atoi(line.substr(0, line.find('\t')).c_str()) == SNPList[i].genomicLoc);
        prev_loc = line.substr(0, line.find('\t'));
        line = line.substr(line.find('\t') + 1);//position
        assert(line[0] == SNPList[i].refNT);//make sure we're referring to the same SNP
        line = line.substr(line.find('\t') + 1);//refBase
        assert(line[0] == SNPList[i].mutNT);
        line = line.substr(line.find('\t') + 1);//sampleGT
        line = line.substr(line.find('\t') + 1);//sampleAlleles
        line = line.substr(line.find('\t') + 1);//dbSNP alleles
        transcript = line.substr(0, line.find('\t'));
        iprime = i;
        while(iprime < (int)SNPList.size() && SNPList[iprime].genomicLoc == SNPList[i].genomicLoc && SNPList[iprime].transcriptID != transcript)
        {
            iprime++;
        }
        if(SNPList[iprime].transcriptID == transcript)
        {
            SNPList[iprime].totAA = line.substr(line.find_last_of("\t")).size() - 1;//1 = tab, 1 = * at the end
        }
        line = line.substr(line.find('\t') + 1);//transcript acc #
        line = line.substr(0, line.find_last_of("\t"));//prot seq
        line = line.substr(0, line.find_last_of("\t"));//miRNAs
        distToSplice = atoi(line.substr(line.find_last_of("\t") + 1).c_str());
        if(((int)line.substr(0, line.find('\t')).find("splice") != -1 || (((int)line.substr(0, line.find('\t')).find("intron") != -1 && SNPList[i].spliceDist == 0 && distToSplice <= 4) ||
            ((int)line.substr(0, line.find('\t')).find("missense") == -1 && (int)line.substr(0, line.find('\t')).find("nonsense") == -1 && SNPList[i].spliceDist == 0 && distToSplice <= 2))))
        {//if it should be a splice, change it
            SNPList[i].is_splice = true;
        }
        line = line.substr(line.find('\t') + 1);//GVS Function
        line = line.substr(line.find('\t') + 1);//dbSNP Function
        if(line.substr(0, line.find('\t')) != "0")
        {//if we found a dbSNP ID, attach it
            if(SNPList[i].dbSNP == ".")
                 SNPList[i].dbSNP = "rs" + line.substr(0, line.find('\t'));
        }
        line = line.substr(line.find('\t') + 1);//rsID
        line = line.substr(line.find(',') + 1);//ref AA
        line = line.substr(line.find('\t') + 1);//mut AA
        line = line.substr(line.find('\t') + 1);//protein position
        SNPList[i].polyPhen = line.substr(0, line.find('\t'));//input polyPhen
        line = line.substr(line.find('\t') + 1);//polyPhen
        SNPList[i].phastCons = atof(line.substr(0, line.find('\t')).c_str());
        line = line.substr(line.find('\t') + 1);//scorePhastCons
        if((atof(line.substr(0, line.find('\t') + 1).c_str()) < GERP || SNPList[i].mutAA == SNPList[i].refAA) && SNPList[i].mutAA != '*' && (!SNPList[i].is_splice && !SNPList[i].is_utr))//GERP + Nonsyn
        //if((SNPList[i].mutAA == SNPList[i].refAA) && !SNPList[i].is_splice && SNPList[i].mutAA != '*' && !SNPList[i].is_utr)//GERP + Nonsyn
        {//if GERP is negative or it's a synonymous coding mutation, don't skip
            continue;
        }
        SNPList[i].GERP = atof(line.substr(0, line.find('\t')).c_str());//else add GERP info
        newSNPList.push_back(SNPList[i]);//reinsert into list
    }
    SNPList = newSNPList;
    cout << "SNP List Size: " << SNPList.size() << "\n";
}

//function classifyAA()
//Outputs the charge/polarity info for a given AA
//Input: 1 letter code for amino acids
//Output: Charge/polarity corresopnding to this code
string classifyAA(const char & AA)
{
    switch(AA)
    {
        case 'A': return "np-neu";
        case 'C': return "np-neu";
        case 'D': return "p-neg";
        case 'E': return "p-neg";
        case 'F': return "np-neu";
        case 'G': return "np-neu";
        case 'H': return "p-neu";
        case 'I': return "np-neu";
        case 'K': return "p-pos";
        case 'L': return "np-neu";
        case 'M': return "np-neu";
        case 'N': return "p-neu";
        case 'P': return "np-neu";
        case 'Q': return "p-neu";
        case 'R': return "p-pos";
        case 'S': return "p-neu";
        case 'T': return "p-neu";
        case 'V': return "np-neu";
        case 'W': return "np-neu";
        case 'Y': return "p-neu";
        case '*': return "stop";
        default: return "N/A";
    };
}

//function ThreeToOne()
//Converts GVS's 3 letter AA code into 1 letter AA code
//Input: 3 letter AA code
//Output: 1 letter AA code
char ThreeToOne(const string & AA)
{
    if(AA == "ALA")
        return 'A';
    if(AA == "CYS")
        return 'C';
    if(AA == "ASP")
        return 'D';
    if(AA == "GLU")
        return 'E';
    if(AA == "PHE")
        return 'F';
    if(AA == "GLY")
        return 'G';
    if(AA == "HIS")
        return 'H';
    if(AA == "ILE")
        return 'I';
    if(AA == "LYS")
        return 'K';
    if(AA == "LEU")
        return 'L';
    if(AA == "MET")
        return 'M';
    if(AA == "ASN")
        return 'N';
    if(AA == "PRO")
        return 'P';
    if(AA == "GLN")
        return 'Q';
    if(AA == "ARG")
        return 'R';
    if(AA == "SER")
        return 'S';
    if(AA == "THR")
        return 'T';
    if(AA == "VAL")
        return 'V';
    if(AA == "TRP")
        return 'W';
    if(AA == "TYR")
        return 'Y';
    if(AA == "stop")
        return '*';
    cout << "Unknown string " << AA << "\n";
    return 'X';
}
