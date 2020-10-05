/************************************************************************
 * CASSI, version 2.51
 * Copyright 2012-2017,
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of CASSI, the SNP interaction analysis program.
 *
 * CASSI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CASSI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CASSI.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Data.h
    \brief This file contains classes for manipulating SNP data.
    
*/

#ifndef __DATA
#define __DATA

#include <set>
#include <map>
#include <list>
#include <stack>
#include <fstream>
#include <string>

using namespace std;


//! Stores SNP data for all subjects for a given SNP.
struct SNPData
{
	
	list<unsigned char> noMinorAllelesAllSubjects; //ordered list of the number of minor alleles for the subjects	
	
	SNPData() : noMinorAllelesAllSubjects() {};
	
	~SNPData()
	{


	};
	
};

//! Joint Genotype Counts
struct QuantitiveTraits
{
	map<unsigned int, double> values; //subject no, quantitive trait

	QuantitiveTraits() : values() {};

	void setupQuantitiveTraits(string & filename, double & missingQTValue);

};

//! Joint Genotype Counts
struct JointGenotypeCounts
{
	unsigned int total;
	unsigned int counts[3][3];

	JointGenotypeCounts()
	{
		resetCounts();	
	};

	void addCount(const unsigned int & count1, const unsigned int & count2);
	void resetCounts();
};

//! Covariate data for all subjects
struct CovariateData
{
	list<list<double> > covariateDataAllSubjects; //ordered individuals, each with ordered list of covariates
	list<unsigned char> caseControls; 

	CovariateData() {};
	CovariateData(string & covariateFilename, string & covariates, string & fname, double & missingQTValue);	

	~CovariateData()
	{

	};

};

//! General class for a SNP window
class SNPWindow
{
protected:

	unsigned int totalNoSubjects;
	unsigned int startSNPNo;
	string filename;
	unsigned int bitCount;
	int one;
	int aBit;	
	char buffer[1];
	ifstream readSNPData;

public:

	SNPWindow(const unsigned int & ts, unsigned int & ssn, string & fname) : totalNoSubjects(ts), startSNPNo(ssn), filename(fname), bitCount(9), one('\1') {};
	

	virtual ~SNPWindow()
	{		
		if(readSNPData.is_open()) readSNPData.close();	
	};

	void advanceSNPData();
	void advanceToFirstWindow();
	unsigned char getNextNoOfMinorAlleles();
	void startNewByte() {bitCount = 9;};
	void reopenBinaryFile();
	void openBinaryFileFirst();
	void closeBinaryFile() {readSNPData.close();};

	virtual void moveToNextSNP() {};
	virtual void startWindowAtStart() {};
	virtual SNPData * getCurrentSNP() {return 0;};
};

//! Stores SNP data for the SNP window in question.
class SNPWindowReadFromFile : public SNPWindow
{
private:

	SNPData * snp;

public:

	SNPWindowReadFromFile(const unsigned int & ts, unsigned int & startSNPNo, string & fname, const bool & readFirstSNP);
	

	~SNPWindowReadFromFile()
	{	
		closeBinaryFile();
		delete snp;		
	};

	void moveToNextSNP();
	void startWindowAtStart();
	SNPData * getCurrentSNP() {return snp;};
};

//! SNP window that stores all the data
class SNPWindowStoreAllData : public SNPWindow
{
private:

	list<SNPData *> allSNPData;
	list<SNPData *>::const_iterator currentSNP;

public:

	SNPWindowStoreAllData(const unsigned int & ts, unsigned int & startSNPNo, string & fname, unsigned int & endSNPNo);
	
	~SNPWindowStoreAllData()
	{		
		for(list<SNPData *>::iterator s = allSNPData.begin(); s != allSNPData.end(); ++s)
		{
			delete *s;
		};
	};

	void moveToNextSNP() {currentSNP++;};
	void startWindowAtStart() {currentSNP = allSNPData.begin();};
	SNPData * getCurrentSNP() {return *currentSNP;};
};


//! SNP window that stores all the data, but in condensed binary format
class SNPWindowStoreAllDataBinary : public SNPWindow
{
private:

	list<char> allSNPData;	
	list<char>::iterator currentByteIt;	

	//store all data somehow
	SNPData * currentSNPData; //update this from all data when necessary, need some sort
	bool firstByte;

public:

	SNPWindowStoreAllDataBinary(const unsigned int & ts, unsigned int & startSNPNo, string & fname, unsigned int & endSNPNo);
	
	~SNPWindowStoreAllDataBinary()
	{	
		delete currentSNPData;
		//for(list<SNPData *>::iterator s = allSNPData.begin(); s != allSNPData.end(); ++s)
		//{
		//	delete *s;
		//};
	};

	void moveToNextSNP() {updateCurrentSNPData();};
	void startWindowAtStart() {currentByteIt = allSNPData.begin(); firstByte = true; updateCurrentSNPData();};
	void updateCurrentSNPData();
	SNPData * getCurrentSNP() {return currentSNPData;};
};

//! Stores info about a SNP.
struct SNPInfo
{
	unsigned int chromosome;
	string name;
	unsigned int bp;

	SNPInfo(unsigned int c, string n, unsigned int b) : chromosome(c), name(n), bp(b) {};

	~SNPInfo() {};
};

//! Names of SNPs
class DescriptionOfSNPs
{
private:
	map<unsigned int, SNPInfo *> snpInfos1; //snpNo, SNPInfo, base-pair positions of all SNPs in order with SNP names
	map<unsigned int, SNPInfo *> snpInfos2;
	bool oneFile;

	map<unsigned int, SNPInfo *>::const_iterator snp1Info;
	map<unsigned int, SNPInfo *>::const_iterator snp2Info;

public:

	DescriptionOfSNPs(string & filename, string & filename2)
	{
		if(filename2 == "") oneFile = true;
		else oneFile = false;

		setUpSNPDesciptionData(filename, false);
		if(filename2 != "") setUpSNPDesciptionData(filename2, true);	
	};

	//methods regarding the SNP descriptions
	void setUpSNPDesciptionData(string & filename, const bool & secondFile);
	void setSNP1DescIterator(unsigned int & snpNo1);
	void setSNP2DescIterator(unsigned int & snpNo2);
	void advanceIterator1() {++snp1Info;};
	void advanceIterator2() {++snp2Info;};
	SNPInfo * getSNPInfo1() {return snp1Info->second;};
	SNPInfo * getSNPInfo2() {return snp2Info->second;};

	unsigned int getNoSNPs(const bool & second) const
	{
		if(second) return snpInfos2.size();
		else return snpInfos1.size();
	};	
};



#endif
