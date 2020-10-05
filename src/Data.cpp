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


/*! \file Data.cpp
    \brief This file contains the source for manipulating SNP data.
    
*/

#include <string>
#include <cstring>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Data.h"
 
//! Adds a count to the joint genotypes between two SNPs.
void JointGenotypeCounts::addCount(const unsigned int & count1, const unsigned int & count2)
{
	if(count1 != 3 && count2 != 3)
	{
		counts[count1][count2]++;
		++total;
	};
};

//! Resets all counts to zero.
void JointGenotypeCounts::resetCounts()
{
	total = 0;
	counts[0][0] = 0;
	counts[0][1] = 0;
	counts[0][2] = 0;
	counts[1][0] = 0;
	counts[1][1] = 0;
	counts[1][2] = 0;
	counts[2][0] = 0;
	counts[2][1] = 0;
	counts[2][2] = 0;	
};

void QuantitiveTraits::setupQuantitiveTraits(string & filename, double & missingQTValue)
{
	//try and find the family file and read in data
	unsigned int length = filename.length();
	string famFilename = filename.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename); outErr("!\n");
		exit(0);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	double phenoType;
	unsigned int noCases = 0;
	unsigned int subjectNo = 1;
	unsigned int noMissingPhenos = 0;

	//loop thro' subjects and store phenotypes
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;
	
		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			
			values[subjectNo] = phenoType;
			subjectNo++;
			if(phenoType == missingQTValue)
			{
				++noMissingPhenos; 				
			};			
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	readFamilyFile.close();
};

//! Sets up SNP data from the .bim file
void DescriptionOfSNPs::setUpSNPDesciptionData(string & filename, const bool & secondFile)
{

	//try and find the binary map file, .bim, and read in data
	unsigned int length = filename.length();
	string mapFilename = filename.substr(0, length-4) + ".bim";

	ifstream readMapFile;
	readMapFile.open(mapFilename.c_str());
	if(!readMapFile.is_open())
	{
		outErr("Cannot read bim file: "); outErr(mapFilename); outErr("!\n");
		exit(0);
	};

	string snpIdentifier, geneticDistance;
	string prevSnpIdentifier = "";
	unsigned int chromosome, basePairPosition;	
	string alleleName1, alleleName2;
	unsigned int snpID = 1;

	//loop thro' subjects and store the cases
	do{
		
		readMapFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >> alleleName1 >> alleleName2;
		
		if(snpIdentifier != prevSnpIdentifier)
		{				
				if(!secondFile) snpInfos1[snpID] = new SNPInfo(chromosome, snpIdentifier, basePairPosition);
				else snpInfos2[snpID] = new SNPInfo(chromosome, snpIdentifier, basePairPosition);
				snpID++;
		};
		
		prevSnpIdentifier = snpIdentifier;				
	}while(!readMapFile.eof());

	readMapFile.close();

	if(!secondFile) out("Data Summary Statistics:\n");
	if(!secondFile && oneFile) {out("Number of SNPs: "); out(snpInfos1.size()); out("\n");}
	else if(!secondFile) {out("Number of SNPs in "); out(filename); out(": "); out(snpInfos1.size()); out("\n");}
	else {out("Number of SNPs in "); out(filename); out(": "); out(snpInfos2.size() ); out("\n");};
	
};


//! Sets up the SNP  desc iterators.
void DescriptionOfSNPs::setSNP1DescIterator(unsigned int & snpNo1)
{	
	snp1Info = snpInfos1.find(snpNo1);
	if(snp1Info == snpInfos1.end())
	{
		outErr("Problem finding SNP no. "); outErr(snpNo1); outErr(" on first SNP window!");
		exit(0);
	};
};

//! Sets up the SNP 2 desc iterators.
void DescriptionOfSNPs::setSNP2DescIterator(unsigned int & snpNo2)
{	
	if(oneFile)
	{
		snp2Info = snpInfos1.find(snpNo2);
		if(snp2Info == snpInfos1.end())
		{
			outErr("Problem finding SNP no. "); outErr(snpNo2); outErr(" on second SNP window!");
			exit(0);
		};
	}
	else
	{
		snp2Info = snpInfos2.find(snpNo2);
		if(snp2Info == snpInfos2.end())
		{
			outErr("Problem finding SNP no. "); outErr(snpNo2); outErr(" on second SNP window!");
			exit(0);
		};
	};
};

void SNPWindow::reopenBinaryFile()
{
	readSNPData.close();
	readSNPData.open(filename.c_str(), ios::binary);
	
	//get past the three special bytes at before the SNP data
	char buffer[3];
	readSNPData.read(buffer, 3);
};

//! Open the binary file for the first time
void SNPWindow::openBinaryFileFirst()
{
	//try and find the binary pedigree file, .bed, and read in data for the first window
	readSNPData.open(filename.c_str(), ios::binary);
	
	if(!readSNPData.is_open())
	{
		outErr("Cannot read binary pedigree file: "); outErr(filename); outErr("!\n");
		exit(0);
	};

	char buffer[3];
	readSNPData.read(buffer, 3);

	//check the plink magic numbers for the file type
	//3rd number indicates format of genotype data, 1 => subjects x SNPs, 0 => SNPs x subjects
	unsigned int magicNumber1 = buffer[0];
	unsigned int magicNumber2 = buffer[1];

	if(magicNumber1 != 108 || magicNumber2 != 27)
	{
		outErr("Detected an old version .bed file!\n");
		outErr("Please use PLINK to update the .bed file.\n");
			
		readSNPData.close();		
		exit(0);
	};

	//determine binary file type
	unsigned int mode = buffer[2];
	if(mode == 0)
	{
		outErr("The binary pedigree file must be in SNP-major mode!\n");
		outErr("Please use PLINK to update the .bed file.\n");
			
		readSNPData.close();		
		exit(0);
	};
};

unsigned char SNPWindow::getNextNoOfMinorAlleles()
{
	int allele1, allele2;
	unsigned char noMinorAlleles = 0;

	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readSNPData.read(buffer, 1);
		if(readSNPData.eof())
		{			
			outErr("Error: reached end of binary SNP file!\n");
			exit(0);
		};
			
		aBit = buffer[0];
			
		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	//if genotype is encoded 1/0 then the genotype is missing so do not add it
	if(allele1 == 1 && allele2 == 1)
	{	
		noMinorAlleles = 0;
	}
	else if(allele1 == 0 && allele2 == 1)
	{	
		noMinorAlleles = 1;
	}
	else if(allele1 == 0 && allele2 == 0)
	{	
		noMinorAlleles = 2;
	}
	else
		noMinorAlleles = 3; //denotes missing genotype

	return noMinorAlleles;

};


//! Creates the SNP window.
SNPWindowReadFromFile::SNPWindowReadFromFile(const unsigned int & ts, unsigned int & ssn, string & fname, const bool & readFirstSNP) : SNPWindow(ts, ssn, fname)
{
	//create SNP data object
	snp = new SNPData();

	//Open the binary file for the first time
	openBinaryFileFirst();

	//put data in for window 1, but not window 2, since window 2 will be reopened
	if(readFirstSNP)
	{
		advanceToFirstWindow();

		//setup the data for the first SNP
		for(unsigned int subjectNo = 1; subjectNo <= ts; ++subjectNo)
		{
			snp->noMinorAllelesAllSubjects.push_back(getNextNoOfMinorAlleles());		
		};
	}
	else
	{
		//setup the data for the first SNP with no counts to overwrite later
		for(unsigned int subjectNo = 1; subjectNo <= ts; ++subjectNo)
		{
			snp->noMinorAllelesAllSubjects.push_back(0);		
		};
	};
};

SNPWindowStoreAllData::SNPWindowStoreAllData(const unsigned int & ts, unsigned int & ssn, string & fname, unsigned int & endSNPNo) : SNPWindow(ts, ssn, fname)
{
	//Open the binary file for the first time
	openBinaryFileFirst();

	//advance thro' the snp data to the first SNP in the window
	advanceToFirstWindow();

	//now loop thro' the remaining SNPs and store data for each SNP until the last SNP
	SNPData * someSNPData;
	unsigned int subjectNo;
	unsigned int snpNo = ssn;
	unsigned int m;

	while(snpNo <= endSNPNo)
	{

		startNewByte();
		someSNPData = new SNPData();

		for(subjectNo = 1; subjectNo <= ts; ++subjectNo)
		{
			m = getNextNoOfMinorAlleles();
			someSNPData->noMinorAllelesAllSubjects.push_back(m);		
		};

		allSNPData.push_back(someSNPData);
		++snpNo;
	};

	currentSNP = allSNPData.begin();

	closeBinaryFile();
};
	
SNPWindowStoreAllDataBinary::SNPWindowStoreAllDataBinary(const unsigned int & ts, unsigned int & ssn, string & fname, unsigned int & endSNPNo) : SNPWindow(ts, ssn, fname)
{
	//Open the binary file for the first time
	openBinaryFileFirst();

	//advance thro' the snp data to the first SNP in the window
	advanceToFirstWindow();

	//now loop thro' the remaining SNPs and store data for each SNP until the last SNP
	unsigned int byteNo;
	char buffer[1];
	unsigned int snpNo = ssn;
	unsigned int noBytesPerSNP = (unsigned int)((double)((double)(ts)/4.0 + 0.8)); //round up if not multiple of 4

	while(snpNo <= endSNPNo)
	{
		for(byteNo = 1; byteNo <= noBytesPerSNP; ++byteNo)
		{
			readSNPData.read(buffer, 1);
			if(readSNPData.eof())
			{			
				outErr("Error: reached end of binary SNP file!\n");
				exit(0);
			};
			allSNPData.push_back(buffer[0]);				
		};
		
		++snpNo;
	};
	
	currentByteIt = allSNPData.begin();
	currentSNPData = new SNPData();

	//setup SNPData with correct number of entries, the number of subjects
	for(unsigned int i = 1; i <= ts; ++i) currentSNPData->noMinorAllelesAllSubjects.push_back(0);

};

void SNPWindowStoreAllDataBinary::updateCurrentSNPData()
{
	bitCount = 9;
	int allele1, allele2;
	unsigned char noMinorAlleles = 0;

	for(list<unsigned char>::iterator sub = currentSNPData->noMinorAllelesAllSubjects.begin(); sub != currentSNPData->noMinorAllelesAllSubjects.end(); ++sub)
	{
		//read in the next piece of data
		if(bitCount == 9)
		{			
			if(firstByte)
			{
				firstByte = false;
			}
			else
			{
				++currentByteIt;

				if(currentByteIt == allSNPData.end())
				{
					outErr("Problem with reading SNP data!\n");						
					exit(0);
				};
			};

			aBit = *currentByteIt;		

			bitCount = 1;
		};		

		allele1 = aBit & one; //read the least significant bit				
		aBit = aBit >> 1; //shift bits to the right
		allele2 = aBit & one; //read the new least significant bit				
		aBit = aBit >> 1; //shift bits to the right for next time

		bitCount += 2;	

		//if genotype is encoded 1/0 then the genotype is missing 
		if(allele1 == 1 && allele2 == 1)
		{	
			noMinorAlleles = 0;
		}
		else if(allele1 == 0 && allele2 == 1)
		{	
			noMinorAlleles = 1;
		}
		else if(allele1 == 0 && allele2 == 0)
		{	
			noMinorAlleles = 2;
		}
		else
			noMinorAlleles = 3; //denotes missing genotype

		*sub = noMinorAlleles;
	};
	
};

//! Sets the SNP to the start of the window
void SNPWindowReadFromFile::startWindowAtStart()
{
	reopenBinaryFile();

	advanceToFirstWindow();

	moveToNextSNP();//read in the data for the first SNP in the window
};

//! Adds the SNP data of the next SNP to the SNP data object
void SNPWindowReadFromFile::moveToNextSNP()
{
	//ensure a new byte is read
	startNewByte();

	for(list<unsigned char>::iterator i = snp->noMinorAllelesAllSubjects.begin(); i != snp->noMinorAllelesAllSubjects.end(); ++i)
	{
		*i = getNextNoOfMinorAlleles();
	};

};

//! Moves through SNP data by one SNP.
void SNPWindow::advanceSNPData()
{
	//a new byte is started after each SNP, 4 subject genotypes per byte,
	// so the no. of bytes is rounded up when divided by 4
	unsigned int bufferSize;
	if(totalNoSubjects%4 == 0) bufferSize = totalNoSubjects/4;
	else bufferSize = (unsigned int)(totalNoSubjects/4) + 1;

	char buffer[1];
	for(unsigned int i = 1; i <= bufferSize; ++i)
	{
		readSNPData.read(buffer, 1);
	};

};

//! Moves SNP window onto the chosen first SNP.
void SNPWindow::advanceToFirstWindow()
{
	unsigned int snpCount = 1;

	do{
		if(snpCount == startSNPNo) break;

		//read SNP data in for the previous SNP that will not be used
		advanceSNPData();

		snpCount++;
	}while(!readSNPData.eof());

};


//! Constucts covariate data and reads in data from file.
CovariateData::CovariateData(string & covariateFilename, string & covariates, string & fname, double & missingQTValue)
{
	//try and find the family file and read in data to create map
	map<string, unsigned int> idToIndivNo; //id name, position of indiv in file
	map<unsigned int, string> indivNoToId; //position of indiv in file, id name

	unsigned int length = fname.length();
	string famFilename = fname.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename); outErr("!\n");
		exit(0);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	double phenoType;	
	unsigned int subjectNo = 1;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + " " + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{		
			idToIndivNo[famIndivID] = subjectNo;
			indivNoToId[subjectNo] = famIndivID;
			subjectNo++;
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	unsigned int totalSubjects = subjectNo - 1;
	readFamilyFile.close();


	ifstream readCovariateFile;
	readCovariateFile.open(covariateFilename.c_str());
	if(!readCovariateFile.is_open())
	{
		outErr("Cannot read covariate file: "); outErr(covariateFilename); outErr("!\n");
		exit(0);
	};

	//determine number of covariates
	unsigned int noCovars = 3;
	ifstream readCovariateFileNoCovar(covariateFilename.c_str());
	string firstLine;
	getline(readCovariateFileNoCovar, firstLine);
	unsigned int noCols = 0;
	bool lastWasSpace = true;	

	for(unsigned int ch = 0; ch < firstLine.length(); ++ch)
	{
		if((firstLine.substr(ch, 1) != " " && firstLine.substr(ch, 1) != "\t") && lastWasSpace)
		{
			++noCols;
			lastWasSpace = false;
		};

		if(firstLine.substr(ch, 1) == " " || firstLine.substr(ch, 1) == "\t" ) lastWasSpace = true;
	};

	if(noCols <= 2)
	{
		outErr("Not enough columns in covariate file: "); outErr(covariateFilename); outErr("!\n");
		exit(0);
	};

	noCovars = noCols - 2;
	readCovariateFileNoCovar.close();

	//determine if header present
	ifstream readCovariateFileHeader(covariateFilename.c_str());
	string aStr;

	readCovariateFileHeader >> aStr >> aStr; //first two are expected strings, even if no header

	bool isANumber = false;
	double aNum;

	//if others are numbers assume there is no header
	for(unsigned int h = 3; h <= noCols; ++h)
	{
		readCovariateFileHeader >> aStr;
		aNum = atof(aStr.c_str());
		if(aNum != 0) isANumber = true;
		else if(aNum == 0 && (aStr.substr(0, 1) == "0" || (aStr.length() >= 2 && aStr.substr(0, 2) == "-0"))) isANumber = true; 
	};

	readCovariateFileHeader.close();

	bool headerExists = !isANumber;

	map<string, unsigned int> headerNames;
	map<string, unsigned int>::const_iterator hn;

	if(headerExists)
	{
		ifstream readCovariateFileHeaderNames(covariateFilename.c_str());

		readCovariateFileHeaderNames >> aStr >> aStr; //fam id, id
		readCovariateFile >> aStr >> aStr; //fam id, id

		for(unsigned int h = 1; h <= noCovars; ++h)
		{
			readCovariateFileHeaderNames >> aStr;
			headerNames[aStr] = h;
			readCovariateFile >> aStr;
		};

		readCovariateFileHeaderNames.close();
	};

	//pick covariates to analyse
	set<unsigned int> covariatesToUse;
	set<unsigned int>::const_iterator ctu;
	set<string> someCovariates;
	unsigned int aCol, bCol;

	if(covariates == "") for(unsigned int co = 1; co <= noCovars; ++co)	covariatesToUse.insert(co);
	else
	{		
		istringstream iss(covariates);
		string token;
		while(getline(iss, token, ','))
		{			
			someCovariates.insert(token);
		};

		unsigned int dashChr = 0;
		string str1, str2;

		//covariates may still be of form 1-3, so check for this
		for(set<string>::const_iterator sc = someCovariates.begin(); sc != someCovariates.end(); ++sc)
		{
			dashChr = 0;
			
			for(unsigned int cch = 0; cch < (*sc).length(); ++cch)
			{
				if((*sc).substr(cch, 1) == "-") dashChr = cch;
			};

			if(dashChr == 0)
			{
				str1 = *sc;
				aCol = atoi(str1.c_str());				

				if(aCol == 0)
				{
					hn = headerNames.find(str1);
					if(hn != headerNames.end()) aCol = hn->second; 
					else
					{
						outErr("Cannot find covariate column: "); outErr(str1); outErr("!\n");
						exit(0);
					};
				};
				
				covariatesToUse.insert(aCol);
			}
			else
			{
				//find covariates a-b, a to b				
				str1 = (*sc).substr(0, dashChr);
				str2 = (*sc).substr((dashChr+1), ((*sc).length()-dashChr-1));
				aCol = atoi(str1.c_str());
				bCol = atoi(str2.c_str());

				if(aCol == 0)
				{
					hn = headerNames.find(str1);
					if(hn != headerNames.end()) aCol = hn->second; 
					else
					{
						outErr("Cannot find covariate column: "); outErr(str1); outErr("!\n");
						exit(0);
					};
				};

				if(bCol == 0)
				{
					hn = headerNames.find(str2);
					if(hn != headerNames.end()) bCol = hn->second; 
					else
					{
						outErr("Cannot find covariate column: "); outErr(str2); outErr("!\n");
						exit(0);
					};
				};				

				//add all covariates from a to b
				for(unsigned int cab = aCol; cab <= bCol; ++cab) covariatesToUse.insert(cab);
			};
		};
	};

	prevFamIndivID = "";
	list<double> covariateValues;
	double aCovarValue;
	unsigned int coToUse;
	map<string, unsigned int>::const_iterator itin;
	map<unsigned int, list<double> > covariateDataSubjects;
	list<double> missingCovars;
	for(unsigned int m = 1; m <= noCovars; ++m)	missingCovars.push_back(missingQTValue);

	//loop thro' subjects and store the covariate values
	do{
		
		readCovariateFile >> famID >> indivID;
		famIndivID = famID + " " + indivID;

		covariateValues.clear();
		coToUse = 1;
		for(unsigned int i = 1; i <= noCovars; ++i)
		{
			readCovariateFile >> aCovarValue;
			ctu = covariatesToUse.find(i);
			if(ctu != covariatesToUse.end())
			{
				covariateValues.push_back(aCovarValue);
			};
		};

		//add covariate values to subject
		if(famIndivID != prevFamIndivID) 
		{		
			itin = idToIndivNo.find(famIndivID);
			if(itin != idToIndivNo.end())
			{					
				covariateDataSubjects[itin->second] = covariateValues;
			}
			else
			{
				outErr("Subject \""); outErr(famID); outErr(" "); outErr(indivID); outErr("\" in covariate file "); outErr(covariateFilename); outErr("\n");
				outErr(" cannot be found in family file "); outErr(famFilename); outErr("!\n");
				exit(0);
			};
		};

		prevFamIndivID = famIndivID;
	}while(!readCovariateFile.eof());

	//add covariate data for use later and set missing
	subjectNo = 1;
	bool existMissingCovar;	
	unsigned int noMissingIndivPhenos = 0;
	
	for(map<unsigned int, list<double> >::const_iterator cds = covariateDataSubjects.begin(); cds != covariateDataSubjects.end(); ++cds)
	{
		if(subjectNo != cds->first)
		{
			//numbers "subject No" to cds->first-1 missing
			for(unsigned int j = subjectNo; j < cds->first; ++j)
			{
				covariateDataAllSubjects.push_back(missingCovars);
				caseControls.push_back(0); //add missing individual
				++noMissingIndivPhenos;
			};
		};

		covariateDataAllSubjects.push_back(cds->second);
		//check for missing covariates, if one is missing add all as missing
		existMissingCovar = false;				
		for(list<double>::const_iterator chkc = cds->second.begin(); chkc != cds->second.end(); ++chkc)
		{
			if(*chkc == missingQTValue) existMissingCovar = true; 
		};

		if(existMissingCovar) caseControls.push_back(0); 
		else caseControls.push_back(1); //add all as controls (or missing) for now will be updated later

		subjectNo = cds->first + 1;
	};

	//add missing indivs at end
	for(unsigned int k = subjectNo; k <= totalSubjects; ++k)
	{
		covariateDataAllSubjects.push_back(missingCovars);
		caseControls.push_back(0); //add missing individual
		++noMissingIndivPhenos;
	};

	readCovariateFile.close();
};

