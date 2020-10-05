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


/*! \file Analysis.cpp
    \brief This file contains the methods the various analyse.
    
*/
#include <iostream>
#include <sstream>

using namespace std; // initiates the "std" or "standard" namespace

#include "Analysis.h"
#include "Statistics.h"
#include "main.h"

//! Runs the chosen analysis.
void Analysis::runAnalysis()
{
	if(snp1EndSNP == 0) snp1EndSNP = descSNPs->getNoSNPs(false);

	if(snp2EndSNP == 0)
	{
			if(filename2 == "") snp2EndSNP = descSNPs->getNoSNPs(false);
			else snp2EndSNP = descSNPs->getNoSNPs(true);
	};

	differentData = true;

	//if the second filename is not given use the first one again, set now to avoid unnec. checks between files for pedigree data 
	if(filename2 == "")
	{
		filename2 = filename;
		differentData = false;
	};

	//setup the SNP windows to analyse, advance thro' data file to before first SNPs in windows
	window1 = new SNPWindowReadFromFile(getTotalNoSubjects(), snp1StartSNP, filename, true);
	if(memoryType == 0) window2 = new SNPWindowReadFromFile(getTotalNoSubjects(), snp2StartSNP, filename2, snp2EndSNP); 
	else if(memoryType == 1) window2 = new SNPWindowStoreAllDataBinary(getTotalNoSubjects(), snp2StartSNP, filename2, snp2EndSNP);
	else window2 = new SNPWindowStoreAllData(getTotalNoSubjects(), snp2StartSNP, filename2, snp2EndSNP);

	//add reference of SNP windows to Statistics to access individual SNP data where nec., eg for LR
	allStatistics->addSNPWindows(window1, window2);

	unsigned int snp1No = snp1StartSNP;
	unsigned int snp2No = snp2StartSNP;
	descSNPs->setSNP1DescIterator(snp1No);
	descSNPs->setSNP2DescIterator(snp2No);
	bool isOK = true; //maximum number of results have not been exceeded
	bool tooCloseForCaseOnly = false;

	//loop through SNP 1
	for( ; ; )
	{

		//start reading the data from the start of SNP window 2
		window2->startWindowAtStart();

		snp2No = snp2StartSNP;
		descSNPs->setSNP2DescIterator(snp2No);

		//loop through SNP 2 to analysis SNP 1 against
		for( ; ; )
		{			
			//determine if too close for a case only test
			getTooCloseForCaseOnly(tooCloseForCaseOnly, snp1No, snp2No);			

			//Analysis SNP 1 against SNP 2, but do not repeat analysis or analyse against the same SNP
			if(differentData || snp2No > snp1No || snp2No < snp1StartSNP || snp1No > snp2EndSNP)
			{
				//calculate the joint genotype counts for SNP1 and SNP2
				updateJointGenotypeCounts();				

				//calculate the statistics and output results if below threshold conditions are met
				if(allStatistics->evaluateStatistics(tooCloseForCaseOnly)) isOK = allStatistics->outputResults(snp1No, snp2No);							
			};

			//move to the next SNP
			++snp2No; descSNPs->advanceIterator2();

			if(snp2No > snp2EndSNP || !isOK) break;

			//read in new SNP data for SNP 2
			window2->moveToNextSNP();
		};

		//move to the next SNP
		++snp1No; descSNPs->advanceIterator1();

		if(snp1No > snp1EndSNP || !isOK) break;

		//read in new SNP data for SNP 1
		window1->moveToNextSNP();
	};

	allStatistics->outputSummaries();
};

//! Sets whether the SNPs are too close for a case only analysis.
void Analysis::getTooCloseForCaseOnly(bool & tooCloseForCaseOnly, unsigned int & snp1No, unsigned int & snp2No)
{
	
	if(differentData) tooCloseForCaseOnly = false;
	else
	{
		int theGap = descSNPs->getSNPInfo1()->bp - descSNPs->getSNPInfo2()->bp;
		if(theGap > caseOnlyGap || theGap < -caseOnlyGap) tooCloseForCaseOnly = false;
		else tooCloseForCaseOnly = true;
	};
	
};

//! Converts an integer to a string
string toString(int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Returns a string of the run time
string getTime(const double & t)
{
	double time = t;
	int days = 0;
	int hours = 0;
	int minutes = 0;
	int seconds = 0;

	string ans = "";
	days = (int) (time / 86400); time -= days*86400;
	hours = (int) (time / 3600); time -= hours*3600;
	minutes = (int) (time / 60); time -= minutes*60;
	seconds = (int) time;

	if(days == 1) ans += "1 day";
	else if(days > 0) ans += toString(days) + " days";

	if(hours > 0)
	{
		if(days != 0)
		{
			if(minutes == 0 && seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(hours == 1) ans += "1 hour";
		else ans += toString(hours) + " hours";
	};

	if(minutes > 0)
	{
		if(ans != "")
		{
			if(seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(minutes == 1) ans += "1 minute";
		else ans += toString(minutes) + " minutes";
	};

	if(seconds > 0)
	{
		if(ans != "")
		{
			ans += " and ";			
		};

		if(seconds == 1) ans += "1 second";
		else ans += toString(seconds) + " seconds";
	};

	if(ans == "") ans = "less than one second";

	return ans;
};

void Analysis::updateJointGenotypeCounts()
{
	if(!calcJointGenotypes) return; //if only performing test not requiring joint counts, eg. LR

	//reset all joint genotype counts to zero in order to calculate new counts for next SNP pair
	caseJointGenotypeCounts->resetCounts();
	controlJointGenotypeCounts->resetCounts();

	SNPData * snp1 = window1->getCurrentSNP();

	//loop thro' all the subjects for each SNP
	list<unsigned char>::const_iterator s1 = snp1->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator s2 = window2->getCurrentSNP()->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator caco = caseControls.begin();

	//count the joint minor alleles for each subject
	do{

		//count cases and controls separately, do not count missing phenos
		if(*caco == 2) caseJointGenotypeCounts->addCount(*s1, *s2);
		else if(*caco == 1) controlJointGenotypeCounts->addCount(*s1, *s2);

		++s1; ++s2; ++caco;
	}while(s1 != snp1->noMinorAllelesAllSubjects.end());

};

//! Sets up family case/control data from the .fam file
void Analysis::setupCaseControls(string & filename)
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
	unsigned int noMissing = 0;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			if(phenoType == 2)
			{
				caseControls.push_back(2);	//case
				++noCases;
			}			
			else if(phenoType == 1) caseControls.push_back(1); //control
			else
			{
				caseControls.push_back(0); //missing pheno
				++noMissing;
			};
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	double totalNoSubjectsNoMiss = caseControls.size() - noMissing;

	double casesPercent = (((double)(noCases))/((double)(totalNoSubjectsNoMiss)))*100;
	double controlsPercent = (((double)(totalNoSubjectsNoMiss - noCases))/((double)(totalNoSubjectsNoMiss)))*100;

	out("Number of subjects: "); out(caseControls.size() ); out("\n");
	out("Number of cases: "); out(noCases); if(totalNoSubjectsNoMiss != 0) {out(" ("); out(casesPercent); out("%)");};  out("\n");
	out("Number of controls: "); out(caseControls.size() - noCases - noMissing); if(totalNoSubjectsNoMiss != 0) {out(" ("); out(controlsPercent); out("%)");};  out("\n");
	out("Number of missing: "); out(noMissing);  out("\n");
	out("\n");

	readFamilyFile.close();
};


//! Checks the family case/control data from the .fam file against that from first filename
void Analysis::checkCaseControls(string & filename2, string & filename)
{
	if(filename2 == "") return;

	//try and find the family file and read in data
	unsigned int length = filename2.length();
	string famFilename2 = filename2.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename2.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename2); outErr("!\n");
		exit(0);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	int phenoType;
	unsigned int noCases = 0;
	list<unsigned char>::const_iterator cc = caseControls.begin();

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
				
				if(((*cc == 1) && phenoType == 2) || ((*cc == 2) && phenoType != 2))
				{
					outErr("\nCase/Control status for "); outErr(filename); outErr(" and "); outErr(filename2); outErr(" do not match!\n (At subject no. "); outErr(noCases); outErr(")\n");
					exit(0);
				};
					
				cc++;
				noCases++;			
		};

		if(cc == caseControls.end()) break;

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	if(noCases != caseControls.size())
	{
		outErr("\nThe number of subjects in "); outErr(filename); outErr(" and "); outErr(filename2); outErr(" do not match!\n ("); outErr(caseControls.size()); outErr(" and "); outErr(noCases); outErr(" respectively)\n");
		exit(0);
	};

	readFamilyFile.close();
};


