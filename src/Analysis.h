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


/*! \file Analysis.h
    \brief This file organises the analysis.
    
*/

#ifndef __ANALYSIS
#define __ANALYSIS

#include <string>

using namespace std; // initiates the "std" or "standard" namespace

#include "Data.h"
#include "main.h"
#include "Statistics.h"

//! Returns a string of the run time
string getTime(const double & t);

//! Organises the correect analysis to perform.
class Analysis
{
private:

	unsigned int snp1StartSNP;
	unsigned int snp1EndSNP;
	unsigned int snp2StartSNP;
	unsigned int snp2EndSNP;
	int caseOnlyGap;

	SNPWindow * window1;
	SNPWindow * window2;
	DescriptionOfSNPs * descSNPs;
	JointGenotypeCounts * caseJointGenotypeCounts;
	JointGenotypeCounts * controlJointGenotypeCounts;
	bool calcJointGenotypes;
	bool differentData;
	list<unsigned char> caseControls; //list of bool values to whether a subject is a case or not
	unsigned int totalNoSubjects;
	unsigned int memoryType;

	string filename;
	string filename2;
	
	AllStatistics * allStatistics;

public:

	Analysis(string & fn, string & fn2, string & ofn, unsigned int & a1, unsigned int & a2, unsigned int & b1, unsigned int & b2, unsigned int & filt, unsigned int & mem, unsigned int & gp, unsigned int & mx, list<Statistic *> & sts) :  filename(fn), filename2(fn2), snp1StartSNP(a1), snp1EndSNP(a2), snp2StartSNP(b1), snp2EndSNP(b2), caseOnlyGap(gp), memoryType(mem)
	{ 
		//setup SNP description info. for both files
		descSNPs = new DescriptionOfSNPs(filename, filename2);
	
		//setup which subject are cases and which are controls
		setupCaseControls(filename);
		checkCaseControls(filename2, filename);
		totalNoSubjects = caseControls.size();

		caseJointGenotypeCounts = new JointGenotypeCounts();
		controlJointGenotypeCounts = new JointGenotypeCounts();
		allStatistics = new AllStatistics(ofn, sts, filt, mx, descSNPs, caseJointGenotypeCounts, controlJointGenotypeCounts, caseControls);
		calcJointGenotypes = allStatistics->getCalcJointGenotypes();
	};

	//! Delete analysis things
	~Analysis()
	{
		delete window1; 
		delete window2;
		delete descSNPs;		
		delete allStatistics;
		delete caseJointGenotypeCounts;
		delete controlJointGenotypeCounts;
	};

	void runAnalysis();

	void updateJointGenotypeCounts();
	void setupCaseControls(string & filename);
	void checkCaseControls(string & filename2, string & filename);
	void getTooCloseForCaseOnly(bool & tooCloseForCaseOnly, unsigned int & snp1No, unsigned int & snp2No);
	unsigned int getTotalNoSubjects() const {return totalNoSubjects;};	
};


#endif
