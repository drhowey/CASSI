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


/*! \file Statistics.cpp
    \brief This file contains the source of models used for logistic regression.
    
*/

#include <map>
#include <iostream>
#include <math.h>

#include <string>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Statistics.h"
#include "Data.h"
#include "Utils.h"

//! Outputs SNP pairs calculated and the number passing the thresholds.
void Statistic::outputTotalCounts()
{
	out("Total SNP pairs calculated: "); out(totalCalculations); out("\n");
	out("Total SNP pair statistics passing threshold: "); out(totalThresholdPassed); out("\n");
};

//! Sets threshold for case/control test.
void Statistic::setThresholdCC(double & th)
{
	if(th >= 1) thresholdCC = 0;
	else thresholdCC = calculateChiSqFromPvalue(th);
};

//! Sets threshold for case only test.
void Statistic::setThresholdCO(double & th)
{
	if(th >= 1) thresholdCC = 0;
	thresholdCO = calculateChiSqFromPvalue(th);
};

//! Determines whether joint genotype counts should be calculated or not.
bool AllStatistics::getCalcJointGenotypes() const
{
	bool ans = false;
	for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
	{		
		ans = ans || (*s)->getCalcJointGenotypes();
	};
	return ans;
};


//! evaluate all the statistics.
bool AllStatistics::evaluateStatistics(bool & tooCloseForCaseOnly)
{
	bool reportResult;

	//reset all test statistics initially not to have passed their thresholds
	for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
	{		
		(*s)->setThresholdPassedFalse();
	};

	if(filterType == 1)
	{
		reportResult = true;
		//do not calculate remaining statistics if one does not pass its threshold
		for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
		{		
			(*s)->evaluateStatistic(tooCloseForCaseOnly);
			if(!(*s)->getThresholdPassed()) {reportResult = false; break;};
		};
	}
	else
	{
		reportResult = false;
		//need to calculate all statistics for these filters
		for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
		{		
			(*s)->evaluateStatistic(tooCloseForCaseOnly);
			if((*s)->getThresholdPassed()) reportResult = true;
		};
	};

	return reportResult;
};

//! Adds reference to SNP Windows.
void AllStatistics::addSNPWindows(SNPWindow * window1, SNPWindow * window2)
{
	for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
	{		
		(*s)->addSNPWindows(window1, window2);
	};
};

//! Output appropriate header for output file.
void AllStatistics::outputHeader()
{	
	resultsFile << "SNP1 CHR1 ID1 BP1 SNP2 CHR2 ID2 BP2";	
	for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
	{
		resultsFile << " ";
		(*s)->outputHeader(resultsFile);
	};
	resultsFile << "\n";
};

//! Outputs a summary of each test statistic used.
void AllStatistics::outputSummaries()
{
	for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
	{		
		(*s)->outputSummaryResults();		
	};	

	outputTotalResults();
};

//! Output calculated results.
bool AllStatistics::outputResults(unsigned int & snp1No, unsigned int & snp2No)
{
	if(totalResults >= maxNoResults && maxNoResults != 0)
	{
		outErr("Total number of results have exceeded the maximum allowed("); outErr(maxNoResults); outErr(")!\n\n");
		return false;
	};

	SNPInfo * snp1 = descSNPs->getSNPInfo1();//descSNPs->getSNPInfo(snp1No, 1);
	SNPInfo * snp2 = descSNPs->getSNPInfo2();//descSNPs->getSNPInfo(snp2No, 2);

	resultsFile << snp1No << " " << snp1->chromosome << " " << snp1->name << " " << snp1->bp << " " << snp2No << " " << snp2->chromosome << " " << snp2->name << " " << snp2->bp;

	
	//output all stat results that were labelled as passing a threshold
	for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
	{
		resultsFile << " ";
		(*s)->outputResults(resultsFile);
	};
	
	resultsFile << "\n";

	totalResults++;

	return true;
};

//! Constructor with thresholds for Joint Effects statistic.
JointEffects::JointEffects(double & thcc, double & thco, unsigned int & cm) : logCache()
{			
	adjJointGenotypeCountsCase = new AdjustedJointGenotypeCounts();
	adjJointGenotypeCountsCon = new AdjustedJointGenotypeCounts();
	thresholdCC = calculateChiSqFromPvalue(thcc);
	thresholdCO = calculateChiSqFromPvalue(thco);
	cellMin = cm;
};

//! Calculate the odds ratio relative risk.
double JointEffects::getOddRatioRelRisk(double & num1, double & num2, double & denom1, double & denom2)
{
	return (num1*num2)/(denom1*denom2);
};

//! Return the log value from the cache, if not add it
double JointEffects::getLog(double v)
{
	map<double, double>::const_iterator i = logCache.find(v);
	if(i == logCache.end())
	{
		double logv = log(v);
		logCache[v] = logv;
		return logv;
	}
	else
	{		
		return i->second;
	};
};

//! Calculates i_** for joint effects test.
void JointEffects::calculateis(double & i22, double & i21, double & i12, double & i11, AdjustedJointGenotypeCounts * ajgc)
{
	i22 = (ajgc->counts[2][2]*ajgc->counts[0][0])/(ajgc->counts[2][0]*ajgc->counts[0][2]);
	i21 = (ajgc->counts[2][1]*ajgc->counts[0][0])/(ajgc->counts[2][0]*ajgc->counts[0][1]);
	i12 = (ajgc->counts[1][2]*ajgc->counts[0][0])/(ajgc->counts[1][0]*ajgc->counts[0][2]);
	i11 = (ajgc->counts[1][1]*ajgc->counts[0][0])/(ajgc->counts[1][0]*ajgc->counts[0][1]);
};

//! Calculates lambda and the variance needed for the joint effect statistic.
void JointEffects::calculateVariables(double & i22, double & i21, double & i12, double & i11, double & totalInverseV, double rowTotals[], AdjustedJointGenotypeCounts * ajc, const bool & useAltStat)
{
	
	//partial diff of f with xi, a diagonal 4 x 4 matrix
	double pdfpdXi1, pdfpdXi2, pdfpdXi3, pdfpdXi4;

	if(useAltStat)
	{
		pdfpdXi1 = sqrt(i22)/(2*sqrt(i22) + 2);
		pdfpdXi2 = i21/(i21 + 1);
		pdfpdXi3 = i12/(i12 + 1);
		pdfpdXi4 = 1.0;
	}
	else
	{
		pdfpdXi1 = 0.5;
		pdfpdXi2 = 1.0;
		pdfpdXi3 = 1.0;
		pdfpdXi4 = 2*i11/(2*i11 - 1);
	};

	double invq22, invq02, invq20, invq00, invq21, invq12, invq01, invq10, invq11;
	
	invq00 = 1.0/ajc->counts[0][0];
	invq01 = 1.0/ajc->counts[0][1];
	invq02 = 1.0/ajc->counts[0][2];
	invq10 = 1.0/ajc->counts[1][0];
	invq11 = 1.0/ajc->counts[1][1];
	invq12 = 1.0/ajc->counts[1][2];
	invq20 = 1.0/ajc->counts[2][0];
	invq21 = 1.0/ajc->counts[2][1];
	invq22 = 1.0/ajc->counts[2][2];

	//Do not add in n^(-1) scalar factor of matrix C (and V) add into w estimates later 
	
	//build matrix V firstly
	list< list<double> > matrixV, inverse;

	list<double> row1;
	row1.push_back((invq22 + invq02 + invq20 + invq00)*pdfpdXi1*pdfpdXi1);
	double v12 = (invq20 + invq00)*pdfpdXi1*pdfpdXi2;
	row1.push_back(v12);
	double v13 = (invq02 + invq00)*pdfpdXi1*pdfpdXi3;
	row1.push_back(v13);
	double v14 = invq00*pdfpdXi1*pdfpdXi4;
	row1.push_back(v14);

	matrixV.push_back(row1);

	list<double> row2;
	row2.push_back(v12);
	row2.push_back((invq21 + invq20 + invq01 + invq00)*pdfpdXi2*pdfpdXi2);
	double v23 = invq00*pdfpdXi2*pdfpdXi3;
	row2.push_back(v23);
	double v24 = (invq01 + invq00)*pdfpdXi2*pdfpdXi4;
	row2.push_back(v24);

	matrixV.push_back(row2);

	list<double> row3;
	row3.push_back(v13);
	row3.push_back(v23);
	row3.push_back((invq12 + invq10 + invq02 + invq00)*pdfpdXi3*pdfpdXi3);
	double v34 = (invq10 + invq00)*pdfpdXi3*pdfpdXi4;
	row3.push_back(v34);

	matrixV.push_back(row3);

	list<double> row4;
	row4.push_back(v14);
	row4.push_back(v24);
	row4.push_back(v34);
	row4.push_back((invq11 + invq10 + invq01 + invq00)*pdfpdXi4*pdfpdXi4);

	matrixV.push_back(row4);

	//set up the inverse of the matrix - perhaps this can be done better within the calculating inverse method
	row1.clear(); row2.clear(); row3.clear(); row4.clear();

	row1.push_back(1); row1.push_back(0); row1.push_back(0); row1.push_back(0);
	row2.push_back(0); row2.push_back(1); row2.push_back(0); row2.push_back(0);
	row3.push_back(0); row3.push_back(0); row3.push_back(1); row3.push_back(0);
	row4.push_back(0); row4.push_back(0); row4.push_back(0); row4.push_back(1);

	inverse.push_back(row1);
	inverse.push_back(row2);
	inverse.push_back(row3);
	inverse.push_back(row4);

	//calculate inverse of V
	getInverseMatrix(matrixV, inverse);

	//sum all elements of inverse of V to give variance (1^T V^(-1) 1)^(-1)
	unsigned int rowCounter = 0;//need to start counting from 0 for array.

	//loop thro' rows
	for(list< list<double> >::const_iterator r = inverse.begin(); r != inverse.end(); ++r)	
	{

		//loop thro' columns
		for(list<double>::const_iterator c = r->begin(); c != r->end(); ++c)		
		{			
			rowTotals[rowCounter] += *c;
		};

		++rowCounter;
	};

	//not including the 1/total scalar
	totalInverseV = 1.0/(rowTotals[0]+rowTotals[1]+rowTotals[2]+rowTotals[3]);

};

//! Adjust the genotype counts, adding a half if req or needed
void AdjustedJointGenotypeCounts::setAdjustedCounts(const bool & addHalf, JointGenotypeCounts * jointGenotypeCounts, const bool & adjustBaseline)
{
	double adj = 0;
	
	//add 0.5 if any cells are 0
	if(addHalf || jointGenotypeCounts->counts[0][0] == 0 || jointGenotypeCounts->counts[0][1] == 0
		|| jointGenotypeCounts->counts[0][2] == 0 || jointGenotypeCounts->counts[1][0] == 0 || jointGenotypeCounts->counts[2][0] == 0
		|| jointGenotypeCounts->counts[1][1] == 0 || jointGenotypeCounts->counts[1][2] == 0 || jointGenotypeCounts->counts[2][1] == 0
		|| jointGenotypeCounts->counts[2][2] == 0)
	{
		halfAdded = true;
		adj = 0.5;		
	}
	else
		halfAdded = false;

	counts[0][0] = adj + (double)(jointGenotypeCounts->counts[0][0]);
	counts[0][1] = adj + (double)(jointGenotypeCounts->counts[0][1]);
	counts[0][2] = adj + (double)(jointGenotypeCounts->counts[0][2]);
	counts[1][0] = adj + (double)(jointGenotypeCounts->counts[1][0]);
	counts[1][1] = adj + (double)(jointGenotypeCounts->counts[1][1]);
	counts[1][2] = adj + (double)(jointGenotypeCounts->counts[1][2]);
	counts[2][0] = adj + (double)(jointGenotypeCounts->counts[2][0]);
	counts[2][1] = adj + (double)(jointGenotypeCounts->counts[2][1]);
	counts[2][2] = adj + (double)(jointGenotypeCounts->counts[2][2]);


	//adjust all the elements if the baseline frequency is less than 0.01
	if(adjustBaseline && counts[0][0]/(double)(jointGenotypeCounts->total) < 0.01)
	{
		double factor = ((double)(jointGenotypeCounts->total))/(0.01*(double)(jointGenotypeCounts->total) + counts[0][1] + counts[0][2] + counts[1][0] +
			             counts[1][1] + counts[1][2] + counts[2][0] + counts[2][1]  + counts[2][2]);

		counts[0][0] = 0.01*factor*(double)(jointGenotypeCounts->total);
		counts[0][1] *= factor;
		counts[0][2] *= factor;
		counts[1][0] *= factor;
		counts[1][1] *= factor;
		counts[1][2] *= factor;
		counts[2][0] *= factor;
		counts[2][1] *= factor;
		counts[2][2] *= factor;
	};

};

//! Check cell counts are at the required minimum or more.
bool checkMinimumCellCounts(const unsigned int & cellMin, JointGenotypeCounts * jointGenotypeCounts)
{
	return (jointGenotypeCounts->counts[2][2] < cellMin || jointGenotypeCounts->counts[1][2] < cellMin || jointGenotypeCounts->counts[2][1] < cellMin || 
		 jointGenotypeCounts->counts[0][2] < cellMin || jointGenotypeCounts->counts[2][0] < cellMin || jointGenotypeCounts->counts[1][1] < cellMin ||
		 jointGenotypeCounts->counts[0][1] < cellMin || jointGenotypeCounts->counts[1][0] < cellMin || jointGenotypeCounts->counts[0][0] < cellMin);
};

//! Calculates the statistic for the joint effects test.
void JointEffects::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	cellCountTooSmallCase = checkMinimumCellCounts(cellMin, caseJointGenotypeCounts);
	if(cellCountTooSmallCase) return;

	double diff;
	bool halfAddedCases2ndTime = false;
	//case variables
	double i22ca, i21ca, i12ca, i11ca;
	double rowTotalsca[4] = {0,0,0,0};
	//control variables
	double i22co, i21co, i12co, i11co;
	double rowTotalsco[4] = {0,0,0,0};

	useAltStatCases = false;
	useAltStatCaseCons = false;
	muAOverTotalInverseVca = 0;

	//adjust the genotype counts for cases to account for zeroes and if the baseline freq. is too low (< 0.01)
	adjJointGenotypeCountsCase->setAdjustedCounts(false, caseJointGenotypeCounts);

	calculateis(i22ca, i21ca, i12ca, i11ca, adjJointGenotypeCountsCase);

	useAltStatCases = (i11ca <= 0.5);

	//calculate variables needed to calculate the case-only stat
	calculateVariables(i22ca, i21ca, i12ca, i11ca, totalInverseVca, rowTotalsca, adjJointGenotypeCountsCase, useAltStatCases);

	//calculate case only statistic 
	if(useAltStatCases)
	{		
		muAOverTotalInverseVca = (rowTotalsca[0]*log((sqrt(i22ca) + 1)*0.5) + rowTotalsca[1]*log((i21ca + 1)*0.5) + rowTotalsca[2]*log((i12ca + 1)*0.5) + rowTotalsca[3]*log(i11ca));

		caseOnlyStat = muAOverTotalInverseVca*muAOverTotalInverseVca*totalInverseVca;
	}
	else
	{
		//use w's to calculate lambda
		lambdaAOverTotalInverseVca = (rowTotalsca[0]*log(i22ca)*0.5 + rowTotalsca[1]*log(i21ca) + rowTotalsca[2]*log(i12ca) + rowTotalsca[3]*log(2*i11ca - 1));

		caseOnlyStat = lambdaAOverTotalInverseVca*lambdaAOverTotalInverseVca*totalInverseVca;	
	};
	
	if(calcCaseControl)
	{
		cellCountTooSmallCon = checkMinimumCellCounts(cellMin, controlJointGenotypeCounts);
		if(!cellCountTooSmallCon)
		{
			//adjust the genotype counts for controls to account for zeroes and if the baseline freq. is too low (< 0.01)
			adjJointGenotypeCountsCon->setAdjustedCounts(adjJointGenotypeCountsCase->halfAdded, controlJointGenotypeCounts);

			//if 0.5 was added to the control counts then we also need to add 0.5 to the case counts if not already done so
			if(adjJointGenotypeCountsCon->halfAdded && !adjJointGenotypeCountsCase->halfAdded)
			{
				halfAddedCases2ndTime = true;		
				adjJointGenotypeCountsCase->setAdjustedCounts(true, caseJointGenotypeCounts);

				calculateis(i22ca, i21ca, i12ca, i11ca, adjJointGenotypeCountsCase);
			};
	
			//calculate i's for controls
			calculateis(i22co, i21co, i12co, i11co, adjJointGenotypeCountsCon);

			useAltStatCaseCons = (i11co <= 0.5) || (i11ca <= 0.5);

			//case variables need to be recalculated if counts were adjusted or the alt stat is to be used where it was not previously (or vice versa)
			if(halfAddedCases2ndTime || (useAltStatCaseCons && !useAltStatCases) || (!useAltStatCaseCons && useAltStatCases))
			{
		
				if(!halfAddedCases2ndTime && ((useAltStatCaseCons && !useAltStatCases) || (!useAltStatCaseCons && useAltStatCases)))
				{
					calculateis(i22ca, i21ca, i12ca, i11ca, adjJointGenotypeCountsCase);
				};

				//reset row totals for recalculation
				rowTotalsca[0] = 0; rowTotalsca[1] = 0; rowTotalsca[2] = 0; rowTotalsca[3] = 0;
				//now recalculate the case variables with the new adjusted counts or for the alt stat for the first time
				calculateVariables(i22ca, i21ca, i12ca, i11ca, totalInverseVca, rowTotalsca, adjJointGenotypeCountsCase, useAltStatCaseCons);
			};

			//calculate variables for controls, either lambda or mu calcs, (and redo cases if nec.)

			//calculate the variables needed to calculate the case-control stat from the control data
			calculateVariables(i22co, i21co, i12co, i11co, totalInverseVco, rowTotalsco, adjJointGenotypeCountsCon, useAltStatCaseCons);
	
			if(useAltStatCaseCons)
			{		
				muNOverTotalInverseVco = (rowTotalsco[0]*log((sqrt(i22co) + 1)*0.5) + rowTotalsco[1]*log((i21co + 1)*0.5) + rowTotalsco[2]*log((i12co + 1)*0.5) + rowTotalsco[3]*log(i11co));

				if(halfAddedCases2ndTime || !useAltStatCases)
				{
					muAOverTotalInverseVca = (rowTotalsca[0]*getLog((sqrt(i22ca) + 1)*0.5) + rowTotalsca[1]*getLog((i21ca + 1)*0.5) + rowTotalsca[2]*getLog((i12ca + 1)*0.5) + rowTotalsca[3]*getLog(i11ca));
				};

				diff = (muAOverTotalInverseVca*totalInverseVca - muNOverTotalInverseVco*totalInverseVco);
				caseControlStat = diff*diff/(totalInverseVca + totalInverseVco);
			}
			else
			{
				//use w's to calculate lambda
				lambdaNOverTotalInverseVco = (rowTotalsco[0]*log(i22co)*0.5 + rowTotalsco[1]*log(i21co) + rowTotalsco[2]*log(i12co) + rowTotalsco[3]*log(2*i11co - 1));

				if(halfAddedCases2ndTime || useAltStatCases)
				{
					lambdaAOverTotalInverseVca = (rowTotalsca[0]*log(i22ca)*0.5 + rowTotalsca[1]*log(i21ca) + rowTotalsca[2]*log(i12ca) + rowTotalsca[3]*log(2*i11ca - 1));
				};

				diff = (lambdaAOverTotalInverseVca*totalInverseVca - lambdaNOverTotalInverseVco*totalInverseVco);
				caseControlStat = diff*diff/(totalInverseVca + totalInverseVco);		
			};
		}
		else caseControlStat = 0;
	};
	

	if((calcCaseControl && caseControlStat >= thresholdCC) ||
	   (calcCaseOnly && !tooCloseForCaseOnly && caseOnlyStat >= thresholdCO))
	{
       thresholdPassed = true;
	   totalThresholdPassed++;
	};

	totalCalculations++;
};

//! Outputs the header for the joint effects test.
void JointEffects::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
		//For testing
		//resultsFile << "CA00 CA01 CA02 CA10 CA11 CA12 CA20 CA21 CA22 ";
		//resultsFile << "CO00 CO01 CO02 CO10 CO11 CO12 CO20 CO21 CO22 ALT_CO ALT_CC ";

		resultsFile << "JE_CASE_LOG_OR JE_CASE_SE ";
		if(calcCaseControl) resultsFile << "JE_CTRL_LOG_OR JE_CTRL_SE JE_CC_CHISQ JE_CC_P JE_CC_ALT";	
		if(calcCaseControl && calcCaseOnly) resultsFile << " ";
		if(calcCaseOnly) resultsFile << "JE_CO_CHISQ JE_CO_P JE_CO_ALT";
	};
};

//!For testing
//void outputCounts(JointGenotypeCounts & jc, ofstream & resultsFile)
//{
//	out("\n");
//	out(jc.counts[0][0]); out(" "); out(jc.counts[0][1]); out(" "); out(jc.counts[0][2]); out("\n");
//	out(jc.counts[1][0]); out(" "); out(jc.counts[1][1]); out(" "); out(jc.counts[1][2]); out("\n");
//	out(jc.counts[2][0]); out(" "); out(jc.counts[2][1]); out(" "); out(jc.counts[2][2]); out("\n\n");
//	resultsFile << jc.counts[0][0] << " "<< jc.counts[0][1] << " " << jc.counts[0][2] << " ";
//	resultsFile << jc.counts[1][0] << " "<< jc.counts[1][1] << " " << jc.counts[1][2] << " ";
//	resultsFile << jc.counts[2][0] << " "<< jc.counts[2][1] << " " << jc.counts[2][2] << " ";
//};

//! Output results of the joint effects test.
void JointEffects::outputResults(ofstream & resultsFile)
{
	//For testing
	//outputCounts(*caseJointGenotypeCounts, resultsFile);
	//outputCounts(*controlJointGenotypeCounts, resultsFile);
	//resultsFile << useAltStatCases << " " << useAltStatCaseCons << " ";

	if(!suppressResults) 
	{
		//output alternative variables if used, alt variance already set
		if(useAltStatCases) lambdaAOverTotalInverseVca = muAOverTotalInverseVca;
		if(useAltStatCaseCons) lambdaNOverTotalInverseVco = muNOverTotalInverseVco;
		
		if(!cellCountTooSmallCase) resultsFile << lambdaAOverTotalInverseVca*totalInverseVca << " " << sqrt(totalInverseVca) << " ";
		else resultsFile << "NA NA ";

		if(calcCaseControl)
		{
			if(!cellCountTooSmallCase && !cellCountTooSmallCon) resultsFile << lambdaNOverTotalInverseVco*totalInverseVco << " " << sqrt(totalInverseVco) << " " << caseControlStat << " " << getPvalueChiSq1DF(caseControlStat);
			else resultsFile << "NA NA NA NA";
			if(useAltStatCaseCons) resultsFile << " Y"; else resultsFile << " N"; 
		};

		if(calcCaseControl && calcCaseOnly) resultsFile << " ";

		if(calcCaseOnly)
		{
			if(!cellCountTooSmallCase) resultsFile << caseOnlyStat << " " << getPvalueChiSq1DF(caseOnlyStat);
			else resultsFile << "NA NA";
			if(useAltStatCases) resultsFile << " Y"; else resultsFile << " N";
		};
	};
};

//! Outputs summary of results for joint effects test.
void JointEffects::outputSummaryResults()
{
	out("Test Statistic: Joint Effects\n");
	if(calcCaseControl) {out("P-value threshold for case/control results: "); out(getPvalueChiSq1DF(thresholdCC)); out("\n");};
	if(calcCaseOnly) {out("P-value threshold for case only results: "); out(getPvalueChiSq1DF(thresholdCO)); out("\n");};
	out("Minimum cell count required for test: "); out(cellMin); out("\n");
	outputTotalCounts();
	out("\n");
};

//! Constructor with thresholds for Fast Epistasis statistic.
AdjustedFastEpistasis::AdjustedFastEpistasis(double & thcc, double & thco)
{		
	adjJointGenotypeCountsCase = new AdjustedJointGenotypeCounts();
	adjJointGenotypeCountsCon = new AdjustedJointGenotypeCounts();
	thresholdCC = calculateChiSqFromPvalue(thcc);
	thresholdCO = calculateChiSqFromPvalue(thco);
};

//! Calculates log odds and adjusted variance for Fast Epistatis test.
void AdjustedFastEpistasis::calculateLogOddsAdjustedVariance(double & logOdds, double & variance, AdjustedJointGenotypeCounts * jointGenotypeCounts)
{
	//Collapsed cell counts 
	double P1 = 4*jointGenotypeCounts->counts[0][0] + 2*jointGenotypeCounts->counts[1][0] + 2*jointGenotypeCounts->counts[0][1] + jointGenotypeCounts->counts[1][1];
	double P2 = 4*jointGenotypeCounts->counts[0][2] + 2*jointGenotypeCounts->counts[1][2] + 2*jointGenotypeCounts->counts[0][1] + jointGenotypeCounts->counts[1][1];
	double P3 = 4*jointGenotypeCounts->counts[2][0] + 2*jointGenotypeCounts->counts[1][0] + 2*jointGenotypeCounts->counts[2][1] + jointGenotypeCounts->counts[1][1];	
	double P4 = 4*jointGenotypeCounts->counts[2][2] + 2*jointGenotypeCounts->counts[1][2] + 2*jointGenotypeCounts->counts[2][1] + jointGenotypeCounts->counts[1][1];

	double rP1 = 1.0/P1;
	double rP2 = 1.0/P2;
	double rP3 = 1.0/P3;
	double rP4 = 1.0/P4;
	
	double bit2 = 2*(rP1 - rP2);
	double bit4 = 2*(rP1 - rP3);
	double bit5 = (rP1 - rP2 - rP3 + rP4);
	double bit6 = 2*(-rP2 + rP4);
	double bit8 = 2*(-rP3 + rP4);

	logOdds = log((P1*P4)/(P2*P3));

	variance = (16*rP1*rP1*jointGenotypeCounts->counts[0][0] + bit2*bit2*jointGenotypeCounts->counts[0][1]
			 + 16*rP2*rP2*jointGenotypeCounts->counts[0][2] + bit4*bit4*jointGenotypeCounts->counts[1][0]
			 + bit5*bit5*jointGenotypeCounts->counts[1][1] + bit6*bit6*jointGenotypeCounts->counts[1][2]
			 + 16*rP3*rP3*jointGenotypeCounts->counts[2][0] + bit8*bit8*jointGenotypeCounts->counts[2][1]
			 + 16*rP4*rP4*jointGenotypeCounts->counts[2][2]);

	//variance = (rP1+rP2+rP3+rP4); //variance for non-adjusted test
};

//! Calculates the statistic for the Adjusted Fast Epistasis test.
void AdjustedFastEpistasis::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	
	//adjust the genotype counts for cases to account for zeroes
	adjJointGenotypeCountsCase->setAdjustedCounts(false, caseJointGenotypeCounts, false);

	//calculate log odds and adjusted variance for cases
	calculateLogOddsAdjustedVariance(caseLogOdds, caseVar, adjJointGenotypeCountsCase);	

	chiSqStatCO = (caseLogOdds*caseLogOdds)/caseVar;

	if(calcCaseControl)
	{
		//adjust the genotype counts for controls to account for zeroes		
		adjJointGenotypeCountsCon->setAdjustedCounts(false, controlJointGenotypeCounts, false);

		//calculate log odds and adjusted variance for controls
		calculateLogOddsAdjustedVariance(conLogOdds, conVar, adjJointGenotypeCountsCon);	
	
		double diff = caseLogOdds - conLogOdds; 
		chiSqStatCC = (diff*diff)/(caseVar + conVar);

	};

	if((calcCaseControl && chiSqStatCC >= thresholdCC) ||
	   (calcCaseOnly && !tooCloseForCaseOnly && chiSqStatCO >= thresholdCO))
	{
		thresholdPassed = true;
		totalThresholdPassed++;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	};
	
	totalCalculations++;
};

//! Outputs the header for the Fast Epistasis test.
void AdjustedFastEpistasis::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
		resultsFile << "AFE_CASE_LOG_OR AFE_CASE_SE ";
		if(calcCaseControl) resultsFile << "AFE_CTRL_LOG_OR AFE_CTRL_SE AFE_CC_CHISQ AFE_CC_P";
		if(calcCaseControl && calcCaseOnly) resultsFile << " ";
		if(calcCaseOnly) resultsFile << "AFE_CO_CHISQ AFE_CO_P"; 
	};
};

//! Output results of the Fast Epistasis test.
void AdjustedFastEpistasis::outputResults(ofstream & resultsFile)
{
	if(!suppressResults)
	{		
		resultsFile << caseLogOdds << " " << sqrt(caseVar) << " ";
		if(calcCaseControl) resultsFile << conLogOdds << " " << sqrt(conVar) << " " << chiSqStatCC << " " << getPvalueChiSq1DF(chiSqStatCC);
		if(calcCaseControl && calcCaseOnly) resultsFile << " ";
		if(calcCaseOnly) resultsFile << chiSqStatCO << " " << getPvalueChiSq1DF(chiSqStatCO);
	};
};

//! Outputs summary of results for Fast Epistasis test.
void AdjustedFastEpistasis::outputSummaryResults()
{
	out("Test Statistic: Adjusted Fast Epistasis\n");
	if(calcCaseControl) {out("P-value threshold for case/control results: "); out(getPvalueChiSq1DF(thresholdCC)); out("\n");};
	if(calcCaseOnly) {out("P-value threshold for case only results: "); out(getPvalueChiSq1DF(thresholdCO)); out("\n");};
	outputTotalCounts();
	out("\n");
};

//! Constructor with thresholds for Adjusted Wu statistic.
AdjustedWu::AdjustedWu(double & thcc, double & thco) 
{			
	thresholdCC = calculateChiSqFromPvalue(thcc);
	thresholdCO = calculateChiSqFromPvalue(thco);
};

//! Calculates the log odds and the variance after haplotype freqs have been estimated.
void AdjustedWu::calculateLogORVar(double & logOdds, double & variance)
{
	//adjust estimated haplotypes if too small
	if(P11 < 1e-8) P11 = 1e-8;
	if(P12 < 1e-8) P12 = 1e-8;
	if(P21 < 1e-8) P21 = 1e-8;
	if(P22 < 1e-8) P22 = 1e-8;

	double sum = P11 + P12 + P21 + P22;
	P11 = P11/sum;	
	P12 = P12/sum;
	P21 = P21/sum;
	P22 = P22/sum;
	
	logOdds = log((P11*P22/(P12*P21)));

	p = P11 + P12;
	u = P11 + P21;
	D = P11 -p*u;
	q = 1 - p;
	v = 1 - u;

	//as defined in Ueki and Cordell Appendix A
	double phi = (P11*P22*(P11 + P22) + P12*P21*(P12 + P21))/((P11*P22 + P12*P21)*P11*P12*P21*P22);

	//define elements of sigma 3 by 3 matrix (without 1/2 constant, to be put in final result)
	double b1, b2, b3, b4, b5, b6, b7, b8, b9;

	b1 = p*q;		b2 = D;			b3 = D*(q-p);
	b4 = D;			b5 = u*v;		b6 = D*(v-u);
	b7 = D*(q-p);	b8 = D*(v-u);

	b9 = (1.0/phi + D*D*(p*q*(v-u)*(v-u) + u*v*(q-p)*(q-p)) - 2*D*D*D*(q-p)*(v-u))/(p*q*u*v - D*D);

	//define elements of vector
	double a1, a2, a3;

	a1 = u/P11 - v/P12 + u/P21 - v/P22;
	a2 = p/P11 + p/P12 - q/P21 - q/P22;
	a3 = 1.0/P11 + 1.0/P12 + 1.0/P21 + 1.0/P22;

	variance = (a1*(b1*a1 + b2*a2 + b3*a3) + a2*(b4*a1 + b5*a2 + b6*a3) + a3*(b7*a1 + b8*a2 + b9*a3))/(2*n);

	//variance = a3/(2*n); //variance for non-adjusted test
};

//! Calculates the statistic for the Adjusted Wu test.
void AdjustedWu::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	//estimate haplotype frequencies P11, P12, P21, P22 for cases
	estimateHaplotypeFreqsEM(caseJointGenotypeCounts, P11, P12, P21, P22);

	//set total number of people to cases
	n = caseJointGenotypeCounts->total;

	//calculate the logs odds and variance for the cases
	calculateLogORVar(lambdaA, varianceCases);
	
	caP11 = P11;
	caP12 = P12;
	caP21 = P21;
	caP22 = P22;

	caseOnlyStat = (lambdaA*lambdaA)/varianceCases;

	if(calcCaseControl)
	{
		//estimate haplotype frequencies P11, P12, P21, P22 for controls
		estimateHaplotypeFreqsEM(controlJointGenotypeCounts, P11, P12, P21, P22);

		//set total number of people to controls
		n = controlJointGenotypeCounts->total;

		//calculate the logs odds and variance for the controls
		calculateLogORVar(lambdaN, varianceControls);

		//set case/control statistic and case only statistic
		double diff = (lambdaA - lambdaN);
		caseControlStat = diff*diff/(varianceCases + varianceControls);
	};
	
	if((calcCaseControl && caseControlStat >= thresholdCC) ||
	   (calcCaseOnly && !tooCloseForCaseOnly && caseOnlyStat >= thresholdCO))
	{
		thresholdPassed = true;	
		totalThresholdPassed++;
	};

	totalCalculations++;
};


//! Outputs the header for the Adjusted Wu test.
void AdjustedWu::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
			resultsFile << "AWU_CASE_LOG_OR AWU_CASE_SE ";
			if(calcCaseControl) resultsFile << "AWU_CTRL_LOG_OR AWU_CTRL_SE AWU_CC_CHISQ AWU_CC_P"; 
			if(calcCaseControl && calcCaseOnly) resultsFile << " ";
			if(calcCaseOnly) resultsFile << "AWU_CO_CHISQ AWU_CO_P";			
	};
};

//! Output results of the Adjusted Wu test.
void AdjustedWu::outputResults(ofstream & resultsFile)
{
	if(!suppressResults)
	{
			if(lambdaA*0 != 0) resultsFile << "NA NA ";
			else resultsFile << lambdaA << " " << sqrt(varianceCases) << " ";

			if(calcCaseControl)
			{
				if(lambdaN*0 != 0) resultsFile << "NA NA NA NA";
				else
				{
						resultsFile << lambdaN << " " << sqrt(varianceControls) << " ";

						if(caseControlStat*0 != 0) resultsFile << "NA NA";
						else resultsFile << caseControlStat << " " << getPvalueChiSq1DF(caseControlStat);
				};
			};

			if(calcCaseControl && calcCaseOnly) resultsFile << " ";
			if(calcCaseOnly)
			{
				if(lambdaA*0 != 0) resultsFile << "NA NA";
				else resultsFile << caseOnlyStat << " " << getPvalueChiSq1DF(caseOnlyStat); 
			};
	};
};

//! Outputs summary of results for Adjusted Wu test.
void AdjustedWu::outputSummaryResults()
{
	out("Test Statistic: Adjusted Wu\n");
	if(calcCaseControl) {out("P-value threshold for case/control results: "); out(getPvalueChiSq1DF(thresholdCC)); out("\n");};
	if(calcCaseOnly) {out("P-value threshold for case only results: "); out(getPvalueChiSq1DF(thresholdCO)); out("\n");};
	outputTotalCounts();
	out("\n");
};

//! Constructor with thresholds for Wellek Ziegler statistic.
WellekZiegler::WellekZiegler(double & thcc, double & thco)
{	
	adjJointGenotypeCountsCase = new AdjustedJointGenotypeCounts();
	adjJointGenotypeCountsCon = new AdjustedJointGenotypeCounts();
	thresholdCC = calculateChiSqFromPvalue(thcc);
	thresholdCO = calculateChiSqFromPvalue(thco);
};

//! Calculates correlation and adjusted variance for Wellek Ziegler test.
void WellekZiegler::calculateRVariance(double & r, double & variance, AdjustedJointGenotypeCounts * jointGenotypeCounts)
{
	
	double total = jointGenotypeCounts->counts[0][0] + jointGenotypeCounts->counts[0][1] + jointGenotypeCounts->counts[0][2] +
					jointGenotypeCounts->counts[1][0] + jointGenotypeCounts->counts[1][1] + jointGenotypeCounts->counts[1][2] +
					jointGenotypeCounts->counts[2][0] + jointGenotypeCounts->counts[2][1] + jointGenotypeCounts->counts[2][2];

	//plug in the frequencies into the formula to calculate the correlation
	double sum1 = (double)(jointGenotypeCounts->counts[1][0] + jointGenotypeCounts->counts[1][1] + jointGenotypeCounts->counts[1][2])/total;
	double sum2 = (double)(jointGenotypeCounts->counts[2][0] + jointGenotypeCounts->counts[2][1] + jointGenotypeCounts->counts[2][2])/total;
	double sum3 = (double)(jointGenotypeCounts->counts[0][1] + jointGenotypeCounts->counts[1][1] + jointGenotypeCounts->counts[2][1])/total;
	double sum4 = (double)(jointGenotypeCounts->counts[0][2] + jointGenotypeCounts->counts[1][2] + jointGenotypeCounts->counts[2][2])/total;

	double e1 = sum1 + 2*sum2;
	double e2 = sum3 + 2*sum4;
	double c = (double)(jointGenotypeCounts->counts[1][1] + 2*(jointGenotypeCounts->counts[1][2] + jointGenotypeCounts->counts[2][1]) + 4*jointGenotypeCounts->counts[2][2])/total - e1*e2;
	double v1 = sum1 + 4*sum2 - e1*e1;
	double v2 = sum3 + 4*sum4 - e2*e2;

	//denomiator
	double d1 = sqrt(v1*v2);

	if(v1 == 0 || v2 == 0) r = 1; //avoid nans
	else r = c/d1; //correlation = r^2

	double g2 = 0.5*r; 

	double i1 = g2*(1-2*e1)/v1;
	double i2 = g2*4*(1-e1)/v1;
	double j1 = g2*(1-2*e2)/v2;
	double j2 = g2*4*(1-e2)/v2;

	double a0 = 1.0/d1;
	double a1 = e1/d1;
	double a2 = e2/d1;

	//double C00 = 0;
	double C01 = -a1                   - j1;
	double C02 = -2*a1                 - j2;
	double C10 = -a2                   - i1;
	double C11 = (a0-a2-a1)     - (i1 + j1);
	double C12 = (2*a0-a2-2*a1) - (i1 + j2);
	double C20 = -2*a2                 - i2;
	double C21 = (2*a0-2*a2-a1) - (i2 + j1);
	double C22 = 2*(2*a0-a2-a1) - (i2 + j2);

	double varianceBit1 = (C01*C01*jointGenotypeCounts->counts[0][1]
			 + C02*C02*jointGenotypeCounts->counts[0][2]
			 + C10*C10*jointGenotypeCounts->counts[1][0]
			 + C11*C11*jointGenotypeCounts->counts[1][1]
			 + C12*C12*jointGenotypeCounts->counts[1][2]
			 + C20*C20*jointGenotypeCounts->counts[2][0]
			 + C21*C21*jointGenotypeCounts->counts[2][1]
			 + C22*C22*jointGenotypeCounts->counts[2][2])/total;

	double varianceBit2 = (C01*jointGenotypeCounts->counts[0][1]
			 + C02*jointGenotypeCounts->counts[0][2]
			 + C10*jointGenotypeCounts->counts[1][0]
			 + C11*jointGenotypeCounts->counts[1][1]
			 + C12*jointGenotypeCounts->counts[1][2]
			 + C20*jointGenotypeCounts->counts[2][0]
			 + C21*jointGenotypeCounts->counts[2][1]
			 + C22*jointGenotypeCounts->counts[2][2])/total;

	variance = 	(varianceBit1 - varianceBit2*varianceBit2)/total; 
};

//! Calculates the statistic for the Wellek Ziegler test.
void WellekZiegler::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	//adjust the genotype counts for cases to account for zeroes
	adjJointGenotypeCountsCase->setAdjustedCounts(false, caseJointGenotypeCounts, false);

	//calculate log odds and adjusted variance for cases
	calculateRVariance(caseR, caseVar, adjJointGenotypeCountsCase);	

	chiSqStatCO = caseR*caseR/caseVar;

	if(calcCaseControl)
	{
		//adjust the genotype counts for controls to account for zeroes		
		adjJointGenotypeCountsCon->setAdjustedCounts(false, controlJointGenotypeCounts, false);

		//calculate log odds and adjusted variance for controls
		calculateRVariance(conR, conVar, adjJointGenotypeCountsCon);		
	
		double diff = caseR - conR; 
		chiSqStatCC = (diff*diff)/(caseVar + conVar);
	};

	if((calcCaseControl && chiSqStatCC >= thresholdCC) ||
	   (calcCaseOnly && !tooCloseForCaseOnly && chiSqStatCO >= thresholdCO))
	{
		thresholdPassed = true;
		totalThresholdPassed++;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	};
	
	totalCalculations++;
};

//! Outputs the header for the Wellek Ziegler test.
void WellekZiegler::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
		resultsFile << "WZ_CASE_R WZ_CASE_VAR ";
		if(calcCaseControl) resultsFile << "WZ_CTRL_R WZ_CTRL_VAR WZ_CC_CHISQ WZ_CC_P";	
		if(calcCaseControl && calcCaseOnly) resultsFile << " ";
		if(calcCaseOnly) resultsFile << "WZ_CO_CHISQ WZ_CO_P"; 
	};
};

//! Output results of the Wellek Ziegler test.
void WellekZiegler::outputResults(ofstream & resultsFile)
{
	if(!suppressResults) 
	{
		resultsFile << caseR << " " << caseVar << " ";
		if(calcCaseControl) resultsFile << conR << " " << conVar << " " << chiSqStatCC << " " << getPvalueChiSq1DF(chiSqStatCC);
		if(calcCaseControl && calcCaseOnly) resultsFile << " ";
		if(calcCaseOnly) resultsFile << chiSqStatCO << " " << getPvalueChiSq1DF(chiSqStatCO);	
	};
};

//! Outputs summary of results for Wellek Ziegler test.
void WellekZiegler::outputSummaryResults()
{
	out("Test Statistic: Wellek Ziegler\n");
	if(calcCaseControl) {out("P-value threshold for case/control results: "); out(getPvalueChiSq1DF(thresholdCC)); out("\n");};
	if(calcCaseOnly) {out("P-value threshold for case only results: "); out(getPvalueChiSq1DF(thresholdCO)); out("\n");};
	outputTotalCounts();
	out("\n");
};

//! Constructor with thresholds for R Squared statistic.
RSquared::RSquared()
{	
	adjJointGenotypeCountsCase = new AdjustedJointGenotypeCounts();
	adjJointGenotypeCountsCon = new AdjustedJointGenotypeCounts();
};

//! Calculates correlation and adjusted variance for Wellek Ziegler test.
void RSquared::calculateRSquared(double & rsq, AdjustedJointGenotypeCounts * jointGenotypeCounts)
{
	
	double total = jointGenotypeCounts->counts[0][0] + jointGenotypeCounts->counts[0][1] + jointGenotypeCounts->counts[0][2] +
					jointGenotypeCounts->counts[1][0] + jointGenotypeCounts->counts[1][1] + jointGenotypeCounts->counts[1][2] +
					jointGenotypeCounts->counts[2][0] + jointGenotypeCounts->counts[2][1] + jointGenotypeCounts->counts[2][2];

	//plug in the frequencies into the formula to calculate the correlation
	double sum1 = (double)(jointGenotypeCounts->counts[1][0] + jointGenotypeCounts->counts[1][1] + jointGenotypeCounts->counts[1][2])/total;
	double sum2 = (double)(jointGenotypeCounts->counts[2][0] + jointGenotypeCounts->counts[2][1] + jointGenotypeCounts->counts[2][2])/total;
	double sum3 = (double)(jointGenotypeCounts->counts[0][1] + jointGenotypeCounts->counts[1][1] + jointGenotypeCounts->counts[2][1])/total;
	double sum4 = (double)(jointGenotypeCounts->counts[0][2] + jointGenotypeCounts->counts[1][2] + jointGenotypeCounts->counts[2][2])/total;

	double e1 = sum1 + 2*sum2;
	double e2 = sum3 + 2*sum4;
	double c = (double)(jointGenotypeCounts->counts[1][1] + 2*(jointGenotypeCounts->counts[1][2] + jointGenotypeCounts->counts[2][1]) + 4*jointGenotypeCounts->counts[2][2])/total - e1*e2;
	double v1 = sum1 + 4*sum2 - e1*e1;
	double v2 = sum3 + 4*sum4 - e2*e2;

	//denomiator
	double d1 = sqrt(v1*v2);

	if(v1 == 0 || v2 == 0) rsq = 1; //avoid nans
	else rsq = c/d1; //correlation = r^2

	rsq = rsq*rsq;
};

//! Calculates the statistic for the Wellek Ziegler test.
void RSquared::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	//adjust the genotype counts for cases to account for zeroes
	adjJointGenotypeCountsCase->setAdjustedCounts(false, caseJointGenotypeCounts, false);

	//calculate log odds and adjusted variance for cases
	calculateRSquared(caseRSq, adjJointGenotypeCountsCase);	

	//adjust the genotype counts for controls to account for zeroes		
	adjJointGenotypeCountsCon->setAdjustedCounts(false, controlJointGenotypeCounts, false);

	//calculate log odds and adjusted variance for controls
	calculateRSquared(conRSq, adjJointGenotypeCountsCon);		
	
	totalCalculations++;
};

//! Outputs the header for the Wellek Ziegler test.
void RSquared::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
		resultsFile << "CASE_RSQ CTRL_RSQ";		
	};
};

//! Output results of the Wellek Ziegler test.
void RSquared::outputResults(ofstream & resultsFile)
{
	if(!suppressResults) 
	{
		resultsFile << caseRSq << " " << conRSq;	
	};
};

//! Outputs summary of results for Wellek Ziegler test.
void RSquared::outputSummaryResults()
{
	out("Statistic: R^2\n");
	outputTotalCounts();
	out("\n");
};

//! Outputs SNP pairs calculated only.
void RSquared::outputTotalCounts()
{
	out("Total SNP pairs calculated: "); out(totalCalculations); out("\n");	
};

//! Calculates the statistic for the Adjusted Wu test.
void DPrime::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	//estimate haplotype frequencies P11, P12, P21, P22 for cases
	estimateHaplotypeFreqsEM(caseJointGenotypeCounts, P11, P12, P21, P22);

	calculateDPrime(caseDPrime);
	
	//estimate haplotype frequencies P11, P12, P21, P22 for controls
	estimateHaplotypeFreqsEM(controlJointGenotypeCounts, P11, P12, P21, P22);

	calculateDPrime(controlDPrime);

	totalCalculations++;
};

//! Calculates the D Prime for the set estimated haplotypes.
void DPrime::calculateDPrime(double & aDPrime)
{
	D = P11*P22 - P12*P21;

	if(D < 0)
	{
		Dmax = (P11+P12)*(P11+P21);
		Dmax2 = (P21+P22)*(P12+P22);
		if(Dmax2 < Dmax) Dmax = Dmax2;
	}
	else
	{
		Dmax = (P11+P12)*(P12+P22);
		Dmax2 = (P11+P21)*(P21+P22);
		if(Dmax2 < Dmax) Dmax = Dmax2;
	};

	aDPrime = D/Dmax;
	if(aDPrime < 0) aDPrime = -aDPrime; //take absolute value of D', |D'|
};

//! Outputs the header for the Adjusted Wu test.
void DPrime::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
			resultsFile << "CASE_DPRIME CTRL_DPRIME";					
	};
};

//! Output results of the Adjusted Wu test.
void DPrime::outputResults(ofstream & resultsFile)
{
	if(!suppressResults)
	{
			if(caseDPrime*0 != 0) resultsFile << "NA ";
			else resultsFile << caseDPrime << " ";

			if(controlDPrime*0 != 0) resultsFile << "NA";
			else resultsFile << controlDPrime;
	};
};

//! Outputs summary of results for Adjusted Wu test.
void DPrime::outputSummaryResults()
{
	out("Statistic: D'\n");
	outputTotalCounts();
	out("\n");
};

//! Outputs SNP pairs calculated only.
void DPrime::outputTotalCounts()
{
	out("Total SNP pairs calculated: "); out(totalCalculations); out("\n");	
};

LogisticRegression::LogisticRegression(double & th) : LogLinearRegression()
{
	thresholdCC = calculateChiSqFromPvalue(th);
	model2SNPs = new TwoSNPLogRegModel();	
	useCovariates = false;
	covariateFile = "";
	covariatesToUseStr = "";
	missingQTValue = -9;
};

LinearRegression::LinearRegression(double & th) : LogLinearRegression()
{
	setThreshold(th);
	QuantitiveTraits * quantitiveTraits = new QuantitiveTraits();
	model2SNPs = new TwoSNPLinearRegModel(quantitiveTraits, missingQTValue);	
	useCovariates = false;
	covariateFile = "";
	covariatesToUseStr = "";
	missingQTValue = -9;
		
};

//! Sets threshold for linear regression.
void LinearRegression::setThreshold(double & th)
{
	if(th >= 1) pValThreshold = 1;
	else if(th < 0) pValThreshold = 0;
	else pValThreshold = th;
};

//! Calculates the statistic for the Logistic or Linear Regression test.
void LogLinearRegression::evaluateStatistic(bool & tooCloseForCaseOnly)
{
	fitModels();
	
	totalCalculations++;
};

//! Outputs the header for the Logistic Regression test.
void LogisticRegression::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
		if(useCovariates) resultsFile << "LR_COVAR_LOG_OR LR_COVAR_SE LR_COVAR_CHISQ LR_COVAR_P";
		else resultsFile << "LR_LOG_OR LR_SE LR_CHISQ LR_P";
	};
};

//! Output results of the Logistic Regression test.
void LogisticRegression::outputResults(ofstream & resultsFile)
{
	if(!suppressResults) 
	{
		if(fitStatusM3 == "Y" && fitStatusM4 =="Y")
		{		
			double beta = model2SNPs->getParameter(4);
			double SE;

			resultsFile << model2SNPs->getParameter(4) << " ";

			if(stat > 0)
			{
				if(beta > 0) SE = beta/sqrt(stat);
				else SE = -beta/sqrt(stat);

				resultsFile << SE << " ";
			}
			else
				resultsFile << "NA ";

			resultsFile << stat << " " << getPvalueChiSq1DF(stat);
		}
		//else if(fitStatusM3 == "D" || fitStatusM4 == "D") resultsFile << " NA NA NA";
		else resultsFile << "NA NA NA NA";		
	};
};

//! Outputs summary of results for Logistic Regression test.
void LogisticRegression::outputSummaryResults()
{
	out("Test Statistic: Logistic Regression\n");
	if(covariateFile != "")
	{
			out("Covariate file: "); out(covariateFile);
			if(covariatesToUseStr != "")
			{
				out(" (using covariate(s) \""); out(covariatesToUseStr); out("\" with missing value "); out(missingQTValue); out(")");
			}
			else
			{				
				out(" (using all covariates, with missing value "); out(missingQTValue); out(")");
			};
			out("\n");
	};
	out("P-value threshold for case/control results: "); out(getPvalueChiSq1DF(thresholdCC)); out("\n");
	outputTotalCounts();
	out("\n");
};

//! Calculates the SE for the interaction parameter.
double LinearRegression::calculateSE()
{
	//S.E. is given by sigma^ sqrt(g_{33}) where g_{33} is (3, 3) entry in (X^T X)^{-1}
	double sigma2 = rss0_1_2_3/(double)nMinuspf;

	double g_44;
	bool isOK = model2SNPs->calc_g44(g_44, useCovariates);

	if(isOK) return sqrt(sigma2*g_44);
	else return -1;
};

//! Outputs the header for the Linear Regression test.
void LinearRegression::outputHeader(ofstream & resultsFile)
{
	if(!suppressResults)
	{
		if(useCovariates) resultsFile << "LIN_COVAR_BETA LIN_COVAR_BETA_SE LIN_COVAR_FSTAT LIN_COVAR_P";
		else resultsFile << "LIN_BETA LIN_BETA_SE LIN_FSTAT LIN_P";
	};
};

//! Output results of the Linear Regression test.
void LinearRegression::outputResults(ofstream & resultsFile)
{
	if(!suppressResults) 
	{
		//Output F-stat for interaction test result
		//FSTAT = ((RSSr - RSSf)/(pf - pr))/RSSf/(n-pf)
		//RSSr residual sum of sqs reduced model, RSSf full model
		//pr, no of parameters reduced model, pf full model, FSTAT ~ F_{pf-pr, n-pf}
		if(fitStatusM3 == "Y" && fitStatusM4 =="Y")
		{
			if(rss0_1_2_3 == 0) if(rss0_1_2 > 0) resultsFile << "NA NA NA 0"; else resultsFile << "NA NA NA 1";
			else
			{	
				resultsFile << model2SNPs->getParameter(4) << " ";

				double SE = calculateSE();
				if(SE >= 0 ) resultsFile << SE << " "; else resultsFile << "NA ";

				resultsFile << fstat << " " << calculatedPval;
			};
		}
		//else if(fitStatusM3 == "D" || fitStatusM4 == "D") resultsFile << " NA NA NA NA";
		else resultsFile << "NA NA NA NA";
	};
};

//! Outputs summary of results for Linear Regression test.
void LinearRegression::outputSummaryResults()
{
	out("Test Statistic: Linear Regression\n");
	if(covariateFile != "")
	{
			out("Covariate file: "); out(covariateFile);
			if(covariatesToUseStr != "")
			{				
				out(" (using covariate(s) \""); out(covariatesToUseStr); out("\" with missing value "); out(missingQTValue); out(")");
			}
			else
			{				
				out(" (using all covariates, with missing value "); out(missingQTValue); out(")");
			};
			out("\n");
	};

	out("P-value threshold for results: "); out(pValThreshold); out("\n");
	outputTotalCounts();
	out("\n");	
};

//! Sets up initial parameters model M3.
void LogisticRegression::updateM3InitialParameters()
{
	if(useCovariates) updateCovarInitialParametersM3();
		
	initialParasM3[1] = threeParameters(prevM3Beta0, prevM3Beta1, prevM3Beta2);	
};

//! Sets up initial parameters model M4.
void LogisticRegression::updateM4InitialParameters()
{
	if(useCovariates) updateCovarInitialParametersM4();

	initialParasM4[1] = fourParameters(prevM3Beta0, prevM3Beta1, prevM3Beta2, prevM4Beta3);	
	initialParasM4[2] = fourParameters(prevM3Beta0, prevM3Beta1, prevM3Beta2, 0);
	initialParasM4[3] = fourParameters(prevM4Beta0, prevM4Beta1, prevM4Beta2, prevM4Beta3);
};

//! Sets up initial parameters model M3.
void LogisticRegression::updateCovarInitialParametersM3()
{
	initialParasCovarM3[1] = prevCovarM3;
};

//! Sets up initial parameters model M4.
void LogisticRegression::updateCovarInitialParametersM4()
{
	initialParasCovarM4[1] = prevCovarM3;
	initialParasCovarM4[2] = prevCovarM3;
	initialParasCovarM4[3] = prevCovarM4;	
};

//! Updates the previous covariate parameters for M3.
void LogisticRegression::updatePrevCovarParametersM3()
{
	unsigned int uc = 1; 
	for(list<double>::iterator pc = prevCovarM3.begin(); pc != prevCovarM3.end(); ++pc, ++uc) *pc = model2SNPs->getParameter(4+uc); 
};

//! Updates the previous covariate parameters for M4.
void LogisticRegression::updatePrevCovarParametersM4()
{
	unsigned int uc = 1; 
	for(list<double>::iterator pc = prevCovarM4.begin(); pc != prevCovarM4.end(); ++pc, ++uc) *pc = model2SNPs->getParameter(4+uc); 
};

//! Updates covariate parameters.
void LogisticRegression::updateCovarParameters(const unsigned int & parametersNo, const unsigned int & modelNo)
{
	unsigned int uc = 1;
	map<unsigned int, list<double> > initialParasCovar;
	map<unsigned int, list<double> >::const_iterator covars;

	if(modelNo == 3) initialParasCovar = initialParasCovarM3;
	else initialParasCovar = initialParasCovarM4;

	covars = initialParasCovar.find(parametersNo);
	
	if(covars != initialParasCovar.end())
	{
		for(list<double>::const_iterator ip = covars->second.begin(); ip != covars->second.end(); ++ip, ++uc)
		{
			model2SNPs->setParameter(4+uc, *ip);			
		};
	}
	else
	{
		covars = initialParasCovar.begin();

		//set covariates paras to 0 otherwise
		for(list<double>::const_iterator ip = covars->second.begin(); ip != covars->second.end(); ++ip, ++uc)
		{
			model2SNPs->setParameter(4+uc, 0);			
		};
	};			
};

//! Sets up initial parameters.
void LogisticRegression::setUpInitialParameters()
{
	//set initial start values of parameters
	//setup values for fitting models 
	prevM3Beta0 = 0; prevM3Beta1 = 0; prevM3Beta2 = 0;
	prevM4Beta0 = 0; prevM4Beta1 = 0; prevM4Beta2 = 0; prevM4Beta3 = 0;

	//setup starting values for fitting model M3	
	initialParasM3[2] = threeParameters(0, 0, 0);

	//setup starting values for fitting model M3	
	initialParasM4[4] = fourParameters(0, 0, 0, 0);

	//set covariate parameters
	if(useCovariates)
	{
		unsigned int noCovars = model2SNPs->getNoCovariates();
		for(unsigned int uc = 1; uc <= noCovars; ++uc)
		{					
			prevCovarM3.push_back(0);
			prevCovarM4.push_back(0);
		};
	};

};

//! Sets case-control status for LR and other initial stuff.
void LogisticRegression::setCaseControls(const list<unsigned char> & ccs)
{
	model2SNPs->setCaseControlsCovar(ccs);
	if(useCovariates) setupCovariates();
	setUpInitialParameters();
};

//! Sets case-control status for LR and other initial stuff.
void LinearRegression::setCaseControls(const list<unsigned char> & ccs)
{
	if(useCovariates) setupCovariates();
};

//! Sets up covariate data for LR.
void LogLinearRegression::setupCovariates()
{	
	CovariateData * covariateData = new CovariateData(covariateFile, covariatesToUseStr, filename, missingQTValue);
	model2SNPs->setCovariateData(covariateData, useCovariates);
	model2SNPs->updateCaseControlCovarWithMissing(covariateData->caseControls);
};

//! Fits logistic regression models.
void LogisticRegression::fitModels()
{
	//setup data with LR model
	if(useCovariates)
	{
		//update SNP data for fitting covariate models
		model2SNPs->setSNPData(window1->getCurrentSNP(), window2->getCurrentSNP());	
	}
	else
	{
		model2SNPs->setJointGenoTypeCounts(caseJointGenotypeCounts, controlJointGenotypeCounts);
	};

	FindFit findFit(model2SNPs);

	double negLogBeta0_1_2 = 0, negLogBeta0_1_2_3 = 0;
	map<unsigned int, double> parameters;
	set<unsigned int> parasToFit;

	//setup parameters
	parameters[1] = 0; //beta0
	parameters[2] = 0; //beta1	
	parameters[3] = 0; //beta2
	parameters[4] = 0; //beta3

	//set up covariate parameters, parameter 5 onwards
	if(useCovariates)
	{
		unsigned int noCovars = model2SNPs->getNoCovariates();
		for(unsigned int uc = 1; uc <= noCovars; ++uc)
		{
			parameters[4+uc] = 0; 
			parasToFit.insert(4+uc); 
		};
	};

	model2SNPs->setNewParameters(parameters);
	
	//check data is ok for fitting further models	
	bool anchorDataOK = true;
	bool partnerDataOK = true;
	bool anchorPartnerDiffOK = true;

	if(!useCovariates)
	{
		anchorDataOK = model2SNPs->checkAnchorSNPData();
		partnerDataOK = model2SNPs->checkPartnerSNPData();
		anchorPartnerDiffOK = model2SNPs->checkAnchorDiffPartnerSNPData();
	};

	parasToFit.insert(1); //beta0	
	parasToFit.insert(2); //beta1		
	parasToFit.insert(3); //beta2		

	fitStatusM3 = "Y";
	fitStatusM4 = "Y";

	bool fittedOK = false;

	//fit model M3
	fittedOK = false;

	if(anchorDataOK && partnerDataOK && anchorPartnerDiffOK)
	{
		updateM3InitialParameters();		

		//try fitting model M3 with different initial values
		for(map<unsigned int, threeParameters>::const_iterator par = initialParasM3.begin(); !fittedOK && par != initialParasM3.end(); ++par)
		{
			//set parameters
			model2SNPs->setParameter(1, par->second.beta0);			
			model2SNPs->setParameter(2, par->second.beta1);	
			model2SNPs->setParameter(3, par->second.beta2);

			//set covariate parameters			
			if(useCovariates) updateCovarParameters(par->first, 3);

			//fit the model	
			fittedOK = findFit.newtonsMethod(negLogBeta0_1_2, parasToFit);
			fittedOK = fittedOK && negLogBeta0_1_2*0 == 0;

		};		

		if(!fittedOK)
		{
			fitStatusM3 = "N";
		};		
	}
	else
	{
		fitStatusM3 = "D";
	};

	//if data is poor then a fit may not be found, and there will be no sign. improvement
	if(fittedOK) 
	{
		prevM3Beta0 = model2SNPs->getParameter(1);
		prevM3Beta1 = model2SNPs->getParameter(2);
		prevM3Beta2 = model2SNPs->getParameter(3);
		if(useCovariates) updatePrevCovarParametersM3();
	};

	if(!fittedOK)
	{
		fitStatusM4 = fitStatusM3;
	}
	else if(anchorDataOK && partnerDataOK && anchorPartnerDiffOK && (useCovariates || model2SNPs->checkInteractionSNPData()))
	{

		//fit model M4
		fittedOK = false;

		updateM4InitialParameters();
		
		parasToFit.insert(4); //beta3

		//try fitting model M4 with different initial values
		for(map<unsigned int, fourParameters >::const_iterator par = initialParasM4.begin(); !fittedOK && par != initialParasM4.end(); ++par)
		{
			//set parameters
			model2SNPs->setParameter(1, par->second.beta0);			
			model2SNPs->setParameter(2, par->second.beta1);	
			model2SNPs->setParameter(3, par->second.beta2);
			model2SNPs->setParameter(4, par->second.beta3);

			//set covariate parameters			
			if(useCovariates) updateCovarParameters(par->first, 4);

			//fit the model	
			fittedOK = findFit.newtonsMethod(negLogBeta0_1_2_3, parasToFit, 1e-6, 100);			
			fittedOK = fittedOK && negLogBeta0_1_2_3*0 == 0 && negLogBeta0_1_2_3 <= (negLogBeta0_1_2 + 1e-12);
		};	
	
		if(!fittedOK)
		{
			fitStatusM4 = "N";
		};

		//if data is poor then a fit may not be found, and there will be no sign. improvement
		if(!fittedOK) negLogBeta0_1_2_3 = negLogBeta0_1_2;
		else
		{
			prevM4Beta0 = model2SNPs->getParameter(1);
			prevM4Beta1 = model2SNPs->getParameter(2);
			prevM4Beta2 = model2SNPs->getParameter(3);	
			prevM4Beta3 = model2SNPs->getParameter(4);
			if(useCovariates) updatePrevCovarParametersM4();
		};
				
	}
	else
	{
		//no minor allele interaction data
		negLogBeta0_1_2_3 = negLogBeta0_1_2;
		fitStatusM4 = "D";
	};

	if(negLogBeta0_1_2_3 > negLogBeta0_1_2) negLogBeta0_1_2_3 = negLogBeta0_1_2;

	if(fitStatusM3 == "Y" && fitStatusM4 =="Y")
	{
		stat = 2*(negLogBeta0_1_2 - negLogBeta0_1_2_3);		
	}
	else stat = 0;

	if(calcCaseControl && stat >= thresholdCC)
	{
		thresholdPassed = true;
		totalThresholdPassed++;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	};
};

//! Fits linear regression models.
void LinearRegression::fitModels()
{
	//update SNP data
	model2SNPs->setSNPData(window1->getCurrentSNP(), window2->getCurrentSNP());

	fitStatusM3 = "Y";
	fitStatusM4 = "Y";	
	
	//fit M3 model
	if(!model2SNPs->fitModel(rss0_1_2, 1, 1, 0)) fitStatusM3 = "D";

	//fit M4 model
	if(!model2SNPs->fitModel(rss0_1_2_3, 1, 1, 1)) fitStatusM4 = "D";
	
	//check valid degrees of freedom
	unsigned int noMissSubjects = model2SNPs->getNoNonMissingSubjects();
	if(useCovariates)
	{
		unsigned int noCovariates = model2SNPs->getNoCovariates();
		if(noMissSubjects > 4 + noCovariates) nMinuspf = noMissSubjects - 4 - noCovariates;
		else fitStatusM4 = "N";
	}
	else
	{
		if(noMissSubjects > 4) nMinuspf = model2SNPs->getNoNonMissingSubjects() - 4;
		else fitStatusM4 = "N";
	};

	calculatedPval = 1;

	//calulate F-stat and p-value
	if(fitStatusM4 == "Y" && fitStatusM3 == "Y")
	{
		fstat = (rss0_1_2 - rss0_1_2_3)/(rss0_1_2_3/((double)nMinuspf));
		calculatedPval = getPvalueFStat(fstat, 1, nMinuspf);
	};

	//check if threshold passed
	if(calculatedPval <= pValThreshold || pValThreshold == 1)
	{
		thresholdPassed = true;
		totalThresholdPassed++;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	};
};
