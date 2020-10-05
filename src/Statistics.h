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


/*! \file Statistics.h
    \brief This file contains the models used for logistic regression.
    
*/

#ifndef __STATISTICS
#define __STATISTICS


#include <map>

using namespace std;

#include "Data.h"
#include "Model.h"
#include "ModelFitting.h"

// General Class for Statistics for interaction between two SNPs.
class Statistic
{
protected:

	unsigned int totalCalculations;
	unsigned int totalThresholdPassed;
	double thresholdCC;
	double thresholdCO;	
	bool thresholdPassed;
	bool calcCaseControl;
	bool calcCaseOnly;
	bool suppressResults;

	//Joint genotype data for calculating statistics
	JointGenotypeCounts * caseJointGenotypeCounts;
	JointGenotypeCounts * controlJointGenotypeCounts;
	
	//Pointers to SNP windows if individual SNP data is needed, eg for LR
	SNPWindow * window1;
	SNPWindow * window2;

public:

	Statistic()	: totalCalculations(0), totalThresholdPassed(0), thresholdCC(0), thresholdCO(0),  calcCaseControl(true), calcCaseOnly(true), suppressResults(false) 	
	{	
				
	};
	
	virtual ~Statistic()
	{
		
	};

	bool getThresholdPassed() const {return thresholdPassed;};
	virtual void setThresholdPassedFalse() {thresholdPassed = false;};
	virtual void outputTotalCounts();
	void setJointCounts(JointGenotypeCounts * cajc, JointGenotypeCounts * cojc) {caseJointGenotypeCounts = cajc; controlJointGenotypeCounts = cojc;};
	void setCalcStatsCCandCO(const bool & cc, const bool & co) {calcCaseControl = cc; calcCaseOnly = co;};
	void setSuppressResults(const bool & b) {suppressResults = b;};
	void addSNPWindows(SNPWindow * w1, SNPWindow * w2) {window1 = w1; window2 = w2;};

	//Required Methods For Statistics
	virtual void setCaseControls(const list<unsigned char> & ccs) {}; //statistic may need this if subjects need to bo considered indiv., eg LR with covariates
	virtual bool getCalcJointGenotypes() const {return true;};

	virtual void setThresholdCC(double & th); //by default set to a chi sq distribution with 1 degree of freedom
	virtual void setThresholdCO(double & th); //by default set to a chi sq distribution with 1 degree of freedom
	virtual void evaluateStatistic(bool & tooCloseForCaseOnly) {};
	virtual void outputHeader(ofstream & resultsFile) {};
	virtual void outputResults(ofstream & resultsFile) {};	
	virtual void outputSummaryResults() {};
};

//! Class that organises all the tests together.
class AllStatistics
{
private:
	
	ofstream resultsFile;
	list<Statistic *> statistics; //ordered list of test statistics to calculate
	unsigned int filterType; //which kind of filter to apply when there are multiple statistics to calculate
							 //1 = all, report all stats if all pass their thresholds; 
							 //2 = any, report all stats if any passes their threshold;							 
	unsigned int maxNoResults;
	unsigned int totalResults;
	DescriptionOfSNPs * descSNPs;

public:

	AllStatistics(string & resultsFilename, list<Statistic *> & sts, unsigned int & filt, unsigned int & mx, DescriptionOfSNPs * ds, JointGenotypeCounts * cajc, JointGenotypeCounts * cojc, list<unsigned char> & caco) : statistics(sts), filterType(filt), maxNoResults(mx), totalResults(0), descSNPs(ds)		
	{	
		for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
		{			
			(*s)->setJointCounts(cajc, cojc);
			(*s)->setCaseControls(caco);
		};
		//open results file
		resultsFile.open(resultsFilename.c_str());
		outputHeader();		
	};
	
	~AllStatistics()
	{
		resultsFile.close();
		for(list<Statistic *>::const_iterator s = statistics.begin(); s != statistics.end(); ++s)
		{			
			delete (*s);
		};
	};

	void outputTotalResults() {out("Number of SNP pairs with results: "); out(totalResults); out("\n");};

	void addSNPWindows(SNPWindow * window1, SNPWindow * window2);
	bool getCalcJointGenotypes() const;
	bool evaluateStatistics(bool & tooCloseForCaseOnly);
	void outputHeader();
	bool outputResults(unsigned int & snp1No, unsigned int & snp2No);
	void outputSummaries();
};

//! Adjusted Joint Genotype Counts
struct AdjustedJointGenotypeCounts
{	
	double counts[3][3];
	bool halfAdded;

	AdjustedJointGenotypeCounts()
	{
			
	};

	~AdjustedJointGenotypeCounts()
	{
			
	};

	void setAdjustedCounts(const bool & addHalf, JointGenotypeCounts * jointGenotypeCounts, const bool & adjustBaseline = true);
};

//! Joint Effects statistic by Ueki and Cordell.
class JointEffects : public Statistic
{
private:

	map<double, double> logCache;	

	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCase;
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCon;

	double lambdaAOverTotalInverseVca; //cases
	double muAOverTotalInverseVca;
	double lambdaNOverTotalInverseVco, muNOverTotalInverseVco; //controls
	double totalInverseVca, totalInverseVco;
	double caseControlStat;
	double caseOnlyStat;
	bool useAltStatCases;
	bool useAltStatCaseCons;
	unsigned int cellMin;
	bool cellCountTooSmallCase;
	bool cellCountTooSmallCon;

public:

	JointEffects() {};

	JointEffects(double & thjecc, double & thjeco, unsigned int & cm);
	
	~JointEffects()
	{
		delete adjJointGenotypeCountsCase;
		delete adjJointGenotypeCountsCon;	
	};

	//General Statistic Methods
	void evaluateStatistic(bool & tooCloseForCaseOnly);
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	
	//Methods for Joint Effects Test	
	void calculateis(double & i22, double & i21, double & i12, double & i11, AdjustedJointGenotypeCounts * ajgc);
	void calculateVariables(double & i22, double & i21, double & i12, double & i11, double & totalInverseV, double rowTotals[], AdjustedJointGenotypeCounts * ajc, const bool & useAltStat);
	double getOddRatioRelRisk(double & num1, double & num2, double & denom1, double & denom2);
	double getLog(double v);	
	void setCellMin(unsigned int & cm) {cellMin = cm;};
};

//! Fast Epistasis statistic as Ueki and Cordell.
class AdjustedFastEpistasis : public Statistic
{
private:
	
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCase;
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCon;

	double chiSqStatCC, chiSqStatCO;	

	double caseVar, conVar;
	double caseLogOdds, conLogOdds;

public:

	AdjustedFastEpistasis() {};

	AdjustedFastEpistasis(double & thcc, double & thco);
	
	~AdjustedFastEpistasis()
	{
		delete adjJointGenotypeCountsCase;
		delete adjJointGenotypeCountsCon;		
	};

	//General Statistic Methods
	void evaluateStatistic(bool & tooCloseForCaseOnly);
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	
	//Methods for Fast Epistasis Test	
	void calculateLogOddsAdjustedVariance(double & logOdds, double & variance,  AdjustedJointGenotypeCounts * jointGenotypeCounts);		
};

//! Adjusted Wu statistic as descibed by Ueki and Cordell
class AdjustedWu : public Statistic
{
private:

	double P11, P12, P21, P22, D, p, q, u, v, n;
	double caP11, caP12, caP21, caP22;
	double lambdaA, lambdaN; //log odds cases and controls
	double varianceCases, varianceControls;

	double caseControlStat;
	double caseOnlyStat;
	
public:

	AdjustedWu() {};

	AdjustedWu(double & thcc, double & thco);
	
	~AdjustedWu()
	{
			
	};

	//General Statistic Methods
	void evaluateStatistic(bool & tooCloseForCaseOnly);
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	
	//Methods for Adjusted Wu Test
	void calculateLogORVar(double & logOdds, double & variance);		
};


//! Wellek Ziegler inspired statistic based on correlation coefficient as decribed by Ueki and Cordell.
class WellekZiegler : public Statistic
{
private:
	
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCase;
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCon;

	double chiSqStatCC, chiSqStatCO;	

	double caseVar, conVar;
	double caseR, conR;


public:

	WellekZiegler() {};

	WellekZiegler(double & thcc, double & thco);
	
	~WellekZiegler()
	{
		delete adjJointGenotypeCountsCase;
		delete adjJointGenotypeCountsCon;	
	};

	//General Statistic Methods
	void evaluateStatistic(bool & tooCloseForCaseOnly);
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	
	//Methods for Wellek Ziegler Test	
	void calculateRVariance(double & r, double & variance, AdjustedJointGenotypeCounts * jointGenotypeCounts);
		
};

//! Calculate Dprime for the cases and controls
class DPrime : public Statistic
{
private:

	double P11, P12, P21, P22, D, Dmax, Dmax2;
	
	double caseDPrime;
	double controlDPrime;
	
public:

	DPrime() {};
	
	~DPrime()
	{
			
	};

	//General Statistic Methods
	void outputTotalCounts();
	void evaluateStatistic(bool & tooCloseForCaseOnly);
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	
	//Methods for Dprime
	void calculateDPrime(double & aDPrime);
	void setThresholdPassed(const bool & p) {thresholdPassed = p;};
	void setThresholdPassedFalse() {}; //do not reset this, as not realy a test
};

//! Outputs the Rsquared for cases and controls.
class RSquared : public Statistic
{
private:
	
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCase;
	AdjustedJointGenotypeCounts * adjJointGenotypeCountsCon;

	double caseRSq, conRSq;

public:

	RSquared();
	
	~RSquared()
	{
		delete adjJointGenotypeCountsCase;
		delete adjJointGenotypeCountsCon;	
	};

	//General Statistic Methods
	void outputTotalCounts();
	void evaluateStatistic(bool & tooCloseForCaseOnly);
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	
	//Methods for RSquared Test	
	void calculateRSquared(double & rsq, AdjustedJointGenotypeCounts * jointGenotypeCounts);
	void setThresholdPassed(const bool & p) {thresholdPassed = p;};	
	void setThresholdPassedFalse() {}; //do not reset this, as not realy a test
};

//! Class for storing three parameters for LR model fitting.
struct threeParameters
{
	double beta0, beta1, beta2;

	threeParameters() : beta0(0), beta1(0), beta2(0) {};
	threeParameters(const double & v1, const double & v2, const double & v3) : beta0(v1), beta1(v2), beta2(v3) {};
};

//! Class for storing four parameters for LR model fitting.
struct fourParameters
{
	double beta0, beta1, beta2, beta3;

	fourParameters() : beta0(0), beta1(0), beta2(0), beta3(0) {};
	fourParameters(const double & v1, const double & v2, const double & v3, const double & v4) : beta0(v1), beta1(v2), beta2(v3), beta3(v4) {};
};

//! Logistic or Linear Regression.
class LogLinearRegression : public Statistic
{
protected:
	
	double stat;	//chi Sq Stat for logistic and F stat for linear

	bool useCovariates;
	string covariateFile;
	string filename;
	string covariatesToUseStr;
	double missingQTValue;
	Model2SNPs * model2SNPs;

	string fitStatusM3;
	string fitStatusM4;

public:

	LogLinearRegression() : Statistic() {};

	~LogLinearRegression()
	{
		delete model2SNPs;	
	};

	//General Statistic Methods
	void evaluateStatistic(bool & tooCloseForCaseOnly);

	//modified general method
	virtual void setCaseControls(const list<unsigned char> & ccs) {};
	bool getCalcJointGenotypes() const {return !useCovariates;};

	//Methods for Logistic or Linear Regression Test	
	virtual void fitModels() {};
	virtual void setFilename(string & f) {filename = f;};
	void setupCovariates();
	void setCovariateFile(string & cf) {covariateFile = cf; useCovariates = true;};
	void setCovariatesToUseStr(string & ctu) {covariatesToUseStr = ctu;};
	void setMissingQTValue(double & q) {missingQTValue = q;};
};

//! Logistic Regression.
class LogisticRegression : public LogLinearRegression
{
private:
	
	//used for fitting models, previous fitted parameters used as starting values
	double prevM3Beta0, prevM3Beta1, prevM3Beta2;
	double prevM4Beta0, prevM4Beta1, prevM4Beta2, prevM4Beta3;
	list<double> prevCovarM3;
	list<double> prevCovarM4;
	map<unsigned int, threeParameters > initialParasM3;
	map<unsigned int, fourParameters > initialParasM4;
	map<unsigned int, list<double> > initialParasCovarM3;
	map<unsigned int, list<double> > initialParasCovarM4;

public:

	LogisticRegression() {};

	LogisticRegression(double & th);
	
	~LogisticRegression()
	{
		
	};

	//General Statistic Methods
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	void setCaseControls(const list<unsigned char> & ccs);

	//Methods for Logistic Regression Test	
	void fitModels();
	void updateM3InitialParameters();	
	void updateCovarInitialParametersM3();
	void updateM4InitialParameters();	
	void setUpInitialParameters();
	void updateCovarInitialParametersM4();
	void updateCovarParameters(const unsigned int & parametersNo, const unsigned int & modelNo);
	void updatePrevCovarParametersM3();
	void updatePrevCovarParametersM4();
};

//! Linear Regression.
class LinearRegression : public LogLinearRegression
{
private:
	double pValThreshold;
	double fstat;
	double calculatedPval;
	double rss0_1_2, rss0_1_2_3;
	unsigned int nMinuspf;

public:

	LinearRegression() {};

	LinearRegression(double & th);
	
	~LinearRegression()
	{
		
	};

	//General Statistic Methods
	void outputHeader(ofstream & resultsFile);
	void outputResults(ofstream & resultsFile);	
	void outputSummaryResults();
	void setCaseControls(const list<unsigned char> & ccs);

	//Methods for Linear Regression Test	
	void fitModels();
	void setFilename(string & f) {filename = f; model2SNPs->setupQuantitiveTraits(filename);};
	double calculateSE();
	void setThreshold(double & th);
};

#endif
