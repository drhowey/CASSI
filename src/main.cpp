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


/*! \file main.cpp
    \brief This file reads in the initial input files and options.
    
    This file also outputs usage instructions and program details.
*/

#include <iostream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <string>

using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "Analysis.h"
#include "Statistics.h"

bool outputToScreen = true; 

ofstream logFile;
stringstream stringLogFile;

//! Output program title to screen
void header()
{
	out("\nCASSI: SNP interaction analysis software, v2.51\n");
	out("-----------------------------------------------\n");
	out("Copyright 2012-2017 Richard Howey, GNU General Public License, v3\n");
	out("Institute of Genetic Medicine, Newcastle University\n\n");
};

//! Output program usage to screen
void usage()
{
		header();
	 	
		out("Usage:\n\t ./cassi [options] file.bed\n");
		out(" or ./cassi parameterfile.pf [file.bed]\n\n");

		out("Options:\n");
		out("  -snp1 a1 a2   -- first SNP window, a1 = Start SNP no., a2 = End SNP no.\n");
		out("  -snp2 b1 b2   -- second SNP window, b1 = Start SNP no., b2 = End SNP no.\n");
		out("  -i file.bed   -- input file\n");
		out("  -i2 file2.bed -- second (optional) input file for second SNP window\n");
		out("  -o file.out   -- output file\n");
		out("  -log file.log -- output log file\n");	
		out("  -max m        -- maximum no. of results\n");			
		out("  -filter-all   -- all statistic thresholds must be met (as ordered)\n");	
		out("  -filter-any   -- any statistic threshold may be met\n");
		out("  -gap          -- gap in bp that SNPs need to be for case only tests\n");
		out("  -mem0         -- (slow) do not store SNPs in memory \n");
		out("  -mem1         -- (faster) store SNPs in memory, binary\n");
		out("  -mem2         -- (fastest) store SNPs in memory, integer\n");
		out("  -rsq          -- report R squared statistic between SNP pairs for cases and controls\n");
		out("  -dprime       -- report D' statistic between SNP pairs for cases and controls\n");
		out("  -so           -- suppress output to screen\n\n");

		out("Test Statistic Options:\n");
		out("  -je           -- do joint effects test\n");
		out("  -awu          -- do adjusted Wu test\n");
		out("  -wz           -- do Wellek Ziegler test\n");
		out("  -afe          -- do adjusted fast epistasis test\n");
		out("  -lr           -- do logistic regression test\n");
		out("  -lin          -- do linear regression test\n\n");

		out("Individual Test Statistic Options:\n");
		out("  -**-thcc t           -- p-value threshold for case/control test\n");
		out("  -**-thco t           -- p-value threshold for case only test\n");
		out("  -**-th t             -- p-value threshold for either test\n");
		out("  -**-cc-only          -- perform case/control test only\n");
		out("  -**-co-only          -- perform case only test only\n");
		out("  -**-so               -- suppress results of test\n");
		out("  -je-cellmin m        -- indiv. cell count minimum to do test\n");
		out("  -lr-covar covars.dat -- set covariates\n");		
		out("  -lr-covar-number no  -- covariate numbers, no\n");
		out("  -lr-covar-name na    -- covariate names, na\n");
		out("  -lr-covar-miss x     -- missing covariate value, x\n");
		out("  -lin-covar covars.dat-- set covariates\n");		
		out("  -lin-covar-number no -- covariate numbers, no\n");
		out("  -lin-covar-name na   -- covariate names, na\n");
		out("  -lin-covar-miss x    -- missing covariate value, x\n\n");

		out("Default Options in Effect:\n");	
		out("  -o filename cassi.out\n");
		out("  -je\n");	
		out("  -je-th 1e-4\n");
		out("  -je-cellmin 5\n");
		out("  -gap 1000\n");
		out("  -mem1\n");
		out("  -max 1e6\n\n");	
		
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(unsigned int & anUnInt, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{	
	if(useParaFile)
	{
		if(readParaFile.eof()) return;
		readParaFile >> anUnInt;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;				
		anUnInt = atoi(argv[argcount]);		
	};
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(double & aDouble, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aDouble;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;			
		aDouble = atof(argv[argcount]);		
	};
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(string & aString, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aString;
	}
	else
	{
		argcount++; if(argcount >= argc) return;		
		aString = argv[argcount];
	};
};


//! Gets the log filename from the results file name
string getDefaultLogFileName(string & outFileName)
{
	unsigned int filenameLength = outFileName.length();
	string logFileName;

	//find extension
	unsigned int a = filenameLength - 1;
	while(a > 0)
	{
		if(outFileName.substr(a, 1) == ".") break;
		a--;
	};

	if(a > 0) logFileName = outFileName.substr(0, a) + ".log";
	else  logFileName = outFileName + ".log";

	return logFileName;
};

//! The start of the program
int main(int argc, char * argv[])
{
	time_t start,end;
	double dif;
	time(&start);

	int argcount = 1;
	string option = "";
	string filename = "";
	string filename2 = "";
	string paraFilename = "";
	string outputFileName = "cassi.out";
	string logFilename = "";

	unsigned int snp1StartSNP = 1;
	unsigned int snp1EndSNP = 0;
	unsigned int snp2StartSNP = 1;
	unsigned int snp2EndSNP = 0;
	unsigned int caseOnlyGap = 1000;
	outputToScreen = true;	
	unsigned int maxNoResults = 1000000;
	bool useParaFile = false;
	ifstream readParaFile;
	list<Statistic *> statistics;	
	unsigned int filterType = 1;
	unsigned int memoryType = 1;

	//Joint Effects Parameters	
	double thresholdJECC = 1e-4;
	double thresholdJECO = 1e-4;
	double thresholdJE = -1;
	JointEffects * aJointEffects = new JointEffects();
	bool jointEffectsExist = false;
	unsigned int cellMin = 5;

	//Adjusted Wu Parameters	
	double thresholdAWUCC = 1e-4;
	double thresholdAWUCO = 1e-4;
	double thresholdAWU = -1;
	AdjustedWu * anAdjustedWu = new AdjustedWu();
	bool adjustedWuExist = false;

	//Wellek Ziegler Parameters	
	double thresholdWZCC = 1e-4;
	double thresholdWZCO = 1e-4;
	double thresholdWZ = -1;
	WellekZiegler * aWellekZiegler = new WellekZiegler();
	bool wellekZieglerExist = false;

	//Fast Epistasis Parameters
	double thresholdFECC = 1e-4;
	double thresholdFECO = 1e-4;
	double thresholdFE = 1e-4;	
	AdjustedFastEpistasis * aFastEpistasis = new AdjustedFastEpistasis();
	bool fastEpistasisExist = false;

	//Logistic Regression Parameters
	double thresholdLR = 1e-4;		
	LogisticRegression * aLogisticRegression = new LogisticRegression();
	bool logisticRegressionExist = false;
	string covariateFile, covariates;
	double missingQTValue = -9;

	//Linear Regression Parameters
	double thresholdLin = 1e-4;		
	LinearRegression * aLinearRegression = new LinearRegression();
	bool linearRegressionExist = false;
	string covariateFileLin, covariatesLin;
	double missingQTValueLin = -9;

	//R Squared Parameters	
	RSquared * aRSquared = new RSquared();
	bool rSquaredExist = false;

	//D Prime Parameters	
	DPrime * aDPrime = new DPrime();
	bool dPrimeExist = false;

	if(argcount < argc) option = argv[argcount];

	//deal with parameter file
	if(option == "-pf" || (option.length() >= 3 && option.substr(option.length()-3) == ".pf"))
	{
		if(option == "-pf")
		{
			argcount++; 
			if(argcount < argc) paraFilename = argv[argcount];
		}
		else paraFilename = option;

		//open parameter file		
		readParaFile.open(paraFilename.c_str());
		if(!readParaFile.is_open())
		{
			header();
			outErr("Cannot read parameter file: "); outErr(paraFilename); outErr("!\n");
			exit(0);
		};

		argcount++; 
		useParaFile = true;
	};

	
	//set given options
	while((!useParaFile && argcount < argc && argv[argcount][0] == '-') || (useParaFile && !readParaFile.eof()))
	{

		if(useParaFile)
		{
			//find the start of the next command
			do{
				readParaFile >> option;
				if(option.length() >= 2 && option.substr(0,1) == "-") break;
			}while(!readParaFile.eof());
		}
		else
		{
			option = argv[argcount];
		};

		if(useParaFile && readParaFile.eof()) break;

		if(option ==  "-snp1")
		{	
			getOptionValue(snp1StartSNP, useParaFile, argcount, argc, argv, readParaFile);
			getOptionValue(snp1EndSNP, useParaFile, argcount, argc, argv, readParaFile);						
		}
		else if(option ==  "-snp2")
		{	
			getOptionValue(snp2StartSNP, useParaFile, argcount, argc, argv, readParaFile);
			getOptionValue(snp2EndSNP, useParaFile, argcount, argc, argv, readParaFile);			
		}
		else if(option == "-i") getOptionValue(filename, useParaFile, argcount, argc, argv, readParaFile);					
		else if(option == "-i2") getOptionValue(filename2, useParaFile, argcount, argc, argv, readParaFile);		
		else if(option == "-o") getOptionValue(outputFileName, useParaFile, argcount, argc, argv, readParaFile);					
		else if(option == "-log") getOptionValue(logFilename, useParaFile, argcount, argc, argv, readParaFile);							
		else if(option == "-max") getOptionValue(maxNoResults, useParaFile, argcount, argc, argv, readParaFile);
		else if(option == "-gap") getOptionValue(caseOnlyGap, useParaFile, argcount, argc, argv, readParaFile);	
		else if(option == "-so") outputToScreen = false;		
		else if(option == "-filter-any") filterType = 2;
		else if(option == "-filter-all") filterType = 1;
		else if(option == "-mem0") memoryType = 0;
		else if(option == "-mem1") memoryType = 1;
		else if(option == "-mem2") memoryType = 2;
		else if(option == "-je") {aJointEffects = new JointEffects(thresholdJECC, thresholdJECO, cellMin); statistics.push_back(aJointEffects); jointEffectsExist = true;}
		else if(option == "-je-thcc") {getOptionValue(thresholdJECC, useParaFile, argcount, argc, argv, readParaFile); if(jointEffectsExist) aJointEffects->setThresholdCC(thresholdJECC);}	
		else if(option == "-je-thco") {getOptionValue(thresholdJECO, useParaFile, argcount, argc, argv, readParaFile); if(jointEffectsExist) aJointEffects->setThresholdCO(thresholdJECO);}	
		else if(option == "-je-th") {getOptionValue(thresholdJE, useParaFile, argcount, argc, argv, readParaFile); if(jointEffectsExist) {aJointEffects->setThresholdCC(thresholdJE); aJointEffects->setThresholdCO(thresholdJE);};}	
		else if(option == "-je-cc-only") {if(jointEffectsExist) aJointEffects->setCalcStatsCCandCO(true, false);}
		else if(option == "-je-co-only") {if(jointEffectsExist) aJointEffects->setCalcStatsCCandCO(false, true);}
		else if(option == "-je-so") {if(jointEffectsExist) aJointEffects->setSuppressResults(true);}	
		else if(option == "-je-cellmin") {getOptionValue(cellMin, useParaFile, argcount, argc, argv, readParaFile); if(jointEffectsExist) aJointEffects->setCellMin(cellMin);}
		else if(option == "-awu") {anAdjustedWu = new AdjustedWu(thresholdAWUCC, thresholdAWUCO); statistics.push_back(anAdjustedWu); adjustedWuExist = true;}				
		else if(option == "-awu-thcc") {getOptionValue(thresholdAWUCC, useParaFile, argcount, argc, argv, readParaFile); if(adjustedWuExist) anAdjustedWu->setThresholdCC(thresholdAWUCC);}	
		else if(option == "-awu-thco") {getOptionValue(thresholdAWUCO, useParaFile, argcount, argc, argv, readParaFile); if(adjustedWuExist) anAdjustedWu->setThresholdCO(thresholdAWUCO);}	
		else if(option == "-awu-th") {getOptionValue(thresholdAWU, useParaFile, argcount, argc, argv, readParaFile); if(adjustedWuExist) {anAdjustedWu->setThresholdCC(thresholdAWU); anAdjustedWu->setThresholdCO(thresholdAWU);};}
		else if(option == "-awu-cc-only") {if(adjustedWuExist) anAdjustedWu->setCalcStatsCCandCO(true, false);}
		else if(option == "-awu-co-only") {if(adjustedWuExist) anAdjustedWu->setCalcStatsCCandCO(false, true);}
		else if(option == "-awu-so") {if(adjustedWuExist) anAdjustedWu->setSuppressResults(true);}			
		else if(option == "-afe") {aFastEpistasis = new AdjustedFastEpistasis(thresholdFECC, thresholdFECO); statistics.push_back(aFastEpistasis); fastEpistasisExist = true;}				
		else if(option == "-afe-thcc") {getOptionValue(thresholdFECC, useParaFile, argcount, argc, argv, readParaFile); if(fastEpistasisExist) aFastEpistasis->setThresholdCC(thresholdFECC);}	
		else if(option == "-afe-thco") {getOptionValue(thresholdFECO, useParaFile, argcount, argc, argv, readParaFile); if(fastEpistasisExist) aFastEpistasis->setThresholdCO(thresholdFECO);}	
		else if(option == "-afe-th") {getOptionValue(thresholdFE, useParaFile, argcount, argc, argv, readParaFile); if(fastEpistasisExist) {aFastEpistasis->setThresholdCC(thresholdFE); aFastEpistasis->setThresholdCO(thresholdFE);};}
		else if(option == "-afe-cc-only") {if(fastEpistasisExist) aFastEpistasis->setCalcStatsCCandCO(true, false);}
		else if(option == "-afe-co-only") {if(fastEpistasisExist) aFastEpistasis->setCalcStatsCCandCO(false, true);}
		else if(option == "-afe-so") {if(fastEpistasisExist) aFastEpistasis->setSuppressResults(true);}
		else if(option == "-wz") {aWellekZiegler = new WellekZiegler(thresholdWZCC, thresholdWZCO); statistics.push_back(aWellekZiegler); wellekZieglerExist = true;}				
		else if(option == "-wz-thcc") {getOptionValue(thresholdWZCC, useParaFile, argcount, argc, argv, readParaFile); if(wellekZieglerExist) aWellekZiegler->setThresholdCC(thresholdWZCC);}	
		else if(option == "-wz-thco") {getOptionValue(thresholdWZCO, useParaFile, argcount, argc, argv, readParaFile); if(wellekZieglerExist) aWellekZiegler->setThresholdCO(thresholdWZCO);}	
		else if(option == "-wz-th") {getOptionValue(thresholdWZ, useParaFile, argcount, argc, argv, readParaFile); if(wellekZieglerExist) {aWellekZiegler->setThresholdCC(thresholdWZ); aWellekZiegler->setThresholdCO(thresholdWZ);};}
		else if(option == "-wz-cc-only") {if(wellekZieglerExist) aWellekZiegler->setCalcStatsCCandCO(true, false);}
		else if(option == "-wz-co-only") {if(wellekZieglerExist) aWellekZiegler->setCalcStatsCCandCO(false, true);}
		else if(option == "-wz-so") {if(wellekZieglerExist) aWellekZiegler->setSuppressResults(true);}
		else if(option == "-lr") {aLogisticRegression = new LogisticRegression(thresholdLR); statistics.push_back(aLogisticRegression); logisticRegressionExist = true;}				
		else if(option == "-lr-th" || option == "-lr-thcc") {getOptionValue(thresholdLR, useParaFile, argcount, argc, argv, readParaFile); if(logisticRegressionExist) aLogisticRegression->setThresholdCC(thresholdLR);}	
		else if(option == "-lr-covar") {getOptionValue(covariateFile, useParaFile, argcount, argc, argv, readParaFile); if(logisticRegressionExist) aLogisticRegression->setCovariateFile(covariateFile);}	
		else if(option == "-lr-covar-name" || option == "-lr-covar-number") {getOptionValue(covariates, useParaFile, argcount, argc, argv, readParaFile); if(logisticRegressionExist) aLogisticRegression->setCovariatesToUseStr(covariates);}
		else if(option == "-lr-covar-miss") {getOptionValue(missingQTValue, useParaFile, argcount, argc, argv, readParaFile); if(logisticRegressionExist) aLogisticRegression->setMissingQTValue(missingQTValue);}	
		else if(option == "-lr-so") {if(aLogisticRegression) aLogisticRegression->setSuppressResults(true);}
		else if(option == "-lin") {aLinearRegression = new LinearRegression(thresholdLin); statistics.push_back(aLinearRegression); linearRegressionExist = true;}				
		else if(option == "-lin-th" || option == "-lin-thcc") {getOptionValue(thresholdLin, useParaFile, argcount, argc, argv, readParaFile); if(linearRegressionExist) aLinearRegression->setThreshold(thresholdLin);}	
		else if(option == "-lin-covar") {getOptionValue(covariateFileLin, useParaFile, argcount, argc, argv, readParaFile); if(linearRegressionExist) aLinearRegression->setCovariateFile(covariateFileLin);}	
		else if(option == "-lin-covar-name" || option == "-lin-covar-number") {getOptionValue(covariatesLin, useParaFile, argcount, argc, argv, readParaFile); if(linearRegressionExist) aLinearRegression->setCovariatesToUseStr(covariatesLin);}
		else if(option == "-lin-covar-miss") {getOptionValue(missingQTValueLin, useParaFile, argcount, argc, argv, readParaFile); if(linearRegressionExist) aLinearRegression->setMissingQTValue(missingQTValueLin);}	
		else if(option == "-lin-so") {if(aLinearRegression) aLinearRegression->setSuppressResults(true);}
		else if(option == "-rsq") {aRSquared = new RSquared(); statistics.push_back(aRSquared); rSquaredExist = true;}
		else if(option == "-dprime") {aDPrime = new DPrime(); statistics.push_back(aDPrime); dPrimeExist = true;}
		else
		{
			if(logFilename != "") logFile.open(logFilename.c_str());
			else logFile.open(getDefaultLogFileName(outputFileName).c_str());

			header();
    		if(useParaFile) {outErr("Unrecognised option: "); outErr(option); outErr("\n\n");}
			else {outErr("Unrecognised command line switch: "); outErr(option); outErr("\n\n");};

			logFile << stringLogFile.str();
			logFile.close();
    		exit(0);			
		};

		if(!useParaFile) argcount++;
	};

	//setup the input files if not given by the "-i" options
	if(argcount < argc)
	{
		filename = argv[argcount];
		argcount++;
		if(argcount < argc) filename2 = argv[argcount];
	};

	//set up log filename
	if(logFilename == "") logFilename = getDefaultLogFileName(outputFileName);

	//set default test to perform as the joint test statistic
	if(statistics.size() == 0)
	{
		if(thresholdJE != -1) //over ride thresholds if overall threshold set
		{
			thresholdJECC = thresholdJE;
			thresholdJECO = thresholdJE;
		};

		aJointEffects = new JointEffects(thresholdJECC, thresholdJECO, cellMin);		
		statistics.push_back(aJointEffects);
		jointEffectsExist = true;
	};

	if(argc==1)
	{
		usage();
		exit(0);
	}
	else if(filename == "")
	{
		outErr("\nInput file not set!\n\n");
		usage();
		exit(0);
	}
	else if(filename.length() >= 4 && filename.substr(filename.length()-4) != ".bed" ||
		 filename2.length() >= 4 && filename2.substr(filename2.length()-4) != ".bed")
	{
		header();
		outErr("A binary pedigree file (.bed) is required for CASSI!\n\n");
		outErr("Try using PLINK with\n");
		outErr("plink --noweb --file mydata --make-bed\n\n");
		exit(0);
	};

	//setup filename for LR, needed to setup covariate data
	if(logisticRegressionExist) aLogisticRegression->setFilename(filename);
	if(linearRegressionExist) aLinearRegression->setFilename(filename);

	//set R square and D Prime whether to be included or not
	//if no test statistics then always include otherwise with any passed results
	if(dPrimeExist)
	{
		if(statistics.size() == 1) aDPrime->setThresholdPassed(true);
		else if(statistics.size() == 2 && rSquaredExist) aDPrime->setThresholdPassed(true);
		else {if(filterType == 1) aDPrime->setThresholdPassed(true); else aDPrime->setThresholdPassed(false);};
	};

	if(rSquaredExist)
	{
		if(statistics.size() == 1) aRSquared->setThresholdPassed(true);
		else if(statistics.size() == 2 && dPrimeExist) aRSquared->setThresholdPassed(true);
		else {if(filterType == 1) aRSquared->setThresholdPassed(true); else aRSquared->setThresholdPassed(false);};
	};

	//output options to screen
	header();

	if(filename == filename2)
	{
		out("WARNING: input file 1 is the same as input file 2, switching to one input file mode.\n\n");
		filename2 = "";
	};

	out("Parameters:\n");
	out("Input file: "); out(filename); out("\n");
	if(filename2 != "") {out("Second input file: "); out(filename2); out("\n");};
	out("Output file: "); out(outputFileName); out("\n");
	out("Log file: "); out(logFilename); out("\n");

	out("Start SNP of first SNP window: "); out(snp1StartSNP); out("\n");
	if(snp1EndSNP != 0) {out("End SNP of first SNP window: "); out(snp1EndSNP); out("\n");}
	else out("End SNP of first window is the last SNP\n");

	out("Start SNP of second SNP window: "); out(snp2StartSNP); out("\n");
	if(snp2EndSNP != 0) {out("End SNP of second SNP window: "); out(snp2EndSNP); out("\n");}
	else out("End SNP of second window is the last SNP\n");

	if(jointEffectsExist || fastEpistasisExist || wellekZieglerExist || adjustedWuExist) {out("Case only base pair gap: "); out(caseOnlyGap); out("\n");};

	if(memoryType == 0) out("Memory: SNPs not stored\n");
	else if(memoryType == 1) out("Memory: SNPs stored in binary format\n");
	else if(memoryType == 2) out("Memory: SNPs stored\n");

	if(maxNoResults == 0) out("Maximum no. of results: no maximum\n");
	else {out("Maximum no. of results: "); out(maxNoResults); out("\n");};
	
	if(statistics.size() > 1)
	{
		if(filterType == 1) out("Filter: all statistic thresholds must pass\n");
		else if(filterType == 2) out("Filter: any statistic threshold may pass\n");		
	};
	
	out("\n");
	
	//create analysis option and run analysis
	Analysis anAnalysis(filename, filename2, outputFileName, snp1StartSNP, snp1EndSNP, snp2StartSNP, snp2EndSNP, filterType, memoryType, caseOnlyGap, maxNoResults, statistics);

	anAnalysis.runAnalysis();
	
	time(&end);
	dif = difftime(end, start);
	out("\nRun time: "); out(getTime(dif)); out("\n\n");

	logFile.open(logFilename.c_str());
	logFile << stringLogFile.str();			
	logFile.close();
};

