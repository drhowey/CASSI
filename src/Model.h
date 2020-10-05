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


/*! \file Model.h
    \brief This file contains the models used for logistic regression.
    
*/

#ifndef __MODEL
#define __MODEL


#include <map>

using namespace std;

struct JointGenotypeCounts;
struct SNPData;
struct QuantitiveTraits;
struct CovariateData;

#include "Data.h"

//! General Class for Models.
class Model
{
protected:

	map<unsigned int, double> parameters;
	
	//used for fitting covariates
	bool useCovariates;
	CovariateData * covariateData;	
	SNPData * anchorSNP;
	SNPData * partnerSNP;
	list<unsigned char> caseControlsCovar;

public:

	Model() : parameters(), useCovariates(false) {};
	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
	Model(map<unsigned int, double> paras) : parameters(paras), useCovariates(false) {};
	
	virtual ~Model() {};

	double getParameter(const unsigned int & no) const {return parameters.find(no)->second;};
	map<unsigned int, double> getParameters() const {return parameters;};
	void setNewParameters(map<unsigned int, double> & paras);
	void setCovariateData(CovariateData * co, bool & uc) {covariateData = co; useCovariates = uc;};
	void setSNPData(SNPData * s1, SNPData * s2) {anchorSNP = s1; partnerSNP = s2;};	
	void setCaseControlsCovar(const list<unsigned char> & ccs) {caseControlsCovar = ccs;};
	void updateCaseControlCovarWithMissing(list<unsigned char> & caseControls);
	unsigned int getNoCovariates() {return covariateData->covariateDataAllSubjects.begin()->size();};

	virtual void setParameter(const unsigned int & no, const double & value) {parameters[no] = value;};	
	virtual double negLogLikelihood() {return 0;};
	virtual void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta3) const {};
	virtual void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta3) const {};

};

//! General Class for Models with 2 SNPs.
class Model2SNPs : public Model
{
protected:
	
	double beta0, beta1, beta2, beta3; // base, SNP1 effect, SNP2 effect, SNP1xSNP2 effect
	
public:

	Model2SNPs() {};
	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
		
	virtual ~Model2SNPs() {};

	void setParameter(const unsigned int & no, const double & value);	
	void setPartnerData(SNPData * par) {partnerSNP = par;};	

	virtual void setJointGenoTypeCounts(JointGenotypeCounts * gjcca, JointGenotypeCounts * gjcco) {};
	
	virtual unsigned int getNoNonMissingSubjects() const {return 0;};
	virtual bool checkAnchorSNPData() const {return true;};
	virtual bool checkPartnerSNPData() const {return true;};
	virtual bool checkInteractionSNPData() const {return true;};
	virtual bool checkAnchorDiffPartnerSNPData() const {return true;};
	virtual bool fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & fitInter) {return 0;};
	virtual void setupQuantitiveTraits(string & filename) {};
	virtual bool calc_g44(double & g44, bool & useCovariates) {return true;};
};

//! Two SNP logistic regression model.
class TwoSNPLogRegModel : public Model2SNPs
{
private:

	JointGenotypeCounts * partnerJointGenotypeCountsCases; //points to the same partner object in the SNP window object
	JointGenotypeCounts * partnerJointGenotypeCountsControls;

protected:

	//varibles used for fitting the model
	mutable double effects;
	mutable double expEffects;
	mutable double noCases, noControls;
	mutable double aNumber;	
	mutable list<double> expEffectsAllIndiv;

public:

	TwoSNPLogRegModel() {};
	
	~TwoSNPLogRegModel() {};
	
	void setJointGenoTypeCounts(JointGenotypeCounts * gjcca, JointGenotypeCounts * gjcco) {partnerJointGenotypeCountsCases = gjcca; partnerJointGenotypeCountsControls = gjcco;};

	double negLogLikelihood();
	
	void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta3) const;
	void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta3) const;

	bool checkAnchorSNPData() const;
	bool checkPartnerSNPData() const;
	bool checkAnchorDiffPartnerSNPData() const;
	bool checkInteractionSNPData() const;

	double negLogLikelihoodCovar();
	void getGradientVectorCovar(map<unsigned int, double> & gradientVector, const bool & fitBeta3) const;
	void getHessianMatrixCovar(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta3) const;
};

//! Two SNP linear regression model.
class TwoSNPLinearRegModel : public Model2SNPs
{
private:

	QuantitiveTraits * quantitiveTraits;
	double missingQTValue;
	unsigned int totalNotMissingBetweenSNPs;

public:

	TwoSNPLinearRegModel() {};
	
	TwoSNPLinearRegModel(QuantitiveTraits * qt, double & mQTV) : quantitiveTraits(qt), missingQTValue(mQTV) {};
	
	~TwoSNPLinearRegModel()
	{
		delete quantitiveTraits;
	};

	unsigned int getNoNonMissingSubjects() const {return totalNotMissingBetweenSNPs;};
	
	void setupQuantitiveTraits(string & filename) {quantitiveTraits->setupQuantitiveTraits(filename, missingQTValue);};
	bool fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & fitInter); 
	bool fitModelCovar(double & rss, const bool & fitAnc, const bool & fitPart, const bool & fitInter);
	bool calc_g44(double & g44, bool & useCovariates);
	bool calc_g44Covar(double & g44);
};

#endif
