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



/*! \file Model.cpp
    \brief This file contains the source of models used for logistic regression.
    
*/

#include <map>
#include <iostream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Model.h"
#include "Data.h"
#include "ModelFitting.h"
#include "Utils.h"

//! Sets the values of all of the parameters.
void Model::setNewParameters(map<unsigned int, double> & paras)
{
	for(map<unsigned int, double>::const_iterator p = paras.begin(); p != paras.end(); ++p)
	{
		setParameter(p->first, p->second);
	};
};

//!Updates caseControl status with more missing indivs (with missing covariates)
void Model::updateCaseControlCovarWithMissing(list<unsigned char> & otherCaseControls)
{
	list<unsigned char>::iterator cc = caseControlsCovar.begin();

	if(otherCaseControls.size() > 0)
	{
		for(list<unsigned char>::const_iterator oc = otherCaseControls.begin(); oc != otherCaseControls.end(); ++oc, ++cc)
		{	
			if(*oc == 0) *cc = 0;			
		};
	};
};

//! Returns the value of the negative log likelihood for the given parameters and joint genotype counts.
double TwoSNPLogRegModel::negLogLikelihood()
{
	if(useCovariates) return negLogLikelihoodCovar();
	double ans = 0;
	
	//loop thro' the number of minor alleles for SNP1, i.e. start SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
		//loop thro' the number of minor alleles for SNP2, i.e. partner SNP
		for(unsigned int snp2 = 0; snp2 <= 2; ++snp2)
		{			
			effects = beta0 + snp1*beta1 + snp2*beta2 + snp1*snp2*beta3;
			noCases = (double)partnerJointGenotypeCountsCases->counts[snp1][snp2];
			noControls = (double)partnerJointGenotypeCountsControls->counts[snp1][snp2];
		
			ans += -noCases*effects + (noCases+noControls)*log(1 + exp(effects));
		};
	};
	
	return ans;
};

//! Returns the value of the negative log likelihood for the given parameters and joint genotype counts.
double TwoSNPLogRegModel::negLogLikelihoodCovar()
{
	double ans = 0;	

	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControlsCovar.begin();

	unsigned int covariateParameterNo;

	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++part, ++covar, ++cc)
	{
		if(*anc != 3 && *part != 3 && *cc != 0) //skip missing data
		{
			effects = beta0 + (*anc)*beta1 + (*part)*beta2 + (*anc)*(*part)*beta3;

			covariateParameterNo = 5;
			//add covariates effects
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			{
				effects += (*cv)*parameters[covariateParameterNo];
				covariateParameterNo++;
			};

			if(*cc == 2) ans -= effects;
			ans += log(1 + exp(effects));
		};
	};

	return ans;
};

//! Sets the value of a given parameter for the model.
void Model2SNPs::setParameter(const unsigned int & no, const double & value)
{
	parameters[no] = value;

	switch(no)
	{
		case 1:
			beta0 = value; return;
		case 2:
			beta1 = value; return;
		case 3:
			beta2 = value; return;
		case 4:
			beta3 = value; return;			
	};
};	

//! Checks if there are any minor alleles for the anchor SNP. 
bool TwoSNPLogRegModel::checkAnchorSNPData() const
{
	return (partnerJointGenotypeCountsCases->counts[1][0] > 0 || partnerJointGenotypeCountsCases->counts[1][1] > 0 || partnerJointGenotypeCountsCases->counts[1][2] > 0 ||
	partnerJointGenotypeCountsCases->counts[2][0] > 0 || partnerJointGenotypeCountsCases->counts[2][1] > 0 || partnerJointGenotypeCountsCases->counts[2][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[1][0] > 0 || partnerJointGenotypeCountsControls->counts[1][1] > 0 || partnerJointGenotypeCountsControls->counts[1][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[2][0] > 0 || partnerJointGenotypeCountsControls->counts[2][1] > 0 || partnerJointGenotypeCountsControls->counts[2][2] > 0 );
};

//! Checks if there are any minor alleles for the partner SNP. 
bool TwoSNPLogRegModel::checkPartnerSNPData() const
{
	 return (partnerJointGenotypeCountsCases->counts[0][1] > 0 || partnerJointGenotypeCountsCases->counts[0][2] > 0 ||
		 partnerJointGenotypeCountsCases->counts[1][1] > 0 || partnerJointGenotypeCountsCases->counts[1][2] > 0 ||
		 partnerJointGenotypeCountsCases->counts[2][1] > 0 || partnerJointGenotypeCountsCases->counts[2][2] > 0 ||
		 partnerJointGenotypeCountsControls->counts[0][1] > 0 || partnerJointGenotypeCountsControls->counts[0][2] > 0 ||
		 partnerJointGenotypeCountsControls->counts[1][1] > 0 || partnerJointGenotypeCountsControls->counts[1][2] > 0 ||
		 partnerJointGenotypeCountsControls->counts[2][1] > 0 || partnerJointGenotypeCountsControls->counts[2][2] > 0);
};

//! Checks if there are any minor allele data for the anchor and partner is different.
bool TwoSNPLogRegModel::checkAnchorDiffPartnerSNPData() const
{	
	return (partnerJointGenotypeCountsCases->counts[0][1] > 0 || partnerJointGenotypeCountsCases->counts[0][2] > 0 ||
	partnerJointGenotypeCountsCases->counts[1][0] > 0 || partnerJointGenotypeCountsCases->counts[1][2] > 0 ||
	partnerJointGenotypeCountsCases->counts[2][0] > 0 || partnerJointGenotypeCountsCases->counts[2][1] > 0 ||
	partnerJointGenotypeCountsControls->counts[0][1] > 0 || partnerJointGenotypeCountsControls->counts[0][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[1][0] > 0 || partnerJointGenotypeCountsControls->counts[1][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[2][0] > 0 || partnerJointGenotypeCountsControls->counts[2][1] > 0);
};

//! Checks if there are any minor allele interactions for the anchor partner SNP.
bool TwoSNPLogRegModel::checkInteractionSNPData() const
{	
	return (partnerJointGenotypeCountsCases->counts[1][1] + partnerJointGenotypeCountsCases->counts[1][2] +
		partnerJointGenotypeCountsCases->counts[2][1] + partnerJointGenotypeCountsCases->counts[2][2] +
		partnerJointGenotypeCountsControls->counts[1][1] + partnerJointGenotypeCountsControls->counts[1][2] +
		partnerJointGenotypeCountsControls->counts[2][1] + partnerJointGenotypeCountsControls->counts[2][2] != 0);			
};

//! Calculates the gradient vector used for Newton's Method.
void TwoSNPLogRegModel::getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta3)  const
{
	if(useCovariates) {getGradientVectorCovar(gradientVector, fitBeta3); return;};

	double ans[4] = {0, 0, 0 ,0};

	//loop thro' the number of minor alleles for SNP1, i.e. start SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
		//loop thro' the number of minor alleles for SNP2, i.e. partner SNP
		for(unsigned int snp2 = 0; snp2 <= 2; ++snp2)
		{			
			effects = beta0 + snp1*beta1 + snp2*beta2 + snp1*snp2*beta3;
			noCases = (double)partnerJointGenotypeCountsCases->counts[snp1][snp2];
			noControls = (double)partnerJointGenotypeCountsControls->counts[snp1][snp2];
			aNumber = -noCases + (noCases + noControls)/(exp(-effects) + 1);
			
			ans[0] += aNumber;			
			ans[1] += snp1*aNumber;
			ans[2] += snp2*aNumber;			
			if(fitBeta3) ans[3] += snp1*snp2*aNumber;				
		};

	};

	gradientVector[1] = ans[0];
	gradientVector[2] = ans[1];
	gradientVector[3] = ans[2];	
	if(fitBeta3) gradientVector[4] = ans[3];	
};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
void TwoSNPLogRegModel::getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta3) const
{
	if(useCovariates) {getHessianMatrixCovar(hessianMatrix, fitBeta3); return;};

	double ans[4][4] = {{0, 0, 0 ,0}, {0, 0, 0 ,0}, {0, 0, 0 ,0}, {0, 0, 0 ,0}}; //col, row
	
	//loop thro' the number of minor alleles for SNP1, i.e. anchor SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
		//loop thro' the number of minor alleles for SNP2, i.e. partner SNP
		for(unsigned int snp2 = 0; snp2 <= 2; ++snp2)
		{			
			effects = beta0 + snp1*beta1 + snp2*beta2 + snp1*snp2*beta3;
			expEffects = exp(effects);
			noCases = (double)partnerJointGenotypeCountsCases->counts[snp1][snp2];
			noControls = (double)partnerJointGenotypeCountsControls->counts[snp1][snp2];
			aNumber = ((noCases+noControls)*expEffects)/((expEffects + 1)*(expEffects + 1));

			//only calculate for the parameters we are fitting over
			ans[0][0] += aNumber;			
			ans[0][1] += snp1*aNumber;	
			ans[1][1] += snp1*snp1*aNumber;			
			ans[0][2] += snp2*aNumber;
			ans[2][2] += snp2*snp2*aNumber;
			ans[1][2] += snp1*snp2*aNumber;	
			
			if(fitBeta3)
			{
				ans[0][3] += snp1*snp2*aNumber;								
				ans[3][3] += snp1*snp2*snp1*snp2*aNumber;
				ans[1][3] += snp1*snp1*snp2*aNumber;			
				ans[2][3] += snp2*snp1*snp2*aNumber;	
			};
			
		};

	};

	//The Hessian matrix is symetric so only calculate half and then copy
	ans[1][0] = ans[0][1];		
	ans[2][0] = ans[0][2];
	ans[2][1] = ans[1][2];
	
	if(fitBeta3)
	{
		ans[3][0] = ans[0][3];
		ans[3][1] = ans[1][3];			
		ans[3][2] = ans[2][3];	
	};

	//setup the matrix with calculated values
	map<unsigned int, double> aCol;
	unsigned int noParas = 3;

	if(fitBeta3) noParas = 4;
	
	for(unsigned int col = 0; col < noParas; ++col)
	{		
		for(unsigned int row = 0; row < noParas; ++row)
		{
			aCol[row + 1] = ans[col][row];	
		};
		hessianMatrix[col + 1] = aCol;
	};

};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood with covariates.
void TwoSNPLogRegModel::getHessianMatrixCovar(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta3) const
{
	double ans[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}; //col, row
	map<pair<unsigned int, unsigned int>, double> covarAns; //starting count from 0, covariates start after other parameters, either 3 or 4

	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControlsCovar.begin();
	list<double>::const_iterator cv1, cv2;

	//unsigned int covariateParameterNo;
	unsigned int cp = 3; //first parameter no. for covariate in covarAns matrix
	if(fitBeta3) cp++;
	
	double val1;

	unsigned int noParas = cp + covar->size();

	list<double>::const_iterator expEffs = expEffectsAllIndiv.begin();//set up in Vector

	//loop thro' individuals
	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++part, ++covar, ++cc)
	{
		if(*anc != 3 && *part != 3 && *cc != 0) //skip missing data
		{						

			aNumber = (*expEffs)/((*expEffs + 1)*(*expEffs + 1));
			++expEffs;

			//only calculate for the parameters we are fitting over
			ans[0][0] += aNumber;			
			ans[0][1] += (*anc)*aNumber;	
			ans[1][1] += (*anc)*(*anc)*aNumber;				
			ans[0][2] += (*part)*aNumber;
			ans[2][2] += (*part)*(*part)*aNumber;
			ans[1][2] += (*anc)*(*part)*aNumber;

			if(fitBeta3)
			{
				val1 = (*anc)*(*part);
				ans[0][3] += val1*aNumber;								
				ans[3][3] += val1*val1*aNumber;
				ans[1][3] += (*anc)*val1*aNumber;			
				ans[2][3] += (*part)*val1*aNumber;	
			};			

			cv1 = covar->begin();

			//add covariate 2nd derivs
			for(unsigned int i1 = 0; i1 < noParas; )
			{
				cv2 = covar->begin();

				for(unsigned int i2 = 0; i2 < noParas; )
				{					
					//update any bit of matrix that is part of covariate in top half of matrix 
					if(i2 >= cp && i2 >= i1)
					{
						if(i1 < cp) //covariate and snp parameter 2nd deriv
						{
							if(i1 == 0) val1 = 1;
							else if(i1 == 1)
							{
								val1 = (*anc);								
							}
							else if(i1 == 2)
							{							
								val1 = (*part);
							}
							else if(i1 == 3) //interaction
							{
								val1 = (*anc)*(*part);	
							};

							covarAns[make_pair(i1, i2)] += val1*(*cv2)*aNumber;
						}
						else //covariate and covariate 2nd deriv
						{
							covarAns[make_pair(i1, i2)] += (*cv1)*(*cv2)*aNumber;	
						};
					};

					++i2;
					if(i2 > cp) ++cv2; //move to next covariate value
				};

				++i1;
				if(i1 > cp) ++cv1; //move to next covariate value
			};

		};//end of skipping missing data

	};//end of subject loop

	//The Hessian matrix is symetric so only calculate half and then copy
	ans[1][0] = ans[0][1];		
	ans[2][0] = ans[0][2];
	ans[2][1] = ans[1][2];

	if(fitBeta3)
	{
		ans[3][0] = ans[0][3];
		ans[3][1] = ans[1][3];			
		ans[3][2] = ans[2][3];	
	};	
	
	//setup the matrix with calculated values
	map<unsigned int, double> aCol;

	//no parameters not including covariates is given by cp
	for(unsigned int col0 = 0; col0 < cp; ++col0)
	{		
		for(unsigned int row0 = 0; row0 < cp; ++row0)
		{
			aCol[row0 + 1] = ans[col0][row0];	
		};
		hessianMatrix[col0 + 1] = aCol;
	};

	//add calculated covariate 2nd derivs
	for(unsigned int col = 0; col < noParas; ++col)
	{	
		if(col < cp) //get predefined column
		{
			aCol = hessianMatrix[col + 1];
		};

		for(unsigned int row = 0; row < noParas; ++row)
		{
			if(col >= cp && col >= row)	//top right half - except for non-covariate
			{
				aCol[row + 1] = covarAns[make_pair(row, col)];
			}
			else if(row >= cp && col < row)	//bottom left half - except for non-covariate
			{
				aCol[row + 1] = covarAns[make_pair(col, row)];
			};
		};

		hessianMatrix[col + 1] = aCol;
	};

};

//! Calculates the inverse of (X^T X) and returns the (4, 4) element needed for SE of interaction parameter.
bool TwoSNPLinearRegModel::calc_g44(double & g44, bool & useCovariates)
{
	if(useCovariates) return calc_g44Covar(g44);

	//construct matrix (X^T X), assumes all parameters are being fitted here
	list< list<double> > matrixXTX;
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	
	double sumAnc = 0, sumPart = 0, sumAncPart = 0;
	double sumAncSq = 0, sumPartSq = 0, sumAncPartSq = 0;
	double sumAncAncPart = 0, sumAncPartPart = 0;

	unsigned int totalNotMissing = 0;
	
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i)
	{
		if((*anc) != 3 && (*part) != 3 && i->second != missingQTValue) // miss out indivs where there is missing SNP data or missing QT data, for standard LR not intersted in missing partner SNP data
		{
			totalNotMissing++;			

			sumAnc += *anc;
			sumAncSq += (*anc)*(*anc);
			
			sumPart += *part;
			sumPartSq += (*part)*(*part);
			sumAncPart += (*anc)*(*part);					
	
			sumAncAncPart += (*anc)*(*anc)*(*part);
			sumAncPartPart += (*anc)*(*part)*(*part);
			sumAncPartSq += (*anc)*(*part)*(*anc)*(*part);

		};
		++anc;
		++part;
	};

	//construct matrix XTX
	//row 1
	list<double> rowOne;
	rowOne.push_back(totalNotMissing);
	rowOne.push_back(sumAnc);
	rowOne.push_back(sumPart);
	rowOne.push_back(sumAncPart);
	matrixXTX.push_back(rowOne);
	
	//row 2
	list<double> rowTwo;
	rowTwo.push_back(sumAnc);
	rowTwo.push_back(sumAncSq);
	rowTwo.push_back(sumAncPart);	
	rowTwo.push_back(sumAncAncPart);
	matrixXTX.push_back(rowTwo);
	
	//row 3
	list<double> rowThree;
	rowThree.push_back(sumPart);
	rowThree.push_back(sumAncPart);
	rowThree.push_back(sumPartSq);
	rowThree.push_back(sumAncPartPart);
	matrixXTX.push_back(rowThree);
	
	//row 4
	list<double> rowFour;
	rowFour.push_back(sumAncPart);
	rowFour.push_back(sumAncAncPart);
	rowFour.push_back(sumAncPartPart);
	rowFour.push_back(sumAncPartSq);
	matrixXTX.push_back(rowFour);
	
	//setup inverse matrix
	list< list<double> > inverse;
	list<double> row1, row2, row3, row4;
	row1.push_back(1); row1.push_back(0); row1.push_back(0); row1.push_back(0);
	row2.push_back(0); row2.push_back(1); row2.push_back(0); row2.push_back(0);
	row3.push_back(0); row3.push_back(0); row3.push_back(1); row3.push_back(0);
	row4.push_back(0); row4.push_back(0); row4.push_back(0); row4.push_back(1);

	inverse.push_back(row1);
	inverse.push_back(row2);
	inverse.push_back(row3);
	inverse.push_back(row4);

	//calculate inverse of XTX
	getInverseMatrix(matrixXTX, inverse);

	//get (4, 4) element of inverse
	list< list<double> >::const_iterator i = inverse.begin();
	i++; i++; i++;
	list<double>::const_iterator j = (*i).begin();
	j++; j++; j++;
	g44 = *j;

	return g44*0 == 0;
};

//! Calculates the inverse of (X^T X) and returns the (4, 4) element needed for SE of interaction parameter.
bool TwoSNPLinearRegModel::calc_g44Covar(double & g44)
{
	
	//construct matrix (X^T X), assumes all parameters are being fitted here
	list< list<double> > matrixXTX;
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator cc = covariateData->caseControls.begin();

	list< list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<double>::const_iterator cv, cv2;
	list<double> sumCV, covarAnc, covarPart, covarAncPart;
	list<double>::iterator sCV, cvAn, cvPa, cvAnPa;

	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		sumCV.push_back(0);
		covarAnc.push_back(0);
		covarPart.push_back(0);	
		covarAncPart.push_back(0);
	};

	list< list<double> > covarCovarSums;
	list<double> aCovarCovarSum;
	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		aCovarCovarSum.push_back(0); //set initial totals to 0			
	};

	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		covarCovarSums.push_back(aCovarCovarSum);
	};

	list< list<double> >::iterator cvcvSums;
	list<double>::iterator cvcvSums2;

	double sumAnc = 0, sumPart = 0, sumAncPart = 0;
	double sumAncSq = 0, sumPartSq = 0, sumAncPartSq = 0;
	double sumAncAncPart = 0, sumAncPartPart = 0;

	unsigned int totalNotMissing = 0;
	
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i, ++anc, ++part, ++cc, ++covar)
	{
		if((*anc) != 3 && (*part) != 3 && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data, for standard LR not intersted in missing partner SNP data
		{
			totalNotMissing++;			

			sumAnc += *anc;
			sumAncSq += (*anc)*(*anc);
			
			sumPart += *part;
			sumPartSq += (*part)*(*part);
			sumAncPart += (*anc)*(*part);					
	
			sumAncAncPart += (*anc)*(*anc)*(*part);
			sumAncPartPart += (*anc)*(*part)*(*part);
			sumAncPartSq += (*anc)*(*part)*(*anc)*(*part);

			//do covariate sums			
			sCV = sumCV.begin();
			cvAn = covarAnc.begin();
			cvPa = covarPart.begin();
			cvAnPa = covarAncPart.begin();
			cvcvSums = covarCovarSums.begin();	

			for(cv = covar->begin(); cv != covar->end(); ++cv)
			{				
				*sCV += (*cv);
				*cvAn += (*anc)*(*cv);
				*cvPa += (*part)*(*cv);
				*cvAnPa += (*anc)*(*part)*(*cv);

				//loop thro' covar*covar sums
				cvcvSums2 = cvcvSums->begin();
				for(list<double>::const_iterator cv2 = covar->begin(); cv2 != covar->end(); ++cv2, ++cvcvSums2)
				{
					*cvcvSums2 += (*cv)*(*cv2);
				};

				++sCV; ++cvAn; ++cvPa;  ++cvAnPa;				
				++cvcvSums;
			};
		};
	};

	//construct matrix XTX
	//row 1
	list<double> rowOne;
	rowOne.push_back(totalNotMissing);
	rowOne.push_back(sumAnc);
	rowOne.push_back(sumPart);
	rowOne.push_back(sumAncPart);

	//do covariate elements, row 5 onwards
	for(sCV = sumCV.begin(); sCV != sumCV.end(); ++sCV)
	{
		rowOne.push_back(*sCV);		
	};

	matrixXTX.push_back(rowOne);
	
	//row 2
	list<double> rowTwo;
	rowTwo.push_back(sumAnc);
	rowTwo.push_back(sumAncSq);
	rowTwo.push_back(sumAncPart);	
	rowTwo.push_back(sumAncAncPart);

	//do covariate elements, row 5 onwards
	for(cvAn = covarAnc.begin(); cvAn != covarAnc.end(); ++cvAn)
	{
		rowTwo.push_back(*cvAn);		
	};

	matrixXTX.push_back(rowTwo);
	
	//row 3
	list<double> rowThree;
	rowThree.push_back(sumPart);
	rowThree.push_back(sumAncPart);
	rowThree.push_back(sumPartSq);
	rowThree.push_back(sumAncPartPart);
	
	//do covariate elements, row 5 onwards
	for(cvPa = covarPart.begin(); cvPa != covarPart.end(); ++cvPa)
	{
		rowThree.push_back(*cvPa);		
	};

	matrixXTX.push_back(rowThree);
	
	//row 4
	list<double> rowFour;
	rowFour.push_back(sumAncPart);
	rowFour.push_back(sumAncAncPart);
	rowFour.push_back(sumAncPartPart);
	rowFour.push_back(sumAncPartSq);

	//do covariate elements, row 5 onwards
	for(cvAnPa = covarAncPart.begin(); cvAnPa != covarAncPart.end(); ++cvAnPa)
	{
		rowFour.push_back(*cvAnPa);		
	};

	matrixXTX.push_back(rowFour);
	
	//do covariate rows, row 5 onwards
	sCV = sumCV.begin();
	cvAn = covarAnc.begin();
	cvPa = covarPart.begin();
	cvAnPa = covarAncPart.begin();
	list<double> rowCovariate;
	

	for(cvcvSums = covarCovarSums.begin(); sCV != sumCV.end(); ++cvcvSums)
	{
		rowCovariate.push_back(*sCV); 
		rowCovariate.push_back(*cvAn); 
		rowCovariate.push_back(*cvPa);
		rowCovariate.push_back(*cvAnPa);
		
		for(cvcvSums2 = cvcvSums->begin(); cvcvSums2 != cvcvSums->end(); ++cvcvSums2)
		{
			rowCovariate.push_back(*cvcvSums2);
		};

		matrixXTX.push_back(rowCovariate);

		rowCovariate.clear();
		++sCV; ++cvAn; ++cvPa; ++cvAnPa;
	};


	//setup inverse matrix
	list< list<double> > inverse;
	unsigned int noRows = matrixXTX.size();
	list<double> aRow;

	for(unsigned int row = 1; row <= noRows; ++row)
	{
		for(unsigned int col = 1; col <= noRows; ++col)
		{
			if(row != col) aRow.push_back(0); else aRow.push_back(1);
		};

		inverse.push_back(aRow);

		aRow.clear();
	};

	//calculate inverse of XTX
	getInverseMatrix(matrixXTX, inverse);

	//get (4, 4) element of inverse
	list< list<double> >::const_iterator i = inverse.begin();
	i++; i++; i++;
	list<double>::const_iterator j = (*i).begin();
	j++; j++; j++;
	g44 = *j;

	if(g44 < 0) g44 = 0; //aviod accuracy issues giving negative answer

	return g44*0 == 0;
};

//! Fits 2 SNP model to quantitive traits.
bool TwoSNPLinearRegModel::fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & fitInter)
{
	if(useCovariates) return fitModelCovar(rss, fitAnc, fitPart, fitInter);

	//Solve equn X^T y = (X^T X)betaHat, to find betaHat, where X is the design matrix
	
	//construct vector X^T y and matrix (X^T X)
	map<unsigned int, double> vectorXTy;
	map<unsigned int, map<unsigned int, double> > matrixXTX;
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	
	double totalQT = 0, ancQT = 0, partQT = 0, ancPartQT = 0;
	double sumAnc = 0, sumPart = 0, sumAncPart = 0;
	double sumAncSq = 0, sumPartSq = 0, sumAncPartSq = 0;
	double sumAncAncPart = 0, sumAncPartPart = 0;

	unsigned int totalNotMissing = 0;
	
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i)
	{
		if((*anc) != 3 && (*part) != 3 && i->second != missingQTValue) // miss out indivs where there is missing SNP data or missing QT data, for standard LR not intersted in missing partner SNP data
		{
			totalQT += i->second;

			totalNotMissing++;			

			if(fitAnc)
			{
				ancQT += (i->second)*(*anc);
				sumAnc += *anc;
				sumAncSq += (*anc)*(*anc);
			};

			if(fitPart)
			{
				partQT += (i->second)*(*part);
				sumPart += *part;
				sumPartSq += (*part)*(*part);
				if(fitAnc) sumAncPart += (*anc)*(*part);			
			};			

			if(fitInter)
			{
				ancPartQT += (i->second)*(*anc)*(*part);
				sumAncAncPart += (*anc)*(*anc)*(*part);
				sumAncPartPart += (*anc)*(*part)*(*part);
				sumAncPartSq += (*anc)*(*part)*(*anc)*(*part);
			};

		};
		++anc;
		++part;
	};

	//construct vector XTy
	vectorXTy[1] = totalQT;
	if(fitAnc) vectorXTy[2] = ancQT;
	if(fitPart) vectorXTy[3] = partQT;	
	if(fitInter) vectorXTy[4] = ancPartQT;

	//construct matrix XTX
	//row 1
	map<unsigned int, double> rowOne;
	rowOne[1] = totalNotMissing;
	if(fitAnc) rowOne[2] = sumAnc;
	if(fitPart) rowOne[3] = sumPart;
	if(fitInter) rowOne[4] = sumAncPart;
	matrixXTX[1] = rowOne;

	if(fitAnc)
	{
		//row 2
		map<unsigned int, double> rowTwo;
		rowTwo[1] = sumAnc;
		rowTwo[2] = sumAncSq;
		if(fitPart) rowTwo[3] = sumAncPart;	
		if(fitInter) rowTwo[4] = sumAncAncPart;
		matrixXTX[2] = rowTwo;
	};

	if(fitPart)
	{
		//row 3
		map<unsigned int, double> rowThree;
		rowThree[1] = sumPart;
		if(fitAnc) rowThree[2] = sumAncPart;
		rowThree[3] = sumPartSq;
		if(fitInter) rowThree[4] = sumAncPartPart;
		matrixXTX[3] = rowThree;
	};

	if(fitInter)
	{
		//row 4
		map<unsigned int, double> rowFour;
		rowFour[1] = sumAncPart;
		rowFour[2] = sumAncAncPart;
		rowFour[3] = sumAncPartPart;
		rowFour[4] = sumAncPartSq;
		matrixXTX[4] = rowFour;
	};


	//solve matrix equation to get betaHat
	map<unsigned int, double> ans = getSolnMatrixEqun(matrixXTX, vectorXTy);

	beta0 = 0; beta1= 0; beta2 = 0;

	//set parameters
	for(map<unsigned int, double>::const_iterator par = ans.begin(); par != ans.end(); ++par)
	{
		setParameter(par->first, par->second);	
	};

	//set residual sum of squares
	anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	part = partnerSNP->noMinorAllelesAllSubjects.begin();

	double diff;
	double xiTbetaHat; //model fit

	rss = 0;

	//loop thro' subjects and calc difference of model fit with the observed data.
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i)
	{	
		if((*anc) != 3 && (*part) != 3 && i->second != missingQTValue) // miss out indivs where there is missing SNP data or missing QT data
		{
			xiTbetaHat = beta0;
			if(fitAnc) xiTbetaHat += beta1*(*anc);
			if(fitPart) xiTbetaHat += beta2*(*part);		
			if(fitInter) xiTbetaHat += beta3*(*anc)*(*part);

			diff = i->second - xiTbetaHat;
			rss += diff*diff;
		};
		++anc;
		++part;
	};

	totalNotMissingBetweenSNPs = totalNotMissing;

	return rss*0 == 0;
};

//! Fits 2 SNP model to quantitive traits.
bool TwoSNPLinearRegModel::fitModelCovar(double & rss, const bool & fitAnc, const bool & fitPart, const bool & fitInter)
{
	//Solve equn X^T y = (X^T X)betaHat, to find betaHat, where X is the design matrix
	
	//construct vector X^T y and matrix (X^T X)
	map<unsigned int, double> vectorXTy;
	map<unsigned int, map<unsigned int, double> > matrixXTX;
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator cc = covariateData->caseControls.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<double>::const_iterator cv, cv2;

	double totalQT = 0, ancQT = 0, partQT = 0, ancPartQT = 0;
	list<double> covarQT, sumCV, covarAnc, covarPart, covarAncPart;
	list<double>::iterator cvQT, sCV, cvAn, cvPa, cvAnPa;
	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		covarQT.push_back(0); //set initial totals to 0
		sumCV.push_back(0);
		covarAnc.push_back(0);
		covarPart.push_back(0);	
		covarAncPart.push_back(0);
	};

	list<list<double> > covarCovarSums;
	list<double> aCovarCovarSum;
	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		aCovarCovarSum.push_back(0); //set initial totals to 0			
	};

	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		covarCovarSums.push_back(aCovarCovarSum);
	};

	list<list<double> >::iterator cvcvSums;
	list<double>::iterator cvcvSums2;
	
	double sumAnc = 0, sumPart = 0, sumAncPart = 0;
	double sumAncSq = 0, sumPartSq = 0, sumAncPartSq = 0;
	double sumAncAncPart = 0, sumAncPartPart = 0;

	unsigned int covariateParameterNo;

	unsigned int totalNotMissing = 0;

	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i, ++anc, ++part, ++cc, ++covar)
	{
		if((*anc) != 3 && (*part) != 3 && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data
		{
			totalQT += i->second;

			totalNotMissing++;			

			if(fitAnc)
			{
				ancQT += (i->second)*(*anc);
				sumAnc += *anc;
				sumAncSq += (*anc)*(*anc);
			};

			if(fitPart)
			{
				partQT += (i->second)*(*part);
				sumPart += *part;
				sumPartSq += (*part)*(*part);
				if(fitAnc) sumAncPart += (*anc)*(*part);			
			};			
		
			if(fitInter)
			{
				ancPartQT += (i->second)*(*anc)*(*part);
				sumAncAncPart += (*anc)*(*anc)*(*part);
				sumAncPartPart += (*anc)*(*part)*(*part);
				sumAncPartSq += (*anc)*(*part)*(*anc)*(*part);
			};

			//do covariate sums		
			cvQT = covarQT.begin();
			sCV = sumCV.begin();
			cvAn = covarAnc.begin();
			cvPa = covarPart.begin();
			if(fitInter) cvAnPa = covarAncPart.begin();
			cvcvSums = covarCovarSums.begin();	
			for(cv = covar->begin(); cv != covar->end(); ++cv)
			{
				*cvQT += (i->second)*(*cv);
				*sCV += (*cv);
				*cvAn += (*anc)*(*cv);
				*cvPa += (*part)*(*cv);
				if(fitInter) *cvAnPa += (*anc)*(*part)*(*cv);

				//loop thro' covar*covar sums
				cvcvSums2 = cvcvSums->begin();
				for(list<double>::const_iterator cv2 = covar->begin(); cv2 != covar->end(); ++cv2, ++cvcvSums2)
				{
					*cvcvSums2 += (*cv)*(*cv2);
				};

				++cvQT; ++sCV; ++cvAn; ++cvPa; 
				if(fitInter) ++cvAnPa;
				++cvcvSums;
			};

		};		
	};

	//construct vector XTy
	vectorXTy[1] = totalQT;
	if(fitAnc) vectorXTy[2] = ancQT;
	if(fitPart) vectorXTy[3] = partQT;	
	if(fitInter) vectorXTy[4] = ancPartQT;

	//do covariate elements, row 5 onwards
	for(covariateParameterNo = 5, cvQT = covarQT.begin(); cvQT != covarQT.end(); ++cvQT, ++covariateParameterNo)
	{
		vectorXTy[covariateParameterNo] = *cvQT;		
	};

	//construct matrix XTX
	//row 1
	map<unsigned int, double> rowOne;
	rowOne[1] = totalNotMissing;
	if(fitAnc) rowOne[2] = sumAnc;
	if(fitPart) rowOne[3] = sumPart;
	if(fitInter) rowOne[4] = sumAncPart;

	//do covariate elements, row 5 onwards
	for(covariateParameterNo = 5, sCV = sumCV.begin(); sCV != sumCV.end(); ++sCV, ++covariateParameterNo)
	{
		rowOne[covariateParameterNo] = *sCV;		
	};

	matrixXTX[1] = rowOne;

	if(fitAnc)
	{
		//row 2
		map<unsigned int, double> rowTwo;
		rowTwo[1] = sumAnc;
		rowTwo[2] = sumAncSq;
		if(fitPart) rowTwo[3] = sumAncPart;	
		if(fitInter) rowTwo[4] = sumAncAncPart;
		//do covariate elements, row 5 onwards
		for(covariateParameterNo = 5, cvAn = covarAnc.begin(); cvAn != covarAnc.end(); ++cvAn, ++covariateParameterNo)
		{
			rowTwo[covariateParameterNo] = *cvAn;		
		};
		matrixXTX[2] = rowTwo;
	};

	if(fitPart)
	{
		//row 3
		map<unsigned int, double> rowThree;
		rowThree[1] = sumPart;
		if(fitAnc) rowThree[2] = sumAncPart;
		rowThree[3] = sumPartSq;
		if(fitInter) rowThree[4] = sumAncPartPart;
		//do covariate elements, row 5 onwards
		for(covariateParameterNo = 5, cvPa = covarPart.begin(); cvPa != covarPart.end(); ++cvPa, ++covariateParameterNo)
		{
			rowThree[covariateParameterNo] = *cvPa;		
		};
		matrixXTX[3] = rowThree;
	};

	if(fitInter)
	{
		//row 4
		map<unsigned int, double> rowFour;
		rowFour[1] = sumAncPart;
		rowFour[2] = sumAncAncPart;
		rowFour[3] = sumAncPartPart;
		rowFour[4] = sumAncPartSq;
		//do covariate elements, row 5 onwards
		for(covariateParameterNo = 5, cvAnPa = covarAncPart.begin(); cvAnPa != covarAncPart.end(); ++cvAnPa, ++covariateParameterNo)
		{
			rowFour[covariateParameterNo] = *cvAnPa;		
		};
		matrixXTX[4] = rowFour;
	};

	//do covariate rows, row 5 onwards
	sCV = sumCV.begin();
	cvAn = covarAnc.begin();
	cvPa = covarPart.begin();
	cvAnPa = covarAncPart.begin();
	map<unsigned int, double> rowCovariate;
	unsigned int covariateParameterNo2;
	cvcvSums = covarCovarSums.begin();

	for(covariateParameterNo = 5; sCV != sumCV.end(); ++covariateParameterNo, ++cvcvSums)
	{
		rowCovariate[1] = *sCV; 
		if(fitAnc) rowCovariate[2] = *cvAn; 
		if(fitPart) rowCovariate[3] = *cvPa;
		if(fitInter) rowCovariate[4] = *cvAnPa;
		
		for(covariateParameterNo2 = 5, cvcvSums2 = cvcvSums->begin(); cvcvSums2 != cvcvSums->end(); ++cvcvSums2, ++covariateParameterNo2)
		{
			rowCovariate[covariateParameterNo2] = *cvcvSums2;
		};

		matrixXTX[covariateParameterNo] = rowCovariate;
		++sCV; ++cvAn; ++cvPa; ++cvAnPa;
	};

	//solve matrix equation to get betaHat
	map<unsigned int, double> ans = getSolnMatrixEqun(matrixXTX, vectorXTy);

	beta0 = 0; beta1= 0; beta2 = 0;

	//set parameters
	for(map<unsigned int, double>::const_iterator par = ans.begin(); par != ans.end(); ++par)
	{
		setParameter(par->first, par->second);	
	};

	//set residual sum of squares
	anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	part = partnerSNP->noMinorAllelesAllSubjects.begin();
	cc = covariateData->caseControls.begin();
	covar = covariateData->covariateDataAllSubjects.begin();

	double diff;
	double xiTbetaHat; //model fit

	rss = 0;

	//loop thro' subjects and calc difference of model fit with the observed data.
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i, ++anc, ++part, ++cc, ++covar)
	{	
		if((*anc) != 3 && (*part) != 3 && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data
		{
			xiTbetaHat = beta0;
			if(fitAnc) xiTbetaHat += beta1*(*anc);
			if(fitPart) xiTbetaHat += beta2*(*part);
			if(fitInter) xiTbetaHat += beta3*(*anc)*(*part);
			covariateParameterNo = 5;
			for(cv = covar->begin(); cv != covar->end(); ++cv)
			{
				xiTbetaHat += (*cv)*(getParameter(covariateParameterNo));
				covariateParameterNo++;
			};

			diff = i->second - xiTbetaHat;
			rss += diff*diff;
		};
		
	};

	totalNotMissingBetweenSNPs = totalNotMissing;

	return rss*0 == 0;
};

//! Gradient vector with covariates.
void TwoSNPLogRegModel::getGradientVectorCovar(map<unsigned int, double> & gradientVector,  const bool & fitBeta3)  const
{
	double ans[4] = {0, 0, 0, 0};
		
	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControlsCovar.begin();

	unsigned int covariateParameterNo;
	unsigned int cp = 4; //first parameter no. for covariate in gradientVector vector
	if(fitBeta3) ++cp;
	unsigned int cpi = cp;

	for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv, ++cpi) gradientVector[cpi] = 0;	

	expEffectsAllIndiv.clear();

	//add covariates parameter derivs
	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++part, ++covar, ++cc)
	{
		if(*anc != 3 && *part != 3 && *cc != 0) //skip missing data
		{
			effects = beta0 + (*anc)*beta1 + (*part)*beta2 + (*anc)*(*part)*beta3;

			covariateParameterNo = 5;
			//add covariates effects
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			{
				effects += (*cv)*getParameter(covariateParameterNo);
				covariateParameterNo++;
			};

			//cache exp(effects) for use with the hessian matrix
			expEffects = exp(effects); 
			expEffectsAllIndiv.push_back(expEffects);

			if(*cc == 2) aNumber = -1 + expEffects/(expEffects + 1);
			else aNumber = expEffects/(expEffects + 1);		

			ans[0] += aNumber;
			ans[1] += (*anc)*aNumber;
			ans[2] += (*part)*aNumber;
			if(fitBeta3) ans[3] += (*anc)*(*part)*aNumber;

			cpi = cp;
			//add covariates parameter derivs
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv, ++cpi)
			{
				gradientVector[cpi] += (*cv)*aNumber;			
			};
		};
	};

	gradientVector[1] = ans[0];
	gradientVector[2] = ans[1];
	gradientVector[3] = ans[2];
	if(fitBeta3) gradientVector[4] = ans[3];	
};

