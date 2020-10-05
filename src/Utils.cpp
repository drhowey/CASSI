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


/*! \file Utils.cpp
    \brief This file contains the source some useful functions for calculating statistics. .
    
*/

#include <map>
#include <list>
#include <math.h>
#include <string>

using namespace std; // initiates the "std" or "standard" namespace

#include "Data.h"
#include "Utils.h"
#include "cdflib.h"

//! Inverts a matrix. the inverse should be set to the identity when given to this function
void getInverseMatrix(list< list<double> > & matrix, list< list<double> > & inverse)
{

	//loop thro' rows of matrix, m, and the inverse, i
	list< list<double> >::iterator mrow = matrix.begin();
	list< list<double> >::iterator irow = inverse.begin();
	list< list<double> >::iterator mrow2;
	list< list<double> >::iterator irow2;
	list<double>::iterator mcol;
	list<double>::iterator icol;
	list<double>::iterator mcol2;
	list<double>::iterator icol2;
	double factor;
	unsigned int rowNo = 1;
	unsigned int colNo = 1;

	for( ; mrow != matrix.end(); ++mrow, ++irow)
	{
		//set column to the first column in the row
		mcol = mrow->begin();
		icol = irow->begin();
		colNo = 1;

		//advance the column until the the row no. is equal to the column no.
		while(colNo != rowNo)
		{
			mcol++;			
			++colNo;
		};

		//divide the row in (m and i) by the value in the matrix, m, at (rowNo, colNo)
		factor = 1.0/(*mcol); //divide all elements instead?
		*mcol = 1;		
		mcol++;		

		//scale the remaining elements in the row - if there are any
		while(mcol != mrow->end())
		{
			*mcol *= factor;			
			mcol++;			
		};

		//scale all of the elements in the inverse for this row
		while(icol != irow->end())
		{			
			*icol *= factor;			
			icol++;
		};

		//subtract the row in question away from the remaining rows scaled  s.t. column = mrow will be zero below this row in matrix m
		mrow2 = mrow;
		irow2 = irow;
		mrow2++;
		irow2++;
		//loop thro' remaining rows
		while(mrow2 != matrix.end())
		{
			//set column iterators to begining of the rows
			mcol2 = mrow2->begin();
			icol2 = irow2->begin();
			mcol = mrow->begin();
			icol = irow->begin();

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = 1;
			while(colNo != rowNo)
			{
				mcol++;
				mcol2++;				
				++colNo;
			};

			factor = *mcol2; //factor to multiple row, rowNo by to take away from subseq. rows
			*mcol2 -= (*mcol)*factor;//0;
			mcol++;
			mcol2++;

			//subtract scaled row for the rest of the matrix, m
			while(mcol2 != mrow2->end())
			{
				*mcol2 -= (*mcol)*factor;				
				mcol++;
				mcol2++;				
			};

			//now perform the same row operation on the inverse matrix, i
			while(icol2 != irow2->end())
			{
				*icol2 -= (*icol)*factor;
				icol++;
				icol2++;				
			};

			mrow2++;
			irow2++;
		};//end of performing row operations to set column (=rowNo) to zeroes below row = rowNo

		++rowNo;
	};//end of performing row operations to set lower left of matrix, m, to zeroes

	//Now reduce the upper right of matrix, m, to zero
	list< list<double> >::reverse_iterator mrowre = matrix.rbegin();
	list< list<double> >::reverse_iterator irowre = inverse.rbegin();
	list< list<double> >::reverse_iterator mrowre2 = matrix.rbegin();
	list< list<double> >::reverse_iterator irowre2 = inverse.rbegin();
	list<double>::reverse_iterator mcolre2;
	list<double>::reverse_iterator mcolre;

	rowNo = matrix.size();

	for( ; mrowre != matrix.rend(); ++mrowre, ++irowre)
	{

		mrowre2 = mrowre;
		irowre2 = irowre;
		mrowre2++;
		irowre2++;

		//loop tho' the remaining rows backwards - if there are any
		while(mrowre2 != matrix.rend())
		{			
			//set column iterators to begining of the rows
			mcolre2 = mrowre2->rbegin();
			icol2 = irowre2->begin();
			mcolre = mrowre->rbegin();
			icol = irowre->begin();

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = mrowre2->size();//size will be 4
			while(colNo != rowNo)
			{
				mcolre++;
				mcolre2++;				
				--colNo;
			};

			factor = *mcolre2; //factor to multiple row, rowNo by to take away from subseq. rows
			*mcolre2 -= (*mcolre)*factor;//0;
			mcolre++;
			mcolre2++;

			//subtract scaled row from the rest of the matrix, m
			while(mcolre2 != mrowre2->rend())//could stop at when col < row
			{
				*mcolre2 -= (*mcolre)*factor;				
				mcolre++;
				mcolre2++;				
			};

			//now perform the same row operation on the inverse matrix, i
			while(icol2 != irowre2->end())
			{
				*icol2 -= (*icol)*factor;
				icol++;
				icol2++;				
			};			

			mrowre2++;
			irowre2++;
		};

		--rowNo;
	};

};

//! Calculates the p-value from a Chi square value with 1 df.
double getPvalueChiSq1DF(const double & chisq)
{	
	if(chisq*0 != 0) return 1.0;
	double a = sqrt(chisq)*oneOverSqRoot2;
	int ind = 0;
	return erfc1(&ind, &a);
};

//! Calculates the p-value from a f-statistic with d1 and d2 dfs.
double getPvalueFStat(double & fstat, const unsigned int & d1, const unsigned int & d2)
{	
	if(fstat <= 0) return 1.0;
	double x = (double)(d1*fstat)/(double)(d1*fstat + d2);
	double y = 1 - x;
	double a = (double)(d1)/2.0;
	double b = (double)(d2)/2.0;
	double w, w1;
	int ierr;

	bratio(&a, &b, &x, &y, &w, &w1, &ierr);

	return w1;
};

//! Calculates the p-value for a given standard normal z.
double getPvalueZSqd(double & zsqd)
{
	double integral;
	double compIntegral; // = 1 - integral	
	double z = sqrt(zsqd);

	cumnor(&z, &integral, &compIntegral); //computes integral of -inf to x of standard normal distribution

	return 2*compIntegral; //2 sided z-score p-value
};

//! Calculate the Z value corresponding to a p-value.
double calculateZSqdFromPvalue(const double & pval)
{
	double p = pval*0.5;
	double q = 1 - p;
	double z = dinvnr(&p, &q);
	return z*z;
};

// Calculate the Chi sq value (with 1 df) corresponding to a p-value using normal distribution.
double calculateChiSqFromPvalue(const double & pval)
{
	return calculateZSqdFromPvalue(pval);
};

//! Estimates the haplotype frequencies between two SNPs using the Expectation-Maximization Algorithm.
void estimateHaplotypeFreqsEM(JointGenotypeCounts * jointGenotypeCounts, double & P11, double & P12, double & P21, double & P22)
{	
	//total count of haplotypes
	double total = 2*(double)(jointGenotypeCounts->total);
	
	//total count of known haplotypes
	double knownHapCount11 = 2*jointGenotypeCounts->counts[0][0]
							+ jointGenotypeCounts->counts[0][1]
							+ jointGenotypeCounts->counts[1][0];

	double knownHapCount12 = jointGenotypeCounts->counts[0][1]
							+ 2*jointGenotypeCounts->counts[0][2]							
							+ jointGenotypeCounts->counts[1][2];

	double knownHapCount21 = jointGenotypeCounts->counts[1][0]						
							+ 2*jointGenotypeCounts->counts[2][0]
							+ jointGenotypeCounts->counts[2][1];

	double knownHapCount22 = jointGenotypeCounts->counts[1][2]
							+ jointGenotypeCounts->counts[2][1]
							+ 2*jointGenotypeCounts->counts[2][2];

	double unknownCount = jointGenotypeCounts->counts[1][1]*0.5;

	double prevEstHapFreq11, prevEstHapFreq12, prevEstHapFreq21, prevEstHapFreq22;
	double genoTypeCount11 =  jointGenotypeCounts->counts[1][1];
	double estFreqHaps11and22;
	double estFreqHaps12and21;
	double diff;
	double accuracy = 1e-8;
	unsigned int maxIts = 1000;
	unsigned int iterations = 0;

	//"2" is the minor allele, "1" is the major allele
	P11 =  (knownHapCount11 + unknownCount)/total;
	P12 =  (knownHapCount12 + unknownCount)/total;
	P21 =  (knownHapCount21 + unknownCount)/total;
	P22 =  (knownHapCount22 + unknownCount)/total;

	do{

		prevEstHapFreq11 = P11;
		prevEstHapFreq12 = P12;
		prevEstHapFreq21 = P21;
		prevEstHapFreq22 = P22;

		//expectation step
		estFreqHaps11and22 = (P11*P22)/(P11*P22 + P12*P21);
		estFreqHaps12and21 = 1 - estFreqHaps11and22;

		//maximization step
		P11 = (genoTypeCount11*estFreqHaps11and22 + knownHapCount11)/total;
		P22 = (genoTypeCount11*estFreqHaps11and22 + knownHapCount22)/total;

		P12 = (genoTypeCount11*estFreqHaps12and21 + knownHapCount12)/total;
		P21 = (genoTypeCount11*estFreqHaps12and21 + knownHapCount21)/total;

		diff = fabs(P11 - prevEstHapFreq11) + fabs(P12 - prevEstHapFreq12) + fabs(P21 - prevEstHapFreq21) + fabs(P22 - prevEstHapFreq22);

		iterations++;

	}while(diff > accuracy && iterations < maxIts);

};

