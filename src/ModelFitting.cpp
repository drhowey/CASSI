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

/*! \file ModelFitting.cpp
    \brief This file contains the source of the model fitting methods.
    
*/

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "ModelFitting.h"
#include "Model.h"
#include "main.h"

//! Copies values from one map to another.
void setParameterValues(map<unsigned int, double> & parametersToSet, map<unsigned int, double> & values)
{
	for(map<unsigned int, double>::const_iterator j = values.begin();  j != values.end(); ++j)
	{			
		parametersToSet[j->first] = j->second;
	};
};

//! Solves matrix equation.
map<unsigned int, double> getSolnMatrixEqun(map<unsigned int, map<unsigned int, double> > & matrix, map<unsigned int, double> & vect)
{
	map<unsigned int, double> ans = vect;

	//loop thro' rows of matrix, m, and the inverse, i
	map<unsigned int, map<unsigned int, double> >::iterator mrow = matrix.begin();	
	map<unsigned int, map<unsigned int, double> >::iterator mrow2;	
	map<unsigned int, double>::iterator mcol;	
	map<unsigned int, double>::iterator mcol2;	
	map<unsigned int, double>::iterator ansIt = ans.begin();
	map<unsigned int, double>::iterator ansIt2;
	
	double factor;
	unsigned int rowNo = 1;
	unsigned int colNo = 1;

	for( ; mrow != matrix.end(); ++mrow, ++ansIt)
	{
		//set column to the first column in the row
		mcol = mrow->second.begin();		
		colNo = 1;

		//advance the column until the the row no. is equal to the column no.
		while(colNo != rowNo)
		{
			mcol++;			
			++colNo;
		};

		//divide the row in (m and i) by the value in the matrix, m, at (rowNo, colNo)
		factor = 1.0/(mcol->second); 
		mcol->second = 1;		
		mcol++;		

		//scale the remaining elements in the row - if there are any
		while(mcol != mrow->second.end())
		{
			mcol->second *= factor;			
			mcol++;			
		};

		//scale answer vector by the same factor
		ansIt->second *= factor;

		//subtract the row in question away from the remaining rows scaled  s.t. column = mrow will be zero below this row in matrix m
		mrow2 = mrow;
		ansIt2 = ansIt;		
		mrow2++;
		ansIt2++;
		
		//loop thro' remaining rows
		while(mrow2 != matrix.end())
		{
			//set column iterators to begining of the rows
			mcol2 = mrow2->second.begin();			
			mcol = mrow->second.begin();			

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = 1;
			while(colNo != rowNo)
			{
				mcol++;
				mcol2++;				
				++colNo;
			};

			factor = mcol2->second; //factor to multiple row, rowNo by to take away from subseq. rows
			mcol2->second -= (mcol->second)*factor;//0;
			mcol++;
			mcol2++;

			//subtract scaled row for the rest of the matrix, m
			while(mcol2 != mrow2->second.end())
			{
				mcol2->second -= (mcol->second)*factor;				
				mcol++;
				mcol2++;				
			};

			//subtract scaled value from ans vector
			ansIt2->second -= (ansIt->second)*factor;			

			mrow2++;
			ansIt2++;			
		};//end of performing row operations to set column (=rowNo) to zeroes below row = rowNo

		++rowNo;	
	};//end of performing row operations to set lower left of matrix, m, to zeroes

	//Now reduce the upper right of matrix, m, to zero
	map<unsigned int, map<unsigned int, double> >::reverse_iterator mrowre = matrix.rbegin();	
	map<unsigned int, map<unsigned int, double> >::reverse_iterator mrowre2 = matrix.rbegin();	
	map<unsigned int, double>::reverse_iterator ansItre = ans.rbegin();
	map<unsigned int, double>::reverse_iterator ansItre2;
	map<unsigned int, double>::reverse_iterator mcolre2;
	map<unsigned int, double>::reverse_iterator mcolre;

	rowNo = matrix.size();

	for( ; mrowre != matrix.rend(); ++mrowre, ++ansItre)
	{

		mrowre2 = mrowre;		
		ansItre2 = ansItre;
		mrowre2++;
		ansItre2++;		

		//loop tho' the remaining rows backwards - if there are any
		while(mrowre2 != matrix.rend())
		{			
			//set column iterators to begining of the rows
			mcolre2 = mrowre2->second.rbegin();			
			mcolre = mrowre->second.rbegin();		

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = mrowre2->second.size();//size will be 4
			while(colNo != rowNo)
			{
				mcolre++;
				mcolre2++;				
				--colNo;
			};

			factor = mcolre2->second; //factor to multiple row, rowNo by to take away from subseq. rows
			mcolre2->second -= (mcolre->second)*factor;//0;
			mcolre++;
			mcolre2++;

			//subtract scaled row from the rest of the matrix, m
			while(mcolre2 != mrowre2->second.rend())//could stop at when col < row
			{
				mcolre2->second -= (mcolre->second)*factor;				
				mcolre++;
				mcolre2++;				
			};

			//subtract scaled value from ans vector
			ansItre2->second -= (ansItre->second)*factor;			

			mrowre2++;
			ansItre2++;			
		};

		--rowNo;
	};

	return ans;
};

//! Function used in Newton's Method.
//double getDistanceVectors(map<unsigned int, double> & p1, map<unsigned int, double> & p2)
//{
//	double ans = 0;
//	map<unsigned int, double>::const_iterator v1 = p1.begin();
//	map<unsigned int, double>::const_iterator v2 = p2.begin();
//
//	while(v1 != p1.end())	
//	{
//		ans += (v1->second - v2->second)*(v1->second - v2->second);
//		++v1; ++v2;
//	};
//
//	return ans;
//};

//! Finds the minimum negative log likelihood.
bool FindFit::newtonsMethod(double & eval, set<unsigned int> & parasToFit, double accuracy, unsigned int maxIterations) const
{
	double prevEval = evaluate();

	if(prevEval == 0)
	{
			return false;
	};

	double distance = 0;
	double diff;
	bool fitBeta3 = false;
	
	set<unsigned int>::const_iterator fb = parasToFit.find(4);
	if(fb != parasToFit.end()) fitBeta3 = true;

	map<unsigned int, double> prevPoint;	
	map<unsigned int, double> nextPoint;

	unsigned int noIterations = 1;

	//setup prevPoint for the initial point	
	map<unsigned int, double> modelParameters = modelToFit->getParameters(); 
 
	for(map<unsigned int, double>::const_iterator mp = modelParameters.begin(); mp != modelParameters.end(); ++mp)
	{
		if(parasToFit.find(mp->first) != parasToFit.end())
		{
			prevPoint[mp->first] = mp->second;
			nextPoint[mp->first] = mp->second;
		};
	};

	do{

		getNextPoint(nextPoint, fitBeta3);

		modelToFit->setNewParameters(nextPoint);

		//check for convergence compare last two points and/or loglikelihoods
		//distance = getDistanceVectors(prevPoint, nextPoint);
		//if(distance*0 != 0) return false;
		//else if(distance < accuracy) break;

		//use neg log likelihood to monitor convergence
		eval = evaluate();
		diff = prevEval - eval;
		if(diff < 0) diff = -diff;

		if(diff*0 != 0) return false;
		else if(diff < accuracy) break;

		setParameterValues(prevPoint, nextPoint);
		prevEval = eval;
		noIterations++;
	}while(noIterations <= maxIterations);

	//eval = evaluate();

	return true;
};

//! Calculates the next point for one iteration of Newton's method for optimization.
void FindFit::getNextPoint(map<unsigned int, double> & point, const bool & fitBeta3) const
{
	map<unsigned int, double> pointChange = getPointChange(fitBeta3);
	map<unsigned int, double>::const_iterator pc = pointChange.begin();
	map<unsigned int, double>::const_iterator pp = point.begin();

	for( ; pp != point.end(); ++pp, ++pc)
	{
		point[pp->first] -= pc->second;
	};
};

//! Returns vector on how to change the best fit parameters when using Newton's Method.
map<unsigned int, double> FindFit::getPointChange(const bool & fitBeta3) const
{	
	map<unsigned int, double> ans;

	map<unsigned int, map<unsigned int, double> > hessianMatrix; // col, row
	map<unsigned int, double> gradientVector;
	map<unsigned int, double> aVector;
	
	//get the gradient vector, 1st derivatives with repect to each parameters
	modelToFit->getGradientVector(gradientVector, fitBeta3);

	//get the Hessian matrix, 2nd derivatives with repect to each pair of parameters
	modelToFit->getHessianMatrix(hessianMatrix, fitBeta3);

	//solve Hess*ans = gradVect
	ans = getSolnMatrixEqun(hessianMatrix, gradientVector);

	return ans;
};

