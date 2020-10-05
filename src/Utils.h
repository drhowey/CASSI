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


/*! \file utils.h
    \brief This file defines some useful functions for calculating statistics. 
    
*/

#ifndef __UTILS
#define __UTILS

//! Inverts a matrix.
void getInverseMatrix(list< list<double> > & matrix, list< list<double> > & inverse);

//! Constant used in chi square calculations.
const double oneOverSqRoot2 = 0.7071067811865475244008443621048490392848359376884740365883;

//! Calculates the p-value for a given chi square value with df=1.
double getPvalueChiSq1DF(const double & chisq);

//! Calculates the p-value from a f-statistic with d1 and d2 dfs.
double getPvalueFStat(double & fstat, const unsigned int & d1, const unsigned int & d2);

//! Calculates the p-value for a given standard normal z.
double getPvalueZSqd(double & zsqd);

//! Calculate the Z value corresponding to a p-value.
double calculateZSqdFromPvalue(const double & pval);

//! Calculate the Chi sq value (with 1 df) corresponding to a p-value.
double calculateChiSqFromPvalue(const double & pval);

//! Estimates haplotype frequencies from joint genotype counts.
void estimateHaplotypeFreqsEM(JointGenotypeCounts * jointGenotypeCounts, double & P11, double & P12, double & P21, double & P22);

#endif
