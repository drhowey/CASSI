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


/*! \file ModelFitting.h
    \brief This file contains model fitting methods.
    
*/

#ifndef __MODELFITTING
#define __MODELFITTING

#include <map>

using namespace std;

#include "Model.h"

map<unsigned int, double> getSolnMatrixEqun(map<unsigned int, map<unsigned int, double> > & matrix, map<unsigned int, double> & vect);
	
//class to find a good fit of the parameters to the data 
class FindFit
{
private:

	Model * modelToFit; //model
	
public:

	FindFit(Model * m) : modelToFit(m) { };

	~FindFit()
	{
		
	};

	bool newtonsMethod(double & eval, set<unsigned int> & parasToFit, double accuracy = 1e-6, unsigned int maxIterations = 50) const;

	double evaluate() const {return modelToFit->negLogLikelihood();};
	
private:

	//Newton's method for optimization methods	
	void getNextPoint(map<unsigned int, double> & point, const bool & fitBeta3) const;
	map<unsigned int, double> getPointChange(const bool & fitBeta3) const;
};


#endif
