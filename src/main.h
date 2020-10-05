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


/*! \file main.h
    \brief This file defines a global variable to suppress output to screen. 
    
*/

#ifndef __MAIN
#define __MAIN

#include <cstdlib>
#include <sstream>

extern bool outputToScreen; 

extern ofstream logFile; 
extern stringstream stringLogFile; //avoid logfile going "missing", bus errors for long jobs with virtual memory pages disappearing in linux

template<typename T>

//! Outputs message to screen and log file
void out(const T & text)
{	
	if(outputToScreen) cout << text;
	stringLogFile << text;
};

template<typename T>

//! Outputs error message to screen and log file
void outErr(const T & text)
{
	cerr << text;
	stringLogFile << text;
};

#endif
