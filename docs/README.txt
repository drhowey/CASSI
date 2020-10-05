
/************************************************************************
 * JointEffects, version 1.0-Beta
 * Copyright 2012,
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of JointEffects, the SNP interaction analysis program.
 *
 * JointEffects is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * JointEffects is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JointEffects.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


REFERENCE
---------

Ueki M, Cordell HJ (2012)
Improved statistics for genome-wide interaction analysis. 
PLoS Genet. 2012; 8(4):e1002625.
PMID:22496670



INSTRUCTIONS 14/05/2012
-----------------------


The current version is designed to run under Linux only. To see a list of 
command-line options, type

./je

This will produce a listing of the possible options (which is output both to 
the screen and to the file je.log)

To run on some data, type:

./je -i file.bed

where file.bed is a PLINK format .bed file (requires PLINK file.bim and 
file.fam files to exist in same directory)

If both SNPS/windows of SNPs to be analysed exist within this file, then 
no additional input files are required. Alternatively, you can specify 
that the first SNP/window of SNPs should come from this file and the 
second SNP/window of SNPs should come from a different (second) file, using:

./je -i file.bed -i2 file2.bed




If the executable file (je) does not run on your system, you may need 
to recompile using:

g++ -O3 *.cpp -o je2

and then run using 

./je2

