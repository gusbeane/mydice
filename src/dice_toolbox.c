/*-----------------------------------------------------------------------------
/
/ Filename: dice_init.c
/ Author: Valentin Perret
/ Author's email: perret.valentin@gmail.com
/ Description: DICE creates galaxies. 
/
/	       DICE uses the GNU Scientific Library (GSL). You can
/	       download the GSL source code from:
/
/		http://www.gnu.org/software/gsl
/
/	       or replace it with another math library.
/
/ Copyright Information:
/
/ Copyright (c) 2014       Valentin Perret
/
/ This program is free software; you can redistribute it and/or modify
/ it under the terms of the GNU General Public License as published by
/ the Free Software Foundation; either version 2 of the License, or
/ (at your option) any later version.
/
/ This program is distributed in the hope that it will be useful,
/ but WITHOUT ANY WARRANTY; without even the implied warranty of
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/ GNU General Public License for more details.
/
/ You should have received a copy of the GNU General Public License
/ along with this program; if not, write to the Free Software
/ Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
/
/ The license is also available at:
/
/		http://www.gnu.org/copyleft/gpl.html .
/
/ Date: October 2014
/
*///---------------------------------------------------------------------------

#include "dice.h"

// Find the minimum between two values
double min(double a, double b) {
	if(a<b) return a;
	else    return b;
}

// Find the maximum between two values
double max(double a, double b) {
	if(a>b) return a;
	else    return b;
}

// Derivate a function using a 4-point central scheme
double deriv_central(galaxy *gal, double x, double h, function_to_derivate F){
	double new_x,f1,f2,f3,f4,derivative;

    new_x = x+2.0*h;
    f1 = F(new_x,gal);
    new_x = x+h;
    f2 = F(new_x,gal);
    new_x = x-h;
    f3 = F(new_x,gal);
    new_x = x-2.0*h;
    f4 = F(new_x,gal);
    derivative = (-f1 + 8.0*f2 - 8.0*f3 + f4)/(12.0*h*kpc);

	return derivative;
}

// Derivate a function using a 2-point forward scheme
double deriv_forward(galaxy *gal, double x, double h, function_to_derivate F){
	double new_x,f1,f2,f3,f4,derivative;

    f1 = F(x,gal);
    new_x = x+h;
    f2 = F(new_x,gal);
    
    derivative = (f2-f1)/(h*kpc);
	return derivative;
}
