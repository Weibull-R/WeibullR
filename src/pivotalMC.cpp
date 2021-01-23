/* pivotalMC.cpp					
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * These functions are distributed in the hope that they will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 * Author: David J. Silkworth
 *     Copyright (C) 2013-2021 OpenReliability.org
 */
 
#include "WeibullR.h"	
#include "Pivotal.h"
			
    using namespace Rcpp ;
	
SEXP pivotalMC(SEXP arg1) {
// Construct an LSLRmodel with this input argument
	std::unique_ptr<Pivotal> PIV(new Pivotal(arg1));
	return PIV->Execute();
}