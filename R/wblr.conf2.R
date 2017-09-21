## wblr.conf.R
## refactored from code originally authored by Jurgen Symynck, April 2014
## Copyright 2014-2017 OpenReliability.org
#
# For more info, visit http://www.openreliability.org/
#
# For the latest version of this file, check the Subversion repository at
# http://r-forge.r-project.org/projects/Weibull-R/
##-------------------------------------------------------------------------------
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
wblr.conf <- function(x,...){
    # x is a single wblr object
        if(!(class(x)=="wblr")){
            stop('Argument \"x\" is not of class \"wblr\" ')
        }
## using findMaxDataRange from plot.wblr, which takes a list of wblr objects
## so simply convert this x to a single item list
        dr <- findMaxDataRange(list(x))

            if(!is.null(x$fit)){	
				x$fit[[length(x$fit)]]<- calculateSingleConf(
					x$fit[[length(x$fit)]],
					x$data, opadata=x$options,datarange=dr,...
				)	
            }	
            x	
        		
	

}	


