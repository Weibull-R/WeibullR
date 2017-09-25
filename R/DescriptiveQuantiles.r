# DescriptiveQuantiles.r
 ##
 ## Author: David J. Silkworth
 ##   (c)2017 OpenReliability.org
##
#-------------------------------------------------------------------------------
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
## Descriptive quantiles are the percentile positions at which points on a [potentially] curved line are 
## calculated, so that  a smoothed curve can be plotted by either spline or linear interpolation.
## It can be desired to set these to specific values so that tabular output of calculated points
## can be reported for comparison with other software. This function will return a vector of quantiles
## based on a named set that can be stored on the options.wblr list.
## It is expected that this function will most often be called by its aleas, DQ.

DescriptiveQuantiles<-function(dqlabel)  {

	if(tolower(dqlabel)=="minitab") {
	## these descriptive quantiles match Minitab unchangeable defaults (27 values)
		dq<-c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))
	}
	
	if(tolower(dqlabel)=="supersmith") {
		## descriptive quantiles for comparison with SuperSMITH (limit of 15 values)	
		dq<-c(.01, .02, .05, .10, .15, .20, .30, .40, .50,  .60, .70, .80, .90, .95, .99)
	}
	
	if(tolower(dqlabel)=="abrem")  {
	## this is the original default by Jurgen Symynck for package abrem
	## it produces 25 evenly spaced points across the y limits of a weibull canvas, including ends
	#F0(seq(F0inv(1e-3), F0inv(0.999),length.out=25))
	dq<-1-exp(-exp(seq(log(qweibull((1e-3),1,1)), log(qweibull((0.999),1,1)),length.out=25)))	
	}
	
	dq
}

DQ<-DescriptiveQuantiles