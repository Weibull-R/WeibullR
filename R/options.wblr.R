## options.wblr.R
## Refactored from options.abrem.R originally authored by Jurgen Symynck
## (c) 2014-2017 OpenReliability.org
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


options.wblr<- function(...){
    # function to handle the many options of the Weibull-R functions
    # the option list should only be manipulated through this function!
    
    # TODO: WARNING: partial matching is in effect!
    # options.wblr()$ylim will return options.wblr()$ylim.default if
    # $ylim was set to NULL!

    single <- FALSE
    args <- list(...)

    if(!exists(as.character(substitute(options_wblr))))
        # if the globally accessible variable was not defined yet, then
        # create it here with default values OR reset to default values
        # message ("Resetting Weibull-R options to default values...")
        options_wblr <<- list(
            dist="weibull",
            method.fit=c("rr","xony"),
##            conf.what="blives",
##            conf.blives.sides="double",
##            unrel.n=25,
            method.conf="mcpivotals",
            pp="median",#,"benard","hazen","mean", "kaplan.meier", "blom"),
            S=1e4,
            pivotals=FALSE,
##            cl=0.9,
			ci=0.9,
##            unrel=c(0.1,0.05,0.01),
			dq="abrem",
            qblives=c(0.1,0.05,0.01),
            mar=c(5.1,4.1,5.1,2.1),
            main="Probability Plot",
            main.contour="Contour Plot",
                # a default title for Contour plots
            sub=NULL,
            sub.contour=NULL,
            xlim=NULL,
            ylim=NULL,
            xlab="Time To Failure",
            ylab="Unreliability [%]",
            log="x", # maybe this should be removed in favor of "canvas"
            #canvas="weibull",
            coordinate.text.size=0.7,
            signif=4,
            pch=1,
            lwd=2,
            lwd.points=2,
            cex.points=1,
            lty=1,
            col="black",
            col.grid="gray",
			interval.col="black",
			interval.lty="dashed",
			interval.lwd=1,

            is.plot.grid=TRUE,
            is.plot.fit=TRUE,
            is.plot.pp=TRUE,
            is.plot.ppcoordinates=FALSE,
            is.plot.legend=TRUE,

            #         legend.position="bottomright",
            legend.text.size=0.7,
            label="",
            in.legend=TRUE,
            in.legend.blives=TRUE,
            in.legend.gof=TRUE,
            is.plot.cb = TRUE,
            persistent=TRUE)
            
    if (!length(args))
        args <- options_wblr
           # return the current option list
    else {
        if (all(unlist(lapply(args, is.character))))
            # if all items in the args are characters, then
            # treat them as the names of the options.
            args <- as.list(unlist(args))
        if (length(args) == 1) {
            if (is.list(args[[1L]]) | is.null(args[[1L]]))
                args <- args[[1L]]
                # unlist the first (and only) argument to a string
            else if(is.null(names(args)))
                # if there is no name to args then
                # the arg itself is the name (?)
            single <- TRUE
        }
    }
    value <- args
    if(options_wblr$persistent){
        options_wblr <<-modifyList(options_wblr, value)
    }
    if(!is.null(args$persistent)){
        value <- args
        if(args$persistent){
            options_wblr <<-modifyList(options_wblr, value)
        }
    }
    # make the options stick between calls of options.wblr()
    if(is.null(names(args)))
        value <- options_wblr[match(args,names(options_wblr))]
    if(single) value <- value[[1L]]
    value
}
# TODO :options that are NULL are not shown in the printout


splitargs <- function(...){
    arg         <- list(...)
    argnames    <- names(arg)
    parplot     <- plot_default_args()
    ret         <- list()
    opanames    <- names(options.wblr())
    ret$opa     <- arg[tolower(argnames) %in% tolower(opanames)]
        # ret$opa can be an emply list, which is compatible with modifyList()
    ret$rem     <- arg[!(tolower(argnames) %in% tolower(opanames))]
    ret
}

plot_default_args <- function(){
    paronly <- c("ask","fig", "fin","lheight","mai", "mar", "mex", "mfcol",
        "mfrow", "mfg","new","oma", "omd", "omi","pin", "plt", "ps", "pty",
        "usr","xlog", "ylog","ylbias")
        # parameters that can only be set using par()
        # see $par() for the origin of this list
    parreadonly <- c("xlog", "ylog", "adj", "ann", "ask", 
        "bg", "bty", "cex", "cex.axis", "cex.lab", 
        "cex.main", "cex.sub", "col", "col.axis", "col.lab", 
        "col.main", "col.sub", "crt", "err", "family", 
        "fg", "fig", "fin", "font", "font.axis",
        "font.lab", "font.main", "font.sub", "lab", "las", 
        "lend", "lheight", "ljoin", "lmitre", "lty", 
        "lwd", "mai", "mar", "mex", "mfcol", 
        "mfg", "mfrow", "mgp", "mkh", "new", 
        "oma", "omd", "omi", "pch", "pin", 
        "plt", "ps", "pty", "smo", "srt", 
        "tck", "tcl", "usr", "xaxp", "xaxs", 
        "xaxt", "xpd", "yaxp", "yaxs", "yaxt", 
        "ylbias")
        # par() parameter that can be set
        # par(no.readonly=TRUE)
    parplot <- unique(sort(c(parreadonly[!(parreadonly %in% paronly)],
        "type","xlim","ylim","log","main","sub","xlab","ylab",
          "ann","axes","frame.plot","panel.first","panel.last","asp")))
          # all valid (?) graphical parameters that can be supplied
          # to plot.default
    parplot
}

