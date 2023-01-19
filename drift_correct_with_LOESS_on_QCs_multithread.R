drift_correct_with_LOESS_on_QCs_vec <- function(data.input, is.QC, running.order, loess.span.vals = NULL, loess.cv.folds = NULL){
        require( package="bisoreg" ) # Required for "loess.wrapper".
        require( package="caret" ) # Required for "R2".
        
        
        # extract training data
        data.training <- data.input[is.QC]
        running.order.training <- running.order[is.QC]
        
        # Use default parameters if none supplied
        if(is.null( loess.span.vals )) loess.span.vals <- seq(0.25, 1, by = 0.05)
        if(is.null( loess.cv.folds ))  loess.cv.folds  <- 5
        
        
        # Cross-validate the LOESS model for the feature
        model <- loess.wrapper( x=running.order.training, y=data.training, span.vals=loess.span.vals, folds=loess.cv.folds ) 
        
        
        # Predict all observations for the feature using the optimal model.
        predictions <- predict( object=model, newdata=running.order )
        
        
        # Correct the observations that have a non-missing prediction.
        # (Extrapolations are not possible.)
        not.na <- !is.na( predictions )
        
        
        # make the output list
        output <- list()
        output$x.corrected <- output$x.predicted <- output$x.original <- data.input
        
        
        # Denoise the signal by subtracting the predicted technical signal from the observations for the feature
        output$x.corrected[ not.na ] <- data.input[ not.na ] - predictions[ not.na ] + mean( data.input, na.rm=TRUE )
        
        
        # Store the predictions.
        output$x.predicted <- predictions
        
        
        # Compute the coefficient of determination for all observations.
        output$R.sq <- rep( x=NA, times=length( output$x.corrected ) )
        output$R.sq <- R2( pred=predictions[ not.na ], obs=data.input[ not.na ] )
        
        
        return(output)
}










drift_correct_with_LOESS_on_QCs <- function( x, running.order, is.QC, loess.cv.folds=NULL, loess.span.vals=NULL, ncores = NULL ) {
        
        #
        # Arguments:
        #
        # x: Matrix of observations by metabolites.
        # running.order: Numeric vector of running order.
        # is.QC: Boolean vector informing about the QC samples.
        # loess.span.vals: (Optional) values of the tuning parameter for the LOESS regression.
        # loess.cv.folds: (Optional) number of cross-validation folds.
        #
        # Author:
        #
        # Tommi Suvitaival
        # tsuv0001@regionh.dk
        #
        # 9.5.2017
        #
        
        # Check that the number of observations is the same in all the inputs.
        
        if ( length( running.order ) != nrow( x ) | length( is.QC ) != nrow( x ) ) {
                
                stop()
                
        }
        
        require( package="rlist" ) # for getting the results into matrixes easily
        require( package="parallel" ) # for parallel apply
        
        
        x.list = as.list(as.data.frame(x))
        
        if(is.null( ncores )) ncores = min(detectCores(), ncol(x))
        cl <- makeCluster(ncores)
        
        
        output.l <- parLapply(cl, # cluster
                              x.list, #input
                              drift_correct_with_LOESS_on_QCs_vec, #fun
                              is.QC = is.QC, 
                              running.order = running.order, 
                              loess.span.vals = NULL, 
                              loess.cv.folds = NULL
        )
        
        
        stopCluster(cl)
                                            
        output <- list()
        output$x.original  <-  unname(list.cbind(unlist(list.select(output.l, x.original ),  recursive=FALSE)))
        output$x.predicted <-  unname(list.cbind(unlist(list.select(output.l, x.predicted),  recursive=FALSE)))
        output$x.corrected  <- unname(list.cbind(unlist(list.select(output.l, x.corrected ), recursive=FALSE)))
        output$R.sq <- unname(unlist(list.select(output.l, R.sq)))
        
        output$running.order <- running.order
        output$is.QC <- is.QC
        
        # put the original names
        names( output$R.sq ) <- colnames( x )
        names( output$running.order ) <- rownames( running.order )
        names( output$is.QC ) <- rownames( is.QC )
        
        
        colnames( output$x.original ) <- colnames( x )
        colnames( output$x.predicted ) <- colnames( x )
        colnames( output$x.corrected ) <- colnames( x )
        rownames( output$x.original ) <- rownames( x )
        rownames( output$x.predicted ) <- rownames( x )
        rownames( output$x.corrected ) <- rownames( x )
        
        return( output )
        
}



plot_LOESS_correction <- function( x, index=NULL ) {
        
        #
        # Arguments:
        # x: value from the "drift_correct_with_LOESS_on_QCs" function.
        #
        # Author:
        #
        # Tommi Suvitaival
        # tsuv0001@regionh.dk
        #
        # 9.5.2017
        #
        
        if ( !is.null( index ) ) {
                
                x$x.original <- x$x.original[ , index, drop=FALSE ]
                x$x.predicted <- x$x.predicted[ , index, drop=FALSE ]
                x$x.corrected <- x$x.corrected[ , index, drop=FALSE ]
                x$R.sq <- x$R.sq[ index ]
                
        }
        
        N.features <- ncol( x$x.corrected )
        
        N.features.per.plot <- 10
        N.plots <- ceiling( N.features / N.features.per.plot )
        
        cols <- 1:N.features.per.plot
        pchs <- paste( ( 1:N.features.per.plot )-1 )
        
        running.number.matrix <- matrix( data=x$running.order, nrow=nrow( x$x.original ), ncol=N.features.per.plot, byrow=FALSE )
        
        for ( i in 1:N.plots ) {
                
                if ( i < N.plots ) {
                        idx.i <- ( ( i-1 )*( N.features.per.plot )+1 ):( i*N.features.per.plot )
                } else {
                        idx.i <- ( ( i-1 )*( N.features.per.plot )+1 ):( ncol( x$x.corrected ) )
                }
                
                matplot( x=running.number.matrix[ , 1:length( idx.i ) ], y=x$x.original[ , idx.i ], col=cols, pch=pchs, xlab="Sample running order", main="Original" )
                matplot( x=running.number.matrix[ which( x$is.QC ), 1:length( idx.i ) ], y=x$x.original[ x$is.QC, idx.i ], col=cols, pch=15, add=TRUE )
                
                matplot( x=running.number.matrix[ , 1:length( idx.i ) ], y=x$x.predicted[ , idx.i ], col=1, pch=pchs, type="l", add=TRUE )
                legend( x="top", legend=colnames( x$x.original )[ idx.i ], col=cols, pch=pchs )
                
                matplot( x=running.number.matrix[ , 1:length( idx.i ) ], y=x$x.corrected[ , idx.i ], col=cols, pch=pchs, xlab="Sample running order", main="Corrected" )
                matplot( x=running.number.matrix[ which( x$is.QC ), 1:length( idx.i ) ], y=x$x.corrected[ x$is.QC, idx.i ], col=cols, pch=15, add=TRUE )
                legend( x="top", legend=paste( colnames( x$x.corrected )[ idx.i ], "R.sq=", signif( x=x$R.sq[ idx.i ], digits=2 ) ), col=cols, pch=pchs )
                
        }
        
        return()
        
}
