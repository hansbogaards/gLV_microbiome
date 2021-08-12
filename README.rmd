---
title: R code "Correlations cannot distinguish interaction types in microbial networks"
author: "Susanne Pinto, Elisa Benincà, Egbert H. van Nes, Marten Scheffer and Johannes A. Bogaards"
output: 
  github_document:
    toc: true
    toc_depth: 3
editor_options:
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE, comment = NA, cache = TRUE,
  fig.width = 7, fig.height = 7)
```


## Introduction
R code belonging to the paper: “Correlations cannot distinguish interaction types in microbial networks” by Susanne Pinto, Elisa Benincà, Egbert H. van Nes, Marten Scheffer and Johannes A. Bogaards. This document shows the script used to perform the base case study of the paper.

We assessed how correlations between bacterial groups are shaped by their ecological interactions, and how correlations can be employed to signify these interactions. We specifically investigate how inference of microbial interactions is affected by interindividual variation in process parameters, and mainly focus on the necessary model assumptions and the meaning of the obtained correlation-based networks. We use the generalized Lotka-Volterra (gLV) model to demonstrate the performance of correlation-based network reconstruction.


## Load the packages
Load the packages in R.

```{r load the packages, echo = TRUE, results = FALSE}
library( tidyverse )
library( deSolve )
library( seqtime ) 
library( purrr )
library( reshape2 )
library( truncnorm )
library( Matrix )
library( ggplot2 )
```



## The generalized host-specific Lotka-Volterra model
We simulated bacterial communities by means of the generalized Lotka-Volterra (gLV) model. We simulated the population growth of species **i** with intrinsic growth rates (**r**), carrying capacities (**K**) and interactions with other species (**a**) that are specific for species **i** as well as for host **m** (see methods).


```{r the model, echo = TRUE}
gLV.model <- function( time, y, parms ) { 
  r_im = parms$r_im %>% as.numeric()
  K_im = parms$K_im %>% as.numeric() 
  a_ijm = parms$a_ijm %>% as.matrix() 
  
  # The function
  dydt <- r_im * y * ( 1 - ( 1/K_im ) * y + ( a_ijm %*% y ))
  
  # Store data in a list
  dydt %>% list() 
} 
```



## Choose your settings and draw random parameters
First set some standard settings. Optionally you can choose a different number of species and samples. 

```{r settings, echo = TRUE}
# Set date and time of the simulation  
Date <- Sys.time() %>% format( ., format="Date %d-%m-%Y_Time %H.%M" )

# Choose the number of simulations
# In the paper, 1000 simulations were used, here we choose 10 simulations, because it is more time efficient
n.simulations = 10

# Size of the population / total number of species
n.species = 10

# Number of samples (hosts)
n.samples = 300
# Later the samples with species that go extinct are removed, therefor simulate more samples than preferred 
n.samples.total = 400

# Timescale
tstart = 0
tend = 150
tstep = 1
times <- seq( from = tstart, to = tend, by = tstep ) 
```

Set your seeds for reproducibility.

```{r set seeds, echo = TRUE}
# Sample integers for seeds
seed <- 100000 %>% sample.int( ., n.simulations, replace = FALSE )
```

With the following script the growth rates, carrying capacities and the interaction matrix are computed randomly from a chosen distribution.

```{r random parameters, echo = TRUE}
# Make a list of n.simulations
parms <- n.simulations %>% vector( mode = "list", length = . )

# Do the simulation n.simulations * n.samples times
for ( j in 1:length( parms ) ) {
  # Make a list in the list of length N.samples
  parms[[j]] <- n.samples.total %>% vector( mode = "list", length = . )
  
  for ( k in 1:n.samples.total ) {
    
    # Use the seeds
    seed[j] %>% set.seed()
    
    # Set random parameters

    # initial abundances per species
    parms[[j]][[k]]$yini = n.species %>% runif( ., 0, 1 ) 
    # Intrinsic growth rates per species
    # Keep them positive to make all the species able to grow on their own
    parms[[j]][[k]]$r_i = n.species %>% runif( ., 0.05, 0.1 ) 
    # carrying capacities per species
    parms[[j]][[k]]$K_i = n.species %>% runif( ., 0, 1 ) 
    
    # First draw the interactions (a_ij) than fill the interaction matrix
    # Interactions are drawn from two truncated gaussian mixture distributions (which results in a bimodal distribution)
    negative_interactions <- ( ( n.species * n.species )/ 2 ) %>% 
      rtruncnorm( n = ., a = -0.5, b = 0.15, mean = -0.25, sd = 0.1 )
    positive_interactions <- ( ( n.species * n.species )/ 2 ) %>% 
      rtruncnorm( n = ., a = -0.15, b = 0.5, mean = 0.25, sd = 0.1 )
    
    # Determine the sparcity
    # The lower the denominator the sparcer the A matrix will be
    zero_index.negative_interactions <- sample( length( negative_interactions ) )[ 1:round(
      length( negative_interactions ) / 1.25 ) ]
    zero_index.positive_interactions <- sample( length( positive_interactions ) )[ 1:round(
      length( positive_interactions ) / 1.25 ) ]
    
    # Remove the values determined by the zero_index
    negative_interactions[ zero_index.negative_interactions ] <- 0
    positive_interactions[ zero_index.positive_interactions ] <- 0
    
    # Remove n.species/2 number of 0s because of the empty diagonal in the A matrix
    negative_interactions <- negative_interactions[ -sample( which( negative_interactions == 0 ), n.species/ 2 ) ]
    positive_interactions <- positive_interactions[ -sample( which( positive_interactions == 0), n.species/ 2 ) ]
    
    # Combine positive and negative interactions into a gaussian mixture distribution
    data.all <- c(negative_interactions, positive_interactions)
    # Randomize order
    data.all <- data.all %>% sample
    
    # Make an empty A matrix with diagonal 0 
   A.matrix <- matrix( NA, ncol = n.species, nrow = n.species )
    diag(A.matrix ) = 0
    
    # Fill in the A matrix (but not the diagonal)
   A.matrix <- data.all %>% ifelse( is.na(A.matrix ), ., 0 )

    # Set the row- and column names
    rownames(A.matrix ) <- "Species." %>% paste0( ., 1:n.species )
    colnames(A.matrix ) <- "Species." %>% paste0( ., 1:n.species )
    
    parms[[j]][[k]]$A.matrix <-A.matrix
    
  }
}
```

Visualize one of the interaction matrices

```{r visualize matrix, echo = TRUE}
# Visualize one interaction matrix
A.matrix_sparse <- parms[[1]][[1]]$A.matrix %>% Matrix(., sparse=TRUE)
A.matrix_sparse %>% image(main = "Interaction matrix (A matrix)", 
                          xlab = "Bacterial species", ylab = "Bacterial species")
```

We added multiplicative noise to the parameters to make them host specific. Here only the interspecies interactions are host specific, but it is possible to make the growth rates and carrying capacities also host-specific.

```{r multiplicative noise, echo = TRUE}
for ( l in 1:length( parms ) ) {
  for ( m in 1:n.samples.total ) {
    # Add noise to the interaction matrix A
    noise.on.a_ij <- parms[[l]][[m]]$A.matrix %>% length() %>% rnorm( ., 0, 0.25 )
    # Only add noise when the value in the A.matrix is not 0
    parms[[l]][[m]]$a_ijm <- ifelse( parms[[l]][[m]]$A.matrix != 0, parms[[l]][[m]]$A.matrix * exp( noise.on.a_ij ), 0 )
    
    # Add noise to the growth rates (in the base case this noise is 0)
    noise.on.r_i <- n.species %>% rnorm( ., 0, 0 )
    parms[[l]][[m]]$r_im <- parms[[l]][[m]]$r_i * exp( noise.on.r_i )
    
    # Add noise to the carrying capacities (in the base case this noise is 0)
    noise.on.K_i <- n.species %>% rnorm( ., 0, 0 )
    parms[[l]][[m]]$K_im <- parms[[l]][[m]]$K_i * exp( noise.on.K_i )
  }
}
```



## Solve the model
The gLV model was solved with the lsoda function (from the deSolve package version 1.24). Because, most gut microbiota studies are limited to only a few samples in time, presenting mere ‘snapshots’ of the intestinal ecosystem, we only 'took' one sample of the species abundances in the equilibrium. 


```{r solving the model, echo = TRUE, results = "hide" }
# Make a list of n.simulations to store the time series
community.times <- n.simulations %>% vector( mode = "list", length = . )
# Make a list of n.simulations to store the samples taken during equilibrium
community.samples <- n.simulations %>% vector( mode = "list", length = . )

# Do for every simulation
for ( n in 1:length( community.times ) ) {
  # Make a list in the list of length N.samples
    community.times[[n]] <- n.samples.total %>% vector( mode = "list", length = . )
    community.samples[[n]] <- matrix( NA, nrow = n.species, ncol = n.samples.total )
  
  tryCatch( {
    
    # Do for every sample / individual
    for ( o in 1:length( community.times[[n]] ) ) {
      
      # Random yini than yini = N.species.total %>% runif()
      community.time = lsoda( y = parms[[n]][[o]][ c( "yini" ) ] %>% unlist %>% as.numeric, 
                              times = times, 
                              func = gLV.model, 
                              parms = parms[[n]][[o]] )
      
      # First column of community time dataset is time
      time = community.time[,1] 
      # Second column till the last species is results oda
      community.time = community.time[, 2:ncol( community.time ) ] 
      # Transpose table
      community.time = t( community.time ) 
      # Use times as columnnames
      colnames( community.time ) = times
      
      # Change the data to a data frame
      abundances <- community.time %>% as.data.frame()
      # Get the name of the last column
      keep <- c( tend + 1 )
      # Keep only the abundances at the last time point (in the equilibrium)
      abundances <- abundances[keep] %>% simplify2array()
      
      # Store the abundances at the last time point in the list
      community.samples[[n]][,o] <- abundances
    }
    # Do not stop the simulation because of an error or warning (but do not return any results)
  }, warning = function( e ){ return( NULL ) } )
  # Print number of simulation, to see progress
  n %>% print()
}

```



Choose a community and plot to see the time series.

```{r plot the community, echo = TRUE}
# Show just one figure
par( mfrow = c( 1, 1 ) )
# Make a time plot (from the package seqtime)
community.time %>% tsplot( ., legend = T, title = "A host specific timeseries" )
```


## Create the dataset
We added additive noise, drawn from a uniform distribution, to the abundances sampled in equilibrium to represent uncertainty in measurements.

```{r measurement noise, echo = TRUE  }
# Make a list of n.simulations
community.samples.withnoise <- n.simulations %>% vector( mode = "list", length = . )

# Do for every simulation
for (s in 1:length( community.samples.withnoise )) {
  # Create the level of noise
  measurement.noise <- community.samples[[s]] %>% length() %>% 
    runif( ., -0.01, 0.01 ) %>% 
    matrix( ., nrow = n.species, ncol = n.samples.total )
  # Add additive noise to the dataset
  community.samples.withnoise[[s]] <- community.samples[[s]] + measurement.noise
  
  # Create column- and rownames
  nms.of.samples = "Sample." %>% paste0( ., 1:n.samples.total )
  nms.of.species = parms[[1]][[1]]$A.matrix %>% colnames()
  # Name columns and rows
  colnames( community.samples.withnoise[[s]] ) = nms.of.samples 
  rownames( community.samples.withnoise[[s]] ) = nms.of.species
  }
```

If species did not co-exist (when the abundance of a species drops below 0.001), we rejected the simulation. 

```{r remove species that got extinct, echo = TRUE  }
# Make a list of n.simulations
community.samples.survival <- n.simulations %>% vector( mode = "list", length = . )

for ( t in 1:length( community.samples.survival ) ) {
  # Transpose the dataset
  community.samples.transpose <- community.samples.withnoise[[t]] %>% t()
  # Change the values of species with an abundance below 0.01 to NA
  community.samples.survive <- ifelse( community.samples.transpose < 0.001, NA, 
                                       community.samples.transpose )
  community.samples.survive <- community.samples.survive %>%  as.data.frame()
  # Remove columns with all NAs
  community.samples.survive <- Filter( function( x ) !all( is.na( x ) ), community.samples.survive )
  # Remove rows with some NAs
  community.samples.survive <- community.samples.survive[complete.cases( community.samples.survive ), ] %>% data.frame()
  
  # Number of species left
  species.left <- ncol( community.samples.survive ) 
  # Number of samples left
  samples.left <- nrow( community.samples.survive ) 
  
  survival <- c( species.left, samples.left ) %>% print()
  
  # Remove the simulation if not enough species and samples are left
  if ( species.left < n.species  | samples.left < n.samples ) {
    print( "not enough species survived" )
  } else {
    community.samples.survival[[t]]$data <- community.samples.survive %>% as.matrix()
  }
}
```

```{r remove samples with not enough species, echo = TRUE }
# Check if there are simulations with not enough species or samples
removed.simulations.survival_error <- n.simulations %>% character()
for ( u in 1:length( community.samples.survival ) ) {
  # If errors occured, print number
  if( length( unlist( community.samples.survival[[u]] ) ) == 0 ) {
    # Store the numbers of the simulations with errors
    removed.simulations.survival_error[[u]] <- print( u )
  } else {
    removed.simulations.survival_error[[u]] <-print( NA )
  }
}

# Remove empty list elements
community.samples.survival <- community.samples.survival %>% compact()
```


The samples that are removed in the previous step, also needs to be removed from the parms data.

```{r remove samples and species, echo = TRUE, error = TRUE  }

# Removed simulations after extinction of species
removed.simulations <- removed.simulations.survival_error %>% as.numeric()

# Keep only the parms that were not removed (== NA) 
for (v in 1:length( parms )) {
  if ( removed.simulations[[v]]  %>% is.na() ) {
    parms[[v]] <- parms[[v]]
  }   else {
    parms[[v]] <- list( NULL )
  }
}

# Remove the simulations where errors or not enough samples/species survived have occurred
parms <- Filter( function( x ) any( unlist( x ) != 0 ), parms )
```

We choose a higher number of samples than the preferred number we want to end up with. Therefore, randomly remove some of the samples to end up with the preferred number of samples in which all species coexist.

```{r update parms, echo = TRUE}
# Make a list of n.simulations
data <- length( community.samples.survival ) %>% vector( mode = "list", length = . )

# Keep randomly the preferred number of samples
for ( w in 1:length( data ) ) {
  # Number of samples left
  samples.left <- nrow( community.samples.survival[[w]][["data"]] ) 
  # Get the number of samples that needs to be removed
  n.samples.to.remove = samples.left - n.samples
  # Transpose the data
  data.samples <- community.samples.survival[[w]][["data"]] %>% t() %>% as.data.frame()
  # Choose the samples that will be removed
  samples.to.remove <- data.samples %>% 
    sample( ., n.samples.to.remove )
  # Get the names of the samples to be removed
  names.samples.to.remove <- samples.to.remove %>% colnames()
  
  # Select the samples to keep
  data[[w]]$data <- data.samples %>% 
    select( ., -one_of( names.samples.to.remove ) )
}

# Transpose the data
for ( x in 1:length( data ) ) {
  data[[x]]$data <- data[[x]][["data"]] %>% t() %>% as.data.frame()
}

# Number of species
data[[1]][["data"]] %>% ncol()
# Number of samples
data[[1]][["data"]] %>% nrow()

# Make a list of n.simulations
names <- length( data ) %>% vector( mode = "list", length = . )

# Do for every sample / individual
for ( y in 1:length( data ) ) {
  names[[y]][["Samples"]] <- rownames( data[[y]][["data"]] )
  names[[y]][["Species"]] <- colnames( data[[y]][["data"]] )
}

```


## Calculate the partial correlations
First, we calculate the partial correlations.

```{r partial correlations, echo = TRUE  }
# Make a list of length  data
results.pcor <- length( data ) %>% vector( mode = "list", length = . )

for ( i in 1:length( data ) ) {
  # Calculate the covariance matrices from the abundance tables
  cx <- data[[i]][["data"]] %>% as.data.frame() %>% cov()
  
  # Calculate the inverse of the covariance matrix
  dx <- cx %>% solve( ., sparse = TRUE,  tol = 1e-150 )
  # Calculate the partial correlation matrix
  pcor <- -cov2cor( dx )
  
  # Change the diagonal to 1 (instead of -1)
  diag( pcor ) <- 1 
  
  # Store the results in the list
  results.pcor[[i]]$pcor <- pcor %>% as.matrix()
}
```

Hereafter we should keep only the significant partial correlations (corrected with the Benjamini Hochberg procedure). These significant correlation matrices are transformed to a matrix where, negative correlations are assigned a "-1", positive correlations are assigned a "1" and if a correlations was not significant than the value "0" is implemented.

```{r significant correlations, echo = TRUE  }
# Do for all the pcors
pcor.pvalues.significant <- length( results.pcor ) %>% vector( mode = "list", length = . )

for ( i in 1:length( results.pcor ) ) {
  # Number of samples
  pcor.n <- dim( data[[i]][["data"]] )[1]
  pcor.gp <- dim( data[[i]][["data"]] )[2]-2
  
  # Do the statistics
  pcor.statistic <- results.pcor[[i]][["pcor"]] * 
    sqrt( ( pcor.n-2-pcor.gp ) / ( 1-results.pcor[[i]][["pcor"]]^2 ) )
  # Get the p values
  pcor.pvalues <- 2 * pnorm( -abs( pcor.statistic ) )
  
  # Adjust the p-values with Benjamini Hochberg procedure
  pcor.pvalues.adjusted <- array(NA, dim=dim(pcor.pvalues))
  names <- colnames(pcor.pvalues)
  colnames(pcor.pvalues.adjusted) <- names
  rownames(pcor.pvalues.adjusted) <- names
  for (j in 1:dim(pcor.pvalues)[1])
    pcor.pvalues.adjusted[j,] <- p.adjust(pcor.pvalues[j,], method = "BH")
  
  # Keep the significant results
  pcor.pvalues.sign <- pcor.pvalues.adjusted < 0.05
  # Change true in 1 and false to 0
  pcor.pvalues.sign <- 1 * pcor.pvalues.sign
  
  # Change the matrices to negative (-1), no (0) and positive interaction (-1) metrices
  pcor.pvalues.sign <- ifelse( pcor.pvalues.sign == 0, 0,
                               ifelse( results.pcor[[i]][["pcor"]] >= 0, 1, -1 ) )
  
  # Remove the diagonal
  diag( pcor.pvalues.sign ) <- NA
  
  # Store in the results
  pcor.pvalues.significant[[i]] <- pcor.pvalues.sign
}
```



## Calculate precision, recall and the F1-scores
These measures are used to compare the original interaction matrix from the model with the inferred significant partial correlations. 

```{r prepare the data, echo = TRUE}
# Make a list 
A.matrix.analysis <- list() 

# Subtract all the interaction matrices
for ( i in 1:length( parms ) ) {
  # Make a list per simulation
 A.matrix.analysis[[i]] <- length( parms[[i]] ) %>% vector( mode = "list", length = . )
  for ( j in 1:length( parms[[i]] ) ) { 
    # Here the a matrix without noise is used
   A.matrix.analysis[[i]]<- parms[[i]][[j]][ c( "A.matrix" ) ] 
  }
}

# Flatten the list of lists
A.matrix.analysis <-A.matrix.analysis %>% flatten()

```

After selecting the variables of interest, the interaction matrix is transformed in a similar way as the significant partial correlation matrix, with the -1, 0 and 1 for respectively negative, none and positive interactions.

```{r transform the interaction matrix, echo = TRUE}
# Do for the a.matrices
A.matrix.round <- length(A.matrix.analysis ) %>% vector( mode = "list", length = . )

for ( i in 1:length(A.matrix.round ) ) {
  # change the matrices to negative (-1), no (0) and positive interaction (-1) metrices
 A.matrix.analysis.01 <- ifelse(A.matrix.analysis[[i]] == 0, 0,
                                  ifelse(A.matrix.analysis[[i]] >= 0, 1, -1 ) )
  
  # Remove the diagonal, because that has nothing to do with interactions
  diag(A.matrix.analysis.01 ) <- NA
  
  # Store in the results
 A.matrix.round[[i]] <-A.matrix.analysis.01
}
```

Calculate the confusion matrix.

```{r confusion matrix, echo = TRUE}
# Make a new list
confusion.matrix <- length(A.matrix.round) %>% vector( mode = "list", length = . )

# Do for every simulation
for ( i in 1:length(A.matrix.round ) ) {
  
  # Make an actual and predicted matrix
  actual <-A.matrix.round[[i]]
  predicted <- pcor.pvalues.significant[[i]]
  
  # Change the lower half of the predicted matrix to na (because it is symmetric)
  predicted[ lower.tri( predicted ) ] <- NA
  
  # Make a empty confusion matrix to store the results
  # Actual is rows, predicted is columns
  results <- matrix( 0, 3, 3 )
  colnames( results ) <- c( "-1","0", "1" )
  rownames( results ) <- c( "-1","0", "1" )
  
  # Go through every row and columns
  for ( r in 1:nrow( predicted ) ) {   
    for ( c in 1:ncol( predicted ) ) {
      
      # If the value is na, than it was a diagonal and do nothing
      if( is.na( predicted[r, c] ) ) {
        results[1, 1] = results[1, 1] + 0
        
        # If cor is 0 and interactions are both 0
      } else if ( predicted[r, c] == 0 & ( actual[r, c] == 0 & actual[c, r] == 0 ) ) {
        # Than a true negative
        results[2, 2] = results[2, 2] + 1
        # If cor is 0 and  one of the interactions interactions are 1
      } else if ( predicted[r, c] == 0 & ( actual[r, c] == 1 | actual[c, r] == 1 ) ) {
        # Than a false negative
        results[3, 2] = results[3, 2] + 1
        # If cor is 0 and  one of the interactions interactions are -1
      } else if ( predicted[r, c] == 0 & ( actual[r, c] == -1 | actual[c, r] == -1 ) ) {
        # Than a false negative
        results[1, 2] = results[1, 2] + 1
        
        # If cor is 1 and one of the interactions interactions are 1
      } else if ( predicted[r, c] == 1 & ( actual[r, c] == 1 | actual[c, r] == 1 ) ) {
        # Than a true positive
        results[3, 3] = results[3, 3] + 1
        # If cor is 1 and interactions are both 0
      } else if ( predicted[r, c] == 1 & ( actual[r, c] == 0 & actual[c, r] == 0 ) ) {
        # Than a false positive
        results[2, 3] = results[2, 3] + 1
        # If cor is 1 and interactions are both -1
      } else if ( predicted[r, c] == 1 & ( actual[r, c] == -1 & actual[c, r] == -1 ) ) {
        # Than a false positive (wrong sign)
        results[1, 3] = results[1, 3] + 1
        
        # If cor is -1 and one of the interactions are  -1
      } else if ( predicted[r, c] == -1 & ( actual[r, c] == -1 | actual[c, r] == -1 ) ) {
        # Than a true positive
        results[1, 1] = results[1, 1] + 1
        # If cor is -1 and interactions are both 0
      } else if ( predicted[r, c] == -1 & ( actual[r, c] == 0 & actual[c, r] == 0 ) ) {
        # Than a false positive
        results[2, 1] = results[2, 1] + 1
        # If cor is -1 and interactions are both 1
      } else if( predicted[r, c] == -1 & ( actual[r, c] == 1 & actual[c, r] == 1 ) ) {
        # Than a false positive (wrong sign)
        results[3, 1] = results[3, 1] + 1
      }
      # Store the results in the list
      confusion.matrix[[i]] <- results
    }
  }
}
```

Calculate the precision, recall and F1-scores

```{r scores,  , echo = TRUE}
# Make a list to store the results
F1.score <- length( confusion.matrix ) %>% vector( mode = "list", length = . )

# Calculate F1_score with paying attention to the sign
for ( i in 1:length( confusion.matrix ) ) {
  confusion <- confusion.matrix[[i]]
  
  # Calculate the precision
  precision = ( confusion[1, 1] + confusion[3, 3] ) / ( ( confusion[1, 1] + confusion[3, 3] ) + 
                                                          ( confusion[2, 1] + confusion[3, 1] + confusion[1, 3] + confusion[2, 3] ) )
  F1.score[[i]]$precision <- precision
  
  # Calculate the recall
  recall = ( confusion[1,1] + confusion[3,3] ) / ( ( confusion[1,1] + confusion[3,3] ) + 
                                                     ( confusion[1,2] + confusion[3,2] ) )
  F1.score[[i]]$recall <- recall
  
  # Calculate the F1-scores
  F1score <- 2 * ( ( precision * recall ) / ( precision + recall ) )
  F1.score[[i]]$F1score <- F1score
}
```

Plot the results 
```{r plot the results, echo = TRUE}
# Make an empty dataframe
F1.score.data <- matrix( 0, nrow = length( F1.score ), ncol = 3 )

# Do for every simulation
# Make a dataframe with the results
for ( i in 1:length( F1.score ) ) {
  F1.score.data[i, 1] <- F1.score[[i]][[ "precision" ]]
  F1.score.data[i, 2] <- F1.score[[i]][[ "recall" ]]
  F1.score.data[i, 3] <- F1.score[[i]][[ "F1score" ]]
}

# Change the NA's in the data to 0
F1.score.data <- F1.score.data %>% replace( ., is.na( F1.score.data ), 0 ) 
# Name the columns
colnames( F1.score.data ) <- c( "Precision", " Recall", " F1.score" )

F1.score.data.2 <- F1.score.data %>% melt( ., id.vars = "F1.score" )
str( F1.score.data.2 )
F1.score.data.2$value <- as.numeric( as.character( F1.score.data.2$value ) )

# Make the boxplot
boxplot <- ggplot( data = F1.score.data.2, 
                   mapping = aes( x = Var2, y = value ) ) +
  geom_boxplot() +
  # Name the labels
  labs( x = "Score type", y = "Value" ) +
  # Text over an angle
  theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )  + 
  # Color the boxes
  geom_boxplot(
    aes( fill = Var2 ),
    position = position_dodge( 0.9 ) ) + 
  ggtitle( paste( "Results", sep = ",", collapse = NULL) )  +
  guides(fill=guide_legend(title="Score type")) +
  scale_y_continuous( limits = c( 0, 1 ) )
# Plot the boxplot
boxplot
```
