##### Functions Used #####

rss_ds <- function(yi.less, yi.larg, y.less.avg, y.larg.avg) {
    
    # Function to find the training RSS of the labels y = {yi.less, yi.larg} on DS
    #
    # Input: yi.less (vector of yi for x < s), 
    #        yi.larg (vector of yi for x >= s),  
    #        y.less.avg (avg training y with x < s), 
    #        y.larg.avg (avg training y with x >= s)
    #
    # Output: rss
    
    rss <- sum((yi.less - y.less.avg)^2) + sum((yi.larg - y.larg.avg)^2)
    
    return(rss)
}


ds_alg <- function(X, y) {
    
    # Function to find the optimal decision stump using a basic grid-search
    # method between the minimum and maximum values of exactly two attributes
    #
    # Input: X (training data in an nx2 matrix; row i, col j), y (training labels)
    # 
    # Output: c(attrib, s, y.less, y.larg) (vector containing attribute j, threshold s, 
    # the average y for training observations with xij < s (y.less), 
    # and similar but with xij >= s (y.larg). These define the decision stump)
    
    # Create the grid to search over
    attrib1.arr <- seq(min(X[,1]), max(X[,1])+0.1, by = 0.1)
    attrib2.arr <- seq(min(X[,2]), max(X[,2])+0.1, by = 0.1)
    
    # Initialise
    attrib <- 0
    smallest.rss <- 0
    
    # Attribute 1
    for (trial.s in attrib1.arr) {
        # Find y.less, y.larg
        yi.less <- y[which(X[,1] < trial.s)]
        yi.larg <- y[which(X[,1] >= trial.s)]
        y.less.avg <- mean(yi.less)
        y.larg.avg <- mean(yi.larg)
        
        # Find training RSS and keep it with values for the DS if RSS 
        # is smaller than the current smallest RSS (or if it's the first run)
        train.rss <- rss_ds(yi.less, yi.larg, y.less.avg, y.larg.avg)
        if (train.rss < smallest.rss | attrib == 0) {
            smallest.rss <- train.rss
            attrib <- 1
            s <- trial.s
            y.less <- y.less.avg
            y.larg <- y.larg.avg
        }
    }
    
    # Attribute 2
    for (trial.s in attrib2.arr) {
        # Find y.less, y.larg
        yi.less <- y[which(X[,2] < trial.s)]
        yi.larg <- y[which(X[,2] >= trial.s)]
        y.less.avg <- mean(yi.less)
        y.larg.avg <- mean(yi.larg)
        
        # Find training RSS and keep it with values for the DS if RSS 
        # is smaller than the current smallest RSS
        train.rss <- rss_ds(yi.less, yi.larg, y.less.avg, y.larg.avg)
        if (train.rss < smallest.rss) {
            smallest.rss <- train.rss
            attrib <- 2
            s <- trial.s
            y.less <- y.less.avg
            y.larg <- y.larg.avg
        }
    }
    
    return(c(attrib, s, y.less, y.larg))
}

ds_pred <- function(X, ds.params) {
    
    # Function to predict labels y.hat from data X using the input ds.params
    #
    # Input: X (data; only two attributes), 
    #        ds.params (parameters c(attrib, s, y.less, y.larg)
    #                   defining the decision stump)
    #
    # Output: y.hat (predicted labels)
    
    # Unpack ds.params
    attrib <- ds.params[1]
    s <- ds.params[2]
    y.less <- ds.params[3]
    y.larg <- ds.params[4]
    
    # Init a double vector
    y.hat <- double(nrow(X))
    
    # Set predictions yi.hat (depending on the branch used, or on attrib.i) 
    # to the avg training y at that node
    y.hat[which(X[,attrib] < s)] <- y.less
    y.hat[which(X[,attrib] >= s)] <- y.larg
    
    return(y.hat)
}

bds_alg <- function(X, y, learn.rate, B) {
    
    # Function to apply boosting to the decision stump algorithm
    #
    # Input: X (training data; 2 attributes only), y (training labels), 
    #        learn.rate (learning rate), B (the number of trees used)
    #
    # Output: fb.arr (array of DSs (fb), which define the BDS prediction rule)
    
    # Initialise fb.arr and set r to y for the first iteration (f.hat = 0)
    fb.arr <- matrix(0, nrow=B, ncol=4)
    r <- y
    
    # Repeatedly find DSs for predicting the residuals
    # The residuals are updated to give a new residual to predict
    for (b in 1:B) {
        fb.arr[b,] <- ds_alg(X, r)
        r <- r - learn.rate*ds_pred(X, fb.arr[b,])
    }
    
    return(fb.arr)
}

bds_pred <- function(X, fb.arr, learn.rate, B) {
    
    # Function to predict labels y.hat from data X using the input fb.arr
    # and learn.rate
    #
    # Input: X (data; only two attributes),
    #        fb.arr (array of parameters defining B DSs),
    #        learn.rate (learning rate),
    #        B (number of trees)
    #
    # Output: y.hat (predicted labels)
    
    # Find the sum of the predictions for the residuals
    summ.fb.hat <- double(nrow(X))  # f.hat = 0
    for (b in 1:B) {
        summ.fb.hat <- summ.fb.hat + ds_pred(X, fb.arr[b,])
    }
    
    # Multiply by learning rate to get the overall prediction rule
    y.hat <- learn.rate*summ.fb.hat
    
    return(y.hat)
}

mse <- function(y, y.hat) {
    # Function to find the out.mse given y, y.hat
    #
    # Input: y (actual labels), y.hat (predicted labels)
    #
    # Output: out.mse (output MSE)
    
    out.mse <- sum((y - y.hat)^2)/length(y)
    
    return(out.mse)
}

############################################################################

##### Code to answer the tasks #####

# Set-up:

library(MASS)
data(Boston)
attach(Boston)

# Get X (lstat and rm only)
X <- matrix(c(lstat, rm), nrow = nrow(Boston), ncol = 2)

# Split into equally sized training and test sets
set.seed(0912)
train <- sample(1:nrow(Boston), nrow(Boston)%/%2)  # %/% = integer division
X.train <- X[train,]
X.test <- X[-train,]
y.train <- medv[train]
y.test <- medv[-train]


# Task 1:

# Fit using training set
ds.out <- ds_alg(X.train, y.train)
# Predict on test set
y.hat <- ds_pred(X.test, ds.out)
# Test MSE
mse.out <- mse(y.test, y.hat)


# Task 2: 

# Init
learn.rate <- 0.01
B <- 1000

# Fit using training set
fb.arr <- bds_alg(X.train, y.train, learn.rate, B)
# Predict on test set
y.hat <- bds_pred(X.test, fb.arr, learn.rate, B)
# Test MSE
mse.out <- mse(y.test, y.hat)


# Task 3:

# Init
learn.rate <- 0.01
B.arr <- c(1, seq(10, 100, by=10), seq(200, 1000, by=100), seq(2000, 10000, by=1000), 100000)
mse.out <- double(length(B.arr))

# Loop over trees
for (i in 1:length(B.arr)) {
    # Fit Model
    fb.arr <- bds_alg(X.train, y.train, learn.rate, B.arr[i])
    # Predict
    y.hat <- bds_pred(X.test, fb.arr, learn.rate, B.arr[i])
    # Find Test MSE
    mse.out[i] <- mse(y.test, y.hat)
}

# Show Arrays
mse.out
B.arr

# Plot
plot(B.arr[12:29], mse.out[12:29], main='Test MSE v Number of Trees', xlab='Number of Trees (B)', ylab='Test MSE', type='b', col='blue')



