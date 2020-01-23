# Q1 Code:
log_F <- function(x) {
    
    # Function to calculate sigma(x), the logistic function applied to x
    #
    # Input: x (double)
    #
    # Output: sigma(x) (double)
    
    sigma_x = exp(x) / (1 + exp(x))
    
    return(sigma_x)
}

logP_F <- function(x) {
    
    # Function to calculate sigma'(x), the logistic function differentiated 
    # with respect to x, applied to input x
    #
    # Input: x (double)
    #
    # Output: sigma'(x) (double)
    
    sigmaP_x = exp(x) / (1 + exp(x))^2
    
    return(sigmaP_x)
}

gradDes_logReg <- function(X, y, B, iter, n) {
    
    # Function to apply gradient decent to logistic regression
    #
    # Input: X = training set (matrix), y = training labels (vector)
    #        B = initial random weights (vector), iter = training steps (int), 
    #        n = learning rate (double)
    #
    # Output: Optimised weights, B
    
    numSamples = nrow(X)
    numAttribs = ncol(X)
    p_X = double(numSamples)  # vector of doubles, length = numSamples
    updateC = double(numSamples)
    for (i in 1:iter) {
        # FORWARD: Find p(X_j) using current weights beta_j
        for (j in 1:numSamples) {
            p_X[j] = log_F(B[1] + sum(B[-1] * X[j,]))
        }
        
        # BACKWARD (1): Calculate 'update constant' and 'learning rate constant'
        for (j in 1:numSamples) {
            updateC[j] = logP_F(B[1] + sum(B[-1] * X[j,]))*(y[j] - p_X[j])
        }
        learnRateC = 2*n/numSamples
        
        # BACKWARD (2): Update wieghts beta_j
        summ = 0
        for (k in 1:numSamples) {
            summ = summ + updateC[k]
        }
        B[1] = B[1] + learnRateC*summ
        for (j in 2:numAttribs+1) {
            summ = 0
            for (k in 1:numSamples) {
                summ = summ + updateC[k]*X[k,j-1]
            }
            B[j] = B[j] + learnRateC*summ
        }
    }
    
    return(B)
}

mse <- function(X, y, B) {
    
    # Function to find the MSE of the predictions (for probability of 1)
    # 
    # Input: X = set of samples (matrix), y = labels (vector), 
    #        B = weights (vector)
    # 
    # Output: mse_out (double)
    
    numSamples = nrow(X)
    summ = 0
    for (i in 1:numSamples) {
        summ = summ + (y[i] - log_F(B[1] + sum(B[-1] * X[i,])))^2
    }
    mse_out = summ/numSamples
    
    return(mse_out)
}

pred_logReg <- function(X, B, minYes) {
    
    # Function to predict 1/0 based on the predicted probabilities of 1
    # (not needed, but it's nice to see these results)
    # 
    # Input: X = set of samples (matrix), B = weights (vector), 
    #        minYes = the minimum for deciding a prediction of yes [0,1]
    # 
    # Output: y.pred = vector of predictions
    
    # Get probabilities
    numSamples = nrow(X)
    p_X = double(numSamples)
    for (i in 1:numSamples) {
        p_X[i] = log_F(B[1] + sum(B[-1] * X[i,]))
    }
    
    # Create predictions
    y.pred <- rep(0, numSamples)
    y.pred[p_X > minYes] <- 1
    
    return(y.pred)
}

# Set Up; loading of data, standardizing, splitting, etc. (not shown here)

# Q4 Code: 
testMat = matrix(0, nrow = 5, ncol = 5)
trainMat = matrix(0, nrow = 5, ncol = 5)
i = 1
for (iter in c(200, 1000, 5000, 10000, 20000)) {  # cols
    for (n in c(0.01, 0.05, 0.1, 0.3, 1.0)) {  # rows
        B = gradDes_logReg(X.train, y.train, B_init, iter, n)
        testMat[i] = mse(X.test, y.test, B)
        trainMat[i] = mse(X.train, y.train, B)
        i = i + 1
    }
}

# Labels
colnames(testMat) <- c('TS = 200', 'TS = 1,000', 'TS = 5,000', 'TS = 10,000', 'TS = 20,000')
rownames(testMat) <- c('LR = 0.01', 'LR = 0.05', 'LR = 0.1', 'LR = 0.3', 'LR = 1.0')
colnames(trainMat) <- c('TS = 200', 'TS = 1,000', 'TS = 5,000', 'TS = 10,000', 'TS = 20,000')
rownames(trainMat) <- c('LR = 0.01', 'LR = 0.05', 'LR = 0.1', 'LR = 0.3', 'LR = 1.0')

# Print Tables
testMat.table <- as.table(testMat)
trainMat.table <- as.table(trainMat)
print(testMat.table)
print(trainMat.table)

# Q6 Code: 
iter = 5000
n = 0.1
mseVec = double(100)
for (i in 1:100) {
    B_init = runif(7, -0.7, 0.7)
    B = gradDes_logReg(X.train, y.train, B_init, iter, n)
    mseVec[i] = mse(X.test, y.test, B)
}

boxplot(mseVec, ylab = 'Test MSE', main = 'Test MSEs For Different Initial Weights')

