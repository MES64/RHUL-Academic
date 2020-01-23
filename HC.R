euclid_dist <- function(x, y) {
    # Function to find the Euclidean distance between two vectors x and y
    #
    # Input: x, y (vectors)
    #
    # Output: d (distance)
    
    vec_len <- length(x)
    summ <- 0
    for (i in 1:vec_len) {
        summ <- summ + (x[i] - y[i])^2
    }
    d <- sqrt(summ)
    
    return(d)
}

hierarch_agglom_clust <- function(X, linkage) {
    # Function to perform heirarchical agglomerative clustering given input 
    # data X, and returning the dendrogram and a list listing the clusters for
    # every number of clusters k = 1,...,n. 
    #
    # Input: X (matrix of data),
    #        linkage (str from {'single', 'complete', 'average', 'centroid'}; 
    #                 how to measure dissimilarity)
    #
    # Ouput: dend (dendrogram),
    #        clusters (list of cluster decompositions, with index k=1,...n
    #                  giving the number of clusters. Each observation is 
    #                  represented as the row index in X). 
    
    
    # LISTS:
    # List of lists: list(list())
    # Unlist() to use functions on arrays
    # Coerce to a list with as.list()
    # Merge two lists into one using c(list1, list2)
    # Create a list of lists using list(list1, list2)
    # Use X[[i]] to access actual element (not a list like in X[i])
    # X[[1]][[2]] = X[[c(1,2)]]
    # Use str(X[i]) to output elements nicely
    
    # STEPS:
    # 0: Pre-calculate a distance matrix between all points (unless centroid)
    # 1: Loop through each pair of clusters (element of overall list)
    # 2: Find distances, using the indexes of the points to look them up
    # 3: Get a matrix of distances for these two clusters
    # 4: Indetify pair of clusters with least dissim. (min func. or loop)
    # 5: Merge these clusters (e.g. 1,3) by: X <- list(1,2,3,4)
    #                                        X[[1]] <- c(X[1],X[3])
    #                                        X <- X[1:4 != 3]
    # 6: Repeat until all clusters are merged (i = n, n-1, ..., 1)
    
    # Create object of type dendrogram; a list of lists ... with attributes
    # Use plot(dend, center = TRUE); which puts the branches at the midpoint
    # dend = list of list of list of ... as above. Add attributes per list:
    # attr(dend[[i]], 'members') <- num of members in that sub-tree (all)
    # attr(dend[[i]], 'height') <- measure of dissim. (all)
    # attr(dend[[i]], 'leaf') <- TRUE (only add if it is a leaf)
    # attr(dend[[i]], 'class') <- "dendrogram" (all)
    
    # - Start with n clusters with attribs. members = 1, height = 0, 
    # leaf = TRUE and class = "dendrogram"
    # - With every merge, add attribs. members = sum of members of two lists,
    # height = dissim., and class = "dendrogram". 
    # - End up with a dendrogram object to plot. 
    
    # Also, create a list of vectors for each merge (not a dendrogram object,
    # just a list of vectors per cluster). This can produce scatter plots like:
    # plot(x1, y1, col=1, pch=20, xlim=c(-1, 4), ylim=c(-1,4))
    # points(x2, y2, col=2, pch=20)
    # points(x3, y3, col=3, pch=20)
    # ... or tables comparing the clusters to the true labels if they exist
    
    # Initialise the lists
    # Clusters
    X_rows <- nrow(X)
    clusters <- list(double(X_rows))
    clusters[[X_rows]] <- as.list(1:X_rows)
    # Dendrogram
    dend <- as.list(1:X_rows)
    for (i in 1:X_rows) {
        attr(dend[[i]], 'members') <- 1
        attr(dend[[i]], 'height') <- 0
        attr(dend[[i]], 'leaf') <- TRUE
        attr(dend[[i]], 'class') <- "dendrogram"
    }
    
    # Pre-calculate distances
    if (linkage != 'centroid') {
        X_dist <- matrix(0, nrow=X_rows, ncol=X_rows)
        for (j in 1:(X_rows-1)) {
            for (i in (j+1):X_rows) {
                X_dist[i,j] <- euclid_dist(X[i,], X[j,])
                X_dist[j,i] <- X_dist[i,j]
            }
        }
    } else {  # If centroid linkage, initialise centroids vector
        centroids <- X  # row = centroid; col = coords
    }
    
    # Loop until all the clusters are merged
    for (k in (X_rows-1):1) {
        prev_clusters <- clusters[[k+1]]
        num_clust <- length(prev_clusters)
        
        # Create dissim. matrix
        dissim <- matrix(Inf, nrow=num_clust, ncol=num_clust)
        # Loop through each of the cluster pairs
        for (j in 1:(num_clust-1)) {
            for (i in (j+1):num_clust) {
                if (linkage != 'centroid') {
                    # Find dissim. per cluster-pair combo
                    clust1 <- prev_clusters[[i]]
                    clust2 <- prev_clusters[[j]]
                    clust_dist <- X_dist[clust1, clust2]
                    if (linkage == 'single') {
                        dissim[i,j] <- min(clust_dist)
                    } else if (linkage == 'complete') {
                        dissim[i,j] <- max(clust_dist)
                    } else if (linkage == 'average') { 
                        dissim[i,j] <- mean(clust_dist)
                    }
                } else {  # If centroid linkage
                    cent1 <- centroids[i,]
                    cent2 <- centroids[j,]
                    dissim[i,j] <- euclid_dist(cent1, cent2)
                }
            }
        }
        # Find the clusters with the min. dissim.
        min_dissim <- min(dissim)
        dissim_indx <- which(dissim == min_dissim, arr.ind = TRUE)
        indx <- dissim_indx[1]
        jndx <- dissim_indx[2]
        
        # Merge these two clusters
        # Centroids
        if (linkage == 'centroid') {
            clust1 <- prev_clusters[[indx]]
            clust2 <- prev_clusters[[jndx]]
            len_clust1 <- length(clust1)
            len_clust2 <- length(clust2)
            if (len_clust1 > 1) {
                clust1_sum <- colSums(X[clust1,])
            } else {
                clust1_sum <- X[clust1,]
            }
            if (len_clust2 > 1) {
                clust2_sum <- colSums(X[clust2,])
            } else {
                clust2_sum <- X[clust2,]
            }
            centroids[indx,] <- (clust1_sum + clust2_sum)/(len_clust1+len_clust2)
            centroids <- centroids[1:num_clust != jndx,]
        }
        # Clusters
        prev_clusters[[indx]] <- c(prev_clusters[[indx]], prev_clusters[[jndx]])
        clusters[[k]] <- prev_clusters[1:num_clust != jndx]
        # Dendrogram
        dend[[indx]] <- list(dend[[indx]], dend[[jndx]])
        attr(dend[[indx]], 'members') <- attr(dend[[indx]][[1]], 'members') + attr(dend[[indx]][[2]], 'members')
        attr(dend[[indx]], 'height') <- min_dissim
        attr(dend[[indx]], 'class') <- "dendrogram"
        dend <- dend[1:num_clust != jndx]
    }
    # Finish dendrogram
    dend <- dend[[1]]
    
    return(list(clusters, dend))
}