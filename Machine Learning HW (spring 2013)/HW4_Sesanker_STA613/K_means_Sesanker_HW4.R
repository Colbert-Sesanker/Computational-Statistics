# Colbert Sesanker: HW 4, k-means clustering, source file for demo, Feb 2013
pts_mat <- as.matrix(read.table("points_hw4.txt", sep=" "))
points  <- split(pts_mat, c(row(pts_mat))) # split matrix into a mapable list of points
means   <- list(c(3,3), c(4, .7))# list(c(2,3), c(1.5, 3.5)) is an example where it converges to local minimum
classes <- 1:length(means)
tol     <- .00001 # l2 distance beetween updated centroids as stopping criteria

l2 <- function(x,y) sqrt(sum((x-y)^2))  # l2 distance between x and y
# assigns a point to closest mean (returns the index (class) of the mean closest to the point)
assign_cluster <- function(means, point)  which.min(Map(partial(l2, point), means)) 

vec_sum <- function(x) Reduce("+", x) # vector sum of a list of vectors
update  <- function(means) {
    assigns  <- Map(partial(assign_cluster, means), points) # list of mean assignments for each point
    in_class <- function(class, assign, point) if (class==assign) point else 0  # returns point if point is in "class", 0 otherwise
    new_mean <- function(class) {# returns new mean using all points assigned to old mean
        class_points <- Map(partial(in_class, class), assigns, points) # zeros all points of points not in specified "class"
        return(vec_sum(class_points)/sum(assigns == class)) # sums points in given class and divides by the number of points 
    }
    return(list("means"=Map(new_mean, classes), "assignments"=assigns))
}
# defined recursively, breaks on tol or count
k_means <- function(means, count=1) {    
    step   <- update(means)
    change <- l2(vec_sum(step$means), vec_sum(means)) # l2 distance between sum of previous centriods to current centroids
    if (count >= 1000 || change <= tol) return(c(step,count=count)) # count added to the list to see when it converges
    k_means(step$means, count=count+1)
}
# Graphing
km <- k_means(means)
x_min <- min(pts_mat[,1]); x_max <- max(pts_mat[,1]);  y_min <- min(pts_mat[,2]);  y_max <- max(pts_mat[,2]); 
#png(file="K-means_clusters", width=512, height=512,res=72)
plot(1, type="n", xlim=c(x_min,x_max), ylim=c(y_min,y_max), xlab="x", ylab="y",
     main=paste("K-means with initial means", means[1], means[2]))
plot_tuple <- function(x, color) points(x[1],x[2], col=color)
plot_tuple(km$means[[1]], "red"); plot_tuple(km$means[[2]], "red"); 
for (p in 1:length(points)) { # color points to respective clusters
  if (km$assignments[[p]]==2) {
    plot_tuple(points[[p]], "blue")
  }
  else
    plot_tuple(points[[p]], "green")
}
#dev.off()

# Tends to converge after about eight iterations

# For partial application of functions
# partial(f, arg1) returns a function with arg1 bound to the first argument of function f
partial <- function(f, ...) {
  capture <- list(...)
  function(...) {
    do.call(f, c(capture, list(...)))
  }
}
