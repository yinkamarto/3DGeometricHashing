#Project Idea
#----------------------------------------------------------------------
#create 2 loops to go through each pair and create hash table
#The first loop goes through each point, the second loop goes
#through every other point after it
#1. Create a hash table by using all possible number of pairs 
# as basis of the input
#2. For rotation, getting the dot product produces the angle of freedom
#3. All points can be displayed using the plot function

#required packages for download
#import the misc3d library - dependency for plot3D
library("misc3d", lib.loc="library/misc3d_0.8-4/")
#import the plot3D library
library("plot3D", lib.loc="library/plot3D_1.1.1/")
#imports the hash library
library("hash", lib.loc="library/3.2/")

#normalizes a vector
mynorm<-function(v) {
  magn=sqrt(sum(v*v))
  result<-v/magn
  return (result)
}

# translates a set of matrix so that a certain pair is the origin
mymatorig<-function(set,basis){
  v<-colMeans(basis)
  newmat<-sweep(set,MARGIN=2,v,FUN="-")
  return(newmat)
}

#Function to compare if 2 matrices have equal values
compareMatrices <- function(m1,m2){
  isEqual <- FALSE
  if(is.matrix(m1) && is.matrix(m2) && (dim(m1)==dim(m2))){
    if(all(as.numeric(m1) == as.numeric(m2))){
      isEqual = TRUE
    }
  }
  return(isEqual)
}
#function that returns a set of bounds(upper, lower, left, right)
getBounds <- function(m){
  topbound <- ceiling(m[,2])
  bottombound <- floor(m[,2])
  leftbound <- floor(m[,1])
  rightbound <- ceiling(m[,1])
  frontbound <- floor(m[,3])
  backbound <- ceiling(m[,3])
  return(c(topbound, bottombound, leftbound, rightbound, frontbound, backbound))
}
#Function to scale and scatter the plot for the bin
scatterPlot <- function(m, scale){
  centroidx <- (m[,3] + m[,4])/2 + scale/10 #i/10 used to scatter the data
  centroidy <- (m[,1] + m[,2])/2 
  centroidz <- (m[,5] + m[,6])/2 
  return(cbind(centroidx, centroidy, centroidz))
}

#Compares two bins to see if they are equal
compareBin<-function(bin1, bin2) {
  if(compareMatrices(bin1, bin2)){
    cat("All match\n")
    print(bin1)
    print(bin2)
    return(TRUE);
  }else{
    return(FALSE);
  }
}

###########################################################
#######Recognition stage ##################################
###########################################################
cat("###################################################\n")
cat("###############Recognition Phase###################\n")
cat("###################################################\n")

#can usefile.choose() to ask user for data
#image points
image <- as.matrix(read.table("object.txt"))
#This plots the data using a polygon3D function
polygon3D(data[c(1,2,3,4,5),1],data[c(1,2,3,4,5),2],data[c(1,2,3,4,5),3], phi = 40, theta = 40,
          col = gg.col(n = 8, alpha = 0.8), NAcol = "white", breaks = NULL,
          colkey = NULL, panel.first = NULL,
          clim = NULL, clab = NULL,
          bty = "b", CI = NULL, surf = NULL,
          xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",
          add = FALSE, plot = TRUE, lwd = 3, border = "black", zlim = c(-0.2,1))

size <- nrow(image)

#important variables
#This specifies what basis to rotate the new image 
i <- 1
j <- 2
k <- 3


#rotated points from file
fileInput <- as.data.frame(read.csv("rotatedPoints.csv", row.names = 1, as.is=TRUE))
hashnames <- as.vector(as.data.frame(read.csv("hashnames.csv")))
h <- hash(keys=hashnames[,1])
allPoints <- matrix(as.matrix(fileInput), nrow=((ncol(fileInput))*nrow(fileInput)), ncol=1)

#convert file input to list of matrices
hashcount <- 1
mholder <- matrix(nrow=5, ncol=3, byrow=TRUE)
newM <- matrix(ncol=3, byrow=TRUE)
#Removes empty first row from data frame
newM <- newM[-1,]
fileList <- list()
contentSize <- length(image)

for(ind in seq(1,nrow(allPoints),15)){
  mholder <- cbind(allPoints[ind:(ind + (size-1))], allPoints[(ind+size):(ind + (size*2)-1)], 
                   allPoints[(ind+(size*2)):(ind + (size*3)-1)])
  fileList[[paste(hashnames[hashcount,])]] <- mholder
  newM <- rbind(newM, mholder)
  hashcount <- hashcount + 1
}

cat("Using ", i, ", ", j , " and ", k, " as basis \n")

#calculate basis to translate all points to this new origin
cat("calculating basis from points\n")
basis <- matrix(c((image[i,1] + image[j,1] + image[k,1])/3, (image[i,2] + image[j,2] + image[k,2])/3,
                  (image[i,3] + image[j,3] + image[k,3])/3), ncol=3, byrow=TRUE)

#a new matrix that holds translated points
cat("Translating points to midpoint\n")
newmat <- mymatorig(image, basis)

#Using points i, j and k
p <- as.vector(data[i,])
q <- as.vector(data[j,])
r <- as.vector(data[k,])

#Create new vector formed by the three points
print("Calculating new vector from selected points")
nx <- mynorm(q-p)
#line of qr
lineqr <- mynorm(r-q)
nz <- mycross(nx, lineqr)
ny <- mycross(nx, nz)

#list to hold new values
M <- matrix( nrow = size, ncol = 3, byrow=TRUE)
#Variable to iterate within matrix
a <- 1
#rotate each points
print("Rotating each points")
newmatSize <- nrow(newmat)

#Using the 3D rotation matrix formula
mBind <- rbind(nx,ny,nz)

#list to hold new values
M <- matrix(nrow = size, ncol = 3, byrow=TRUE)
#Variable to iterate within matrix
a <- 1
#rotate each points
cat("Rotating each points by alpha \n")
newmatSize <- nrow(newmat)

for(ind in 1:newmatSize){
  #Carry out matrix multiplication of the matrix and the selected points
  M[a, ] <- mBind %*% newmat[ind, ]
  #Increment count variable
  a<-a+1
}

mykey <- paste('M',i, j, k, sep = ",")
polygon3D(M[c(1,2,3,4,5),1],M[c(1,2,3,4,5),2],M[c(1,2,3,4,5),3], phi = 90, theta = 30,
          col = gg.col(n = 8, alpha = 0.8), NAcol = "white", breaks = NULL,
          colkey = NULL, panel.first = NULL,
          clim = NULL, clab = NULL,
          bty = "b", CI = NULL, surf = NULL, main = paste("Object after rotation using ", mykey), 
          xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",
          add = FALSE, plot = TRUE, lwd = 3, border = "black", xlim = c(min(M[,1]), max(M[,1])),
          ylim = c(min(M[,2]), max(M[,2])), zlim = c(min(M[,3]), max(M[,3])))

cat("Checking if there is a rotation of ", mykey, " in hash table \n")
mydatakey <- paste('M',i, j, sep = ".")

cat("Checking if hash table contains points ",i, ",", j, " and ", k," as basis \n")
if(mykey %in% names(fileList)){
  cat("Key found in hash\n")
  cat("--------hashed image------------\n")
  #verifying values
  print(fileList[mykey])
  #convert list item to matrix
  hashvalues <- matrix(unlist(fileList[mykey]), ncol=3, byrow=FALSE)
}else{
  cat("Could not find the specified key in hash table \n")
}

#bin should have top, bottom, left, right, front and back bounds
#Use floor and ceil to get the bounds and the average to center the object
hashbin <- list()
myM <- matrix(ncol=6, byrow=TRUE)
for(i in 1:length(fileList)){
  hashname <- names(fileList[i])
  list <- fileList[[hashname]]
  
  mym <- getBounds(list)
  #insert into hashbin
  hashbin[[hashname]] <- matrix(mym, ncol=6, byrow=FALSE)
}

#hashbin data for new image
myI <- matrix(ncol=6, byrow=TRUE)
newbin <-matrix()
for(i in 1:nrow(M)){
  myI <- getBounds(M)
  #insert into hashbin
  newbin <- matrix(myI, ncol=6, byrow=FALSE)
}

#plot hashbin
centroid <- matrix(ncol=4)
centroid <- centroid[-1,]
hashtable <- list()
scale <- -3
for(i in 1:length(hashbin)){
  myM <- matrix(unlist(hashbin[i]), ncol=6, byrow=FALSE)
  centroid <- rbind(centroid, cbind(scatterPlot(myM, scale), names(hashbin[i])))
  hashtable[[names(hashbin[i])]] <- scatterPlot(myM,scale)
  
  scale <- scale + 1
  if(scale > 3)
    scale <- -3
}
#new image
scale <- -3
abin <- matrix(ncol=3)
abin <- abin[-1,]

for(i in 1:nrow(newbin)){
  myI <- matrix(newbin[i,], ncol=6, byrow=FALSE)
  abin <- rbind(abin, scatterPlot(myI, scale))
}

#########################################
###This fixes 0 points in y##############
###Check that code doesn't need in z#####
########################################
#used to fix zero points from plot
for(ind in nrow(centroid):1){
  if(centroid[ind,2] == 0)
    centroid[ind,2] <- 0.5
}
for(ind in nrow(abin):1){
  if(abin[ind,2] == 0)
    abin[ind,2] <- 0.5
}

#######################################################
####Should be rechecked after completion###############
####to see if there is need to customize size##########
#######################################################
# Create a scatter plot
minX <- min(as.numeric(centroid[,1]))
maxX <- max(as.numeric(centroid[,1]))

minY <- min(as.numeric(centroid[,2]))
maxY <- max(as.numeric(centroid[,2]))

minZ <- min(as.numeric(centroid[,3]))
maxZ <- max(as.numeric(centroid[,3]))

#Using absolute helped to fix bug with minus minus in the plot's boundaries
#This generates the hash table in their respective bins
scatter3D(as.numeric(centroid[,1]),as.numeric(centroid[,2]),as.numeric(centroid[,3]), 
          phi = 30, bty = "b2", pch = 1, 
          xlim=c(minX - abs(minX * 0.1), maxX + abs(maxX * 0.1)),
          ylim=c(minY - abs(minY * 0.1), maxY + abs(maxY * 0.1)),
          zlim=c(minZ - abs(minZ * 0.1), maxZ + abs(maxZ * 0.1)),
          cex = 1, ticktype = "detailed")#add = TRUE)

#This adds the new image into it's respective hash bin. It is represented by a square or diamond
scatter3D(as.numeric(abin[,1]),as.numeric(abin[,2]),as.numeric(abin[,3]), 
          phi = 30, bty = "b2", pch = 5, 
          xlim=c(minX - abs(minX * 0.1), maxX + abs(maxX * 0.1)),
          ylim=c(minY - abs(minY * 0.1), maxY + abs(maxY * 0.1)),
          zlim=c(minZ - abs(minZ * 0.1), maxZ + abs(maxZ * 0.1)),
          cex = 1, ticktype = "detailed", add = TRUE)


cat("Picking random image bins to match:\n")
rand <- sample(1:nrow(abin),1)

count <- 0
possibleValues <- matrix(ncol=4,byrow=TRUE)
possibleValues <- possibleValues[-1,]
for(ind in 1:nrow(centroid)){
  if(ceiling(as.numeric(centroid[ind,1]))==ceiling(abin[rand,1]) &&
     ceiling(as.numeric(centroid[ind,2]))==ceiling(abin[rand,2]) &&
     ceiling(as.numeric(centroid[ind,3]))==ceiling(abin[rand,3])){
    possibleValues <- rbind(possibleValues, centroid[ind,])
    count <- count + 1
  }
}
cat(paste("Possible matches:", count, "\n"))
if(nrow(possibleValues)>0){
  cat("Possible basis match:\n")
  cat(possibleValues[,4], ";\n")
}

if(nrow(possibleValues)==1){
  cat("No other competing bin values\n")
  cat("Checking if all other bins with ", possibleValues[,4], " basis match: ")
  actualbin <- hashbin[[possibleValues[,4]]]
  if(compareBin(actualbin, newbin))
    cat("Object rotated on:", possibleValues[,4], " basis")
  else
    cat("Doesn't fully match")
}else if(nrow(possibleValues)==0){
  cat("No match found")
}else{
  cat("Getting actual match:\n")
  count <- 0
  for(ind in 1:nrow(possibleValues))  {
    cat("Checking if all other bins with ", possibleValues[ind,4], " basis match: ")
    actualbin <- hashbin[[possibleValues[ind,4]]]
    
    if(compareMatrices(actualbin, newbin)){
      cat("Match found\n")
      print(actualbin)
      print(newbin)
      cat("Object rotated on:", possibleValues[ind,4], " basis")
      break
    }else{
      cat("Doesn't fully match\n")
    }
  }
}
