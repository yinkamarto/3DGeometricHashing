#Project Idea
#----------------------------------------------------------------------
#create 3 loops to go through each pair and create hash table
#The first loop goes through each point, the second loop goes
#through every other point after it and the third goes through
#each point after the 2nd loop.
#1. Create a hash table by using all possible number of pairs 
# as basis of the input
#2. For rotation, getting the dot product produces the angle of freedom
#3. All points can be displayed using the polygon3D function

#required packages for download
#import the misc3d library - dependency for plot3D
library("misc3d", lib.loc="library/misc3d_0.8-4/")
#import the plot3D library
library("plot3D", lib.loc="library/plot3D_1.1.1/")
#imports the hash library
library("hash", lib.loc="library/hash/3.2/")

#functions for program

#normalizes a vector
mynorm<-function(v) {
  magn=sqrt(sum(v*v))
  result<-v/magn
  return (result)
}

# translates a set of matrix so that a certain pair is the origin
mymatorig<-function(set,basis){
  v<-colMeans(basis)
  newmat<-sweep(set,MARGIN=2,v,FUN="-") #check changing margin to 2
  return(newmat)
}

#get angle between points
getAngle<-function(v1,v2){
  vdot <-v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3]
  v1norm=sqrt(v1[1]*v1[1] + v1[2]*v1[2] + v1[3]*v1[3])
  v2norm=sqrt(v2[1]*v2[1] + v2[2]*v2[2] + v2[3]*v2[3])
  m=vdot/(v1norm * v2norm)
  angle = acos(m) 
  angleInRadians <- angle
  return(angleInRadians)
}

# Cross product of two vectors
mycross<-function(v1,v2){
  xcomp<-v1[2]*v2[3]-v1[3]*v2[2]
  ycomp<-v1[3]*v2[1]-v2[3]*v1[1]
  zcomp<-v1[1]*v2[2]-v1[2]*v2[1]
  vec<-cbind(xcomp,ycomp,zcomp)
  return(vec)
}


#Get degrees from radians
convertToDegrees<-function(angleInRadians){
  angleInDegrees <- angleInRadians * 180.0/3.1415
  return(angleInDegrees)
}

#Using static data to test
data <- matrix(c(1, 3.85, 0.4, 1.37, 6.9, 0.8,  4.74, 5.24, 0.4, 5.43, 0.87, 0.2, 2.56, 0.37, 0.2), 
               ncol=3, byrow=TRUE)
print(data)
#The size is the number of rows in the data. This is used to know the number of possible pairs
size <- nrow(data)

#This plots the data using a polygon3D function
polygon3D(data[c(1,2,3,4,5),1],data[c(1,2,3,4,5),2],data[c(1,2,3,4,5),3], phi = 40, theta = 40,
          col = gg.col(n = 8, alpha = 0.8), NAcol = "white", breaks = NULL,
          colkey = NULL, panel.first = NULL,
          clim = NULL, clab = NULL,
          bty = "b", CI = NULL, surf = NULL,
          xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",
          add = FALSE, plot = TRUE, lwd = 3, border = "black", zlim = c(-0.2,1))

#Rotates the perspective view, not the actual object by the specified theta
#This is for visualization purpose. Can be commented out
#for (angle in seq(0, 360, by = 45)){
#  plotdev(theta = angle, main=paste("Rotation by angle", angle))
#}

#Rotates the perspective view, not the actual object by the specified phi
#This is for visualization purpose. Can be commented out
#for (angle in seq(0, 360, by = 45)){
#  plotdev(phi = angle, theta=40)
#}

#Hashing
hashtable <- list()
h <- hash()

count <- 1
colorIndex <- 1
for(i in 1:(size-2)) {
  for(j in (i+1):(size-1)){
    for(k in (j+1):(size)){
      cat("######################################################\n")
      cat("######### Using points ", i, ", ", j, " and ", k, " as basis######\n")
      cat("######################################################\n")
      
      #midpoint for a 3 dimensional object over 3 points will be the midpoint of 2 points of that object 
      #followed by finding the equation of the line frm that point to the opposite corner
      
      #Using points i, j and k
      p <- as.vector(data[i,])
      q <- as.vector(data[j,])
      r <- as.vector(data[k,])
      
      #midpoint of p, q and r
      basis <- matrix(c((p[1] + q[1] + r[1])/3.0, (p[2] + q[2] + r[2])/3.0, (p[3] + q[3] + r[3])/3.0),ncol=3, byrow=TRUE)
      
      cat("midpoint of points ", i, ", ", j, " and ", k, ": \n")
      print(basis)
      
      ######translation: making the midpoint the origin
      #a new matrix that holds translated points
      print("Translating points to midpoint")
      newmat <- mymatorig(data, basis)
      print("Translated points")
      print(newmat)
      
      #object after translation
      polygon3D(newmat[c(1,2,3,4,5),1],newmat[c(1,2,3,4,5),2],newmat[c(1,2,3,4,5),3], phi = 0, theta = 0,
                col = gg.col(n = 8, alpha = 0.8), NAcol = "white", breaks = NULL,
                colkey = NULL, panel.first = NULL,
                clim = NULL, clab = NULL,
                bty = "b", CI = NULL, surf = NULL, main = "Object after translation", 
                xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",
                add = FALSE, plot = TRUE, lwd = 3, border = "black", xlim = c(min(newmat[,1]), max(newmat[,1])),
                ylim = c(min(newmat[,2]), max(newmat[,2])), zlim = c(min(newmat[,3]), max(newmat[,3])))
      
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
      
      for(ind in 1:newmatSize){
        #Carry out matrix multiplication of the matrix and the selected points
        M[a, ] <- mBind %*% newmat[ind, ]
        #Increment count variable
        a<-a+1
      }
      #Add the matrix to the hashtable
      print("Adding new set to hash table")
      mykey <- paste('M',i, j, k, sep = ",")
      h[mykey] <- count
      hashtable[[mykey]] <- M
      print("Plotting graph of previous points and rotated points")
      #object after rotation
      polygon3D(M[c(1,2,3,4,5),1],M[c(1,2,3,4,5),2],M[c(1,2,3,4,5),3], phi = 90, theta = 30,
                col = gg.col(n = 8, alpha = 0.8), NAcol = "white", breaks = NULL,
                colkey = NULL, panel.first = NULL,
                clim = NULL, clab = NULL,
                bty = "b", CI = NULL, surf = NULL, main = paste("Object after rotation using ", mykey), 
                xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",
                add = FALSE, plot = TRUE, lwd = 3, border = "black", xlim = c(min(M[,1]), max(M[,1])),
                ylim = c(min(M[,2]), max(M[,2])), zlim = c(min(M[,3]), max(M[,3])))
      
    }
    print("Printing newly translated and rotated points")
    print(M)
  }
}

#checking value of hash
print("printing Hash table")
print(hashtable)

print("Writing data to file")
write.csv(hashtable,"rotatedPoints.csv", row.names = TRUE)
print("writing hash names to file")
write.csv(keys(h),"hashnames.csv", row.names = FALSE)
print("Complete")

