#a) K-means clustering 
# some words on this partitioning algorithm:
 # this method consists of partinioning a dataset D of n objects into a set of k clusters such that #observations in the same group/cluster present high similarity and they are very dissimilar            #from the objects partitioned in other groups. #the weakness of this method:
 # - it may be applicable only when the mean is defined --> problem with categorical data -->              # a possible solution is given by using k-modes algorithm: replacing means of clusters with k modes. # - unable to handle noisy data and outliers  - not suitable to discover clusters with non-convex                    # shapes  the algorithm:  considering a range of k=1,...,10 groups 
WGSS <- rep(0,10) #we now want to record the WGSS for each clustering solution where WGSS is  #the sum of squares for each k=1,2,...,10 clustering solution 
#note that the aim of K-means is to minimize the within cluster sum of squares. we are going to #consider the standardized data set prior created (X.std) n = nrow(X.std) WGSS[1]=(n-1)*sum(apply(X.std,2,var)) # we used apply because we want to perform the variance on the #numerical matrix by column(2) 
for(k in 2:10) {
 WGSS[k] <- sum(kmeans(X.std, centers = k)$withinss) } #note that the more clusters we fit and the #smaller is the SS, however I cannot base my choice only looking for the smallest SS.  Thereby, we #may plot the number of clusters against SS and look at the elbow in the curve.
 plot(1:10, WGSS, type="b", xlab="k", ylab="Within group sum of squares") #we decide to take k=3 #after observing the elbow of the graph between k=3 and k=4. # the reason why we chose k=3 is #based on the elbow method. This method looks at the percentage of variance explained as a #function of the number of clusters. One should choose a number of clusters so that adding #another cluster does not give much better #information on the data. Note that the first #cluster/group explain a lot of variance (information). #computing a cross-tabulation of clustering #solution
 K <- 3 cl <- kmeans(X.std, center=K) table(cl$cluster)
 #number of chemical compositions in each cluster # cluster 1 =10, cluster2 = 14, cluster3=21                  #total chemical compositions grouped 45 --> missing 0 #since there are no clusters with very few #objects, we may be satisfied of the solution. #drawing a scatter plot of data objects grouped in #three clusters 
plot(X.std, col= cl$cluster) points(cl$centers, col=1:K, pch=8)
#b) hierarchical clustering 
# some words on this clustering algorithm:  it consists of grouping a dataset into a hierarchy or a tree # of clusters. Unlike partitioning algorithm, this algorithm does not require the number of clusters in #advance,  but it needs terminal conditions. We may define the hierarchical clustering as a bottom-# up strategy and #its results are plotted though the dendogram.
 #the algorithm: 
cl.average = hclust(dist(X.std), method="average") 
#automatically Euclidan Distance has been used to calculate the dissimilarity matrix. #method="average" is the "average linkage" through which the distance between two clusters is #measured as the mean distance between an observation in one cluster and an observation in the #other cluster. Unlike other linkage methods (simple and complete linkage), the average linkage #method uses a more central measure of location. 
plot(cl.average) #dendogram #we decide to take k=3
 # from right to left:
 #group1 ---> 13-|16 no.observations:21
 #group2 ---> 24-|22 no.observations:14 
#group3 ---> 44-|36 no.observations:10 
#observations no: 16,31,18 may look to be outliers. if we have a look at the data (X.std), we observe #that the observations no.16,18,31 present values very far different from the others for Na2O, CaO #and MgO respectively. 
#cut dendogram 
hcl= cutree(hclust(dist(X.std)), 3) #cut the dendogram by setting k=3 table(hcl) 
#cluster1= 21, cluster2= 14, cluster3=10
#c) compare two methods: The Rand Index tab=table(hcl,cl$cluster) # in order to compare two #clustering solutions, we refer to the "classAgreement" #function
library(e1071)
 #call the library where the function is located 
values=classAgreement(tab)
 values$rand #1=strong agreement. We refer to the rand index to establish the agreement between #the two solutions. It is a number with minimum value 0 (little agreement) and maximum value 1 #(strong agreement). #Rand Index= A/(A+D) where A = Agreement and D= Disagreement. In some #cases, It is better to use the #Adjusted Rand Index to take into account the agreement due to the #chance. 
library(mclust) adjustedRandIndex(hcl,cl$cluster) #1
#d) clustering solutions and the kiln variable
 table(hcl,potterydata[,10]) #cross-tabulation of hierarchical clustering solution and the location       #) at which the pottery was found 
table(cl$cluster,potterydata[,10]) #cross-tabulation of k-means solution and the location (the kiln) #we may say that in both cases, the cluster 1 groups chemicals spotted in location 1, cluster 2 groups #observations found in location:2,3 and finally the cluster 3 summarises 10 observations belonging #to the location 4,5. #we may be quite satisfied of our solutions. 

#Question 3
 #3A) Does the new population member suffer from diabetes? Explain how you arrived at your #answer  load the dataset Pima 
pima=read.csv("Pima.csv")
 summary(pima) #The variables present different scales but we do not need to be worried about that #because  the QDA results on standardized and non-standardized features are going to be exactly #the same Eigenvalues, standardized coefficients, structure correlations, discriminant scores - #everything will be the same. 
#Only eigenvectors will differ. The reason why there is no effect of standardization on the main #results in QDA is that QDA decomposes ratio of Between-to-Within covariances, and not the #covariance itself #having its magnitude (as PCA does). Thereby, applying directly Quadratic #discriminant analysis. 
library(MASS) #we need to call the packaged MASS where the function QDA is contained qda.res=qda(type~npreg+glu+bp+skin+bmi+ped+age,data=pima)  #groups~x1+x2+....+x7 attributes(qda.res) # lists of the attributes that I can call 
qda.res$prior # percentage of women with diabetes (Yes) and without (No) ~ Prior probabilities of #groups basically is the same as: summary(pima$type)/nrow(pima)
qda.res$means # in other words we get the estimates of the group prior probabilities ~ group means 
# so now we need to estimate the covariance matrix.  This is simply the weighted average of the #group specific covariance matrix estimates. 
N <- nrow(pima) #number of observations 
pima.yes=subset(pima,type=="Yes")  #subset of women with diabetes
pima.no=subset(pima,type=="No")  #subset of women without diabetes
#computing the matrix covariance for each subset excluding the type variable cov.pima.yes=cov(pima.yes[,1:7]) 
cov.pima.no=cov(pima.no[,1:7]) 
#we want to write the quadratic discriminant function for any observation x for each group g. 
qdf = function(x, prior, mu, covar)  {  x = matrix(as.numeric(x), ncol=1)  log(prior) -(0.5*log(det(covar)))-0.5*t(x-mu)%*%solve(covar)%*%(x-mu) }
#Note: Unlike the Linear discriminant Function, we do not have the assumption that the covariance #of each class is identical ? we have a different discriminant function computed for each class.
#  calculating the value of quadratic discriminant function for the new patient: 
id = c(7,187,50,33,33.9,0.826,30)
 dfyes= qdf(id, qda.res$prior[2], qda.res$mean[2,], cov.pima.yes ) 
dfno=qdf(id, qda.res$prior[1], qda.res$mean[1,], cov.pima.no )
 #so we get two different values referring to the type: Yes or No. 
# the largest one indicates the belonging to the group/ Yes or No in that our initial aim was to #maximize the quadratic discriminant function max(dfyes,dfno)
 # -17.85404 ~ Yes # YES. The new population member suffers of diabetes
#3B) Using the model ?tted in part (a), calculate the misclassi?cation rate of your quadratic #discriminant  analysis model.
 #adding the test data to the training data
 pima[nrow(pima)+1,]=c(2,88,58,26,28.4,0.766,22,"No")
 pima[nrow(pima)+1,] =c(9,170,74,31,44,0.403,43,"Yes") 
pima[nrow(pima)+1,] =c(10,101,76,48,32.9,0.171,63,"No") 
pima[nrow(pima)+1,] =c(5,121,72,23,26.2,0.245,30,"No")
 pima[nrow(pima)+1,] =c(1,93,70,31,30.4,0.315,23,"No")
pimatest=pima[528:532,] #considering only the test data
 res=rep(0,5) 
dfyes=rep(0,5)
 dfno=rep(0,5)
 for (i in 1:5) { dfyes[i]= qdf(pimatest[i,1:7], qda.res$prior[2], qda.res$mean[2,], cov.pima.yes ) dfno[i]= qdf(pimatest[i,1:7], qda.res$prior[1], qda.res$mean[1,], cov.pima.no ) res[i]=max(dfyes[i],dfno[i]) } 
class=rep(0,5)
 for (i in 1:5) { if (res[i]==dfyes[i]){  class[i]="yes" } else {  class[i]="no" }}
tab=table(class,pimatest[,8])
 mis.class.rate=tab[2,1]/nrow(pimatest) #0.2
#3c Why might a user choose to use linear discriminant analysis over quadratic discriminant #analysis? 
library(lattice)
#In order to display the two covariance matrices as a surface, we refer to #a level plot using coloured #regions to identify regions of different heights 
levelplot(cov.pima.yes)
 levelplot(cov.pima.no) # As we can see from the levelplot, the two classes present similar #covariance --> the assumption of equal covariance matrices might be valid=LDA

