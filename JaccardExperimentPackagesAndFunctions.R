#Functions to support the ClusterExperiment_v2 Rmd file.

## ---- packagesAndFunctions ---------

set.seed(816) #Always set the random seed for replicatability of random processes.

library(dplyr)
library(cluster)
library(ggplot2)

########################################################################
##To print numbers,x, with a consistent number of digits,d.
nicePrint<-function(x,d=3) {
  d<-paste("%1.",d,"f",sep="")
  return(sprintf(d, x))
}

########################################################################
##The basic jaccard similarity function, my own simple version
jaccard<-function(S1,S2) {
  x<-length(intersect(S1,S2))/length(union(S1,S2))
  return(x)
}

########################################################################
##The next two functions are for our data experiments.
##They generate event spaces and sequences, 
##which both consist of letters of the alphabet in our data experiments.

#Define the function that creates new event spaces
#The event space has a size (length) and an offset from the first letter, A.
genEventSpace<-function(offset=0,size=10) {
  start=1+offset
  end = start+size-1
  eventSpace<-toupper(letters[start:end])
  return(eventSpace)
}

#Define the function that creates sequences
#Either select a sequence of length that varies randomly between low and high
#or if Low==High, a uniform length.
genOneSequence<-function(eventSpace,Low,High) {
  if(Low==High) {
    s<-sample(eventSpace,Low)
  }
  else{
    s<-sample(eventSpace,sample(Low:High,1))  
  }
  return(s)
}

########################################################################
##This is the basic engine for generating our list of simulated event sequences.
##It samples sequences from N event spaces, also generated from the function. 
##For now it assumes a fixed offset and size, but these could become vectors.

sequenceGenerator<-function(offset=4,size=10,Low=2,High=8,trialsPerEventSpace=200,eventSpaces=2) {
  
  eventSequences<-list()
  
  for(i in 1:eventSpaces) {
    E<-genEventSpace( ( (i-1)*offset ) ,size)
    eventSequences[[i]]<-lapply(1:trialsPerEventSpace, function(i) genOneSequence(E,Low,High))
  }
  eventSequences<-unlist(eventSequences,recursive=FALSE)
  return(eventSequences)
}

######################################################################## 
##Sets up the pairwise operations on the sequences list
##Here we use the jaccard function to calculate the similarity between each pair of sequences
##Then turn that into a distance by subtracting the similarity from 1.
pairwiseDist<-function(SEQS) {
  
  #get the length of the list
  LEN<-length(SEQS) 
  #create a sequence similarity matrix of size NxN, initially populated by 0s
  SIM<-matrix(data=0,nrow=LEN,ncol=LEN)
  #set the diagonals to 1
  diag(SIM)<-1
  #create logical matrices so that the vector of similarities can be put into matrix form
  simLower<-lower.tri(SIM)
  simUpper<-upper.tri(SIM)
  
  #create a vector of all pairwise combos of sequences, comput the jaccard similarity for each pair.
  combos<-combn(LEN,2)
  #this function goes across the columns of combos, reads two rows of the SEQS matrix and computes similarity.
  s<-sapply(1:length(combos[1,]), function(i) jaccard(SEQS[[combos[1,i]]],SEQS[[combos[2,i]]]) )
  #put the vector of similarities into the SIM matrix
  SIM[simLower]<-s  
  SIM[simUpper]<-s
  #Return a distance matrix
  DIST<-1-SIM
  return(DIST)
}

########################################################################
###This similarity function compares event sequences to event spaces.
###It is called from the function clusterIteration, we repeats clustering
###multiple times in order to extract event spaces.

eventJaccard<-function(SEQS,lim=0,boost=0,e1,e2,e3=NULL,e4=NULL) {
  #Plays a role similar to pairwiseDist.  This is used in our 3 stage process, when we have created clusters
  #using the longest sequences and recovered the event spaces.  
  #Then, we just walk through the original list and compare each string to the recovered
  #events spaces, e1 and e2.
  LEN<-length(SEQS)
  s<-NULL
  
  #If we only have 2 event spaces, calcuate the similarity of an event sequence to an event space this way.
  if(is.null(e3)){

    
  #If the event sequence is more similar to event space e1, return the value 1.  If it's more similar to e2, return 2.  Otherwise, 
  #return 0, indicating it's euqually similar to e1 and e2.
  s<-sapply(1:LEN, function(i) 
    ifelse( ( jaccard(SEQS[[i]],e1)-jaccard(SEQS[[i]],e2) )/(boost+jaccard(SEQS[[i]],e1)) > (lim),1,
            ifelse( ( jaccard(SEQS[[i]],e2)-jaccard(SEQS[[i]],e1) )/(boost+jaccard(SEQS[[i]],e2)) > (lim),2,0))
  )
} 

#If we have 4 event spaces, calculate similarity this way.  
else{
  for (i in 1:LEN) {
    x<-NULL
    y<-NULL
    x<-c(-1,-1,-1,-1)
    x[1]<-jaccard(SEQS[[i]],e1)
    x[2]<-jaccard(SEQS[[i]],e2)
    x[3]<-jaccard(SEQS[[i]],e3)
    x[4]<-jaccard(SEQS[[i]],e4)
    m1<-which.max(x)
    y<-x[-m1]
    m2<-which.max(y)
    crit<-(x[m1]-y[m2])/x[m1]
    s<-append(s,ifelse( crit>0,m1,0 ))
  }
  
}
  
  return(s)
} 

########################################################################
###Repeats clustering several times and extracts event spaces.
### o = offset, h and l are high and low values for J(Si,Ei).

clusterIteration<-function(o,h,l,clusters=2,reps=20,cutoff=0.8,divisive=TRUE) {

#generates our initial set of event sequences for this experiment  
L<-sequenceGenerator(offset=o,High=h,Low=l,eventSpaces = clusters)
#get the length of each sequence
q<-sapply(1:length(L), function(i) length(L[[i]]))
#an integer that is the 80th percentile cutoff value
x<-as.numeric(quantile(q,cutoff))
#Generate a set of data that are the longest event sequences
Q<-L[q>=x]
#Run the function kCluster, defined below, which will extract the event sequences.
kClusters(clusters=clusters,Q,L,divisive)
}

#The 'divisive' parameter indicates whether we're using diana or pam clustering.
kClusters<-function(clusters=2,Q,L,divisive) {

events1<-NULL
events2<-NULL
events3<-NULL
events4<-NULL

for(z in 1:reps) {
  ind<-sort(sample(1:length(Q),floor(0.8*length(Q))))
  H<-Q[ind]
  
#We'll go ahead and run both types of clustering, but only return results based on the 'divisive' parameter.
  dv<-pairwiseDist(H) %>%
      diana(diss=TRUE)

  pm2<-pairwiseDist(H) %>% pam(diss=TRUE,k=clusters)

  
  cut<-cutree(as.hclust(dv),k=clusters)
  cutp<-pm2$clustering

if (divisive) {
  events1<-append(unique(unlist(H[cut==1])),events1)
  events2<-append(unique(unlist(H[cut==2])),events2)
  if(clusters==4) {
    events3<-append(unique(unlist(H[cut==3])),events3)
    events4<-append(unique(unlist(H[cut==4])),events4)
  }
}
else {
  events1<-append(unique(unlist(H[cutp==1])),events1)
  events2<-append(unique(unlist(H[cutp==2])),events2)
  if(clusters==4) {
    events3<-append(unique(unlist(H[cutp==3])),events3)
    events4<-append(unique(unlist(H[cutp==4])),events4)
  }
}
  
  
}

zz<-tbl_df(as.factor(events1))

e1<-names(table(zz)[table(zz) > 0.9*reps])

zz<-tbl_df(as.factor(events2))

e2<-names(table(zz)[table(zz) > 0.9*reps])

if(clusters==4) {
  zz<-tbl_df(as.factor(events3))
  e3<-names(table(zz)[table(zz) > 0.9*reps])
  
  zz<-tbl_df(as.factor(events4))
  e4<-names(table(zz)[table(zz) > 0.9*reps])
}

if(clusters==2) {
S<-eventJaccard(L,e1,e2)
J<-list(S,sort(e1),sort(e2),divisive)
names(J)<-c("classification","event_space1","event_space2")
}

if(clusters==4) {
  S<-eventJaccard(L,e1,e2,e3,e4)
  J<-list(S,sort(e1),sort(e2),sort(e3),sort(e4))
  names(J)<-c("classification","event_space1","event_space2","event_space3","event_space4")
}

return(J)

}
