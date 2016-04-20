##  Copyright (C) 2009 Andrés Baselga   <andres.baselga@usc.es>
##                    
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

# Three functions are included:

# beta.sor computes a distance matrix using Sorensen pair-wise dissimilarity 
# beta.sim computes a distance matrix using Simpson pair-wise dissimilarity 
# beta.nes computes a distance matrix using Nestedness-resultant pair-wise 
# dissimilarity 


beta.sor<-function(x){
## x must be a data frame
matr<-as.matrix(x);
nr<-dim(matr)[1];
result<-matrix(nrow=nr,ncol=nr);
rownames(result)<-rownames(matr);
colnames(result)<-rownames(matr);
for(i in 1:nr) {
    for(j in i:nr) {
bij<-sum(matr[i,]&(!matr[j,]));
cij<-sum(matr[j,]&(!matr[i,]));
aij<-sum(matr[i,]&(matr[j,]));
Sorensen.ij<-(bij+cij)/((2*aij)+bij+cij);
result[i,j]<-Sorensen.ij;
result[j,i]<-Sorensen.ij;
}
}
d<-as.dist(result);
d
}

beta.sim<-function(x){
## x must be a data frame
matr<-as.matrix(x);
nr<-dim(matr)[1];
result<-matrix(nrow=nr,ncol=nr);
rownames(result)<-rownames(matr);
colnames(result)<-rownames(matr);
for(i in 1:nr) {
    for(j in i:nr) {
bij<-sum(matr[i,]&(!matr[j,]));
cij<-sum(matr[j,]&(!matr[i,]));
aij<-sum(matr[i,]&(matr[j,]));
Simpson.ij<-min(bij,cij)/(min(bij,cij)+aij);
result[i,j]<-Simpson.ij;
result[j,i]<-Simpson.ij;
}
}
d<-as.dist(result);
d
}

beta.nes<-function(x){
## x must be a data frame
matr<-as.matrix(x);
nr<-dim(matr)[1];
result<-matrix(nrow=nr,ncol=nr);
rownames(result)<-rownames(matr);
colnames(result)<-rownames(matr);
for(i in 1:nr) {
    for(j in i:nr) {
bij<-sum(matr[i,]&(!matr[j,]));
cij<-sum(matr[j,]&(!matr[i,]));
aij<-sum(matr[i,]&(matr[j,]));
nestedness.ij<-((max(bij,cij)-min(bij,cij))/((2*aij)+bij+cij))*(aij/(min(bij,cij)+aij));
result[i,j]<-nestedness.ij;
result[j,i]<-nestedness.ij;
}
}
d<-as.dist(result);
d
}