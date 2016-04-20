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

# beta.SOR computes the Sorensen-based multiple-site dissimilarity 
# beta.SIM computes the Simpson-based multiple-site dissimilarity 
# beta.NES computes the Nestedness-based multiple-site dissimilarity 

# Some previous functions are needed

maxbibj<-function(matrix) {
  ## nr is the number of rows
  
  nr<-dim(matrix)[1];

  ## This variable contains the sum of 'max(bi,bj)' on row pairs.

  sum.maxbibj<-0.;

  ## We 'loop' on every different row pairs.
  
  for(i in 1:nr-1) {
    for(j in (i+1):nr) {
      
      ## 'bi' and 'bj' are respectively the number of species appearing
      ## only in the ith and the jth site (row).

      bi<-sum(matrix[i,]&(!matrix[j,]));
      bj<-sum(matrix[j,]&(!matrix[i,]));
      
      ## We sum up 'max(bi,bj)' on every different pair of rows.

      sum.maxbibj<-sum.maxbibj+max(bi,bj);
    }
  }

  ## We return the sum

  sum.maxbibj;
}


minbibj<-function(matrix) {
  ## nr is the number of rows
  
  nr<-dim(matrix)[1];

  ## This variable contains the sum of 'min(bi,bj)' on row pairs.

  sum.minbibj<-0.;

  ## We 'loop' on every different row pairs.
  
  for(i in 1:nr-1) {
    for(j in (i+1):nr) {
      
      ## 'bi' and 'bj' are respectively the number of species appearing
      ## only in the ith and the jth site (row).

      bi<-sum(matrix[i,]&(!matrix[j,]));
      bj<-sum(matrix[j,]&(!matrix[i,]));
      
      ## We sum up 'min(bi,bj)' on every different pair of rows.

      sum.minbibj<-sum.minbibj+min(bi,bj);
    }
  }

  ## We return the sum

  sum.minbibj;
}


zeroColumn<-function(matrix) {
  sum<-0;
  nc<-dim(matrix)[2];
  for(i in 1:nc) if(!sum(matrix[,i])) sum<-sum+1;
  if(sum!=0)
    warning(sum," missing species. This fact does not affect the result since total number of species is calculated excluding these columns.",call.=FALSE,immediate.=TRUE);
  sum
}

beta.SOR<-function(x) {

  ## x must be a data frame
  
  matrix<-as.matrix(x);
  
  ## Si is the number of species present in the ith site. We sum all
  ## these values.

  sumSi<-sum(matrix);

  ## St is the total number of species in all sites, thus the number of columns, excepting if a given
  ## species is not present at all (in any site), so we must not count this column
  ## out.
  
  St<-ncol(matrix)-zeroColumn(matrix);

  ## 'a' is the number of species shared for at least two sites

  a<-sumSi-St;
  index<-(minbibj(matrix)+maxbibj(matrix))/(minbibj(matrix)+maxbibj(matrix)+(2*a));
  if(any(matrix>1))
    warning("The table contains values >1: the result may be meaningless.",call.=FALSE,immediate. = TRUE);
  if(any(matrix<0))
    warning("The table contains values <0: the result may be meaningless.",call.=FALSE,immediate. = TRUE);  
  index;

}

beta.SIM<-function(x) {

  ## x must be a data frame
  
  matrix<-as.matrix(x);
  
  ## Si is the number of species present in the ith site. We sum all
  ## these values.

  sumSi<-sum(matrix);

  ## St is the total number of species in all sites, thus the number of columns, excepting if a given
  ## species is not present at all (in any site), so we must not count this column
  ## out.
  
  St<-ncol(matrix)-zeroColumn(matrix);

  ## 'a' is the number of species shared for at least two sites

  a<-sumSi-St;
  index<-minbibj(matrix)/(minbibj(matrix)+a);
  if(any(matrix>1))
    warning("The table contains values >1: the result may be meaningless.",call.=FALSE,immediate. = TRUE);
  if(any(matrix<0))
    warning("The table contains values <0: the result may be meaningless.",call.=FALSE,immediate. = TRUE);  
  index;

}

beta.NES<-function(x) {

  ## x must be a data frame
  
  matrix<-as.matrix(x);
  
  ## Si is the number of species present in the ith site. We sum all
  ## these values.

  sumSi<-sum(matrix);

  ## St is the total number of species in all sites, thus the number of columns, excepting if a given
  ## species is not present at all (in any site), so we must not count this column
  ## out.
  
  St<-ncol(matrix)-zeroColumn(matrix);

  ## 'a' is the number of species shared for at least two sites

  a<-sumSi-St;
  index<-(a/(minbibj(matrix)+a))*((maxbibj(matrix)-minbibj(matrix))/((2*a)+maxbibj(matrix)+minbibj(matrix)));
  if(any(matrix>1))
    warning("The table contains values >1: the result may be meaningless.",call.=FALSE,immediate. = TRUE);
  if(any(matrix<0))
    warning("The table contains values <0: the result may be meaningless.",call.=FALSE,immediate. = TRUE);  
  index;

}
