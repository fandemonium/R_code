## in R:
 df<-read.delim("../../../rpf_trees/rpf_w_refs/node_tips/ref_dropped_clade_trait_table.txt", sep="\t", header=F)
library(reshape2)
 df.wide<-dcast(df, oid ~ clades, value.var="counts")
> head(df.wide)
        oid 1805 2164 2244 2710 2922 2993 3111 3228 3441
1 637000001    8   NA   NA   NA   NA   NA   NA   NA   NA
2 637000012   NA   NA   NA   NA   NA   NA   NA    1   NA
3 637000014   NA   NA   NA   NA   NA   NA   NA   NA    1
4 637000015   NA   NA   NA   NA   NA   NA   NA   NA    1
5 637000020   NA   NA   NA   NA   NA   NA   NA   NA    1
6 637000060   NA   NA    1   NA   NA   NA   NA   NA   NA
> 
 df.wide<-dcast(df, oid ~ clades, value.var="counts")
> head(df.wide)
        oid 1805 2164 2244 2710 2922 2993 3111 3228 3441
1 637000001    8   NA   NA   NA   NA   NA   NA   NA   NA
2 637000012   NA   NA   NA   NA   NA   NA   NA    1   NA
3 637000014   NA   NA   NA   NA   NA   NA   NA   NA    1
4 637000015   NA   NA   NA   NA   NA   NA   NA   NA    1
5 637000020   NA   NA   NA   NA   NA   NA   NA   NA    1
6 637000060   NA   NA    1   NA   NA   NA   NA   NA   NA
> 
> head(ref_clades.trait.table)
  data.1037.16s.ref.dropped.traits.oid X1805 X2164 X2244 X2710 X2922 X2993
1                           2265129005     0     0     0     0     0     0
2                           2501939606     0     0     0     0     0     0
3                           2502545039     0     0     0     0     0     0
4                           2502790015     0     0     0     0     0     0
5                           2503283023     0     0     0     0     0     0
6                           2503508007     0     0     0     0     0     0
  X3111 X3228 X3441
1     0     0     0
2     0     0     0
3     0     0     1
4     0     0     0
5     0     0     0
6     0     0     0
