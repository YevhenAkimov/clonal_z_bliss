# clonal_z_bliss

clonal_z_bliss example usage

```r
source("https://raw.githubusercontent.com/YevhenAkimov/clonal_z_bliss/main/clonal_z_bliss.R")
gr_rates1= c(1,7,1,4,3,1,2,1,2,5,5)/10
gr_rates2= c(5,1,4,1,1,3,5,5,5,1,1)/10
gr_control=c(9,7,5,4,8,5,4,6,8,5,7)/10
sizes= rep(1,length(gr_rates1))
sizes=sizes/sum(sizes)
clonal_z_bliss(gr_rates1,gr_rates2,gr_control,sizes,t=5,n_perm=100)
```
