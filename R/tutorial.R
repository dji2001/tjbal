devtools::install_github('xuyiqing/panelView')
library(tjbal)
#> ## Syntax has changed since v.0.4.0.
#> ## See http://bit.ly/tjbal4r for more info.
#> ## Comments and suggestions -> yiqingxu@stanford.edu.
data(tjbal)
ls()
#> [1] "germany" "npc"
head(npc)
library(panelView)
#> ## See bit.ly/panelview4r for more info.
#> ## Report bugs -> yiqingxu@stanford.edu.
panelview(roa ~ treat + so_portion + rev2007, data = npc, show.id = c(1:50), 
          index = c("gvkey","fyear"), xlab = "Year", ylab = "Firm ID",
          axis.lab.gap = c(0,1), by.timing = TRUE, 
          display.all = TRUE, axis.lab = "time")
with(npc, table(fyear, treat))
out.did <- tjbal(roa ~ treat, data = npc,   
                 index = c("gvkey","fyear"), Y.match.npre = 0,estimator="meanfirst",
                 demean = TRUE, vce = "boot", nsims = 200)
print(out.did)
plot(out.did, ylab = "ATT", ylim = c(-0.06, 0.06))
out.mbal <- tjbal(roa ~ treat + so_portion + rev2007, data = npc,
                  index = c("gvkey","fyear"), demean = FALSE, estimator = "meanfirst",
                  vce = "jackknife")
print(out.mbal)
out.mbal <- tjbal(data = npc, Y = "roa", D = "treat", X = c("so_portion","rev2007"), 
                  index = c("gvkey","fyear"), demean = FALSE, estimator = "mean",
                  vce = "jackknife", nsims = 200)
plot(out.mbal, ylim = c(-0.04, 0.06))
begin.time<-Sys.time() 
out.kbal <- tjbal(roa ~ treat + so_portion + rev2007, data = npc,  
                  index = c("gvkey","fyear"), estimator = "meanfirst", demean = FALSE, vce = "jackknife")
print(out.kbal)
panelview(data = germany, gdp ~ treat, index = c("country","year"), 
          xlab = "Year", ylab = "Country", by.timing = TRUE,
          axis.adjust = TRUE, axis.lab.gap = c(1,0))
out2.mbal <- tjbal(data = germany, Y = "gdp", D = "treat", Y.match.time = c(1960:1990), 
                   X = c("trade","infrate","industry", "schooling", "invest80"), 
                   X.avg.time = list(c(1981:1990),c(1981:1990),c(1981:1990),c(1980:1985),c(1980)),
                   index = c("country","year"), demean = TRUE)
plot(out2.mbal, ylim = c(-5000,5000))
