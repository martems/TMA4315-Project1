library(car)
library(GGally)
data(SLID,package = "car")
SLID = SLID[complete.cases(SLID), ]
ds = SLID
colnames(ds)
dim(ds)
summary(ds)
levels(ds$sex)
levels(ds$language)
ggpairs(ds)
lm = lm(formula = wages ~ education+age+sex+language,data=SLID)
summary(lm)
