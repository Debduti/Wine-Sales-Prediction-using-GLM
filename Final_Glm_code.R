#Part 0: Load & Prepare Data

library(readr)
library(dplyr)
library(zoo)
library(psych)
library(ROCR)
library(corrplot)
library(car)
library(InformationValue)
library(rJava)
library(pbkrtest)
library(car)
library(leaps)
library(MASS)
library(corrplot)
library(glm2)
library(aod)
library(mice)
library(Hmisc)
library(xlsxjars)
library(xlsx)

library(VIM)
library(pROC)
library(pscl) # For "counting" models (e.g., Poisson and Negative Binomial)
library(ggplot2) # For graphical tools
library(readr)
library(corrplot)


#setwd

# Read the wine dataset

wine=read.csv("Wine_Training.csv",header=T)

head(wine,1)

# Part 1: Data Exploration

#Data Quality Check
str(wine)
summary(wine)

library(Hmisc)
describe(wine)

#TARGET
par(mfrow=c(1,2))
hist(wine$TARGET, col = "#A71990", xlab = "Target ", main = "histogram of wine sales")
boxplot(wine$TARGET, col = "#A71990", main = "boxplot of wine sales")
par(mfrow = c(1,1))
#Chemistry 

# FixedAcidity and VolatileAcidity
dev.off()
par("mar")
par(mar=c(3,1,1,1))
par(mfrow=c(2,2))
hist(wine$FixedAcidity, col = "#A71930", xlab ="FixedAcidity", main = "Histogram of FixedAcidity")
hist(wine$VolatileAcidity, col = "#A71990", xlab = "VolatileAcidity", main = "Histogram of VolatileAcidity")
boxplot(wine$FixedAcidity, col = "#A71930", main = "Boxplot of FixedAcidity")
boxplot(wine$VolatileAcidity, col = "#A71990", main = "Boxplot of VolatileAcidity")
par(mfrow=c(1,1))

# CitricAcid and ResidualSugar
par(mfrow=c(2,2))
hist(wine$CitricAcid, col = "#D77730", xlab = "CitricAcid", main = "Histogram of CitricAcid")
hist(wine$ResidualSugar, col = "#ABCEAC", xlab = "ResidualSugar ", main = "Histogram of ResidualSugar")
boxplot(wine$CitricAcid, col = "#D77730", main = "Boxplot of CitricAcid")
boxplot(wine$ResidualSugar, col = "#ABCEAC", main = "Boxplot of ResidualSugar")
par(mfrow=c(1,1))

#Chlorides and FreeSulfur Dioxide
par(mfrow=c(2,2))
hist(wine$Chlorides, col = "#A71930", xlab = "Chlorides", main = "Histogram of Chlorides")
hist(wine$FreeSulfurDioxide, col = "#DBCEAC", xlab = "FreeSulfurDioxide ", main = "Histogram of FreeSulfurDioxide")
boxplot(wine$Chlorides, col = "#A71930", main = "Boxplot of Chlorides")
boxplot(wine$FreeSulfurDioxide, col = "#DBCEAC", main = "Boxplot of FreeSulfurDioxide")
par(mfrow=c(1,1))

#TotalSulfurDioxide and Density
par(mfrow=c(1,1))
hist(wine$TotalSulfurDioxide, col = "#D77730", xlab = "TotalSulfurDioxide", main = "Histogram of TotalSulfurDioxide")
hist(wine$Density, col = "#DBCEAC", xlab = "Density", main = "Histogram of Density")
boxplot(wine$TotalSulfurDioxide, col = "#D77730", main = "Boxplot of TotalSulfurDioxide")
boxplot(wine$Density, col = "#DBCEAC", main = "Boxplot of Density")
par(mfrow=c(1,1))

#pH and Sulphates
par(mfrow=c(2,2))
hist(wine$pH, col = "#A71930", xlab = "pH", main = "Histogram of pH")
hist(wine$Sulphates, col = "#09ADAD", xlab = "Sulphates", main = "Histograms of Sulphates")
boxplot(wine$pH, col = "#A71930", main = "Boxplot of pH")
boxplot(wine$Sulphates, col = "#09ADAD", main = "Boxplot of Sulphates")
par(mfrow=c(1,1))


#Alcohol and Acid Index
par(mfrow=c(2,2))
hist(wine$Alcohol, col = "#A71930", xlab = "Alcohol", main = "Histogram of Alcohol")
hist(wine$AcidIndex, col = "#DBCEAC", xlab = "AcidIndex", main = "Histograms of AcidIndex")
boxplot(wine$Alcohol, col = "#A71930", main = "Boxplot of Alcohol")
boxplot(wine$AcidIndex, col = "#DBCEAC", main = "Boxplot of AcidIndex")
par(mfrow=c(1,1))

#Label Appeal and STARS
par(mfrow=c(2,2))
hist(wine$LabelAppeal, col = "#DBCEAC", xlab = "LabelAppeal", main = "Histogram of LabelAppeal ")
hist(wine$STARS, col = "#09ADAD", xlab = "STARS", main = "Histogram of STARS")
boxplot(wine$LabelAppeal, col = "#DBCEAC", main = "Boxplot of LabelAppeal")
boxplot(wine$STARS, col = "#09ADAD", main = "Boxplot of STARS")
par(mfrow=c(1,1))



#############################################################################################
##Part 2: Data Preparation

#Check missing data percentage
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(wine,2,pMiss)

# Outliers, defined as $Q1 - (1.5*IQR)$ and $Q3 + (1.5*IQR)$ are removed from the data set and flag
#variables are appended in a similar fashion to the missing data flags. Finally, 
#all missing and deleted variables are reimputed using regression trees.

wine$NoResidualSugar <- 0
wine$NoResidualSugar [is.na(wine$ResidualSugar)] <- 1

wine$NoChlorides  <- 0
wine$NoChlorides [is.na(wine$Chlorides)] <- 1

wine$NoFreeSulfurDioxide <- 0
wine$NoFreeSulfurDioxide[is.na(wine$FreeSulfurDioxide)] <- 1

wine$NoTotalSulfurDioxide <- 0
wine$NoTotalSulfurDioxide[is.na(wine$TotalSulfurDioxide)] <- 1

wine$NopH <- 0
wine$NopH[is.na(wine$pH)] <- 1

wine$NoSulphates <- 0
wine$NoSulphates [is.na(wine$Sulphates)] <- 1

wine$NoResidualSugar <- 0
wine$NoResidualSugar [is.na(wine$ResidualSugar)] <- 1

wine$NoAlcohol <- 0
wine$NoAlcohol [is.na(wine$Alcohol)] <- 1

wine$NoSTARS<- 0
wine$NoSTARS [is.na(wine$STARS)] <- 1

str(wine)



colnames(wine) <- tolower(colnames(wine))

train.flagged <- wine%>%
  mutate_all(funs(na.flag = ifelse(is.na(.),1,0)))
int_df <- train.flagged%>%
  dplyr::select(-index, -target, -index_na.flag, -target_na.flag)%>%
  dplyr::select_if(is.numeric)
# md.pattern(int_df)
cleaned_cols <- list()

head(train.flagged,10)

for(c in colnames(wine%>%
                  dplyr::select(-index, -target)%>%
                  dplyr::select_if(is.numeric))){
  column <- train.flagged%>%select_(col = c)
  iqr <- quantile(column$col, na.rm = T)[4] - quantile(column$col, na.rm = T)[2]
  low <- quantile(column$col, na.rm = T)[2] - iqr
  high <- quantile(column$col, na.rm = T)[4] + iqr
  
  vals <- c()
  for(i in seq(1:nrow(int_df))){
    ifelse(between(column$col[i], low - (1.5*iqr), high + (1.5*iqr)),
           vals[i] <- column$col[i], 
           ifelse(is.na(column$col[i]), vals[i] <- NA, vals[i] <- NA))
  }
  
  ifelse(length(vals) == nrow(int_df),
         cleaned_cols[[c]] <- vals, 
         cleaned_cols[[c]] <- c(vals,NA))
}


df2 <- bind_cols(
  bind_cols(cleaned_cols)%>%
    scale(center = TRUE)%>%
    data.frame(),
  train.flagged%>%
    dplyr::select(ends_with('na.flag'))%>%
    dplyr::select(-index_na.flag, -target_na.flag)
)



df3 <- df2%>%
  mutate(
    fixedacidity_out.flag = ifelse(is.na(fixedacidity) & fixedacidity_na.flag ==0,1,0),
    volatileacidity_out.flag = ifelse(is.na(volatileacidity) & volatileacidity_na.flag ==0,1,0),
    citricacid_out.flag = ifelse(is.na(citricacid) & citricacid_na.flag ==0,1,0),
    residualsugar_out.flag = ifelse(is.na(residualsugar) & residualsugar_na.flag ==0,1,0),
    chlorides_out.flag = ifelse(is.na(chlorides) & chlorides_na.flag ==0,1,0),
    freesulfurdioxide_out.flag = ifelse(is.na(freesulfurdioxide) & freesulfurdioxide_na.flag ==0,1,0),
    totalsulfurdioxide_out.flag = ifelse(is.na(totalsulfurdioxide) & totalsulfurdioxide_na.flag ==0,1,0),
    density_out.flag = ifelse(is.na(density) & density_na.flag ==0,1,0),
    ph_out.flag = ifelse(is.na(ph) & ph_na.flag ==0,1,0),
    sulphates_out.flag = ifelse(is.na(sulphates) & sulphates_na.flag ==0,1,0),
    alcohol_out.flag = ifelse(is.na(alcohol) & alcohol_na.flag ==0,1,0),
    labelappeal_out.flag = ifelse(is.na(labelappeal) & labelappeal_na.flag ==0,1,0),
    acidindex_out.flag = ifelse(is.na(acidindex) & acidindex_na.flag ==0,1,0),
    stars_out.flag = ifelse(is.na(stars) & stars_na.flag ==0,1,0)
  )

?mice

library(mice)
temp_df <- mice(df3, method = 'cart', maxit = 1)
train <- complete(temp_df)%>%
  bind_cols(wine%>%dplyr::select(index, target))%>%
  dplyr::select(-stars_out.flag, -labelappeal_out.flag, -density_na.flag,
                -labelappeal_na.flag, -acidindex_na.flag, -fixedacidity_na.flag,
                -volatileacidity_na.flag, -citricacid_na.flag)


summary(train)

names(train)



################################################################################################


#Correlation Matrix
subdatnum22 <- subset(train, select=c(
  'fixedacidity',
  'volatileacidity',
  'citricacid',
  'residualsugar',
  'chlorides',
  'freesulfurdioxide',
  'totalsulfurdioxide',
  'density',
  'ph',
  'sulphates',
  'alcohol',
  'labelappeal',
  'acidindex',
  'stars',
  'residualsugar_out.flag',
  'chlorides_na.flag',
  'freesulfurdioxide_out.flag',
  'totalsulfurdioxide_out.flag',
  'ph_out.flag',
  'sulphates_out.flag',
  'stars_na.flag',
  'alcohol_out.flag',
  'target'))

require(corrplot)
mcor <- cor(subdatnum22)
corrplot(mcor, method="number", shade.col=NA, tl.col="black",tl.cex=0.8)
par(mfrow=c(1,1)) 




########################################################################################################################
# Part3: Model Building

#Model 1:Poisson 

library(MASS)



### mlr###

stepwisemodel <- lm(formula = target ~ ., data = train)
stepwise <- stepAIC(stepwisemodel, direction = "both",trace = 0)
summary(stepwise)
anova(stepwise)
summary(stepwise)
#####

base_poission <- glm(target ~ ., family="poisson", data=train)
summary(base_poission)
#Using AIC Stepwise for variable selection
poisson.back <- stepAIC(base_poission, direction = 'backward',trace=0)
summary(poisson.back)
poiss.model <- glm(formula = target ~ volatileacidity + +freesulfurdioxide+ totalsulfurdioxide + 
                     ph + sulphates + chlorides+ alcohol + labelappeal + acidindex + 
                     stars + stars_na.flag + acidindex_out.flag, 
                   family = "poisson", data = train)

poisson.coeffs <- data.frame(var = names(poiss.model$coefficients),
                                  coefficient = poiss.model$coefficients)%>%
  mutate(method = 'Poisson')

anova(poiss.model, test="Chisq")

with(poiss.model, cbind(res.deviance = deviance, df = df.residual,
                        p = pchisq(deviance, df.residual, lower.tail=FALSE)))


poiss.model1 <- glm(formula = target ~ volatileacidity + totalsulfurdioxide + 
                     sulphates+ labelappeal + acidindex + 
                     stars + stars_na.flag + volatileacidity_out.flag + acidindex_out.flag, 
                   family = "poisson", data = train)

poisson.coeffs1 <- data.frame(var = names(poiss.model1$coefficients),
                             coefficient = poiss.model1$coefficients)%>%
  mutate(method = 'Poisson')

anova(poiss.model1, test="Chisq")

with(poiss.model1, cbind(res.deviance = deviance, df = df.residual,
                        p = pchisq(deviance, df.residual, lower.tail=FALSE)))


library(AER)
deviance(poiss.model)/poiss.model$df.residual
dispersiontest(poiss.model)

#what type of dispersion does sample have?
mean(train$target)
var(train$target)

library(car)
influencePlot(poiss.model)
res <- residuals(poiss.model, type="deviance")
plot(log(predict(poiss.model)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)

########################################################################################################

#Model2: Negative Binomial

base_nb <- glm.nb(target ~ ., data=train)
#Using AIC Stepwise for variable selection
nb.back <- stepAIC(base_nb, direction = 'backward',trace=0)

summary(nb.back)

negbinomial.mod <- glm.nb(formula = target ~ volatileacidity + totalsulfurdioxide +
                            ph + sulphates + alcohol + labelappeal + acidindex + 
                            stars + stars_na.flag + volatileacidity_out.flag + acidindex_out.flag, 
                          data = train)

negbinomial.coeffs <- data.frame(var = names(negbinomial.mod$coefficients),
                                 coefficient = negbinomial.mod$coefficients)%>%
  mutate(method = 'Negative Binomial')
odTest(negbinomial.mod)

summary(negbinomial.mod)

library(ggplot2)

ggplot(bind_rows(negbinomial.coeffs, poisson.coeffs),
       aes(x = reorder(var, coefficient), y = coefficient, fill = method))+
  geom_col(position = 'dodge')+
  coord_flip()+
  labs(y = 'Coefficient',
       x = element_blank(),
       title = 'Regression Coefficients',
       subtitle = 'Features selected via backwards variable selection')



#theme_gray()

#############################################

# Zero Inflated Regression

# Zero Inflated Poisson

zinp.mod <- pscl::zeroinfl(formula = target ~ volatileacidity + totalsulfurdioxide + 
                             ph + sulphates + alcohol + labelappeal + acidindex + 
                             stars + stars_na.flag + volatileacidity_out.flag + acidindex_out.flag, 
                           data = train)
zinp.coeffs <- data.frame(var = names(zinp.mod$coefficients$zero),
                          coefficient = zinp.mod$coefficients$zero)%>%
  mutate(method = 'Zero-Inflated Poisson')
summary(zinp.mod)

vuong(poiss.model, zinp.mod)

#Zero Inflated Neg Binom

zinng.mod <- pscl::zeroinfl(formula = target ~ volatileacidity + totalsulfurdioxide + 
                              ph + sulphates + alcohol + labelappeal + acidindex + 
                              stars + stars_na.flag + volatileacidity_out.flag + acidindex_out.flag, 
                            data = train, dist = "negbin")
zinng.coeffs <- data.frame(var = names(zinng.mod$coefficients$zero),
                           coefficient = zinng.mod$coefficients$zero)%>%
  mutate(method = 'Zero-Inflated Negative Binomial')

summary(zinng.mod)

vuong(negbinomial.mod, zinng.mod)

# Comparing Coefficients of Zero Inflated Models

ggplot(bind_rows(zinp.coeffs, zinng.coeffs),
       aes(x = reorder(var, coefficient), y = coefficient, fill = method))+
  geom_col(position = 'dodge')+
  coord_flip()+
  labs(y = 'Coefficient',
       x = element_blank(),
       title = 'Regression Coefficients',
       subtitle = 'Features selected via backwards variable selection')+
  
  theme_gray()


####################################################################################
# Model Evaluation
# Calculating mae and rmse on full train data.We also calculate AIC and BIC score s of all
# the models

library(ModelMetrics)

columns <- c('Poisson', 'Negative Binomial','Zero-Inflated Poisson','Zero-Inflated Negative Binomial')

poiss.mae <- mae(train$target, predict(poiss.model))
poiss.rmse <- rmse(train$target, predict(poiss.model))
AIC.poiss <- AIC(poiss.model)
BIC.poiss <- BIC(poiss.model)

negbin.mae <- mae(train$target, predict(negbinomial.mod))
negbin.rmse <- rmse(train$target, predict(negbinomial.mod))
AIC.nbr <- AIC(negbinomial.mod)
BIC.nbr <- BIC(negbinomial.mod)

zinp.mae <- mae(train$target, predict(zinp.mod))
zinp.rmse <- rmse(train$target, predict(zinp.mod))
AIC.zinp <- AIC(zinp.mod)
BIC.zinp <- BIC(zinp.mod)

zinng.mae <- mae(train$target, predict(zinng.mod))
zinng.rmse <- rmse(train$target, predict(zinng.mod))
AIC.zinng <- AIC(zinng.mod)
BIC.zinng <- BIC(zinng.mod)

data.frame(
  columns,
  rmse = c(poiss.rmse, negbin.rmse, zinp.rmse, zinng.rmse),
  mae = c(poiss.mae, negbin.mae, zinp.mae, zinng.mae),
  AIC = c(AIC.poiss,AIC.nbr,AIC.zinp,AIC.zinng),
  BIC = c(BIC.poiss,BIC.nbr,BIC.zinp,BIC.zinng)
)



