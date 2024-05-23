# Objective
In this project, our aim is to analyze the Wine dataset and build a model to predict the number
of wine cases sold, based on its features, most of which are chemical properties.  

This dataset comes from a Kaggle competition for wine sales prediction and the dataset can be found on
Github.  

We will cover the following topics in this project:  

**a**. We will do an initial analysis to find out nulls in the dataset and process them   
**b**. We will do exploratory analysis on all the variables and produce relevant boxplots, 
scatter plots and histograms     
**c**. We will then do some feature engineering that will help with building our models   
**d**. We will build four different models- **Poisson Regression Model, Negative Binomial Model, Zero Inflated Poisson Model and Zero Inflated Negative Binomial Model**.     
We will compare their performances based on metrics such as Root Mean Squared Error, 
**Mean Absolute Error, AIC and BIC** and select a final best model   

# Dataset

The dataset Wine.csv contains 12795 records and 14 predictors (excluding INDEX). Our response variable is TARGET. The dataset is obtained from a Kaggle competition for prediction of wine sales based on its attributes.

![image](https://github.com/Debduti/Wine-Sales-Prediction-using-GLM/assets/58540839/8335dd71-05ca-4ea1-989b-1e4f6ccca622)  

# Data Description  
Attached below is the data description of the Wine.csv dataset.  
![image](https://github.com/Debduti/Wine-Sales-Prediction-using-GLM/assets/58540839/4cc83aae-daf4-411d-a828-ed1352324d63)  

The following is a summary of the dataset.
![image](https://github.com/Debduti/Wine-Sales-Prediction-using-GLM/assets/58540839/7fbae8db-47fd-49df-b740-46a632db397d)  

The above summary statistics reveals that there are **NULLS** in **ResidualSugar, Chlorides,
FreeSulfurDioxide, TotalSulfurDioxide, pH, Sulphates, Alcohol, and Stars**. We also seem to
have outliers in some of the columns. We will take a closer look into these in our exploratory
data analysis and clean our data accordingly.

# Methods
We conducted the following steps for our project:  
**Exploratory Data Analysis:** We checked the distribution of all the variables and created
histograms and box plots for them to visualize them. We also plotted the distribution of
percentage of missing values in the dataset. We also created a correlation matrix and scatter
plots to understand the correlation of the features with the TARGET.  

**Data Preparation:** Outliers, as observed from the histogram of predictor variables are
removed from the dataset. Binary flags are added to each variable to indicate whether they
have NULLS. The missing values are re-imputed using predictive mean matching.  

**Data Analysis:** We fitted four different models- Poisson Regression , Negative Binomial,
Zero Inflated Poisson Regression and Zero Inflated Negative Binomial Regression models.
We ran a full model with all the variables, and then we did a stepwise backward selection
using AIC to select the significant predictors. This process was done for both Poisson and
Negative Binomial Models. The goodness of fit was checked using Pearson Chi-Square test
for the models. Since there was slight overdispersion and there was a large number of zeros
in our response variable, we proceeded to run zero-inflated regression models for both
Poisson and Negative Binomial Regression Models. We ran the Vuong test between our base
model and zero inflated model to determine whether the zero inflated model had any
difference with the base model.  

Finally, we compared all four models based on different metrics like their **Root Mean
Squared Error, Mean Absolute Error, AIC, and BIC** and chose the best performing model.

# Results
Detailed findings and results of the modeling has been discussed in the attached report  

**a. Regression Coefficients of the Models**
In Figure 17 we see the plot depicting the regression coefficients of the Poisson vs
Negative Binomial Models. Both models assign almost similar coefficients to the
predictors. Missing values for stars and acidindex values significantly harm sales
while labelappeal and present higher values for stars lead to higher sales.  

**![image](https://github.com/Debduti/Wine-Sales-Prediction-using-GLM/assets/58540839/2938bb50-a41c-408b-964d-89fe3b1fa249)
**  
Both sets of models use the same predictors. Some coefficients have been flipped in the
zero-inflated models, Missing Stars and acidity values are strong positive influences on sales.
Label appeal is similarly positive as in previous models, but actual star ratings are negative.
Our assumption is that higher star ratings also come with a higher price tag, which can
explain the negative coefficient for STARS.  

**b. Model Evaluation**
The four models have been evaluated based on Root Mean Square Error(rmse), Mean
Absolute Error(MAE), AIC, and BIC.
![image](https://github.com/Debduti/Wine-Sales-Prediction-using-GLM/assets/58540839/b5a7c30c-f050-4f55-9b3b-93bf031a625c)  
The error rates for the Zero Inflated models are lower than the regular models. Although the
AIC and BIC values are slightly better for the Zero Inflated Poisson Models, our inference is
to choose the Zero Inflated Negative Binomial Model because of the slight overdispersion
present in the response variable and large number of zeros in the response variable.  

# References  

[1] Cameron A. C. and Trivedi P. K., (2013). Regression Analysis of Count Data, Second
Edition, Econometric Society Monograph No. 53, Cambridge University Press, Cambridge.  
[2] Dunn, P.K.; Smyth, G.K. (2018). Generalized Linear Models with Examples in R. New
York: Springer. doi:10.1007/978-1-4419-0118-7.  
[3] Kida, Y. (2019). Generalized Linear Models - Introduction to advanced statistical
modeling. https://towardsdatascience.com/generalized-linear-models-9cbf848bb8ab.  







