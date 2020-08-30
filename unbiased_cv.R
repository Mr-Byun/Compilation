setwd("C:/Users/britz/Desktop/Research/°£¾Ï/")

source("source.R")

df <- read(exc.clinical=T, split.hosp=F)

############
# Untuned  #
rf <- iter_tune(df, iter=100,
                repeats=1, folds=10, # if method="none" is used, repeats and folds are not applied for parameter tuning
                search="none", model="ranger",
                stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                elimFunc=NULL)

# reduced_rf <- iter_tune(df, iter=100,
#                         repeats=1, folds=10,
#                         search="none", model="ranger",
#                         stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
#                         elimFunc=NULL, p=.5)

svm <- iter_tune(df, iter=100,
                 repeats=1, folds=10,
                 search="none", model="svmRadial",
                 stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                 elimFunc=NULL)

glm <- iter_tune(df, iter=100,
                 repeats=1, folds=10,
                 search="none", model="glm",
                 stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                 elimFunc=NULL)

#########
# Tuned #
tuned_rf <- iter_tune(df, iter=100,
                      repeats=1, folds=10,
                      search="repeatedcv", model="ranger",
                      stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                      elimFunc=NULL)

tuned_svm <- iter_tune(df, iter=100,
                       repeats=1, folds=10,
                       search="repeatedcv", model="svmRadial",
                       stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                       elimFunc=NULL)

glmnet <- iter_tune(df, iter=100,
                    repeats=1, folds=10,
                    search="repeatedcv", model="glmnet",
                    stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=glmnet.tuning,
                    elimFunc=NULL)

###########################
# Lasso Feature Selection #
# For Lasso feature selection or RFE, pass ellipsis as (fe.search, fe.repeats, fe.folds)
# Available elimFunc are (glmnet, rf, lr)
lasso_rf <- iter_tune(df, iter=100,
                      repeats=3, folds=10,
                      search="none", model="ranger",
                      stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                      elimFunc="glmnet", fe.search="repeatedcv", fe.repeats=1, fe.folds=10)

lasso_svm <- iter_tune(df, iter=100,
                       repeats=3, folds=10,
                       search="none", model="svmRadial",
                       stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                       elimFunc="glmnet", fe.search="repeatedcv", fe.repeats=1, fe.folds=10)

lasso_glm <- iter_tune(df, iter=100,
                       repeats=3, folds=10,
                       search="none", model="glm",
                       stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=glm.tuning,
                       elimFunc="glmnet", fe.search="repeatedcv", fe.repeats=1, fe.folds=10)


#######
# RFE #
rfe_rf <- iter_tune(df, iter=100,
                    repeats=3, folds=10,
                    search="none", model="ranger",
                    stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                    elimFunc="rf", fe.search="repeatedcv", fe.repeats=1, fe.folds=10)
save.image('temp.RData')

rfe_svm <- iter_tune(df, iter=100,
                     repeats=3, folds=10,
                     search="none", model="svmRadial",
                     stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                     elimFunc="rf", fe.search="repeatedcv", fe.repeats=1, fe.folds=10)
save.image('temp.RData')

rfe_glm <- iter_tune(df, iter=100,
                     repeats=3, folds=10,
                     search="none", model="glm",
                     stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=glm.tuning,
                     elimFunc=rfe_rf, fe.search="repeatedcv", fe.repeats=1, fe.folds=10)
save.image('temp.RData')

####################
# Get Final Result #
final <- as.one.df(order=c("glm", "lasso_glm", "glmnet", "rfe_glm",
                           "rf", "lasso_rf", "tuned_rf", "rfe_rf",
                           "svm", "lasso_svm", "tuned_svm", "rfe_svm"),
                   glm, lasso_glm, glmnet, rfe_glm,
                   rf, lasso_rf, tuned_rf, rfe_rf,
                   svm, lasso_svm, tuned_svm, rfe_svm)

means <- aggregate(AUC ~ Group, final$Test, mean)

ggplot(data = final$Test, aes(x=Group, y=AUC, fill=lGroup)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape=1, size=2, color="red", fill="red") +
  stat_summary(fun=mean, colour="red", geom="text",  show.legend=F,
               vjust=-0.7, aes( label=round(..y.., digits=3))) +
  scale_x_discrete(labels=c("LR\n(182 predictors)", "LR w/ Lasso Selection", "Tuned Lasso\n(182 predictors)", "LR w/ RFE",
                            "RF\n(182 predictors)", "RF w/ Lasso Selection", "Tuned RF\n(182 predictors)", "RF w/ RFE",
                            "SVM\n(182 predictors)", "SVM w/ Lasso Selection", "Tuned SVM\n(182 predictors)", "SVM w/ RFE")) +
  scale_fill_discrete(name = "Model Type", labels = c("Logistic Regression", "Random Forest", "Support Vector Machine"))


#######################
# Univariate Analysis #
variables <- var.analysis(rfe_rf, plot=F)

uni <- uni.analysis(df, variables$optVariables,
                    iter=50,
                    repeats=1, folds=10,
                    search="none", model="glm",
                    stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=tuning,
                    elimFunc=NULL)

uni.all <- rbind(interpolate_AUC(uni[["Ser"]], "Ser"),
           rbind(interpolate_AUC(uni[["C10.2"]], "C10.2"),
           rbind(interpolate_AUC(uni[["AFP"]], "AFP"),
           rbind(interpolate_AUC(uni[["C10"]], "C10"),
           interpolate_AUC(uni[["Ser"]], "Asp")))))

plot_sdAUC(uni.all)
