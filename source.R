library(caret); library(caretEnsemble); library(doParallel);
library(tidyr); library(stringr); library(ROCR); library(gridExtra);

options(scipen=999) # Disable scientific notation

read <- function(exc.clinical=F, split.hosp=F)
{
  df <- read.csv("./data/181_metabolite_clinical.csv")
  df <- df[, -1]
  exc.col <- colnames(df)[unlist(lapply(df, function(x) any(is.na(x))))]
  exc.col <- c(exc.col, c("disease", "Group", "hospital", "hospital.binary", "age.binary", "setID", "Azelaic.acid"))
  if(split.hosp) exc.col <- exc.col[!exc.col == "hospital"]
  
  df <- df[, !colnames(df) %in% exc.col]
  
  df$HCC <- factor(ifelse(df$HCC == 1, "HCC", "Cirrhosis"))
  df$sex <- factor(ifelse(df$sex == 1, "M", "F"))
  df$etiology <- factor(ifelse(df$etiology == 1, "one", ifelse(df$etiology == 2, "two", "three")))
  
  # rownames(df) <- df$Name
  
  if(exc.clinical) df <- subset(df, select=c(-sex, -etiology, -age))
  if(split.hosp)
  {
    df <- lapply(dplyr::group_split(df, hospital), as.data.frame)
    df <- lapply(df, function(x) { x["hospital"] <- NULL; x})
  }
  
  return(df)
}

## model  = "ranger/rf" / "svmLinear" / "svmPoly" / "svmRadial"
tuning <- function(df, search, model, repeats=3, folds=10, stratified=T, parallel=T, verbose=T)
{
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  # on.exit(registerDoSEQ())
  
  if(stratified) stratified.index <- createMultiFolds(df$HCC, k=folds, times=repeats)
  else stratified.index <- NULL
  
  results <- NULL
  if(search != "none")
  {
    if(model == "ranger") importance <- "impurity"
    else if(model == "rf") importance <- T
    else importance <- NULL
    
    trnControl <- trainControl(method = search,
                               number = folds,
                               repeats = repeats,
                               index = stratified.index,
                               search = "random",
                               classProbs = T,
                               summaryFunction = twoClassSummary,
                               returnResamp= 'final',
                               savePredictions = 'final',
                               allowParallel = parallel,
                               verbose=verbose)
    
    results <- train(HCC ~ ., data = df,
                     method = model,
                     metric = "ROC",
                     # tuneLength = 20,
                     importance = importance, # only for RF
                     trControl = trnControl,
                     num.threads = 4)
  }
  else
  { # If method="none", caluclate default parameter for RF and SVM
    param <- estimate_param(model, df)
    
    trnControl <- trainControl(method = search,
                               classProbs = T,
                               summaryFunction = twoClassSummary,
                               verbose=verbose)
    
    results <- train(HCC ~ ., data = df,
                     method = model,
                     tuneGrid = param,
                     metric = "ROC",
                     trControl = trnControl)
  }
  
  if(is.null(results)) stop("Incorrect method")
  return(results)
}

glm.tuning <- function(df, search, model, repeats=3, folds=10, stratified=T, parallel=F, verbose=T)
{
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  # on.exit(registerDoSEQ())
  
  if(stratified)
    stratified.index <- createMultiFolds(df$HCC, k=folds, times=repeats)
  else
    stratified.index <- NULL
  

  trnControl <- trainControl(method = search,
                             number = folds,
                             repeats = repeats,
                             index = stratified.index,
                             search = "random",
                             classProbs = T,
                             # savePredictions = T,
                             summaryFunction = twoClassSummary,
                             returnResamp = 'final',
                             savePredictions = 'final',
                             verbose=verbose)
    
  results <- train(HCC ~ ., data = df,
                   method = model,
                   family = 'binomial',
                   metric = "ROC",
                   trControl = trnControl)
  
  return(results)
}

stacking <- function(df, search, model, repeats=3, folds=10, stratified=T, parallel=T, verbose=T)
{
  df$HCC <- factor(ifelse(df$HCC == 1, "HCC", "Normal"))
  
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  if(any(str_detect(base_models, "glmnet")))
  {
    lambdas <- seq(1,0,-0.01)
    grid <- expand.grid(alpha = 1, lambda = lambdas)
    glm.tune <- list(
      glmnet=caretModelSpec(method="glmnet", tuneGrid=grid)
    )
    base_models <- base_models[!base_models %in% grep(paste0("glmnet", collapse = "|"), base_models, value = T)]
  }
  else
    glm.tune <- NULL
  
  
  stratified.index <- createMultiFolds(df$HCC, k=folds, times=repeats)
  baseCtrl <- trainControl(method = "repeatedcv",
                           number = folds,
                           repeats = repeats,
                           index = stratified.index,
                           summaryFunction = twoClassSummary,
                           savePredictions = "final",
                           classProbs = T,
                           allowParallel = T,
                           verbose=T)
  
  metaCtrl <- trainControl(method = "repeatedcv",
                           number = folds,
                           repeats = repeats,
                           summaryFunction = twoClassSummary,
                           savePredictions = "final",
                           classProbs = T,
                           allowParallel = parallel,
                           verbose=T)
  
  if(verbose) print("Training Base-models")
  models <- caretList(HCC ~ .,
                      data = df,
                      trControl = baseCtrl,
                      metric = "ROC",
                      methodList = base_models,
                      tuneList = glm.tune)
  stopCluster(cl)
  
  # results <- resamples(models)
  # summary(results)
  
  if(verbose) print("Training Meta-learner")
  
  registerDoParallel(cl)
  stack <- caretStack(models,
                      method = meta_learner,
                      metric = "ROC",
                      trControl = metaCtrl)
  
  return(stack)
}

voting <- function(stack)
{
  df <- data.frame(row.names=1:nrow(stack$models[[1]]$pred))
  for(i in 1:length(stack$models))
  {
    temp <- stack$models[[i]]$pred
    temp <- temp[order(temp$Resample, temp$rowIndex), ]
    temp$Resample <- factor(temp$Resample)
    temp <- subset(temp, select=c(obs, HCC, rowIndex, Resample))
    temp$obs <- factor(ifelse(temp$obs == "HCC", 1, 0))
    colnames(temp) <- paste0(names(stack$models)[i], "_", colnames(temp))
    
    df <- cbind(df, temp)
  }
  
  temp <- df[, str_detect(colnames(df), "rowIndex")]
  for(i in 1:(length(stack$models)-1))
  {
    if(!identical(temp[, i], temp[, i+1]))
      stop("Unidentical index between resamples")
  }
  
  temp <- df[, str_detect(colnames(df), "HCC")]
  df$Voting <- rowMeans(temp)
  
  df <- data.frame(obs=df[, 1],
                   Voting=df$Voting,
                   Resample=df[, 4])
  
  res <- c()
  for(i in levels(df$Resample))
  {
    temp <- df[str_detect(df$Resample, i), ]
    
    res <- c(res, MLmetrics::AUC(temp$Voting, temp$obs))
  }
  
  return(res)
}

glm.coef <- function(glm, s=NULL)
{
  if(is.null(s)) s <- glm$bestTune$lambda
    
  coef_f <- coef(glm$finalModel,s=s)
  coef <- data.frame(coef.name = dimnames(coef_f)[[1]], 
                     coef.value = matrix(coef_f))
  
  coef <- coef[-1,]
  coef$coef.name <- as.character(coef$coef.name)
  
  coef <- coef[!(coef$coef.value == 0), ]
  coef <- dplyr::arrange(coef,-coef.value)
  
  return(coef)
}

recursive.elim <- function(df, funcs="rf", fe.search, fe.repeats=1, fe.folds, parallel=T)
{
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  if(funcs == "rf")
  {
    funcs <- rfFuncs
    funcs$fit <- function (x, y, first, last, ...){
      loadNamespace("randomForest")
      randomForest::randomForest(x, y, ...)}
  }
  else if(funcs == "lr") funcs <- lrFuncs
  funcs$summary <- twoClassSummary
  
  ctrl <- rfeControl(functions = funcs,
                     method = fe.search,
                     number = fe.folds,
                     repeats = fe.repeats,
                     rerank = T,
                     saveDetails = T,
                     allowParallel = parallel,
                     verbose = T)
  
  Profile <- rfe(df[,-1], df$HCC,
                 sizes = c(1:100, ncol(df[, -1])),
                 rfeControl = ctrl,
                 metric="ROC",
                 num.threads = 4)
  
  Profile <- list(optsize = Profile$optsize,
                  optVariables = Profile$optVariables,
                  variables = Profile$variables,
                  results = Profile$results)
  
  return(Profile)
}

glmnet.tuning <- function(df, fe.search, model, fe.repeats=1, fe.folds, stratified=T, parallel=T, verbose=T, alpha = 1, returnVar=F)
{
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  # on.exit(registerDoSEQ())
  
  lambdas <- seq(0.0001, 1, length = 100)
  grid <- expand.grid(alpha = alpha, lambda = lambdas)
  if(alpha >= 2)
    grid <- NULL
  
  if(stratified)
    stratified.index <- createMultiFolds(df$HCC, k=fe.folds, times=fe.repeats)
  else
    stratified.index <- NULL
  
  trnControl <- trainControl(method = fe.search,
                             number = fe.folds,
                             repeats = fe.repeats,
                             index = stratified.index,
                             search = "random",
                             classProbs = T,
                             summaryFunction = twoClassSummary,
                             returnResamp = 'final',
                             savePredictions = 'final',
                             allowParallel = parallel,
                             verbose=verbose)
  
  final_model <- train(HCC ~ ., data = df,
                       method = model,
                       family = 'binomial',
                       metric = "ROC",
                       # maximize=F,
                       trControl = trnControl,
                       tuneGrid = grid,
                       num.threads = 4)
  
  if(returnVar)
  {
    coef.cv <- glm.coef(final_model)
    
    results <- list(optsize = length(coef.cv$coef.name),
                    optVariables = coef.cv$coef.name,
                    variables = coef.cv$coef.value,
                    results = final_model$results)
  }
  else results <- final_model
  
  return(results)
}

# search = "none" to disable tuning for train set
iter_tune <- function(df, iter=100,
                      search, model,
                      repeats=3, folds=10,
                      stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=NULL,
                      elimFunc=NULL, p=.66, ...)
{
  result <- list()
  
  for(i in 1:iter)
  {
    if(stratified)
    {
      train.index <- createDataPartition(df$HCC, p = .66, list = FALSE)
      train <- df[ train.index, ]
      test  <- df[-train.index, ]
    }
    else # random sampling
    {
      train.index <- sample(1:nrow(df), round(nrow(df) * .66))
      train <- df[ train.index, ]
      test  <- df[-train.index, ]
    }
    
    if(!is.null(norm))
    {
      if(norm == "scale")
      {
        preProcValues <- preProcess(train, method = c("center", "scale"))
        train <- predict(preProcValues, train)
        
        preProcValues <- preProcess(test, method = c("center", "scale"))
        test <- predict(preProcValues, test)
      }
      else if(norm == "range")
      {
        preProcValues <- preProcess(train, method = "range")
        train <- predict(preProcValues, train)
        
        preProcValues <- preProcess(test, method = "range")
        test <- predict(preProcValues, test)
      }
      else if(norm == "log2")
      {
        train[, sapply(df, is.numeric)] <- log2(train[, sapply(df, is.numeric)] + 1)
        test[, sapply(df, is.numeric)] <- log2(test[, sapply(df, is.numeric)] + 1)
      }
    }
    
    if(!is.null(elimFunc))
    {
      if(elimFunc == "glmnet")
      {
        if(verbose) print("Lasso feature selection on train set")
        
        elim <- glmnet.tuning(train, model=elimFunc, ...,
                              alpha = 1, returnVar=T, stratified=stratified, parallel=parallel, verbose=verbose)
      }
      else if(elimFunc == "rf" | elimFunc == "lr")
      {
        if(verbose) print("RFE on train set")
        
        elim <- recursive.elim(train, funcs=elimFunc, ...,
                               parallel=parallel)
      }
      else stop("Incorrect elimFunc specified")
      
      if(verbose) print(paste(elim$optsize, "variables selected"))
      train <- train[, colnames(train) %in% c("HCC", elim$optVariables)]
      test <- test[, colnames(test) %in% c("HCC", elim$optVariables)]
    }
    else
    {
      elim <- list(optsize = NULL,
                   optVariables = NULL,
                   variables = NULL)
    }
    
    if(search != "none")
      if(verbose) print("Tuning on train set")
    else
      if(verbose) print("Training on train set")
    final_model <- tuneFunc(train, search, model, repeats, folds, stratified, parallel, verbose)
    
    trn_prob <- predict(final_model, train, type="prob")[, 2]; trn_AUC <- calc_AUC(trn_prob, train$HCC)
    test_prob <- predict(final_model, test, type="prob")[, 2]; test_AUC <- calc_AUC(test_prob, test$HCC)
    
    res <- list(Train = trn_AUC,
                Test = test_AUC,
                Tune = final_model$bestTune,
                pred = data.frame(prob = test_prob,
                                  obs = test$HCC),
                optsize = elim$optsize,
                optVariables = elim$optVariables,
                variables = elim$variables,
                results = elim$results)
    
    if(verbose) print(paste(i, round(res$Train, 4), round(res$Test, 4)))
    result <- c(result, list(res))
  }
  
  print(paste(mean(extract_AUC(result)$Train),
              mean(extract_AUC(result)$Test)))
  
  return(result)
}

single_dataset <- function(train, test,
                           search, model,
                           repeats=3, folds=10,
                           stratified=T, norm=NULL, parallel=T, verbose=T, tuneFunc=NULL,
                           elimFunc=NULL, ...)
{
  if(!is.null(norm))
  {
    if(norm == "scale")
    {
      preProcValues <- preProcess(train, method = c("center", "scale"))
      train <- predict(preProcValues, train)
      
      preProcValues <- preProcess(test, method = c("center", "scale"))
      test <- predict(preProcValues, test)
    }
    else if(norm == "range")
    {
      preProcValues <- preProcess(train, method = "range")
      train <- predict(preProcValues, train)
      
      preProcValues <- preProcess(test, method = "range")
      test <- predict(preProcValues, test)
    }
    else if(norm == "log2")
    {
      train[, sapply(df, is.numeric)] <- log2(train[, sapply(df, is.numeric)] + 1)
      test[, sapply(df, is.numeric)] <- log2(test[, sapply(df, is.numeric)] + 1)
    }
  }
  
  if(!is.null(elimFunc))
  {
    if(elimFunc == "glmnet")
    {
      if(verbose) print("Lasso feature selection on train set")
      
      elim <- glmnet.tuning(train, model=elimFunc, ...,
                            alpha = 1, returnVar=T, stratified=stratified, parallel=parallel, verbose=verbose)
    }
    else if(elimFunc == "rf" | elimFunc == "lr")
    {
      if(verbose) print("RFE on train set")
      
      elim <- recursive.elim(train, funcs=elimFunc, ...,
                             parallel=parallel)
    }
    else stop("Incorrect elimFunc specified")
    
    if(verbose) print(paste(elim$optsize, "variables selected"))
    train <- train[, colnames(train) %in% c("HCC", elim$optVariables)]
    test <- test[, colnames(test) %in% c("HCC", elim$optVariables)]
  }
  else
  {
    elim <- list(optsize = NULL,
                 optVariables = NULL,
                 variables = NULL)
  }
  
  if(search != "none")
    if(verbose) print("Tuning on train set")
  else 
    if(verbose) print("Training on train set")
  final_model <- tuneFunc(train, search, model, repeats, folds, stratified, parallel, verbose)
  
  trn_prob <- predict(final_model, train, type="prob")[, 2]; trn_AUC <- calc_AUC(trn_prob, train$HCC)
  test_prob <- predict(final_model, test, type="prob")[, 2]; test_AUC <- calc_AUC(test_prob, test$HCC)
  
  res <- list(Train = trn_AUC,
              Test = test_AUC,
              Tune = final_model$bestTune,
              pred = data.frame(prob = test_prob,
                                obs = test$HCC),
              optsize = elim$optsize,
              optVariables = elim$optVariables,
              variables = elim$variables,
              result = elim$results)
  
  if(verbose) print(paste(round(res$Train, 4), round(res$Test, 4)))
  
  return(res)
}

# Helper functions
estimate_param <- function(model, df)
{
  param <- NULL
  
  if(model == "ranger")
  {
    param <- expand.grid(mtry = round(sqrt(ncol(df))),
                         splitrule = "gini",
                         min.node.size = 1)
  }
  if(model == "svmRadial")
  {
    param <- expand.grid(sigma = kernlab::sigest(HCC ~ ., data=df)[2],
                         C = 1)
  }
  
  return(param)
}

calc_AUC <- function(prob, target)
{
  target <- factor(ifelse(target == "HCC", 1, 0))
  pred_ROCR <- prediction(prob, target)
  auc_ROCR <- performance(pred_ROCR, measure = "auc")
  auc_ROCR <- auc_ROCR@y.values[[1]]
  
  return(auc_ROCR)
}

extract_AUC <- function(result)
{
  Train <- Test <- c()
  
  for(i in 1:length(result))
  {
    Train <- c(Train, result[[i]]$Train)
    Test <- c(Test, result[[i]]$Test)
  }
  
  return(list(Train = Train, Test = Test))
}

mean_AUC <- function(result)
{
  Train <- Test <- c()
  
  for(i in 1:length(result))
  {
    Train <- c(Train, result[[i]]$Train)
    Test <- c(Test, result[[i]]$Test)
  }
  
  return(data.frame(Train = mean(Train), Test = mean(Test)))
}

as.one.df <- function(order, ...)
{
  object <- as.list(substitute(list(...)))[-1L]
  mrn <- is.null(row.names)
  x <- list(...)
  names(x) <- order
  
  x <- lapply(x, extract_AUC)
  
  Train <- Test <- list()
  lapply(1:length(x), function(i){
    Train <<- c(Train, list(data.frame(AUC = x[[i]]$Train,
                                      Group = order[i])))
    
    Test <<- c(Test, list(data.frame(AUC = x[[i]]$Test,
                                     Group = order[i])))
  })
  
  Train <- data.table::rbindlist(Train)
  Test <- data.table::rbindlist(Test)
  
  Train$lGroup <- ifelse(str_detect(Train$Group, "glm"), "glm", ifelse(str_detect(Train$Group, "svm"), "svm", ifelse(str_detect(Train$Group, "rf"), "rf", NA)))
  Test$lGroup <- ifelse(str_detect(Test$Group, "glm"), "glm", ifelse(str_detect(Test$Group, "svm"), "svm", ifelse(str_detect(Test$Group, "rf"), "rf", NA)))
  
  return(list(Train = Train, Test = Test))
}

interpolate_AUC <- function(pred_list, group)
{
  mean_fpr <- seq(0, 1, length=100)
  tpr <- list()
  
  lapply(pred_list, function(x){
    roc <- pROC::roc(response = x$pred$obs, x$pred$prob, positive="HCC")
    se_sp <- pROC::coords(roc, seq(0, 1, length=100), transpose = T)[-1, ]
    
    # tpr <- se_sp[2, ]
    # fpr <- 1-se_sp[1, ]
    
    interp_tpr <- approx(1-se_sp[1, ], se_sp[2, ], mean_fpr)[[2]]
    interp_tpr[1] = 0
    
    tpr <<- c(tpr, list(interp_tpr))
  })
  tpr <- do.call("rbind", tpr)
  mean_tpr <- colMeans(tpr)
  mean_tpr[length(mean_tpr)] <- 1
  
  std_tpr <- apply(tpr, 2, sd)
  tpr_upper <- pmin(mean_tpr + std_tpr, 1)
  tpr_lower <- pmax(mean_tpr - std_tpr, 0)
  
  data <- data.frame(mean_tpr = mean_tpr,
                     mean_fpr = mean_fpr,
                     lower = tpr_lower,
                     upper = tpr_upper,
                     Groups = group)
  
  return(data)
}

plot_sdAUC <- function(df)
{
  df$Groups <- as.factor(df$Groups)
  
  plot <- ggplot(df, aes(mean_fpr, mean_tpr, group=Groups, colour=Groups, fill=Groups)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.1,
                colour = NA) +
    geom_abline(intercept = 0, slope = 1, col="red", size=1) +
    xlab("False Positive Rate (FPR)") + ylab("True Positive Rate (TPR)")
  
  print(plot)
}

# Univariate analysis
var.analysis <- function(res, plot = T)
{
  variable.table <- c()
  score <- list()
  optNum <- c()
  
  lapply(res, function(df){
    variable.table <<- c(variable.table, df$optVariables)
    score <<- c(score, list(df$variables))
    optNum <<- c(optNum, df$optsize)
  })
  score <- do.call(rbind, score)
  score.plot <- aggregate(Overall ~ var, score, mean)
  score.plot <- score.plot[order(score.plot$Overall, decreasing = T), ]
  score.plot <- head(score.plot, median(optNum))
  score.plot <- score[score$var %in% score.plot$var, ]
  
  variable.table <- as.data.frame(table(variable.table) / length(res))
  colnames(variable.table) <- c("variable", "freq")
  variable.table <- variable.table[order(variable.table$freq, decreasing = T), ]
  table.plot <- head(variable.table, median(optNum))
  
  result <- list(variable.table = as.data.frame(variable.table),
                 optVariables = factor(table.plot$variable),
                 score = score,
                 optNum = round(median(optNum)))
  # result$variable.table <- result$variable.table[order(result$variable.table$freq, decreasing = T), ]
  
  if(plot)
  {
    a <- ggplot(table.plot, aes(reorder(variable, freq), freq)) +
      geom_bar(stat='identity') +
      coord_flip(ylim=c(0.25, 1)) +
      xlab("Variables") + ylab("Probability")
    
    # ordering <- with(score.plot, reorder(var, Overall, median, order = TRUE))
    # score.plot$var <- factor(score.plot$var, levels = levels(ordering))
    # 
    # b <- ggplot(score.plot, aes(var, Overall)) +
    #   geom_boxplot() +
    #   coord_flip() +
    #   xlab("") + ylab("Overall Score")
    
    # grid.arrange(a, b, nrow=1)
    return(a)
    # print(b)
  }
  
  return(result)
}

uni.analysis <- function(df, variables, ...)
{
  result <- list()
  for(i in variables)
  {
    print(paste(i, paste()))
    
    uni <- df[, c("HCC", i)]
    glm <- iter_tune(uni, ...)
    
    result <- c(result, list(glm))
  }
  names(result) <- variables
  
  return(result)
}
