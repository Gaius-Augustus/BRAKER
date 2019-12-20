setwd("~/")
dat = read.table("prothint_eval.lst", header = T) # result of analyse_prothint.pl

# normalize mult by median
dat$mult_norm = dat$mult/median(dat$mult)

# ratio TRUE/ALL
sum(dat$label==TRUE)/length(dat$label)
# 85%

# split into training and test set
entry_index = 1:length(dat$mult)
testset_index = sample(entry_index, 500)
testset = dat[testset_index,]
train_index = setdiff(entry_index, testset_index)
trainset = dat[train_index,]

glm.fit <- glm(label ~ mult_norm + al_score, data = trainset, family=binomial(link='logit'))
summary(glm.fit)
#            Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.00529    0.04935  -81.16   <2e-16 ***
# mult_norm    4.73909    0.04662  101.66   <2e-16 ***
# al_score     9.09026    0.14741   61.67   <2e-16 ***


fitted.results <- predict(glm.fit ,newdata=testset, type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != testset$label)
print(paste('Accuracy',1-misClasificError))
# "Accuracy 0.93"

library(ROCR)
p <- predict(glm.fit, newdata=testset, type="response")
pr <- prediction(p, testset$label)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
# 0.94

# apply assumed model to testset

manual_prediction = -4.00529 + 4.73909 * testset$mult_norm + 9.09026 * testset$al_score
predicted_label = rep(FALSE, length(manual_prediction))
predicted_label[manual_prediction >= 0.85] = TRUE # logit von 0.7 is 0.85, statt 0.5
misClasificError <- mean(predicted_label != testset$label)
print(paste('Accuracy',1-misClasificError))
# 0.93
# count number of correct examples in testset
sum(testset$label==TRUE)/(sum(testset$label==TRUE) + sum(testset$label==FALSE))
# 0.802



# repeat logistic regression with Arabidopsis species excluded and check whether that makes a difference
setwd("~/")
dat = read.table("prothint_eval.lst", header = T)

# normalize mult by median
dat$mult_norm = dat$mult/median(dat$mult)

# split into training and test set
entry_index = 1:length(dat$mult)
testset_index = sample(entry_index, 500)
testset = dat[testset_index,]
train_index = setdiff(entry_index, testset_index)
trainset = dat[train_index,]
manual_prediction = -4.00529 + 4.73909 * testset$mult_norm + 9.09026 * testset$al_score
predicted_label = rep(FALSE, length(manual_prediction))
predicted_label[manual_prediction >= 0.85] = TRUE # logit von 0.7 is 0.85, statt 0.5
misClasificError <- mean(predicted_label != testset$label)
print(paste('Accuracy',1-misClasificError))
# 95% - > ok. building the model on this data set does not change much.