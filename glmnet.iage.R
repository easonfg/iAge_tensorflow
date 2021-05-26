library(dplyr)
library(tibble)
library(glmnet)
iage.data = read.csv('../data/combine_031517.csv')
dim(iage.data)
iage.data = iage.data[!is.na(iage.data$BCELLS),]
iage.data = subset(iage.data, select = -c(CHEX1, CHEX2, CHEX3, CHEX4))
dim(iage.data)
sum(iage.data$GENDER == 0)

## just use scale
#make_standardize = function(x){
#    #n,p = dim(x)
#    sizes = dim(x)
#    X = x
#    mus = c()
#    stds = c()
#    for (i in 1:sizes[2]){
#        x = X[,i]
#        mu = mean(x)
#        st = sd(x)
#        x = x - mu
#        x = x / st
#        X[,i] = x
#        mus = c(mus, mu)
#        stds = c(stds, st)
#    }
#    # browser()
#    return(list('X' = X,'mus' = mus,'stds' = stds))
#}


remove_cyto_outliers = function(lD_cyto){
  X = lD_cyto
  standardized.data = scale(X)
  upper = attr(standardized.data, 'scaled:center')+attr(standardized.data, 'scaled:scale')*3
  lower = attr(standardized.data, 'scaled:center')-attr(standardized.data, 'scaled:scale')*3
  for (i in 1:dim(X)[1]){
      for (j in 1:dim(X)[2]){
          if (X[i,j]>upper[j]){ X[i,j] = upper[j]}
          if (X[i,j]<lower[j]){ X[i,j] = lower[j]}
      }
  }
  Xs = scale(log(X+0.00001))
  return(Xs)
}

# cytof.data = scale(iage.data[,9:ncol(iage.data)])
# cytof.data = scale(iage.data[,grep('CD40L', colnames(iage.data)):ncol(iage.data)])
# attr(cytof.data, 'scaled:scale')

# make_standardize(iage.data[,grep('CD40L', colnames(iage.data)):ncol(iage.data)])
rm.scaled = remove_cyto_outliers(iage.data[,grep('CD40L', colnames(iage.data)):ncol(iage.data)])
rm.scaled
# sum(rm.scaled - cytof.data)
#### LASSO ####
res.coef = c()
max.auc = c()
res.coef.values = c()
for (all.rep in 1:10){
  for (j in 1:10){
    print(paste('I: ', all.rep, 'J: ', j))
    md3cv <- cv.glmnet(as.matrix(rm.scaled),
                       # iage.data$AGE)
                       scale(iage.data$AGE))
    # plot(md3cv)
    
    return_features = function(coeff){
      top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
      top_features = top_features[order(top_features$coefficient),]
      # print(top_features)
      return(top_features)
    }
    
    tvec = return_features(coef(md3cv, s = 'lambda.min'))$coefficient
    names(tvec) = return_features(coef(md3cv, s = 'lambda.min'))$name
    
    res.coef.values = c(res.coef.values, tvec)
    
    res.coef = c(res.coef, return_features(coef(md3cv, s = 'lambda.min'))$name)
    max.auc = c(max.auc, max(md3cv$cvm))
    # print(res.coef)
    # print(max.auc)
    # res.coef$name
  }
  print(mean(max.auc))
  print(sd(max.auc))
  # jpeg(paste('pred_within_subtypes_res/Moderate/lasso/Moderate_binary_glmnet_example.jpeg', sep = ''),
  #      units="in", width=10, height=10, res=500)
  # plot(md3cv)
  # dev.off()
}

df.res.coef.values = data.frame(res.coef.values)
df.res.coef.values$names = names(res.coef.values)


## count from the 100 reps
genes.count = table(df.res.coef.values$names)
df.genes.count = data.frame(genes.count)
df.genes.count = df.genes.count %>% column_to_rownames(var = 'Var1')
hist(genes.count[order(genes.count, decreasing = T)])
(genes.count[order(genes.count, decreasing = T)])

coef.means = df.res.coef.values %>% group_by(names) %>%
  summarize(coef_means = mean(res.coef.values, na.rm = TRUE))

## weighted coef base on number of appearances in the 100 reps
coef.means = coef.means %>% column_to_rownames(var = "names") 
coef.means
# coef.means * 
weighted.coef = coef.means * df.genes.count[rownames(coef.means),]
weighted.coef = weighted.coef[-1,, drop = FALSE]
head(weighted.coef)
weighted.coef = data.frame(scale(weighted.coef))
weighted.coef = weighted.coef[order(weighted.coef, decreasing = T),,drop = FALSE ]
weighted.coef

# write.csv((genes.count[order(genes.count, decreasing = T)])[-1], 'res/all_proteins/binary_ordered_genes_count.csv')
# write.csv(weighted.coef, 'res/all_proteins/binary_weighted_coef.csv')
# write.csv(coef.means, 'res/all_proteins/binary_mean_coef.csv')

graph.data = coef.means
graph.data$freq = genes.count[rownames(graph.data)]
graph.data = graph.data[-1,]

graph.data$names = sapply(strsplit(rownames(graph.data),'\\.'), function(x) x[1])

library(ggrepel)
# print(
ggplot(graph.data, aes(coef_means, freq, label = names)) +
  # ggplot(graph.data, aes(freq, coef_means)) +
  geom_point(aes(color = abs(coef_means))) +
  # geom_text(aes(label=names),size = 3, hjust=-0.5, vjust=0.5) +
  # geom_text(aes(label=ifelse(freq>75,as.character(names),'')),size = 3, hjust=-0.5, vjust=0.5) +
  ## best alpha
  # geom_text_repel(data          = subset(graph.data, freq > 75 | abs(coef_means)>0.05)) +
  ## lasso
  geom_text_repel(data          = subset(graph.data, freq > 90 | abs(coef_means)>0.025), max.overlaps = 1000) +
  # geom_text(position=position_jitter(width=1,height=1)) +
  # nudge_y       = 32 - subset(graph.data, freq > 25)$freq,
  # size          = 4,
  # box.padding   = 1.5,
  # point.padding = 0.5,
  # force         = 100,
  # segment.size  = 0.2,
  # segment.color = "grey50",
  # direction     = "x") +
  xlab("Average Coefficients") + 
  ylab("Frequencies") +
  labs(color='Absolute Average\ 
       Coefficients') + ## legend title
  ggtitle('Binary Prediction (Severe A)')+
  theme_bw()
# )

graph.data[order(graph.data$coef_means),]
graph.data[order(abs(graph.data$coef_means)),]

## predict age
as.matrix(rm.scaled) %*% as.vector(coef.means[-1])
dim(as.vector(coef.means[-1]))
as.matrix(rm.scaled) %>%dim()
rm.scaled%>%dim()
rm.scaled %*% rep(3,50)
pred_age = rm.scaled[, rownames(coef.means[-1,,drop = F])] %*% coef.means[-1,]
library(Metrics)
rmse(iage.data$AGE, pred_age)
