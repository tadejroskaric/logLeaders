### All functions need to be in global environment

set.seed(2024)
# Mixed data (columns 2 and 3 are categorical)
mix_data = data.frame(ggplot2::diamonds[c(1,2,3,5,6,7)])
mix_data = mix_data[sample(nrow(mix_data),100),]
clustering_mix = logLeaders(data=mix_data,noClusters=3,factors=c(2,3),niter=20,BIC=F,leaders=NULL)

# Only categorical data
cat_data = data.frame(ggplot2::diamonds[c(2,3,4)])
cat_data = cat_data[sample(nrow(cat_data),100),]
clustering_cat = logLeaders(data=cat_data,noClusters=3,factors=1:3,niter=20,BIC=F,leaders=NULL)

# Only numeric data
num_data = data.frame(ggplot2::diamonds[c(1,5,6,7)])
num_data = num_data[sample(nrow(num_data),100),]
clustering_num = logLeaders(data=num_data,noClusters=3,factors=NULL,niter=20,BIC=F,leaders=NULL)

############################

# Choosing initial leaders and BIC=T (mixed data)
df = data.frame(ggplot2::diamonds[c(1,2,3,5,6,7)])
df = df[sample(nrow(df),100),]
clustering = logLeaders(data=df,noClusters=3,factors=c(2,3),niter=20,BIC=T,leaders=c(2,33,60))

clustering$output
clustering$BIC

