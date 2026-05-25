library(psych)
options(scipen = 999)

latent_int_5rep <- read.table("~/GitHub/latent_interaction/pre_review_5reps/latent_int_5rep.dat", quote="\"", comment.char="")

colnames(latent_int_5rep) <- c("categories", "group_prob", "rsq_prod", "N", "loading", 
                                "n_items", "model_type", "param.id", "est", "true", 
                                "bias", "rel.bias", "squared_bias","pval", "sig", "ci.cov")


# Relative Bias
relbias_mean <- describeBy(rel.bias ~ categories + group_prob + rsq_prod + N + loading + n_items + model_type + param.id, 
                           data = latent_int_5rep, mat = T)
relbias_mean <- relbias_mean[complete.cases(relbias_mean),]
relbias_mean <- cbind.data.frame(relbias_mean[,2:9],relbias_mean[,12:13],relbias_mean[,17:18])
relbias_mean <- relbias_mean[!is.infinite(relbias_mean$mean), ]
colnames(relbias_mean) <- c("categories", "group_prob", "rsq_prod", "N", "loading", "n_items", 
                            "model_type", "param.id", "mean","sd","min","max")
relbias_mean


# Significance
sig_mean <- describeBy(sig ~ categories + group_prob + rsq_prod + N + loading + n_items + model_type + param.id, 
                       data = latent_int_5rep, mat = T)
sig_mean <- sig_mean[complete.cases(sig_mean),]
sig_mean <- cbind.data.frame(sig_mean[,2:9],sig_mean[,12:13],sig_mean[,17:18])
sig_mean <- sig_mean[!is.infinite(sig_mean$mean), ]
colnames(sig_mean) <- c("categories", "group_prob", "rsq_prod", "N", "loading", "n_items", 
                        "model_type", "param.id", "mean","sd","min","max")
sig_mean


# Confidence Interval Coverage
cicov_mean <- describeBy(ci.cov ~ categories + group_prob + rsq_prod + N + loading + n_items + model_type + param.id, 
                         data = latent_int_5rep, mat = T)
cicov_mean <- cicov_mean[complete.cases(cicov_mean),]
cicov_mean <- cbind.data.frame(cicov_mean[,2:9],cicov_mean[,12:13],cicov_mean[,17:18])
cicov_mean <- cicov_mean[!is.infinite(cicov_mean$mean), ]
colnames(cicov_mean) <- c("categories", "group_prob", "rsq_prod", "N", "loading", "n_items", 
                          "model_type", "param.id", "mean","sd","min","max")
cicov_mean


est.true <- aggregate(latent_int_5rep[,c("est","true")],
                      by = list(categories = latent_int_5rep$categories,
                                group_prob = latent_int_5rep$group_prob,
                                rsq_prod = latent_int_5rep$rsq_prod,
                                N = latent_int_5rep$N,
                                loading = latent_int_5rep$loading,
                                n_items = latent_int_5rep$n_items,
                                model_type = latent_int_5rep$model_type,
                                param.id = latent_int_5rep$param.id),
                      FUN = mean)



cicov_mean <- cicov_mean[cicov_mean$param.id != "beta_0",]
relbias_mean <- relbias_mean[relbias_mean$param.id != "beta_0",]
sig_mean <- sig_mean[sig_mean$param.id != "beta_0",]

est.true <- est.true[est.true$param.id != "beta_0" & est.true$param.id != "beta_G2" & est.true$param.id != "beta_G3",]
subset <- est.true[est.true$param.id == "beta_XG2",]
