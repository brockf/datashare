library(data.table)
library(dplyr)
library(lme4)
library(beepr)
library(ggplot2)
library(reshape2)

alpha <- .05
sims <- 50000

sim_data <- data.table(
                  simulation = 1:sims,
                  alpha = 0,
                  beta = 0,
                  sample_sizes = 0,
                  num_trials = 0,
                  condition_1 = 0,
                  condition_2 = 0,
                  condition_1_sd = 0,
                  condition_2_sd = 0,
                  pop_sd = 0,
                  sig_logit = 0,
                  sig_elogit = 0,
                  sig_arcsin = 0,
                  sig_prop = 0
              )

for (i in 1:sims) {
  # sample an alpha and beta for our population
  # samples will be sampled from same population
  beta_alpha <- runif(1,1,1000)
  beta_beta <- runif(1,1,1000)
  
  # uniformly sample sample sizes between 10-100 participants
  sample_sizes <- round(runif(1, 10, 100),0)
 
  # generate study data
  # 10 trials per participant
  num_trials <- round(runif(1, 2, 100),0)
  
  # sample subject means from beta distributions
  condition_1_subjects <- rbeta(sample_sizes, beta_alpha, beta_beta)
  condition_2_subjects <- rbeta(sample_sizes, beta_alpha, beta_beta)
  
  # set to zero to trigger sampling below
  sample_1 <- 0
  sample_2 <- 0
  
  # check samples for no variance (which will cause issues, so we re-sample)
  while (length(sample_1) == 1 | (var(sample_1) == 0 | var(sample_2) == 0)) {
    sample_1 <- rbinom(num_trials*sample_sizes, 1, rep(condition_1_subjects,each=num_trials))
    sample_2 <- rbinom(num_trials*sample_sizes, 1, rep(condition_2_subjects,each=num_trials))
  }
  
  study_data <- data.frame(
                      subject = rep(paste0('subject',1:(sample_sizes*2)), each=num_trials),
                      condition = rep(c('condition1','condition2'), each=(sample_sizes*num_trials)),
                      trial = rep(1:num_trials, times=(sample_sizes*2)),
                      response = c(sample_1,sample_2)
                  )
  
  study_data_agg <- study_data %>%
                    group_by(subject, condition) %>%
                    summarise(y=sum(response), N=length(response)) %>%
                    ungroup() %>%
                    mutate(
                      Prop = y / N,
                      elog = log( (y + .5) / (N - y + .5) ), # empirical logit
                      Arcsin = asin(sqrt(Prop)) # arcsin-sqrt 
                    )
  
  condition_1 <- with(subset(study_data_agg, condition == 'condition1'), mean(Prop))
  condition_2 <- with(subset(study_data_agg, condition == 'condition2'), mean(Prop))
  condition_1_sd <- with(subset(study_data_agg, condition == 'condition1'), sd(Prop))
  condition_2_sd <- with(subset(study_data_agg, condition == 'condition2'), sd(Prop))
  pop_sd <- sd(study_data_agg$Prop)
  
  # logistic regression
  model <- glmer(response ~ condition + (1 | subject), data=study_data, family="binomial")
  zs <- fixef(model)/sqrt(diag(vcov(model)))
  sig_logit <- ifelse(2*pnorm(abs(zs),lower.tail=FALSE)[2] <= alpha, TRUE, FALSE)
    
  # empirical logit
  sig_elogit <- ifelse(t.test(elog ~ condition, data=study_data_agg, var.equal=T)$p.value <= alpha, TRUE, FALSE)
  
  # arcsin-root
  sig_arcsin <- ifelse(t.test(Arcsin ~ condition, data=study_data_agg, var.equal=T)$p.value <= alpha, TRUE, FALSE)
  
  # proportions
  sig_prop <- ifelse(t.test(Prop ~ condition, data=study_data_agg, var.equal=T)$p.value <= alpha, TRUE, FALSE)
  
  # enter into data
  data <- c(
    i,
    beta_alpha,
    beta_beta,
    sample_sizes,
    num_trials,
    condition_1,
    condition_2,
    condition_1_sd,
    condition_2_sd,
    pop_sd,
    sig_logit,
    sig_elogit,
    sig_arcsin,
    sig_prop
  )
  data <- as.vector(data)
  
  sim_data[i, names(sim_data) := as.list(data)]
}

# play the Mario sound when it's done
beep(8)

# calculate grand_mean for each sim
sim_data$grand_mean <- (sim_data$condition_1 + sim_data$condition_2) / 2

# calculate diff for each sim
sim_data$diff <- abs(sim_data$condition_1 - sim_data$condition_2)

sim_data <- data.frame(sim_data)

# overall Type I error rates
colMeans(sim_data[, c('sig_logit','sig_elogit','sig_arcsin','sig_prop')])

######### Analysis 1: type 1 error rates by sample size

# melt and visualize
sim_data_melted <- melt(sim_data, id.vars=c('simulation','sample_sizes','num_trials','grand_mean'), measure.vars = c('sig_logit','sig_elogit','sig_arcsin','sig_prop'), variable.name="test", value.name="sig")

# aggregate by sample sizes, and bin by N=10
agg_by_sample <- sim_data_melted %>%
                 filter(sample_sizes < 99) %>%
                 mutate(BinSample = sample_sizes %/% 5) %>%
                 group_by(BinSample,test) %>%
                 summarise(mean_sig = mean(sig), N=length(sig)) %>%
                 ungroup()

p <- ggplot(agg_by_sample, aes(x=BinSample, y=mean_sig, color=test)) +
                                stat_smooth(method="loess", se=FALSE, size=1) +
                                geom_hline(yint=.05, linetype="dashed", alpha=.5) +
                                scale_color_discrete(name="", labels=c('Logistic Regression','Empirical Logit','Arcsin-Root','Raw Proportions')) +
                                scale_y_continuous(name="Type I Error Rate") +
                                scale_x_continuous(name="Sample Size", labels=seq(0,100,by=25)) +
                                theme(legend.position="top")

png(filename = "figures/type1-by-N.png",
    width = 3000, height = 2500, units = "px", pointsize = 16,
    bg = "white", res = 400)
p
dev.off()
p

######### Analysis 2: type 1 error rates by grand mean

# aggregate by mean, and bin by .10
agg_by_mean <- sim_data_melted %>%
                  filter(grand_mean > .1 & grand_mean < .9) %>%
                  mutate(BinMean = round(grand_mean,2)) %>%
                  group_by(BinMean,test) %>%
                  summarise(mean_sig = mean(sig), N=length(sig)) %>%
                  ungroup()

p <- ggplot(agg_by_mean, aes(x=BinMean, y=mean_sig, color=test)) +
                                stat_smooth(method="loess", se=FALSE, size=1) +
                                geom_hline(yint=.05, linetype="dashed", alpha=.5) +
                                scale_color_discrete(name="", labels=c('Logistic Regression','Empirical Logit','Arcsin-Root','Raw Proportions')) +
                                scale_y_continuous(name="Type I Error Rate") +
                                scale_x_continuous(name="Grand Mean", breaks=seq(0,1.0,by=.1)) +
                                theme(legend.position="top")

png(filename = "figures/type1-by-grand-mean.png",
    width = 3000, height = 2500, units = "px", pointsize = 16,
    bg = "white", res = 400)
p
dev.off()
p

######### Analysis 3: when do these tests not agree? (should be worse at the tails)

sim_data <- sim_data %>%
  mutate(SumSig = sig_logit + sig_elogit + sig_arcsin + sig_prop,
         Agree = ifelse(SumSig == 4 | SumSig == 0, TRUE, FALSE))

agreement_agg_by_mean <- sim_data %>%
                    filter(grand_mean > .1 & grand_mean < .9) %>%
                    mutate(BinMean = round(grand_mean,2)) %>%
                    group_by(BinMean) %>%
                    summarise(mean_agreement = mean(Agree), N=length(Agree)) %>%
                    ungroup()

p <- ggplot(agreement_agg_by_mean, aes(x=BinMean, y=mean_agreement)) +
                  stat_smooth(method="loess", se=FALSE, size=1) +
                  scale_y_continuous(name="Agreement between Tests") +
                  scale_x_continuous(name="Grand Mean", breaks=seq(0,1.0,by=.10)) +
                  theme(legend.position="top")

png(filename = "figures/null-agreement-by-grand-mean.png",
    width = 3000, height = 2500, units = "px", pointsize = 16,
    bg = "white", res = 400)
p
dev.off()
p

# correlate significances
cor(sim_data[, c('sig_logit','sig_elogit','sig_arcsin','sig_prop')])