library(stringr);
library(ggplot2);
library(grid);
library(reshape2);
library(dplyr);
library(plyr);
library(parallel);
library(epitools)
rpar <- function(n, xm, a) {
  v <- runif(n)
  xm / v^(1.0/a)
} ## power-law distribution

if(FALSE) { ## test a bit about the power-law distribution
  data.frame(x=rpar(10000,10,2)) %>%
    ggplot(aes(x)) +
    geom_histogram() +
    scale_x_log10() + 
    scale_y_log10()
}

num_generations <- 10000;
num_genes <- 500; ## only consider the top 100 genes. Other lower expressed genes are of no consequence in terms of translation error
min_expr <- 100000;
powerlaw_param <- 2;
sum(rpar(num_genes,min_expr,powerlaw_param)); ## this should be at around the 1e8 level (total number of proteins in Ecoli at about 1e6, yeast at about 1e8, mammalian at about 1e10, according to Bionumbers)
pop_size <- 1000;
Lnt <- 1000; ## nucleotide length
Lp <- 300; ## protein length
e <- 100;##epistasis
shape_est_TR <- 5.86
scale_est_TR <- 0.173
shape_est_TL <- 5.44
scale_est_TL <- 0.197
min_correct_frac <- 0.5

ancPop <- data.frame(
  ind = rep(1:pop_size,each=num_genes),
  beta1=rep(runif(num_genes, 10^-5, 10^-5),pop_size),
  beta2=rep(runif(num_genes, 10^-4, 10^-4),pop_size),
  N = rep(rpar(num_genes,min_expr,powerlaw_param),pop_size) ## using power-law distribution for expression
);
pop <- ancPop;

dfParamSearch <- expand.grid(
  mut_rate = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01),
  param_c =  c(1e-5,1e-6,1e-7,1e-8)
);

trackFitness <-list();
trackCor <- list();
relativeFitness <- list();
trackFitness_data <- list()

mclapply(
  c(1:nrow(dfParamSearch)),
  mc.cores = 30,
  mc.preschedule = F,
  function(thisParam){
    mut_rate = dfParamSearch[thisParam,"mut_rate"];
    param_c = dfParamSearch[thisParam,"param_c"];
    filename = paste0("01.sim_",mut_rate,"_",param_c,".RData");

    trackFitness[[filename]] <- rep(NA,num_generations);
    trackCor[[filename]] <- rep(NA,num_generations);

    for(thisGen in c(1:num_generations)) {
      curPop <- pop %>%  ## current population
        select(beta1,beta2,N) %>%
        mutate(ind = rep(1:pop_size,each=num_genes)); ## rename the individuals
      nRow <- nrow(curPop);
      mut_effect1 <- rgamma(nRow, shape=shape_est_TR, scale=scale_est_TR);
      curPop$beta1 <- curPop$beta1 * ifelse(runif(nRow) < mut_rate, mut_effect1,1);
      mut_effect2 <- rgamma(nRow, shape=shape_est_TL, scale=scale_est_TL);
      curPop$beta2 <- curPop$beta2 * ifelse(runif(nRow) < mut_rate, mut_effect2,1);
      curPop$err1 <- 1 - (1-curPop$beta1)^Lnt; ## prob of at least one transcription error
      curPop$err2 <- 1 - (1-curPop$beta2)^Lp; ## prob of at least one translation error
      curPop$totalCorrect <- (1-curPop$err1)*(1-curPop$err2);
      curPop$totalErr <- 1 - curPop$totalCorrect; ## prob of at least one error
      curPop$bothErr <- curPop$err1 * curPop$err2;
      ind2fitness <- curPop %>%
        group_by(ind) %>%
        dplyr::summarise(
          totalErrCnt = sum(N * err1*(1-err2)+N*err2*(1-err1)+N*e*err1*err2),  ## toxicity from error molecules
          enoughCorrect = prod(totalCorrect>min_correct_frac)## sufficient correct molecules
        ) %>%
        mutate(fitness = enoughCorrect *exp(-param_c * totalErrCnt));
      prolif <- sample.int(pop_size,prob=ind2fitness$fitness,replace=T);
      prolif_row <- lapply(prolif,function(x){seq((x-1) * num_genes + 1,x*num_genes,by=1)}) %>% unlist();
      pop <- curPop[prolif_row,]; ## proliferation
      
      trackFitness[thisGen] <- mean(ind2fitness$fitness);
      trackCor[thisGen] <- cor.test(pop$beta1[1:num_genes],pop$beta2[1:num_genes],method="pearson")$estimate;
      if(thisGen %% 200 == 0) {
        save(file=filename,list=c("filename","pop","trackFitness","trackCor","mut_rate","param_c","thisGen"));
        cat(paste0(filename," Generation ",thisGen," . Average fitness ",mean(ind2fitness$fitness)," correlation at ",trackCor[thisGen]," Prob(cor<0): ", mean(trackCor<0,na.rm=T),"\n"));
      }
      
    
    }
  }
)

