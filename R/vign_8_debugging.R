# Vignette 8: Debugging


lib.loc <- NULL
#lib.loc <- "c:/work/wham/old_packages/lab"
library("wham", lib.loc = lib.loc)
wham.dir <- find.package("wham")
asap3 <- read_asap3_dat(file.path(wham.dir,"extdata","ex7_SNEMAYT.dat"))
input <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                              NAA_re = list(sigma="rec+1", cor="iid"))

mod <- fit_wham(input, do.osa = F, do.retro = F)

# using grep to check component nlls for NAs
mod <- fit_wham(input, do.fit = F)
mod$fn()
rep <- mod$rep
nlls <- grep("nll", names(rep), value = TRUE)
nlls
sapply(nlls, function(x) sum(rep[[x]]))

rep$nll_agg_catch
# Year 42 is the problem!

input$data$agg_catch[42,]

input$data$agg_catch_sigma[42,]
input$par$log_catch_sig_scale

rep$pred_stock_CAA[1,1,42,]
rep$waa_catch[1,42,] 
#It's the weight at age 1!

# checking selectivity parameters
# 	check selAA[[]]
# fit_8 from day 2 session 1 has a high selpar estimated.
# fit_0 from day 2 session 2 has a hig selpar


# checking variance parameters/random effects
sort(mod$opt$par)
mod$parList

# age comp and index cv inputs needs to be 
# Dirichlet multinomial input Neff should be an upper bound
# Index CVs from ASAP3 file are often adjusted from the CVs from the bottom trawl surveys 



# internal calculation of reference points needs a different intial value
# round(fit$rep$log_SPR0 + log(0.4) - fit$rep$log_SPR_FXSPR, 4)
# fit$rep$log_SPR0_static + log(0.4) - fit$rep$log_SPR_FXSPR_static
# exp(fit$rep$log_FXSPR)
# input$data$FXSPR_static_init
# nofit$env$data$FXSPR_static_init

# grep("static", names(nofit_0$rep), value = TRUE)
# plot_wham_output(nofit_0, dir.main = tmp.dir)


#no features for estimation in phases 

#have to do this by hand

#using parList from previous fit
#most important for fixed effects