# Vignette 5: Multiple regions, stocks, and movement

lib.loc <- NULL
#lib.loc <- "c:/work/wham/old_packages/lab"
library("wham", lib.loc = lib.loc)

path_to_examples <- system.file("extdata", package="wham")

two_stocks_asap <- read_asap3_dat(file.path(path_to_examples,c("ex1_SNEMAYT.dat","ex1_SNEMAYT.dat")))
ini_2_stocks <- prepare_wham_input(two_stocks_asap)
fit_2_stocks <- fit_wham(ini_2_stocks, do.osa = FALSE, do.retro= FALSE, do.sdrep = FALSE)
fit_2_stocks_RDS <- file.path(res_dir,"fit_2_stocks.RDS")
saveRDS(fit_2_stocks, fit_2_stocks_RDS)

fit_2_stocks$rep$SSB


#########################################################
#different stocks
diff_stocks_asap <- read_asap3_dat(file.path(getwd(),"data",c("north.dat","south.dat")))
selectivity <- list(model = rep(c("logistic", "age-specific"),c(8,4)), n_selblocks = 12,
	fix_pars = c(rep(list(NULL),8), list(2:8,3:8,3:8,2:8)),
	initial_pars = c(rep(list(c(2,0.2)),8),list(rep(c(0.5,1),c(1,7)), rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(1,7)))))
diff_stocks_input <- prepare_wham_input(diff_stocks_asap, selectivity = selectivity)
fit_diff_stocks <- fit_wham(diff_stocks_input, do.osa = FALSE, do.retro= FALSE, do.sdrep = FALSE)



#########################################################
#not using asap3 dat file
input <- diff_stocks_input

catch <- input$data$agg_catch
catch_cv <- cbind(sqrt(exp(input$data$agg_catch_sigma^2) - 1))
catch_Neff <- cbind(input$data$catch_Neff)
catch_paa <- input$data$catch_paa
use_catch_paa <- input$data$use_catch_paa
selblock_pointer_fleets <- input$data$selblock_pointer_fleets

index <- input$data$agg_indices
index_cv <- cbind(sqrt(exp(input$data$agg_index_sigma^2) - 1))
index_Neff <- cbind(input$data$index_Neff)
index_fracyr <- input$data$fracyr_indices
units_indices <- input$data$units_indices
units_index_paa <- input$data$units_index_paa
use_indices <- input$data$use_indices
use_index_paa <- input$data$use_index_paa
index_paa <- input$data$index_paa
selblock_pointer_indices <- input$data$selblock_pointer_indices

waa <- input$data$waa
mat <- input$data$mature
waa_pointer_indices <- input$data$waa_pointer_indices
waa_pointer_fleets <- input$data$waa_pointer_fleets
waa_pointer_ssb <- input$data$waa_pointer_ssb
fracyr_ssb <-  input$data$fracyr_SSB

basic_info <- list(
	n_stocks = 2L,
	n_regions = 2L,
	ages = 1:input$data$n_ages,
	n_seasons = 1L,
	n_fleets = input$data$n_fleets,
	fracyr_SSB = fracyr_ssb,
	maturity = mat,
	years = as.integer(input$years),
	waa = waa,
	waa_pointer_ssb = waa_pointer_ssb,
	spawn_regions = 1:2,
	spawn_seasons = c(1,1)
)

# x <- input$par$logit_selpars
# x[] <- input$map$logit_selpars

catch_info <- list(
	n_fleets = NCOL(catch),
	agg_catch = catch,
	agg_catch_cv = catch_cv,
	catch_paa = catch_paa,
	use_catch_paa = use_catch_paa,
	catch_Neff = catch_Neff,
	selblock_pointer_fleets = selblock_pointer_fleets,
	waa_pointer_fleets = waa_pointer_fleets,
	fleet_regions = rep(c(1:2),c(2,2))
)

index_info <- list(
	n_indices = NCOL(index),
	agg_indices = index,
	units_indices = units_indices,
	units_index_paa = units_index_paa,
	agg_index_cv = index_cv,
	fracyr_indices = index_fracyr,
	use_indices = use_indices,
	use_index_paa = use_index_paa,
	index_paa = index_paa,
	index_Neff = index_Neff,
	selblock_pointer_indices = selblock_pointer_indices,
	waa_pointer_indices = waa_pointer_indices,
	index_regions = rep(c(1:2),c(2,2))
)

MAA <- exp(input$par$M_re)
for(i in 1:2) MAA[i,i,,] <- MAA[i,i,,]*exp(matrix(input$par$Mpars[i,i,], length(basic_info$years), length(basic_info$ages), byrow = TRUE))

M_in <- list(initial_MAA = MAA)

F_in <- list(F = matrix(0.3, length(basic_info$years), catch_info$n_fleets))
q_in <- list(initial_q = rep(1e-6, index_info$n_indices))


#all at once 
input_all <- prepare_wham_input(
	basic_info = basic_info,
	NAA_re = list(sigma = "rec+1"),
	selectivity = selectivity, 
	catch_info = catch_info, 
	index_info = index_info, 
	M = M_in, 
	F = F_in, 
	catchability = q_in)



#nlls <- grep("nll", names(nofit_all$rep))
#sapply(nofit_all$rep[nlls], sum)

#piece by piece
input_seq <- prepare_wham_input(basic_info = basic_info)
input_seq <- set_NAA(input_seq, list(sigma = "rec+1"))
input_seq <- set_M(input_seq, M = M_in)
input_seq <- set_catch(input_seq, catch_info = catch_info)
input_seq <- set_F(input_seq, F_in)
input_seq <- set_indices(input_seq, index_info = index_info)
input_seq <- set_q(input_seq, q_in)
input_seq <- set_selectivity(input_seq, selectivity = selectivity)

input_asap <- prepare_wham_input(diff_stocks_asap, selectivity = selectivity, NAA_re = list(sigma = "rec+1"))

#compare 
nofit_all <- fit_wham(input_all, do.fit = FALSE)
nofit_seq <- fit_wham(input_seq, do.fit = FALSE)
nofit_asap <- fit_wham(input_asap, do.fit = FALSE)
cbind(nofit_all$par, nofit_asap$par)

nofit_seq$fn() - nofit_all$fn() #equal
length(nofit_seq$par) - length(nofit_all$par) #equal

length(nofit_seq$par) - length(nofit_asap$par) #equal

nofit_all$fn() - nofit_asap$fn() #different

nofit_asap$par-nofit_all$par #different

nofit_all$fn() - nofit_asap$fn(nofit_all$par) #equal

fit_seq <- fit_wham(input_seq, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
fit_all <- fit_wham(input_all, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
fit_asap <- fit_wham(input_asap, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
fit_seq$opt$obj - fit_all$opt$obj # equal
fit_asap$opt$obj - fit_all$opt$obj # equal

fit_all$fn(fit_asap$opt$par)
res_dir <- file.path(getwd(),"temp")


saveRDS(fit_asap, file.path(res_dir,"fit_2_stocks_asap.RDS"))
saveRDS(fit_all, file.path(res_dir,"fit_2_stocks_all.RDS"))
saveRDS(fit_seq, file.path(res_dir,"fit_2_stocks_seq.RDS"))



fit_asap <- do_reference_points(fit_asap, do.sdrep = TRUE)
fit_asap$peels <- retro(fit_asap)
fit_asap <- make_osa_residuals(fit_asap)

tmp.dir <- tempdir(check=TRUE)
plot_wham_output(fit_asap, dir.main = tmp.dir)

fit_RDS <- file.path(res_dir,"fit_2_stocks.RDS")
saveRDS(fit_asap, fit_RDS)

x <- jitter_wham(fit_RDS = fit_RDS, n_jitter = 10, res_dir = res_dir, do_parallel = FALSE)
sapply(x[[1]], function(y) y$obj) #nlls


