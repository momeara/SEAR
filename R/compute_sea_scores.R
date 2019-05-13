# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(fdrtool)

#' Compute target vs. target scores for all targets in the library
#' This sets up the sge cluster job that is then submitted from the command line:
#'
#' run_base <- <runs_dir>/<library_id>_vs_<library_id>
#'
#'   cd <run_base>/logs
#'   qsub -t 1-<n_targets> <path/to/sea-wrapper.sh> <run_base> <library_fname> <fp_type>
#'
#' then run SEAR::collect_compute_sea_network(...) to generate scores file
#'
#' @param library_fname .sea file
#'    generated with SEAR::pack.library(...)
#'    and then the background fit with plot_fits.library(...) and then set_fit.library(...)
#' @return run_base
#' @export
submit_compute_sea_network <- function(
	library_fname,
	library_id,
	fp_type,
	runs_dir=".",
	config_files=NULL,
	target_range="all",
	verbose=TRUE,
	dry_run=FALSE
){

	if(!dir.exists(runs_dir)){
		if(verbose){
			cat("Creating '", runs_dir, "' ...\n", sep="")
		}
		dir.create(runs_dir)
	}

	sea_library <- unpack.library(library_fname, verbose=T)
	sea_library <- sea_library$targets %>%
		left_join(
			sea_library$molecules,
			by=c("compound"))

	run_base <- paste0(
		runs_dir, "/",
		library_id, "_vs_", library_id)

	if(verbose){
		cat("Creating run_base: '", run_base, "' ...\n", sep="")
	}
	unlink(run_base, recursive=T, force=T)
	dir.create(run_base)

	if(verbose){
		cat(
			"Computing SEA target vs. target for all targets in library:\n",
			"\tlibrary_id: ", library_id, "\n",
			"\tfp_type: ", fp_type, "\n",
			"\tconfig_files: ", paste(names(config_files), config_files, collapse=", ", sep=": "), "\n",
			"\trun_base: ", run_base, "\n",
			"\tsea_library_fname: ", library_fname, "\n",
			"\tdry_run: ", ifelse(dry_run, "yes", "no"), "\n",
			sep="")
	}

	n_targets <- sea_library %>% distinct(target) %>% nrow
	if(target_range == 'all'){
		target_range <- paste0("1-", n_targets)
	}

	if(verbose){
		cat(
			"\tn_targets: ", n_targets, "\n",
			"\ttarget_range: ", target_range, "\n",
			sep="")
	}

	#### prepare inputs
	if(verbose){
		cat("Preparing run directories and data...\n")
	}

	if(!is.null(config_files)){
		do.call(SEAR::write.config_files, c(dir=run_base, config_files))
	}

	outputs_path <- paste0(run_base, "/outputs")
	unlink(outputs_path, recursive=T, force=T)
	dir.create(outputs_path)

	logs_path <-paste0(run_base, "/logs")
	unlink(logs_path, recursive=T, force=T)
	dir.create(logs_path)

	inputs_path <- paste0(run_base, "/inputs")
	unlink(inputs_path, recursive=T, force=T)
	dir.create(inputs_path, recursive=T)
	sea_library %>%
		plyr::d_ply(c("target"), function(df){
			target <- df$target[1]
			target_fname <- paste0(inputs_path, "/", target, ".csv")
			df %>%
				dplyr::select(compound, smiles) %>%
				write.table(
					target_fname,
					quote=F,
					sep=",",
					row.names=F,
					col.names=T)
		})

	sea_wrapper_fname <- paste0(path.package("SEAR"), "/scripts/sea-wrapper.sh")

	cmd <- paste0(
		sea_wrapper_fname, " ",
		run_base, " ",
		library_fname, " ",
		fp_type)

	script <- paste0("
	cd ", logs_path, "
	qsub -t ", target_range, " ", cmd, "\n")

  if(verbose){
		cat("script:\n", script, sep="")
	}

	if(!dry_run){
		system(paste0(script, "\n"))
		if(verbose){
			cat("Run submitted, to check if it is done do 'qstat'\n\n", sep="")
		}
	} else {
		if(verbose){
			cat("To submit run copy and paste in shell, then to check if it is done do 'qstat'\n\n", sep="")
		}
	}
	run_base
}


#' Collect target vs. target scores computed with submit_compute_sea_network
#' Add Q-values using the fdrtool package
#'
#' @param run_base returned by submit_compute_sea_network(...)
#' @export
collect_compute_sea_network <- function(
	run_base,
	verbose
){
	if(verbose){
		cat("Collecting scores from '", paste0(run_base, "/outputs' ...\n"), sep="")
	}
	scores <- list.files(paste0(run_base, "/outputs")) %>%
		plyr::ldply(function(target){
			df <- readr::read_csv(
				paste0(run_base, "/outputs/", target, "/", target, ".csv.out.csv"),
				col_types=cols(
					`Query ID` = col_character(),
					`Target ID` = col_character(),
					`Affinity Threshold (nM)` = col_integer(),
					`P-Value` = col_double(),
					`Cut Sum` = col_double(),
					`Max Tc` = col_double(),
					`Z-Score` = col_double(),
					Name = col_character(),
					Description = col_character())) %>%
				dplyr::transmute(
					target1 = `Target ID`,
					target2 = `Query ID`,
					affinity = `Affinity Threshold (nM)`,
					Pvalue = `P-Value`,
					CutSum = `Cut Sum`,
					MaxTC = `Max Tc`,
					Zscore = `Z-Score`)
		})

	if(verbose){
		cat("Computing Q-values ...\n")
	}
	ref_targets <- scores %>% dplyr::distinct(target1) %>% magrittr::extract2("target1")
	query_targets <- scores %>% dplyr::distinct(target2) %>% magrittr::extract2("target2")
	fdr <- expand.grid(
		target1 = ref_targets,
		target2 = query_targets) %>%
		dplyr::left_join(
			scores %>% dplyr::select(target1, target2, Pvalue),
			by=c("target1", "target2")) %>%
		dplyr::mutate(
			Pvalue = ifelse(is.na(Pvalue), runif(n()), Pvalue))

	n_ref_target <- length(ref_targets)
	n_query_target <- length(query_targets)
	n_shared_targets <- length(ref_targets[ref_targets %in% query_targets])
	
	a <- fdrtool::fdrtool(fdr$Pvalue, statistic="pvalue")
	fdr <- fdr %>%
		dplyr::mutate(
			Qvalue=a$qval,
			Evalue=Pvalue * ((n_ref_target * n_query_target)/2 - n_shared_targets)) %>%
		dplyr::select(-Pvalue)

	# this is slow because there are lot of fdr values
	scores <- scores %>%
		left_join(fdr, by=c("target1", "target2"))

	s <- scores
}

