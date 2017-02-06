# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#' Config files should be a list <config_fname>=<config_file_contents>
#' for one of the sea config files.
#' these are written to the current working directory so they can be parsed
#' so as to take presidence over other specified configuration options
#' @export
write.config_files <- function(fitcore=NULL, fpcore=NULL, seacore=NULL, sea_interface=NULL){
  library(readr)
  config_files = as.list(match.call())
  for( config_name in names(config_files)[-1] ){
      config_fname <-paste( getwd(), "/", config_name, ".cfg", sep="")
      if(file.exists(config_fname)){
          existing_config_file <- readr::read_file(config_fname)
          if(existing_config_file != config_files[[config_name]]){
              stop(paste("Failed to write '", config_fname, "' because the file already exists and is different from the provided one, please resolve this discrepancy before trying again.", sep=""))
          }
      }
      cat(config_files[[config_name]], file=config_fname, sep="", append=FALSE)
  }
}

#' @export
cleanup.config_files <- function(
    fitcore=TRUE, fpcore=TRUE, seacore=TRUE, sea_interface=TRUE){

    config_files = as.list(match.call())
    for( config_name in names(config_files) ){
        config_fname <-paste( getwd(), "/", config_name, ".cfg", sep="")
        if(!file.exists(config_fname)){
            continue
        }
        tryCatch({
            file.remove(config_fname)
        }, error=function(e){
            print(paste("Unable to clean up config '", config_fname, "'.", sep=""))
        })
    }
}

#' check that the smiles and compound columns are unique
#' @export
check_duplicates.smi <-function(smi, smiles="smiles", compound="compound") {
	duplicates <- smi %>%
		dplyr::group_by_(compound) %>%
			dplyr::tally() %>%
			dplyr::ungroup() %>%
		dplyr::filter(n > 1) %>%
		dplyr:::select(-n) %>%
		dplyr::left_join(smi, by=compound)
	if(nrow(duplicates) > 1) {
		cat("WARNING: Multiple smiles strings for the same compound were found:\n")
		print(duplicates %>% head(20) %>% as.data.frame )
		cat("WARNING: You can removed duplicates by running clean_duplicates.smi\n")
	}

	duplicates <- smi %>%
		dplyr::group_by_(smiles) %>%
			dplyr::tally %>%
			dplyr::ungroup %>%
		dplyr::filter(n > 1) %>%
		dplyr:::select(-n) %>%
		dplyr::left_join(smi, by=smiles)
	if(nrow(duplicates) > 1) {
		cat("WARNING: Multiple compounds having the same smiles string were found:\n")
		print(duplicates %>% head(20) %>% as.data.frame)
		cat("WARNING: You can removed duplicates by running clean_duplicates.smi\n")
	}
}

#' @export
clean_duplicates.smi <- function(smi, smiles="smiles", compound="compound"){
	smi %>%
		dplyr::group_by_(compound) %>%
			dplyr::slice(1) %>%
		dplyr::ungroup %>%
		dplyr::group_by_(smiles) %>%
			dplyr::slice(1) %>%
		dplyr::ungroup
}


#' returns data.frame with columns [smiles, compound]
#' @export
read.smi <- function(fname, smiles="smiles", compound="compound"){
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .smi file, '", fname, "', does not exist.", sep=""))
	}
	z <- data.table::fread(fname, sep=";", header=F, data.table=F)
	data.table::setnames(z, c(smiles, compound))
	check_duplicates.smi(z)
	z
}


#' combine compounds together, removing duplicates -- warn if there are
#' different smiles strings for the same compound
#'
#' INPUT:
#'   smis: named list of smi data.frames with at least columns [smiles, compound]
#'   smiles: name of smiles column in smi data.frames
#'   compound: name of compound column in smi data.frame
#'
#' OUTPUT:
#'   data.frame: columns of smis and column <set>
#' @export
combine.smi <- function(smis, smiles="smiles", compound="compound") {
	z <- plyr::ldply(smis, function(df) df, .id="set")
	check_duplicates.smi(z)
	z
}



#' input data.frame with columns [smiles, compound]
#' write .smi file to given fname
#'
#' filter unique compound identifiers
#' @export
write.smi <- function(
	smi,
	fname,
	smiles="smiles",
	compound="compound",
	seaware_format=T,
	verbose=F){

	if(verbose){
		cat("Writing smiles data frame to '", fname, "':\n", sep="")
		cat("\tsmiles column: ", smiles, "\n", sep="")
		cat("\tcompound column: ", compound, "\n", sep="")
		if(seaware_format){
			cat("\tusing the seaware format <cid>,<smiles>\n")
		} else{
			cat("\tusing the legacy format <smiles>;<cid>\n")
		}
	}

	smi <- as.data.frame(smi)

	if(seaware_format){
		if(!stringr::str_detect(fname, ".csv$")){
			stop(paste0("WARNING: requested to write to output smiles file '", fname, "'. For SEA calculations it should have a '.csv' extension."))
		}
	}else{
		if(!stringr::str_detect(fname, ".smi$")){
			stop(paste0("WARNING: requested to write to output .smi file '", fname, "'. For SEA calculations it should have a '.smi' extension."))
		}
	}

	if(!file.exists(dirname(fname))){
		stop(paste0("ERROR: Cannot write '", fname, "', because the directory '", dirname(fname), "' does not exist."))
	}

	check_column <- function(col_type, col_name){
		if( !(col_name %in% names(smi))){
			if(seware_format){
				stop(paste0(
					"ERROR: cannot write smiles file '", fname, "' because the input data.frame does not have column ", col_type, " having name '", col_name, "'\n",
					"ERROR:    input data.frame has columns: [", paste(names(smi), collapse=", "), "]\n",
					"ERROR:    required columns for seaware smiles format: [", compound, ", ", smiles, "]"))
			} else {
				stop(paste0(
					"ERROR: cannot write .smi file '", fname, "' because the input data.frame does not have column ", col_type, " having name '", col_name, "'\n",
					"ERROR:    input data.frame has columns: [", paste(names(smi), collapse=", "), "]\n",
					"ERROR:    required columns for old-sea .smi format: [", smiles, ", ", compound, "]"))
			}
		}
		if(any(is.na(smi[[col_name]]))){
			stop(paste0("ERROR: ", sum(is.na(smi[[col_name]])), " ", col_type, " entries are NA."))
		}

		if(!seaware_format){
			m <- stringr::str_detect(smi[[col_name]], fixed(";"))
			if(any(m, na.rm=T)){
				stop(paste0("ERROR: ", sum(m), " ", col_type, " values in column '", target, "' contain ';', which is used as the separator for SEA input files.\n"))
			}
		}
	}

	check_column("compound", compound)
	check_column("smiles", smiles)
	check_duplicates.smi(smi=smi, smiles=smiles, compound=compound)

	if(seaware_format){
		write.table(
			as.data.frame(smi)[,c(compound, smiles)],
			fname,
			quote=T,
			sep=",",
			row.names=F,
			col.names=c("molecule id", "smiles"))
	} else {
		write.table(
			as.data.frame(smi)[, c(smiles,compound)],
			fname,
			quote=F,
			sep=";",
			row.names=F,
			col.names=F)
	}
	fname
}



#' @export
write.singleton.set <- function(smi, fname, target="compound", name="compound", compound="compound"){
	if( !(target %in% names(smi))){
		stop(paste0("ERROR: cannot write .set file '", fname, "' because the input data.frame does not have column '", target, "'\n\t data.frame columns: [", paste(names(smi), collapse=", "), "\n\trequired columns: [", paste(unique(target, name, compound), collapse=", "), "]"))
	}
	if( !(name %in% names(smi))){
		stop(paste0("ERROR: cannot write .set file '", fname, "' because the input data.frame does not have column '", name, "'\n\t data.frame columns: [", paste(names(smi), collapse=", "), "]\n\trequired columns: [", paste(unique(target, name, compound), collapse=", "), "]"))
	}
	if( !(compound %in% names(smi))){
		stop(paste0("ERROR: cannot write .set file '", fname, "' because the input data.frame does not have column '", compound, "'\n\t data.frame columns: [", paste(names(smi), collapse=", "), "]\n\trequired columns: [", paste(unique(target, name, compound), collapse=", "), "]"))
	}

	if(!stringr::str_detect(fname, ".set$")){
		stop(paste0("WARNING: requested to write to output .set file '", fname, "'. For SEA calculations it should have a '.set' extension."))
	}

	if(!file.exists(dirname(fname))){
		stop(paste0("ERROR: cannot write '", fname, "', because the directory '", dirname(fname), "' does not exist."))
	}

	write.table(
		as.data.frame(smi)[,c(target, name, compound)],
		fname,
		quote=F,
		sep=";",
		row.names=F,
		col.names=F)
}

pack.library <- function(
	molecules,
	targets,
	library_fname = NULL,
	fingerprint_type="rdkit_ecfp4",
	name=NULL,
	config_files = NULL,
	verbose=F){

	tmp_base <- tempfile()
	molecules_fname <- paste0(tmp_base, ".csv")
	targets_fname <- paste0(tmp_base, ".set")
	if(is.null(library_fname)){
		library_fname <- paste0(tmp_base, ".sea")
	}

	write.smi(molecules, molecules_fname)
	write.set(targets, targets_fname)

	do.call(write.config_files, config_files)

	cmd <- paste(
		"SeaShell library pack",
		"--generate-fingerprint", fingerprint_type,
		ifelse(!is.null(name), paste("--name", shQuote(name)), ""),
		library_fname, molecules_fname, targets_fname)
	if(verbose){
		cat(cmd, "\n", sep="")
	}

	system(cmd)
#	file.remove(molecules_fname)
#	file.remove(targets_fname)
#	cleanup.config_files()
	library_fname
}


#' unpack a library
#' returns: list(molecules, targets) where
#'   molecules is a data.frame with columns: compound, smiles, fingerprint
#'   targets is a data.frame with columns: target, name, affinity, compound, description
#'
#' Note: the molecules and targets tables are joined on the compound in the molecules and targets is the join column
#' @export
unpack.library <- function(
	fname,
	fingerprint_format="sea_native",
	molecules_fname = NULL,
	targets_fname = NULL,
	unpack_compounds=TRUE,
	verbose=FALSE){

	if(!is.character(fname) | !file.exists(fname)){
		stop(paste0("ERROR: The SEA .smi file, '", fname, "', does not exist."))
	}

	if(is.null(molecules_fname)){
		molecules_fname <- paste0(tempfile(), ".csv")
	}

	if(is.null(targets_fname)){
		targets_fname <- paste0(tempfile(), ".set")
	}

	cmd <- paste("SeaShell library unpack --yes --fingerprint-format", fingerprint_format, fname, molecules_fname, targets_fname)
	if(verbose){
			cat("Unpacking targets and molecules form '", fname, "'\n", sep="")
			cat(cmd, "\n", sep="")
	}
	system(cmd)

	molecules <- data.table::fread(molecules_fname, sep=",", header=T, data.table=F, skip=1) %>%
		dplyr::select(
			compound=`molecule id`,
			smiles,
			fingerprint)

	targets <- data.table::fread(targets_fname, sep=",", header=T, data.table=F) %>%
		dplyr::select(
			target=`target id`,
			name,
			affinity,
			compound=molecules,
			description)

	if(unpack_compounds){
		targets <- targets %>%
			transform(compound = stingr::str_split(compound, ":")) %>%
			tidyr::unnest(compound)
	}
	list(molecules=molecules, targets=targets)
}


plot_fits.library <- function(
	library_fname,
	diagnostics_fname=library_fname,
	verbose=F){

	if(verbose){
		cat("Making background ...\n")
	}
	cmd <- paste(
		"SeaShell library fit --yes",
		"--plot", diagnostics_fname,
		library_fname)

	if(verbose){
		cat(cmd, "\n", sep="")
	}

	system(cmd)

	if(verbose){
		cat("check '", diagnostics_fname, ".dist' and '", diagnostics_fname, ".png'\n", sep="")
		cat("looking up best fit ... \n")
	}
}

set_fit.library <- function(
	library_fname,
	threshold=.28,
	diagnostics_fname=library_fname,
	verbose=F){

	cmd <- paste(
		"SeaShell library fit --yes",
		"--select", threshold, diagnostics_fname,
		library_fname)
	if(verbose){
		cat(cmd, "\n", sep="")
	}
	system(cmd)
}


#' returns data.frame with columns [target, name, compound]
#' if unpack_compounds, put each compound on separate row
#' @export
read.set <- function(
	fname,
	target="target",
	name="name",
	compound="compound",
	unpack_compounds=T){
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .set file, '", fname, "', does not exist.", sep=""))
	}
	sets <- data.table::fread(fname, sep=";", header=F, data.table=F)
	sets %>% data.table::setnames(c(target, name, compound))
	if(unpack_compounds) {
		sets <- sets %>%
			transform(compound = stringr::str_split(compound, ":")) %>%
			tidyr::unnest(compound)
	}
	sets
}

#' INPUT:
#'   set: data.frame with columns [target, name, affinity, description, compound]
#'        if affinity and description are not present, assume it is the old style
#'        format and don't write them out.
#'
#'        warning note the order of the columns in the input idfferes in the order of the output
#'   set_fname: .set file name
#' OUTPUT:
#'   creates <set_fname>
#' @export
write.set <- function(
	set,
	fname,
	target="target",
	name="name",
	compound="compound",
	affinity="affinity",
	description="description",
	seaware_format=T){

	if(!stringr::str_detect(fname, ".set$")){
		stop(paste0("WARNING: requested to write to output .set file '", fname, "'. For SEA calculations it should have a '.set' extension."))
	}

	if(!file.exists(dirname(fname))){
		stop(paste0("ERROR: cannot write '", fname, "', because the directory '", dirname(fname), "' does not exist."))
	}

	if( "affinity" %in% names(set) || "description" %in% names(set) ){
		cat("INFO: Using seaware .set format: <target>;<name>;<affinity>;<description>;<compounds>\n")
		if(!seaware_format){
			stop("Please specify seaware_format=True to use seaware .set format")
		}
		seaware_format <- TRUE
	} else {
		cat("INFO: Using old-sea .set format: <starget>;<name>;<compounds>\n")
		if(seaware_format){
			stop("please specify seaware_format=False to use old-sea .set format")
		}
		seaware_format <- FALSE
	}

	check_column <- function(col_type, col_name){
		if( !(col_name %in% names(set))){
			if(seware_format){
				stop(paste0(
					"ERROR: cannot write .set file '", fname, "' because the input data.frame does not have column ", col_type, " having name '", col_name, "'\n",
					"ERROR:    input data.frame has columns: [", paste(names(set), collapse=", "), "]\n",
					"ERROR:    required columns for seaware .set format: [", target, ", ", name, ", ", affinity, ",", description, ",", compound, "]"))
			} else {
				stop(paste0(
					"ERROR: cannot write .set file '", fname, "' because the input data.frame does not have column ", col_type, " having name '", col_name, "'\n",
					"ERROR:    input data.frame has columns: [", paste(names(set), collapse=", "), "]\n",
					"ERROR:    required columns for old-sea .set format: [", target, ", ", name, ", ", compound, "]"))
			}
		}
		if(any(is.na(set[[col_name]]))){
			stop(paste0("ERROR: ", sum(is.na(set[[col_name]])), " ", col_type, " entries are NA."))
		}

		if(!seaware_format){
			m <- stringr::str_detect(set[[col_name]], fixed(";"))
			if(any(m, na.rm=T)){
				stop(paste0("ERROR: ", sum(m), " ", col_type, " values in column '", target, "' contain ';', which is used as the separator for SEA input files.\n"))
			}
		}
	}

	check_column("target", target)
	check_column("name", name)
	check_column("compound", compound)
	check_column("affinity", affinity)
	check_column("description", description)

	if(seaware_format){
		set %>%
			group_by_("target") %>%
				dplyr:::summarize(
					name=name[1],
					affinity=affinity[1],
					compounds=paste(compound, collapse=":"),
					description=description[1]) %>%
			write.table(fname, quote=T, sep=",", row.names=F, col.names=T)
	} else {
		set %>%
			dplyr::group_by_("target") %>%
				dplyr:::summarize(
					name[1],
					paste(compound, collapse=":")) %>%
			write.table(fname, quote=F, sep=";", row.names=F, col.names=F)
	}
}


#' create a fingerprint file for the given file or data frame with columns [compound, smiles]
#' return based on the read parameter
#' @export
create.fp <- function(
	smi,
	compound="compound",
	smiles="smiles",
	output_fname=NULL,
	fp_format=c('sea_native', 'bitstring'),
	fp_type=c('rdkit_ecfp', 'rdkit_path', 'rdkit_ecfp_mask'),
	standardize=T,
	verbose=F){

	if(is.character(smi)){
		smi_fname <- smi
	} else if(is.data.frame(smi)){
		smi_fname <- paste0(tempfile(), ".smi")
		write.smi(smi, smi_fname, compound=compound, smiles=smiles, seaware_format=T, verbose)
		cat("DONE\n")
	} else{
		cat("ERROR: unrecognized smi class: ", paste0(class, collapse=", "), "\n", sep="")
	}

	if(is.null(output_fname)){
		output_fname <- tempfile()
	}
	fp_format <- match.args(fp_format)
	fp_type <- match.args(fp_type)

	if(verbose){
		cat("Writing a ", fp_format, " ", fp_type, " fingerprint file to '", output_fname, "'...\n", sep="")
	}

	cmd <- paste0(
		"SeaShell fingerprint ",
		"--to-format ", fp_format, " ",
		ifelse(standardize, '', '--skip-standardization '),
		"--generate-fingerprint ", fp_type, " ",
		"'", smi_fname, "' ",
		"'", output_fname, "'")

	cat(cmd, "\n")
	system(cmd)
}

#' @export
read.fp <- function(fname, fingerprint="fingerprint", compound="compound"){
	require(sqldf)
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .fp file, '", fname, "', does not exist.", sep=""))
	}
	fp <- readr::read_delim(fname, colnames=, delim=";")
	names(fp) <- c(fingerprint, compound)
	fp
}

#' @export
read.fp.compound <- function(fname, compound="compound"){
	require(sqldf)
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .fp file, '", fname, "', does not exist.", sep=""))
	}
	fp <- read.csv.sql(fname, "SELECT V2 FROM file", header=F, sep=";")
	names(fp) <- c(compound)
	fp
}

#' @export
write.fp <- function(fp, fname, fingerprint="fingerprint", compound="compound"){
	if( !(fingerprint %in% names(fp))){
		stop(paste0("ERROR: cannot write .fp file '", fname, "' because the input data.frame does not have column '", fingerprint, "'\n\t data.frame columns: [", paste(names(fp), collapse=", "), "]\n\trequired columns: [", fingerprint, ", ", compound, "]"))
	}
	if( !(compound %in% names(fp))){
		stop(paste0("ERROR: cannot write .fp file '", fname, "' because the input data.frame does not have column '", compound, "'\n\t data.frame columns: [", paste(names(fp), collapse=", "), "]\n\trequired columns: [", fingerprint, ", ", compound, "]"))
	}

	if(!stringr::str_detect(fname, ".fp$")){
		stop(paste0("WARNING: requested to write to output fingerprint file '", fname, "'. For SEA calculations it should have a '.fp' extension."))
	}

	if(!file.exists(dirname(fname))){
		stop(paste0("ERROR: cannot write '", fname, "', because the directory '", dirname(fname), "' does not exist."))
	}

	write.table(
		fp[,c(fingerprint, compound)],
		fname,
		quote=F,
		sep=";",
		row.names=F,
		col.names=F)
}

#' @export
read.scores <- function(fname, use_cache=T){
	if(!is.character(fname) || !file.exists(fname)){
		stop(paste("ERROR: The SEA .scores.csv file, '", fname, "', does not exist.", sep=""))
	}

	cache_fname <- str_replace(fname, ".csv", ".Rdata")
	if(use_cache && file.exists(cache_fname)){
		envir = environment()
		load(cache_fname, envir=envir)
		return(envir$scores)
	}

	scores <- data.table::fread(fname, sep=",", header=T, data.table=F)

	names_map <- c(
			"EValue"="EValue",
			"EValue"="e-value",
			"MaxTC"="MaxTC",
			"MaxTC"="max tc",
			"ZScore"="ZScore",
			"ZScore"="zscore",
			"target1"="target1",
			"target2"="target2",
			"target1_Ki"="target1_Ki",
			"target2_Ki"="target2_ki")

	data.table::setnames(
		scores,
		old=names_map[names_map %in% names(scores)],
		new=names(names_map)[names_map %in% names(scores)])

	if(use_cache) {
		save(scores, file=cache_fname)
	}

	scores
}

#' @export
write.scores <- function(scores, fname){
	write.table(
		scores,
		fname,
		quote=F,
		sep=",",
		row.names=F,
		col.names=F)
}

# run SEA
#
# inputs:
#   background_fname: filename of .fit file
#   ref_set: compound sets for reference targets
#     path to .set file or
#     data.frame with columns [target, name, compound]
#   ref_smi: compounds for reference targets
#     path to .smi file or
#     data.frame with columns [compound, smiles]
#   query_set: compound sets for query targets
#     path to .set file or
#     data.frame with columns [target, name, compound]
#   query_smi: compounds for query targets
#     path to .smi file or
#     data.frame with columns [compound, smiles]
#   run_tag: a tag given help identify the sea analysis in the SEA viewer
#     default NULL
#   output_base: path and file prefix where the files should be written
#     default to a temporary file

# Output:
#   data.frame with columns [target1, target2, EValue, MaxTC] for SEA associations
run.sea <- function(
	background_fname,
	ref_set,
	ref_smi,
	query_set,
	query_smi,
	run_tag=NULL,
	output_base=NULL,
	similarity_measure="tanimoto",
	tversky_alpha=0.5,
	tversky_beta=0.5
){
	if(!stringr::str_detect(background_fname, ".fit$")){
		cat("WARNING: background .fit file '", background_fname, "', does not end in '.fit'\n", sep="")
	}

	if(is.null(output_base)){
		output_base <- tempfile()
	}
	input_base <- tempfile()

	# the run tag is helpful for looking at the results in the sea viewer
	if(!is.null(run_tag)){
		input_base <- paste(input_base, run_tag, sep=".")
	}

	if(similarity_measure=="tversky"){
		similarity_measure <- paste0("\"tversky -a ", tversky_alpha, " -b ", tversky_beta, "\"")
	}

	if(is.character(ref_set)){
		if(!stringr::str_detect(ref_set, ".set$")){
			cat("WARNING: reference .set file '", ref_set, "', does not end in '.set'\n", sep="")
		}
		ref_set_fname <- ref_set
	} else {
		ref_set_fname <- paste0(input_base, ".ref.set")
		cat("writing ref_set to -> '", ref_set_fname, "' ... ", sep="")
		write.set(ref_set, ref_set_fname)
		cat("DONE\n")
	}

	if(is.character(query_set)){
		if(!stingr::str_detect(query_set, ".set$")){
			cat("WARNING: query .set file '", query_set, "', does not end in '.set'\n", sep="")
		}
		query_set_fname <- query_set
	} else {
		query_set_fname <- paste0(input_base, ".query.set")
		cat("writing query_set to -> '", query_set_fname, "' ... ", sep="")
		write.set(query_set, query_set_fname)
		cat("DONE\n")
	}

	if(is.character(ref_smi)){
		if(!stringr::str_detect(ref_smi, ".smi$")){
			cat("WARNING: referene .smi file '", ref_smi, "', does not end in '.smi'\n", sep="")
		}
		ref_smi_fname <- ref_smi
	} else {
		ref_smi_fname <- paste0(input_base, ".ref.smi")
		cat("writing ref_smi to -> '", ref_smi_fname, "' ... ", sep="")
		write.smi(ref_smi, ref_smi_fname)
		cat("DONE\n")
	}

	if(is.character(query_smi)){
		if(!stringr::str_detect(query_smi, ".smi$")){
			cat("WARNING: query .smi file '", query_smi, "', does not end in '.smi'\n", sep="")
		}
		query_smi_fname <- query_smi
	} else {
		query_smi_fname <- paste0(input_base, ".query.smi")
		cat("writing query_smi to -> '", query_smi_fname, "' ... ", sep="")
		write.smi(query_smi, query_smi_fname)
		cat("DONE\n")
	}


	cmd <- paste(
			paste0(sea_scripts_dir, "/run_sea.sh"),
			background_fname,
			similarity_measure,
			ref_set_fname,
			query_set_fname,
			output_base)
	cat(cmd, "\n", sep="")
	system(cmd)
	scores <- read.scores(paste0(output_base, "/sea.scores.csv"))
}




# build redirect urls into the sea viewer a list of target pairs
#
#  sea_locator_url: url to the locate_association.php script
#    you may need to copy this file into the sea viewer base directory.
#  sea_analsis_id: each sea run is assigned an analysis id
#  query_targets, reference_targets: vectors of target codes
#  flip: if true put the query target downstairs
#  return:
#     a character vector of urls into the SEA viewer
sea_viewer_locator <- function(
	query_targets,
	query_activity_thresholds,
	reference_targets,
	reference_activity_thresholds,
	sea_locator_url="http://sea16.ucsf.bkslab.org/custom/search/target_vs_target"){
	if(length(query_targets) != length(reference_targets)){
		stop(paste0(
			"target1 and target2 should be two vectors of target codes the same length, ",
			"but the length of query targets is ", length(query_targets), ", ",
			"and the length of the reference target is '", length(reference_targets), "'"))
	}

	paste0(
		sea_locator_url,
		"?library=chembl_10uM",
		"&target_id=", paste(reference_targets, reference_activity_thresholds, sep="+"),
		"&query_id=", paste(query_targets, query_activity_thresholds, sep="+"),
		"&rank_cutoff=75")
}

#' scores is a data.frame with columns target1, target2 as uniprot_entries
#' @export
hit.report <- function(
	hits
){
	library(BioChemPantry)
	pantry <- get_pantry()

	targets <- pantry %>% schema_tbl("hgnc_151120.genes")
	sea_scores <- pantry %>% schema_tbl("sea_chembl20.scores")
	sequence_similarity <- pantry %>% schema_tbl("sea_chembl20.blastp_target_vs_target")
	co_expression <- pantry %>% schema_tbl("coexpnet.human")
	phenotype_associations <- pantry %>%	schema_tbl("phenocarta_151014.associations")
	ppi <- pantry %>% schema_tbl("biogrid_3_2_121.human_silver")

	uniprot_chordate_taxa <- Zr::uniprot_taxa(
		user_agent=httr::user_agent("httr mattjomeara@gmail.com"),
		ancestor=7711)

	check_column <- function(col_name){
		if( !(col_name %in% names(hits))){
			stop(paste0(
				"ERROR: cannot make a hit report file '", fname, "' because the input data.frame does not have column '", col_name, ",\n",
				"ERROR:    input data.frame has columns: [", paste(names(hits), collapse=", "), "]\n",
				"ERROR:    required columns: [ target1, target2, ...]"))
		}
	}
	check_column("target1")
	check_column("target2")
	additional_cols <- names(hits)
	additional_cols <- additional_cols[which(!(additional_cols %in% c("target1", "target2")))]

	hits <- hits %>%
		dplyr::mutate(
			mnemonic1 = stringr::str_extract(target1, "[A-Z0-9]+$"),
			mnemonic2 = stringr::str_extract(target2, "[A-Z0-9]+$"),
			is_human1 = mnemonic1 == "HUMAN",
			is_human2 = mnemonic2 == "HUMAN",
			viewer_link = sea_viewer_locator(target1, 5, target2, 5)) %>%
		dplyr::left_join(
			uniprot_chordate_taxa %>% dplyr::select(mnemonic1=Mnemonic, taxon1 = Taxon, common_name1 = `Common name`),
			by="mnemonic1") %>%
		dplyr::left_join(
			uniprot_chordate_taxa %>% dplyr::select(mnemonic2=Mnemonic, taxon2 = Taxon, common_name2 = `Common name`),
			by="mnemonic2") %>%
		dplyr::mutate(
			is_chordate1 = !is.na(taxon1),
			is_chordate2 = !is.na(taxon2)) %>%
		dplyr::select(-taxon1, -taxon2, -mnemonic1, -mnemonic2)

	hits_db <- dplyr::copy_to(pantry, df=hits, name="hits", temporary=T)

	hits_p1 <- hits_db %>%
		dplyr::inner_join(
			sequence_similarity %>%
				dplyr::filter(EValue < sql("1e-50")) %>%
				dplyr::select(target1, target1p = target2, p1_blast_Evalue=EValue),
			by=c("target1")) %>%
		dplyr::inner_join(
			sea_scores %>%
				dplyr::filter(MaxTC == 1) %>%
				dplyr::select(
					target1p = target1, target2, p1_sea_Qvalue = Qvalue),
			by=c("target1p", "target2")) %>%
		dplyr::collapse() %>%
		dplyr::group_by(target1, target2) %>%
			dplyr::arrange(p1_sea_Qvalue, p1_blast_Evalue) %>%
			dplyr::filter(row_number() <= 5) %>%
			dplyr::collapse() %>% # https://github.com/hadley/dplyr/issues/1680
			dplyr::summarize(
				p1_blast_Evalue = paste(sql("to_char(\"p1_blast_Evalue\", '9D99EEEE')"), collapse=";"),
				p1_sea_Qvalue = paste(sql("to_char(\"p1_sea_Qvalue\", '9D99EEEE')"), collapse=";"),
				target1p = paste(target1p, collapse = ";")) %>%
		dplyr::ungroup()

	hits_p2 <- hits_db %>%
		dplyr::inner_join(
			sequence_similarity %>%
				dplyr::filter(EValue < sql("1e-50")) %>%
				dplyr::select(target2p = target1, target2, p2_blast_Evalue=EValue),
			by=c("target2")) %>%
		dplyr::inner_join(
			sea_scores %>%
				dplyr::filter(MaxTC == 1) %>%
				dplyr::select(target1, target2p = target2, p2_sea_Qvalue = Qvalue),
			by=c("target1", "target2p")) %>%
		dplyr::collapse() %>%
		group_by(target1, target2) %>%
			dplyr::arrange(p2_sea_Qvalue, p2_blast_Evalue) %>%
			dplyr::filter(row_number() <= 5) %>%
			dplyr::collapse() %>% # https://github.com/hadley/dplyr/issues/1680
			dplyr::summarize(
				p2_blast_Evalue = paste(sql("to_char(\"p2_blast_Evalue\", '9D99EEEE')"), collapse=";"),
				p2_sea_Qvalue = paste(sql("to_char(\"p2_sea_Qvalue\", '9D99EEEE')"), collapse=";"),
				target2p = paste(target2p, collapse = ";")) %>%
		dplyr::ungroup()

	hits <- hits_db %>%
		dplyr::left_join(
			targets %>%
				dplyr::select(target1 = uniprot_entry, gene1=symbol),
			by="target1") %>%
		dplyr::left_join(
			targets %>%
				dplyr::select(target2 = uniprot_entry, gene2=symbol),
			by="target2") %>%
		dplyr::left_join(
			sea_scores %>%
				dplyr::mutate(
					sea_Qvalue = Qvalue,
					sea_MaxTC = MaxTC) %>%
				dplyr::select(
					target1, target2,
					description1, description2,
					target1_class_1, target1_class_2,
					target2_class_1, target2_class_2,
					sea_Qvalue, sea_MaxTC),
			by=c("target1", "target2")) %>%
		dplyr::left_join(
			sequence_similarity %>%
				dplyr::mutate(
					blast_Evalue = EValue) %>%
				dplyr::select(target1, target2, blast_Evalue),
				by=c("target1", "target2")) %>%
		dplyr::left_join(
			phenotype_associations %>%
				dplyr::select(target1, target2, phenotype_name),
			by=c("target1", "target2")) %>%
		dplyr::left_join(
			co_expression %>%
				dplyr::mutate(
					coexp_rank_fraction = rank_fraction) %>%
				dplyr::select(target1, target2, coexp_rank_fraction),
			by=c("target1", "target2")) %>%
		dplyr::left_join(
			ppi %>%
				dplyr::mutate(BioGRID_silver=1) %>%
				dplyr::select(target1, target2, BioGRID_silver),
			by=c("target1", "target2")) %>%
		dplyr::left_join(hits_p1, by=c("target1", "target2")) %>%
		dplyr::left_join(hits_p2, by=c("target1", "target2")) %>%
		dplyr::select_(.dots=c(
			"target1", "target2", "description1", "description2",
			"is_human1", "is_human2", "is_chordate1", "is_chordate2",
			"target1_class_1", "target1_class_2", "target2_class_1", "target2_class_2",
			additional_cols,
			"sea_MaxTC", "sea_Qvalue", "viewer_link",
			"blast_Evalue", "phenotype_name", "coexp_rank_fraction", "BioGRID_silver",
			"target1p", "p1_blast_Evalue", "p1_sea_Qvalue",
			"target2p", "p2_blast_Evalue", "p2_sea_Qvalue") %>%
			lapply(as.name)) %>% # https://github.com/hadley/dplyr/issues/1392
		dplyr::collect()
}
