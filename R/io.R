# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#' returns data.frame with columns [smiles, compound]
#' @export
read.smi <- function(fname, smiles="smiles", compound="compound"){
	require(data.table)
	require(plyr)
	require(dplyr)
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .smi file, '", fname, "', does not exist.", sep=""))
	}
	z <- fread(fname, sep=";", header=F, data.table=F)
	setnames(z, c(smiles, compound))
	check_duplicates.smi(z)
	z
}


#' input data.frame with columns [smiles, compound]
#' write .smi file to given fname
#'
#' filter unique compound identifiers
#' @export
write.smi <- function(smi, fname, smiles="smiles", compound="compound"){
	require(plyr)
	require(dplyr)
	require(data.table)
	smi <- as.data.frame(smi)
	if( !(smiles %in% names(smi))){
		stop(paste0("ERROR: cannot write .smi file '", fname, "' because the input data.frame does not have column '", smiles, "'\n\t data.frame columns: [", paste(names(smi), collapse=", "), "]\n\trequired columns: [", smiles, ", ", compound, "]"))
	}
	if( !(compound %in% names(smi))){
		stop(paste0("ERROR: cannot write .smi file '", fname, "' because the input data.frame does not have column '", compound, "'\n\t data.frame columns: [", paste(names(smi), collapse=", "), "]\n\trequired columns: [", smiles, ", ", compound, "]"))
	}

	if(!str_detect(fname, ".smi$")){
		stop(paste0("WARNING: requested to write to output .smi file '", fname, "'. For SEA calculations it should have a '.smi' extension."))
	}

	if(!file.exists(dirname(fname))){
		stop(paste0("ERROR: Cannot write '", fname, "', because the directory '", dirname(fname), "' does not exist."))
	}

	if(sum(is.na(smi[[compound]]) > 0)){
		stop(paste0("ERROR: Cannot write '", fname, "', because the compound column contains '", sum(is.na(smi[[compound]])), "' NA values."))
	}

	if(sum(is.na(smi[[smiles]]) > 0)){
		stop(paste0("ERROR: Cannot write '", fname, "', because the smiles column contains '", sum(is.na(smi[[smiles]])), "' NA values."))
	}

	check_duplicates.smi(smi)

	write.table(
		as.data.frame(smi)[, c(smiles,compound)],
		fname,
		quote=F,
		sep=";",
		row.names=F,
		col.names=F)
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

	if(!str_detect(fname, ".set$")){
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

#' returns data.frame with columns [target, name, compound]
#' if unpack_compounds, put each compound on separate row
#' @export
read.set <- function(
	fname,
	target="target",
	name="name",
	compound="compound",
	unpack_compounds=T){
	require(data.table)
	require(dplyr)
	require(stringr)
	require(tidyr)
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .set file, '", fname, "', does not exist.", sep=""))
	}
	sets <- fread(fname, sep=";", header=F, data.table=F)
	sets %>% setnames(c(target, name, compound))
	if(unpack_compounds) {
		sets <- sets %>%
			transform(compound = str_split(compound, ":")) %>%
			unnest(compound)
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
	description="description"){
	require(plyr)
	require(dplyr)
	require(stringr)

	if(!str_detect(fname, ".set$")){
		stop(paste0("WARNING: requested to write to output .set file '", fname, "'. For SEA calculations it should have a '.set' extension."))
	}

	if(!file.exists(dirname(fname))){
		stop(paste0("ERROR: cannot write '", fname, "', because the directory '", dirname(fname), "' does not exist."))
	}

	if( "affinity" %in% names(set) || "description" %in% names(set) ){
		cat("INFO: Using seaware .set format: <target>;<name>;<affinity>;<description>;<compounds>\n")
		seware_format <- TRUE
	} else {
		cat("INFO: Using old-sea .set format: <starget>;<name>;<compounds>\n")
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

		m <- str_detect(set[[col_name]], fixed(";"))
		if(any(m, na.rm=T)){
			stop(paste0("ERROR: ", sum(m), " ", col_type, " values in column '", target, "' contain ';', which is used as the separator for SEA input files.\n"))
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
					name[1],
					activity[1],
					description[1],
					paste(compound, collapse=":")) %>%
			write.table(fname, quote=F, sep=";", row.names=F, col.names=F)
	} else {
		set %>%
			group_by_("target") %>%
				dplyr:::summarize(
					name[1],
					paste(compound, collapse=":")) %>%
			write.table(fname, quote=F, sep=";", row.names=F, col.names=F)
	}
}


#' create a fingerprint file for the given file or data frame with columns [compound, smiles]
#' return based on the read parameter
#'   read == 'none' -> return NULL
#'   read == 'fname' -> return the filename of the generated fingerprint file
#'   read == 'compound' -> return just the list of compounds in the fingerprint file
#'   read == 'full' -> read in the full fingerprint file
#' @export
create.fp <- function(
	smi,
	read=c("none", "fname", "compound", "full"),
	compound="compound", smiles="smiles",
	output_fname=NULL){

	read <- match.arg(read)

	if(is.character(smi)){
		if(!str_detect(smi, ".smi$")){
			cat("WARNING: .smi file '", smi, "', does not end in '.smi'\n")
		}
		smi_fname <- smi
	} else if(is.data.frame(smi)){
		smi_fname <- paste0(tempfile(), ".smi")
		cat("writing .smi to -> '", smi_fname, "' ... ", sep="")
		write.smi(smi, smi_fname, compound=compound, smiles=smiles)
		cat("DONE\n")
	} else{
		cat("ERROR: unrecognized smi class: ", paste0(class, collapse=", "), "\n")
	}

	system(
		paste(
			"sea-molecule-fingerprint",
			smi_fname))

	if(is.null(output_fname)){
		fp_fname <- str_replace(smi_fname, fixed(".smi"), ".fp")
	} else {
		fp_fname = output_fname
	}

	if(read == "none"){
		return()
	} else if(read=="fname"){
		return(fp_fname)
	} else if(read=="compound"){
		return(read.fp.compound(fp_fname, compound=compound))
	} else if(read=="full") {
		return(read.fp(fp_fname, compound=compound))
	} else {
		stop("ERROR: read type not recognized.\n")
	}
}

#' @export
read.fp <- function(fname, fingerprint="fingerprint", compound="compound"){
	require(sqldf)
	if(!is.character(fname) | !file.exists(fname)){
		stop(paste("ERROR: The SEA .fp file, '", fname, "', does not exist.", sep=""))
	}
	fp <- read.csv.sql(fname, "SELECT * FROM file", header=F, sep=";")
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

	if(!str_detect(fname, ".fp$")){
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

	library(data.table)
	scores <- fread(fname, sep=",", header=T, data.table=F)

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

	setnames(
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
