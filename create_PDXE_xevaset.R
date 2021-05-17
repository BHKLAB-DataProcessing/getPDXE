library(Xeva)
library(Biobase)
library(biocompute)

model <- read.csv("/pfs/getPDXE/model_info.csv")
experiment <- read.csv("/pfs/getPDXE/expriment.csv")
drugInfo <- read.csv("/pfs/getPDXE/drug_info.csv")

control.drug <- "untreated"
drugs <- unique(model$drug)
drugs <- drugs[drugs!=control.drug]
##------create expriment design list

expDesign <- list()
for(p in unique(model$patient.id))
{
  for(d in drugs)
  {
    bt <- list(batch.name = sprintf("%s.%s", p,d),
               control = model$model.id[model$patient.id==p & model$drug==control.drug],
               treatment=model$model.id[model$patient.id==p & model$drug==d])
    if(length(bt$control)>0 | length(bt$treatment)>0)
    { expDesign[[bt$batch.name]] <- bt }
  }
}

###----------------
modToBiobaseMap <- read.csv("/pfs/getPDXE/modToBiobaseMap.csv")

##-----read mol data ---
RNASeq <- readRDS("/pfs/getPDXE/molProf_RNASeq.rds")
cnv <- readRDS("/pfs/getPDXE/molProf_cnv.rds")
mutation <- readRDS("/pfs/getPDXE/molProf_mutation.rds")

##-----check if ids are matching ----

RNASeq <- RNASeq[, sampleNames(RNASeq)%in%modToBiobaseMap$biobase.id]
cnv <- cnv[, sampleNames(cnv)%in%modToBiobaseMap$biobase.id]
mutation <- mutation[, sampleNames(mutation)%in%modToBiobaseMap$biobase.id]


##==== create XevaSet ====
pdxe = createXevaSet(name="PDXE xevaSet",
                     model = model,
                            drug = drugInfo,
                            experiment = experiment,
                            expDesign = expDesign,
        molecularProfiles=list(RNASeq = RNASeq,
                               mutation=mutation,
                               cnv=cnv),
                            modToBiobaseMap = modToBiobaseMap)

##===== set response =====
for(res in c("mRECIST", "slope", "AUC", "angle", "abc", "TGI"))
{
  pdxe <- setResponse(pdxe, res.measure = res, verbose=TRUE)
}

saveRDS(pdxe, "/pfs/out/Xeva_PDXE.rds")

### CREATE BIOCOMPUTE OBJECT ###

library(biocompute)
###########################
#####Provenance Domain#####
###########################

#Created and modified dates
#Sys.setenv(TZ = "EST")
created <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
modified <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")

#Contributions
contributors <- data.frame(
  "name" = c("Anthony Mammoliti", "Petr Smirnov", "Benjamin Haibe-Kains"),
  "affiliation" = c(rep("University Health Network", 3)),
  "email" = c("anthony.mammoliti@uhnresearch.ca", "petr.smirnov@utoronto.ca", "Benjamin.Haibe-Kains@uhnresearch.ca"),
  "contribution" = c("createdBy","createdBy","authoredBy"),
  "orcid" = c(NA,NA,"https://orcid.org/0000-0002-7684-0079"),
  stringsAsFactors = FALSE
)

#License
license <- "https://opensource.org/licenses/Apache-2.0"

#Name of biocompute object
name <- "PDXE"

#Version of biocompute object
bio_version <- "1.0.0"

#Embargo (none)
embargo <- c()

#Derived from and obsolete after (none)
derived_from <- c()
obsolete_after <- c()

#reviewers (none)
review <- c()

#compile domain
provenance <- compose_provenance_v1.3.0(
  name, bio_version, review, derived_from, obsolete_after,
  embargo, created, modified, contributors, license
)
provenance %>% convert_json()


############################
#####Description Domain#####
############################
times_rnaseq <- as.POSIXct("2020-02-19T4:10:32", format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
#Keywords and platform info
keywords <- c("Biomedical", "Xenograph", "Pharmacogenomics")
platform <- c("Pachyderm", "ORCESTRA (orcestra.ca)", "Linux/Ubuntu")

#Metadata for each pipeline step
pipeline_meta <- data.frame(
  "step_number" = c("1"),
  "name" = c("Build data object"),
  "description" = c("(1) Compiled processed RNA-seq, cnv, mutation data; (2) Download of appropriate sample and treatment identifiers; (3) Compile drug sensitivity data; (4) Build data object"),
  "version" = c(1.0),
  stringsAsFactors = FALSE
)

#Inputs for each pipeline step
pipeline_input <- data.frame(
  "step_number" = c("1","1","1","1","1","1","1"),
  "filename" = c("Processed RNA-seq data",
                 "Processed cnv data",
                 "Processed mutation data",
                 "Sample annotation data",
                 "Treatment annotation data",
                 "Processed drug sensitivity data",
                 "Script for data object generation"),
  "uri" = c(
    "/pfs/getPDXE/molProf_RNASeq.rds",
    "/pfs/getPDXE/molProf_cnv.rds",
    "/pfs/getPDXE/molProf_mutation.rds",
    "/pfs/getPDXE/model_info.csv",
    "/pfs/getPDXE/drug_info.csv",
    "/pfs/getPDXE/expriment.csv",
    "https://github.com/BHKLAB-Pachyderm/getPDXE/blob/main/create_PDXE_xevaset.R"
  ),
  "access_time" = c(times_rnaseq,created,created,created,created,created,created),
  stringsAsFactors = FALSE
)


#Outputs for each pipeline step
pipeline_output <- data.frame(
  "step_number" = c("1"),
  "filename" = c("Data object"),
  "uri" = c(
    "/pfs/out/Xeva_PDXE.rds"
  ),
  "access_time" = c(created),
  stringsAsFactors = FALSE
)

#xref (none)
xref <- c()

#pipeline prereq (none)
pipeline_prerequisite <- c()

#compile domain
description <- compose_description_v1.3.0(
  keywords, xref, platform,
  pipeline_meta, pipeline_prerequisite, pipeline_input, pipeline_output
)
description %>% convert_json()


############################
######Execution Domain######
############################

script <- c()
script_driver <- c()

#software/tools and its versions used for data object creation
software_prerequisites <- data.frame(
  "name" = c("Pachyderm", "Docker Image"),
  "version" = c("1.9.3", "v2"),
  "uri" = c(
    "https://www.pachyderm.com", "https://hub.docker.com/r/bhklab/xeva"
  ),
  stringsAsFactors = FALSE
)

software_prerequisites[,"access_time"] <- rep(NA, length(software_prerequisites$name))
software_prerequisites[,"sha1_chksum"] <- rep(NA, length(software_prerequisites$name))

external_data_endpoints <- c()
environment_variables <- c()

execution <- compose_execution_v1.3.0(
  script, script_driver, software_prerequisites, external_data_endpoints, environment_variables
)
execution %>% convert_json()


############################
######Extension Domain######
############################

#repo of scripts/data used
scm_repository <- data.frame("extension_schema"= c("https://github.com/BHKLAB-Pachyderm"))
scm_type <- "git"
scm_commit <- c()
scm_path <- c()
scm_preview <- c()

scm <- compose_scm(scm_repository, scm_type, scm_commit, scm_path, scm_preview)
scm %>% convert_json()

extension <- compose_extension_v1.3.0(scm)
extension %>% convert_json()

############################
######Parametric Domain#####
############################

df_parametric <- data.frame(
  "param" = c(),
  "value" = c(),
  "step" = c(),
  stringsAsFactors = FALSE
)

parametric <- compose_parametric_v1.3.0(df_parametric)
parametric %>% convert_json()



############################
######Usability Domain######
############################

#usability of our data objects
text <- c(
  "Pipeline for creating PDXE data object through ORCESTRA (orcestra.ca), a platform for the reproducible and transparent processing, sharing, and analysis of biomedical data."
)

usability <- compose_usability_v1.3.0(text)
usability %>% convert_json()


######################
######I/O Domain######
######################

input_subdomain <- data.frame(
  c("Processed RNA-seq data",
    "Processed cnv data",
    "Processed mutation data",
    "Sample annotation data",
    "Treatment annotation data",
    "Processed drug sensitivity data",
    "Script for data object generation"),
  "uri" = c(
    "/pfs/getPDXE/molProf_RNASeq.rds",
    "/pfs/getPDXE/molProf_cnv.rds",
    "/pfs/getPDXE/molProf_mutation.rds",
    "/pfs/getPDXE/model_info.csv",
    "/pfs/getPDXE/drug_info.csv",
    "/pfs/getPDXE/expriment.csv",
    "https://github.com/BHKLAB-Pachyderm/getPDXE/blob/main/create_PDXE_xevaset.R"
  ),
  "access_time" = c(times_rnaseq,created,created,created,created,created,created),
  stringsAsFactors = FALSE
)

output_subdomain <- data.frame(
  "mediatype" = c("RDS"),
  "uri" = c(
    "/pfs/out/Xeva_PDXE.rds"
  ),
  "access_time" = c(created,created,created),
  stringsAsFactors = FALSE
)

io <- compose_io_v1.3.0(input_subdomain, output_subdomain)
io %>% convert_json()


########################
######Error Domain######
########################

empirical <- c()
algorithmic <- c()

error <- compose_error(empirical, algorithmic)
error %>% convert_json()


####Retrieve Top Level Fields####
tlf <- compose_tlf_v1.3.0(
  provenance, usability, extension, description,
  execution, parametric, io, error
)
tlf %>% convert_json()


####Complete BCO####

bco <- biocompute::compose_v1.3.0(
  tlf, provenance, usability, extension, description,
  execution, parametric, io, error
)
bco %>% convert_json() %>% export_json("/pfs/out/PDXE_BCO.json") %>% validate_checksum()

