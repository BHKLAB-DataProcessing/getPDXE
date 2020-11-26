library(Xeva)
library(Biobase)

setwd("~/Desktop/PDXE_xevaSet-main")
model <- read.csv("data/raw_data/model_info.csv")
experiment <- read.csv("data/raw_data/expriment.csv")
drugInfo <- read.csv("data/raw_data/drug_info.csv")

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
modToBiobaseMap <- read.csv("data/raw_data/modToBiobaseMap.csv")

##-----read mol data ---
RNASeq <- readRDS("data/raw_data/molProf_RNASeq.rds")
cnv <- readRDS("data/raw_data/molProf_cnv.rds")
mutation <- readRDS("data/raw_data/molProf_mutation.rds")

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


##=========subset by tissue type ====
mi <- modelInfo(pdxe)
for(tissue in c("Breast Cancer", "Colorectal Cancer", "Cutaneous Melanoma",
                "Gastric Cancer", "Non-small Cell Lung Carcinoma",
                "Pancreatic Ductal Carcinoma"))
{
  pdxe.tissue <- subsetXeva(pdxe, ids= mi$model.id[mi$tissue.name==tissue],
                            id.name = "model.id")
  saveRDS(pdxe.tissue, sprintf("data/PDXE_xevaSet/PDXE_%s_XevaSet.rds", gsub(" ", "_", tissue)))
}

