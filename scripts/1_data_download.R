#################################### Load Libraries ####################################
library(TCGAbiolinks)
library(DESeq2)
library(DT)
library(data.table)

#################################### Global Variables #################################### 
input_path <- "/home/jesusprada2/cancer/data"
output_path <- "/home/jesusprada2/cancer/data"
cancer_types <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-DLBC","TCGA-ESCA","TCGA-GBM",
                  "TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC",
                  "TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD",
                  "TCGA-TGCT","TCGA-THCA","TCGA-THYM","TCGA-UCEC","TCGA-UCS","TCGA-UVM")

dat <- data.table()
for (cancer in cancer_types){
  
  #################################### Download Data #################################### 
  print(sprintf("Download %s data...", cancer))
  query <- GDCquery(project = cancer,
                    legacy = T,
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq",
                    file.type = "results",
                    experimental.strategy  = "RNA-Seq",
                    sample.type = c("Primary Tumor","Solid Tissue Normal"),
                    access = "open")
  GDCdownload(query, 
              method = "api", 
              files.per.chunk = 1000, 
              directory = input_path)
  
  #################################### Read Data #################################### 
  data <- GDCprepare(query, directory =  input_path)
  
  
  #################################### Create Dataset ####################################
  # Get information about dataset
  info_dat <- as.data.frame(colData(data))
  
  # Extract gene count matrix
  subdat <- data.frame(assay(data))
  
  # Scale by patient
  row_names <- rownames(subdat)
  subdat <- data.frame(sapply(subdat, function(x){x_scaled <- x / sum(x, na.rm = T); return(x_scaled)}))
  rownames(subdat) <- row_names 
  
  # Format dataset into ML format
  patients <- colnames(subdat)
  subdat <- data.table(t(subdat))
  rownames(subdat) <- patients
  subdat$barcode <- gsub("\\.", "-", patients)
  
  # Merge gene expression with patient information
  subdat <- merge(subdat, info_dat, by = "barcode")
  subdat$barcode <- NULL
  
  # Combine results
  dat <- rbind(dat, subdat, fill = TRUE)
}


#################################### EDA #################################### 
# Check dataset composition
table(dat$project_id)

table(dat$sample_type)
table(dat$project_id, dat$sample_type)

table(dat$vital_status)
table(dat$project_id, dat$vital_status)

# Check for differences in gene expression
healthy_summary <- sapply(dat[sample_type == "Solid Tissue Normal"], mean)
tumor_summary <- sapply(dat[target != "Solid Tissue Normal"], mean)
summary <- data.table(gene = names(healthy_summary), healthy = healthy_summary, tumor = tumor_summary)
summary$ratio <- summary$tumor / summary$healthy
summary$ratio[summary$ratio == Inf] <- NA
summary$diff <- abs(summary$tumor - summary$healthy)


#################################### Save Output #################################### 
saveRDS(dat, file.path(output_path, "raw_data.rds"))
write.table(dat, file.path(output_path, "raw_data.csv"), sep = "|", col.names = T, row.names = F)
