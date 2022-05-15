#################################### Load Libraries ####################################
library(data.table);
library(treemap);
library(foreach)
library(doParallel)
library(caret)


#################################### Global Paths #################################### 
scripts_path <- "/bigdata/utiles_jesusprada/TCGA/scripts/R/"
input_path <- "/bigdata/utiles_jesusprada/TCGA/data/"
output_path <- "/bigdata/utiles_jesusprada/TCGA/data/"


#################################### Source Auxiliary Scripts #################################### 
source(file.path(scripts_path, "aux.R"))
source(file.path(scripts_path, "clean_data.R"))

#################################### Global Variables #################################### 
perc_splitting <- c(0.7, 0.15, 0.15)
years_window <- 3
fixed_var <- c("target", "dataset", "index")




#################################### Read Data #################################### 
dat <- readRDS(file.path(input_path, "raw_data.rds"))
dat$index <- 1:nrow(dat);
dicc_kegg_brite_id <- fread(file.path(input_path, "dicc", "1_kegg_brite_id.csv"), header = F, colClasses = "character")
dicc_kegg_gene_id <- fread(file.path(input_path, "dicc", "2_kegg_gene_id.csv"), header = F)
dicc_gene_to_hugo <- fread(file.path(input_path, "dicc", "3_gene_to_hugo.csv"), header = F)
dicc_ensemble_to_hugo <- fread(file.path(input_path, "dicc", "4_ensemble_to_hugo.csv"))


#################################### Filter cancer samples ###############################
dat <- dat[dat$definition == "Primary solid Tumor"]
dat <- dat[, .SD[1], by = "patient"]

#################################### Set target #################################### 
dat$days_to_collection[is.na(dat$days_to_collection)] <- 0;
dat$days_to_last_follow_up[is.na(dat$days_to_last_follow_up)] <- dat$days_to_collection[is.na(dat$days_to_last_follow_up)] + years_window*365;
dat$days_to_death <- dat$days_to_death - dat$days_to_collection
dat$days_to_last_follow_up <- dat$days_to_last_follow_up - dat$days_to_collection
dat$target <- as.numeric(dat$days_to_death > 0 & dat$days_to_death <= years_window*365)
dat$target[is.na(dat$target)] <- 0;
dat$days_to_death[is.na(dat$days_to_death)] <- 0;
dat <- dat[!(target == 0 & days_to_last_follow_up < years_window*365)]
dat <- dat[!(days_to_death < 0)]
dat <- dat[!(target == 0 & vital_status == "Dead" & dat$days_to_death == 0)]
table(dat$target)
table(dat$target)/nrow(dat)



#################################### Data Splitting ####################################
new_dat <- data.table();
for (cancer in unique(dat$project_id)){
  for (target_value in unique(dat$target)){
    subdat <- dat[project_id == cancer & target == target_value]
    if (length(subdat) > 0){
      train_index <- sample(1:nrow(subdat), 
                            floor(nrow(subdat)*perc_splitting[1]))
      index <- setdiff(1:nrow(subdat), train_index)
      val_index <- sample(index, 
                          floor(nrow(subdat)*perc_splitting[2]))
      test_index <- setdiff(index, val_index)
      subdat$dataset <- "train"
      subdat$dataset[val_index] <- "val"
      subdat$dataset[test_index] <- "test"
      new_dat <- rbind(new_dat, subdat);
    }
  }
}
dat <- new_dat

# Free memory
rm(new_dat)
gc()


#################################### Remove unuseful variables ####################################
unuseful_variables <- c(colnames(dat)[grepl("paper_", colnames(dat))],
                        "patient",
                        "sample",
                        "shortLetterCode",
                        "definition",
                        "sample_submitter_id",
                        "sample_type_id",
                        "tumor_descriptor",
                        "sample_id",
                        "sample_type",
                        "composition",
                        "state",
                        "pathology_report_uuid",
                        "submitter_id",
                        "is_ffpe",
                        "tissue_type",
                        "days_to_diagnosis",
                        "last_known_disease_status",
                        "classification_of_tumor",
                        "diagnosis_id",
                        "tumor_grade",
                        "progression_or_recurrence",
                        "exposure_id",
                        "year_of_birth",
                        "demographic_id",
                        "bcr_patient_barcode",
                        "releasable",
                        "released",
                        "preservation_method",
                        "sample.aux",
                        "igcccg_stage"
                        )
dat <- dat[, setdiff(colnames(dat), unuseful_variables), with =F]

# Free memory
gc();

#################################### Remove cheating variables ####################################
cheating_variables <- c("treatments",
                        "days_to_sample_procurement",
                        "days_to_last_follow_up",
                        "days_to_death",
                        "year_of_death",
                        "dead_in_year",
                        "vital_status")
dat <- dat[, setdiff(colnames(dat), cheating_variables), with =F]

# Free memory
gc();

################################# Set classes #########################################
dat$primary_site <- sapply(dat$primary_site, '[', 1)
dat$disease_type <- sapply(dat$disease_type, '[', 1)
dat$year_of_diagnosis <- as.factor(dat$year_of_diagnosis)

binary_vars <- colnames(dat)[sapply(dat, function(x){length(unique(x)) <= 2}) & sapply(dat, class) == "character"]
dat[, binary_vars] <- data.table(sapply(dat[, binary_vars, with = F], function(x){as.numeric(x == unique(x)[1])}))



# Character to factor
dat[, which(sapply(dat, class) == "character")] <- data.table(sapply(dat[, sapply(dat, class) == "character", with = F], as.factor),
                                                              stringsAsFactors = TRUE);

# Integer to numeric
dat[, which(sapply(dat, class) == "integer")] <- data.table(sapply(dat[, sapply(dat, class) == "integer", with = F], as.numeric));

#################################### Clean Data ####################################
dat_original <- dat
ret_clean <- clean_data(dat,
                        selected_variables = c(),
                        factor_datos = 0,
                        over_factor = 0,
                        k = 5,
                        normalization = "tip",
                        fill = "mean", 
                        threshold_variance = 0,
                        redundant_threshold = 0.99,
                        irrelevant_threshold = 0.01,
                        n_perc_varimp = 0,
                        seed=TRUE,
                        threshold_nas = 0.8, 
                        fixed_var);
dat <- ret_clean$dat;

# Free memory
ret_clean[["dat"]] <- NULL;
gc();

# Get top 1000 most relevant genes
top_1000_genes <- ret_clean$varimp$variable[grepl("\\|", ret_clean$varimp$variable)][1:1000]
gen_dat <- dat[,top_1000_genes, with = F]
dat <- cbind(dat[, colnames(dat)[!grepl("\\|",colnames(dat))] , with = F], gen_dat)

# Save Output 
saveRDS(dat, file.path(output_path, "clean_data.rds"))
write.table(dat, file.path(output_path, "clean_data.csv"), sep = ";", col.names = T, row.names = F)
saveRDS(ret_clean, file.path(output_path, "ret_clean.rds"))

#################################### Cast to Images ####################################
indexes <- dat$index
# Get gen matrix
top_1000_genes <- ret_clean$varimp$variable[grepl("\\|", ret_clean$varimp$variable)][1:1000]
gen_dat <- dat[,top_1000_genes, with = F]

dicc_kegg_brite_id$V1 <- NULL
dicc_kegg_brite_id$V1 <- dicc_kegg_brite_id$V2
dicc_kegg_brite_id$V2 <- NULL
dicc_kegg_gene_id$V1 <- gsub(":", "", gsub("[a-z]","", dicc_kegg_gene_id$V1))
dicc_kegg <- merge(dicc_kegg_gene_id, dicc_kegg_brite_id, by = "V1")
setnames(dicc_gene_to_hugo, c("V1", "V2"), c("V2", "V6"))
dicc_kegg <- merge(dicc_kegg, dicc_gene_to_hugo, by = "V2")
dicc_kegg$name <- sapply(strsplit(dicc_kegg$V6, "; "), '[', 2)

registerDoParallel(cores = detectCores())

foreach (i = 1:nrow(gen_dat))%dopar%{
  index <- indexes[i]
  patient <- gen_dat[i,]
  patient <- data.table(symbol = colnames(patient), value = as.numeric(patient))
  patient$symbol <- sapply(strsplit(patient$symbol, "\\|"), '[', 1)
  patient <- merge(patient, dicc_ensemble_to_hugo[, c("symbol", "name"), with = F], by = "symbol")
  patient <- merge(patient, dicc_kegg[, c("name", "V3", "V4", "V5"), with = F], by = "name")
  setnames(patient, c("V3", "V4", "V5", "name"), c("level_1", "level_2", "level_3", "level_4"))
  patient$size <- 1;
  
  #Create image
  png(filename=file.path(output_path, "images", sprintf("%s.png", index)),width=1024, height=1024)
  treemap(patient, 
          index=c("level_1", "level_2", "level_3", "level_4"), 
          vSize="size", 
          vColor="value", 
          type="value", 
          palette="RdYlBu",
          fontsize.labels= 0,title = "", title.legend = "", fontsize.title = 0, fontsize.legend = 0, position.legend = "none", 
          )
  dev.off()
}
