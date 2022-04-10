#################################### Load Libraries ####################################
library(data.table);
library(treemap);


#################################### Source Auxiliary Scripts #################################### 
source("/home/jesusprada2/cancer/scripts/aux.R")
source("/home/jesusprada2/cancer/scripts/clean_data.R")

#################################### Global Variables #################################### 
input_path <- "/home/jesusprada2/cancer/data/"
output_path <- "/home/jesusprada2/cancer/data/"
perc_splitting <- c(0.7, 0.15, 0.15)

fixed_var <- "target"
#################################### Read Data #################################### 
dat <- readRDS(file.path(input_path, "raw_data.rds"))
dicc_kegg_brite_id <- fread(file.path(input_path, "dicc", "1_kegg_brite_id.csv"), header = F, colClasses = "character")
dicc_kegg_gene_id <- fread(file.path(input_path, "dicc", "2_kegg_gene_id.csv"), header = F)
dicc_gene_to_hugo <- fread(file.path(input_path, "dicc", "3_gene_to_hugo.csv"), header = F)
dicc_ensemble_to_hugo <- fread(file.path(input_path, "dicc", "4_ensemble_to_hugo.csv"))

#################################### Set target #################################### 
dat$target <- as.numeric(dat$sample_type != "Solid Tissue Normal")
dat$sample_type <- NULL;

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

#################################### Cast to Images ####################################

# Get gen matrix
last_gen <- which(colnames(dat) == "patient") - 1
gen_dat <- dat[, colnames(dat)[1:last_gen], with = F]


# Fill missing values with 0
gen_dat[is.na(dat)] <- 0

# Remove constant columns
gen_dat <- gen_dat[, sapply(dat, sd) != 0, with = F]

dicc_kegg_brite_id$V1 < NULL
dicc_kegg_brite_id$V1 <- dicc_kegg_brite_id$V2
dicc_kegg_brite_id$V2 <- NULL
dicc_kegg_gene_id$V1 <- gsub(":", "", gsub("[a-z]","", dicc_kegg_gene_id$V1))
dicc_kegg <- merge(dicc_kegg_gene_id, dicc_kegg_brite_id, by = "V1")
setnames(dicc_gene_to_hugo, c("V1", "V2"), c("V2", "V6"))
dicc_kegg <- merge(dicc_kegg, dicc_gene_to_hugo, by = "V2")
dicc_kegg$name <- sapply(strsplit(dicc_kegg$V6, "; "), '[', 2)

for (i in 1:nrow(dat)){
  barcode <- dat$barcode[i]
  patient <- gen_dat[i,]
  patient <- data.table(symbol = colnames(patient), value = as.numeric(patient))
  patient$symbol <- sapply(strsplit(patient$symbol, "\\|"), '[', 1)
  patient <- merge(patient, dicc_ensemble_to_hugo[, c("symbol", "name"), with = F], by = "symbol")
  patient <- merge(patient, dicc_kegg[, c("name", "V3", "V4", "V5"), with = F], by = "name")
  setnames(patient, c("V3", "V4", "V5", "name"), c("level_1", "level_2", "level_3", "level_4"))
  patient$size <- 1;
  
  #Create image
  png(filename=file.path(output_path, "images", sprintf("%s.png", barcode)),width=1024, height=1024)
  treemap(patient, 
          index=c("level_1", "level_2", "level_3", "level_4"), 
          vSize="size", 
          vColor="value", 
          type="value", 
          palette="RdYlBu",
          fontsize.labels= 0)
  dev.off()
}



#################################### Clean Data ####################################
ret_clean <- clean_data(dat,
                        unlist(best_outer$selected_variables),best_outer$factor_datos,best_outer$threshold_variables,
                        best_outer$over_factor,best_outer$k,best_outer$normalization,
                        best_outer$fill, best_outer$anomaly_factor,best_outer$threshold_variance,
                        best_outer$redundant_threshold,best_outer$irrelevant_threshold,
                        best_outer$n_perc_varimp,
                        seed=TRUE,threshold_nas, 
                        global_variables$fixed_var, global_variables$diag_vars, global_variables$proc_vars,
                        target, country);
dat <- ret_clean$dat;

# Free memory
ret_clean[["dat"]] <- NULL;
gc();

# Save Output 
saveRDS(dat, file.path(output_path, "clean_data.rds"))
write.table(dat, file.path(output_path, "clean_data.csv"), sep = "|", col.names = T, row.names = F)



