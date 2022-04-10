#################################### SELECT VARIABLES #################################

select_variables <- function(dat, selected_variables, fixed_var){
  # Selected variables for the model
  if (length(selected_variables) > 0){
    dat <- dat[, c(fixed_var, selected_variables[selected_variables %in% colnames(dat)]), with=F];
  }
  return(list(dat=dat));
}

######################### REMOVE CONSTANT ##########################################
constant_variables <- function(dat){
  unique_values <- apply(dat,2,function(x){length(unique(x))})
  constant_variables <- names(unique_values)[unique_values == 1];
  return(constant_variables);
}

###############################  REMOVE NOT INFORMED VARIABLES ############################## 

not_informed <- function(dat, threshold_nas=0.4){
  ratio_nas <- apply(dat,2,function(x){sum(is.na(x) | x=="NA")/length(x)});
  not_informed_var <- names(ratio_nas)[ratio_nas>threshold_nas];
  return(not_informed_var);
}

######################### FILL_VALUES ##########################################
# Cast categorical missing values to "NA"
missing_to_NA_category <- function(x){
  if (sum(is.na(x)) > 0){
    levels(x) <- c(levels(x), "NA");
    x[is.na(x)] <- "NA"
  } 
  return(x);
}

# Fill missing values
fill_values <- function(dat, fill){
  no_categorical <- colnames(dat)[sapply(dat, class) == "numeric"];
  if (sum(sapply(dat[, no_categorical, with = F],function(x){sum(is.na(x))})) > 0){
    # Get original order;
    dat$ord<-1:nrow(dat);
    
    # Generalize dataset
    if (sum(colnames(dat)=="dataset")==0){
      dat$dataset <- "dummy";
    }
    
    # Fill missing values with mean
    if (fill=="mean"){
      m <- apply(dat[dataset!="original_train", which(sapply(dat, class) == "numeric"), with = F],
                 2, mean, na.rm=TRUE);
      m[is.na(m)] <- 0;
      dat[, no_categorical] <- data.table(apply(rbind(t(m), 
                                                      dat[, no_categorical, with = F]), 
                                                2, replace));
      
      # Fill missing using mice package
    } else if (fill=="mice"){
      n_m<-5;
      for (iter in 1:2){
        try(imputed<-mice(dat[,-"target",with=F],maxit=5,m=n_m));
        if (exists("imputed")){
          res<-complete(imputed,1)[,no_categorical];
          for (i in 2:n_m){
            res<-res+complete(imputed,i)[,no_categorical];
          }
          dat[,no_categorical]<-res/n_m;
          break;
        }
      }
      
      # Fill persisting missing values with mean
      m<-apply(dat[,no_categorical,with=F],2,mean,na.rm=TRUE);
      dat[,no_categorical] <- data.table(apply(rbind(t(m), 
                                                     dat[, no_categorical, with = F]), 
                                               2, replace));
      
      # Fill missing using Amelia package
    } else if (fill=="amelia"){
      nominal<-setdiff(colnames(dat),c("target","index",no_categorical))
      idvars<-c("index","order",colnames(dat)[sapply(dat,class)=="character" & sapply(dat,function(x){sum(is.na(x))}) > 0]);
      noms<-setdiff(colnames(dat)[sapply(dat,class)=="character"],idvars)
      
      dat[,no_categorical]<-amelia(dat[,-"target",with=F],noms=noms,
                                   idvars=idvars,m = 1,parallel="multicore",cl=4)$imputations[[1]][,no_categorical];
      
      # Fill persisting missing values with mean
      m<-apply(dat[,no_categorical,with=F],2,mean,na.rm=TRUE);
      dat[,no_categorical]<-apply(rbind(t(m),dat[,no_categorical,with=F]),2,replace);
    }
    dat<-dat[order(ord)];
    dat<-dat[,-"ord",with=F];
  }
  
  # Generalize dataset
  if (sum(dat$dataset=="dummy")==nrow(dat)){
    dat <- dat[,-"dataset",with=F]
  }
  
  return(dat);
}

# Compute "untouchable" categories (the ones with anomaly over the target)
untouchable_values_f <- function (name, x, target, min_transactions=2, anomaly_factor=1){
  # Format values
  values <- as.character(x);
  
  # List of categories
  unique_values <- unique(values);
  
  # Get mean of target per category
  dt <- data.table(target = target, value = values);
  mean_target <- dt[, mean(target), by = "value"];
  var_names <- mean_target$value;
  mean_target <- mean_target$V1;
  names(mean_target) <- var_names;
  
  # Get global frequency (classification) or target value (regression)
  global_value <- dt[, .N, by = "value"];
  var_names <- global_value$value;
  global_value <- global_value$N;
  names(global_value) <- var_names;
  
  # Get statistics of normal behaviour
  m <- mean(mean_target);
  s <- sd(mean_target);
  
  # Set untouchable variables as the one with anomaly fraud frequency
  untouchable <- names(mean_target)[((mean_target > (m + anomaly_factor*s)) | (mean_target < (m - anomaly_factor*s))) & global_value >= min_transactions];
  
  if (length(untouchable)>0){
    return(data.table(name=name,untouchable=list(untouchable)))
  } else {
    return(data.table(name=name,untouchable=NA));
  }
}

######################### GROUP CATEGORIES LIST ##########################################
# Get list of categories to group
group_categories_list <- function (name, dat_train, threshold_variables,min_size_others=2,max_n_categories=10,
                                   untouchable_values=data.table(), diag_vars, proc_vars){
  
  grouped_values<-NA;
  if (name %in% diag_vars){
    x <- c()
    for (diag_name in diag_vars){
      x <- c(x, as.character(dat_train[[diag_name]]))
    }
  } else if (name %in% proc_vars){
    x <- c()
    for (proc_name in proc_vars){
      x <- c(x, as.character(dat_train[[proc_name]]))
    }
  }  else {
    x <- dat_train[[name]]
  }
  x <- as.character(x);
  ratios <- table(x)/length(x)
  ind_others<- (ratios <  threshold_variables);
  if (length(ind_others) >= max_n_categories & sum(ind_others) > min_size_others){
    grouped_values<-names(ratios)[ind_others];
    if(length(untouchable_values)>0){
      grouped_values<-grouped_values[!(grouped_values %in% unlist(untouchable_values$untouchable[untouchable_values$name==name]))]
    }
  }
  
  if (!is.na(grouped_values)[1] & length(grouped_values)>0){
    return(data.table(name=name,grouped_values=list(grouped_values)))
  } else {
    return(data.table(name=name,grouped_values=NA));
  }
}


######################### GROUP CATEGORIES ##########################################
# Group values with few instances
group_values <- function (name,values_to_group,x, diag_vars, proc_vars){
  
  # Reformat values
  x<-as.character(x);
  
  # Group categories
  if (name %in% c(diag_vars, proc_vars)){
    x[x%in% values_to_group]<-substr(x[x%in% values_to_group], 1, 3);
  } else {
    x[x%in% values_to_group]<-"others";
  }
  
  
  return(x);
}

######################### NEW TO OTHERS ##########################################
new_to_others<- function(name, dat, dat_train, diag_vars, proc_vars){
  # Reformat values
  x<-as.character(dat[[name]]);
  
  if (name %in% diag_vars){
    x_train <- c()
    for (diag_name in diag_vars){
      x_train <- c(x_train, as.character(dat_train[[diag_name]]))
    }
  } else if (name %in% proc_vars){
    x_train <- c()
    for (proc_name in proc_vars){
      x_train <- c(x_train, as.character(dat_train[[proc_name]]))
    }
  }  else {
    x_train <- as.character(dat_train[[name]])
  }
  
  train_values<-unique(x_train);
  if (name %in% c(diag_vars, proc_vars)){
    x[!(x %in% train_values)] <- substr(x[!(x %in% train_values)], 1, 3);
  } else {
    x[!(x %in% train_values)] <- "others";
  }
  return(x);
}

######################### NEW TO OTHERS ##########################################
new_to_others_list <- function(name, dat, dat_train, diag_vars, proc_vars){
  
  # Reformat values
  x<-as.character(dat[[name]]);
  
  if (name %in% diag_vars){
    x_train <- c()
    for (diag_name in diag_vars){
      x_train <- c(x_train, as.character(dat_train[[diag_name]]))
    }
  } else if (name %in% proc_vars){
    x_train <- c()
    for (proc_name in proc_vars){
      x_train <- c(x_train, as.character(dat_train[[proc_name]]))
    }
  }  else {
    x_train <- as.character(dat_train[[name]])
  }
  
  train_values<-unique(x_train);
  
  new_values <- unique(x[!(x %in% train_values)])
  
  return(data.table(new_values = list(new_values)));
}

######################### GROUP CATEGORIES ##########################################
# Group categories
group_categories <- function(categorical, dat, dat_train, threshold_variables, min_size_others, max_n_categories, untouchable_values, 
                             diag_vars, proc_vars){
  
  # Get list of categories to group
  grouped_values_list<-data.table(t(sapply(categorical, 
                                           function(x){group_categories_list(x, dat_train, threshold_variables,
                                                                             min_size_others,max_n_categories,untouchable_values,
                                                                             diag_vars[-1], proc_vars[-1])})),
                                  stringsAsFactors = TRUE);
  # Apply categories grouping                                                                                                                                                      min_size_others,max_n_categories,untouchable_values)})));
  dat[,categorical]<-data.table(sapply(categorical, 
                                       function(x){group_values(name = x,
                                                                values_to_group = unlist(grouped_values_list[name==x]$grouped_values),
                                                                x = dat[[x]],
                                                                diag_vars, 
                                                                proc_vars)}),
                                stringsAsFactors = TRUE);
  
  # # Get values of val and test not in train
  new_values_list<-data.table(sapply(categorical,function(x){new_to_others_list(x, dat, dat_train, diag_vars, proc_vars)}),
                              stringsAsFactors = TRUE);
  
  # # Set values of val and test not in train to others
  dat[,categorical]<-data.table(sapply(categorical,function(x){new_to_others(x, dat, dat_train, diag_vars, proc_vars)}),
                                stringsAsFactors = TRUE);
  
  # Merge lists
  for (i in 1:nrow(grouped_values_list)){
    grouped_values_list$grouped_values[i] <- list(c(unlist(grouped_values_list$grouped_values[i]), unlist(new_values_list$new_values[i])))
  }
  
  
  return(list(dat=dat,grouped_values_list=grouped_values_list));
  
}


######################### ONE HOT DIAG PROC ##########################################

one_hit_diag_proc <- function(dat, diag_vars, proc_vars){
  
  res <- foreach (vars = list(diag_vars, proc_vars)) %dopar%{
    # Set categories
    categories <- c()
    for (var in vars){
      categories <- c(categories, as.character(unique(dat[[var]])))
    }
    categories <- unique(categories)
    
    # Free memory
    gc();
    
    # Fill diag data
    dat_flags <- matrix(0, nrow = nrow(dat), ncol = length(categories))
    dat_ranks <- matrix(0, nrow = nrow(dat), ncol = length(categories))
    colnames(dat_flags) <- categories
    colnames(dat_ranks) <- categories
    for (var in vars){
      x <- dat[[var]]
      for (value in unique(x)){
        dat_flags[, value] <- dat_flags[, value] + as.numeric(x == value)
        dat_ranks[x == value, value] <- unlist(sapply(dat_ranks[x == value, value], function(x){x[x!=0] <- sapply(x[x!=0], min, 1 + which(vars == var)); x[x == 0] <- 1 + which(vars == var); return(x)}))
      }
    }
    dat_flags[dat_flags == 0] <- length(vars) + 2;
    dat_ranks[dat_ranks == 0] <- length(vars) + 2;
    colnames(dat_flags) <- paste0(gsub("[0-9]", "", vars[1]), "sec_flag_", colnames(dat_flags))
    colnames(dat_ranks) <- paste0(gsub("[0-9]", "", vars[1]), "sec_rank_", colnames(dat_ranks))
    
    # Free memory
    gc();
    
    # Set colnames
    dat_flags <- data.table(dat_flags)
    dat_ranks <- data.table(dat_ranks)
    
    # Return output
    list(dat_flags = dat_flags, dat_ranks = dat_ranks)
  }
  
  # Replace diag data
  dat <- dat[, setdiff(colnames(dat), c(diag_vars, proc_vars)), with = F]
  dat <- cbind(dat, 
               res[[1]]$dat_flags, res[[1]]$dat_ranks,
               res[[2]]$dat_flags, res[[2]]$dat_ranks)
  
  # Free memory
  gc();
  
  
  return(dat);
  
}


######################### ONE HOT ENCODING ##########################################
onehot <- function(x){
  values<-unique(x);
  
  ret <- matrix(0,nrow = length(x),ncol=length(values));
  
  for (i in 1:length(values)){
    ret[, i] = as.numeric(x==values[i]);
  }
  
  colnames(ret)<-values;
  return(ret);
}

one_hot_encoding <- function(dat,categorical,no_categorical,fixed_var){
  
  cat_dat <- as.matrix(dat[, categorical, with = F]);
  vars <- sample(categorical, length(categorical));
  ret <- foreach (c = vars, .combine = cbind) %dopar%{
    inner_ret <- onehot(cat_dat[,c]);
    colnames(inner_ret) <- paste0(c, "_", colnames(inner_ret));
    inner_ret;
  }
  ret <- data.table(ret);
  
  if (!is.na(fixed_var)[1]){
    dat <- data.table(dat[,fixed_var,with=F],ret, dat[,no_categorical,with=F]);
  } else {
    dat <- data.table(ret, dat[,no_categorical,with=F]);
  }
  return(dat);
}

################################ GET STATISTIC #################################################
# Get statistic of a vector after removing tail values
get_statistic <- function(x, fun, na.rm = FALSE){
  if (length(unique(remove_tails(x))) > 1){
    x <- remove_tails(x);
  }
  ret <- fun(x, na.rm = na.rm);
  return(ret);
}


################################ Remove outliers #################################################
# Remove tail values ('outliers') from a vector
remove_tails <- function(x, thresholds = c(0.01, 0.99)){
  percentiles <- quantile(x, thresholds);
  x[x < percentiles[1]] <- percentiles[1];
  x[x > percentiles[2]] <- percentiles[2];
  return(x);
}


#################################### TIPIFY COLUMN #################################
tipify_column <- function(column, m, s){
  return((column-m)/s);
}

#################################### TIPIFY #################################
tipify <- function(dat,m,s){
  
  # Convert to data.table
  if (class(dat)[1]!="data.table"){
    dat<-data.table(dat);
  }
  
  # Remove constant variables
  constant_var<-which(s==0);
  if (length(constant_var)>0){
    dat<-dat[,-constant_var,with=F];
    m<-m[-constant_var];
    s<-s[-constant_var];
  }
  print(sprintf("%s constant variables removed",length(constant_var)));
  
  #Free memory
  gc();
  
  # Cast to matrix
  dat <- as.matrix(dat)
  
  #Free memory
  gc();
  
  # Tipify
  for (i in 1:ncol(dat)){
    dat[, i] <- (dat[, i] - m[i])/s[i];
  }
  
  #Free memory
  gc();
  
  # Cast back to data.table
  dat <- as.data.table(dat)
  
  #Free memory
  gc();
  
  # Set back column names
  colnames(dat) <- names(m)
  
  #Free memory
  gc();
  
  return(list(dat=dat,m=m,s=s));
}

#################################### UNTIPIFY #################################
untipify <- function(dat,m,s){
  
  # Convert to data.table
  if (class(dat)[1]!="data.table"){
    dat<-data.table(dat);
  }
  
  # Remove constant variables
  constant_var<-which(s==0);
  if (length(constant_var)>0){
    dat<-dat[,-constant_var,with=F];
    m<-m[-constant_var];
    s<-s[-constant_var];
  }
  print(sprintf("%s constant variables removed",length(constant_var)));
  
  dat <- as.matrix(dat);
  for (i in 1:ncol(dat)){
    dat[,i] <-  dat[,i]*s[i] + m[i];
  }
  dat <- data.table(dat)
  
  
  colnames(dat)<-names(m);
  
  return(list(dat=dat,m=m,s=s));
}

#################################### IRRELEVANT VARIABLES #################################

remove_irrelevant<-function(correlations,irrelevant_threshold){
  # Irrelevant variables
  relevance<-correlations[1,-1];
  irrelevant_variables<-names(relevance)[is.na(relevance) | relevance<irrelevant_threshold];
  return(irrelevant_variables);
}

#################################### REDUNDANT VARIABLES #################################
remove_redundant <- function(correlations,redundant_threshold){
  redundancy<-apply(correlations,2,function(x){which(x>redundant_threshold)});
  redundancy<-redundancy[which(sapply(redundancy,length)>1)]
  
  redundant_variables<-c();
  for (i in redundancy){
    imp<-sort(correlations[1,i],decreasing = TRUE);
    redundant_variables<-c(redundant_variables,names(imp)[2:length(i)])
  }
  redundant_variables<-unique(redundant_variables);
} 

#################################### IMPORTANT VARIABLES #################################
select_important<-function(dat,n_perc_varimp,priority_class=2,y=NA){
  if (is.na(y[1])){
    varimp<-filterVarImp(x=dat[,-fixed_var,with=F],y=dat$target,nonpara=TRUE);
  } else {
    varimp<-filterVarImp(x=dat,y=y,nonpara=TRUE);
  }
  priority_class <- min(priority_class, ncol(varimp))
  varimp<-data.table(variable=rownames(varimp),imp=varimp[,priority_class])
  varimp<-varimp[order(-imp)];
  
  if (n_perc_varimp > 0){
    selected<-varimp$variable[1:ceiling(n_perc_varimp*nrow(varimp))];
  } else {
    selected <-NA
  }
  
  return(list(selected = selected, varimp = varimp));
}

#################################### PCA #################################
pca <- function(dat,dat_train,threshold_variance, fixed_var){
  # PCA
  prin_comp <- prcomp(dat_train[,-fixed_var,with=F],center = FALSE,scale. = FALSE); 
  
  # choose number of variables
  prop_varex <- summary(prin_comp)$importance[3,];
  number_of_variables <- min(which(prop_varex > threshold_variance))
  
  # Apply PCA
  dat<-data.table(dat[,fixed_var,with=F],
                  data.table(predict(prin_comp, newdata = dat[,-fixed_var,with=F]))[,1:number_of_variables,with=F]);
  
  # Get mean and standard deviation for train data (excluding target) 
  m_tip<-apply(dat[dataset=="train",-fixed_var,with=F],2,mean);
  s_tip<-apply(dat[dataset=="train",-fixed_var,with=F],2,sd);
  
  # Tipify data
  ret<-tipify(dat[,-fixed_var,with=F],m_tip,s_tip);
  dat<-data.table(dat[,fixed_var,with=F],ret$dat);
  m_tip<-ret$m;
  s_tip<-ret$s;
  
  return(list(dat=dat,prin_comp=prin_comp,m_tip=m_tip,s_tip=s_tip));
}