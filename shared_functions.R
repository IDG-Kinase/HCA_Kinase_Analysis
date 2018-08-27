read_h5_file_to_tidy = function(h5_file) {
  gene_set = data.frame(
    gene_names = h5read(h5_file,"GRCh38/gene_names"),
    ensembl_gene_id = h5read(h5_file,"GRCh38/genes")) %>% 
    mutate(gene_index = 0:(n()-1))
  
  barcodes = data.frame(
    bar_str = h5read(h5_file,"GRCh38/barcodes"))
  
  data_to_barcode = h5read(h5_file,"GRCh38/indptr")
  barcode_read_counts = c(data_to_barcode,NA)- c(NA,data_to_barcode)
  barcode_read_counts = barcode_read_counts[2:(length(barcode_read_counts)-1)]
  
  
  counts_and_gene_index = data.frame(
    counts = h5read(h5_file,"GRCh38/data"),
    gene_index = h5read(h5_file,"GRCh38/indices"),
    barcode = rep(barcodes[,1],barcode_read_counts)
  ) %>% left_join(gene_set)
  
  return(counts_and_gene_index)
}

filter_single_cell = function(tidy_10X, read_total_range = c(-Inf,Inf), unique_genes_range = c(-Inf,Inf)) {
  
  passed_cells = tidy_10X %>% group_by(barcode) %>%
    summarise(read_counts = sum(counts),
              unique_genes = n()) %>%
    filter(between(read_counts, read_total_range[1], read_total_range[2]),
           between(unique_genes, unique_genes_range[1], unique_genes_range[2]))
  
  tidy_10X = tidy_10X %>% filter(barcode %in% passed_cells$barcode)
}

correlate_single_cell_read_counts = function(tidy_10X, min_shared_cells = 100) {
  correlation_sets = list()
  
  gene_list = unique(tidy_10X$gene_names)
  for (gene_1_index in 1:(length(gene_list)-1)) {
    gene_1 = gene_list[gene_1_index]
    gene_1_set = tidy_10X %>% 
      select(gene_names,barcode,counts,class) %>% 
      filter(gene_names == gene_1)
    
    if (dim(gene_1_set)[1] < min_shared_cells) {
      next;
    }
    
    for (gene_2_index in (gene_1_index+1):length(gene_list)) {
      gene_2 = gene_list[gene_2_index]
      if (gene_1 == gene_2) {
        next;
      }
      gene_2_set = tidy_10X %>% 
        select(gene_names,barcode,counts,class) %>% 
        filter(gene_names == gene_2)
      
      overlap_set = inner_join(gene_1_set, gene_2_set, by="barcode")
      
      if (dim(overlap_set)[1] <= min_shared_cells) {
        next;
      }
      
      if (sd(overlap_set$counts.x) == 0 | sd(overlap_set$counts.y) == 0) {
        next;
      }
      
      correlation_sets$gene_1 = c(correlation_sets$gene_1,gene_1)
      correlation_sets$gene_2 = c(correlation_sets$gene_2,gene_2)
      correlation_sets$cell_count = c(correlation_sets$cell_count,dim(overlap_set)[1])
      correlation_sets$value = c(correlation_sets$value,cor(overlap_set$counts.x,overlap_set$counts.y))
    }
  }
  
  correlation_sets = data.frame(correlation_sets)
}

gather_gene_sets = function(tidy_10X, min_shared_cells = 100) {
  
  total_cells = length(unique(tidy_10X$barcode))
  
  barcode_sets = list()
  percent_cells = list()
  for (gene in unique(tidy_10X$symbol)) {
    temp = tidy_10X %>%
      filter(symbol == gene) %>%
      select(barcode)
    
    if (dim(temp)[1] >= min_shared_cells) {
      barcode_sets[[gene]] = temp
      percent_cells[[gene]] = length(barcode_sets[[gene]]$barcode)/total_cells
    }
  }
  
  ##Doublets
  kinase_combinations = combn(names(barcode_sets),2)
  
  passed_doublets = list()
  
  library(progress)
  
  pb <- progress_bar$new(total = dim(kinase_combinations)[2])
  
  for (index in 1:dim(kinase_combinations)[2]) {
    pb$tick()
    gene_1 = as.character(kinase_combinations[1,index])
    gene_2 = as.character(kinase_combinations[2,index])
    
    overlap_set = intersect(barcode_sets[[gene_1]]$barcode,barcode_sets[[gene_2]]$barcode)
    
    if (length(overlap_set) < min_shared_cells) {
      next;
    } else {
      gene_str = glue::collapse(sort(c(gene_1,gene_2)),sep="|")
      passed_doublets[[gene_str]] = overlap_set
    }
    
  }
  
  ##Triplets
  kinase_combinations = combn(names(barcode_sets),3)
  
  passed_sets = list()
  
  
  
  pb <- progress_bar$new(total = dim(kinase_combinations)[2])
  count = 0;
  for (index in 1:dim(kinase_combinations)[2]) {
    pb$tick()
    
    gene_1 = as.character(kinase_combinations[1,index])
    gene_2 = as.character(kinase_combinations[2,index])
    
    gene_pair_1 = glue::collapse(sort(c(gene_1,gene_2)),sep="|")
    if (is.null(passed_doublets[[gene_pair_1]])) {
      next
    }
    
    gene_3 = as.character(kinase_combinations[3,index])
    
    gene_pair_2 = glue::collapse(sort(c(gene_2,gene_3)),sep="|")
    if (is.null(passed_doublets[[gene_pair_2]])) {
      next
    }
    
    gene_pair_3 = glue::collapse(sort(c(gene_1,gene_3)),sep="|")
    if (is.null(passed_doublets[[gene_pair_3]])) {
      next
    }
    
    overlap_count = length(intersect(passed_doublets[[gene_pair_1]],passed_doublets[[gene_pair_2]]))
    
    if (overlap_count < min_shared_cells) {
      next;
    }
    
    overlap_expect = (percent_cells[[gene_1]]*
                        percent_cells[[gene_2]]*
                        percent_cells[[gene_3]])*total_cells
    
    if (overlap_count > overlap_expect) {
      passed_sets$gene_1 = c(passed_sets$gene_1,gene_1)
      passed_sets$gene_2 = c(passed_sets$gene_2,gene_2)
      passed_sets$gene_3 = c(passed_sets$gene_3,gene_3)
      passed_sets$overlap_expect = c(passed_sets$overlap_expect,overlap_expect)
      passed_sets$overlap_observe = c(passed_sets$overlap_observe,overlap_count)
    }
  }
  
  passed_sets = as.data.frame(passed_sets)
  
  passed_sets = passed_sets %>% mutate(overlap_diff = overlap_observe - overlap_expect) %>%
    arrange(desc(overlap_diff))
  
  passed_sets
}

gather_gene_sets_parallel = function(tidy_10X, min_shared_cells = 100) {
  total_cells = length(unique(tidy_10X$barcode))
  
  library(doParallel)
  
  cl <- makeCluster(10)
  
  registerDoParallel(cl)
  
  min_shared_cells = 100
  
  gene_list = unique(bone_marrow_kinases$symbol)
  barcode_sets = foreach(gene = gene_list, 
                         .packages = 'tidyverse') %dopar% {
                           temp = bone_marrow_kinases %>%
                             filter(symbol == gene) %>%
                             select(barcode)
                           temp
                         }
  barcode_sets = setNames(barcode_sets,gene_list)
  
  for (gene in names(barcode_sets)) {
    if (length(barcode_sets[[gene]]) < min_shared_cells) {
      barcode_sets[[gene]] <- NULL
    }
  }
  
  gene_list = names(barcode_sets)
  percent_cells = foreach (gene = gene_list) %dopar% {
    length(barcode_sets[[gene]]$barcode)/total_cells
  }
  percent_cells = setNames(percent_cells,gene_list)
  
  ##Doublets
  kinase_combinations = combn(names(barcode_sets),2)
  
  passed_doublets = list()
  
  library(progress)
  
  pb <- progress_bar$new(total = dim(kinase_combinations)[2])
  
  passed_doublets = foreach (index = 1:dim(kinase_combinations)[2]) %dopar% {
    library(tidyverse)
    pb$tick()
    gene_1 = as.character(kinase_combinations[1,index])
    gene_2 = as.character(kinase_combinations[2,index])
    
    overlap_set = intersect(barcode_sets[[gene_1]]$barcode,barcode_sets[[gene_2]]$barcode)
    
    gene_str = glue::collapse(sort(c(gene_1,gene_2)),sep="|")
    passed_doublets[[gene_str]] = overlap_set
  }
  
  ##Triplets
  kinase_combinations = combn(names(barcode_sets),3)
  
  passed_sets = list()
  
  
  pb <- progress_bar$new(total = dim(kinase_combinations)[2])
  count = 0;
  foreach (index = 1:dim(kinase_combinations)[2]) %dopar% {
    library(tidyverse)
    pb$tick()
    
    gene_1 = as.character(kinase_combinations[1,index])
    gene_2 = as.character(kinase_combinations[2,index])
    
    gene_pair_1 = glue::collapse(sort(c(gene_1,gene_2)),sep="|")
    if (is.null(passed_doublets[[gene_pair_1]])) {
      next
    }
    
    gene_3 = as.character(kinase_combinations[3,index])
    
    gene_pair_2 = glue::collapse(sort(c(gene_2,gene_3)),sep="|")
    if (is.null(passed_doublets[[gene_pair_2]])) {
      next
    }
    
    gene_pair_3 = glue::collapse(sort(c(gene_1,gene_3)),sep="|")
    if (is.null(passed_doublets[[gene_pair_3]])) {
      next
    }
    
    overlap_count = length(intersect(passed_doublets[[gene_pair_1]],passed_doublets[[gene_pair_2]]))
    
    if (overlap_count < min_shared_cells) {
      next;
    }
    
    overlap_expect = (percent_cells[[gene_1]]*
                        percent_cells[[gene_2]]*
                        percent_cells[[gene_3]])*total_cells
    
    if (overlap_count > overlap_expect) {
      passed_sets$gene_1 = c(passed_sets$gene_1,gene_1)
      passed_sets$gene_2 = c(passed_sets$gene_2,gene_2)
      passed_sets$gene_3 = c(passed_sets$gene_3,gene_3)
      passed_sets$gene_names = glue::collapse(sort(c(gene_1,gene_2,gene_3)),sep="|")
      passed_sets$overlap
      passed_sets$overlap_expect = c(passed_sets$overlap_expect,overlap_expect)
      passed_sets$overlap_observe = c(passed_sets$overlap_observe,overlap_count)
    }
  }
  
  passed_sets = as.data.frame(passed_sets)
  
  passed_sets = passed_sets %>% mutate(overlap_diff = overlap_observe - overlap_expect) %>%
    arrange(desc(overlap_diff))
  
  passed_sets
}

gather_gene_sets_general = function(tidy_10X, min_shared_cells = 100, min_percent_cells = NA,
                                    max_cluster_size = 5) {
  
  total_cells = length(unique(tidy_10X$barcode))
  print(paste0('Found ', total_cells, ' cells in the data set.'))
  
  if (! is.na(min_percent_cells)) {
    min_shared_cells = round(total_cells*min_percent_cells);
  }
  
  barcode_sets = list()
  barcode_sets[[1]] = list()
  
  percent_cells = list()
  library(progress)
  print('Working on gathering cell ID barcodes for each gene.')
  pb <- progress_bar$new(total = length(unique(tidy_10X$symbol)))
  for (gene in unique(tidy_10X$symbol)) {
    pb$tick()
    temp = tidy_10X %>%
      filter(symbol == gene) %>%
      select(barcode)
    
    if (dim(temp)[1] >= min_shared_cells) {
      barcode_sets[[1]][[gene]] = temp$barcode
      percent_cells[[gene]] = length(barcode_sets[[1]][[gene]])/total_cells
    }
  }
  print(paste0('Found ',length(barcode_sets[[1]]), ' genes present in at least ',
               min_shared_cells, ' cells.'))

  overrepresented = list()
  cluster_num = 1
  
  while (length(names(barcode_sets[[cluster_num]])) > 0 & cluster_num < max_cluster_size) {
    cluster_num = cluster_num + 1
    
    prev_cluster_size = cluster_num - 1
    
    previous_kinase_combinations = names(barcode_sets[[prev_cluster_size]])
    
    barcode_sets[[cluster_num]] = list()
    overrepresented[[cluster_num]] = list()
    
    print(paste0('Working on size ', cluster_num, " clusters."))
    pb <- progress_bar$new(total = length(previous_kinase_combinations))
    
    tested = list()
    for (this_combination in previous_kinase_combinations) {
      pb$tick()
      
      starting_gene_set = strsplit(this_combination,'\\|')[[1]]
      
      search_set = setdiff(names(barcode_sets[[1]]),starting_gene_set)
      
      for (this_gene in search_set) {
        gene_set = c(starting_gene_set,this_gene)
        full_gene_set = glue::collapse(sort(c(gene_set)),sep="|");
        if (! is.null(tested[[full_gene_set]])) {
          next;
        }
        
        overlap_set = intersect(barcode_sets[[prev_cluster_size]][[this_combination]],
                                barcode_sets[[1]][[this_gene]])
        
        tested[[full_gene_set]] = 1
        if (length(overlap_set) > min_shared_cells) {
          barcode_sets[[cluster_num]][[full_gene_set]] = overlap_set;
        } else {
          next;
        }

        expected_overlap = percent_cells[[gene_set[1]]]*total_cells;
        for (gene_name in gene_set[2:length(gene_set)]) {
          expected_overlap = expected_overlap*percent_cells[[gene_name]]
        }
        
        if (length(overlap_set) > expected_overlap) {
          overrepresented[[cluster_num]]$gene_set = c(overrepresented[[cluster_num]]$gene_set,
                                                     glue::collapse(sort(gene_set),sep="|"))
          overrepresented[[cluster_num]]$overlap_expect = c(overrepresented[[cluster_num]]$overlap_expect,
                                                           expected_overlap)
          overrepresented[[cluster_num]]$overlap_observe = c(overrepresented[[cluster_num]]$overlap_observe,
                                                            length(overlap_set))
        }
      }
    }
    # overrepresented[[cluster_num]] = as.data.frame(over_represented[[cluster_num]])
    print(paste0('Found ', length(barcode_sets[[cluster_num]]), 
                 ' sets, of which ', length(overrepresented[[cluster_num]]$gene_set), 
                 ' are over-represented.'))
  }
  overrepresented
}

