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

filter_single_cell = function(tidy_10X, read_total_range = c(-Inf,Inf),
                              unique_genes_range = c(-Inf,Inf), 
                              mt_range = c(-Inf, Inf),
                              mt_mad_range = c(-Inf, Inf)) {
  
  passed_cells = tidy_10X %>% group_by(barcode) %>%
    summarise(read_counts = sum(counts),
              unique_genes = n(),
              mt_perc = sum(counts[grep('MT', gene_names)])/sum(counts)*100) %>%
    filter(between(read_counts, read_total_range[1], read_total_range[2]),
           between(unique_genes, unique_genes_range[1], unique_genes_range[2]),
           between(mt_perc, mt_range[1], mt_range[2]),
           between(mt_perc, median(mt_perc) + mad(mt_perc)*mt_mad_range[1], 
                   median(mt_perc) + mad(mt_perc)*mt_mad_range[2]))
  
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

###############################################################################
# Gene Set Clustering
###############################################################################
gather_gene_sets = function(tidy_10X, 
                            min_shared_cells = 100, 
                            min_percent_cells = NA,
                            max_cluster_size = Inf) {
  
  total_cells = length(unique(tidy_10X$barcode))
  print(paste0('Found ', total_cells, ' cells in the data set.'))
  
  if (! is.na(min_percent_cells)) {
    min_shared_cells = round(total_cells*min_percent_cells);
  }
  
  tidy_10X = tidy_10X %>% left_join(
    data.frame(barcode = unique(tidy_10X$barcode)) %>% mutate(barcode_num = 1:n())
  )
  
  barcode_sets = list()
  barcode_sets[[1]] = list()
  
  percent_cells = list()
  library(progress)
  print('Working on gathering cell ID barcodes for each gene.')
  pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                         total = length(unique(tidy_10X$symbol)))
  for (gene in unique(tidy_10X$symbol)) {
    pb$tick()
    temp = tidy_10X %>%
      filter(symbol == gene) %>%
      select(barcode_num)
    
    if (dim(temp)[1] >= min_shared_cells) {
      barcode_sets[[1]][[gene]] = temp$barcode_num
      percent_cells[[gene]] = length(barcode_sets[[1]][[gene]])/total_cells
    }
  }
  print(paste0('Found ',length(barcode_sets[[1]]), ' genes present in at least ',
               min_shared_cells, ' cells.'))
  
  cluster_props = list()
  cluster_num = 1
  
  while (length(names(barcode_sets[[cluster_num]])) > 0 & cluster_num < max_cluster_size) {
    cluster_num = cluster_num + 1
    
    prev_cluster_size = cluster_num - 1
    
    previous_kinase_combinations = names(barcode_sets[[prev_cluster_size]])
    
    barcode_sets[[cluster_num]] = list()
    cluster_props[[cluster_num]] = list()
    
    print(paste0('Working on size ', cluster_num, " clusters."))
    pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                           total = length(previous_kinase_combinations))
    
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
        
        cluster_props[[cluster_num]]$gene_set = c(cluster_props[[cluster_num]]$gene_set,
                                                  glue::collapse(sort(gene_set),sep="|"))
        cluster_props[[cluster_num]]$overlap_expect = c(cluster_props[[cluster_num]]$overlap_expect,
                                                        expected_overlap)
        cluster_props[[cluster_num]]$overlap_observe = c(cluster_props[[cluster_num]]$overlap_observe,
                                                         length(overlap_set))
      }
    }
    
    if (length(cluster_props[[cluster_num]]) != 0) {
      cluster_props[[cluster_num]] = as.data.frame(cluster_props[[cluster_num]]) %>%
        mutate(overlap_diff = overlap_observe - overlap_expect) %>%
        mutate(overlap_diff_percent = overlap_diff/total_cells)
    }
    
    print(paste0('Found ', length(barcode_sets[[cluster_num]]), 
                 ' sets, tested ', length(tested), ' sets.'))
  }
  cluster_props
}

overlap_gene_sets <- function(gene_set,barcode_sets,cluster_num,min_shared_cells) {
  overlap_set = list()
  for (this_index in 1:dim(gene_set)[1]) {
    previous_set = as.character(gene_set[this_index,]$previous_set);
    added_gene = as.character(gene_set[this_index,]$added_gene);
    barcode_overlap = intersect(barcode_sets[[cluster_num]][[previous_set]],
                                barcode_sets[[1]][[added_gene]])
    
    if (length(barcode_overlap) > min_shared_cells) {
      split_previous_set = strsplit(previous_set,'\\|')[[1]]
      full_gene_set = glue::collapse(sort(c(previous_set,added_gene)),sep="|");
      overlap_set[[full_gene_set]] = barcode_overlap;
    }
  }
  overlap_set
}


gather_gene_set_parallel <- function(tidy_10X, 
                                     min_shared_cells = 100, 
                                     min_percent_cells = NA,
                                     max_cluster_size = Inf,
                                     num_cores = NA,
                                     parallel_process_set_size = 100000) {
  
  total_cells = length(unique(tidy_10X$barcode))
  print(paste0('Found ', total_cells, ' cells in the data set.'))
  
  if (! is.na(min_percent_cells)) {
    min_shared_cells = round(total_cells*min_percent_cells);
  }
  
  if (is.na(num_cores)) {
    num_cores = parallel::detectCores()-1
  }
  
  tidy_10X = tidy_10X %>% left_join(
    data.frame(barcode = unique(tidy_10X$barcode)) %>% mutate(barcode_num = 1:n())
  )
  
  barcode_sets = list()
  barcode_sets[[1]] = list()
  
  percent_cells = list()
  library(progress)
  print('Working on gathering cell ID barcodes for each gene.')
  pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                         total = length(unique(tidy_10X$symbol)))
  for (gene in unique(tidy_10X$symbol)) {
    pb$tick()
    temp = tidy_10X %>%
      filter(symbol == gene) %>%
      select(barcode_num)
    
    if (dim(temp)[1] >= min_shared_cells) {
      barcode_sets[[1]][[gene]] = temp$barcode_num
      percent_cells[[gene]] = length(barcode_sets[[1]][[gene]])/total_cells
    }
  }
  print(paste0('Found ',length(barcode_sets[[1]]), ' genes present in at least ',
               min_shared_cells, ' cells.'))
  
  print("Working on building possible 2 clusters.")  
  to_test = as.data.frame(t(combn(names(barcode_sets[[1]]),2)))
  names(to_test) <- c('previous_set','added_gene')
  to_test$previous_set = as.character(to_test$previous_set)
  to_test$added_gene = as.character(to_test$added_gene)
  
  split_doublets = split(to_test,
                         rep(1:(parallel::detectCores()-1),
                             length.out = dim(to_test)[1]))
  
  overlap_gene_sets <- function(gene_set,barcode_sets,cluster_num,min_shared_cells) {
    overlap_set = list()
    for (this_index in 1:dim(gene_set)[1]) {
      previous_set = as.character(gene_set[this_index,]$previous_set);
      added_gene = as.character(gene_set[this_index,]$added_gene);
      barcode_overlap = intersect(barcode_sets[[cluster_num]][[previous_set]],
                                  barcode_sets[[1]][[added_gene]])
      
      if (length(barcode_overlap) > min_shared_cells) {
        split_previous_set = strsplit(previous_set,'\\|')[[1]]
        full_gene_set = glue::collapse(sort(c(previous_set,added_gene)),sep="|");
        overlap_set[[full_gene_set]] = barcode_overlap;
      }
    }
    overlap_set
  }
  
  overlap_sets_split = mclapply(split_doublets,
                                function(x) { overlap_gene_sets(x,barcode_sets,1,min_shared_cells); },
                                mc.cores=num_cores);
  
  barcode_sets[[2]] = list()
  for (this_overlap_set in overlap_sets_split) {
    for (this_combo in names(this_overlap_set)) {
      barcode_sets[[2]][[this_combo]] = this_overlap_set[[this_combo]]
    }
  }
  print(paste0('Found ', length(barcode_sets[[2]]), " sets, tested ", dim(to_test)[1], " sets."))
  
  cluster_props = list()
  cluster_num = 2
  
  while (length(names(barcode_sets[[cluster_num]])) > 0 & cluster_num < max_cluster_size) {
    cluster_num = cluster_num + 1
    
    prev_cluster_size = cluster_num - 1
    
    previous_kinase_combinations = sort(names(barcode_sets[[prev_cluster_size]]))
    
    barcode_sets[[cluster_num]] = list()
    cluster_props[[cluster_num]] = list()
    
    print(paste0('Working on building possible ', cluster_num, " clusters."))
    pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                           total = length(previous_kinase_combinations))
    
    to_test = data.frame(previous_set=as.character(),added_gene=as.character())
    tested_count = 0;
    for (this_combination in previous_kinase_combinations) {
      pb$tick()
      
      starting_gene_set = strsplit(this_combination,'\\|')[[1]]
      
      search_set = glue::collapse(starting_gene_set[2:length(starting_gene_set)],sep="\\|")
      
      matching_clusters = grep(paste0('^',search_set),previous_kinase_combinations, value=TRUE)
      last_gene = sapply(matching_clusters,function(x) { tail(strsplit(x,'\\|')[[1]],1) })
      if (length(last_gene) == 0) {
        next;
      }
      to_test = to_test %>% add_row(previous_set = this_combination,added_gene = last_gene)
      
      if (dim(to_test)[1] > parallel_process_set_size) {
        split_sets = split(to_test,
                           rep(1:(parallel::detectCores()-1),
                               length.out = dim(to_test)[1]))
        
        overlap_sets_split = mclapply(split_sets,
                                      function(x) { overlap_gene_sets(x,barcode_sets,prev_cluster_size,min_shared_cells); },
                                      mc.cores=num_cores);
        
        for (this_overlap_set in overlap_sets_split) {
          for (this_combo in names(this_overlap_set)) {
            barcode_sets[[cluster_num]][[this_combo]] = this_overlap_set[[this_combo]]
          }
        }
        tested_count = tested_count + dim(to_test)[1]
        to_test = data.frame(previous_set=as.character(),added_gene=as.character())
      }
    }
    
    split_sets = split(to_test,
                       rep(1:(parallel::detectCores()-1),
                           length.out = dim(to_test)[1]))
    
    overlap_sets_split = mclapply(split_sets,
                                  function(x) { overlap_gene_sets(x,barcode_sets,prev_cluster_size,min_shared_cells); },
                                  mc.cores=num_cores);
    
    for (this_overlap_set in overlap_sets_split) {
      for (this_combo in names(this_overlap_set)) {
        barcode_sets[[cluster_num]][[this_combo]] = this_overlap_set[[this_combo]]
      }
    }
    tested_count = tested_count + dim(to_test)[1]
    
    print(paste0('Found ', length(barcode_sets[[cluster_num]]), " sets, tested ", tested_count, " sets."))
    
  }
  
  cluster_props = list();
  for (cluster_num in 2:(length(barcode_sets)-1)) {
    cluster_props[[cluster_num]] = list()
    for (combo_set in names(barcode_sets[[cluster_num]])) {
      gene_set = strsplit(combo_set,'\\|')[[1]]
      
      expected_overlap = percent_cells[[gene_set[1]]]*total_cells;
      for (gene_name in gene_set[2:length(gene_set)]) {
        expected_overlap = expected_overlap*percent_cells[[gene_name]]
      }
      
      cluster_props[[cluster_num]]$gene_set = c(cluster_props[[cluster_num]]$gene_set,
                                                glue::collapse(sort(gene_set),sep="|"))
      cluster_props[[cluster_num]]$overlap_expect = c(cluster_props[[cluster_num]]$overlap_expect,
                                                      expected_overlap)
      cluster_props[[cluster_num]]$overlap_observe = c(cluster_props[[cluster_num]]$overlap_observe,
                                                       length(barcode_sets[[cluster_num]][[combo_set]]))
      
    }
    
    if (length(cluster_props[[cluster_num]]) != 0) {
      cluster_props[[cluster_num]] = as.data.frame(cluster_props[[cluster_num]]) %>%
        mutate(overlap_diff = overlap_observe - overlap_expect) %>%
        mutate(overlap_diff_percent = overlap_diff/total_cells)
    }
  }
  cluster_props
}