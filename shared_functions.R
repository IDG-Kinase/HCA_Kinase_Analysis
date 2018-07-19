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