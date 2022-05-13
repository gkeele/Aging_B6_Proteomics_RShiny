
########### Functions

library(tidyverse)
library(ggbeeswarm)
library(ggrepel)
#library(ensimplR)
library(plotly)

protein_all_tissues_plotting_dat <- read.csv("./data/protein_all_tissues_plotting.csv", header = TRUE) %>%
  mutate(mouse.id = factor(mouse.id),Sex = factor(Sex, levels = c("M", "F")), Age = factor(Age, levels = c("Y", "O")))

protein_all_tissues_plotting_dat$Age = plyr::revalue(protein_all_tissues_plotting_dat$Age, c("O" = "Old","Y" = "Young"))
protein_all_tissues_plotting_dat$Sex = plyr::revalue(protein_all_tissues_plotting_dat$Sex, c("M" = "Male","F" = "Female"))
  
## Read in tissue-specific results and Convert to long data frame 
tissue_specific_results <- readRDS("./data/tissue_specific_results.rds")
tissue_specific_results_long <- plyr::ldply(tissue_specific_results, data.frame) %>%
  dplyr::select(-.id) %>% group_by(tissue) %>%
  mutate(sex_qval = p.adjust(sex_pval, method = "BH"),
         age_qval = p.adjust(age_pval, method = "BH"))

## GO gene sets		 
gene_sets_dat <- readRDS("./data/all_gene_sets.RDS")

### New dataset loading
fgsea_age_tables <- readRDS("./data/fgsea_age_tables.rds")
fgsea_sex_tables <- readRDS("./data/fgsea_sex_tables.rds")

## Effects
#tissue_specific_results <- readRDS("./data/GO/tissue_specific_results.rds")

## Grab protein.id, gene.id, and symbol information
protein_to_symbol <- do.call("rbind", tissue_specific_results) %>%
  dplyr::select(protein.id, gene.id, symbol) %>%
  distinct

  

		 
### 1. Sex & Age effect function for Box plots

### Sex by Age only
Sex_Age_Effect_1 <- function(pdata = NULL) {
  
  # Plot sex_age effects in all tissues
  S_A_plot = ggplot(data = pdata, aes(y = Intensity, x = Sex)) +
    geom_boxplot(aes(fill=Age),alpha = 0,outlier.shape = NA) +			### aes(fill=Age),  alpha = 0,
    #geom_jitter(col = "gray") +
    geom_point(aes(color=Age),position=position_jitterdodge())+		### aes(group=Age),
    facet_grid(toupper(symbol)~tissue,scales = "free_x")+
	#facet_wrap(toupper(symbol)	~tissue,ncol=5)+			### toupper(symbol)
	scale_x_discrete(drop = FALSE)+
	#theme(text=element_text(size=21))+
	#theme(strip.text.x = element_text(size = 30)) +
    theme_bw()
	
  #ggplotly(S_A_plot) 
  #S_A_plot = ggplotly(S_A_plot) %>% layout(boxmode='group')	
  return(list(S_A_plot = S_A_plot))
     
}

Sex_Age_Effect_2 <- function(pdata = NULL, type="SA") {
  
  # Plot sex_age effects in all tissues
  if(type=="SA"){
	  S_A_plot = ggplot(data = pdata, aes(y = Intensity, x = Sex)) +
		geom_boxplot(aes(fill=Age),alpha = 0,outlier.shape = NA) +			### aes(fill=Age),  alpha = 0,
		#geom_jitter(col = "gray") +
		geom_point(aes(color=Age),position=position_jitterdodge())+		### aes(group=Age),
		facet_grid(toupper(symbol)~tissue,scales = "free_x")+
		#facet_wrap(toupper(symbol)	~tissue,ncol=5)+			### toupper(symbol)
		scale_x_discrete(drop = FALSE)+
		#theme(text=element_text(size=21))+
		#theme(strip.text.x = element_text(size = 30)) +
		theme_bw()
	}
	
	if(type =="sex"){
		S_A_plot = ggplot(data = pdata, aes(y = Intensity, x = Sex)) +
		geom_boxplot(alpha = 0,outlier.shape = NA) +			### aes(fill=Age),  alpha = 0,
		#geom_jitter(col = "gray") +
		geom_point(aes(color=Sex),position=position_jitterdodge())+		### aes(group=Age),
		facet_grid(toupper(symbol)~tissue,scales = "free_x")+
		scale_x_discrete(drop = FALSE)+
		theme_bw()		
	}
	
	if(type =="age"){
		S_A_plot = ggplot(data = pdata, aes(y = Intensity, x = Age)) +
		geom_boxplot(alpha = 0,outlier.shape = NA) +			### aes(fill=Age),  alpha = 0,
		#geom_jitter(col = "gray") +
		geom_point(aes(color=Age),position=position_jitterdodge())+		### aes(group=Age),
		facet_grid(toupper(symbol)~tissue,scales = "free_x")+
		scale_x_discrete(drop = FALSE)+
		theme_bw()		
	}

	
  return(list(S_A_plot = S_A_plot))
     
}

####### Function to capitalize first letter of word
single_simple_cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
simple_cap <- function(x) {
  sapply(1:length(x), function(i) single_simple_cap(x[i]))
}



### 3. Sex & Age effect comparison plots

tissue_fgsea_compare_pathwaylevel <- function(tissue1,
                                              tissue2,
                                              fgsea_list) {
  
  tissue1_fgsea_dat <- fgsea_list[[tissue1]]
  tissue2_fgsea_dat <- fgsea_list[[tissue2]]
  
  ## Merge data
  fgsea_dat <- inner_join(tissue1_fgsea_dat %>%
                            dplyr::select(pathway, NES, padj) %>%
                            dplyr::rename(tissue1_effect = NES,
                                          tissue1_padj = padj),
                          tissue2_fgsea_dat %>%
                            dplyr::select(pathway, NES, padj) %>%
                            dplyr::rename(tissue2_effect = NES,
                                          tissue2_padj = padj)) %>%
    mutate(consistent = sign(tissue1_effect) == sign(tissue2_effect)) %>%
    dplyr::rowwise() %>%
    mutate(mean_padj = mean(tissue1_padj, tissue2_padj)) %>%
    dplyr::ungroup() %>%
    arrange(mean_padj) %>%
    dplyr::select(-mean_padj)
  
  fgsea_dat
}

## Plot to compare effects
tissue_fgsea_compare_genelevel_plot <- function(tissue1,
                                                tissue2,
                                                effects_list,
                                                fgsea_list,
                                                fdr = 0.1,
                                                include_thresh = NULL,
                                                gene_sets_dat,
                                                effect_type = c("age", "sex"),
                                                highlight_pathway,
                                                highlight_col = "dodgerblue3") {
  
  include_thresh <- ifelse(is.null(include_thresh), fdr, include_thresh)
  tissue1_fgsea_dat <- fgsea_list[[tissue1]]
  tissue1_effects_dat <- effects_list[[tissue1]]
  tissue2_fgsea_dat <- fgsea_list[[tissue2]]
  tissue2_effects_dat <- effects_list[[tissue2]]
  effect_type <- effect_type[1]
  
  ## Expand fgsea data to gene-level
  highlight_genes <- gene_sets_dat %>%
    filter(gs_name == highlight_pathway) %>% 
    pull(ensembl_gene) %>%
    unique
  
  genes_fgsea_dat <- inner_join(tissue1_effects_dat %>%
                                  dplyr::select(gene.id, protein.id, symbol, paste(effect_type, c("effect", "effect_se", "pval"), sep = "_")) %>%
                                  dplyr::rename(tissue1_effect = paste(effect_type, "effect", sep = "_"),
                                                tissue1_effect_se = paste(effect_type, "effect_se", sep = "_"),
                                                tissue1_pval = paste(effect_type, "pval", sep = "_")) %>%
                                  mutate(tissue1_padj = p.adjust(tissue1_pval, method = "BH")) %>%
                                  dplyr::select(-tissue1_pval), 
                                tissue2_effects_dat %>%
                                  dplyr::select(gene.id, protein.id, symbol, paste(effect_type, c("effect", "effect_se", "pval"), sep = "_")) %>%
                                  dplyr::rename(tissue2_effect = paste(effect_type, "effect", sep = "_"),
                                                tissue2_effect_se = paste(effect_type, "effect_se", sep = "_"),
                                                tissue2_pval = paste(effect_type, "pval", sep = "_")) %>%
                                  mutate(tissue2_padj = p.adjust(tissue2_pval, method = "BH")) %>%
                                  dplyr::select(-tissue2_pval)) %>%
    mutate(fdr_status = case_when(tissue1_padj > fdr & tissue2_padj > fdr ~ "Neither tissue",
                                  (tissue1_padj < fdr & tissue2_padj > fdr) | (tissue1_padj > fdr & tissue2_padj < fdr) ~ "One tissue",
                                  tissue1_padj < fdr & tissue2_padj < fdr ~ "Both tissues")) %>%
    mutate(fdr_status = factor(fdr_status, levels = c("Neither tissue", "One tissue", "Both tissues"))) %>%
    mutate(highlight = gene.id %in% highlight_genes)

  p <- ggplot(data = genes_fgsea_dat %>%
                arrange(highlight) %>%
                filter(tissue1_padj < include_thresh | tissue2_padj < include_thresh),
              aes(y = tissue1_effect/tissue1_effect_se,
                  x = tissue2_effect/tissue2_effect_se,
                  size = fdr_status,
                  col = highlight)) +
    geom_point() +
    scale_color_manual(values = c("gray", highlight_col), drop = FALSE) +
    scale_size_manual(values = c(1, 2, 3), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    ylab(paste("Standardized", effect_type, "effect in", tissue1)) + xlab(paste("Standardized", effect_type, "effect in", tissue2)) +
    ggtitle(paste("Comparison of", effect_type, "effects on proteins between", tissue1, "and", tissue2),
            subtitle = paste("GO category:", highlight_pathway)) +
    theme_bw() +
    theme(legend.position = c(0.2, 0.85),
          legend.text = element_text(size = 8)) +
    guides(size = guide_legend(title = paste("FDR <", fdr)),
           col = guide_legend(title = "GO category protein"))
  
  p
}

## Return gene-level table
tissue_fgsea_compare_genelevel_table <- function(tissue1,
                                                 tissue2,
                                                 effects_list,
                                                 annot_dat,
                                                 effect_type = c("age", "sex")) {
  
  tissue1_effects_dat <- effects_list[[tissue1]]
  tissue2_effects_dat <- effects_list[[tissue2]]
  effect_type <- effect_type[1]
  
  ## Grab gene-level intersection
  genes_dat <- inner_join(tissue1_effects_dat %>%
                            dplyr::select(gene.id, protein.id, paste(effect_type, c("effect", "effect_se", "pval"), sep = "_")) %>%
                            dplyr::rename(tissue1_effect = paste(effect_type, "effect", sep = "_"),
                                          tissue1_effect_se = paste(effect_type, "effect_se", sep = "_"),
                                          tissue1_pval = paste(effect_type, "pval", sep = "_")) %>%
                            mutate(tissue1_padj = p.adjust(tissue1_pval, method = "BH")) %>%
                            dplyr::select(-tissue1_pval), 
                          tissue2_effects_dat %>%
                            dplyr::select(gene.id, protein.id, paste(effect_type, c("effect", "effect_se", "pval"), sep = "_")) %>%
                            dplyr::rename(tissue2_effect = paste(effect_type, "effect", sep = "_"),
                                          tissue2_effect_se = paste(effect_type, "effect_se", sep = "_"),
                                          tissue2_pval = paste(effect_type, "pval", sep = "_")) %>%
                            mutate(tissue2_padj = p.adjust(tissue2_pval, method = "BH")) %>%
                            dplyr::select(-tissue2_pval)) %>%
    left_join(annot_dat) %>%
    arrange(tissue1_padj)
  
  genes_dat
}

## Generate data to compare effects
tissue_fgsea_data <- function(tissue1,
                              tissue2,
                              effects_list,
                              fgsea_list,
                              fdr = 0.1,
							  gene_sets_dat,
                              effect_type = c("age", "sex"),
                              highlight_pathway,
                              highlight_col = "dodgerblue3",
                              include_nonsig = FALSE) {
  
  tissue1_fgsea_dat <- fgsea_list[[tissue1]]
  tissue1_effects_dat <- effects_list[[tissue1]]
  tissue2_fgsea_dat <- fgsea_list[[tissue2]]
  tissue2_effects_dat <- effects_list[[tissue2]]
  effect_type <- effect_type[1]
  
  ## Expand fgsea data to gene-level
  highlight_genes <- gene_sets_dat %>%
    filter(gs_name == highlight_pathway) %>% 
    pull(ensembl_gene) %>%
    unique
   
  genes_fgsea_dat <- inner_join(tissue1_effects_dat %>%
                                  dplyr::select(gene.id, protein.id, symbol, paste(effect_type, c("effect", "effect_se", "pval"), sep = "_")) %>%
                                  dplyr::rename(tissue1_effect = paste(effect_type, "effect", sep = "_"),
                                                tissue1_effect_se = paste(effect_type, "effect_se", sep = "_"),
                                                tissue1_pval = paste(effect_type, "pval", sep = "_")) %>%
                                  mutate(tissue1_padj = p.adjust(tissue1_pval, method = "BH")) %>%
                                  dplyr::select(-tissue1_pval), 
                                tissue2_effects_dat %>%
                                  dplyr::select(gene.id, protein.id, symbol, paste(effect_type, c("effect", "effect_se", "pval"), sep = "_")) %>%
                                  dplyr::rename(tissue2_effect = paste(effect_type, "effect", sep = "_"),
                                                tissue2_effect_se = paste(effect_type, "effect_se", sep = "_"),
                                                tissue2_pval = paste(effect_type, "pval", sep = "_")) %>%
                                  mutate(tissue2_padj = p.adjust(tissue2_pval, method = "BH")) %>%
                                  dplyr::select(-tissue2_pval)) %>%
    mutate(fdr_status = case_when(tissue1_padj > fdr & tissue2_padj > fdr ~ "Neither tissue",
                                  (tissue1_padj < fdr & tissue2_padj > fdr) | (tissue1_padj > fdr & tissue2_padj < fdr) ~ "One tissue",
                                  tissue1_padj < fdr & tissue2_padj < fdr ~ "Both tissues")) %>%
    mutate(fdr_status = factor(fdr_status, levels = c("Neither tissue", "One tissue", "Both tissues"))) %>%
    mutate(highlight_T = gene.id %in% highlight_genes) %>% 
	mutate(highlight = ifelse(highlight_T == TRUE, "GO category protein", "Not GO category protein"))
	
	t_level = levels(genes_fgsea_dat$fdr_status)
	genes_fgsea_dat$fdr_status = paste0(paste0('FDR<', fdr, ' in '),genes_fgsea_dat$fdr_status)
	genes_fgsea_dat$fdr_status = factor(genes_fgsea_dat$fdr_status, levels = paste0(paste0('FDR<', fdr, ' in '),t_level))
    genes_fgsea_dat$highlight = factor(genes_fgsea_dat$highlight, levels = c("Not GO category protein","GO category protein"))
      
  return(genes_fgsea_dat)
}


