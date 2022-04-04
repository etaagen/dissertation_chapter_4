##### ASR_input_CV

# get ASR formatted input variables from raw full pop data, using CV matrix 

ASR_input_CV <- function(founder_data, QTL_n){
  ### STEP 0 
  # info column tracks repeating rows (multiple transcripts)
  founder_data$info <- as.factor(founder_data$info)
  # get loci with HIGH SnpEff category
  # if there are multiple transcripts, select 1 randomly (don't want multiple variants of same gene)
  # X maintains the genetic map order 
  high <- filter(founder_data, grepl("HIGH", SnpEff)) %>% group_by(info) %>% sample_n(1) %>% arrange(X)
  non_synon_coding <- filter(founder_data, grepl("NON_SYNONYMOUS_CODING", SnpEff)) %>% group_by(info) %>% sample_n(1) %>% arrange(X)
  # take a random % (not all DV are causal / impact simulated phenotype)
  # we expect the vast majority of `HIGH` category to be deleterious, take 95%
  # group by chromosome to maintain ratios 
  # will later serve as QTL base in _DV approach
  high_sub <- high %>% group_by(Chr) %>%  slice_sample(prop = 0.90) %>% arrange(X)
  # and subset to 25% of non_synon_coding (2827) within the 20 cM increased recombination target
  chr1A <- non_synon_coding %>% filter(Chr == "1A", Pos.cM. > 0.35 & Pos.cM. < 0.55)
  chr1B <- non_synon_coding %>% filter(Chr == "1B", Pos.cM. > 0.42 & Pos.cM. < 0.62)
  chr1D <- non_synon_coding %>% filter(Chr == "1D", Pos.cM. > 0.35 & Pos.cM. < 0.55)
  chr2A <- non_synon_coding %>% filter(Chr == "2A", Pos.cM. > 0.57 & Pos.cM. < 0.77)
  chr2B <- non_synon_coding %>% filter(Chr == "2B", Pos.cM. > 0.53 & Pos.cM. < 0.73)
  chr2D <- non_synon_coding %>% filter(Chr == "2D", Pos.cM. > 0.55 & Pos.cM. < 0.75)
  chr3A <- non_synon_coding %>% filter(Chr == "3A", Pos.cM. > 0.39 & Pos.cM. < 0.59)
  chr3B <- non_synon_coding %>% filter(Chr == "3B", Pos.cM. > 0.42 & Pos.cM. < 0.62)
  chr3D <- non_synon_coding %>% filter(Chr == "3D", Pos.cM. > 0.35 & Pos.cM. < 0.55)
  chr4A <- non_synon_coding %>% filter(Chr == "4A", Pos.cM. > 1 & Pos.cM. < 1.2)
  chr4B <- non_synon_coding %>% filter(Chr == "4B", Pos.cM. > 0.43 & Pos.cM. < 0.63)
  chr4D <- non_synon_coding %>% filter(Chr == "4D") # not enough markers to filter
  chr5A <- non_synon_coding %>% filter(Chr == "5A", Pos.cM. > 0.44 & Pos.cM. < 0.64)
  chr5B <- non_synon_coding %>% filter(Chr == "5B", Pos.cM. > 0.33 & Pos.cM. < 0.53)
  chr5D <- non_synon_coding %>% filter(Chr == "5D", Pos.cM. > 0.48 & Pos.cM. < 0.68)
  chr6A <- non_synon_coding %>% filter(Chr == "6A", Pos.cM. > 0.42 & Pos.cM. < 0.62)
  chr6B <- non_synon_coding %>% filter(Chr == "6B", Pos.cM. > 0.39 & Pos.cM. < 0.59)
  chr6D <- non_synon_coding %>% filter(Chr == "6D", Pos.cM. > 0.4 & Pos.cM. < 0.6)
  chr7A <- non_synon_coding %>% filter(Chr == "7A", Pos.cM. > 0.55 & Pos.cM. < 0.75)
  chr7B <- non_synon_coding %>% filter(Chr == "7B", Pos.cM. > 0.45 & Pos.cM. < 0.65)
  chr7D <- non_synon_coding %>% filter(Chr == "7D", Pos.cM. > 0.74 & Pos.cM. < 0.94)
  
  non_synon_coding_sub <- rbind(chr1A, chr1B, chr1D,
                                chr2A, chr2B, chr2D,
                                chr3A, chr3B, chr3D,
                                chr4A, chr4B, chr4D,
                                chr5A, chr5B, chr5D,
                                chr6A, chr6B, chr6D,
                                chr7A, chr7B, chr7D)  %>% arrange (X)
  non_synon_coding_sub <- non_synon_coding_sub %>% ungroup() %>% slice_sample(n = 2827) %>% arrange(X)
  DV_sub <- rbind(high_sub, non_synon_coding_sub) 
  # save loci for SnpChip column in `ASR_input_info`
  DV_loci <- DV_sub$X %>% as.data.frame()
  colnames(DV_loci) <- "X"
  DV_loci$Type <-  rep("DV_SNP", nrow(DV_sub))
  
  ### first selects what loci get genotyped (identical for traditional and DV approach)
  ### next selects what loci are QTL (random for traditional approach and DV ratio for DV approach)
  # select a single transcript to represent variant impact
  founder_data_transcript_filter <- founder_data %>%  group_by(info) %>% sample_n(1) %>% arrange(X)
  
  #### set up the loci for SnpChip (addSnpChip())
  # select 3 markers from each cM bin (on average 300 SNPs per chr)
  # founder_data_cM_filter <- founder_data_transcript_filter %>% group_by(Chr, Pos.cM.) %>% slice_sample(n=3) %>% arrange(X)
  # # save loci for SnpChip column in `ASR_input_info`
  # SnpChip_loci <- founder_data_cM_filter$X %>% as.data.frame()
  # colnames(SnpChip_loci) <- "X"
  # SnpChip_loci$Type <-  rep("SnpChip", nrow(founder_data_cM_filter))
  
  ###############################################################################
  
  #### set up the loci for QTL (addTraitA())
  # first for random QTL, output has `_R` suffix
  # next for DV QTL, output has `_DV` suffix
  
  ###############################################################################
  ##### Random QTL ##### 
    # select QTL_n markers per chromosome 
    QTL_data_cM_filter <- founder_data_transcript_filter %>% group_by(Chr) %>% slice_sample(n=QTL_n) %>% arrange(X)
    # save loci for QTL column in `ASR_input_info`
    QTL_loci <- QTL_data_cM_filter$X %>% as.data.frame()
    colnames(QTL_loci) <- "X"
    QTL_loci$Type <-  rep("QTL", nrow(QTL_data_cM_filter)) 
    
    ### combine founder_data_cM_filter and QTL_data_cM_filter for create genMap and haplotype ASR input
    ASR_input_data <- QTL_data_cM_filter #rbind(founder_data_cM_filter, QTL_data_cM_filter, DV_sub) %>% arrange(X)
    # there may be some loci shared between SnpChip and QTL - reduce duplicates from ASR_input_data
    ASR_input_data <- ASR_input_data %>% distinct()
    
    # combine the loci info
    # drop SnpChip from ASR_loci to maintain number of QTL / chr 
    # drop the DV_loci from both SnpChip and QTL to maintain random assignment
    # filter SnpChip_loci for X not present in QTL_loci
    # SnpChip_loci <- SnpChip_loci[!SnpChip_loci$X %in% QTL_loci$X, ] 
    # DV_loci_R <- DV_loci[!DV_loci$X %in% SnpChip_loci$X, ]
    # DV_loci_R <- DV_loci_R[!DV_loci_R$X %in% QTL_loci$X, ] 
    ASR_loci <- QTL_loci # rbind(SnpChip_loci, QTL_loci, DV_loci_R) %>% arrange(X)

    # genetic map is a list for each chromosome of the "Pos.cM." column
    ASR_input_data_split <- split(ASR_input_data, ASR_input_data$Chr)
    
    # jitters positions 
    jitter_Pos.cM. <- list()
    for(i in 1:length(ASR_input_data_split)){
      # add slight jitters to genetic map
      jitter_map <- jitter(ASR_input_data_split[[i]]$Pos.cM., factor = 1, amount = 0.00001) %>% as.data.frame()
      # make sure order is ascending 
      jitter_map <- jitter_map[order(jitter_map$.),] %>% as.data.frame()
      # make sure first site is 0, and not negative 
      jitter_map[1,1] <- 0
      # and if multiple 0 cM sites then remainder could be negative, prevent this with abs()
      jitter_map <- abs(jitter_map)
      # and now those values might not be ascending
      jitter_map <- jitter_map[order(jitter_map$.),] %>% as.data.frame()
      jitter_Pos.cM.[[i]] <- jitter_map
    }
    # get cM pos for genMap
    genMap <- list()
    for(i in 1:length(jitter_Pos.cM.)){
      genMap[[i]] <- jitter_Pos.cM.[[i]]$.
    }
    # output genMap variable
    genMap_input_R <<- genMap
    
    
    #### haplotype list of df where each row represents cultivar, column represent segregating sites 
    # store the marker info (snp or qtl, or DV) for assigning in ASR
    # add marker_pos column, because haplotypes will just be a numbered list 1:last marker that can correspond to the X and info column
    snp_info <- list()
    for(i in 1:length(ASR_input_data_split)){
      marker_pos <- 1:nrow(ASR_input_data_split[[i]]) %>%  as.data.frame()
      new_info <- cbind(marker_pos, ASR_input_data_split[[i]])
      new_info <- rename(new_info, "marker_pos" = `.`)
      snp_info[[i]] <- new_info
    }
    # convert to df
    founder_data_filter_df <- do.call(rbind, lapply(snp_info, as.data.frame))
    # select col "marker_pos", "X", "info", "Chr", "Pos.cM." and "SnpEff" 
    snp_info_snp_filter <- cbind(founder_data_filter_df[,1:5], founder_data_filter_df[,16]) 
    # add the "Type" col from ASR_loci
    snp_info_snp_filter <- left_join(snp_info_snp_filter, ASR_loci, by = "X")
    info_ASR_input_R <<- snp_info_snp_filter %>% rename(SnpEff = `founder_data_filter_df[, 16]`)
    
    #### make haplotype file  
    # gather just chromosome (column 4) and genotype columns (17:42)
    genotype <- cbind(founder_data_filter_df[,4], founder_data_filter_df[,17:42]) 
    split_Chr_haplotype <- split(genotype, genotype$`founder_data_filter_df[, 4]`)
    # remove chr column, transpose (each column is part of the genetic map) 
    haplotype <- list()
    for(i in 1:length(split_Chr_haplotype)){
      # remove chr column
      snps <- split_Chr_haplotype[[i]][,2:27]
      # transpose
      snps_t <- t(snps) # %>% as.data.frame()
      # set the snp position names as their order on each chr. (corresponds to snp_info_snp_filter$marker_pos)
      marker_pos_col_name <- c(1:ncol(snps_t))
      colnames(snps_t) <- marker_pos_col_name
      haplotype[[i]] <- snps_t
    }
    # save cultivar column order for reference 
    cultivar_info <<- rownames(haplotype[[1]])
    haplotype_input_R <<- haplotype
    
    ###############################################################################
    ##### Deleterious variant QTL ##### 
    # assign DV_sub as QTL sites 
    
    DV_sub_QTL <- DV_sub %>%
      arrange(X) %>% group_by(Chr) %>% 
      slice_sample(n=QTL_n) %>% 
      arrange(X)
    
    # not every chr has QTL_n variants ... too few DV on some, esp. D genome --> need to tell ASR how many QTL / chr
    
    # save loci for QTL column in `ASR_input_info`
    QTL_loci_dv <- DV_sub_QTL$X %>% as.data.frame()
    colnames(QTL_loci_dv) <- "X"
    QTL_loci_dv$Type <-  rep("QTL", nrow(DV_sub_QTL)) 
    
    
    # need to use the same snpchip genetic map positions (which came from jitter)
    consensus_pos <- do.call(rbind, lapply(jitter_Pos.cM., as.data.frame))
    colnames(consensus_pos) <- "Pos.jitter"
    consensus_SnpChip <- cbind(info_ASR_input_R, consensus_pos) %>% relocate("Pos.jitter", .before = SnpEff) %>% filter(Type == "SnpChip") 
  
    ### combine founder_data_cM_filter and DV_sub_QTL for create genMap and haplotype ASR input
    # founder_data_cM_filter needs to match consensus_SnpChip
    #founder_data_cM_filter_dv <- founder_data_cM_filter[founder_data_cM_filter$X %in% consensus_SnpChip$X, ] 
    # assign the Pos.jitter from consensus_SnpChip
    #founder_data_cM_filter_dv$Pos.cM. <- as.double(consensus_SnpChip$Pos.jitter)
    # combine with DV QTL, and extra DV sites from DV_sub
    # first remove DV_sub present in DV sites
    #high_sub_DV <- DV_sub[!DV_sub$X %in% DV_sub_QTL$X, ] 
    ASR_input_data_dv <- DV_sub_QTL 
    # there may be some loci shared between SnpChip and QTL - reduce duplicates from ASR_input_data_dv
    # slightly different from _R approach because jitter pos. all rows are not identical, only X column
    ASR_input_data_dv = ASR_input_data_dv[!duplicated(ASR_input_data_dv$X),] 
    
    # combine the loci info
    # drop QTL from ASR_loci to maintain same SnpChip for _R and _DV
    # filter QTL_loci  for X not present in SnpChip_loci
    # DO THIS FOR HIGH SNP FILTERED OUT DUE TO nQTL LIMIT
    # QTL_loci_dv <- QTL_loci_dv[!QTL_loci_dv$X %in% SnpChip_loci$X, ] 
    # DV_loci_dv <- DV_loci[!DV_loci$X %in% QTL_loci_dv$X, ] 
    # DV_loci_dv <- DV_loci_dv[!DV_loci_dv$X %in% SnpChip_loci$X, ]
    ASR_loci_dv <- QTL_loci_dv #rbind(SnpChip_loci, QTL_loci_dv, DV_loci_dv) %>% arrange(X)
    
    # genetic map is a list for each chromosome of the "Pos.cM." column
    ASR_input_data_split_dv <- split(ASR_input_data_dv, ASR_input_data_dv$Chr)
    
    #### CURRENTLY CAN HAVE THE SAME SnpChip MAP BUT DON'T HAVE SAME # OF QTL / CHR on all chr (drop some on _R)
    
    # get cM pos for genMap
    genMap <- list()
    for(i in 1:length(ASR_input_data_split_dv)){
      genMap[[i]] <- ASR_input_data_split_dv[[i]]$Pos.cM.
    }
    # output genMap variable
    genMap_input_DV <<- genMap
    
    #### haplotype list of df where each row represents cultivar, column represent segregating sites 
    # store the marker info (snp or qtl) for assigning in ASR
    # add marker_pos column, because haplotypes will just be a numbered list 1:last marker that can correspond to the X and info column
    snp_info <- list()
    for(i in 1:length(ASR_input_data_split_dv)){
      marker_pos <- 1:nrow(ASR_input_data_split_dv[[i]]) %>%  as.data.frame()
      new_info <- cbind(marker_pos, ASR_input_data_split_dv[[i]])
      new_info <- rename(new_info, "marker_pos" = `.`)
      snp_info[[i]] <- new_info
    }
    # convert to df
    founder_data_filter_df <- do.call(rbind, lapply(snp_info, as.data.frame))
    # select col "marker_pos", "X", "info", "Chr", "Pos.cM." and "SnpEff" 
    snp_info_snp_filter <- cbind(founder_data_filter_df[,1:5], founder_data_filter_df[,16]) 
    # add the "Type" col from ASR_loci_dv
    snp_info_snp_filter <- left_join(snp_info_snp_filter, ASR_loci_dv, by = "X")
    info_ASR_input_DV <<- snp_info_snp_filter %>% rename(SnpEff = `founder_data_filter_df[, 16]`)
    
    #### make haplotype file  
    # gather just chromosome (column 4) and genotype columns (17:42)
    genotype <- cbind(founder_data_filter_df[,4], founder_data_filter_df[,17:42]) 
    split_Chr_haplotype <- split(genotype, genotype$`founder_data_filter_df[, 4]`)
    # remove chr column, transpose (each column is part of the genetic map) 
    haplotype <- list()
    for(i in 1:length(split_Chr_haplotype)){
      # remove chr column
      snps <- split_Chr_haplotype[[i]][,2:27]
      # transpose
      snps_t <- t(snps) # %>% as.data.frame()
      # set the snp position names as their order on each chr. (corresponds to snp_info_snp_filter$marker_pos)
      marker_pos_col_name <- c(1:ncol(snps_t))
      colnames(snps_t) <- marker_pos_col_name
      haplotype[[i]] <- snps_t
    }
    # save cultivar column order for reference 
    #cultivar_info <<- rownames(haplotype[[1]]) same for both _R and _DV
    haplotype_input_DV <<- haplotype
    
}

# make sure QTL count per chr matches DV and R
# note which segregating sites in R are dropped and get a new "type" 
# ok to keep in geno and haplo because they aren't used for QTL or SnpChip - neutral sites
ASR_QTL_per_chr_CV <- function(info_ASR_input_DV, info_ASR_input_R){
  # check info_ASR_input_DV for true QTL/chr count 
  QTL_per_chr <<- info_ASR_input_DV %>% filter(Type == "QTL") %>% dplyr::select(Chr) %>% table() %>% as.data.frame() 
  # filter for just R QTL
  QTL_info_ASR_input_R <- info_ASR_input_R %>% filter(Type == "QTL")
  # loop over each chr to match freq in QTL_per_chr
  QTL_info_ASR_input_R_split <- split(QTL_info_ASR_input_R, QTL_info_ASR_input_R$Chr)
  
  freq <- QTL_per_chr$Freq
  store_slice <- list()
  for(i in 1:length(QTL_info_ASR_input_R_split)){
    split_slice <- QTL_info_ASR_input_R_split[[i]] %>% slice_sample(n = freq[i]) %>% arrange(X)
    store_slice[[i]] <- split_slice
  }
  QTL_info_ASR_input_R_filter <- do.call(rbind, lapply(store_slice, as.data.frame))
  #neutral_info_ASR_input_R <- QTL_info_ASR_input_R[!QTL_info_ASR_input_R$X %in% QTL_info_ASR_input_R_filter$X, ]
  #neutral_info_ASR_input_R$Type <- "Neutral"
  # could add filter where check if any are HIGH, and then get "DV_SNP" anno. 
  # join back with SnpChip_info
  #SNP_info_ASR_input_R <- info_ASR_input_R %>% filter(Type == "SnpChip")
  #DV_info_ASR_input_R <- info_ASR_input_R %>% filter(Type == "DV_SNP")
  info_ASR_input_R <<- QTL_info_ASR_input_R_filter 
  #rbind(QTL_info_ASR_input_R_filter, SNP_info_ASR_input_R,DV_info_ASR_input_R, neutral_info_ASR_input_R) %>% arrange(X)
}


#### This approach attempted to take the extra SNPs that were QTL in _R and get dropped from the info 
# and also drop them from the geno and haplo input - however the output as input in ASR ran into bugs
# specifically when setting the QTL loci based on marker postition, which is not in perfect numerical order 
# after removing QTL sites --- 

# ASR_QTL_per_chr <- function(info_ASR_input_DV, info_ASR_input_R, genMap_input_R, haplotype_input_R){
#   # because _DV dictates how many QTL per chr need to reduce get the # of segregating sites to match in _R
#   # this is important for assigning QTL sites, as well as keeping track of the DV allele frequency (which considers all segrgating sites)
#   # check info_ASR_input_DV for true QTL/chr count 
#   QTL_per_chr <<- info_ASR_input_DV %>% filter(Type == "QTL") %>% select(Chr) %>% table() %>% as.data.frame() 
#   # filter for just _R QTL sites
#   QTL_info_ASR_input_R <- info_ASR_input_R %>% filter(Type == "QTL")
#   # loop over each chr to match freq in QTL_per_chr
#   QTL_info_ASR_input_R_split <- split(QTL_info_ASR_input_R, QTL_info_ASR_input_R$Chr)
# 
#   freq <- QTL_per_chr$Freq
#   store_slice <- list()
#   for(i in 1:length(QTL_info_ASR_input_R_split)){
#     split_slice <- QTL_info_ASR_input_R_split[[i]] %>% slice_sample(n = freq[i]) %>% arrange(X)
#     store_slice[[i]] <- split_slice
#   }
#   QTL_info_ASR_input_R_filter <- do.call(rbind, lapply(store_slice, as.data.frame))
#   # join back with SnpChip_info
#   SNP_info_ASR_input_R <- info_ASR_input_R %>% filter(Type == "SnpChip")
#   DV_info_ASR_input_R <- info_ASR_input_R %>% filter(Type == "DV_SNP")
#   info_ASR_input_R <<- rbind(QTL_info_ASR_input_R_filter, SNP_info_ASR_input_R, DV_info_ASR_input_R) %>% arrange(X)
#   
#   # remove the sites from genMap_input_R and haplotype_input_R
#   remove_sites <- QTL_info_ASR_input_R[!QTL_info_ASR_input_R$X %in% QTL_info_ASR_input_R_filter$X, ]
#   # add the numerical chr order 
#   remove_sites_Chr <- remove_sites$Chr
#   remove_sites$ASR_Chr <- remove_sites_Chr %>% gsub("1A", 1, .) %>% 
#     gsub("1B", 2, .) %>% 
#     gsub("1D", 3, .) %>% 
#     gsub("2A", 4, .) %>%  
#     gsub("2B", 5, .) %>%
#     gsub("2D", 6, .) %>%
#     gsub("3A", 7, .) %>% 
#     gsub("3B", 8, .) %>%
#     gsub("3D", 9, .) %>%
#     gsub("4A", 10, .) %>% 
#     gsub("4B", 11, .) %>%
#     gsub("4D", 12, .) %>%
#     gsub("5A", 13, .) %>% 
#     gsub("5B", 14, .) %>%
#     gsub("5D", 15, .) %>%
#     gsub("6A", 16, .) %>% 
#     gsub("6B", 17, .) %>%
#     gsub("6D", 18, .) %>%
#     gsub("7A", 19, .) %>% 
#     gsub("7B", 20, .) %>%
#     gsub("7D", 21, .) 
#   
#   # remove_sites$marker_pos corresponds to haplotype_input column position to drop and genMap row
#   genMap_1 <- genMap_input_R[[1]] #%>% as.data.frame()
#   haplo_1 <- haplotype_input_R[[1]] #%>% as.data.frame()
#   remove_1 <- remove_sites %>% filter(ASR_Chr == 1)
#   if(nrow(remove_1) > 0){
#     genMap_1 <- genMap_1[-remove_1$marker_pos] #%>% as.data.frame()
#     haplo_1 <- haplo_1[,-remove_1$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_2 <- genMap_input_R[[2]] #%>% as.data.frame()
#   haplo_2 <- haplotype_input_R[[2]] #%>% as.data.frame()
#   remove_2 <- remove_sites %>% filter(ASR_Chr == 2)
#   if(nrow(remove_2) > 0){
#     genMap_2 <- genMap_2[-remove_2$marker_pos] #%>% as.data.frame()
#     haplo_2 <- haplo_2[,-remove_2$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_3 <- genMap_input_R[[3]] #%>% as.data.frame()
#   haplo_3 <- haplotype_input_R[[3]] #%>% as.data.frame()
#   remove_3 <- remove_sites %>% filter(ASR_Chr == 3)
#   if(nrow(remove_3) > 0){
#     genMap_3 <- genMap_3[-remove_3$marker_pos] #%>% as.data.frame()
#     haplo_3 <- haplo_3[,-remove_3$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_4 <- genMap_input_R[[4]] #%>% as.data.frame()
#   haplo_4 <- haplotype_input_R[[4]] #%>% as.data.frame()
#   remove_4 <- remove_sites %>% filter(ASR_Chr == 4)
#   if(nrow(remove_4) > 0){
#     genMap_4 <- genMap_4[-remove_4$marker_pos] #%>% as.data.frame()
#     haplo_4 <- haplo_4[,-remove_4$marker_pos]#%>% as.data.frame()
#   }
# 
#   
#   genMap_5 <- genMap_input_R[[5]] #%>% as.data.frame()
#   haplo_5 <- haplotype_input_R[[5]] #%>% as.data.frame()
#   remove_5 <- remove_sites %>% filter(ASR_Chr == 5)
#   if(nrow(remove_5) > 0){
#     genMap_5 <- genMap_5[-remove_5$marker_pos] #%>% as.data.frame()
#     haplo_5 <- haplo_5[,-remove_5$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_6 <- genMap_input_R[[6]] #%>% as.data.frame()
#   haplo_6 <- haplotype_input_R[[6]] #%>% as.data.frame()
#   remove_6 <- remove_sites %>% filter(ASR_Chr == 6)
#   if(nrow(remove_6) > 0){
#     genMap_6 <- genMap_6[-remove_6$marker_pos] #%>% as.data.frame()
#     haplo_6 <- haplo_6[,-remove_6$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_7 <- genMap_input_R[[7]] #%>% as.data.frame()
#   haplo_7 <- haplotype_input_R[[7]] #%>% as.data.frame()
#   remove_7 <- remove_sites %>% filter(ASR_Chr == 7)
#   if(nrow(remove_7) > 0){
#     genMap_7 <- genMap_7[-remove_7$marker_pos] #%>% as.data.frame()
#     haplo_7 <- haplo_7[,-remove_7$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_8 <- genMap_input_R[[8]] #%>% as.data.frame()
#   haplo_8 <- haplotype_input_R[[8]] #%>% as.data.frame()
#   remove_8 <- remove_sites %>% filter(ASR_Chr == 8)
#   if(nrow(remove_8) > 0){
#     genMap_8 <- genMap_8[-remove_8$marker_pos] #%>% as.data.frame()
#     haplo_8 <- haplo_8[,-remove_8$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_9 <- genMap_input_R[[9]] #%>% as.data.frame()
#   haplo_9 <- haplotype_input_R[[9]] #%>% as.data.frame()
#   remove_9 <- remove_sites %>% filter(ASR_Chr == 9)
#   if(nrow(remove_9) > 0){
#     genMap_9 <- genMap_9[-remove_9$marker_pos] #%>% as.data.frame()
#     haplo_9 <- haplo_9[,-remove_9$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_10 <- genMap_input_R[[10]] #%>% as.data.frame()
#   haplo_10 <- haplotype_input_R[[10]] #%>% as.data.frame()
#   remove_10 <- remove_sites %>% filter(ASR_Chr == 10)
#   if(nrow(remove_10) > 0){
#     genMap_10 <- genMap_10[-remove_10$marker_pos] #%>% as.data.frame()
#     haplo_10 <- haplo_10[,-remove_10$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_11 <- genMap_input_R[[11]] #%>% as.data.frame()
#   haplo_11 <- haplotype_input_R[[11]] #%>% as.data.frame()
#   remove_11 <- remove_sites %>% filter(ASR_Chr == 11)
#   if(nrow(remove_11) > 0){
#     genMap_11 <- genMap_11[-remove_11$marker_pos] #%>% as.data.frame()
#     haplo_11 <- haplo_11[,-remove_11$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_12 <- genMap_input_R[[12]] #%>% as.data.frame()
#   haplo_12 <- haplotype_input_R[[12]] #%>% as.data.frame()
#   remove_12 <- remove_sites %>% filter(ASR_Chr == 12)
#   if(nrow(remove_12) > 0){
#     genMap_12 <- genMap_12[-remove_12$marker_pos] #%>% as.data.frame()
#     haplo_12 <- haplo_12[,-remove_12$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_13 <- genMap_input_R[[13]] #%>% as.data.frame()
#   haplo_13 <- haplotype_input_R[[13]] #%>% as.data.frame()
#   remove_13 <- remove_sites %>% filter(ASR_Chr == 13)
#   if(nrow(remove_13) > 0){
#     genMap_13 <- genMap_13[-remove_13$marker_pos] #%>% as.data.frame()
#     haplo_13 <- haplo_13[,-remove_13$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_14 <- genMap_input_R[[14]] #%>% as.data.frame()
#   haplo_14 <- haplotype_input_R[[14]] #%>% as.data.frame()
#   remove_14 <- remove_sites %>% filter(ASR_Chr == 14)
#   if(nrow(remove_14) > 0){
#     genMap_14 <- genMap_14[-remove_14$marker_pos] #%>% as.data.frame()
#     haplo_14 <- haplo_14[,-remove_14$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_15 <- genMap_input_R[[15]] #%>% as.data.frame()
#   haplo_15 <- haplotype_input_R[[15]] #%>% as.data.frame()
#   remove_15 <- remove_sites %>% filter(ASR_Chr == 15)
#   if(nrow(remove_15) > 0){
#     genMap_15 <- genMap_15[-remove_15$marker_pos] #%>% as.data.frame()
#     haplo_15 <- haplo_15[,-remove_15$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_16 <- genMap_input_R[[16]] #%>% as.data.frame()
#   haplo_16 <- haplotype_input_R[[16]] #%>% as.data.frame()
#   remove_16 <- remove_sites %>% filter(ASR_Chr == 16)
#   if(nrow(remove_16) > 0){
#     genMap_16 <- genMap_16[-remove_16$marker_pos] #%>% as.data.frame()
#     haplo_16 <- haplo_16[,-remove_16$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_17 <- genMap_input_R[[17]] #%>% as.data.frame()
#   haplo_17 <- haplotype_input_R[[17]] #%>% as.data.frame()
#   remove_17 <- remove_sites %>% filter(ASR_Chr == 17)
#   if(nrow(remove_17) > 0){
#     genMap_17 <- genMap_17[-remove_17$marker_pos] #%>% as.data.frame()
#     haplo_17 <- haplo_17[,-remove_17$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_18 <- genMap_input_R[[18]] #%>% as.data.frame()
#   haplo_18 <- haplotype_input_R[[18]] #%>% as.data.frame()
#   remove_18 <- remove_sites %>% filter(ASR_Chr == 18)
#   if(nrow(remove_18) > 0){
#     genMap_18 <- genMap_18[-remove_18$marker_pos] #%>% as.data.frame()
#     haplo_18 <- haplo_18[,-remove_18$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_19 <- genMap_input_R[[19]] #%>% as.data.frame()
#   haplo_19 <- haplotype_input_R[[19]] #%>% as.data.frame()
#   remove_19 <- remove_sites %>% filter(ASR_Chr == 19)
#   if(nrow(remove_19) > 0){
#     genMap_19 <- genMap_19[-remove_19$marker_pos] #%>% as.data.frame()
#     haplo_19 <- haplo_19[,-remove_19$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_20 <- genMap_input_R[[20]] #%>% as.data.frame()
#   haplo_20 <- haplotype_input_R[[20]] #%>% as.data.frame()
#   remove_20 <- remove_sites %>% filter(ASR_Chr == 20)
#   if(nrow(remove_20) > 0){
#     genMap_20 <- genMap_20[-remove_20$marker_pos] #%>% as.data.frame()
#     haplo_20 <- haplo_20[,-remove_20$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_21 <- genMap_input_R[[21]] #%>% as.data.frame()
#   haplo_21 <- haplotype_input_R[[21]] #%>% as.data.frame()
#   remove_21 <- remove_sites %>% filter(ASR_Chr == 21)
#   if(nrow(remove_21) > 0){
#     genMap_21 <- genMap_21[-remove_21$marker_pos] #%>% as.data.frame()
#     haplo_21 <- haplo_21[,-remove_21$marker_pos]#%>% as.data.frame()
#   }
#   
#   genMap_input_R <<- list(genMap_1, genMap_2, genMap_3, genMap_4, genMap_5, genMap_6, genMap_7, genMap_8, genMap_9, genMap_10,
#                genMap_11, genMap_12, genMap_13, genMap_14, genMap_15, genMap_16, genMap_17, genMap_18, genMap_19, genMap_20, genMap_21)
#   
#   haplotype_input_R <<- list(haplo_1, haplo_2, haplo_3, haplo_4, haplo_5, haplo_6, haplo_7, haplo_8, haplo_9, haplo_10,
#                           haplo_11, haplo_12, haplo_13, haplo_14, haplo_15, haplo_16, haplo_17, haplo_18, haplo_19, haplo_20, haplo_21)
#   
#   
#   # it is OK that info_ASR_input_R and info_ASR_input_DV are different lengths because table($Type) shows QTL count and SnpChip are the same
#   # the DV_SNP is driving the difference, which is nececarry for tracking DV AF at all segregating sites  
#   
#   
#   # # lazy, but noticing here that info_ASR_input_R$pos.cM is true cM, and info_ASR_input_DV$pos.cM is jitter ONLY for SnpChip pos
#   # # easiest to get the jitter pos for info_ASR_input_R from genMap_input_R
#   # pos.Jitter <- do.call(rbind, lapply(genMap_input_R, as.data.frame))
#   # 
#   # # easiest to get the true cM pos for info_ASR_input_DV from SNP_info_ASR_input_R and QTL_info_ASR_input_DV
#   # 
#   # # info_ASR_input_R needs SNP_info_ASR_input_DV jitter pos combined with QTL_info_ASR_input_R cM positions 
#   # and QTL_info_ASR_input_DV positions 
#   # #add the get the pos.cM and pos.Jitter positions for info_ASR_input_R and info_ASR_input_DV
#   # # info_ASR_input_R has pos.cM, needs pos.Jitter
#   # info_ASR_input_R$pos.Jitter <- 
#   # info_ASR_input_DV
# }
