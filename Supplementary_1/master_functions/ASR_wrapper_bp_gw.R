##### 1.20.22
# update: the QTL effects are kept constant across the different R and DV studies
##### Wrapper function for creating ASR founder haplotypes, setting simulation parameters, and modeling the breeding program

ASR_wrapper_bp <- function(type, map_scale, peri_scale, genMap_input, haplotypes_input, info_ASR_input, dv_loci, QTL_effects, QTL_effect_ratio,
                        heritability, legend){
  ######### Wrapper function for creating ASR founder haplotypes, setting simulation parameters, and modeling the breeding program
  # note QTL # is selected prior to this function input
  # type is "R" (random SNP assigned as QTL) or "DV" (deleterious SnpEff assigned as QTL)
  # map_scale is numeric, by how much do you want to change map size? 1 is WT
  # peri_scale is numeric, by how much do you want to scale DV 20 cM region? 1 is WT
  # genMap_input is the ASR formatted input list of SNP positions, an output of ASR_input()
  # haplotypes_input is the ASR formatted input list of genotypes, an output of ASR_input()
  # QTL_effects is a list of QTL effects to keep constant among the different studies
  # info_ASR_input relates to _R or _DV, for choosing specific SNP and QTL sites, an output of ASR_input()
  # QTL_effect_ratio is formatted as "alternate" or c(X, Y) where
  ##### alternate has QTL effect go +/-/+/- etc. 
  ##### X is ratio of + QTL effect (select against minor allele)
  ##### Y is ratio of - QTL effect (select against major allele)
  ##### Values between 0 and 1, ie. c(0.5, 0.5) is half + and half -
  # heritability is a numeric, ie. 0.5
  # legend provides DF name info column for identifying parameters, used in plotting / filtering
  
  ######### Create the founder populations with newMapPops() ########## 
  founderPop = newMapPop(genMap = genMap_input, haplotypes = haplotypes_input, inbred = TRUE, ploidy = 2L)
  
  # set the simulation parameters
  SP = SimParam$new(founderPop)
  
  ########## manually add trait (QTL) and SnpChip (genotype) positions based on ASR_input() results  
  
  # set up additive trait QTL (yield)
  # assign QTL_n loci for QTl per chr
  SP$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
  # Extract the trait
  trait = SP$traits[[1]]
  
  # loci we know we want to assign to QTL 
  QTL_loci <- info_ASR_input %>% filter(Type == "QTL") %>% select(marker_pos)
  # change loci
  trait@lociLoc = QTL_loci$marker_pos
  # change the QTL_effects to be constant among studies
  trait@addEff = QTL_effects 
  # change the additive effects signs based on QTL_effect_ratio
  # "alternate" means 1/-1/1/-1 etc.
  if(QTL_effect_ratio[1] == "alternate"){
    # make everything absolute value first
    trait@addEff = abs(trait@addEff)
    change_effect = length(trait@addEff)
    # new_effect_ratio alternates +/- for each QTL
    new_effect_ratio <- c(rep(x = c(1, -1), length.out = change_effect)) 
    # length(new_effect_ratio) may not equal change_effect
    if(length(new_effect_ratio) == change_effect){
      trait@addEff = trait@addEff*new_effect_ratio
    } 
    # or instead of "alternate", could be a ratio that's not fully random, and has some ratio of selecting against minor allele
  } else if(QTL_effect_ratio[1] != "alternate"){ # if not fully random
    # QTL_effect_ratio[1] is ratio that's +
    # QTL_effect_ratio[2] is ratio that -
    # make everything absolute value first
    trait@addEff = abs(trait@addEff)
    change_effect = length(trait@addEff)
    # new_effect_ratio takes half of QTL_effect_ratio[1] as + and -, and all of QTL_effect_ratio[2] is +
    new_effect_ratio <- c(rep.int(1, times = round((QTL_effect_ratio[1])*change_effect)), 
                          rep.int(-1, times = round(QTL_effect_ratio[2]*change_effect)))
    new_effect = sample(new_effect_ratio, size = length(new_effect_ratio), replace = FALSE)
    if(length(new_effect) == change_effect){
      trait@addEff = trait@addEff*new_effect
    } else if (length(new_effect) - change_effect < 0) { # if new_effect too small, add remaining QTL_effects
      difference = length(new_effect) - change_effect
      new_effect = c(new_effect, 
                     sample(new_effect_ratio, size = abs(difference), replace = TRUE))
      trait@addEff = trait@addEff*new_effect
    } else if (length(new_effect) - change_effect > 0) { # if new_effect too large, remove QTL_effects
      difference = length(new_effect) - change_effect 
      new_effect <- new_effect[1:(length(new_effect)-abs(difference))]
      trait@addEff = trait@addEff*new_effect
    }
  }
  
  # Replace trait loci
  SP$switchTrait(traitPos=1, lociMap=trait)
  # NOTE WE COME BACK AND RE-ASSIGN
  # need the variance of the founder pop to correct for different QTL effect dist. under different conditions
  
  # newPop(founderPop) takes into account SP settings
  parents = newPop(founderPop, simParam = SP) 
  parents_Var <- varA(pop = parents, simParam = SP)
  
  # corrected QTL effects to set the same variance among conditons
  corrected_addEff <- trait@addEff / as.vector(sqrt(parents_Var))
  # Replace QTL additive effects 
  trait@addEff <- corrected_addEff
  # Replace trait loci
  SP$resetPed()
  SP$switchTrait(traitPos=1, lociMap=trait)
  
  # set up SnpChip (genotypes)
  # assign SnpChip_n loci for SNPs per chr 
  SnpChip_n <- table(info_ASR_input$Chr, info_ASR_input$Type) %>% as.data.frame() %>% filter(Var2 == "SnpChip") %>% select(Freq)
  SP$addSnpChip(SnpChip_n$Freq)
  # Extract the SnpChip
  genotype = SP$snpChips[[1]]
  # loci we know we want to assign to genotype
  SnpChip_loci <- info_ASR_input %>% filter(Type == "SnpChip") %>% select(marker_pos)
  # change loci
  genotype@lociLoc = SnpChip_loci$marker_pos
  # Replace genotype loci
  SP$snpChips[[1]]@lociLoc = genotype@lociLoc
  
  # set error variance with heritability ~ to yield
  SP$setVarE(H2 = heritability)
  
  # keep track of recombination 
  # SP$setTrackRec(TRUE) # not sure of utility yet ... 
  
  assign("SP", SP, envir = .GlobalEnv)
  
  # ---------------------------------------------------
  # cycles of phenotypic selection burnin
  # ---------------------------------------------------
  # newPop(founderPop) takes into account SP settings
  parents = newPop(founderPop, simParam = SP) 
  parents_Var <- varA(pop = parents, simParam = SP)
  parents_BE <- varA(pop = parents, simParam = SP)/genicVarA(pop = parents, simParam = SP)
  #parents_LD <- QTL_SNP_LD(pop = parents, simParam = SP, info_ASR_input = info_ASR_input)
  # QTL fixation ratio
  QTL_fixation(trait = trait, pop = parents, simParam = SP, pop_name = "founder")

  # measure the allele frequency of parents (26 cultivars)
  # parents_AF = pullSegSiteGeno(parents, simParam = SP) %>% as.data.table()
  
  #Want to create a larger training population, create F1 DHs 
  
  #nBurnInCyc <- 1 #number of burnin cycles
  nCrs_burn <- 400 #number of crosses
  nProgPerCrs <- 1 #number of progeny per cross
  nDHperProg <- 1 #number of DH from each F1
  nSelect <- 100 #number of progeny to select and advance

  burnin <- parents
    # make F1s
    burnin = randCross(pop = burnin, nCrosses = nCrs_burn,
                       nProgeny = nProgPerCrs, simParam = SP)
    # make DHs
    burnin = makeDH(pop = burnin, nDH = nDHperProg, simParam = SP)
    # measure the DH phenotypes
    burnin = setPheno(pop = burnin, H2 = heritability, simParam = SP)
    # if (cycle < nBurnInCyc){
    #   # phenotypic selection for next generation of burnin
    #   burnin = selectInd(pop = burnin, nInd = nSelect, use = "pheno")
    # }

  burnin_Var <- varA(pop = burnin, simParam = SP)
  # re-standardize variance to 1 with burnin_Var
  # corrected QTL effects to set the same variance among conditons
  corrected_addEff <- trait@addEff / as.vector(sqrt(burnin_Var))
  # Replace QTL additive effects 
  trait@addEff <- corrected_addEff
  # Replace trait loci
  SP$resetPed()
  SP$switchTrait(traitPos=1, lociMap=trait)
  
  assign("SP", SP, envir = .GlobalEnv)
  # 
  burnin_Var <- varA(pop = burnin, simParam = SP)
  burnin_BE <- varA(pop = burnin, simParam = SP)/genicVarA(pop = burnin, simParam = SP)
  #burnin_LD <- QTL_SNP_LD(pop = burnin, simParam = SP, info_ASR_input = info_ASR_input)
  # fixed QTL ratio
  QTL_fixation(trait = trait, pop = burnin, simParam = SP, pop_name = "burnin")
  # measure the allele frequency (used for DV AF later)
  #burnin_AF = pullSegSiteGeno(burnin, simParam = SP) %>% as.data.table()
  
  # trying to re-set the genetic map, burnin becomes a new founder pop ...
  # the genMap_input should be the same - positions haven't changed from WT 
  # haplotypes come from burnin 
  burnin_haplo <- pullSegSiteGeno(burnin) %>% as.data.frame()
  # need to turn columns for each chr into lists
  chr_length <- NULL
  for(i in 1:length(haplotypes_input)){
    SNP_count <- ncol(haplotypes_input[[i]])
    chr_length <- c(chr_length, SNP_count)
  }
  chr_length <- cumsum(chr_length) + c(0, rep(1, 19), 0)
  # pull chr columns for list object
  chr1A <- burnin_haplo[, 1:chr_length[1]]
  chr1B <- burnin_haplo[, (chr_length[1]+1):(chr_length[2]-1)]
  chr1D <- burnin_haplo[, (chr_length[2]):(chr_length[3]-1)]
  chr2A <- burnin_haplo[, (chr_length[3]):(chr_length[4]-1)]
  chr2B <- burnin_haplo[, (chr_length[4]):(chr_length[5]-1)]
  chr2D <- burnin_haplo[, (chr_length[5]):(chr_length[6]-1)]
  chr3A <- burnin_haplo[, (chr_length[6]):(chr_length[7]-1)]
  chr3B <- burnin_haplo[, (chr_length[7]):(chr_length[8]-1)]
  chr3D <- burnin_haplo[, (chr_length[8]):(chr_length[9]-1)]
  chr4A <- burnin_haplo[, (chr_length[9]):(chr_length[10]-1)]
  chr4B <- burnin_haplo[, (chr_length[10]):(chr_length[11]-1)]
  chr4D <- burnin_haplo[, (chr_length[11]):(chr_length[12]-1)]
  chr5A <- burnin_haplo[, (chr_length[12]):(chr_length[13]-1)]
  chr5B <- burnin_haplo[, (chr_length[13]):(chr_length[14]-1)]
  chr5D <- burnin_haplo[, (chr_length[14]):(chr_length[15]-1)]
  chr6A <- burnin_haplo[, (chr_length[15]):(chr_length[16]-1)]
  chr6B <- burnin_haplo[, (chr_length[16]):(chr_length[17]-1)]
  chr6D <- burnin_haplo[, (chr_length[17]):(chr_length[18]-1)]
  chr7A <- burnin_haplo[, (chr_length[18]):(chr_length[19]-1)]
  chr7B <- burnin_haplo[, (chr_length[19]):(chr_length[20]-1)]
  chr7D <- burnin_haplo[, (chr_length[20]):(chr_length[21])]

  burnin_haplo <- list(chr1A, chr1B, chr1D, 
          chr2A, chr2B, chr2D,
          chr3A, chr3B, chr3D,
          chr4A, chr4B, chr4D,
          chr5A, chr5B, chr5D,
          chr6A, chr6B, chr6D,
          chr7A, chr7B, chr7D) 
  
  founderPop2 = newMapPop(genMap = genMap_input, haplotypes = burnin_haplo, inbred = TRUE, ploidy = 2L)
  
  ########## re-set genetic map size 
  if(map_scale > 1){
    genMap_change = founderPop2@genMap
    
    # Multiply genetic length by map_scale
    my_genMap_change = lapply(genMap_change, function(x) map_scale*x)
    
    # Replace genetic map
    founderPop2@genMap = my_genMap_change
  }
  
  if(peri_scale > 1){
    my_genMap_change = founderPop2@genMap
    
    # call DV_scale_map
    DV_scale_map(genMap_change = my_genMap_change, DV_scale = peri_scale, BP = TRUE)
    
    # Replace genetic map
    founderPop2@genMap = genMap_change_output
  }
  
  ######## remake SP as SP2
  SP2 = SimParam$new(founderPop2)
  
  ########## manually add trait (QTL) and SnpChip (genotype) positions based on ASR_input() results  
  
  # set up additive trait QTL (yield)
  # assign QTL_n loci for QTl per chr
  SP2$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
  # Extract the trait
  trait = SP2$traits[[1]]
  
  # loci we know we want to assign to QTL 
  QTL_loci <- info_ASR_input %>% filter(Type == "QTL") %>% select(marker_pos)
  # change loci
  trait@lociLoc = QTL_loci$marker_pos
  # change the QTL_effects to be constant among studies
  trait@addEff = QTL_effects 
  # change the additive effects signs based on QTL_effect_ratio
  # "alternate" means 1/-1/1/-1 etc.
  if(QTL_effect_ratio[1] == "alternate"){
    # make everything absolute value first
    trait@addEff = abs(trait@addEff)
    change_effect = length(trait@addEff)
    # new_effect_ratio alternates +/- for each QTL
    new_effect_ratio <- c(rep(x = c(1, -1), length.out = change_effect)) 
    # length(new_effect_ratio) may not equal change_effect
    if(length(new_effect_ratio) == change_effect){
      trait@addEff = trait@addEff*new_effect_ratio
    } 
    # or instead of "alternate", could be a ratio that's not fully random, and has some ratio of selecting against minor allele
  } else if(QTL_effect_ratio[1] != "alternate"){ # if not fully random
    # QTL_effect_ratio[1] is ratio that's +
    # QTL_effect_ratio[2] is ratio that -
    # make everything absolute value first
    trait@addEff = abs(trait@addEff)
    change_effect = length(trait@addEff)
    # new_effect_ratio takes half of QTL_effect_ratio[1] as + and -, and all of QTL_effect_ratio[2] is +
    new_effect_ratio <- c(rep.int(1, times = round((QTL_effect_ratio[1])*change_effect)), 
                          rep.int(-1, times = round(QTL_effect_ratio[2]*change_effect)))
    new_effect = sample(new_effect_ratio, size = length(new_effect_ratio), replace = FALSE)
    if(length(new_effect) == change_effect){
      trait@addEff = trait@addEff*new_effect
    } else if (length(new_effect) - change_effect < 0) { # if new_effect too small, add remaining QTL_effects
      difference = length(new_effect) - change_effect
      new_effect = c(new_effect, 
                     sample(new_effect_ratio, size = abs(difference), replace = TRUE))
      trait@addEff = trait@addEff*new_effect
    } else if (length(new_effect) - change_effect > 0) { # if new_effect too large, remove QTL_effects
      difference = length(new_effect) - change_effect 
      new_effect <- new_effect[1:(length(new_effect)-abs(difference))]
      trait@addEff = trait@addEff*new_effect
    }
  }
  
  # Replace trait loci
  SP2$switchTrait(traitPos=1, lociMap=trait)
  # NOTE WE COME BACK AND RE-ASSIGN
  # need the variance of the founder pop to correct for different QTL effect dist. under different conditions
  
  # newPop(founderPop2) takes into account SP2 settings
  parents2 = newPop(founderPop2, simParam = SP2) 
  parents2_Var <- varA(pop = parents2, simParam = SP2)
  
  # corrected QTL effects to set the same variance among conditons
  corrected_addEff2 <- trait@addEff / as.vector(sqrt(parents2_Var))
  # Replace QTL additive effects 
  trait@addEff <- corrected_addEff2
  # Replace trait loci
  SP2$resetPed()
  SP2$switchTrait(traitPos=1, lociMap=trait)
  
  # set up SnpChip (genotypes)
  # assign SnpChip_n loci for SNPs per chr 
  SnpChip_n <- table(info_ASR_input$Chr, info_ASR_input$Type) %>% as.data.frame() %>% filter(Var2 == "SnpChip") %>% select(Freq)
  SP2$addSnpChip(SnpChip_n$Freq)
  # Extract the SnpChip
  genotype = SP2$snpChips[[1]]
  # loci we know we want to assign to genotype
  SnpChip_loci <- info_ASR_input %>% filter(Type == "SnpChip") %>% select(marker_pos)
  # change loci
  genotype@lociLoc = SnpChip_loci$marker_pos
  # Replace genotype loci
  SP2$snpChips[[1]]@lociLoc = genotype@lociLoc
  
  # set error variance with heritability ~ to yield
  SP2$setVarE(H2 = heritability)
  
  # keep track of recombination 
  # SP2$setTrackRec(TRUE) # not sure of utility yet ... 
  
  assign("SP2", SP2, envir = .GlobalEnv)
  
 
  # ---------------------------------------------------
  # set TP and GS parameters (geno and pheno from last generation of burnin)
  # ---------------------------------------------------
  # newPop(founderPop) takes into account SP settings
  # founderPop2 = burnin with genetic map changes
  TP = newPop(founderPop2, simParam = SP2) 
  # set phenotypes
  TP = setPheno(pop = TP, H2 = heritability, simParam = SP2)
  
  gsPop <- TP 
  #founder_AF = pullSegSiteGeno(gsPop, simParam = SP) %>% as.data.table() # used for tracking DV AF beginning to end
  
  nGScyc <- 10
  nCrs_GS <- 400
  nToUpdateTP <- 160 # top 40% return to TP
  nToDropTP <- 80 # drop bottom 20% of TP
  nToCross <- 20 # top 5% selection intensity 
  
  # track population mean and variation
  allMeans <- NULL
  allVars <- NULL
  allAcuracy <- NULL 
  allBulmer <- NULL
  #allLD <- NULL
  
  for (cycle in 1:nGScyc){
    F1 = randCross(pop = gsPop, nCrosses = nCrs_GS, 
                   nProgeny = nProgPerCrs, simParam = SP2) 
    allMeans = c(allMeans, meanG(F1))
    allVars = c(allVars, varA(F1, simParam=SP2))
    # Make DHs
    DH = makeDH(pop = F1, nDH = nDHperProg, simParam = SP2) 
    ans = RRBLUP(TP, simParam=SP2)
    DH = setEBV(DH, ans, simParam=SP2)
    #Evaluate accuracy
    allAcuracy <- c(allAcuracy, cor(gv(DH), ebv(DH)))
    # Evaluate Bulmer effect, LD and fixed QTL
    be <- varA(pop = DH, simParam = SP2)/genicVarA(pop = DH, simParam = SP2)
    allBulmer <- c(allBulmer, be)
    # LD <- QTL_SNP_LD(pop = DH, simParam = SP2, info_ASR_input = info_ASR_input)
    # allLD <- c(allLD, LD)
    QTL_fixation(trait = trait, pop = DH, simParam = SP2, pop_name = cycle)
    ###JLJ update the TP here (there could be different ways of doing it)
    updateTPpop <- selectInd(DH, nInd = nToUpdateTP, trait = 1,
                              use = "ebv", simParam = SP2)
    # # fully replace the TP with new generation
    #updateTPpop <- DH
    updateTPpop <- setPheno(updateTPpop, H2 = heritability, simParam = SP2)
    ###JLJ the TP will get big fast: you may want to remove the older DHs
    TP <- c(TP, updateTPpop)
    # and drop the bottom 20%
    TP <- selectInd(TP, nInd = (TP@nInd - nToDropTP), trait = 1,
                             use = "pheno", simParam = SP2)
    # fully replace the TP
    #TP <- updateTPpop
    ###JLJ now select the DHs that will be crossed
    gsPop = selectInd(DH, nInd = nToCross, trait = 1, use = "ebv", simParam = SP2)
  }
  
  # compare AF to Founder_AF, used for tracking DV AF beginning to end
  #end_AF = pullSegSiteGeno(gsPop, simParam = SP2) %>% as.data.table()
  
  cycle = 1:nGScyc
  genetic_gain = cbind(cycle, allMeans) %>% as.data.frame()
  # note using what's called cycle0 instead of burnin gv
  
  # for each cycle, genetic gain is difference between mean gv @ cycle 0 and cycle X
  genetic_gain$genetic_gain = genetic_gain$allMeans - genetic_gain[1,2]
  genetic_gain$legend = legend
  colnames(genetic_gain) = c("cycle", "pop_mean_gv", "genetic_gain", "legend")
  assign(paste0("GG_", legend), genetic_gain, envir = .GlobalEnv)
  
  genetic_variance = cbind(cycle, allVars) %>% as.data.frame()
  # add parent and burnin variance
  parents_Var = data.frame("founder", parents_Var)
  colnames(parents_Var) = c("cycle", "allVars")
  burnin_Var = data.frame("burnin", burnin_Var)
  colnames(burnin_Var) = c("cycle", "allVars")
  genetic_variance = rbind(parents_Var, burnin_Var, genetic_variance)
  genetic_variance$legend = legend
  colnames(genetic_variance) = c("cycle", "additive_variance", "legend")
  assign(paste0("GV_", legend), genetic_variance, envir = .GlobalEnv)
  
  prediction_accuracy = cbind(cycle, allAcuracy) %>% as.data.frame()
  prediction_accuracy$legend = legend
  colnames(prediction_accuracy) = c("cycle", "prediction_accuracy", "legend")
  assign(paste0("PA_", legend), prediction_accuracy, envir = .GlobalEnv)
  
  bulmer_effect = cbind(cycle, allBulmer) %>% as.data.frame()
  # add parent and burnin BE
  parents_BE = data.frame("founder", parents_BE)
  colnames(parents_BE) = c("cycle", "allBulmer")
  burnin_BE = data.frame("burnin", burnin_BE)
  colnames(burnin_BE) = c("cycle", "allBulmer")
  bulmer_effect = rbind(parents_BE, burnin_BE, bulmer_effect)
  bulmer_effect$legend = legend
  colnames(bulmer_effect) = c("cycle", "bulmer_effect", "legend")
  assign(paste0("BE_", legend), bulmer_effect, envir = .GlobalEnv)
  
  # LD_QTL_SNP = cbind(cycle, allLD) %>% as.data.frame()
  # # add parent and burnin variance
  # parents_LD = data.frame("founder", parents_LD)
  # colnames(parents_LD) = c("cycle", "allLD")
  # burnin_LD = data.frame("burnin", burnin_LD)
  # colnames(burnin_LD) = c("cycle", "allLD")
  # LD_QTL_SNP = rbind(parents_LD, burnin_LD, LD_QTL_SNP)
  # LD_QTL_SNP$legend = legend
  # colnames(LD_QTL_SNP) = c("cycle", "LD_QTL_SNP", "legend")
  # assign(paste0("LD_", legend), LD_QTL_SNP, envir = .GlobalEnv)
  
  allQTL <- rbind(founder_QTL, burnin_QTL, `1_QTL`, `2_QTL`, `3_QTL`, `4_QTL`, `5_QTL`, `6_QTL`,
                  `7_QTL`, `8_QTL`, `9_QTL`, `10_QTL`) %>% as.data.frame()
  allQTL$legend = legend
  colnames(allQTL) = c("better_QTL_fix_ratio", "worse_QTL_fix_ratio", "cycle", "legend")
  assign(paste0("QTL_", legend), allQTL, envir = .GlobalEnv)
  
  # # QTL effect size AF change
  # QTL_effect_size_AF_change(trait = trait, pop_1 = burnin, pop_2 = gsPop, simParam = SP2)
  # QTL_ES_AF_change$legend = legend
  # assign(paste0("QTL_ES_", legend), QTL_ES_AF_change, envir = .GlobalEnv)
  
  # QTL AF change
  my_QTL_AF_change(trait = trait, pop_1 = burnin, pop_2 = gsPop, simParam = SP2)
  QTL_AF_change$legend = legend
  assign(paste0("QTL_AF_", legend), QTL_AF_change, envir = .GlobalEnv)
  
}

