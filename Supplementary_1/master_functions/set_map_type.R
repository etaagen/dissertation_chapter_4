# function for 10X genetic map in 20 cM with most DV

DV_scale_map <- function(genMap_change, DV_scale, BP){
    genMap_1 <- genMap_change[[1]]
    # get left distal positions
    l_distal_1 <- genMap_1[genMap_1 < 0.35] 
    # get pericentromere positions
    peri_1 <- genMap_1[genMap_1 > 0.35 & genMap_1 < 0.55] 
    # multiply by DV_scale 
    peri_1 <- peri_1*DV_scale
    # get right distal positions
    r_distal_1 <- genMap_1[genMap_1 > 0.55] 
    # adjust for difference between peri_1[length(peri_1)] and r_distal_1[1], and distance from l_distal
    r_distal_1 <- r_distal_1 + peri_1[length(peri_1)] + peri_1[1]
    # combine
    if(length(l_distal_1) == 0 | length(peri_1) == 0 | length(r_distal_1) == 0){
        genMap_1 <- genMap_1
    } else{
        genMap_1 <- c(l_distal_1, peri_1, r_distal_1)
    }
    
    # chr 2
    genMap_2 <- genMap_change[[2]] #%>% as.data.frame()
    # get left distal positions
    l_distal_2 <- genMap_2[genMap_2 < 0.42] 
    # get pericentromere positions
    peri_2 <- genMap_2[genMap_2 > 0.42 & genMap_2 < 0.62] 
    # multiply by DV_scale 
    peri_2 <- peri_2*DV_scale
    # get right distal positions
    r_distal_2 <- genMap_2[genMap_2 > 0.62] 
    # adjust for difference between peri_2[length(peri_2)] and r_distal_2[2]
    r_distal_2 <- r_distal_2 + peri_2[length(peri_2)]  + peri_2[1]
    # combine
    if(length(l_distal_2) == 0 | length(peri_2) == 0 | length(r_distal_2) == 0){
        genMap_2 <- genMap_2
    } else{
        genMap_2 <- c(l_distal_2, peri_2, r_distal_2)
    }
    
    # chr 3
    genMap_3 <- genMap_change[[3]] #%>% as.data.frame()
    # get left distal positions
    l_distal_3 <- genMap_3[genMap_3 < 0.35] 
    # get pericentromere positions
    peri_3 <- genMap_3[genMap_3 > 0.35 & genMap_3 < 0.55] 
    # multiply by DV_scale 
    peri_3 <- peri_3*DV_scale
    # get right distal positions
    r_distal_3 <- genMap_3[genMap_3 > 0.55] 
    # adjust for difference between peri_3[length(peri_3)] and r_distal_3[3]
    r_distal_3 <- r_distal_3 + peri_3[length(peri_3)] + peri_3[1]
    # combine
    if(length(l_distal_3) == 0 | length(peri_3) == 0 | length(r_distal_3) == 0){
        genMap_3 <- genMap_3
    } else{
        genMap_3 <- c(l_distal_3, peri_3, r_distal_3)
    }
    
    # chr 4
    genMap_4 <- genMap_change[[4]] #%>% as.data.frame()
    # get left distal positions
    l_distal_4 <- genMap_4[genMap_4 < 0.57] 
    # get pericentromere positions
    peri_4 <- genMap_4[genMap_4 > 0.57 & genMap_4 < 0.77] 
    # multiply by DV_scale 
    peri_4 <- peri_4*DV_scale
    # get right distal positions
    r_distal_4 <- genMap_4[genMap_4 > 0.77] 
    # adjust for difference between peri_4[length(peri_4)] and r_distal_4[4]
    r_distal_4 <- r_distal_4 + peri_4[length(peri_4)]+ peri_4[1]
    # combine
    if(length(l_distal_4) == 0 | length(peri_4) == 0 | length(r_distal_4) == 0){
        genMap_4 <- genMap_4
    } else{
        genMap_4 <- c(l_distal_4, peri_4, r_distal_4)
    }
    
    # chr 5
    genMap_5 <- genMap_change[[5]] #%>% as.data.frame()
    # get left distal positions
    l_distal_5 <- genMap_5[genMap_5 < 0.53] 
    # get pericentromere positions
    peri_5 <- genMap_5[genMap_5 > 0.53 & genMap_5 < 0.73] 
    # multiply by DV_scale 
    peri_5 <- peri_5*DV_scale
    # get right distal positions
    r_distal_5 <- genMap_5[genMap_5 > 0.73] 
    # adjust for difference between peri_5[length(peri_5)] and r_distal_5[5]
    r_distal_5 <- r_distal_5 + peri_5[length(peri_5)] + peri_5[1]
    # combine
    if(length(l_distal_5) == 0 | length(peri_5) == 0 | length(r_distal_5) == 0){
        genMap_5 <- genMap_5
    } else{
        genMap_5 <- c(l_distal_5, peri_5, r_distal_5)
    }
    
    # chr 6
    genMap_6 <- genMap_change[[6]] #%>% as.data.frame()
    # get left distal positions
    l_distal_6 <- genMap_6[genMap_6 < 0.55] 
    # get pericentromere positions
    peri_6 <- genMap_6[genMap_6 > 0.55 & genMap_6 < 0.75] 
    # multiply by DV_scale 
    peri_6 <- peri_6*DV_scale
    # get right distal positions
    r_distal_6 <- genMap_6[genMap_6 > 0.75] 
    # adjust for difference between peri_6[length(peri_6)] and r_distal_6[6]
    r_distal_6 <- r_distal_6 + peri_6[length(peri_6)] + peri_6[1]
    # combine
    if(length(l_distal_6) == 0 | length(peri_6) == 0 | length(r_distal_6) == 0){
        genMap_6 <- genMap_6
    } else{
        genMap_6 <- c(l_distal_6, peri_6, r_distal_6)
    }
    
    # chr 7
    genMap_7 <- genMap_change[[7]] #%>% as.data.frame()
    # get left distal positions
    l_distal_7 <- genMap_7[genMap_7 < 0.39] 
    # get pericentromere positions
    peri_7 <- genMap_7[genMap_7 > 0.39 & genMap_7 < 0.59] 
    # multiply by DV_scale 
    peri_7 <- peri_7*DV_scale
    # get right distal positions
    r_distal_7 <- genMap_7[genMap_7 > 0.59] 
    # adjust for difference between peri_7[length(peri_7)] and r_distal_7[7]
    r_distal_7 <- r_distal_7 + peri_7[length(peri_7)] + peri_7[1]
    # combine
    if(length(l_distal_7) == 0 | length(peri_7) == 0 | length(r_distal_7) == 0){
        genMap_7 <- genMap_7
    } else{
        genMap_7 <- c(l_distal_7, peri_7, r_distal_7)
    }
    
    # chr 8
    genMap_8 <- genMap_change[[8]] #%>% as.data.frame()
    # get left distal positions
    l_distal_8 <- genMap_8[genMap_8 < 0.42] 
    # get pericentromere positions
    peri_8 <- genMap_8[genMap_8 > 0.42 & genMap_8 < 0.62] 
    # multiply by DV_scale 
    peri_8 <- peri_8*DV_scale
    # get right distal positions
    r_distal_8 <- genMap_8[genMap_8 > 0.62] 
    # adjust for difference between peri_8[length(peri_8)] and r_distal_8[8]
    r_distal_8 <- r_distal_8 + peri_8[length(peri_8)]+ peri_8[1]
    # combine
    if(length(l_distal_8) == 0 | length(peri_8) == 0 | length(r_distal_8) == 0){
        genMap_8 <- genMap_8
    } else{
        genMap_8 <- c(l_distal_8, peri_8, r_distal_8)
    }
    
    # chr 9
    if(BP == TRUE){
        genMap_9 <- genMap_change[[9]] #%>% as.data.frame()
        # get left distal positions
        l_distal_9 <- genMap_9[genMap_9 < 0.8] 
        # get pericentromere positions
        peri_9 <- genMap_9[genMap_9 > 0.8 & genMap_9 < 1] 
        # multiply by DV_scale 
        peri_9 <- peri_9*DV_scale
        # get right distal positions
        r_distal_9 <- genMap_9[genMap_9 > 1] 
        # adjust for difference between peri_9[length(peri_9)] and r_distal_9[9]
        r_distal_9 <- r_distal_9 + peri_9[length(peri_9)]+ peri_9[1]
        # combine
        if(length(l_distal_9) == 0 | length(peri_9) == 0 | length(r_distal_9) == 0){
            genMap_9 <- genMap_9
        } else{
            genMap_9 <- c(l_distal_9, peri_9, r_distal_9)
        }
    }else{
        genMap_9 <- genMap_change[[9]] #%>% as.data.frame()
        # get left distal positions
        l_distal_9 <- genMap_9[genMap_9 < 0.35] 
        # get pericentromere positions
        peri_9 <- genMap_9[genMap_9 > 0.35 & genMap_9 < 0.55] 
        # multiply by DV_scale 
        peri_9 <- peri_9*DV_scale
        # get right distal positions
        r_distal_9 <- genMap_9[genMap_9 > 0.55] 
        # adjust for difference between peri_9[length(peri_9)] and r_distal_9[9]
        r_distal_9 <- r_distal_9 + peri_9[length(peri_9)]+ peri_9[1]
        # combine
        if(length(l_distal_9) == 0 | length(peri_9) == 0 | length(r_distal_9) == 0){
            genMap_9 <- genMap_9
        } else{
            genMap_9 <- c(l_distal_9, peri_9, r_distal_9)
        }
    }
 
    # chr 10
    genMap_10 <- genMap_change[[10]] #%>% as.data.frame()
    # get left distal positions
    l_distal_10 <- genMap_10[genMap_10 < 1] 
    # get pericentromere positions
    peri_10 <- genMap_10[genMap_10 > 1 & genMap_10 < 1.2] 
    # multiply by DV_scale 
    peri_10 <- peri_10*DV_scale
    # get right distal positions
    r_distal_10 <- genMap_10[genMap_10 > 1.2] 
    # adjust for difference between peri_10[length(peri_10)] and r_distal_10[10]
    r_distal_10 <- r_distal_10 + peri_10[length(peri_10)] + peri_10[1]
    # combine
    if(length(l_distal_10) == 0 | length(peri_10) == 0 | length(r_distal_10) == 0){
        genMap_10 <- genMap_10
    } else{
        genMap_10 <- c(l_distal_10, peri_10, r_distal_10)
    }
    
    # chr 11
    genMap_11 <- genMap_change[[11]] #%>% as.data.frame()
    # get left distal positions
    l_distal_11 <- genMap_11[genMap_11 < 0.43] 
    # get pericentromere positions
    peri_11 <- genMap_11[genMap_11 > 0.43 & genMap_11 < 0.63] 
    # multiply by DV_scale 
    peri_11 <- peri_11*DV_scale
    # get right distal positions
    r_distal_11 <- genMap_11[genMap_11 > 0.63] 
    # adjust for difference between peri_11[length(peri_11)] and r_distal_11[11]
    r_distal_11 <- r_distal_11 + peri_11[length(peri_11)] + peri_11[1]
    # combine
    if(length(l_distal_11) == 0 | length(peri_11) == 0 | length(r_distal_11) == 0){
        genMap_11 <- genMap_11
    } else{
        genMap_11 <- c(l_distal_11, peri_11, r_distal_11)
    }
    
    # chr 12
    if(BP == TRUE){
        genMap_12 <- genMap_change[[12]] #%>% as.data.frame()
        # get left distal positions
        l_distal_12 <- genMap_12[genMap_12 < 0.1] 
        # get pericentromere positions
        peri_12 <- genMap_12[genMap_12 > 0.1 & genMap_12 < 0.3] 
        # multiply by DV_scale 
        peri_12 <- peri_12*DV_scale
        # get right distal positions
        r_distal_12 <- genMap_12[genMap_12 > 0.3] 
        # adjust for difference between peri_12[length(peri_12)] and r_distal_12[12]
        r_distal_12 <- r_distal_12 + peri_12[length(peri_12)]+ peri_12[1]
        # combine
        if(length(l_distal_12) == 0 | length(peri_12) == 0 | length(r_distal_12) == 0){
            genMap_12 <- genMap_12
        } else{
            genMap_12 <- c(l_distal_12, peri_12, r_distal_12)
        }
    } else{
        genMap_12 <- genMap_change[[12]] #%>% as.data.frame()
        # get left distal positions
        l_distal_12 <- genMap_12[genMap_12 < 0.45] 
        # get pericentromere positions
        peri_12 <- genMap_12[genMap_12 > 0.45 & genMap_12 < 0.65] 
        # multiply by DV_scale 
        peri_12 <- peri_12*DV_scale
        # get right distal positions
        r_distal_12 <- genMap_12[genMap_12 > 0.65] 
        # adjust for difference between peri_12[length(peri_12)] and r_distal_12[12]
        r_distal_12 <- r_distal_12 + peri_12[length(peri_12)]+ peri_12[1]
        # combine
        if(length(l_distal_12) == 0 | length(peri_12) == 0 | length(r_distal_12) == 0){
            genMap_12 <- genMap_12
        } else{
            genMap_12 <- c(l_distal_12, peri_12, r_distal_12)
        }
    }
   
    
    # chr 13
    genMap_13 <- genMap_change[[13]] #%>% as.data.frame()
    # get left distal positions
    l_distal_13 <- genMap_13[genMap_13 < 0.44] 
    # get pericentromere positions
    peri_13 <- genMap_13[genMap_13 > 0.44 & genMap_13 < 0.64] 
    # multiply by DV_scale 
    peri_13 <- peri_13*DV_scale
    # get right distal positions
    r_distal_13 <- genMap_13[genMap_13 > 0.64] 
    # adjust for difference between peri_13[length(peri_13)] and r_distal_13[13]
    r_distal_13 <- r_distal_13 + peri_13[length(peri_13)]+ peri_13[1]
    # combine
    if(length(l_distal_13) == 0 | length(peri_13) == 0 | length(r_distal_13) == 0){
        genMap_13 <- genMap_13
    } else{
        genMap_13 <- c(l_distal_13, peri_13, r_distal_13)
    }
    
    # chr 14
    genMap_14 <- genMap_change[[14]] #%>% as.data.frame()
    # get left distal positions
    l_distal_14 <- genMap_14[genMap_14 < 0.33] 
    # get pericentromere positions
    peri_14 <- genMap_14[genMap_14 > 0.33 & genMap_14 < 0.53] 
    # multiply by DV_scale 
    peri_14 <- peri_14*DV_scale
    # get right distal positions
    r_distal_14 <- genMap_14[genMap_14 > 0.53] 
    # adjust for difference between peri_14[length(peri_14)] and r_distal_14[14]
    r_distal_14 <- r_distal_14 + peri_14[length(peri_14)]+ peri_14[1]
    # combine
    if(length(l_distal_14) == 0 | length(peri_14) == 0 | length(r_distal_14) == 0){
        genMap_14 <- genMap_14
    } else{
        genMap_14 <- c(l_distal_14, peri_14, r_distal_14)
    }
    
    # chr 15
    if(BP == TRUE){
        genMap_15 <- genMap_change[[15]] #%>% as.data.frame()
        # get left distal positions
        l_distal_15 <- genMap_15[genMap_15 < 1.1] 
        # get pericentromere positions
        peri_15 <- genMap_15[genMap_15 > 1.1 & genMap_15 < 1.3] 
        # multiply by DV_scale 
        peri_15 <- peri_15*DV_scale
        # get right distal positions
        r_distal_15 <- genMap_15[genMap_15 > 1.3] 
        # adjust for difference between peri_15[length(peri_15)] and r_distal_15[15]
        r_distal_15 <- r_distal_15 + peri_15[length(peri_15)]+ peri_15[1]
        # combine
        if(length(l_distal_15) == 0 | length(peri_15) == 0 | length(r_distal_15) == 0){
            genMap_15 <- genMap_15
        } else{
            genMap_15 <- c(l_distal_15, peri_15, r_distal_15)
        }
    }else{
        genMap_15 <- genMap_change[[15]] #%>% as.data.frame()
        # get left distal positions
        l_distal_15 <- genMap_15[genMap_15 < 0.48] 
        # get pericentromere positions
        peri_15 <- genMap_15[genMap_15 > 0.48 & genMap_15 < 0.68] 
        # multiply by DV_scale 
        peri_15 <- peri_15*DV_scale
        # get right distal positions
        r_distal_15 <- genMap_15[genMap_15 > 0.68] 
        # adjust for difference between peri_15[length(peri_15)] and r_distal_15[15]
        r_distal_15 <- r_distal_15 + peri_15[length(peri_15)]+ peri_15[1]
        # combine
        if(length(l_distal_15) == 0 | length(peri_15) == 0 | length(r_distal_15) == 0){
            genMap_15 <- genMap_15
        } else{
            genMap_15 <- c(l_distal_15, peri_15, r_distal_15)
        }
    }
   
    
    # chr 16
    genMap_16 <- genMap_change[[16]] #%>% as.data.frame()
    # get left distal positions
    l_distal_16 <- genMap_16[genMap_16 < 0.42] 
    # get pericentromere positions
    peri_16 <- genMap_16[genMap_16 > 0.42 & genMap_16 < 0.62] 
    # multiply by DV_scale 
    peri_16 <- peri_16*DV_scale
    # get right distal positions
    r_distal_16 <- genMap_16[genMap_16 > 0.62] 
    # adjust for difference between peri_16[length(peri_16)] and r_distal_16[16]
    r_distal_16 <- r_distal_16 + peri_16[length(peri_16)]+ peri_16[1]
    # combine
    if(length(l_distal_16) == 0 | length(peri_16) == 0 | length(r_distal_16) == 0){
        genMap_16 <- genMap_16
    } else{
        genMap_16 <- c(l_distal_16, peri_16, r_distal_16)
    }
    
    # chr 17
    genMap_17 <- genMap_change[[17]] #%>% as.data.frame()
    # get left distal positions
    l_distal_17 <- genMap_17[genMap_17 < 0.39] 
    # get pericentromere positions
    peri_17 <- genMap_17[genMap_17 > 0.39 & genMap_17 < 0.59] 
    # multiply by DV_scale 
    peri_17 <- peri_17*DV_scale
    # get right distal positions
    r_distal_17 <- genMap_17[genMap_17 > 0.59] 
    # adjust for difference between peri_17[length(peri_17)] and r_distal_17[17]
    r_distal_17 <- r_distal_17 + peri_17[length(peri_17)]+ peri_17[1]
    # combine
    if(length(l_distal_17) == 0 | length(peri_17) == 0 | length(r_distal_17) == 0){
        genMap_17 <- genMap_17
    } else{
        genMap_17 <- c(l_distal_17, peri_17, r_distal_17)
    }
    
    # chr 18
    genMap_18 <- genMap_change[[18]] #%>% as.data.frame()
    # get left distal positions
    l_distal_18 <- genMap_18[genMap_18 < 0.4] 
    # get pericentromere positions
    peri_18 <- genMap_18[genMap_18 > 0.4 & genMap_18 < 0.6] 
    # multiply by DV_scale 
    peri_18 <- peri_18*DV_scale
    # get right distal positions
    r_distal_18 <- genMap_18[genMap_18 > 0.6] 
    # adjust for difference between peri_18[length(peri_18)] and r_distal_18[18]
    r_distal_18 <- r_distal_18 + peri_18[length(peri_18)]+ peri_18[1]
    # combine
    if(length(l_distal_18) == 0 | length(peri_18) == 0 | length(r_distal_18) == 0){
        genMap_18 <- genMap_18
    } else{
        genMap_18 <- c(l_distal_18, peri_18, r_distal_18)
    }
    
    # chr 19
    genMap_19 <- genMap_change[[19]] #%>% as.data.frame()
    # get left distal positions
    l_distal_19 <- genMap_19[genMap_19 < 0.55] 
    # get pericentromere positions
    peri_19 <- genMap_19[genMap_19 > 0.55 & genMap_19 < 0.75] 
    # multiply by DV_scale 
    peri_19 <- peri_19*DV_scale
    # get right distal positions
    r_distal_19 <- genMap_19[genMap_19 > 0.75] 
    # adjust for difference between peri_19[length(peri_19)] and r_distal_19[19]
    r_distal_19 <- r_distal_19 + peri_19[length(peri_19)]+ peri_19[1]
    # combine
    if(length(l_distal_19) == 0 | length(peri_19) == 0 | length(r_distal_19) == 0){
        genMap_19 <- genMap_19
    } else{
        genMap_19 <- c(l_distal_19, peri_19, r_distal_19)
    }
    
    # chr 20
    genMap_20 <- genMap_change[[20]] #%>% as.data.frame()
    # get left distal positions
    l_distal_20 <- genMap_20[genMap_20 < 0.45] 
    # get pericentromere positions
    peri_20 <- genMap_20[genMap_20 > 0.45 & genMap_20 < 0.65] 
    # multiply by DV_scale 
    peri_20 <- peri_20*DV_scale
    # get right distal positions
    r_distal_20 <- genMap_20[genMap_20 > 0.65] 
    # adjust for difference between peri_20[length(peri_20)] and r_distal_20[20]
    r_distal_20 <- r_distal_20 + peri_20[length(peri_20)]+ peri_20[1]
    # combine
    if(length(l_distal_20) == 0 | length(peri_20) == 0 | length(r_distal_20) == 0){
        genMap_20 <- genMap_20
    } else{
        genMap_20 <- c(l_distal_20, peri_20, r_distal_20)
    }
    
    # chr 21
    genMap_21 <- genMap_change[[21]] #%>% as.data.frame()
    # get left distal positions
    l_distal_21 <- genMap_21[genMap_21 < 0.74] 
    # get pericentromere positions
    peri_21 <- genMap_21[genMap_21 > 0.74 & genMap_21 < 0.94] 
    # multiply by DV_scale 
    peri_21 <- peri_21*DV_scale
    # get right distal positions
    r_distal_21 <- genMap_21[genMap_21 > 0.94] 
    # adjust for difference between peri_21[length(peri_21)] and r_distal_21[21]
    r_distal_21 <- r_distal_21 + peri_21[length(peri_21)] + peri_21[1]
    # combine
    if(length(l_distal_21) == 0 | length(peri_21) == 0 | length(r_distal_21) == 0){
        genMap_21 <- genMap_21
    } else{
        genMap_21 <- c(l_distal_21, peri_21, r_distal_21)
    }
    
    genMap_change_output <- list(genMap_1, genMap_2, genMap_3, genMap_4, genMap_5, genMap_6, genMap_7, genMap_8, genMap_9, genMap_10, genMap_11, 
                          genMap_12, genMap_13, genMap_14, genMap_15, genMap_16, genMap_17, genMap_18, genMap_19, genMap_20, genMap_21)
    genMap_change_output <<- genMap_change_output
}

