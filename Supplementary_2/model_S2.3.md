Full versus BP pop model analysis, cycle 6 & 10
================

<details>
<summary>
Click here: Load packages, format data
</summary>

``` r
library(tidyverse)
library(data.table)
library(kableExtra)
library(ggsci)
library(gt)
library(lme4)
library(emmeans)
library(car)
library(broom.mixed)

gg <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gg.csv") %>% as.data.frame()
gv <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gv.csv") %>% as.data.frame()
pa <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_pa.csv") %>% as.data.frame()
be <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_be.csv") %>% as.data.frame()
# be `value` is 2-fold higher than it should be (ASR does not know how to treat VarA/GenicVarA for DH)
be$value <- be$value/2
qtl <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_qtl.csv") %>% as.data.frame()
qtl_af <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_qtl_af.csv") %>% as.data.frame()


gg_bp <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/bp_gg.csv") %>% as.data.frame()
gv_bp <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/bp_gv.csv") %>% as.data.frame()
pa_bp <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/bp_pa.csv") %>% as.data.frame()
be_bp <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/bp_be.csv") %>% as.data.frame()
# be `value` is 2-fold higher than it should be (ASR does not know how to treat VarA/GenicVarA for DH)
be_bp$value <- be_bp$value/2
qtl_bp <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/bp_qtl.csv") %>% as.data.frame()
qtl_af_bp <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/bp_qtl_af.csv") %>% as.data.frame()

#filter full pop for repulsion conditions 3 and 4 and QTL type random
gg <- gg %>% filter(Repulsion %in% c(3, 4), QTL_type == "R")
gv <- gv %>% filter(Repulsion %in% c(3, 4), QTL_type == "R")
be <- be %>% filter(Repulsion %in% c(3, 4), QTL_type == "R")
pa <- pa %>% filter(Repulsion %in% c(3, 4), QTL_type == "R")
qtl <- qtl %>% filter(Repulsion %in% c(3, 4), QTL_type == "R")
qtl_af <- qtl_af %>% filter(Repulsion %in% c(3, 4), QTL_type == "R")

# combine
gg <- rbind(gg, gg_bp)
gv <- rbind(gv, gv_bp)
be <- rbind(be, be_bp)
pa <- rbind(pa, pa_bp)
qtl <- rbind(qtl, qtl_bp)
qtl_af <- rbind(qtl_af, qtl_af_bp)
```

</details>

### Genetic gain

**Linear model, ANOVA table**

<details>
<summary>
Click here: Model, ANOVA tables cycle 6 & 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + founder population)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   55.8562  2  7.430e-13 ***
    ## Recombination              33.1889  1  8.363e-09 ***
    ## QTL                      4982.5855  1  < 2.2e-16 ***
    ## H2                      17458.8189  1  < 2.2e-16 ***
    ## Repulsion                 644.6819  1  < 2.2e-16 ***
    ## Matrix                   3044.8654  1  < 2.2e-16 ***
    ## Pop                      1241.5490  1  < 2.2e-16 ***
    ## Map_type:Recombination     21.2793  2  2.395e-05 ***
    ## Map_type:QTL                9.4028  2  0.0090824 ** 
    ## Map_type:H2                17.4501  2  0.0001625 ***
    ## Map_type:Repulsion         12.6661  2  0.0017766 ** 
    ## Map_type:Matrix            42.8628  2  4.926e-10 ***
    ## Map_type:Pop                1.0498  2  0.5916214    
    ## Recombination:QTL           2.3307  1  0.1268484    
    ## Recombination:H2            6.7208  1  0.0095293 ** 
    ## Recombination:Repulsion     6.6860  1  0.0097174 ** 
    ## Recombination:Matrix       24.2934  1  8.272e-07 ***
    ## Recombination:Pop           5.0407  1  0.0247585 *  
    ## QTL:H2                   4774.3538  1  < 2.2e-16 ***
    ## QTL:Repulsion            1567.5107  1  < 2.2e-16 ***
    ## QTL:Matrix                467.7522  1  < 2.2e-16 ***
    ## QTL:Pop                  1148.3618  1  < 2.2e-16 ***
    ## H2:Repulsion                0.2058  1  0.6501056    
    ## H2:Matrix                 119.3445  1  < 2.2e-16 ***
    ## H2:Pop                     31.0053  1  2.573e-08 ***
    ## Repulsion:Matrix          177.2795  1  < 2.2e-16 ***
    ## Repulsion:Pop            1516.8281  1  < 2.2e-16 ***
    ## Matrix:Pop                216.5118  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   27.9042  2  8.723e-07 ***
    ## Recombination               3.4472  1  0.0633575 .  
    ## QTL                     31256.3429  1  < 2.2e-16 ***
    ## H2                      22858.5556  1  < 2.2e-16 ***
    ## Repulsion                 153.1343  1  < 2.2e-16 ***
    ## Matrix                   1448.1796  1  < 2.2e-16 ***
    ## Pop                       320.4048  1  < 2.2e-16 ***
    ## Map_type:Recombination      5.0441  2  0.0802956 .  
    ## Map_type:QTL               19.1750  2  6.858e-05 ***
    ## Map_type:H2                57.0826  2  4.024e-13 ***
    ## Map_type:Repulsion         21.0227  2  2.723e-05 ***
    ## Map_type:Matrix            47.5893  2  4.636e-11 ***
    ## Map_type:Pop               25.8656  2  2.417e-06 ***
    ## Recombination:QTL          11.0752  1  0.0008749 ***
    ## Recombination:H2           24.6609  1  6.836e-07 ***
    ## Recombination:Repulsion     8.6642  1  0.0032453 ** 
    ## Recombination:Matrix       30.6566  1  3.080e-08 ***
    ## Recombination:Pop          19.8990  1  8.164e-06 ***
    ## QTL:H2                  12377.0703  1  < 2.2e-16 ***
    ## QTL:Repulsion            1667.2067  1  < 2.2e-16 ***
    ## QTL:Matrix                  1.9020  1  0.1678523    
    ## QTL:Pop                  2823.3491  1  < 2.2e-16 ***
    ## H2:Repulsion               65.6844  1  5.292e-16 ***
    ## H2:Matrix                   8.6043  1  0.0033537 ** 
    ## H2:Pop                    183.4964  1  < 2.2e-16 ***
    ## Repulsion:Matrix          123.9723  1  < 2.2e-16 ***
    ## Repulsion:Pop             904.7389  1  < 2.2e-16 ***
    ## Matrix:Pop                333.7961  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>

### Genetic variance

**Linear model, ANOVA table**
<details>
<summary>
Click here: Model, ANOVA tables cycle 6 & 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + founder population)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  805.4446  2  < 2.2e-16 ***
    ## Recombination             414.5369  1  < 2.2e-16 ***
    ## QTL                     14324.4222  1  < 2.2e-16 ***
    ## H2                       4880.4469  1  < 2.2e-16 ***
    ## Repulsion                 149.8899  1  < 2.2e-16 ***
    ## Matrix                   2315.5318  1  < 2.2e-16 ***
    ## Pop                       736.7711  1  < 2.2e-16 ***
    ## Map_type:Recombination    171.5885  2  < 2.2e-16 ***
    ## Map_type:QTL              269.0285  2  < 2.2e-16 ***
    ## Map_type:H2                 8.2580  2  0.0160991 *  
    ## Map_type:Repulsion         16.4942  2  0.0002620 ***
    ## Map_type:Matrix            46.3123  2  8.778e-11 ***
    ## Map_type:Pop               24.6846  2  4.363e-06 ***
    ## Recombination:QTL         131.8310  1  < 2.2e-16 ***
    ## Recombination:H2            1.4628  1  0.2264850    
    ## Recombination:Repulsion     1.7987  1  0.1798719    
    ## Recombination:Matrix       41.1642  1  1.400e-10 ***
    ## Recombination:Pop           3.1250  1  0.0770994 .  
    ## QTL:H2                     96.7989  1  < 2.2e-16 ***
    ## QTL:Repulsion             272.3814  1  < 2.2e-16 ***
    ## QTL:Matrix                722.2122  1  < 2.2e-16 ***
    ## QTL:Pop                  1351.0845  1  < 2.2e-16 ***
    ## H2:Repulsion               14.7205  1  0.0001247 ***
    ## H2:Matrix                 288.3233  1  < 2.2e-16 ***
    ## H2:Pop                     15.1690  1  9.830e-05 ***
    ## Repulsion:Matrix            2.4946  1  0.1142410    
    ## Repulsion:Pop               2.6989  1  0.1004205    
    ## Matrix:Pop                 27.0769  1  1.955e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                 1079.9839  2  < 2.2e-16 ***
    ## Recombination             534.3883  1  < 2.2e-16 ***
    ## QTL                     19553.1362  1  < 2.2e-16 ***
    ## H2                       2272.5577  1  < 2.2e-16 ***
    ## Repulsion                 203.1311  1  < 2.2e-16 ***
    ## Matrix                   1066.4256  1  < 2.2e-16 ***
    ## Pop                       544.6201  1  < 2.2e-16 ***
    ## Map_type:Recombination    311.5962  2  < 2.2e-16 ***
    ## Map_type:QTL              666.1265  2  < 2.2e-16 ***
    ## Map_type:H2                35.6621  2  1.803e-08 ***
    ## Map_type:Repulsion         35.2940  2  2.168e-08 ***
    ## Map_type:Matrix            34.7542  2  2.839e-08 ***
    ## Map_type:Pop               23.1357  2  9.466e-06 ***
    ## Recombination:QTL         298.2391  1  < 2.2e-16 ***
    ## Recombination:H2           15.1004  1  0.0001019 ***
    ## Recombination:Repulsion     3.4250  1  0.0642150 .  
    ## Recombination:Matrix        8.0471  1  0.0045575 ** 
    ## Recombination:Pop           3.7995  1  0.0512680 .  
    ## QTL:H2                    196.4485  1  < 2.2e-16 ***
    ## QTL:Repulsion             145.4550  1  < 2.2e-16 ***
    ## QTL:Matrix                 71.1235  1  < 2.2e-16 ***
    ## QTL:Pop                   917.1399  1  < 2.2e-16 ***
    ## H2:Repulsion                5.1189  1  0.0236660 *  
    ## H2:Matrix                 227.6344  1  < 2.2e-16 ***
    ## H2:Pop                      3.9548  1  0.0467386 *  
    ## Repulsion:Matrix            5.5448  1  0.0185360 *  
    ## Repulsion:Pop              14.6126  1  0.0001320 ***
    ## Matrix:Pop                  2.6278  1  0.1050052    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>

### Bulmer effect

**Linear model, ANOVA table**  
<details>
<summary>
Click here: Model, ANOVA tables cycle 6 & 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + founder population)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  551.4361  2  < 2.2e-16 ***
    ## Recombination             218.0605  1  < 2.2e-16 ***
    ## QTL                     17420.4556  1  < 2.2e-16 ***
    ## H2                       4367.8499  1  < 2.2e-16 ***
    ## Repulsion                6193.3212  1  < 2.2e-16 ***
    ## Matrix                     18.8756  1  1.395e-05 ***
    ## Pop                      4776.1203  1  < 2.2e-16 ***
    ## Map_type:Recombination    153.8638  2  < 2.2e-16 ***
    ## Map_type:QTL               56.8488  2  4.523e-13 ***
    ## Map_type:H2                10.4742  2  0.0053156 ** 
    ## Map_type:Repulsion         11.8652  2  0.0026516 ** 
    ## Map_type:Matrix            30.5412  2  2.334e-07 ***
    ## Map_type:Pop               29.2980  2  4.345e-07 ***
    ## Recombination:QTL          63.5269  1  1.582e-15 ***
    ## Recombination:H2            2.6872  1  0.1011567    
    ## Recombination:Repulsion     3.2213  1  0.0726847 .  
    ## Recombination:Matrix       14.5916  1  0.0001335 ***
    ## Recombination:Pop           4.7810  1  0.0287750 *  
    ## QTL:H2                    467.8708  1  < 2.2e-16 ***
    ## QTL:Repulsion            2829.3072  1  < 2.2e-16 ***
    ## QTL:Matrix                247.9898  1  < 2.2e-16 ***
    ## QTL:Pop                    23.5205  1  1.236e-06 ***
    ## H2:Repulsion              103.1661  1  < 2.2e-16 ***
    ## H2:Matrix                  81.2898  1  < 2.2e-16 ***
    ## H2:Pop                     59.7311  1  1.087e-14 ***
    ## Repulsion:Matrix          121.1705  1  < 2.2e-16 ***
    ## Repulsion:Pop             499.0956  1  < 2.2e-16 ***
    ## Matrix:Pop                151.8188  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                             Chisq Df Pr(>Chisq)    
    ## Map_type                 309.8924  2  < 2.2e-16 ***
    ## Recombination             74.4298  1  < 2.2e-16 ***
    ## QTL                     9964.5197  1  < 2.2e-16 ***
    ## H2                      4599.0218  1  < 2.2e-16 ***
    ## Repulsion               4560.3687  1  < 2.2e-16 ***
    ## Matrix                     6.2291  1   0.012567 *  
    ## Pop                     2396.1822  1  < 2.2e-16 ***
    ## Map_type:Recombination    96.3451  2  < 2.2e-16 ***
    ## Map_type:QTL              29.4262  2  4.076e-07 ***
    ## Map_type:H2               13.7999  2   0.001008 ** 
    ## Map_type:Repulsion         0.7119  2   0.700494    
    ## Map_type:Matrix            3.9958  2   0.135622    
    ## Map_type:Pop              12.7568  2   0.001698 ** 
    ## Recombination:QTL         29.0288  1  7.131e-08 ***
    ## Recombination:H2           2.1586  1   0.141772    
    ## Recombination:Repulsion    1.3471  1   0.245778    
    ## Recombination:Matrix       1.7978  1   0.179976    
    ## Recombination:Pop          0.2787  1   0.597557    
    ## QTL:H2                   522.4084  1  < 2.2e-16 ***
    ## QTL:Repulsion            862.6891  1  < 2.2e-16 ***
    ## QTL:Matrix                 6.6834  1   0.009732 ** 
    ## QTL:Pop                   16.7098  1  4.356e-05 ***
    ## H2:Repulsion             152.4720  1  < 2.2e-16 ***
    ## H2:Matrix                 51.4593  1  7.310e-13 ***
    ## H2:Pop                     4.7450  1   0.029383 *  
    ## Repulsion:Matrix           7.3473  1   0.006716 ** 
    ## Repulsion:Pop            337.1305  1  < 2.2e-16 ***
    ## Matrix:Pop                30.9450  1  2.654e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>

### Prediction accuracy

**Linear model, ANOVA table**  
<details>
<summary>
Click here: Model, ANOVA tables cycle 6 & 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + founder population)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   42.2310  2  6.756e-10 ***
    ## Recombination              22.4956  1  2.106e-06 ***
    ## QTL                      1041.5400  1  < 2.2e-16 ***
    ## H2                      10104.1949  1  < 2.2e-16 ***
    ## Repulsion                   5.4588  1  0.0194704 *  
    ## Matrix                   1582.0466  1  < 2.2e-16 ***
    ## Pop                       100.7561  1  < 2.2e-16 ***
    ## Map_type:Recombination      8.5162  2  0.0141493 *  
    ## Map_type:QTL               39.7344  2  2.354e-09 ***
    ## Map_type:H2                 2.1307  2  0.3445987    
    ## Map_type:Repulsion          0.2354  2  0.8889786    
    ## Map_type:Matrix            15.0641  2  0.0005356 ***
    ## Map_type:Pop                3.2927  2  0.1927566    
    ## Recombination:QTL           1.6069  1  0.2049221    
    ## Recombination:H2            0.1626  1  0.6868016    
    ## Recombination:Repulsion     0.6919  1  0.4055206    
    ## Recombination:Matrix        5.1228  1  0.0236137 *  
    ## Recombination:Pop           1.7228  1  0.1893291    
    ## QTL:H2                    665.6454  1  < 2.2e-16 ***
    ## QTL:Repulsion              11.6672  1  0.0006361 ***
    ## QTL:Matrix                861.0899  1  < 2.2e-16 ***
    ## QTL:Pop                   134.9132  1  < 2.2e-16 ***
    ## H2:Repulsion                1.0301  1  0.3101345    
    ## H2:Matrix                 278.8678  1  < 2.2e-16 ***
    ## H2:Pop                     86.6546  1  < 2.2e-16 ***
    ## Repulsion:Matrix            0.6738  1  0.4117391    
    ## Repulsion:Pop              75.3651  1  < 2.2e-16 ***
    ## Matrix:Pop                106.2746  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                    5.9763  2  0.0503817 .  
    ## Recombination               2.4479  1  0.1176812    
    ## QTL                      7900.6082  1  < 2.2e-16 ***
    ## H2                      14306.9502  1  < 2.2e-16 ***
    ## Repulsion                  27.3222  1  1.722e-07 ***
    ## Matrix                    336.4915  1  < 2.2e-16 ***
    ## Pop                        11.7580  1  0.0006058 ***
    ## Map_type:Recombination     17.3796  2  0.0001683 ***
    ## Map_type:QTL               35.5288  2  1.928e-08 ***
    ## Map_type:H2                 0.1595  2  0.9233570    
    ## Map_type:Repulsion          0.1314  2  0.9364028    
    ## Map_type:Matrix            36.8795  2  9.811e-09 ***
    ## Map_type:Pop               12.6268  2  0.0018118 ** 
    ## Recombination:QTL           3.8322  1  0.0502782 .  
    ## Recombination:H2            0.5850  1  0.4443659    
    ## Recombination:Repulsion     1.2114  1  0.2710567    
    ## Recombination:Matrix       21.8053  1  3.018e-06 ***
    ## Recombination:Pop           3.4858  1  0.0618957 .  
    ## QTL:H2                   2726.2398  1  < 2.2e-16 ***
    ## QTL:Repulsion               1.2171  1  0.2699393    
    ## QTL:Matrix                 71.3264  1  < 2.2e-16 ***
    ## QTL:Pop                   152.5548  1  < 2.2e-16 ***
    ## H2:Repulsion               23.0575  1  1.572e-06 ***
    ## H2:Matrix                  47.7022  1  4.961e-12 ***
    ## H2:Pop                      0.3300  1  0.5656678    
    ## Repulsion:Matrix            0.2435  1  0.6217113    
    ## Repulsion:Pop               0.3133  1  0.5756731    
    ## Matrix:Pop                 93.8142  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>

### QTL fixation

**Linear model has `singular fit`, omitted / plots have little
variation**

### QTL allele frequency change

**Linear model, ANOVA table**

<details>
<summary>
Click here: Model, ANOVA table
</summary>

**Note:** total cycle change observations. See .Rmd file for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + founder pop + allele)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')

ANOVA table:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                    2.5062  2   0.285624    
    ## Recombination               0.0006  1   0.980858    
    ## QTL                     96869.9942  1  < 2.2e-16 ***
    ## H2                       2058.8998  1  < 2.2e-16 ***
    ## Repulsion                2628.4210  1  < 2.2e-16 ***
    ## Matrix                   2547.5476  1  < 2.2e-16 ***
    ## Pop                     18758.9025  1  < 2.2e-16 ***
    ## Allele                  17474.5789  2  < 2.2e-16 ***
    ## Map_type:Recombination      0.6637  2   0.717612    
    ## Map_type:QTL                0.5593  2   0.756067    
    ## Map_type:H2                 3.9275  2   0.140334    
    ## Map_type:Repulsion         10.3287  2   0.005717 ** 
    ## Map_type:Matrix             1.6976  2   0.427933    
    ## Map_type:Pop                4.3760  2   0.112138    
    ## Map_type:Allele             2.8114  4   0.589864    
    ## Recombination:QTL           0.3984  1   0.527920    
    ## Recombination:H2            0.2594  1   0.610523    
    ## Recombination:Repulsion     0.2877  1   0.591676    
    ## Recombination:Matrix        2.6531  1   0.103347    
    ## Recombination:Pop           1.7944  1   0.180390    
    ## Recombination:Allele        1.6089  2   0.447331    
    ## QTL:H2                    713.1140  1  < 2.2e-16 ***
    ## QTL:Repulsion             158.1588  1  < 2.2e-16 ***
    ## QTL:Matrix               2242.9123  1  < 2.2e-16 ***
    ## QTL:Pop                 20184.1481  1  < 2.2e-16 ***
    ## QTL:Allele               3283.0215  2  < 2.2e-16 ***
    ## H2:Repulsion                7.9980  1   0.004683 ** 
    ## H2:Matrix                   4.0402  1   0.044428 *  
    ## H2:Pop                   1884.7381  1  < 2.2e-16 ***
    ## H2:Allele                2523.4523  2  < 2.2e-16 ***
    ## Repulsion:Matrix           37.2100  1  1.061e-09 ***
    ## Repulsion:Pop              39.7650  1  2.864e-10 ***
    ## Repulsion:Allele          118.1441  2  < 2.2e-16 ***
    ## Matrix:Pop                376.1254  1  < 2.2e-16 ***
    ## Matrix:Allele             515.1623  2  < 2.2e-16 ***
    ## Pop:Allele              17789.7536  2  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
