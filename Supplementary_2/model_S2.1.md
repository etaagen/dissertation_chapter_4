Controlled recombination simulation, model analysis, cycle 6
================

### Notes:

Each simulation rep draws from 26 founders, 10 cycles of phenotypic
selection burnin 400 F1s -&gt; 400 DHs, phenotype, serves as TP and GS
cycle 1 starting pop, variance standardized to 1 -&gt; change
recombination distribution (WT, 20X Peri or 20X Chr) –&gt; 400 F1s -&gt;
400 DHs, train model, phenotype and update TP with top 160 (40%) DHs and
drop bottom 80 (20%), advance top 20 DHs (5% selection intensity) for
next cycle of 400 new F1s, repeat for 10 cycles

**WT map type does not get a recombination change, WT is always WT while
Peri and Chr are always 2 or 20X**

<details>
<summary>
Click here: simulation parameters
</summary>

#### Variables:

**Recombination:** WT, 2X or 20X  
**Genetic map change:** WT, Pericentromere or Chromosome  
**Repulsion:** 1:5  
**QTL type:** Random or deleterious variant **Heritability:** 0.2 or
0.8  
**QTL per Chr:** 2 or 200  
**Relationship matrix:** genomewide or causal variant

#### Load packages and data

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

# anova and summary plot functions
source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/summary_stats_anova.R")
source("/Users/ellietaagen/Desktop/github/cr_simulation/Results/table_view_function.R")

gg <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_gg.csv") %>% as.data.frame()
gv <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_gv.csv") %>% as.data.frame()
pa <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_pa.csv") %>% as.data.frame()
be <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_be.csv") %>% as.data.frame()
# be `value` is 2-fold higher than it should be (ASR does not know how to treat VarA/GenicVarA for DH)
be$value <- be$value/2
qtl <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_qtl.csv") %>% as.data.frame()
qtl_af <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_qtl_af.csv") %>% as.data.frame()
```

**Data frames**

gg: population’s genetic gain  
gv: additive genetic variance of population  
pa: prediction accuracy of genomic selection  
be: Bulmer effect (varA/genicVarA)  
qtl: positive and negative effect QTL fixation ratio  
qtl\_af: change in QTL negative allele frequency, subset by small (Q1),
medium, or large (Q3) effect size

**Columns**

**rep** is 1:number of reps simulation was run  
**cycle** designates founder / burnin (0), and GS cycle 1:10  
**value** is the raw response variable measurement, (or average after
summary function, grouped by legend and cycle, of all rep)  
**Matrix** is the relationship matrix used in RRBLUP, genomewide or
causal variant  
**Pop** is the full founder set (F) or biparental (bp)  
**Recombination** is the scale of map change, 2 or 20X  
**H2** is the broad sense heritability  
**QTL** is the number of QTL per chromosome  
**QTL\_type** is R or DV, random or deleterious variant  
**Map\_Type** is the WT, Pericentromere, or Chromosome-wide change to
the genetic map, given a Recombination scale  
**Repulsion** is 1:5 representing different coupling and replusion
ratios

-   1: Additive effect signs are positive for all QTL (select against
    minor allele)

-   2: Random 2/3 of additive effect signs are positive and 1/3 are
    negative for QTL

-   3: Random 1/2 of additive effect signs are positive and 1/2 are
    negative for QTL

-   4: 1/2 of additive effect signs are positive and 1/2 are negative
    for QTL, alternating positive or negative each QTL

-   5: Random 1/3 of additive effect signs are positive and 2/3 are
    negative for QTL (most selection for minor allele)

</details>

### Genetic gain

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 6
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  101.0194  2  < 2.2e-16 ***
    ## Recombination             110.5811  1  < 2.2e-16 ***
    ## QTL                     13437.2069  1  < 2.2e-16 ***
    ## H2                      30885.4319  1  < 2.2e-16 ***
    ## Repulsion               44368.8764  4  < 2.2e-16 ***
    ## Matrix                   8913.4296  1  < 2.2e-16 ***
    ## QTL_type                  228.6894  1  < 2.2e-16 ***
    ## Map_type:Recombination     53.8515  2  2.024e-12 ***
    ## Map_type:QTL               49.8383  2  1.506e-11 ***
    ## Map_type:H2                21.1028  2  2.616e-05 ***
    ## Map_type:Repulsion         21.9181  8    0.00507 ** 
    ## Map_type:Matrix            98.8036  2  < 2.2e-16 ***
    ## Map_type:QTL_type          21.6756  2  1.964e-05 ***
    ## Recombination:QTL          16.4759  1  4.927e-05 ***
    ## Recombination:H2            2.7113  1    0.09964 .  
    ## Recombination:Repulsion    29.0321  4  7.701e-06 ***
    ## Recombination:Matrix       56.5239  1  5.552e-14 ***
    ## Recombination:QTL_type      0.4699  1    0.49304    
    ## QTL:H2                  13374.6500  1  < 2.2e-16 ***
    ## QTL:Repulsion           20740.1563  4  < 2.2e-16 ***
    ## QTL:Matrix                303.3262  1  < 2.2e-16 ***
    ## QTL:QTL_type                0.5557  1    0.45598    
    ## H2:Repulsion              961.6636  4  < 2.2e-16 ***
    ## H2:Matrix                  64.7833  1  8.360e-16 ***
    ## H2:QTL_type                28.5104  1  9.320e-08 ***
    ## Repulsion:Matrix         1529.8692  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        144.6101  4  < 2.2e-16 ***
    ## Matrix:QTL_type            19.1637  1  1.200e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 6, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 6, WT map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

#### Cycle 6, Pericentromere map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

#### Cycle 6, Chromosome map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

</details>

### Genetic variance

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 6
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  790.8788  2  < 2.2e-16 ***
    ## Recombination             413.1177  1  < 2.2e-16 ***
    ## QTL                     30294.9966  1  < 2.2e-16 ***
    ## H2                       7559.8438  1  < 2.2e-16 ***
    ## Repulsion                5698.5540  4  < 2.2e-16 ***
    ## Matrix                   6754.3559  1  < 2.2e-16 ***
    ## QTL_type                  426.2978  1  < 2.2e-16 ***
    ## Map_type:Recombination    170.4159  2  < 2.2e-16 ***
    ## Map_type:QTL              288.3043  2  < 2.2e-16 ***
    ## Map_type:H2                15.2764  2  0.0004817 ***
    ## Map_type:Repulsion        146.8844  8  < 2.2e-16 ***
    ## Map_type:Matrix            88.4942  2  < 2.2e-16 ***
    ## Map_type:QTL_type          37.3038  2  7.936e-09 ***
    ## Recombination:QTL         138.6647  1  < 2.2e-16 ***
    ## Recombination:H2            0.3073  1  0.5793625    
    ## Recombination:Repulsion    92.3304  4  < 2.2e-16 ***
    ## Recombination:Matrix       27.4763  1  1.590e-07 ***
    ## Recombination:QTL_type      8.2971  1  0.0039710 ** 
    ## QTL:H2                     44.8762  1  2.099e-11 ***
    ## QTL:Repulsion            4063.3353  4  < 2.2e-16 ***
    ## QTL:Matrix                850.4611  1  < 2.2e-16 ***
    ## QTL:QTL_type              285.4757  1  < 2.2e-16 ***
    ## H2:Repulsion              669.3071  4  < 2.2e-16 ***
    ## H2:Matrix                 908.8516  1  < 2.2e-16 ***
    ## H2:QTL_type                22.2243  1  2.426e-06 ***
    ## Repulsion:Matrix          656.0847  4  < 2.2e-16 ***
    ## Repulsion:QTL_type          7.2796  4  0.1218303    
    ## Matrix:QTL_type            33.4535  1  7.299e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Click to expand: Cycle 6, best and worst simulation
parameter prediction summaries
</summary>

#### Cycle 6, WT map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

#### Cycle 6, Pericentromere map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

#### Cycle 6, Chromosome map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

</details>

### Bulmer effect

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 6
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  238.4673  2  < 2.2e-16 ***
    ## Recombination             119.9017  1  < 2.2e-16 ***
    ## QTL                         5.6399  1  0.0175564 *  
    ## H2                      14321.9523  1  < 2.2e-16 ***
    ## Repulsion               49687.4998  4  < 2.2e-16 ***
    ## Matrix                   1345.8369  1  < 2.2e-16 ***
    ## QTL_type                   10.9079  1  0.0009576 ***
    ## Map_type:Recombination     58.3355  2  2.151e-13 ***
    ## Map_type:QTL               76.6621  2  < 2.2e-16 ***
    ## Map_type:H2                62.0610  2  3.339e-14 ***
    ## Map_type:Repulsion         24.0929  8  0.0022110 ** 
    ## Map_type:Matrix            18.4938  2  9.641e-05 ***
    ## Map_type:QTL_type          18.9998  2  7.486e-05 ***
    ## Recombination:QTL          35.2666  1  2.875e-09 ***
    ## Recombination:H2           26.7614  1  2.302e-07 ***
    ## Recombination:Repulsion    22.9195  4  0.0001314 ***
    ## Recombination:Matrix        5.5149  1  0.0188547 *  
    ## Recombination:QTL_type      5.3654  1  0.0205405 *  
    ## QTL:H2                   2875.9284  1  < 2.2e-16 ***
    ## QTL:Repulsion           18417.7635  4  < 2.2e-16 ***
    ## QTL:Matrix                196.9912  1  < 2.2e-16 ***
    ## QTL:QTL_type               58.6410  1  1.892e-14 ***
    ## H2:Repulsion             8370.6627  4  < 2.2e-16 ***
    ## H2:Matrix                   0.1766  1  0.6743049    
    ## H2:QTL_type               149.4851  1  < 2.2e-16 ***
    ## Repulsion:Matrix          933.5257  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        881.8917  4  < 2.2e-16 ***
    ## Matrix:QTL_type           117.7713  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 6, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 6, WT map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->

#### Cycle 6, Pericentromere map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

#### Cycle 6, Chromosome map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

</details>

### Prediction accuracy

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 6
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  134.9856  2  < 2.2e-16 ***
    ## Recombination              60.4413  1  7.581e-15 ***
    ## QTL                      2841.1749  1  < 2.2e-16 ***
    ## H2                      15828.2522  1  < 2.2e-16 ***
    ## Repulsion                1757.6344  4  < 2.2e-16 ***
    ## Matrix                   2680.1305  1  < 2.2e-16 ***
    ## QTL_type                   23.7940  1  1.072e-06 ***
    ## Map_type:Recombination     52.6813  2  3.634e-12 ***
    ## Map_type:QTL               81.3879  2  < 2.2e-16 ***
    ## Map_type:H2                 0.8367  2   0.658130    
    ## Map_type:Repulsion         11.3497  8   0.182661    
    ## Map_type:Matrix            71.7847  2  2.583e-16 ***
    ## Map_type:QTL_type           4.9238  2   0.085271 .  
    ## Recombination:QTL          14.0391  1   0.000179 ***
    ## Recombination:H2            1.7826  1   0.181830    
    ## Recombination:Repulsion     8.3627  4   0.079161 .  
    ## Recombination:Matrix       24.1059  1  9.118e-07 ***
    ## Recombination:QTL_type      1.8431  1   0.174588    
    ## QTL:H2                   1960.3430  1  < 2.2e-16 ***
    ## QTL:Repulsion             339.6510  4  < 2.2e-16 ***
    ## QTL:Matrix                716.6802  1  < 2.2e-16 ***
    ## QTL:QTL_type                3.2833  1   0.069989 .  
    ## H2:Repulsion               38.6133  4  8.373e-08 ***
    ## H2:Matrix                 547.3578  1  < 2.2e-16 ***
    ## H2:QTL_type                 1.3411  1   0.246837    
    ## Repulsion:Matrix           54.4870  4  4.161e-11 ***
    ## Repulsion:QTL_type          6.1297  4   0.189672    
    ## Matrix:QTL_type            43.9344  1  3.396e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 6, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 6, WT map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-18-4.png)<!-- -->

#### Cycle 6, Pericentromere map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->

#### Cycle 6, Chromosome map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

</details>

### QTL fixation

**Linear model has `singular fit`, omitted**

### QTL positive allele frequency change

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA table
</summary>

**Note:** total cycle change observations. See .Rmd file for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type + allele)^2 + (1|rep)`

ANOVA table:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                    7.5916  2  0.0224654 *  
    ## Recombination              13.8127  1  0.0002020 ***
    ## QTL                     72156.0864  1  < 2.2e-16 ***
    ## H2                          4.8586  1  0.0275086 *  
    ## Repulsion               31640.6715  4  < 2.2e-16 ***
    ## Matrix                   2645.3858  1  < 2.2e-16 ***
    ## QTL_type                  366.9488  1  < 2.2e-16 ***
    ## Allele                   4355.9982  2  < 2.2e-16 ***
    ## Map_type:Recombination      6.6114  2  0.0366740 *  
    ## Map_type:QTL                2.5264  2  0.2827471    
    ## Map_type:H2                13.8616  2  0.0009772 ***
    ## Map_type:Repulsion          2.6898  8  0.9522909    
    ## Map_type:Matrix            19.9313  2  4.699e-05 ***
    ## Map_type:QTL_type           1.7159  2  0.4240226    
    ## Map_type:Allele             2.1492  4  0.7083346    
    ## Recombination:QTL           2.3565  1  0.1247642    
    ## Recombination:H2            1.6002  1  0.2058765    
    ## Recombination:Repulsion     4.3689  4  0.3583720    
    ## Recombination:Matrix        9.9975  1  0.0015675 ** 
    ## Recombination:QTL_type      0.4094  1  0.5222764    
    ## Recombination:Allele        5.4240  2  0.0664053 .  
    ## QTL:H2                   1189.9138  1  < 2.2e-16 ***
    ## QTL:Repulsion           25751.5391  4  < 2.2e-16 ***
    ## QTL:Matrix               1475.7672  1  < 2.2e-16 ***
    ## QTL:QTL_type               88.0456  1  < 2.2e-16 ***
    ## QTL:Allele              20841.6066  2  < 2.2e-16 ***
    ## H2:Repulsion               27.3418  4  1.695e-05 ***
    ## H2:Matrix                  48.6184  1  3.109e-12 ***
    ## H2:QTL_type                23.2883  1  1.394e-06 ***
    ## H2:Allele               15462.7855  2  < 2.2e-16 ***
    ## Repulsion:Matrix         1374.6931  4  < 2.2e-16 ***
    ## Repulsion:QTL_type         31.3299  4  2.622e-06 ***
    ## Repulsion:Allele         2253.8078  8  < 2.2e-16 ***
    ## Matrix:QTL_type             9.8533  1  0.0016953 ** 
    ## Matrix:Allele             183.2536  2  < 2.2e-16 ***
    ## QTL_type:Allele            34.6539  2  2.985e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Best and worst simulation parameter prediction summaries
</summary>

#### WT map

<div id="hozddgfprz" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#hozddgfprz .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#hozddgfprz .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#hozddgfprz .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#hozddgfprz .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#hozddgfprz .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hozddgfprz .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#hozddgfprz .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#hozddgfprz .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#hozddgfprz .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#hozddgfprz .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#hozddgfprz .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#hozddgfprz .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#hozddgfprz .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#hozddgfprz .gt_from_md > :first-child {
  margin-top: 0;
}

#hozddgfprz .gt_from_md > :last-child {
  margin-bottom: 0;
}

#hozddgfprz .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#hozddgfprz .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#hozddgfprz .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hozddgfprz .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#hozddgfprz .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hozddgfprz .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#hozddgfprz .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#hozddgfprz .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hozddgfprz .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#hozddgfprz .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#hozddgfprz .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#hozddgfprz .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#hozddgfprz .gt_left {
  text-align: left;
}

#hozddgfprz .gt_center {
  text-align: center;
}

#hozddgfprz .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#hozddgfprz .gt_font_normal {
  font-weight: normal;
}

#hozddgfprz .gt_font_bold {
  font-weight: bold;
}

#hozddgfprz .gt_font_italic {
  font-style: italic;
}

#hozddgfprz .gt_super {
  font-size: 65%;
}

#hozddgfprz .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  <thead class="gt_header">
    <tr>
      <th colspan="8" class="gt_heading gt_title gt_font_normal" style><strong>Positive QTL AF change estimates, top 10</strong></th>
    </tr>
    <tr>
      <th colspan="8" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style><strong>Map type:</strong> WT, <strong>QTL type:</strong> R</th>
    </tr>
  </thead>
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">QTL</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">H2</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Repulsion</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Matrix</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Allele</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">estimate</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">std.error</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p.value</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2750611</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2738180</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">8</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2515346</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">8</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2511277</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">8</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">small_pos</td>
<td class="gt_row gt_right">0.2396553</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">8</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">small_pos</td>
<td class="gt_row gt_right">0.2387044</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">3</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2384267</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">3</td>
<td class="gt_row gt_left">CV</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2379634</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">GW</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2268367</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
    <tr><td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">2</td>
<td class="gt_row gt_left">5</td>
<td class="gt_row gt_left">GW</td>
<td class="gt_row gt_left">med_pos</td>
<td class="gt_row gt_right">0.2235032</td>
<td class="gt_row gt_right">0.001516022</td>
<td class="gt_row gt_right">0</td></tr>
  </tbody>
  
  
</table>
</div>

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->

#### Pericentromere map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-24-4.png)<!-- -->

#### Chromosome map

![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-25-3.png)<!-- -->![](master_results_ASR_models_cycle6_files/figure-gfm/unnamed-chunk-25-4.png)<!-- -->

</details>
