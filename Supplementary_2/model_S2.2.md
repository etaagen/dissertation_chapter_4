Controlled recombination simulation, model analysis, cycle 10
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

gg <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gg.csv") %>% as.data.frame()
gv <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gv.csv") %>% as.data.frame()
pa <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_pa.csv") %>% as.data.frame()
be <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_be.csv") %>% as.data.frame()
# be `value` is 2-fold higher than it should be (ASR does not know how to treat VarA/GenicVarA for DH)
be$value <- be$value/2
qtl <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_qtl.csv") %>% as.data.frame()
qtl_af <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_qtl_af.csv") %>% as.data.frame()
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
Click here: Model, ANOVA tables cycle 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   17.4836  2  0.0001598 ***
    ## Recombination              27.9115  1  1.270e-07 ***
    ## QTL                     69342.8187  1  < 2.2e-16 ***
    ## H2                      47495.4618  1  < 2.2e-16 ***
    ## Repulsion               53002.8238  4  < 2.2e-16 ***
    ## Matrix                   6608.9818  1  < 2.2e-16 ***
    ## QTL_type                  756.9646  1  < 2.2e-16 ***
    ## Map_type:Recombination     29.2230  2  4.511e-07 ***
    ## Map_type:QTL               29.1559  2  4.665e-07 ***
    ## Map_type:H2                89.4358  2  < 2.2e-16 ***
    ## Map_type:Repulsion         52.8155  8  1.171e-08 ***
    ## Map_type:Matrix           173.8639  2  < 2.2e-16 ***
    ## Map_type:QTL_type          18.6364  2  8.978e-05 ***
    ## Recombination:QTL           0.0692  1  0.7924825    
    ## Recombination:H2           19.3244  1  1.103e-05 ***
    ## Recombination:Repulsion    55.1625  4  3.004e-11 ***
    ## Recombination:Matrix      111.8359  1  < 2.2e-16 ***
    ## Recombination:QTL_type      9.8184  1  0.0017277 ** 
    ## QTL:H2                  33239.7357  1  < 2.2e-16 ***
    ## QTL:Repulsion           20158.2196  4  < 2.2e-16 ***
    ## QTL:Matrix                132.6092  1  < 2.2e-16 ***
    ## QTL:QTL_type              177.4185  1  < 2.2e-16 ***
    ## H2:Repulsion             1238.3575  4  < 2.2e-16 ***
    ## H2:Matrix                  50.8971  1  9.734e-13 ***
    ## H2:QTL_type                84.4262  1  < 2.2e-16 ***
    ## Repulsion:Matrix         1053.9478  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        165.3248  4  < 2.2e-16 ***
    ## Matrix:QTL_type            31.1773  1  2.355e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 10, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 10, WT map

![](model_S2.2_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

#### Cycle 10, Pericentromere map

![](model_S2.2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

#### Cycle 10, Chromosome map

![](model_S2.2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

</details>

### Genetic variance

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                 1787.1225  2  < 2.2e-16 ***
    ## Recombination             723.0802  1  < 2.2e-16 ***
    ## QTL                     43292.4871  1  < 2.2e-16 ***
    ## H2                       4559.8736  1  < 2.2e-16 ***
    ## Repulsion                5749.2764  4  < 2.2e-16 ***
    ## Matrix                   3470.2565  1  < 2.2e-16 ***
    ## QTL_type                  716.8011  1  < 2.2e-16 ***
    ## Map_type:Recombination    355.3263  2  < 2.2e-16 ***
    ## Map_type:QTL             1078.7935  2  < 2.2e-16 ***
    ## Map_type:H2               116.3580  2  < 2.2e-16 ***
    ## Map_type:Repulsion        204.2068  8  < 2.2e-16 ***
    ## Map_type:Matrix           144.6909  2  < 2.2e-16 ***
    ## Map_type:QTL_type          63.3410  2  1.761e-14 ***
    ## Recombination:QTL         361.7722  1  < 2.2e-16 ***
    ## Recombination:H2           11.3351  1  0.0007606 ***
    ## Recombination:Repulsion   123.4045  4  < 2.2e-16 ***
    ## Recombination:Matrix       45.2834  1  1.705e-11 ***
    ## Recombination:QTL_type     11.8542  1  0.0005753 ***
    ## QTL:H2                    407.0931  1  < 2.2e-16 ***
    ## QTL:Repulsion            4555.4160  4  < 2.2e-16 ***
    ## QTL:Matrix                 47.2668  1  6.195e-12 ***
    ## QTL:QTL_type              704.4960  1  < 2.2e-16 ***
    ## H2:Repulsion              265.0598  4  < 2.2e-16 ***
    ## H2:Matrix                 466.3194  1  < 2.2e-16 ***
    ## H2:QTL_type                 3.6266  1  0.0568639 .  
    ## Repulsion:Matrix          178.2702  4  < 2.2e-16 ***
    ## Repulsion:QTL_type         13.0644  4  0.0109652 *  
    ## Matrix:QTL_type            32.1576  1  1.422e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 10, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 10, WT map

![](model_S2.2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

#### Cycle 10, Pericentromere map

![](model_S2.2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

#### Cycle 10, Chromosome map

![](model_S2.2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

</details>

### Bulmer effect

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  283.5006  2  < 2.2e-16 ***
    ## Recombination             112.6348  1  < 2.2e-16 ***
    ## QTL                       804.3940  1  < 2.2e-16 ***
    ## H2                      18014.8549  1  < 2.2e-16 ***
    ## Repulsion               35298.4123  4  < 2.2e-16 ***
    ## Matrix                   1272.5126  1  < 2.2e-16 ***
    ## QTL_type                    1.4461  1    0.22915    
    ## Map_type:Recombination     95.3983  2  < 2.2e-16 ***
    ## Map_type:QTL              108.9090  2  < 2.2e-16 ***
    ## Map_type:H2                57.0505  2  4.089e-13 ***
    ## Map_type:Repulsion         45.8449  8  2.544e-07 ***
    ## Map_type:Matrix            54.7342  2  1.302e-12 ***
    ## Map_type:QTL_type          19.6860  2  5.312e-05 ***
    ## Recombination:QTL          36.5648  1  1.477e-09 ***
    ## Recombination:H2           24.3517  1  8.026e-07 ***
    ## Recombination:Repulsion    12.7153  4    0.01275 *  
    ## Recombination:Matrix       32.0628  1  1.493e-08 ***
    ## Recombination:QTL_type      1.5575  1    0.21203    
    ## QTL:H2                   2200.6151  1  < 2.2e-16 ***
    ## QTL:Repulsion            8544.5213  4  < 2.2e-16 ***
    ## QTL:Matrix                          0               
    ## QTL:QTL_type               24.0027  1  9.620e-07 ***
    ## H2:Repulsion             5604.3104  4  < 2.2e-16 ***
    ## H2:Matrix                   0.0000  1    0.99935    
    ## H2:QTL_type                95.0250  1  < 2.2e-16 ***
    ## Repulsion:Matrix          961.7626  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        564.4787  4  < 2.2e-16 ***
    ## Matrix:QTL_type            63.8338  1  1.354e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 10, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 10, WT map

![](model_S2.2_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->

#### Cycle 10, Pericentromere map

![](model_S2.2_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

#### Cycle 10, Chromosome map

![](model_S2.2_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

</details>

### Prediction accuracy

**Linear model, ANOVA table, best and worst simulation parameters by
genetic map and QTL type**

<details>
<summary>
Click here: Model, ANOVA tables cycle 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

    ## boundary (singular) fit: see help('isSingular')

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   22.5904  2  1.243e-05 ***
    ## Recombination              16.4655  1  4.954e-05 ***
    ## QTL                     12777.1672  1  < 2.2e-16 ***
    ## H2                      41250.5010  1  < 2.2e-16 ***
    ## Repulsion                3650.3679  4  < 2.2e-16 ***
    ## Matrix                   2429.6476  1  < 2.2e-16 ***
    ## QTL_type                  231.3494  1  < 2.2e-16 ***
    ## Map_type:Recombination     28.9159  2  5.260e-07 ***
    ## Map_type:QTL               69.9379  2  6.504e-16 ***
    ## Map_type:H2                 6.7443  2   0.034315 *  
    ## Map_type:Repulsion         14.2120  8   0.076405 .  
    ## Map_type:Matrix           137.9069  2  < 2.2e-16 ***
    ## Map_type:QTL_type           9.4468  2   0.008885 ** 
    ## Recombination:QTL          33.6729  1  6.520e-09 ***
    ## Recombination:H2            0.3249  1   0.568694    
    ## Recombination:Repulsion     0.6902  4   0.952537    
    ## Recombination:Matrix       57.0471  1  4.255e-14 ***
    ## Recombination:QTL_type      5.4186  1   0.019923 *  
    ## QTL:H2                   5755.9031  1  < 2.2e-16 ***
    ## QTL:Repulsion             251.4932  4  < 2.2e-16 ***
    ## QTL:Matrix                          0               
    ## QTL:QTL_type               96.0096  1  < 2.2e-16 ***
    ## H2:Repulsion               62.7645  4  7.606e-13 ***
    ## H2:Matrix                 913.4720  1  < 2.2e-16 ***
    ## H2:QTL_type                 2.6504  1   0.103526    
    ## Repulsion:Matrix          129.3546  4  < 2.2e-16 ***
    ## Repulsion:QTL_type         30.3812  4  4.093e-06 ***
    ## Matrix:QTL_type           105.8726  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
<details>
<summary>
Click here: Cycle 10, best and worst simulation parameter prediction
summaries
</summary>

#### Cycle 10, WT map

![](model_S2.2_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-18-4.png)<!-- -->

#### Cycle 10, Pericentromere map

![](model_S2.2_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->

#### Cycle 10, Chromosome map

![](model_S2.2_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

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

![](model_S2.2_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->

#### Pericentromere map

![](model_S2.2_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-24-4.png)<!-- -->

#### Chromosome map

![](model_S2.2_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-25-3.png)<!-- -->![](model_S2.2_files/figure-gfm/unnamed-chunk-25-4.png)<!-- -->

</details>
