Controlled recombination simulation, model analysis, cycle 6 & 10
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

**Linear model, ANOVA table**

<details>
<summary>
Click here: Model, ANOVA tables cycle 6 & 10
</summary>

**Note:** filtered for cycle 6 or cycle 10 observations. See .Rmd file
for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   89.4420  2  < 2.2e-16 ***
    ## Recombination             125.4176  1  < 2.2e-16 ***
    ## QTL                     13597.3675  1  < 2.2e-16 ***
    ## H2                      30875.0986  1  < 2.2e-16 ***
    ## Repulsion               45217.3694  4  < 2.2e-16 ***
    ## Matrix                   8567.8467  1  < 2.2e-16 ***
    ## QTL_type                  195.1159  1  < 2.2e-16 ***
    ## Map_type:Recombination     54.9519  2  1.168e-12 ***
    ## Map_type:QTL               51.1622  2  7.767e-12 ***
    ## Map_type:H2                17.0646  2   0.000197 ***
    ## Map_type:Repulsion         19.2111  8   0.013770 *  
    ## Map_type:Matrix            90.7872  2  < 2.2e-16 ***
    ## Map_type:QTL_type          10.6679  2   0.004825 ** 
    ## Recombination:QTL          31.9720  1  1.564e-08 ***
    ## Recombination:H2            4.9912  1   0.025477 *  
    ## Recombination:Repulsion    29.5366  4  6.081e-06 ***
    ## Recombination:Matrix       69.0208  1  < 2.2e-16 ***
    ## Recombination:QTL_type      1.2630  1   0.261092    
    ## QTL:H2                  13363.7849  1  < 2.2e-16 ***
    ## QTL:Repulsion           20914.5757  4  < 2.2e-16 ***
    ## QTL:Matrix                328.3382  1  < 2.2e-16 ***
    ## QTL:QTL_type                1.4338  1   0.231147    
    ## H2:Repulsion              921.2661  4  < 2.2e-16 ***
    ## H2:Matrix                  77.4231  1  < 2.2e-16 ***
    ## H2:QTL_type                26.1310  1  3.190e-07 ***
    ## Repulsion:Matrix         1579.5358  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        135.0370  4  < 2.2e-16 ***
    ## Matrix:QTL_type            12.8808  1   0.000332 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   24.5680  2  4.625e-06 ***
    ## Recombination              26.0883  1  3.261e-07 ***
    ## QTL                     69683.0142  1  < 2.2e-16 ***
    ## H2                      47470.8522  1  < 2.2e-16 ***
    ## Repulsion               53692.4669  4  < 2.2e-16 ***
    ## Matrix                   6438.8584  1  < 2.2e-16 ***
    ## QTL_type                  705.9874  1  < 2.2e-16 ***
    ## Map_type:Recombination     30.7548  2  2.097e-07 ***
    ## Map_type:QTL               21.9562  2  1.707e-05 ***
    ## Map_type:H2                77.5940  2  < 2.2e-16 ***
    ## Map_type:Repulsion         47.7166  8  1.119e-07 ***
    ## Map_type:Matrix           173.3465  2  < 2.2e-16 ***
    ## Map_type:QTL_type          14.1630  2  0.0008405 ***
    ## Recombination:QTL           2.6505  1  0.1035156    
    ## Recombination:H2           29.0120  1  7.193e-08 ***
    ## Recombination:Repulsion    50.7409  4  2.529e-10 ***
    ## Recombination:Matrix      114.2558  1  < 2.2e-16 ***
    ## Recombination:QTL_type      8.8331  1  0.0029582 ** 
    ## QTL:H2                  33221.7135  1  < 2.2e-16 ***
    ## QTL:Repulsion           20284.4633  4  < 2.2e-16 ***
    ## QTL:Matrix                116.7378  1  < 2.2e-16 ***
    ## QTL:QTL_type              185.3634  1  < 2.2e-16 ***
    ## H2:Repulsion             1214.0012  4  < 2.2e-16 ***
    ## H2:Matrix                  36.1123  1  1.863e-09 ***
    ## H2:QTL_type                79.6111  1  < 2.2e-16 ***
    ## Repulsion:Matrix         1117.0858  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        161.2456  4  < 2.2e-16 ***
    ## Matrix:QTL_type            26.4919  1  2.646e-07 ***
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

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  696.6703  2  < 2.2e-16 ***
    ## Recombination             499.1506  1  < 2.2e-16 ***
    ## QTL                     30191.1671  1  < 2.2e-16 ***
    ## H2                       7659.5367  1  < 2.2e-16 ***
    ## Repulsion                5681.2691  4  < 2.2e-16 ***
    ## Matrix                   6580.7232  1  < 2.2e-16 ***
    ## QTL_type                  476.0269  1  < 2.2e-16 ***
    ## Map_type:Recombination    210.0695  2  < 2.2e-16 ***
    ## Map_type:QTL              257.0303  2  < 2.2e-16 ***
    ## Map_type:H2                11.6328  2   0.002978 ** 
    ## Map_type:Repulsion        140.6820  8  < 2.2e-16 ***
    ## Map_type:Matrix            59.4582  2  1.227e-13 ***
    ## Map_type:QTL_type          23.8410  2  6.653e-06 ***
    ## Recombination:QTL         143.5848  1  < 2.2e-16 ***
    ## Recombination:H2            0.2182  1   0.640406    
    ## Recombination:Repulsion   100.8170  4  < 2.2e-16 ***
    ## Recombination:Matrix       46.9447  1  7.302e-12 ***
    ## Recombination:QTL_type      8.9435  1   0.002785 ** 
    ## QTL:H2                     42.5655  1  6.836e-11 ***
    ## QTL:Repulsion            3966.3256  4  < 2.2e-16 ***
    ## QTL:Matrix                884.2658  1  < 2.2e-16 ***
    ## QTL:QTL_type              296.1741  1  < 2.2e-16 ***
    ## H2:Repulsion              667.7577  4  < 2.2e-16 ***
    ## H2:Matrix                 861.2709  1  < 2.2e-16 ***
    ## H2:QTL_type                25.1489  1  5.307e-07 ***
    ## Repulsion:Matrix          636.9729  4  < 2.2e-16 ***
    ## Repulsion:QTL_type          2.9301  4   0.569583    
    ## Matrix:QTL_type            20.5368  1  5.849e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                 1617.3552  2  < 2.2e-16 ***
    ## Recombination             939.6183  1  < 2.2e-16 ***
    ## QTL                     41806.4056  1  < 2.2e-16 ***
    ## H2                       4447.5608  1  < 2.2e-16 ***
    ## Repulsion                5513.0116  4  < 2.2e-16 ***
    ## Matrix                   3268.5481  1  < 2.2e-16 ***
    ## QTL_type                  738.4682  1  < 2.2e-16 ***
    ## Map_type:Recombination    467.6696  2  < 2.2e-16 ***
    ## Map_type:QTL              999.2793  2  < 2.2e-16 ***
    ## Map_type:H2                98.1838  2  < 2.2e-16 ***
    ## Map_type:Repulsion        202.9028  8  < 2.2e-16 ***
    ## Map_type:Matrix           111.1330  2  < 2.2e-16 ***
    ## Map_type:QTL_type          52.9782  2  3.133e-12 ***
    ## Recombination:QTL         441.5412  1  < 2.2e-16 ***
    ## Recombination:H2           29.9093  1  4.527e-08 ***
    ## Recombination:Repulsion   115.4242  4  < 2.2e-16 ***
    ## Recombination:Matrix       68.2360  1  < 2.2e-16 ***
    ## Recombination:QTL_type      8.6923  1  0.0031956 ** 
    ## QTL:H2                    389.9247  1  < 2.2e-16 ***
    ## QTL:Repulsion            4264.5478  4  < 2.2e-16 ***
    ## QTL:Matrix                 37.8615  1  7.595e-10 ***
    ## QTL:QTL_type              635.8169  1  < 2.2e-16 ***
    ## H2:Repulsion              314.9633  4  < 2.2e-16 ***
    ## H2:Matrix                 417.5357  1  < 2.2e-16 ***
    ## H2:QTL_type                11.9875  1  0.0005356 ***
    ## Repulsion:Matrix          155.2592  4  < 2.2e-16 ***
    ## Repulsion:QTL_type         11.9006  4  0.0181056 *  
    ## Matrix:QTL_type            22.3161  1  2.313e-06 ***
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

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  216.4596  2  < 2.2e-16 ***
    ## Recombination              98.7899  1  < 2.2e-16 ***
    ## QTL                        15.0194  1  0.0001064 ***
    ## H2                      14514.0053  1  < 2.2e-16 ***
    ## Repulsion               49399.2815  4  < 2.2e-16 ***
    ## Matrix                   1285.0020  1  < 2.2e-16 ***
    ## QTL_type                    7.5597  1  0.0059688 ** 
    ## Map_type:Recombination     48.2993  2  3.250e-11 ***
    ## Map_type:QTL               55.5942  2  8.470e-13 ***
    ## Map_type:H2                57.2069  2  3.782e-13 ***
    ## Map_type:Repulsion         39.1829  8  4.546e-06 ***
    ## Map_type:Matrix            15.3014  2  0.0004757 ***
    ## Map_type:QTL_type          13.6109  2  0.0011077 ** 
    ## Recombination:QTL          41.7281  1  1.049e-10 ***
    ## Recombination:H2           32.5880  1  1.139e-08 ***
    ## Recombination:Repulsion    21.4383  4  0.0002592 ***
    ## Recombination:Matrix        2.2666  1  0.1321929    
    ## Recombination:QTL_type      5.5349  1  0.0186412 *  
    ## QTL:H2                   2884.4732  1  < 2.2e-16 ***
    ## QTL:Repulsion           18277.4764  4  < 2.2e-16 ***
    ## QTL:Matrix                165.8114  1  < 2.2e-16 ***
    ## QTL:QTL_type               47.4673  1  5.593e-12 ***
    ## H2:Repulsion             8880.0660  4  < 2.2e-16 ***
    ## H2:Matrix                   0.2135  1  0.6440710    
    ## H2:QTL_type               122.1529  1  < 2.2e-16 ***
    ## Repulsion:Matrix          888.8018  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        860.3565  4  < 2.2e-16 ***
    ## Matrix:QTL_type            93.5076  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  258.3498  2  < 2.2e-16 ***
    ## Recombination             123.3980  1  < 2.2e-16 ***
    ## QTL                       774.5131  1  < 2.2e-16 ***
    ## H2                      17666.1542  1  < 2.2e-16 ***
    ## Repulsion               33662.0688  4  < 2.2e-16 ***
    ## Matrix                   1183.9237  1  < 2.2e-16 ***
    ## QTL_type                    0.6752  1  0.4112584    
    ## Map_type:Recombination     97.1343  2  < 2.2e-16 ***
    ## Map_type:QTL               99.8406  2  < 2.2e-16 ***
    ## Map_type:H2                62.4210  2  2.789e-14 ***
    ## Map_type:Repulsion         38.2611  8  6.737e-06 ***
    ## Map_type:Matrix            47.3400  2  5.251e-11 ***
    ## Map_type:QTL_type          14.5555  2  0.0006907 ***
    ## Recombination:QTL          33.4263  1  7.402e-09 ***
    ## Recombination:H2           26.2315  1  3.028e-07 ***
    ## Recombination:Repulsion     8.2669  4  0.0822779 .  
    ## Recombination:Matrix       27.3481  1  1.699e-07 ***
    ## Recombination:QTL_type      5.3330  1  0.0209256 *  
    ## QTL:H2                   2076.6451  1  < 2.2e-16 ***
    ## QTL:Repulsion            8108.0174  4  < 2.2e-16 ***
    ## QTL:Matrix                          0               
    ## QTL:QTL_type               30.9882  1  2.596e-08 ***
    ## H2:Repulsion             5365.0462  4  < 2.2e-16 ***
    ## H2:Matrix                   0.2239  1  0.6360683    
    ## H2:QTL_type                71.4104  1  < 2.2e-16 ***
    ## Repulsion:Matrix          897.2405  4  < 2.2e-16 ***
    ## Repulsion:QTL_type        565.1629  4  < 2.2e-16 ***
    ## Matrix:QTL_type            65.8398  1  4.891e-16 ***
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

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type)^2 + (1|rep)`

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

ANOVA table, cycle 6:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                  128.4071  2  < 2.2e-16 ***
    ## Recombination              61.7312  1  3.937e-15 ***
    ## QTL                      2715.0770  1  < 2.2e-16 ***
    ## H2                      16414.8840  1  < 2.2e-16 ***
    ## Repulsion                1806.8348  4  < 2.2e-16 ***
    ## Matrix                   2834.3687  1  < 2.2e-16 ***
    ## QTL_type                   15.4883  1  8.302e-05 ***
    ## Map_type:Recombination     53.5982  2  2.298e-12 ***
    ## Map_type:QTL               83.9894  2  < 2.2e-16 ***
    ## Map_type:H2                 6.5914  2   0.037042 *  
    ## Map_type:Repulsion         12.3429  8   0.136546    
    ## Map_type:Matrix            77.6085  2  < 2.2e-16 ***
    ## Map_type:QTL_type           1.3579  2   0.507157    
    ## Recombination:QTL          27.2896  1  1.751e-07 ***
    ## Recombination:H2            0.0037  1   0.951373    
    ## Recombination:Repulsion     5.1547  4   0.271794    
    ## Recombination:Matrix       33.5975  1  6.778e-09 ***
    ## Recombination:QTL_type      1.3339  1   0.248113    
    ## QTL:H2                   1850.9790  1  < 2.2e-16 ***
    ## QTL:Repulsion             319.9080  4  < 2.2e-16 ***
    ## QTL:Matrix                842.1568  1  < 2.2e-16 ***
    ## QTL:QTL_type                9.5621  1   0.001986 ** 
    ## H2:Repulsion               37.7668  4  1.252e-07 ***
    ## H2:Matrix                 581.9232  1  < 2.2e-16 ***
    ## H2:QTL_type                 1.1168  1   0.290602    
    ## Repulsion:Matrix           63.4173  4  5.543e-13 ***
    ## Repulsion:QTL_type         11.5856  4   0.020714 *  
    ## Matrix:QTL_type            64.1690  1  1.142e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ANOVA table, cycle 10:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                   36.9578  2  9.434e-09 ***
    ## Recombination              11.1375  1   0.000846 ***
    ## QTL                     12450.2629  1  < 2.2e-16 ***
    ## H2                      40222.4041  1  < 2.2e-16 ***
    ## Repulsion                3632.1656  4  < 2.2e-16 ***
    ## Matrix                   2253.1764  1  < 2.2e-16 ***
    ## QTL_type                  222.3027  1  < 2.2e-16 ***
    ## Map_type:Recombination     27.6508  2  9.902e-07 ***
    ## Map_type:QTL               69.6119  2  7.655e-16 ***
    ## Map_type:H2                 1.0211  2   0.600162    
    ## Map_type:Repulsion         13.8676  8   0.085283 .  
    ## Map_type:Matrix           126.7120  2  < 2.2e-16 ***
    ## Map_type:QTL_type           9.8018  2   0.007440 ** 
    ## Recombination:QTL          30.6362  1  3.112e-08 ***
    ## Recombination:H2            0.0827  1   0.773681    
    ## Recombination:Repulsion     1.1745  4   0.882274    
    ## Recombination:Matrix       57.4536  1  3.461e-14 ***
    ## Recombination:QTL_type      5.0253  1   0.024980 *  
    ## QTL:H2                   5586.9541  1  < 2.2e-16 ***
    ## QTL:Repulsion             263.5152  4  < 2.2e-16 ***
    ## QTL:Matrix                          0               
    ## QTL:QTL_type               66.8967  1  2.861e-16 ***
    ## H2:Repulsion               64.2368  4  3.726e-13 ***
    ## H2:Matrix                 921.4778  1  < 2.2e-16 ***
    ## H2:QTL_type                 3.9033  1   0.048190 *  
    ## Repulsion:Matrix          115.4210  4  < 2.2e-16 ***
    ## Repulsion:QTL_type         40.5718  4  3.296e-08 ***
    ## Matrix:QTL_type           109.5343  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>

### QTL fixation

**Linear model has `singular fit`, omitted**

### QTL allele frequency change

**Linear model, ANOVA table**

<details>
<summary>
Click here: Model, ANOVA table
</summary>

**Note:** total cycle change observations. See .Rmd file for code.

`response variable ~ (map type + recombination + QTL per Chr + H2 + repulsion + matrix + QTL type + allele)^2 + (1|rep)`

    ## boundary (singular) fit: see help('isSingular')

ANOVA table:

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: value
    ##                              Chisq Df Pr(>Chisq)    
    ## Map_type                    7.6050  2  0.0223149 *  
    ## Recombination               8.7192  1  0.0031487 ** 
    ## QTL                     72069.0154  1  < 2.2e-16 ***
    ## H2                          1.6475  1  0.1993003    
    ## Repulsion               31755.2339  4  < 2.2e-16 ***
    ## Matrix                   2576.8294  1  < 2.2e-16 ***
    ## QTL_type                  344.4908  1  < 2.2e-16 ***
    ## Allele                   4413.3619  2  < 2.2e-16 ***
    ## Map_type:Recombination      6.9795  2  0.0305088 *  
    ## Map_type:QTL                3.1292  2  0.2091685    
    ## Map_type:H2                11.0767  2  0.0039330 ** 
    ## Map_type:Repulsion          2.6763  8  0.9530029    
    ## Map_type:Matrix            18.2048  2  0.0001114 ***
    ## Map_type:QTL_type           0.4092  2  0.8149710    
    ## Map_type:Allele            10.8113  4  0.0287683 *  
    ## Recombination:QTL           0.3534  1  0.5521776    
    ## Recombination:H2            2.8969  1  0.0887507 .  
    ## Recombination:Repulsion     4.4768  4  0.3453020    
    ## Recombination:Matrix       14.4182  1  0.0001464 ***
    ## Recombination:QTL_type      0.4017  1  0.5262103    
    ## Recombination:Allele        1.3715  2  0.5037215    
    ## QTL:H2                   1233.2556  1  < 2.2e-16 ***
    ## QTL:Repulsion           25862.2095  4  < 2.2e-16 ***
    ## QTL:Matrix               1440.5732  1  < 2.2e-16 ***
    ## QTL:QTL_type               77.3361  1  < 2.2e-16 ***
    ## QTL:Allele              20405.0332  2  < 2.2e-16 ***
    ## H2:Repulsion               18.5441  4  0.0009657 ***
    ## H2:Matrix                  64.7589  1  8.465e-16 ***
    ## H2:QTL_type                13.8161  1  0.0002016 ***
    ## H2:Allele               15549.8127  2  < 2.2e-16 ***
    ## Repulsion:Matrix         1383.6361  4  < 2.2e-16 ***
    ## Repulsion:QTL_type         33.9241  4  7.724e-07 ***
    ## Repulsion:Allele         2314.6917  8  < 2.2e-16 ***
    ## Matrix:QTL_type             9.0609  1  0.0026114 ** 
    ## Matrix:Allele             145.4860  2  < 2.2e-16 ***
    ## QTL_type:Allele            30.7960  2  2.055e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

</details>
