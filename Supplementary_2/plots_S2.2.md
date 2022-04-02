Controlled recombination simulation, biparental plots
================

### Notes:

Each simulation rep draws from P1xP2 (“PI\_220431” and “PI\_185715”,
selected from 26 founders for least genetic map lost when filtering for
polymorphisms and greatest SNP density)

Breeding scheme: P1xP2 -&gt; 400 F1s -&gt; 400 DHs, phenotype, serves as
TP and GS cycle 1 starting pop, variance standardized to 1 -&gt; change
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
**Repulsion:** 3:4 **Heritability:** 0.2 or 0.8  
**QTL per Chr:** 2 or 200  
**Relationship matrix:** genomewide or causal variant

#### Load packages and data

``` r
library(tidyverse)
library(data.table)
library(kableExtra)

# anova and summary plot functions
source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/summary_stats_anova.R")
source("/Users/ellietaagen/Desktop/github/cr_simulation/Results/table_view_function.R")

gg <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/bp_gg.csv") %>% as.data.frame()
gv <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/bp_gv.csv") %>% as.data.frame()
pa <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/bp_pa.csv") %>% as.data.frame()
be <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/bp_be.csv") %>% as.data.frame()
# be `value` is 2-fold higher than it should be (ASR does not know how to treat VarA/GenicVarA for DH)
be$value <- be$value/2
qtl <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/bp_qtl.csv") %>% as.data.frame()
qtl_af <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/bp_qtl_af.csv") %>% as.data.frame()
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
**Map\_Type** is the WT, Pericentromere, or Chromosome-wide change to
the genetic map, given a Recombination scale  
**Repulsion** is 1:5 representing different coupling and replusion
ratios

-   3: Random 1/2 of additive effect signs are positive and 1/2 are
    negative for QTL

-   4: 1/2 of additive effect signs are positive and 1/2 are negative
    for QTL, alternating positive or negative each QTL

</details>

### Genetic gain

**Plots**

<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

</details>

### Genetic variance

**Plots**

<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

</details>

### Bulmer effect

**Plots**
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

</details>

### Prediction accuracy

**Plots**

<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

</details>

### QTL fixation

**Plots**

<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

</details>

### QTL positive allele frequency change

**Plots**

<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 200, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.8
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

</details>
<details>
<summary>
Click to expand: QTL/Chr: 2, H^2: 0.2
</summary>

![](master_results_ASR_BP_plots_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

</details>
