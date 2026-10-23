# VAPOR application results

Rebuilt from scratch 2026-10-20 15:37 from this rerun's outputs only.

VAPOR 0.1.0, git 96d90c9. DP base measure Gamma(as = 1, bs = 1); a_dir = b_dir = 0.1; residual-precision floor 1e-8.
Panel: NA filter on the selected series only; both dimensions train 1996Q1-2020Q4 (100 rows), test 2021Q1-2025Q4 (20 rows).
Partial-correlation summaries: mean of per-draw partial correlations pooled over four chains for every Bayesian method; point estimates for VAR + Glasso and the sample covariance. MSPE uses each method's precision point estimate, unchanged.
Main VAR + Glasso and VAR + PrecFact use the full-rank bayesianVARs posterior-mean first stage; the "(rank 3)" rows are the rank-3-truncated versions.
Kastner, VAR + Glasso (both first stages), the sample covariance and the random walk are reused, not refitted.

Simulation sections are absent: the simulations have not been rerun under the new base measure.


# p = 29

## p = 29 (2a) MSPE per chain and mean (SD)

*Source: p29_results.rds*

```
                  Method n_chains                      per_chain            MSPE
       PrecFact (no VAR)        4 0.6367, 0.6404, 0.6732, 0.6911 0.6603 (0.0262)
            VAR + Glasso        1                         0.5604          0.5604
          VAR + PrecFact        4 0.7777, 0.6792, 0.6525, 0.6894 0.6997 (0.0543)
       Kastner (VAR+FSV)        4 0.5614, 0.5714, 0.5704, 0.5677 0.5677 (0.0045)
       Sample covariance        1                         0.6114          0.6114
             Random walk        1                         0.8913          0.8913
          VAPOR (rank 3)        4 0.5631, 0.5790, 0.5746, 0.5750 0.5729 (0.0069)
   VAR + Glasso (rank 3)        1                         0.5602          0.5602
 VAR + PrecFact (rank 3)        4 0.6105, 0.6297, 0.6511, 0.6430 0.6336 (0.0177)
```

## p = 29 (2a) Diebold-Mariano vs VAPOR

HLN corrected, h = 1, per-quarter losses averaged over chains; positive favours VAPOR.

*Source: p29_results.rds*

```
                  Method Difference    SE    DM     p     Won
       PrecFact (no VAR)      0.087 0.027  3.14 0.005 17 / 20
            VAR + Glasso     -0.013 0.024 -0.51 0.615  9 / 20
          VAR + PrecFact      0.127 0.029  4.23 0.000 18 / 20
       Kastner (VAR+FSV)     -0.005 0.024 -0.22 0.832  9 / 20
       Sample covariance      0.038 0.041  0.90 0.377 13 / 20
             Random walk      0.318 0.079  3.93 0.001 16 / 20
          VAPOR (rank 3)         NA    NA    NA    NA       -
   VAR + Glasso (rank 3)     -0.013 0.034 -0.36 0.722 10 / 20
 VAR + PrecFact (rank 3)      0.061 0.019  3.17 0.005 14 / 20
```

## p = 29 (2a) held-out sets in which VAPOR is lower

Of the 100 held-out sets.

*Source: extra_tables.rds*

```
                  Method sets VAPOR_lower share median_gap
       PrecFact (no VAR)  100         100  1.00     0.0878
          VAR + PrecFact  100         100  1.00     0.1211
 VAR + PrecFact (rank 3)  100          95  0.95     0.0613
       Kastner (VAR+FSV)  100          42  0.42    -0.0063
            VAR + Glasso  100          46  0.46    -0.0094
   VAR + Glasso (rank 3)  100          43  0.43    -0.0101
       Sample covariance  100          65  0.65     0.0455
             Random walk  100          98  0.98     0.3092
```

## p = 29 (2b) VAPOR edges, pooled 95% CI rule

pcor is the pooled mean of the per-draw partial correlations over 4000 draws.

*Source: p29_results.rds*

```
                 edge   pcor     lo     hi  Rhat bulk_ESS tail_ESS chains_selecting
       PAYEMS--UNRATE -0.786 -0.884 -0.501 1.004      858     1300                4
 CPIAUCSL--WPSFD49207  0.721  0.465  0.856 1.031      132      295                4
            GS10--BAA  0.653  0.473  0.781 1.000     2897     2416                4
        GDPC1--GPDIC1  0.649  0.214  0.810 1.003     1089      900                4
      FEDFUNDS--TB3MS  0.627  0.369  0.800 1.002     1416     2000                4
       GDPC1--PCECC96  0.585  0.120  0.781 1.005      881     1177                4
      PCECC96--GPDIC1 -0.551 -0.762 -0.136 1.005     1039     1511                4
           GS1--TB3MS  0.534  0.273  0.738 1.003     1030     2039                4
     IPFINAL--IPCONGD  0.532  0.005  0.823 1.006     1436     1889                2
     IPCONGD--IPBUSEQ -0.412 -0.705 -0.015 1.005     1821     2028                4
        GS10--S.P.500  0.409  0.053  0.639 1.002     2523     2156                4
         BAA--S.P.500 -0.378 -0.608 -0.055 1.001     2377     3257                4
            GS1--GS10  0.300  0.037  0.525 1.000     2601     3061                4
```

## p = 29 (2c) all pairs, ranks 1-10 by |pooled mean|

*Source: pairs_p29.rds*

```
                 pair    mean      lo      hi selected rank
       PAYEMS--UNRATE -0.7863 -0.8842 -0.5010     TRUE    1
 CPIAUCSL--WPSFD49207  0.7212  0.4645  0.8557     TRUE    2
            GS10--BAA  0.6530  0.4733  0.7812     TRUE    3
        GDPC1--GPDIC1  0.6490  0.2139  0.8098     TRUE    4
      FEDFUNDS--TB3MS  0.6267  0.3686  0.7995     TRUE    5
       GDPC1--PCECC96  0.5849  0.1199  0.7815     TRUE    6
      PCECC96--GPDIC1 -0.5508 -0.7615 -0.1362     TRUE    7
           GS1--TB3MS  0.5336  0.2733  0.7383     TRUE    8
     IPFINAL--IPCONGD  0.5323  0.0052  0.8232     TRUE    9
     IPFINAL--IPBUSEQ  0.4370 -0.0393  0.7650    FALSE   10
```

## p = 29 (2d) edge counts at thresholds

*Source: p29_results.rds*

```
 threshold VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
      0.05    78      89          134               297               149            140                   155                     157
      0.10    42      49           75               200                69             67                    90                      60
      0.20    20      23           33                64                27             31                    36                      36
      0.30    14      14           17                34                19             22                    21                      21
      0.40    12      10            9                17                16             18                    12                      15
      0.50     9       9            5                12                11             14                     9                      12
```

## p = 29 (2e) edge recovery, per-edge ranks

*Source: p29_results.rds*

```
                     VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
PAYEMS--UNRATE           1      21            5                 4                 2              5                     2                       3
CPIAUCSL--WPSFD49207     2       8            1                 6                12              1                     1                       1
GS10--BAA                3       4            3                 9               284            231                     5                     232
GDPC1--GPDIC1            4       6            8                 7                 7              2                     8                       4
FEDFUNDS--TB3MS          5       5            4                 8                 6              7                     6                       8
GDPC1--PCECC96           6       7           11                11                 8             11                    11                      10
PCECC96--GPDIC1          7      10           38                10                 9             10                    13                       9
GS1--TB3MS               8       9            2                 5                 4              3                     4                       6
IPFINAL--IPCONGD         9       1            6                 1                 1              6                     3                       2
IPCONGD--IPBUSEQ        11       3           44                 2                 3              4                    16                       5
GS10--S.P.500           12      30           30                18               335            318                    22                     217
BAA--S.P.500            13      27           25                13               301            360                    14                     335
GS1--GS10               14      12            9                12                10              9                    10                       7
```

## p = 29 (2e) edge recovery summary

*Source: recheck_tables.rds*

```
            method median_rank in_top_7 in_top_28 in_top_70
             VAPOR           7        7        13        13
 PrecFact (no VAR)           8        6        10        10
    VAR + PrecFact           7        7        10        10
      VAR + Glasso           8        6        10        13
 Sample covariance           8        6        13        13
           Kastner           8        6        12        13
```

## p = 29 (2f) convergence

Rank-normalized split R-hat and bulk/tail ESS over the four chains.

*Source: p29_results.rds*

```
            quantity  Rhat bulk_ESS tail_ESS
              loglik 1.008     1061     1396
               rho_B 0.999     3672     3827
       PAYEMS_UNRATE 1.004      858     1300
 CPIAUCSL_WPSFD49207 1.031      132      295
            GS10_BAA 1.000     2897     2416
        GDPC1_GPDIC1 1.003     1089      900
      FEDFUNDS_TB3MS 1.002     1416     2000
       GDPC1_PCECC96 1.005      881     1177
      PCECC96_GPDIC1 1.005     1039     1511
           GS1_TB3MS 1.003     1030     2039
     IPFINAL_IPCONGD 1.006     1436     1889
     IPCONGD_IPBUSEQ 1.005     1821     2028
        GS10_S.P.500 1.002     2523     2156
         BAA_S.P.500 1.001     2377     3257
            GS1_GS10 1.000     2601     3061
```

## p = 29 (2g) spectral radius of B

*Source: extra_tables.rds*

```
    chain draws median_rho pr_gt_1
 seed 101  1000      0.866  0.0280
 seed 102  1000      0.866  0.0280
 seed 103  1000      0.863  0.0260
 seed 104  1000      0.863  0.0270
   pooled  4000      0.865  0.0272
```

## p = 29 (2h) one-step-ahead posterior predictive check

*Source: p29_results.rds*

```
    statistic observed rep_mean     lo     hi   ppp
      max_eig   12.892   10.902  9.320 12.046 0.000
         eig2    3.408    3.407  2.742  4.166 0.458
 mean_abs_cor    0.385    0.307  0.260  0.348 0.000
  q90_abs_cor    0.794    0.654  0.568  0.722 0.000
 max_abs_pcor    0.924    0.847  0.752  0.923 0.025
 prop_pcor_02    0.840    0.753  0.672  0.823 0.002
      mean_sd    1.004    1.016  0.967  1.070 0.688
     kurtosis   17.581   10.825  9.284 12.572 0.000
    resid_ac1    0.188   -0.008 -0.051  0.035 0.000
```

## p = 29 (2i) residual-precision floor, per series

*Source: p29_results.rds*

```
 seed series iterations_hit first_iteration
  101 (none)              0              NA
  102 (none)              0              NA
  103 (none)              0              NA
  104 (none)              0              NA
```

## p = 29 (2i) residual-precision floor, iterations with any floor

floor_iters counts iterations, burn-in included.

*Source: p29_results.rds*

```
 seed series_floored series_iteration_hits iterations_with_any_floor earliest_iteration   nrun  burn
  101              0                     0                         0                 NA 100000 20000
  102              0                     0                         0                 NA 100000 20000
  103              0                     0                         0                 NA 100000 20000
  104              0                     0                         0                 NA 100000 20000
```

## p = 29 (2j) sensitivity sweep

One 100,000-iteration chain per setting; B stored for every fit.

*Source: sweeps_p29_final.rds*

```
 block  setting seed        source   MSPE edges median_rho pr_rho_gt_1 minutes series_floored
     k      k=5  101   fitted here 0.6286     6      0.865       0.030     0.9              0
     k      k=8  101   fitted here 0.5761    13      0.884       0.028     1.0              0
     k     k=10  101 p29_final.rds 0.5631    12      0.858       0.020     1.1              0
     k     k=15  101   fitted here 0.5456    11      0.810       0.045     1.3              0
  rank   rank=2  101   fitted here 0.5777    13      0.882       0.043     0.9              0
  rank   rank=3  101 p29_final.rds 0.5631    12      0.858       0.020     1.1              0
  rank   rank=4  101   fitted here 0.5803    12      0.824       0.052     1.4              0
   a_L  a_L=0.3  101   fitted here 0.5637    10      0.869       0.038     1.1              0
   a_L  a_L=0.5  101 p29_final.rds 0.5631    12      0.858       0.020     1.1              0
   a_L  a_L=1.0  101   fitted here 0.5783    15      0.862       0.025     1.1              0
  seed seed=101  101 p29_final.rds 0.5631    12      0.858       0.020     1.1              0
  seed seed=102  102 p29_final.rds 0.5790    13      0.865       0.030     1.1              0
  seed seed=103  103 p29_final.rds 0.5746    12      0.863       0.038     1.1              0
  seed seed=104  104 p29_final.rds 0.5750    13      0.860       0.022     1.1              0
  seed seed=105  105   fitted here 0.5731    12      0.866       0.025     1.1              0
  seed seed=106  106   fitted here 0.5768    13      0.862       0.050     1.1              0
```

## p = 29 (2j) edge stability over distinct fits

The reported configuration is counted once, not three times.

*Source: edge_stability.rds*

```
                 edge n_settings out_of reported
 CPIAUCSL--WPSFD49207          8      8     TRUE
      FEDFUNDS--TB3MS          8      8     TRUE
           GS1--TB3MS          8      8     TRUE
       PAYEMS--UNRATE          8      8     TRUE
         BAA--S.P.500          7      8     TRUE
        GDPC1--GPDIC1          7      8     TRUE
            GS1--GS10          7      8     TRUE
            GS10--BAA          7      8     TRUE
        GS10--S.P.500          7      8     TRUE
       GDPC1--PCECC96          6      8     TRUE
      PCECC96--GPDIC1          6      8     TRUE
     IPCONGD--IPBUSEQ          5      8     TRUE
     IPFINAL--IPCONGD          4      8     TRUE
       BAA--BUSLOANSx          1      8    FALSE
    OILPRICEx--M2REAL          1      8    FALSE
       PAYEMS--M2REAL          1      8    FALSE
     S.P.500--VIXCLSx          1      8    FALSE
```

## p = 29 (2j) edges selected by each of seeds 101-106

*Source: edge_stability.rds*

```
                 edge n_seeds out_of reported
         BAA--S.P.500       6      6     TRUE
 CPIAUCSL--WPSFD49207       6      6     TRUE
      FEDFUNDS--TB3MS       6      6     TRUE
        GDPC1--GPDIC1       6      6     TRUE
       GDPC1--PCECC96       6      6     TRUE
            GS1--GS10       6      6     TRUE
           GS1--TB3MS       6      6     TRUE
            GS10--BAA       6      6     TRUE
        GS10--S.P.500       6      6     TRUE
     IPCONGD--IPBUSEQ       6      6     TRUE
       PAYEMS--UNRATE       6      6     TRUE
      PCECC96--GPDIC1       6      6     TRUE
     IPFINAL--IPCONGD       3      6     TRUE
```

## p = 29 (2k) chain length

Configs B, A and C refitted here; 100000/20000/80 reuses the main fits.

*Source: chain_length_p29.rds*

```
 iterations  nrun_burn_thin pairs_selected max_Rhat median_Rhat min_bulk_ESS pooled_edges            MSPE
       7500     7500/2500/5             13    1.094       1.031           28           13 0.5754 (0.0117)
      20000  20000/10000/10             13    1.084       1.010           61           12 0.5815 (0.0010)
      50000  50000/10000/40             13    1.033       1.003          105           13 0.5754 (0.0054)
     100000 100000/20000/80             13    1.031       1.003          132           13 0.5729 (0.0069)
```

## p = 29 (2l) minutes per chain

*Source: p29_results.rds*

```
                method seed    status minutes error
                 VAPOR  101 completed    1.09  <NA>
                 VAPOR  102 completed    1.09  <NA>
                 VAPOR  103 completed    1.09  <NA>
                 VAPOR  104 completed    1.08  <NA>
              PrecFact  101 completed    0.45  <NA>
              PrecFact  102 completed    0.44  <NA>
              PrecFact  103 completed    0.44  <NA>
              PrecFact  104 completed    0.44  <NA>
          VAR+PrecFact  101 completed    0.44  <NA>
          VAR+PrecFact  102 completed    0.45  <NA>
          VAR+PrecFact  103 completed    0.42  <NA>
          VAR+PrecFact  104 completed    0.42  <NA>
 VAR+PrecFact (rank 3)  101 completed    0.41  <NA>
 VAR+PrecFact (rank 3)  102 completed    0.38  <NA>
 VAR+PrecFact (rank 3)  103 completed    0.38  <NA>
 VAR+PrecFact (rank 3)  104 completed    0.37  <NA>
```

## p = 29 (2m) condition number of the training-sample covariance

*Source: extra_tables.rds*

```
  p T_train rank min_eigenvalue max_eigenvalue condition
 29     100   29       0.002297          12.88   5607.67
```

## p = 29 runs of >= 5 identical consecutive draws

*Source: stuck_runs.rds*

```
     result
 none found
```


# p = 100

## p = 100 (2a) MSPE per chain and mean (SD)

*Source: p100_results.rds*

```
                  Method n_chains                              per_chain                MSPE
       PrecFact (no VAR)        4         0.6778, 0.6576, 0.6797, 0.6494     0.6661 (0.0150)
            VAR + Glasso        1                                 0.5918              0.5918
          VAR + PrecFact        4         0.6668, 0.7279, 0.6667, 0.6477     0.6773 (0.0349)
       Kastner (VAR+FSV)        4 270.7653, 106.9019, 165.1070, 332.3876 218.7904 (101.6639)
       Sample covariance        1                                 1.8869              1.8869
             Random walk        1                                 1.1128              1.1128
          VAPOR (rank 3)        4         0.5649, 0.5610, 0.5595, 0.5676     0.5632 (0.0037)
   VAR + Glasso (rank 3)        1                                 0.5843              0.5843
 VAR + PrecFact (rank 3)        4         0.7436, 0.6704, 0.6697, 0.7191     0.7007 (0.0368)
```

## p = 100 (2a) Diebold-Mariano vs VAPOR

HLN corrected, h = 1, per-quarter losses averaged over chains; positive favours VAPOR.

*Source: p100_results.rds*

```
                  Method Difference     SE   DM     p     Won
       PrecFact (no VAR)      0.103  0.024 4.11 0.001 20 / 20
            VAR + Glasso      0.029  0.018 1.56 0.134 13 / 20
          VAR + PrecFact      0.114  0.027 4.06 0.001 18 / 20
       Kastner (VAR+FSV)    218.227 36.166 5.88 0.000 20 / 20
       Sample covariance      1.324  0.187 6.89 0.000 20 / 20
             Random walk      0.550  0.126 4.26 0.000 19 / 20
          VAPOR (rank 3)         NA     NA   NA    NA       -
   VAR + Glasso (rank 3)      0.021  0.012 1.73 0.101 12 / 20
 VAR + PrecFact (rank 3)      0.137  0.031 4.33 0.000 20 / 20
```

## p = 100 (2a) held-out sets in which VAPOR is lower

Of the 100 held-out sets.

*Source: extra_tables.rds*

```
                  Method sets VAPOR_lower share median_gap
       PrecFact (no VAR)  100         100  1.00     0.1007
          VAR + PrecFact  100         100  1.00     0.1146
 VAR + PrecFact (rank 3)  100         100  1.00     0.1394
       Kastner (VAR+FSV)  100         100  1.00   196.6706
            VAR + Glasso  100          81  0.81     0.0279
   VAR + Glasso (rank 3)  100          75  0.75     0.0212
       Sample covariance  100         100  1.00     1.3257
             Random walk  100         100  1.00     0.5140
```

## p = 100 (2b) VAPOR edges, pooled 95% CI rule

pcor is the pooled mean of the per-draw partial correlations over 4000 draws.

*Source: p100_results.rds*

```
                edge   pcor     lo     hi  Rhat bulk_ESS tail_ESS chains_selecting
    AMDMNOx--ANDENOx  0.807  0.582  0.887 1.005      982      457                4
   UNRATE--UNRATESTx  0.751  0.467  0.889 1.001     2820     3016                4
         TCU--CUMFNS  0.750  0.423  0.879 1.003     2105     2190                4
       INDPRO--IPMAT  0.699  0.331  0.866 1.002     2738     1696                4
    IPFINAL--IPCONGD  0.691  0.431  0.845 1.013      813     1449                4
            AAA--BAA  0.667  0.413  0.832 1.002     2075     2820                4
      M1REAL--M2REAL  0.648  0.123  0.885 1.022      194     1924                4
     FEDFUNDS--TB3MS  0.639  0.310  0.822 1.004      850     1121                4
     MANEMP--DMANEMP  0.632  0.336  0.824 1.000     2480     2226                4
      USGOOD--USCONS  0.630  0.230  0.823 1.004     1276      604                4
 UNRATELTx--UEMP27OV  0.616  0.272  0.813 1.003     1649     2451                4
    IPFINAL--IPBUSEQ  0.605  0.285  0.792 1.005     1813     1597                4
    IPCONGD--IPBUSEQ -0.573 -0.759 -0.286 1.003     2641     1649                4
  CPIAUCSL--CPITRNSL  0.558  0.226  0.797 1.003     1233     1818                4
     PCECC96--PCESVx  0.550  0.090  0.785 1.020      221     1429                4
       USPRIV--USLAH  0.521  0.002  0.844 1.011      469      629                2
       HOUST--HOUSTS  0.512  0.077  0.753 1.002     1165     1183                4
       HOUST--HOUSTW  0.504  0.081  0.742 1.003     1306     1244                4
      PCECC96--PCDGx  0.492  0.053  0.751 1.009      644     1523                4
     PPIACO--WPSID61  0.486  0.053  0.751 1.001     2105     3196                4
          TB3MS--GS1  0.467  0.173  0.735 1.008      356      524                4
           GS5--GS10  0.422  0.042  0.745 1.023      129      526                4
   PCECTPI--CPITRNSL  0.419  0.032  0.698 1.003     1211     1505                4
       PCDGx--PCESVx -0.388 -0.668 -0.014 1.060       44     1246                3
   GS5--MORTGAGE30US  0.355  0.020  0.633 1.002     1071     2040                4
            GS1--GS5  0.304  0.016  0.603 1.020      238      808                4
```

## p = 100 (2c) all pairs, ranks 1-10 by |pooled mean|

*Source: pairs_p100.rds*

```
              pair   mean     lo     hi selected rank
  AMDMNOx--ANDENOx 0.8066 0.5818 0.8866     TRUE    1
 UNRATE--UNRATESTx 0.7507 0.4667 0.8887     TRUE    2
       TCU--CUMFNS 0.7500 0.4226 0.8785     TRUE    3
     INDPRO--IPMAT 0.6988 0.3305 0.8662     TRUE    4
  IPFINAL--IPCONGD 0.6910 0.4314 0.8449     TRUE    5
          AAA--BAA 0.6675 0.4127 0.8319     TRUE    6
    M1REAL--M2REAL 0.6482 0.1235 0.8851     TRUE    7
   FEDFUNDS--TB3MS 0.6387 0.3104 0.8224     TRUE    8
   MANEMP--DMANEMP 0.6319 0.3363 0.8240     TRUE    9
    USGOOD--USCONS 0.6298 0.2296 0.8232     TRUE   10
```

## p = 100 (2d) edge counts at thresholds

*Source: p100_results.rds*

```
 threshold VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
      0.05   409     228          476              4712              1291           1112                   492                    1297
      0.10   148      92          196              4474               416            350                   235                     420
      0.20    64      55           65              4001               137            100                    59                     132
      0.30    41      41           31              3490                79             48                    30                      74
      0.40    27      35           20              2938                57             32                    18                      43
      0.50    19      30           10              2400                38             22                     8                      27
```

## p = 100 (2e) edge recovery, per-edge ranks

*Source: p100_results.rds*

```
                    VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
AMDMNOx--ANDENOx        1    2655            2                23               220            263                     1                     228
UNRATE--UNRATESTx       2      16           31                66                39             27                    26                      45
TCU--CUMFNS             3      35           37               665                 1             52                    24                     722
INDPRO--IPMAT           4      13           23              2716                 7             12                    18                      17
IPFINAL--IPCONGD        5      28           17               601                 6             16                    16                       5
AAA--BAA                6    1490           15               586              1155            164                     8                    3723
M1REAL--M2REAL          7     247           20              3176                50             45                    17                       4
FEDFUNDS--TB3MS         8     754           10               132                16             56                     3                      74
MANEMP--DMANEMP         9      31           12                31                26              5                    15                      15
USGOOD--USCONS         10      17           28                10                 5              9                    36                       9
UNRATELTx--UEMP27OV    11      40            1              3387              2119           4007                     2                    3874
IPFINAL--IPBUSEQ       12      37           55              1595                 8             18                    40                       8
IPCONGD--IPBUSEQ       13      33         1370               867                 3              6                  1270                       2
CPIAUCSL--CPITRNSL     15      53           11              3924                37             78                     9                      33
PCECC96--PCESVx        16       5           26                22                46             51                    20                      36
USPRIV--USLAH          17      93           60              3841                63            199                    66                     129
HOUST--HOUSTS          18      29           14               914               857            483                    11                    1451
HOUST--HOUSTW          19      39           13              3927               455            462                    12                     367
PCECC96--PCDGx         20      12           39               123                32             15                    62                      13
PPIACO--WPSID61        22      54            7              2978                47             77                    10                    1445
TB3MS--GS1             23     760            5              3660                15             62                     6                      73
GS5--GS10              25    1753            8                18                30              7                     4                       6
PCECTPI--CPITRNSL      26     146           43              1931                57             39                    60                      72
PCDGx--PCESVx          28      15         1370                 8                51            128                  1270                      31
GS5--MORTGAGE30US      32    1619           62              2275              3500           2599                    59                    4057
GS1--GS5               39      97           33              1293                74            112                    34                      48
```

## p = 100 (2e) edge recovery summary

*Source: recheck_tables.rds*

```
            method median_rank in_top_25 in_top_100 in_top_250
             VAPOR        14.0        22         26         26
 PrecFact (no VAR)        42.5         8         20         21
    VAR + PrecFact        54.0         8         17         21
      VAR + Glasso        21.5        14         24         24
 Sample covariance       890.5         5          7          9
           Kastner        39.5         6         18         20
```

## p = 100 (2f) convergence

Rank-normalized split R-hat and bulk/tail ESS over the four chains.

*Source: p100_results.rds*

```
           quantity  Rhat bulk_ESS tail_ESS
             loglik 1.011      744     2213
              rho_B 1.059       54       40
    AMDMNOx_ANDENOx 1.005      982      457
   UNRATE_UNRATESTx 1.001     2820     3016
         TCU_CUMFNS 1.003     2105     2190
       INDPRO_IPMAT 1.002     2738     1696
    IPFINAL_IPCONGD 1.013      813     1449
            AAA_BAA 1.002     2075     2820
      M1REAL_M2REAL 1.022      194     1924
     FEDFUNDS_TB3MS 1.004      850     1121
     MANEMP_DMANEMP 1.000     2480     2226
      USGOOD_USCONS 1.004     1276      604
 UNRATELTx_UEMP27OV 1.003     1649     2451
    IPFINAL_IPBUSEQ 1.005     1813     1597
    IPCONGD_IPBUSEQ 1.003     2641     1649
  CPIAUCSL_CPITRNSL 1.003     1233     1818
     PCECC96_PCESVx 1.020      221     1429
       USPRIV_USLAH 1.011      469      629
       HOUST_HOUSTS 1.002     1165     1183
       HOUST_HOUSTW 1.003     1306     1244
      PCECC96_PCDGx 1.009      644     1523
     PPIACO_WPSID61 1.001     2105     3196
          TB3MS_GS1 1.008      356      524
           GS5_GS10 1.023      129      526
   PCECTPI_CPITRNSL 1.003     1211     1505
       PCDGx_PCESVx 1.060       44     1246
   GS5_MORTGAGE30US 1.002     1071     2040
            GS1_GS5 1.020      238      808
```

## p = 100 (2g) spectral radius of B

*Source: extra_tables.rds*

```
    chain draws median_rho pr_gt_1
 seed 101  1000      0.785  0.1390
 seed 102  1000      0.799  0.1580
 seed 103  1000      0.785  0.1300
 seed 104  1000      0.977  0.4660
   pooled  4000      0.818  0.2232
```

## p = 100 (2h) one-step-ahead posterior predictive check

*Source: p100_results.rds*

```
    statistic observed rep_mean     lo     hi   ppp
      max_eig   41.167   39.755 38.032 41.240 0.032
         eig2    9.699    9.072  8.042 10.247 0.145
 mean_abs_cor    0.352    0.327  0.308  0.344 0.005
  q90_abs_cor    0.736    0.715  0.693  0.734 0.010
 max_abs_pcor    0.995    0.968  0.920  0.997 0.050
 prop_pcor_02    0.163    0.250  0.073  0.376 0.855
      mean_sd    1.004    1.015  0.991  1.046 0.812
     kurtosis   19.294   14.636 13.751 15.360 0.000
    resid_ac1    0.208   -0.012 -0.034  0.012 0.000
```

## p = 100 (2i) residual-precision floor, per series

*Source: p100_results.rds*

```
 seed series iterations_hit first_iteration
  101 (none)              0              NA
  102 (none)              0              NA
  103 (none)              0              NA
  104 (none)              0              NA
```

## p = 100 (2i) residual-precision floor, iterations with any floor

floor_iters counts iterations, burn-in included.

*Source: p100_results.rds*

```
 seed series_floored series_iteration_hits iterations_with_any_floor earliest_iteration   nrun  burn
  101              0                     0                         0                 NA 100000 20000
  102              0                     0                         0                 NA 100000 20000
  103              0                     0                         0                 NA 100000 20000
  104              0                     0                         0                 NA 100000 20000
```

## p = 100 (2j) sensitivity sweep

One 100,000-iteration chain per setting; B stored for every fit.

*Source: sweeps_p100_final.rds*

```
 block  setting seed         source   MSPE edges median_rho pr_rho_gt_1 minutes series_floored
     k     k=15  101    fitted here 0.6063    14      0.784       0.155    17.8              0
     k     k=25  101 p100_final.rds 0.5649    29      0.779       0.132    14.0              0
     k     k=35  101    fitted here 0.5648    34      0.851       0.162    20.2              0
  rank   rank=2  101    fitted here 0.5656    37      0.777       0.048    12.9              0
  rank   rank=3  101 p100_final.rds 0.5649    29      0.779       0.132    14.0              0
  rank   rank=4  101    fitted here 0.5636    26      0.916       0.245    28.0              0
   a_L  a_L=0.3  101    fitted here 0.5708    36      0.823       0.185    19.2              0
   a_L  a_L=0.5  101 p100_final.rds 0.5649    29      0.779       0.132    14.0              0
   a_L  a_L=1.0  101    fitted here 0.6073    35      0.824       0.128    19.0              0
  seed seed=101  101 p100_final.rds 0.5649    29      0.779       0.132    14.0              0
  seed seed=102  102 p100_final.rds 0.5610    28      0.792       0.150    14.0              0
  seed seed=103  103 p100_final.rds 0.5595    26      0.793       0.112    14.0              0
  seed seed=104  104 p100_final.rds 0.5676    27      0.983       0.482    14.0              0
  seed seed=105  105    fitted here 0.5609    28      0.956       0.452    19.1              0
  seed seed=106  106    fitted here 0.5649    38      0.834       0.242    12.0              0
```

## p = 100 (2j) edge stability over distinct fits

The reported configuration is counted once, not three times.

*Source: edge_stability.rds*

```
                 edge n_settings out_of reported
   CPIAUCSL--CPITRNSL          7      7     TRUE
      FEDFUNDS--TB3MS          7      7     TRUE
            GS5--GS10          7      7     TRUE
      PPIACO--WPSID61          7      7     TRUE
           TB3MS--GS1          7      7     TRUE
    UNRATE--UNRATESTx          7      7     TRUE
       USGOOD--USCONS          7      7     TRUE
             AAA--BAA          6      7     TRUE
     AMDMNOx--ANDENOx          6      7     TRUE
             GS1--GS5          6      7     TRUE
        HOUST--HOUSTS          6      7     TRUE
        HOUST--HOUSTW          6      7     TRUE
        INDPRO--IPMAT          6      7     TRUE
     IPCONGD--IPBUSEQ          6      7     TRUE
      MANEMP--DMANEMP          6      7     TRUE
        GDPC1--GPDIC1          5      7    FALSE
     IPFINAL--IPBUSEQ          5      7     TRUE
     IPFINAL--IPCONGD          5      7     TRUE
       M1REAL--M2REAL          5      7     TRUE
      PCECC96--PCESVx          5      7     TRUE
    PCECTPI--CPITRNSL          5      7     TRUE
  UNRATELTx--UEMP27OV          5      7     TRUE
    GS5--MORTGAGE30US          4      7     TRUE
       HOUSTS--HOUSTW          4      7    FALSE
       PCECC96--PCDGx          4      7     TRUE
       USTPU--USTRADE          4      7    FALSE
         GDPC1--GCEC1          3      7    FALSE
       GDPC1--IMPGSC1          3      7    FALSE
      GPDIC1--IMPGSC1          3      7    FALSE
    IPDMAT--IPMANSICS          3      7    FALSE
        PCDGx--PCESVx          3      7     TRUE
   S.P.500--NASDAQCOM          3      7    FALSE
          TCU--CUMFNS          3      7     TRUE
    UNRATE--UNRATELTx          3      7    FALSE
 UNRATESTx--UNRATELTx          3      7    FALSE
        USPRIV--USLAH          3      7     TRUE
     EXPGSC1--IMPGSC1          2      7    FALSE
          FPIx--PNFIx          2      7    FALSE
       GDPC1--EXPGSC1          2      7    FALSE
        GDPC1--PCESVx          2      7    FALSE
      GPDIC1--EXPGSC1          2      7    FALSE
        GPDIC1--GCEC1          2      7    FALSE
    IPNMAT--IPMANSICS          2      7    FALSE
       PAYEMS--USGOVT          2      7    FALSE
       PCESVx--GPDIC1          2      7    FALSE
      PCESVx--IMPGSC1          2      7    FALSE
        PCESVx--PCNDx          2      7    FALSE
   USSTHPI--SPCS10RSA          2      7    FALSE
    BUSINVx--ISRATIOx          1      7    FALSE
       GCEC1--IMPGSC1          1      7    FALSE
        HOUST--PERMIT          1      7    FALSE
       PAYEMS--USPRIV          1      7    FALSE
        PCDGx--GPDIC1          1      7    FALSE
         PCDGx--PCNDx          1      7    FALSE
        PCDGx--RSAFSx          1      7    FALSE
    PCECTPI--PCEPILFE          1      7    FALSE
   PPIACO--WPSFD49207          1      7    FALSE
```

## p = 100 (2j) edges selected by each of seeds 101-106

*Source: edge_stability.rds*

```
                 edge n_seeds out_of reported
             AAA--BAA       6      6     TRUE
     AMDMNOx--ANDENOx       6      6     TRUE
   CPIAUCSL--CPITRNSL       6      6     TRUE
      FEDFUNDS--TB3MS       6      6     TRUE
             GS1--GS5       6      6     TRUE
            GS5--GS10       6      6     TRUE
        HOUST--HOUSTS       6      6     TRUE
        HOUST--HOUSTW       6      6     TRUE
        INDPRO--IPMAT       6      6     TRUE
     IPCONGD--IPBUSEQ       6      6     TRUE
     IPFINAL--IPBUSEQ       6      6     TRUE
     IPFINAL--IPCONGD       6      6     TRUE
       M1REAL--M2REAL       6      6     TRUE
      MANEMP--DMANEMP       6      6     TRUE
       PCECC96--PCDGx       6      6     TRUE
      PCECC96--PCESVx       6      6     TRUE
    PCECTPI--CPITRNSL       6      6     TRUE
      PPIACO--WPSID61       6      6     TRUE
           TB3MS--GS1       6      6     TRUE
          TCU--CUMFNS       6      6     TRUE
    UNRATE--UNRATESTx       6      6     TRUE
  UNRATELTx--UEMP27OV       6      6     TRUE
       USGOOD--USCONS       6      6     TRUE
        PCDGx--PCESVx       5      6     TRUE
    GS5--MORTGAGE30US       4      6     TRUE
        USPRIV--USLAH       4      6     TRUE
        GDPC1--GPDIC1       3      6    FALSE
      USGOOD--DMANEMP       3      6    FALSE
    IPDMAT--IPMANSICS       2      6    FALSE
    IPNMAT--IPMANSICS       2      6    FALSE
 UNRATESTx--UNRATELTx       2      6    FALSE
       USTPU--USTRADE       2      6    FALSE
     EXPGSC1--IMPGSC1       1      6    FALSE
       GDPC1--EXPGSC1       1      6    FALSE
         GDPC1--GCEC1       1      6    FALSE
       GDPC1--IMPGSC1       1      6    FALSE
        GDPC1--PCESVx       1      6    FALSE
      GPDIC1--EXPGSC1       1      6    FALSE
      GPDIC1--IMPGSC1       1      6    FALSE
       PCESVx--GPDIC1       1      6    FALSE
      PCESVx--IMPGSC1       1      6    FALSE
        PCESVx--PCNDx       1      6    FALSE
   S.P.500--NASDAQCOM       1      6    FALSE
```

## p = 100 (2k) chain length

Configs B, A and C refitted here; 100000/20000/80 reuses the main fits.

*Source: chain_length_p100.rds*

```
 iterations  nrun_burn_thin pairs_selected max_Rhat median_Rhat min_bulk_ESS pooled_edges            MSPE
       7500     7500/2500/5             48    1.558       1.103            7           27 0.5693 (0.0097)
      20000  20000/10000/10             52    1.575       1.095            7           27 0.5626 (0.0054)
      50000  50000/10000/40             35    1.414       1.015            8           24 0.5667 (0.0031)
     100000 100000/20000/80             32    1.233       1.004           13           26 0.5632 (0.0037)
```

## p = 100 (2l) minutes per chain

*Source: p100_results.rds*

```
                method seed    status minutes error
                 VAPOR  101 completed   13.99  <NA>
                 VAPOR  102 completed   13.98  <NA>
                 VAPOR  103 completed   14.00  <NA>
                 VAPOR  104 completed   14.00  <NA>
              PrecFact  101 completed    3.41  <NA>
              PrecFact  102 completed    3.41  <NA>
              PrecFact  103 completed    3.40  <NA>
              PrecFact  104 completed    3.46  <NA>
          VAR+PrecFact  101 completed    3.47  <NA>
          VAR+PrecFact  102 completed    3.47  <NA>
          VAR+PrecFact  103 completed    3.37  <NA>
          VAR+PrecFact  104 completed    3.37  <NA>
 VAR+PrecFact (rank 3)  101 completed    3.36  <NA>
 VAR+PrecFact (rank 3)  102 completed    3.38  <NA>
 VAR+PrecFact (rank 3)  103 completed    3.40  <NA>
 VAR+PrecFact (rank 3)  104 completed    3.39  <NA>
```

## p = 100 (2m) condition number of the training-sample covariance

*Source: extra_tables.rds*

```
   p T_train rank min_eigenvalue max_eigenvalue                   condition
 100     100   99      7.308e-16          41.13 singular (not of full rank)
```

## p = 100 runs of >= 5 identical consecutive draws

*Source: stuck_runs.rds*

```
     result
 none found
```


# p = 29 only

## p = 29 (3) the six pairs among FEDFUNDS, TB3MS, GS1, GS10

*Source: pairs_p29.rds*

```
            pair    mean      lo     hi selected rank
 FEDFUNDS--TB3MS  0.6267  0.3686 0.7995     TRUE    5
   FEDFUNDS--GS1  0.0894 -0.1325 0.4261    FALSE   45
  FEDFUNDS--GS10 -0.0463 -0.2516 0.1538    FALSE   86
      GS1--TB3MS  0.5336  0.2733 0.7383     TRUE    8
     GS10--TB3MS -0.0394 -0.2607 0.1904    FALSE  102
       GS1--GS10  0.3004  0.0365 0.5252     TRUE   14
```

## p = 29 (3) rank and k table (Table S5)

*Source: ranks_p29_final.rds*

```
 rank  k    pr_expl   med_rho   rep_eig         ppp      mspe    edges
 Full  5 0.72000000 1.0923592  6.343882 0.000000000 0.5832016  5.00000
 Full  8 0.74083333 1.0901096  6.836566 0.000000000 0.5690547  7.00000
 Full 10 0.71583333 1.0779173  7.208582 0.000000000 0.5646292  7.00000
 Full 15 0.68166667 1.0692632  8.132276 0.000000000 0.5580581  7.00000
    2  5 0.02416667 0.8902483 11.210664 0.000000000 0.6394312  7.00000
    2  8 0.01750000 0.8869360 11.260493 0.003333333 0.5766957 13.00000
    2 10 0.02916667 0.8815878 10.884214 0.000000000 0.5776061 13.00000
    2 15 0.04083333 0.8190351 10.719980 0.000000000 0.5309051 12.66667
    3  5 0.03750000 0.8908412 11.357049 0.001666667 0.6285600  6.00000
    3  8 0.02583333 0.8767766 11.141614 0.000000000 0.5762521 13.00000
    3 10 0.02916667 0.8621043 10.867921 0.000000000 0.5722450 12.33333
    3 15 0.04000000 0.8144817 10.551531 0.000000000 0.5478724 11.66667
    4  5 0.03916667 0.8692343 11.288303 0.000000000 0.6330320  7.00000
    4  8 0.04916667 0.8753289 11.110432 0.000000000 0.5806847 12.33333
    4 10 0.04500000 0.8230429 11.021056 0.000000000 0.5786284 12.00000
    4 15 0.04750000 0.8009115 10.586905 0.001666667 0.5508561 11.66667
```


# p = 100 only

## p = 100 (4) rate block

*Source: pairs_p100.rds*

```
            pair    mean      lo     hi selected rank
 FEDFUNDS--TB3MS  0.6387  0.3104 0.8224     TRUE    8
   FEDFUNDS--GS1  0.1151 -0.3353 0.5149    FALSE  115
   FEDFUNDS--GS5 -0.1007 -0.4164 0.1297    FALSE  144
  FEDFUNDS--GS10 -0.0295 -0.2312 0.1357    FALSE  994
      TB3MS--GS1  0.4673  0.1726 0.7348     TRUE   23
      TB3MS--GS5  0.0127 -0.2044 0.2419    FALSE 2513
     TB3MS--GS10 -0.0048 -0.1651 0.1614    FALSE 3896
        GS1--GS5  0.3044  0.0158 0.6030     TRUE   39
       GS1--GS10  0.0227 -0.1680 0.2438    FALSE 1441
       GS5--GS10  0.4217  0.0417 0.7450     TRUE   25
```

## p = 100 (4) selected negative edges

Rank of each within every method's own ordering.

*Source: p100_tables_extra.rds*

```
                 VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
IPCONGD--IPBUSEQ    13      33         1370               867                 3              6                  1270                       2
PCDGx--PCESVx       28      15         1370                 8                51            128                  1270                      31
```

## p = 100 (4) group purity of the top 100 pairs

*Source: p100_tables_extra.rds*

```
                  method top_n same_group purity
                   VAPOR   100         93   0.93
                 Kastner   100         77   0.77
            VAR + Glasso   100         73   0.73
       Sample covariance   100         31   0.31
       PrecFact (no VAR)   100         82   0.82
          VAR + PrecFact   100         82   0.82
   VAR + Glasso (rank 3)   100         65   0.65
 VAR + PrecFact (rank 3)   100         74   0.74
```

## p = 100 (4) VAR + Glasso: nonzero off-diagonals and group membership

*Source: extra_tables.rds*

```
 nonzero_off_diagonal within_group share_within    rho
                 1369          294        0.215 0.0304
```

## p = 100 (4) VAR + Glasso: 25 largest pairs VAPOR does not select

*Source: extra_tables.rds*

```
                    pair glasso_pcor VAPOR_mean VAPOR_lo VAPOR_hi same_group
             FPIx--PNFIx       0.613      0.208   -0.027    0.778       TRUE
            BAA--BAA10YM       0.570      0.056   -0.152    0.348       TRUE
      S.P.500--NASDAQCOM       0.530      0.411   -0.024    0.788       TRUE
          USTPU--USTRADE       0.504      0.569   -0.012    0.855       TRUE
       PCECTPI--CPIAUCSL       0.418      0.195   -0.032    0.453       TRUE
 BOGMBASEREALx--TOTRESNS       0.414      0.064   -0.089    0.281       TRUE
      PCEPILFE--CPILFESL       0.406      0.097   -0.077    0.339       TRUE
           PCDGx--RSAFSx       0.372      0.157   -0.133    0.643      FALSE
        TB3MS--GS10TB3Mx      -0.356     -0.022   -0.289    0.240       TRUE
         CE16OV--CIVPART       0.345      0.077   -0.061    0.296       TRUE
         GS10--GS10TB3Mx       0.336      0.009   -0.236    0.273       TRUE
           HOUST--PERMIT       0.325      0.257   -0.047    0.610       TRUE
       USGOVT--CONSUMERx       0.306      0.038   -0.093    0.195      FALSE
        S.P.500--VIXCLSx      -0.301     -0.176   -0.603    0.094       TRUE
           GDPC1--GPDIC1       0.299      0.491   -0.015    0.868       TRUE
        RSAFSx--ISRATIOx      -0.294     -0.190   -0.744    0.114       TRUE
      OILPRICEx--EXCAUSx      -0.290     -0.041   -0.225    0.122      FALSE
        AWHNONAG--M1REAL       0.280      0.110   -0.125    0.464      FALSE
     WPSFD49207--WPSID61       0.274      0.197   -0.044    0.527       TRUE
        IPFUELS--CLAIMSx      -0.250     -0.051   -0.222    0.084      FALSE
       PCECTPI--PCEPILFE       0.249      0.345   -0.022    0.658       TRUE
               GS10--AAA       0.247      0.145   -0.040    0.395       TRUE
      HWIURATIOx--M1REAL      -0.243     -0.055   -0.390    0.164      FALSE
        AWHMAN--AWHNONAG       0.241      0.013   -0.152    0.185       TRUE
          PAYEMS--USPRIV       0.240      0.281   -0.028    0.716       TRUE
```

## p = 100 (4) ranks of the selected rate edges

*Source: p100_tables_extra.rds*

```
                  VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
FEDFUNDS--TB3MS       8     754           10               132                16             56                     3                      74
TB3MS--GS1           23     760            5              3660                15             62                     6                      73
GS5--GS10            25    1753            8                18                30              7                     4                       6
GS5--MORTGAGE30US    32    1619           62              2275              3500           2599                    59                    4057
GS1--GS5             39      97           33              1293                74            112                    34                      48
```

## p = 100 (4) GS5--MORTGAGE30US, every method

*Source: p100_tables_extra.rds*

```
                  VAPOR Kastner VAR + Glasso Sample covariance PrecFact (no VAR) VAR + PrecFact VAR + Glasso (rank 3) VAR + PrecFact (rank 3)
GS5--MORTGAGE30US    32    1619           62              2275              3500           2599                    59                    4057
```

## p = 100 (4) the p = 29 edges in the pooled p = 100 estimate

Cross-dimension table.

*Source: extra_tables.rds*

```
                 edge in_p100 rank_of_4950 p100_mean p100_lo p100_hi selected_at_p100 p29_mean
       PAYEMS--UNRATE    TRUE         2623    -0.012  -0.217   0.176            FALSE   -0.786
 CPIAUCSL--WPSFD49207    TRUE          113     0.119  -0.091   0.398            FALSE    0.721
            GS10--BAA    TRUE          143     0.101  -0.160   0.430            FALSE    0.653
        GDPC1--GPDIC1    TRUE           21     0.491  -0.015   0.868            FALSE    0.649
      FEDFUNDS--TB3MS    TRUE            8     0.639   0.310   0.822             TRUE    0.627
       GDPC1--PCECC96    TRUE          158     0.098  -0.132   0.461            FALSE    0.585
      PCECC96--GPDIC1    TRUE          195    -0.084  -0.411   0.138            FALSE   -0.551
           GS1--TB3MS    TRUE           23     0.467   0.173   0.735             TRUE    0.534
     IPFINAL--IPCONGD    TRUE            5     0.691   0.431   0.845             TRUE    0.532
     IPCONGD--IPBUSEQ    TRUE           13    -0.573  -0.759  -0.286             TRUE   -0.412
        GS10--S.P.500    TRUE          412     0.050  -0.153   0.357            FALSE    0.409
         BAA--S.P.500    TRUE          243    -0.071  -0.344   0.154            FALSE   -0.378
            GS1--GS10    TRUE         1444     0.023  -0.168   0.244            FALSE    0.300
```

## p = 100 (4) NIPA block, per chain

*Source: nipa_p100.rds*

```
 seed GDPC1--GPDIC1 GDPC1--IMPGSC1 GDPC1--EXPGSC1 GPDIC1--IMPGSC1 frac_GDPC1_GPDIC1_gt_0.5
  101        0.5925        -0.4194         0.3783          0.3857                    0.701
  102        0.3560        -0.2196         0.2058          0.2041                    0.337
  103        0.3451        -0.2092         0.2050          0.1861                    0.329
  104        0.6694        -0.4878         0.4325          0.4451                    0.826
```

## p = 100 (4) NIPA block, pooled R-hat and ESS

*Source: nipa_p100.rds*

```
            pair pooled_mean  Rhat bulk_ESS tail_ESS
   GDPC1--GPDIC1      0.4908 1.233       13      420
  GDPC1--IMPGSC1     -0.3340 1.228       13      477
  GDPC1--EXPGSC1      0.3054 1.205       15      624
 GPDIC1--IMPGSC1      0.3052 1.227       13      548
```

## p = 100 (4) NIPA accounting identity, coefficients

*Source: nipa_p100.rds*

```
        term estimate std_error      t        p
 (Intercept)   0.0000    0.0054   0.00 1.00e+00
     PCECC96   0.7510    0.0099  76.01 3.35e-86
      GPDIC1   0.4802    0.0098  49.18 7.09e-69
       GCEC1   0.1140    0.0059  19.24 2.28e-34
     EXPGSC1   0.2947    0.0111  26.66 1.30e-45
     IMPGSC1  -0.3833    0.0152 -25.19 1.38e-43
```

## p = 100 (4) NIPA accounting identity, fit

*Source: nipa_p100.rds*

```
   n df_residual       R2   adj_R2 residual_SD residual_SD_as_share_of_GDPC1_SD F_stat
 100          94 0.997234 0.997087    0.053975                         0.053975 6777.6
```

## p = 100 (4) NIPA block, smallest eigenvalues of the 6 x 6 correlation matrix

Condition number 3282.04.

*Source: nipa_p100.rds*

```
         rank eigenvalue
     smallest 0.00134867
 2nd smallest 0.12687300
 3rd smallest 0.23936600
```

## p = 100 (4) sample correlations among the five rate series

*Source: nipa_p100.rds*

```
   min median   max
 0.276  0.679 0.954
```

## p = 100 (4) partial-correlation convention, per method

*Source: p100_tables_extra.rds*

```
                  method                              estimate chains draws
                   VAPOR mean of per-draw partial correlations      4  4000
                 Kastner mean of per-draw partial correlations      4  4000
            VAR + Glasso                        point estimate      0    NA
       Sample covariance                        point estimate      0    NA
       PrecFact (no VAR) mean of per-draw partial correlations      4  4000
          VAR + PrecFact mean of per-draw partial correlations      4  4000
   VAR + Glasso (rank 3)                        point estimate      0    NA
 VAR + PrecFact (rank 3) mean of per-draw partial correlations      4  4000
```


# Figures

```
diagnostics_p100_ppc.png         2026-10-20 14:18
diagnostics_p100_trace.png       2026-10-20 14:18
diagnostics_p29_ppc.png          2026-10-20 13:35
diagnostics_p29_trace.png        2026-10-20 13:35
heatmap_p100_r3k25.png           2026-10-20 15:24
heatmap_p29_r3k10.png            2026-10-20 15:24
network_p100_grouped.png         2026-10-20 15:24
nipa_trace_p100.png              2026-10-20 15:25
termstructure_p100.png           2026-10-20 15:24
```


# Simulations

Rerun under the DP base measure Gamma(as = 1, bs = 1), for VAPOR and for PrecFact; a_dir = b_dir = 0.1; residual-precision floor 1e-8. Same seeds, replicates and settings as before.
Paper naming: Simulation E is sim_f.R and Simulation F is sim_g.R. **Simulation F has three settings, not four: p = 200, T = 140 ran 17 h to replicate 19 of 25 and was cancelled.**
Every number below is from this run; nothing is carried over from the Gamma(0.1, 0.1) tables.

## Simulation A, Experiment A (varying rho(B), T = 200): every cell

Mean (sd) over 25 replicates.

*Source: sim_A_results.rds*

```
 Scenario            Method          Sens          Spec          Frob
  rho=0.1 PrecFact (no VAR) 0.555 (0.188) 0.998 (0.004) 0.775 (0.093)
  rho=0.1   Glasso (no VAR) 0.825 (0.149) 0.510 (0.138) 0.943 (0.129)
  rho=0.1      VAR + Glasso 0.833 (0.151) 0.512 (0.140) 0.929 (0.113)
  rho=0.1    VAR + PrecFact 0.555 (0.187) 0.998 (0.004) 0.751 (0.077)
  rho=0.1             VAPOR 0.567 (0.154) 0.999 (0.002) 0.684 (0.074)
  rho=0.1 Kastner (VAR+FSV) 0.547 (0.118) 0.934 (0.021) 1.118 (0.397)
  rho=0.5 PrecFact (no VAR) 0.398 (0.213) 0.970 (0.024) 1.524 (0.482)
  rho=0.5   Glasso (no VAR) 0.675 (0.308) 0.590 (0.216) 1.569 (0.380)
  rho=0.5      VAR + Glasso 0.834 (0.143) 0.527 (0.140) 0.932 (0.113)
  rho=0.5    VAR + PrecFact 0.557 (0.184) 0.999 (0.003) 0.755 (0.076)
  rho=0.5             VAPOR 0.566 (0.157) 1.000 (0.001) 0.695 (0.069)
  rho=0.5 Kastner (VAR+FSV) 0.536 (0.120) 0.933 (0.022) 1.199 (0.379)
 rho=0.75 PrecFact (no VAR) 0.304 (0.147) 0.926 (0.036) 2.363 (0.562)
 rho=0.75   Glasso (no VAR) 0.407 (0.363) 0.757 (0.242) 2.031 (0.411)
 rho=0.75      VAR + Glasso 0.837 (0.145) 0.529 (0.140) 0.935 (0.115)
 rho=0.75    VAR + PrecFact 0.561 (0.191) 0.998 (0.004) 0.761 (0.082)
 rho=0.75             VAPOR 0.553 (0.157) 1.000 (0.002) 0.693 (0.067)
 rho=0.75 Kastner (VAR+FSV) 0.536 (0.117) 0.933 (0.020) 1.196 (0.406)
 rho=0.89 PrecFact (no VAR) 0.236 (0.098) 0.888 (0.050) 2.977 (0.604)
 rho=0.89   Glasso (no VAR) 0.146 (0.294) 0.903 (0.225) 2.177 (0.458)
 rho=0.89      VAR + Glasso 0.838 (0.148) 0.527 (0.146) 0.935 (0.113)
 rho=0.89    VAR + PrecFact 0.575 (0.171) 0.998 (0.004) 0.758 (0.077)
 rho=0.89             VAPOR 0.559 (0.160) 1.000 (0.002) 0.694 (0.067)
 rho=0.89 Kastner (VAR+FSV) 0.529 (0.116) 0.932 (0.021) 1.142 (0.396)
 rho=0.95 PrecFact (no VAR) 0.238 (0.111) 0.864 (0.061) 3.269 (0.656)
 rho=0.95   Glasso (no VAR) 0.053 (0.192) 0.967 (0.139) 2.193 (0.468)
 rho=0.95      VAR + Glasso 0.839 (0.144) 0.525 (0.146) 0.934 (0.114)
 rho=0.95    VAR + PrecFact 0.555 (0.181) 0.998 (0.004) 0.762 (0.082)
 rho=0.95             VAPOR 0.550 (0.169) 1.000 (0.002) 0.689 (0.071)
 rho=0.95 Kastner (VAR+FSV) 0.543 (0.109) 0.934 (0.020) 1.178 (0.399)
 rho=0.99 PrecFact (no VAR) 0.239 (0.116) 0.845 (0.068) 3.586 (0.741)
 rho=0.99   Glasso (no VAR) 0.009 (0.018) 0.994 (0.014) 2.197 (0.470)
 rho=0.99      VAR + Glasso 0.830 (0.146) 0.530 (0.141) 0.931 (0.112)
 rho=0.99    VAR + PrecFact 0.569 (0.184) 0.999 (0.004) 0.756 (0.080)
 rho=0.99             VAPOR 0.557 (0.158) 1.000 (0.001) 0.688 (0.068)
 rho=0.99 Kastner (VAR+FSV) 0.547 (0.097) 0.934 (0.019) 1.169 (0.389)
```

## Simulation A, Experiment B (varying T, rho(B) = 0.95): every cell

*Source: sim_A_results.rds*

```
 Scenario            Method          Sens          Spec          Frob
     T=50 PrecFact (no VAR) 0.095 (0.080) 0.950 (0.030) 3.199 (0.624)
     T=50   Glasso (no VAR) 0.022 (0.047) 0.991 (0.019) 2.196 (0.471)
     T=50      VAR + Glasso 0.602 (0.214) 0.716 (0.168) 1.711 (0.359)
     T=50    VAR + PrecFact 0.192 (0.137) 0.994 (0.011) 1.740 (0.305)
     T=50             VAPOR 0.020 (0.039) 1.000 (0.000) 1.652 (0.357)
     T=50 Kastner (VAR+FSV) 0.396 (0.101) 0.913 (0.021) 1.796 (0.444)
    T=100 PrecFact (no VAR) 0.170 (0.099) 0.900 (0.043) 3.249 (0.663)
    T=100   Glasso (no VAR) 0.035 (0.071) 0.984 (0.040) 2.189 (0.473)
    T=100      VAR + Glasso 0.741 (0.153) 0.587 (0.151) 1.321 (0.192)
    T=100    VAR + PrecFact 0.424 (0.187) 0.996 (0.006) 1.143 (0.116)
    T=100             VAPOR 0.372 (0.145) 1.000 (0.001) 1.038 (0.157)
    T=100 Kastner (VAR+FSV) 0.484 (0.113) 0.925 (0.021) 1.386 (0.369)
    T=200 PrecFact (no VAR) 0.238 (0.111) 0.864 (0.061) 3.269 (0.656)
    T=200   Glasso (no VAR) 0.053 (0.192) 0.967 (0.139) 2.193 (0.468)
    T=200      VAR + Glasso 0.839 (0.144) 0.525 (0.146) 0.934 (0.114)
    T=200    VAR + PrecFact 0.555 (0.181) 0.998 (0.004) 0.762 (0.082)
    T=200             VAPOR 0.550 (0.169) 1.000 (0.002) 0.689 (0.071)
    T=200 Kastner (VAR+FSV) 0.543 (0.109) 0.934 (0.020) 1.178 (0.399)
    T=500 PrecFact (no VAR) 0.339 (0.151) 0.810 (0.076) 3.357 (0.743)
    T=500   Glasso (no VAR) 0.106 (0.267) 0.931 (0.190) 2.185 (0.463)
    T=500      VAR + Glasso 0.915 (0.084) 0.476 (0.123) 0.606 (0.071)
    T=500    VAR + PrecFact 0.697 (0.143) 1.000 (0.000) 0.472 (0.056)
    T=500             VAPOR 0.691 (0.131) 1.000 (0.000) 0.439 (0.049)
    T=500 Kastner (VAR+FSV) 0.613 (0.110) 0.944 (0.019) 1.050 (0.443)
   T=1000 PrecFact (no VAR) 0.384 (0.146) 0.782 (0.071) 3.327 (0.706)
   T=1000   Glasso (no VAR) 0.121 (0.296) 0.923 (0.201) 2.182 (0.461)
   T=1000      VAR + Glasso 0.952 (0.044) 0.438 (0.086) 0.453 (0.036)
   T=1000    VAR + PrecFact 0.738 (0.130) 1.000 (0.002) 0.335 (0.046)
   T=1000             VAPOR 0.736 (0.124) 1.000 (0.001) 0.312 (0.037)
   T=1000 Kastner (VAR+FSV) 0.650 (0.113) 0.949 (0.019) 0.999 (0.454)
```

## Simulation B, Experiment A (varying rho(B), T = 200): every cell

*Source: sim_B_results.rds*

```
 Scenario            Method          Sens          Spec          Frob
  rho=0.1 PrecFact (no VAR) 0.633 (0.120) 0.985 (0.010) 1.069 (0.159)
  rho=0.1   Glasso (no VAR) 0.994 (0.018) 0.801 (0.060) 0.633 (0.100)
  rho=0.1      VAR + Glasso 0.994 (0.018) 0.801 (0.056) 0.633 (0.100)
  rho=0.1    VAR + PrecFact 0.635 (0.129) 0.985 (0.012) 1.068 (0.179)
  rho=0.1             VAPOR 0.613 (0.119) 0.992 (0.008) 1.027 (0.180)
  rho=0.1 Kastner (VAR+FSV) 0.766 (0.120) 0.977 (0.015) 0.894 (0.217)
  rho=0.5 PrecFact (no VAR) 0.454 (0.129) 0.978 (0.016) 1.337 (0.243)
  rho=0.5   Glasso (no VAR) 0.930 (0.060) 0.786 (0.063) 1.025 (0.216)
  rho=0.5      VAR + Glasso 0.991 (0.022) 0.796 (0.060) 0.638 (0.102)
  rho=0.5    VAR + PrecFact 0.636 (0.120) 0.984 (0.012) 1.082 (0.183)
  rho=0.5             VAPOR 0.581 (0.141) 0.991 (0.009) 1.027 (0.176)
  rho=0.5 Kastner (VAR+FSV) 0.769 (0.100) 0.977 (0.013) 0.882 (0.213)
 rho=0.75 PrecFact (no VAR) 0.352 (0.129) 0.949 (0.026) 1.853 (0.337)
 rho=0.75   Glasso (no VAR) 0.652 (0.328) 0.851 (0.093) 1.471 (0.244)
 rho=0.75      VAR + Glasso 0.994 (0.018) 0.806 (0.056) 0.639 (0.099)
 rho=0.75    VAR + PrecFact 0.649 (0.130) 0.984 (0.012) 1.090 (0.183)
 rho=0.75             VAPOR 0.592 (0.139) 0.992 (0.008) 1.021 (0.168)
 rho=0.75 Kastner (VAR+FSV) 0.767 (0.093) 0.977 (0.013) 0.905 (0.191)
 rho=0.89 PrecFact (no VAR) 0.298 (0.126) 0.925 (0.033) 2.308 (0.443)
 rho=0.89   Glasso (no VAR) 0.296 (0.359) 0.931 (0.099) 1.649 (0.188)
 rho=0.89      VAR + Glasso 0.994 (0.018) 0.805 (0.060) 0.639 (0.099)
 rho=0.89    VAR + PrecFact 0.652 (0.132) 0.982 (0.012) 1.082 (0.172)
 rho=0.89             VAPOR 0.585 (0.126) 0.991 (0.009) 1.041 (0.183)
 rho=0.89 Kastner (VAR+FSV) 0.765 (0.114) 0.977 (0.014) 0.898 (0.219)
 rho=0.95 PrecFact (no VAR) 0.263 (0.132) 0.909 (0.033) 2.610 (0.520)
 rho=0.95   Glasso (no VAR) 0.091 (0.246) 0.982 (0.053) 1.682 (0.168)
 rho=0.95      VAR + Glasso 0.994 (0.018) 0.803 (0.059) 0.638 (0.099)
 rho=0.95    VAR + PrecFact 0.630 (0.157) 0.982 (0.015) 1.085 (0.169)
 rho=0.95             VAPOR 0.588 (0.111) 0.992 (0.008) 1.030 (0.172)
 rho=0.95 Kastner (VAR+FSV) 0.774 (0.110) 0.977 (0.014) 0.899 (0.229)
 rho=0.99 PrecFact (no VAR) 0.305 (0.120) 0.877 (0.062) 2.909 (0.679)
 rho=0.99   Glasso (no VAR) 0.087 (0.246) 0.983 (0.050) 1.693 (0.153)
 rho=0.99      VAR + Glasso 0.994 (0.018) 0.802 (0.058) 0.636 (0.100)
 rho=0.99    VAR + PrecFact 0.629 (0.103) 0.983 (0.012) 1.088 (0.180)
 rho=0.99             VAPOR 0.575 (0.150) 0.992 (0.007) 1.037 (0.163)
 rho=0.99 Kastner (VAR+FSV) 0.765 (0.103) 0.977 (0.013) 0.888 (0.210)
```

## Simulation B, Experiment B (varying T, rho(B) = 0.95): every cell

*Source: sim_B_results.rds*

```
 Scenario            Method          Sens          Spec          Frob
     T=50 PrecFact (no VAR) 0.082 (0.086) 0.965 (0.030) 2.534 (0.546)
     T=50   Glasso (no VAR) 0.094 (0.145) 0.976 (0.040) 1.685 (0.147)
     T=50      VAR + Glasso 0.767 (0.135) 0.838 (0.045) 1.284 (0.115)
     T=50    VAR + PrecFact 0.124 (0.101) 0.994 (0.008) 1.770 (0.216)
     T=50             VAPOR 0.000 (0.000) 1.000 (0.000) 1.537 (0.107)
     T=50 Kastner (VAR+FSV) 0.419 (0.092) 0.946 (0.014) 1.572 (0.182)
    T=100 PrecFact (no VAR) 0.166 (0.085) 0.929 (0.044) 2.571 (0.525)
    T=100   Glasso (no VAR) 0.110 (0.247) 0.971 (0.077) 1.688 (0.146)
    T=100      VAR + Glasso 0.948 (0.052) 0.804 (0.055) 0.925 (0.135)
    T=100    VAR + PrecFact 0.313 (0.134) 0.992 (0.014) 1.362 (0.176)
    T=100             VAPOR 0.162 (0.106) 0.998 (0.006) 1.290 (0.156)
    T=100 Kastner (VAR+FSV) 0.620 (0.089) 0.964 (0.015) 1.152 (0.183)
    T=200 PrecFact (no VAR) 0.263 (0.132) 0.909 (0.033) 2.610 (0.520)
    T=200   Glasso (no VAR) 0.091 (0.246) 0.982 (0.053) 1.682 (0.168)
    T=200      VAR + Glasso 0.994 (0.018) 0.803 (0.059) 0.638 (0.099)
    T=200    VAR + PrecFact 0.630 (0.157) 0.982 (0.015) 1.085 (0.169)
    T=200             VAPOR 0.588 (0.111) 0.992 (0.008) 1.030 (0.172)
    T=200 Kastner (VAR+FSV) 0.774 (0.110) 0.977 (0.014) 0.899 (0.229)
    T=500 PrecFact (no VAR) 0.383 (0.146) 0.857 (0.073) 2.593 (0.587)
    T=500   Glasso (no VAR) 0.142 (0.310) 0.975 (0.060) 1.675 (0.170)
    T=500      VAR + Glasso 1.000 (0.000) 0.801 (0.055) 0.405 (0.039)
    T=500    VAR + PrecFact 0.746 (0.114) 0.977 (0.011) 0.940 (0.170)
    T=500             VAPOR 0.752 (0.106) 0.981 (0.013) 0.913 (0.177)
    T=500 Kastner (VAR+FSV) 0.838 (0.096) 0.983 (0.012) 0.731 (0.227)
   T=1000 PrecFact (no VAR) 0.450 (0.137) 0.833 (0.072) 2.623 (0.548)
   T=1000   Glasso (no VAR) 0.149 (0.338) 0.970 (0.069) 1.682 (0.162)
   T=1000      VAR + Glasso 1.000 (0.000) 0.788 (0.055) 0.291 (0.039)
   T=1000    VAR + PrecFact 0.777 (0.097) 0.972 (0.014) 0.898 (0.207)
   T=1000             VAPOR 0.789 (0.111) 0.975 (0.015) 0.857 (0.200)
   T=1000 Kastner (VAR+FSV) 0.839 (0.090) 0.984 (0.012) 0.687 (0.237)
```

## Simulation C (factor stochastic volatility DGP): every cell

Held-out prediction MSPE, mean (sd) over 25 replicates.

*Source: sim_C_results.rds*

```
           Scenario            Method            MSPE
        Sparse (1%) PrecFact (no VAR) 0.9093 (0.1250)
        Sparse (1%)    VAR + PrecFact 0.9041 (0.1295)
        Sparse (1%)             VAPOR 0.9097 (0.1264)
        Sparse (1%) Kastner (VAR+FSV) 0.8520 (0.1492)
        Sparse (1%) Sample covariance 0.9190 (0.1543)
 Intermediate (10%) PrecFact (no VAR) 0.9007 (0.1437)
 Intermediate (10%)    VAR + PrecFact 0.8992 (0.1325)
 Intermediate (10%)             VAPOR 0.9085 (0.1413)
 Intermediate (10%) Kastner (VAR+FSV) 0.8403 (0.1283)
 Intermediate (10%) Sample covariance 0.8985 (0.1474)
        Dense (80%) PrecFact (no VAR) 0.7066 (0.2062)
        Dense (80%)    VAR + PrecFact 0.8531 (0.3504)
        Dense (80%)             VAPOR 0.8527 (0.3393)
        Dense (80%) Kastner (VAR+FSV) 0.7427 (0.3266)
        Dense (80%) Sample covariance 0.4668 (0.1551)
```

## Simulation D (varying signal strength): every cell

*Source: sim_D_results.rds*

```
 Scenario            Method          Sens          Spec          Frob
  val=0.5 PrecFact (no VAR) 0.123 (0.062) 0.888 (0.048) 2.651 (0.708)
  val=0.5   Glasso (no VAR) 0.003 (0.011) 0.999 (0.003) 0.508 (0.127)
  val=0.5      VAR + Glasso 0.154 (0.129) 0.956 (0.038) 0.471 (0.122)
  val=0.5    VAR + PrecFact 0.010 (0.027) 1.000 (0.002) 0.512 (0.124)
  val=0.5             VAPOR 0.015 (0.034) 1.000 (0.000) 0.506 (0.109)
  val=0.5 Kastner (VAR+FSV) 0.270 (0.091) 0.896 (0.016) 0.569 (0.112)
    val=1 PrecFact (no VAR) 0.166 (0.073) 0.884 (0.053) 2.873 (0.689)
    val=1   Glasso (no VAR) 0.006 (0.015) 0.997 (0.007) 1.417 (0.327)
    val=1      VAR + Glasso 0.653 (0.193) 0.706 (0.119) 0.818 (0.135)
    val=1    VAR + PrecFact 0.329 (0.185) 0.998 (0.004) 0.763 (0.102)
    val=1             VAPOR 0.351 (0.177) 0.998 (0.003) 0.681 (0.077)
    val=1 Kastner (VAR+FSV) 0.461 (0.112) 0.923 (0.020) 0.888 (0.188)
  val=1.5 PrecFact (no VAR) 0.238 (0.111) 0.864 (0.061) 3.269 (0.656)
  val=1.5   Glasso (no VAR) 0.053 (0.192) 0.967 (0.139) 2.193 (0.468)
  val=1.5      VAR + Glasso 0.839 (0.144) 0.525 (0.146) 0.934 (0.114)
  val=1.5    VAR + PrecFact 0.555 (0.181) 0.998 (0.004) 0.762 (0.082)
  val=1.5             VAPOR 0.550 (0.169) 1.000 (0.002) 0.689 (0.071)
  val=1.5 Kastner (VAR+FSV) 0.543 (0.109) 0.934 (0.020) 1.178 (0.399)
  val=2.5 PrecFact (no VAR) 0.333 (0.147) 0.849 (0.049) 3.882 (0.777)
  val=2.5   Glasso (no VAR) 0.307 (0.404) 0.779 (0.298) 3.115 (0.569)
  val=2.5      VAR + Glasso 0.942 (0.067) 0.360 (0.105) 0.968 (0.110)
  val=2.5    VAR + PrecFact 0.751 (0.101) 0.999 (0.002) 0.702 (0.066)
  val=2.5             VAPOR 0.718 (0.114) 1.000 (0.001) 0.724 (0.113)
  val=2.5 Kastner (VAR+FSV) 0.599 (0.120) 0.942 (0.021) 1.830 (0.691)
    val=4 PrecFact (no VAR) 0.400 (0.146) 0.835 (0.060) 4.387 (0.859)
    val=4   Glasso (no VAR) 0.435 (0.439) 0.673 (0.343) 3.769 (0.621)
    val=4      VAR + Glasso 0.965 (0.033) 0.297 (0.045) 0.970 (0.105)
    val=4    VAR + PrecFact 0.829 (0.090) 0.997 (0.006) 0.680 (0.080)
    val=4             VAPOR 0.778 (0.102) 0.999 (0.002) 0.887 (0.137)
    val=4 Kastner (VAR+FSV) 0.620 (0.117) 0.945 (0.020) 2.453 (0.774)
```

## Simulation E (sim_f.R, approximately low-rank B, T >= p): every cell

Omega_Frob is ||R_bar - R_0||_F; B_Frob is ||B_hat - B_0||_F/||B_0||_F (VAPOR only); Sens and Spec are at an oracle threshold retaining the true number of edges.

*Source: sim_f_final_results.rds*

```
   setting          Method Omega_Frob Omega_Frob_sd B_Frob B_Frob_sd   Sens Sens_sd   Spec Spec_sd
   (50,50)    VAPOR dl r=3     3.9512        0.4136 0.6981    0.0991 0.2972  0.0795 0.9454  0.0110
   (50,50) VAPOR ridge r=3     3.9559        0.4025 1.0945    0.5356 0.2926  0.0680 0.9450  0.0103
   (50,50)         Kastner     4.6015        0.4570    NaN        NA 0.1321  0.0421 0.9329  0.0090
   (50,50)          Glasso     4.3458        0.4622    NaN        NA 0.9309  0.2391 0.0747  0.2584
   (50,75)    VAPOR dl r=3     3.2543        0.5016 0.5656    0.1452 0.4318  0.0922 0.9560  0.0095
   (50,75) VAPOR ridge r=3     3.2378        0.4863 0.5825    0.3046 0.4414  0.0960 0.9567  0.0098
   (50,75)         Kastner     4.4944        0.4778    NaN        NA 0.1416  0.0354 0.9336  0.0087
   (50,75)          Glasso     4.3243        0.4682    NaN        NA 0.8105  0.3445 0.2253  0.4092
  (50,100)    VAPOR dl r=3     2.6882        0.4332 0.4776    0.1549 0.5550  0.0760 0.9656  0.0076
  (50,100) VAPOR ridge r=3     2.6343        0.4367 0.4370    0.2101 0.5664  0.0813 0.9665  0.0077
  (50,100)         Kastner     4.4752        0.4828    NaN        NA 0.1471  0.0421 0.9340  0.0091
  (50,100)          Glasso     4.2611        0.4585    NaN        NA 0.6235  0.3718 0.4905  0.4810
 (100,100)    VAPOR dl r=3     6.8470        1.2962 0.6398    0.0753 0.4678  0.0850 0.9508  0.0146
 (100,100) VAPOR ridge r=3     6.9339        1.2157 0.9365    0.3273 0.4652  0.0825 0.9505  0.0147
 (100,100)         Kastner     9.5971        1.0114    NaN        NA 0.1079  0.0244 0.9189  0.0134
 (100,100)          Glasso     9.4728        0.9987    NaN        NA 0.9008  0.2743 0.1093  0.3022
 (150,150)    VAPOR dl r=3     8.0995        1.2978 0.5695    0.0728 0.5800  0.0566 0.9611  0.0083
 (150,150) VAPOR ridge r=3     8.1196        1.3674 0.7773    0.2539 0.5778  0.0561 0.9609  0.0083
 (150,150)         Kastner    14.5657        1.4236    NaN        NA 0.0949  0.0181 0.9167  0.0093
 (150,150)          Glasso    14.4578        1.4096    NaN        NA 0.7274  0.4056 0.2927  0.4355
```

## Simulation F (sim_g.R, T < p at T/p = 0.7): every cell

Three settings; p = 200, T = 140 was cancelled.

*Source: sim_g_results.rds*

```
   setting          Method Omega_Frob Omega_Frob_sd B_Frob B_Frob_sd   Sens Sens_sd   Spec Spec_sd
   (50,35)    VAPOR dl r=3     4.3591        0.3613 0.8360    0.0765 0.1971  0.0575 0.9379  0.0094
   (50,35) VAPOR ridge r=3     4.3217        0.3838 0.8984    0.2689 0.1971  0.0591 0.9379  0.0095
   (50,35)         Kastner     4.7484        0.4802    NaN        NA 0.1179  0.0408 0.9319  0.0083
   (50,35)          Glasso     4.3490        0.4557    NaN        NA 0.9646  0.1770 0.0373  0.1865
  (100,70)    VAPOR dl r=3     8.9510        1.0524 0.7932    0.0544 0.2341  0.0444 0.9298  0.0153
  (100,70) VAPOR ridge r=3     8.9771        1.0630 0.7945    0.1540 0.2295  0.0456 0.9294  0.0153
  (100,70)         Kastner     9.6320        0.9973    NaN        NA 0.0940  0.0222 0.9177  0.0134
  (100,70)          Glasso     9.4770        1.0062    NaN        NA 1.0000  0.0000 0.0000  0.0000
 (150,105)    VAPOR dl r=3    13.4741        1.5485 0.7666    0.0491 0.2361  0.0503 0.9294  0.0111
 (150,105) VAPOR ridge r=3    13.5582        1.5374 0.7151    0.1110 0.2298  0.0475 0.9288  0.0110
 (150,105)         Kastner    14.5832        1.4414    NaN        NA 0.0939  0.0140 0.9166  0.0096
 (150,105)          Glasso    14.4836        1.4338    NaN        NA 1.0000  0.0000 0.0000  0.0000
```

## Simulation A, Experiment A: paired differences in Frobenius distance

Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is closer to the truth.

*Source: sim_A_results.rds*

```
 Scenario            Method  Diff    SE     t      p    Wins
  rho=0.1 PrecFact (no VAR) 0.091 0.010  8.81 <1e-04 23 / 25
  rho=0.1   Glasso (no VAR) 0.259 0.021 12.50 <1e-04 25 / 25
  rho=0.1      VAR + Glasso 0.246 0.019 12.77 <1e-04 25 / 25
  rho=0.1    VAR + PrecFact 0.067 0.010  6.48 <1e-04 22 / 25
  rho=0.1 Kastner (VAR+FSV) 0.435 0.076  5.74 <1e-04 23 / 25
  rho=0.5 PrecFact (no VAR) 0.830 0.089  9.36 <1e-04 25 / 25
  rho=0.5   Glasso (no VAR) 0.874 0.069 12.75 <1e-04 25 / 25
  rho=0.5      VAR + Glasso 0.237 0.018 13.52 <1e-04 25 / 25
  rho=0.5    VAR + PrecFact 0.060 0.012  4.92 <1e-04 20 / 25
  rho=0.5 Kastner (VAR+FSV) 0.504 0.072  7.02 <1e-04 25 / 25
 rho=0.75 PrecFact (no VAR) 1.671 0.107 15.56 <1e-04 25 / 25
 rho=0.75   Glasso (no VAR) 1.338 0.076 17.56 <1e-04 25 / 25
 rho=0.75      VAR + Glasso 0.242 0.017 14.16 <1e-04 25 / 25
 rho=0.75    VAR + PrecFact 0.068 0.012  5.54 <1e-04 18 / 25
 rho=0.75 Kastner (VAR+FSV) 0.504 0.077  6.51 <1e-04 25 / 25
 rho=0.89 PrecFact (no VAR) 2.284 0.117 19.59 <1e-04 25 / 25
 rho=0.89   Glasso (no VAR) 1.483 0.088 16.89 <1e-04 25 / 25
 rho=0.89      VAR + Glasso 0.241 0.018 13.19 <1e-04 25 / 25
 rho=0.89    VAR + PrecFact 0.064 0.012  5.45 <1e-04 20 / 25
 rho=0.89 Kastner (VAR+FSV) 0.449 0.075  5.96 <1e-04 24 / 25
 rho=0.95 PrecFact (no VAR) 2.580 0.126 20.44 <1e-04 25 / 25
 rho=0.95   Glasso (no VAR) 1.505 0.091 16.58 <1e-04 25 / 25
 rho=0.95      VAR + Glasso 0.246 0.019 12.93 <1e-04 25 / 25
 rho=0.95    VAR + PrecFact 0.073 0.012  6.18 <1e-04 21 / 25
 rho=0.95 Kastner (VAR+FSV) 0.489 0.076  6.47 <1e-04 24 / 25
 rho=0.99 PrecFact (no VAR) 2.898 0.145 20.01 <1e-04 25 / 25
 rho=0.99   Glasso (no VAR) 1.509 0.091 16.55 <1e-04 25 / 25
 rho=0.99      VAR + Glasso 0.243 0.019 12.54 <1e-04 25 / 25
 rho=0.99    VAR + PrecFact 0.067 0.011  6.01 <1e-04 19 / 25
 rho=0.99 Kastner (VAR+FSV) 0.481 0.072  6.65 <1e-04 24 / 25
```

## Simulation A, Experiment B: paired differences in Frobenius distance

Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is closer to the truth.

*Source: sim_A_results.rds*

```
 Scenario            Method  Diff    SE     t       p    Wins
     T=50 PrecFact (no VAR) 1.547 0.116 13.29  <1e-04 25 / 25
     T=50   Glasso (no VAR) 0.544 0.036 14.96  <1e-04 25 / 25
     T=50      VAR + Glasso 0.059 0.024  2.44   0.022 17 / 25
     T=50    VAR + PrecFact 0.088 0.064  1.37    0.18 18 / 25
     T=50 Kastner (VAR+FSV) 0.144 0.033  4.33 0.00023 21 / 25
    T=100 PrecFact (no VAR) 2.211 0.130 17.02  <1e-04 25 / 25
    T=100   Glasso (no VAR) 1.151 0.089 13.00  <1e-04 25 / 25
    T=100      VAR + Glasso 0.283 0.030  9.58  <1e-04 23 / 25
    T=100    VAR + PrecFact 0.105 0.021  5.06  <1e-04 22 / 25
    T=100 Kastner (VAR+FSV) 0.348 0.067  5.17  <1e-04 22 / 25
    T=200 PrecFact (no VAR) 2.580 0.126 20.44  <1e-04 25 / 25
    T=200   Glasso (no VAR) 1.505 0.091 16.58  <1e-04 25 / 25
    T=200      VAR + Glasso 0.246 0.019 12.93  <1e-04 25 / 25
    T=200    VAR + PrecFact 0.073 0.012  6.18  <1e-04 21 / 25
    T=200 Kastner (VAR+FSV) 0.489 0.076  6.47  <1e-04 24 / 25
    T=500 PrecFact (no VAR) 2.918 0.150 19.50  <1e-04 25 / 25
    T=500   Glasso (no VAR) 1.746 0.092 18.89  <1e-04 25 / 25
    T=500      VAR + Glasso 0.168 0.010 16.49  <1e-04 25 / 25
    T=500    VAR + PrecFact 0.033 0.006  5.82  <1e-04 23 / 25
    T=500 Kastner (VAR+FSV) 0.612 0.089  6.89  <1e-04 24 / 25
   T=1000 PrecFact (no VAR) 3.015 0.141 21.33  <1e-04 25 / 25
   T=1000   Glasso (no VAR) 1.870 0.092 20.29  <1e-04 25 / 25
   T=1000      VAR + Glasso 0.141 0.007 20.23  <1e-04 25 / 25
   T=1000    VAR + PrecFact 0.023 0.004  5.14  <1e-04 22 / 25
   T=1000 Kastner (VAR+FSV) 0.687 0.090  7.62  <1e-04 25 / 25
```

## Simulation B, Experiment A: paired differences in Frobenius distance

Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is closer to the truth.

*Source: sim_B_results.rds*

```
 Scenario            Method   Diff    SE      t       p    Wins
  rho=0.1 PrecFact (no VAR)  0.042 0.012   3.57  0.0016 18 / 25
  rho=0.1   Glasso (no VAR) -0.393 0.030 -12.98  <1e-04  0 / 25
  rho=0.1      VAR + Glasso -0.393 0.030 -13.21  <1e-04  0 / 25
  rho=0.1    VAR + PrecFact  0.042 0.009   4.62 0.00011 21 / 25
  rho=0.1 Kastner (VAR+FSV) -0.133 0.021  -6.38  <1e-04  1 / 25
  rho=0.5 PrecFact (no VAR)  0.311 0.040   7.76  <1e-04 25 / 25
  rho=0.5   Glasso (no VAR) -0.002 0.042  -0.04    0.97 11 / 25
  rho=0.5      VAR + Glasso -0.388 0.029 -13.33  <1e-04  0 / 25
  rho=0.5    VAR + PrecFact  0.055 0.009   6.36  <1e-04 22 / 25
  rho=0.5 Kastner (VAR+FSV) -0.145 0.022  -6.57  <1e-04  4 / 25
 rho=0.75 PrecFact (no VAR)  0.832 0.063  13.26  <1e-04 25 / 25
 rho=0.75   Glasso (no VAR)  0.450 0.047   9.61  <1e-04 25 / 25
 rho=0.75      VAR + Glasso -0.382 0.029 -13.18  <1e-04  0 / 25
 rho=0.75    VAR + PrecFact  0.069 0.011   6.30  <1e-04 23 / 25
 rho=0.75 Kastner (VAR+FSV) -0.116 0.016  -7.14  <1e-04  2 / 25
 rho=0.89 PrecFact (no VAR)  1.268 0.086  14.73  <1e-04 25 / 25
 rho=0.89   Glasso (no VAR)  0.609 0.030  20.21  <1e-04 25 / 25
 rho=0.89      VAR + Glasso -0.402 0.032 -12.53  <1e-04  0 / 25
 rho=0.89    VAR + PrecFact  0.041 0.011   3.57  0.0016 21 / 25
 rho=0.89 Kastner (VAR+FSV) -0.143 0.022  -6.61  <1e-04  1 / 25
 rho=0.95 PrecFact (no VAR)  1.581 0.099  15.90  <1e-04 25 / 25
 rho=0.95   Glasso (no VAR)  0.653 0.024  27.09  <1e-04 25 / 25
 rho=0.95      VAR + Glasso -0.392 0.029 -13.35  <1e-04  0 / 25
 rho=0.95    VAR + PrecFact  0.055 0.008   6.68  <1e-04 22 / 25
 rho=0.95 Kastner (VAR+FSV) -0.131 0.022  -5.82  <1e-04  4 / 25
 rho=0.99 PrecFact (no VAR)  1.872 0.137  13.66  <1e-04 25 / 25
 rho=0.99   Glasso (no VAR)  0.656 0.023  28.07  <1e-04 25 / 25
 rho=0.99      VAR + Glasso -0.400 0.027 -14.57  <1e-04  0 / 25
 rho=0.99    VAR + PrecFact  0.051 0.013   4.03 0.00049 21 / 25
 rho=0.99 Kastner (VAR+FSV) -0.149 0.018  -8.34  <1e-04  1 / 25
```

## Simulation B, Experiment B: paired differences in Frobenius distance

Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is closer to the truth.

*Source: sim_B_results.rds*

```
 Scenario            Method   Diff    SE      t      p    Wins
     T=50 PrecFact (no VAR)  0.997 0.103   9.71 <1e-04 25 / 25
     T=50   Glasso (no VAR)  0.148 0.021   7.17 <1e-04 23 / 25
     T=50      VAR + Glasso -0.253 0.020 -12.41 <1e-04  0 / 25
     T=50    VAR + PrecFact  0.232 0.039   5.96 <1e-04 22 / 25
     T=50 Kastner (VAR+FSV)  0.035 0.034   1.03   0.31 14 / 25
    T=100 PrecFact (no VAR)  1.281 0.105  12.22 <1e-04 25 / 25
    T=100   Glasso (no VAR)  0.398 0.025  15.88 <1e-04 25 / 25
    T=100      VAR + Glasso -0.366 0.028 -13.28 <1e-04  0 / 25
    T=100    VAR + PrecFact  0.072 0.019   3.74  0.001 22 / 25
    T=100 Kastner (VAR+FSV) -0.138 0.028  -5.02 <1e-04  5 / 25
    T=200 PrecFact (no VAR)  1.581 0.099  15.90 <1e-04 25 / 25
    T=200   Glasso (no VAR)  0.653 0.024  27.09 <1e-04 25 / 25
    T=200      VAR + Glasso -0.392 0.029 -13.35 <1e-04  0 / 25
    T=200    VAR + PrecFact  0.055 0.008   6.68 <1e-04 22 / 25
    T=200 Kastner (VAR+FSV) -0.131 0.022  -5.82 <1e-04  4 / 25
    T=500 PrecFact (no VAR)  1.681 0.123  13.69 <1e-04 25 / 25
    T=500   Glasso (no VAR)  0.762 0.031  24.32 <1e-04 25 / 25
    T=500      VAR + Glasso -0.508 0.031 -16.34 <1e-04  0 / 25
    T=500    VAR + PrecFact  0.028 0.014   2.01  0.056 18 / 25
    T=500 Kastner (VAR+FSV) -0.181 0.022  -8.28 <1e-04  1 / 25
   T=1000 PrecFact (no VAR)  1.766 0.116  15.29 <1e-04 25 / 25
   T=1000   Glasso (no VAR)  0.825 0.031  26.88 <1e-04 25 / 25
   T=1000      VAR + Glasso -0.566 0.040 -14.15 <1e-04  0 / 25
   T=1000    VAR + PrecFact  0.041 0.019   2.18   0.04 16 / 25
   T=1000 Kastner (VAR+FSV) -0.170 0.023  -7.32 <1e-04  1 / 25
```

## Simulation C: paired differences in held-out MSPE

Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is predicting more accurately.

*Source: sim_C_results.rds*

```
           Scenario            Method    Diff     SE     t       p    Wins
        Sparse (1%) PrecFact (no VAR) -0.0004 0.0037 -0.10    0.92 11 / 25
        Sparse (1%)    VAR + PrecFact -0.0055 0.0025 -2.23   0.035  9 / 25
        Sparse (1%) Kastner (VAR+FSV) -0.0577 0.0104 -5.57  <1e-04  2 / 25
        Sparse (1%) Sample covariance  0.0093 0.0141  0.66    0.52 15 / 25
 Intermediate (10%) PrecFact (no VAR) -0.0078 0.0031 -2.53   0.019  9 / 25
 Intermediate (10%)    VAR + PrecFact -0.0094 0.0062 -1.50    0.15  9 / 25
 Intermediate (10%) Kastner (VAR+FSV) -0.0683 0.0171 -3.99 0.00053  1 / 25
 Intermediate (10%) Sample covariance -0.0100 0.0212 -0.47    0.64 13 / 25
        Dense (80%) PrecFact (no VAR) -0.1462 0.0330 -4.43 0.00018  3 / 25
        Dense (80%)    VAR + PrecFact  0.0003 0.0071  0.04    0.96 11 / 25
        Dense (80%) Kastner (VAR+FSV) -0.1101 0.0176 -6.25  <1e-04  1 / 25
        Dense (80%) Sample covariance -0.3859 0.0690 -5.59  <1e-04  0 / 25
```

## Simulation D: paired differences in Frobenius distance

Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is closer to the truth.

*Source: sim_D_results.rds*

```
 Scenario            Method   Diff    SE     t       p    Wins
  val=0.5 PrecFact (no VAR)  2.146 0.145 14.84  <1e-04 25 / 25
  val=0.5   Glasso (no VAR)  0.002 0.015  0.15    0.88 13 / 25
  val=0.5      VAR + Glasso -0.035 0.013 -2.73   0.012  7 / 25
  val=0.5    VAR + PrecFact  0.006 0.009  0.71    0.49 12 / 25
  val=0.5 Kastner (VAR+FSV)  0.063 0.015  4.09 0.00042 19 / 25
    val=1 PrecFact (no VAR)  2.192 0.134 16.30  <1e-04 25 / 25
    val=1   Glasso (no VAR)  0.736 0.061 12.14  <1e-04 25 / 25
    val=1      VAR + Glasso  0.136 0.022  6.19  <1e-04 21 / 25
    val=1    VAR + PrecFact  0.082 0.011  7.15  <1e-04 23 / 25
    val=1 Kastner (VAR+FSV)  0.207 0.032  6.44  <1e-04 22 / 25
  val=1.5 PrecFact (no VAR)  2.580 0.126 20.44  <1e-04 25 / 25
  val=1.5   Glasso (no VAR)  1.505 0.091 16.58  <1e-04 25 / 25
  val=1.5      VAR + Glasso  0.246 0.019 12.93  <1e-04 25 / 25
  val=1.5    VAR + PrecFact  0.073 0.012  6.18  <1e-04 21 / 25
  val=1.5 Kastner (VAR+FSV)  0.489 0.076  6.47  <1e-04 24 / 25
  val=2.5 PrecFact (no VAR)  3.158 0.146 21.63  <1e-04 25 / 25
  val=2.5   Glasso (no VAR)  2.391 0.104 23.06  <1e-04 25 / 25
  val=2.5      VAR + Glasso  0.244 0.013 18.19  <1e-04 25 / 25
  val=2.5    VAR + PrecFact -0.022 0.019 -1.15    0.26 12 / 25
  val=2.5 Kastner (VAR+FSV)  1.105 0.130  8.52  <1e-04 25 / 25
    val=4 PrecFact (no VAR)  3.500 0.167 20.90  <1e-04 25 / 25
    val=4   Glasso (no VAR)  2.882 0.117 24.63  <1e-04 25 / 25
    val=4      VAR + Glasso  0.083 0.023  3.59  0.0015 19 / 25
    val=4    VAR + PrecFact -0.207 0.032 -6.46  <1e-04  1 / 25
    val=4 Kastner (VAR+FSV)  1.566 0.148 10.57  <1e-04 24 / 25
```

## Simulation E: paired differences against VAPOR (DL)

Reference is VAPOR dl r=3. The 'VAPOR ridge r=3' row IS the DL-vs-ridge difference: positive favours DL.

*Source: sim_f_final_results.rds*

```
     Scenario          Method   Diff    SE     t      p    Wins
   p=50, T=50 VAPOR ridge r=3  0.005 0.015  0.32   0.75 11 / 25
   p=50, T=50         Kastner  0.650 0.048 13.54 <1e-04 25 / 25
   p=50, T=50          Glasso  0.395 0.051  7.68 <1e-04 25 / 25
   p=50, T=75 VAPOR ridge r=3 -0.016 0.018 -0.89   0.38  7 / 25
   p=50, T=75         Kastner  1.240 0.075 16.56 <1e-04 25 / 25
   p=50, T=75          Glasso  1.070 0.082 13.04 <1e-04 25 / 25
  p=50, T=100 VAPOR ridge r=3 -0.054 0.018 -2.95  0.007  6 / 25
  p=50, T=100         Kastner  1.787 0.071 25.20 <1e-04 25 / 25
  p=50, T=100          Glasso  1.573 0.069 22.75 <1e-04 25 / 25
 p=100, T=100 VAPOR ridge r=3  0.087 0.057  1.54   0.14 16 / 25
 p=100, T=100         Kastner  2.750 0.149 18.43 <1e-04 25 / 25
 p=100, T=100          Glasso  2.626 0.152 17.28 <1e-04 25 / 25
 p=150, T=150 VAPOR ridge r=3  0.020 0.090  0.22   0.83 17 / 25
 p=150, T=150         Kastner  6.466 0.242 26.72 <1e-04 25 / 25
 p=150, T=150          Glasso  6.358 0.237 26.81 <1e-04 25 / 25
```

## Simulation F: paired differences against VAPOR (DL)

Reference is VAPOR dl r=3; the 'VAPOR ridge r=3' row is the DL-vs-ridge difference.

*Source: sim_g_results.rds*

```
     Scenario          Method   Diff    SE     t      p    Wins
   p=50, T=35 VAPOR ridge r=3 -0.037 0.011 -3.41 0.0023  7 / 25
   p=50, T=35         Kastner  0.389 0.047  8.34 <1e-04 23 / 25
   p=50, T=35          Glasso -0.010 0.031 -0.32   0.75 13 / 25
  p=100, T=70 VAPOR ridge r=3  0.026 0.020  1.32    0.2 17 / 25
  p=100, T=70         Kastner  0.681 0.043 15.82 <1e-04 25 / 25
  p=100, T=70          Glasso  0.526 0.046 11.42 <1e-04 25 / 25
 p=150, T=105 VAPOR ridge r=3  0.084 0.046  1.83   0.08 18 / 25
 p=150, T=105         Kastner  1.109 0.097 11.38 <1e-04 25 / 25
 p=150, T=105          Glasso  1.009 0.097 10.38 <1e-04 25 / 25
```

## Simulation E: relative coefficient error of VAPOR

||B_hat - B_0||_F / ||B_0||_F, mean (sd) over 25 replicates.

*Source: sim_f_final_results.rds*

```
   setting          Method rel_coef_error     sd
   (50,50)    VAPOR dl r=3         0.6981 0.0991
   (50,50) VAPOR ridge r=3         1.0945 0.5356
   (50,75)    VAPOR dl r=3         0.5656 0.1452
   (50,75) VAPOR ridge r=3         0.5825 0.3046
  (50,100)    VAPOR dl r=3         0.4776 0.1549
  (50,100) VAPOR ridge r=3         0.4370 0.2101
 (100,100)    VAPOR dl r=3         0.6398 0.0753
 (100,100) VAPOR ridge r=3         0.9365 0.3273
 (150,150)    VAPOR dl r=3         0.5695 0.0728
 (150,150) VAPOR ridge r=3         0.7773 0.2539
```

## Simulation E: edge sensitivity and specificity at the oracle threshold

Glasso's values are reported but are degenerate: its off-diagonal mass is near uniform, so the oracle threshold assigns it almost all or almost no edges.

*Source: sim_f_final_results.rds*

```
   setting          Method   Sens Sens_sd   Spec Spec_sd
   (50,50)    VAPOR dl r=3 0.2972  0.0795 0.9454  0.0110
   (50,50) VAPOR ridge r=3 0.2926  0.0680 0.9450  0.0103
   (50,50)         Kastner 0.1321  0.0421 0.9329  0.0090
   (50,50)          Glasso 0.9309  0.2391 0.0747  0.2584
   (50,75)    VAPOR dl r=3 0.4318  0.0922 0.9560  0.0095
   (50,75) VAPOR ridge r=3 0.4414  0.0960 0.9567  0.0098
   (50,75)         Kastner 0.1416  0.0354 0.9336  0.0087
   (50,75)          Glasso 0.8105  0.3445 0.2253  0.4092
  (50,100)    VAPOR dl r=3 0.5550  0.0760 0.9656  0.0076
  (50,100) VAPOR ridge r=3 0.5664  0.0813 0.9665  0.0077
  (50,100)         Kastner 0.1471  0.0421 0.9340  0.0091
  (50,100)          Glasso 0.6235  0.3718 0.4905  0.4810
 (100,100)    VAPOR dl r=3 0.4678  0.0850 0.9508  0.0146
 (100,100) VAPOR ridge r=3 0.4652  0.0825 0.9505  0.0147
 (100,100)         Kastner 0.1079  0.0244 0.9189  0.0134
 (100,100)          Glasso 0.9008  0.2743 0.1093  0.3022
 (150,150)    VAPOR dl r=3 0.5800  0.0566 0.9611  0.0083
 (150,150) VAPOR ridge r=3 0.5778  0.0561 0.9609  0.0083
 (150,150)         Kastner 0.0949  0.0181 0.9167  0.0093
 (150,150)          Glasso 0.7274  0.4056 0.2927  0.4355
```

## Simulation F: relative coefficient error of VAPOR

*Source: sim_g_results.rds*

```
   setting          Method rel_coef_error     sd
   (50,35)    VAPOR dl r=3         0.8360 0.0765
   (50,35) VAPOR ridge r=3         0.8984 0.2689
  (100,70)    VAPOR dl r=3         0.7932 0.0544
  (100,70) VAPOR ridge r=3         0.7945 0.1540
 (150,105)    VAPOR dl r=3         0.7666 0.0491
 (150,105) VAPOR ridge r=3         0.7151 0.1110
```

## Simulation F: edge sensitivity and specificity at the oracle threshold

*Source: sim_g_results.rds*

```
   setting          Method   Sens Sens_sd   Spec Spec_sd
   (50,35)    VAPOR dl r=3 0.1971  0.0575 0.9379  0.0094
   (50,35) VAPOR ridge r=3 0.1971  0.0591 0.9379  0.0095
   (50,35)         Kastner 0.1179  0.0408 0.9319  0.0083
   (50,35)          Glasso 0.9646  0.1770 0.0373  0.1865
  (100,70)    VAPOR dl r=3 0.2341  0.0444 0.9298  0.0153
  (100,70) VAPOR ridge r=3 0.2295  0.0456 0.9294  0.0153
  (100,70)         Kastner 0.0940  0.0222 0.9177  0.0134
  (100,70)          Glasso 1.0000  0.0000 0.0000  0.0000
 (150,105)    VAPOR dl r=3 0.2361  0.0503 0.9294  0.0111
 (150,105) VAPOR ridge r=3 0.2298  0.0475 0.9288  0.0110
 (150,105)         Kastner 0.0939  0.0140 0.9166  0.0096
 (150,105)          Glasso 1.0000  0.0000 0.0000  0.0000
```

## Cells differing from the archived Gamma(0.1, 0.1) tables: count per table

567 cells differ in all. The two PrecFact rows account for 295 of them: their base measure moved from 0.1 to 1 along with VAPOR's. Kastner and the Glasso rows move only through the paired differences and win counts, which are computed against VAPOR.

*Source: sim_tables_diff.rds*

```
                       table cells_differing
          tab:sim-a-autocorr              78
   tab:sim-a-paired-autocorr              82
 tab:sim-a-paired-samplesize              72
        tab:sim-a-samplesize              62
          tab:sim-b-autocorr              83
        tab:sim-b-samplesize              66
                   tab:sim-c              21
                   tab:sim-d              64
                   tab:sim-e              24
                   tab:sim-f              15
```

## Cells differing: count per method

*Source: sim_tables_diff.rds*

```
            method cells_differing
    VAR + PrecFact             148
 PrecFact (no VAR)             147
             VAPOR              71
 Kastner (VAR+FSV)              57
   Glasso (no VAR)              51
      VAR + Glasso              51
     VAPOR (ridge)              20
        VAPOR (DL)              19
 Sample covariance               3
```

## Every differing cell, archived vs new

Formatting is ignored: <1e-04 and $<10^{-4}$ compare equal.

*Source: sim_tables_diff.rds*

```
                       table scenario            method                    column        archived             new
          tab:sim-a-autocorr      0.1 PrecFact (no VAR)               Sensitivity    0.555(0.192)    0.555(0.188)
          tab:sim-a-autocorr      0.1 PrecFact (no VAR)                 Frobenius    0.784(0.093)    0.775(0.093)
          tab:sim-a-autocorr      0.1 PrecFact (no VAR)              Delta_F (SE)    0.098(0.013)    0.091(0.010)
          tab:sim-a-autocorr      0.1   Glasso (no VAR)              Delta_F (SE)    0.257(0.020)    0.259(0.021)
          tab:sim-a-autocorr      0.1      VAR + Glasso              Delta_F (SE)    0.243(0.019)    0.246(0.019)
          tab:sim-a-autocorr      0.1    VAR + PrecFact               Sensitivity    0.563(0.189)    0.555(0.187)
          tab:sim-a-autocorr      0.1    VAR + PrecFact                 Frobenius    0.767(0.074)    0.751(0.077)
          tab:sim-a-autocorr      0.1    VAR + PrecFact              Delta_F (SE)    0.081(0.012)    0.067(0.010)
          tab:sim-a-autocorr      0.1 Kastner (VAR+FSV)              Delta_F (SE)    0.432(0.075)    0.435(0.076)
          tab:sim-a-autocorr      0.1             VAPOR               Sensitivity    0.563(0.163)    0.567(0.154)
          tab:sim-a-autocorr      0.1             VAPOR               Specificity    1.000(0.002)    0.999(0.002)
          tab:sim-a-autocorr      0.1             VAPOR                 Frobenius    0.686(0.073)    0.684(0.074)
          tab:sim-a-autocorr      0.5 PrecFact (no VAR)               Sensitivity    0.387(0.213)    0.398(0.213)
          tab:sim-a-autocorr      0.5 PrecFact (no VAR)               Specificity    0.973(0.019)    0.970(0.024)
          tab:sim-a-autocorr      0.5 PrecFact (no VAR)                 Frobenius    1.534(0.486)    1.524(0.482)
          tab:sim-a-autocorr      0.5 PrecFact (no VAR)              Delta_F (SE)    0.838(0.090)    0.830(0.089)
          tab:sim-a-autocorr      0.5   Glasso (no VAR)              Delta_F (SE)    0.873(0.068)    0.874(0.069)
          tab:sim-a-autocorr      0.5      VAR + Glasso              Delta_F (SE)    0.236(0.017)    0.237(0.018)
          tab:sim-a-autocorr      0.5    VAR + PrecFact               Sensitivity    0.561(0.178)    0.557(0.184)
          tab:sim-a-autocorr      0.5    VAR + PrecFact               Specificity    0.998(0.003)    0.999(0.003)
          tab:sim-a-autocorr      0.5    VAR + PrecFact                 Frobenius    0.765(0.073)    0.755(0.076)
          tab:sim-a-autocorr      0.5    VAR + PrecFact              Delta_F (SE)    0.069(0.011)    0.060(0.012)
          tab:sim-a-autocorr      0.5 Kastner (VAR+FSV)              Delta_F (SE)    0.503(0.072)    0.504(0.072)
          tab:sim-a-autocorr      0.5             VAPOR               Sensitivity    0.563(0.164)    0.566(0.157)
          tab:sim-a-autocorr      0.5             VAPOR                 Frobenius    0.696(0.071)    0.695(0.069)
          tab:sim-a-autocorr     0.75 PrecFact (no VAR)               Sensitivity    0.286(0.151)    0.304(0.147)
          tab:sim-a-autocorr     0.75 PrecFact (no VAR)               Specificity    0.929(0.041)    0.926(0.036)
          tab:sim-a-autocorr     0.75 PrecFact (no VAR)                 Frobenius    2.359(0.569)    2.363(0.562)
          tab:sim-a-autocorr     0.75 PrecFact (no VAR)              Delta_F (SE)    1.663(0.108)    1.671(0.107)
          tab:sim-a-autocorr     0.75   Glasso (no VAR)              Delta_F (SE)    1.335(0.076)    1.338(0.076)
          tab:sim-a-autocorr     0.75      VAR + Glasso              Delta_F (SE)    0.239(0.017)    0.242(0.017)
          tab:sim-a-autocorr     0.75    VAR + PrecFact               Sensitivity    0.557(0.167)    0.561(0.191)
          tab:sim-a-autocorr     0.75    VAR + PrecFact                 Frobenius    0.771(0.080)    0.761(0.082)
          tab:sim-a-autocorr     0.75    VAR + PrecFact              Delta_F (SE)    0.076(0.012)    0.068(0.012)
          tab:sim-a-autocorr     0.75 Kastner (VAR+FSV)              Delta_F (SE)    0.501(0.077)    0.504(0.077)
          tab:sim-a-autocorr     0.75             VAPOR               Sensitivity    0.554(0.172)    0.553(0.157)
          tab:sim-a-autocorr     0.75             VAPOR               Specificity    0.999(0.004)    1.000(0.002)
          tab:sim-a-autocorr     0.75             VAPOR                 Frobenius    0.695(0.068)    0.693(0.067)
          tab:sim-a-autocorr     0.89 PrecFact (no VAR)               Sensitivity    0.239(0.114)    0.236(0.098)
          tab:sim-a-autocorr     0.89 PrecFact (no VAR)               Specificity    0.885(0.059)    0.888(0.050)
          tab:sim-a-autocorr     0.89 PrecFact (no VAR)                 Frobenius    2.985(0.628)    2.977(0.604)
          tab:sim-a-autocorr     0.89 PrecFact (no VAR)              Delta_F (SE)    2.292(0.121)    2.284(0.117)
          tab:sim-a-autocorr     0.89   Glasso (no VAR)              Delta_F (SE)    1.484(0.088)    1.483(0.088)
          tab:sim-a-autocorr     0.89      VAR + Glasso              Delta_F (SE)    0.242(0.018)    0.241(0.018)
          tab:sim-a-autocorr     0.89    VAR + PrecFact               Sensitivity    0.566(0.180)    0.575(0.171)
          tab:sim-a-autocorr     0.89    VAR + PrecFact               Specificity    0.999(0.003)    0.998(0.004)
          tab:sim-a-autocorr     0.89    VAR + PrecFact                 Frobenius    0.763(0.072)    0.758(0.077)
          tab:sim-a-autocorr     0.89    VAR + PrecFact              Delta_F (SE)    0.071(0.013)    0.064(0.012)
          tab:sim-a-autocorr     0.89 Kastner (VAR+FSV)              Delta_F (SE)    0.450(0.075)    0.449(0.075)
          tab:sim-a-autocorr     0.89             VAPOR               Sensitivity    0.554(0.161)    0.559(0.160)
          tab:sim-a-autocorr     0.89             VAPOR               Specificity    0.999(0.002)    1.000(0.002)
          tab:sim-a-autocorr     0.89             VAPOR                 Frobenius    0.692(0.070)    0.694(0.067)
          tab:sim-a-autocorr     0.95 PrecFact (no VAR)               Sensitivity    0.244(0.110)    0.238(0.111)
          tab:sim-a-autocorr     0.95 PrecFact (no VAR)               Specificity    0.863(0.065)    0.864(0.061)
          tab:sim-a-autocorr     0.95 PrecFact (no VAR)                 Frobenius    3.296(0.689)    3.269(0.656)
          tab:sim-a-autocorr     0.95 PrecFact (no VAR)              Delta_F (SE)    2.605(0.133)    2.580(0.126)
          tab:sim-a-autocorr     0.95   Glasso (no VAR)              Delta_F (SE)    1.503(0.091)    1.505(0.091)
          tab:sim-a-autocorr     0.95      VAR + Glasso              Delta_F (SE)    0.243(0.019)    0.246(0.019)
          tab:sim-a-autocorr     0.95    VAR + PrecFact               Sensitivity    0.573(0.188)    0.555(0.181)
          tab:sim-a-autocorr     0.95    VAR + PrecFact                 Frobenius    0.772(0.082)    0.762(0.082)
          tab:sim-a-autocorr     0.95    VAR + PrecFact              Delta_F (SE)    0.081(0.011)    0.073(0.012)
          tab:sim-a-autocorr     0.95 Kastner (VAR+FSV)              Delta_F (SE)    0.487(0.076)    0.489(0.076)
          tab:sim-a-autocorr     0.95             VAPOR               Sensitivity    0.557(0.162)    0.550(0.169)
          tab:sim-a-autocorr     0.95             VAPOR                 Frobenius    0.691(0.069)    0.689(0.071)
          tab:sim-a-autocorr     0.99 PrecFact (no VAR)               Sensitivity    0.276(0.132)    0.239(0.116)
          tab:sim-a-autocorr     0.99 PrecFact (no VAR)               Specificity    0.829(0.063)    0.845(0.068)
          tab:sim-a-autocorr     0.99 PrecFact (no VAR)                 Frobenius    3.582(0.759)    3.586(0.741)
          tab:sim-a-autocorr     0.99 PrecFact (no VAR)              Delta_F (SE)    2.893(0.149)    2.898(0.145)
          tab:sim-a-autocorr     0.99   Glasso (no VAR)              Delta_F (SE)    1.509(0.092)    1.509(0.091)
          tab:sim-a-autocorr     0.99      VAR + Glasso              Delta_F (SE)    0.243(0.020)    0.243(0.019)
          tab:sim-a-autocorr     0.99    VAR + PrecFact               Sensitivity    0.559(0.173)    0.569(0.184)
          tab:sim-a-autocorr     0.99    VAR + PrecFact               Specificity    0.998(0.005)    0.999(0.004)
          tab:sim-a-autocorr     0.99    VAR + PrecFact                 Frobenius    0.772(0.083)    0.756(0.080)
          tab:sim-a-autocorr     0.99    VAR + PrecFact              Delta_F (SE)    0.084(0.012)    0.067(0.011)
          tab:sim-a-autocorr     0.99 Kastner (VAR+FSV)              Delta_F (SE)    0.480(0.073)    0.481(0.072)
          tab:sim-a-autocorr     0.99             VAPOR               Sensitivity    0.558(0.151)    0.557(0.158)
          tab:sim-a-autocorr     0.99             VAPOR               Specificity    1.000(0.002)    1.000(0.001)
          tab:sim-a-autocorr     0.99             VAPOR                 Frobenius    0.689(0.066)    0.688(0.068)
        tab:sim-a-samplesize       50 PrecFact (no VAR)               Sensitivity    0.099(0.087)    0.095(0.080)
        tab:sim-a-samplesize       50 PrecFact (no VAR)               Specificity    0.941(0.029)    0.950(0.030)
        tab:sim-a-samplesize       50 PrecFact (no VAR)                 Frobenius    3.229(0.638)    3.199(0.624)
        tab:sim-a-samplesize       50 PrecFact (no VAR)              Delta_F (SE)    1.581(0.120)    1.547(0.116)
        tab:sim-a-samplesize       50   Glasso (no VAR)              Delta_F (SE)    0.547(0.039)    0.544(0.036)
        tab:sim-a-samplesize       50      VAR + Glasso              Delta_F (SE)    0.062(0.025)    0.059(0.024)
        tab:sim-a-samplesize       50    VAR + PrecFact               Sensitivity    0.218(0.159)    0.192(0.137)
        tab:sim-a-samplesize       50    VAR + PrecFact               Specificity    0.992(0.011)    0.994(0.011)
        tab:sim-a-samplesize       50    VAR + PrecFact                 Frobenius    1.750(0.311)    1.740(0.305)
        tab:sim-a-samplesize       50    VAR + PrecFact              Delta_F (SE)    0.101(0.065)    0.088(0.064)
        tab:sim-a-samplesize       50 Kastner (VAR+FSV)              Delta_F (SE)    0.147(0.035)    0.144(0.033)
        tab:sim-a-samplesize       50             VAPOR               Sensitivity    0.029(0.041)    0.020(0.039)
        tab:sim-a-samplesize       50             VAPOR                 Frobenius    1.649(0.351)    1.652(0.357)
        tab:sim-a-samplesize      100 PrecFact (no VAR)               Sensitivity    0.164(0.109)    0.170(0.099)
        tab:sim-a-samplesize      100 PrecFact (no VAR)               Specificity    0.894(0.047)    0.900(0.043)
        tab:sim-a-samplesize      100 PrecFact (no VAR)                 Frobenius    3.244(0.646)    3.249(0.663)
        tab:sim-a-samplesize      100 PrecFact (no VAR)              Delta_F (SE)    2.203(0.128)    2.211(0.130)
        tab:sim-a-samplesize      100   Glasso (no VAR)              Delta_F (SE)    1.149(0.090)    1.151(0.089)
        tab:sim-a-samplesize      100      VAR + Glasso              Delta_F (SE)    0.280(0.030)    0.283(0.030)
        tab:sim-a-samplesize      100    VAR + PrecFact               Sensitivity    0.444(0.199)    0.424(0.187)
        tab:sim-a-samplesize      100    VAR + PrecFact               Specificity    0.997(0.005)    0.996(0.006)
        tab:sim-a-samplesize      100    VAR + PrecFact                 Frobenius    1.144(0.120)    1.143(0.116)
        tab:sim-a-samplesize      100    VAR + PrecFact              Delta_F (SE)    0.103(0.021)    0.105(0.021)
        tab:sim-a-samplesize      100 Kastner (VAR+FSV)              Delta_F (SE)    0.345(0.069)    0.348(0.067)
        tab:sim-a-samplesize      100             VAPOR               Sensitivity    0.365(0.159)    0.372(0.145)
        tab:sim-a-samplesize      100             VAPOR                 Frobenius    1.041(0.158)    1.038(0.157)
        tab:sim-a-samplesize      200 PrecFact (no VAR)               Sensitivity    0.244(0.110)    0.238(0.111)
        tab:sim-a-samplesize      200 PrecFact (no VAR)               Specificity    0.863(0.065)    0.864(0.061)
        tab:sim-a-samplesize      200 PrecFact (no VAR)                 Frobenius    3.296(0.689)    3.269(0.656)
        tab:sim-a-samplesize      200 PrecFact (no VAR)              Delta_F (SE)    2.605(0.133)    2.580(0.126)
        tab:sim-a-samplesize      200   Glasso (no VAR)              Delta_F (SE)    1.503(0.091)    1.505(0.091)
        tab:sim-a-samplesize      200      VAR + Glasso              Delta_F (SE)    0.243(0.019)    0.246(0.019)
        tab:sim-a-samplesize      200    VAR + PrecFact               Sensitivity    0.573(0.188)    0.555(0.181)
        tab:sim-a-samplesize      200    VAR + PrecFact                 Frobenius    0.772(0.082)    0.762(0.082)
        tab:sim-a-samplesize      200    VAR + PrecFact              Delta_F (SE)    0.081(0.011)    0.073(0.012)
        tab:sim-a-samplesize      200 Kastner (VAR+FSV)              Delta_F (SE)    0.487(0.076)    0.489(0.076)
        tab:sim-a-samplesize      200             VAPOR               Sensitivity    0.557(0.162)    0.550(0.169)
        tab:sim-a-samplesize      200             VAPOR                 Frobenius    0.691(0.069)    0.689(0.071)
        tab:sim-a-samplesize      500 PrecFact (no VAR)               Sensitivity    0.343(0.166)    0.339(0.151)
        tab:sim-a-samplesize      500 PrecFact (no VAR)               Specificity    0.797(0.085)    0.810(0.076)
        tab:sim-a-samplesize      500 PrecFact (no VAR)                 Frobenius    3.330(0.712)    3.357(0.743)
        tab:sim-a-samplesize      500 PrecFact (no VAR)              Delta_F (SE)    2.889(0.143)    2.918(0.150)
        tab:sim-a-samplesize      500   Glasso (no VAR)              Delta_F (SE)    1.745(0.092)    1.746(0.092)
        tab:sim-a-samplesize      500      VAR + Glasso              Delta_F (SE)    0.166(0.010)    0.168(0.010)
        tab:sim-a-samplesize      500    VAR + PrecFact               Sensitivity    0.691(0.140)    0.697(0.143)
        tab:sim-a-samplesize      500    VAR + PrecFact                 Frobenius    0.485(0.056)    0.472(0.056)
        tab:sim-a-samplesize      500    VAR + PrecFact              Delta_F (SE)    0.044(0.007)    0.033(0.006)
        tab:sim-a-samplesize      500 Kastner (VAR+FSV)              Delta_F (SE)    0.610(0.088)    0.612(0.089)
        tab:sim-a-samplesize      500             VAPOR               Sensitivity    0.690(0.130)    0.691(0.131)
        tab:sim-a-samplesize      500             VAPOR               Specificity    1.000(0.001)    1.000(0.000)
        tab:sim-a-samplesize      500             VAPOR                 Frobenius    0.440(0.050)    0.439(0.049)
        tab:sim-a-samplesize     1000 PrecFact (no VAR)               Sensitivity    0.363(0.119)    0.384(0.146)
        tab:sim-a-samplesize     1000 PrecFact (no VAR)               Specificity    0.781(0.071)    0.782(0.071)
        tab:sim-a-samplesize     1000 PrecFact (no VAR)                 Frobenius    3.319(0.722)    3.327(0.706)
        tab:sim-a-samplesize     1000 PrecFact (no VAR)              Delta_F (SE)    3.008(0.145)    3.015(0.141)
        tab:sim-a-samplesize     1000   Glasso (no VAR)              Delta_F (SE)    1.871(0.092)    1.870(0.092)
        tab:sim-a-samplesize     1000    VAR + PrecFact               Sensitivity    0.733(0.134)    0.738(0.130)
        tab:sim-a-samplesize     1000    VAR + PrecFact               Specificity    0.999(0.002)    1.000(0.002)
        tab:sim-a-samplesize     1000    VAR + PrecFact                 Frobenius    0.344(0.055)    0.335(0.046)
        tab:sim-a-samplesize     1000    VAR + PrecFact              Delta_F (SE)    0.033(0.007)    0.023(0.004)
        tab:sim-a-samplesize     1000             VAPOR               Sensitivity    0.737(0.128)    0.736(0.124)
        tab:sim-a-samplesize     1000             VAPOR                 Frobenius    0.311(0.038)    0.312(0.037)
   tab:sim-a-paired-autocorr      0.1 PrecFact (no VAR)                   Delta_F           0.098           0.091
   tab:sim-a-paired-autocorr      0.1 PrecFact (no VAR)                        SE           0.013           0.010
   tab:sim-a-paired-autocorr      0.1 PrecFact (no VAR)                         t            7.75            8.81
   tab:sim-a-paired-autocorr      0.1 PrecFact (no VAR)                      Wins           24/25           23/25
   tab:sim-a-paired-autocorr      0.1   Glasso (no VAR)                   Delta_F           0.257           0.259
   tab:sim-a-paired-autocorr      0.1   Glasso (no VAR)                        SE           0.020           0.021
   tab:sim-a-paired-autocorr      0.1   Glasso (no VAR)                         t           12.61           12.50
   tab:sim-a-paired-autocorr      0.1      VAR + Glasso                   Delta_F           0.243           0.246
   tab:sim-a-paired-autocorr      0.1      VAR + Glasso                         t           12.90           12.77
   tab:sim-a-paired-autocorr      0.1    VAR + PrecFact                   Delta_F           0.081           0.067
   tab:sim-a-paired-autocorr      0.1    VAR + PrecFact                        SE           0.012           0.010
   tab:sim-a-paired-autocorr      0.1    VAR + PrecFact                         t            6.55            6.48
   tab:sim-a-paired-autocorr      0.1 Kastner (VAR+FSV)                   Delta_F           0.432           0.435
   tab:sim-a-paired-autocorr      0.1 Kastner (VAR+FSV)                        SE           0.075           0.076
   tab:sim-a-paired-autocorr      0.5 PrecFact (no VAR)                   Delta_F           0.838           0.830
   tab:sim-a-paired-autocorr      0.5 PrecFact (no VAR)                        SE           0.090           0.089
   tab:sim-a-paired-autocorr      0.5 PrecFact (no VAR)                         t            9.33            9.36
   tab:sim-a-paired-autocorr      0.5   Glasso (no VAR)                   Delta_F           0.873           0.874
   tab:sim-a-paired-autocorr      0.5   Glasso (no VAR)                        SE           0.068           0.069
   tab:sim-a-paired-autocorr      0.5      VAR + Glasso                   Delta_F           0.236           0.237
   tab:sim-a-paired-autocorr      0.5      VAR + Glasso                        SE           0.017           0.018
   tab:sim-a-paired-autocorr      0.5      VAR + Glasso                         t           13.79           13.52
   tab:sim-a-paired-autocorr      0.5    VAR + PrecFact                   Delta_F           0.069           0.060
   tab:sim-a-paired-autocorr      0.5    VAR + PrecFact                        SE           0.011           0.012
   tab:sim-a-paired-autocorr      0.5    VAR + PrecFact                         t            6.07            4.92
   tab:sim-a-paired-autocorr      0.5    VAR + PrecFact                      Wins           22/25           20/25
   tab:sim-a-paired-autocorr      0.5 Kastner (VAR+FSV)                   Delta_F           0.503           0.504
   tab:sim-a-paired-autocorr      0.5 Kastner (VAR+FSV)                         t            7.03            7.02
   tab:sim-a-paired-autocorr     0.75 PrecFact (no VAR)                   Delta_F           1.663           1.671
   tab:sim-a-paired-autocorr     0.75 PrecFact (no VAR)                        SE           0.108           0.107
   tab:sim-a-paired-autocorr     0.75 PrecFact (no VAR)                         t           15.45           15.56
   tab:sim-a-paired-autocorr     0.75   Glasso (no VAR)                   Delta_F           1.335           1.338
   tab:sim-a-paired-autocorr     0.75   Glasso (no VAR)                         t           17.61           17.56
   tab:sim-a-paired-autocorr     0.75      VAR + Glasso                   Delta_F           0.239           0.242
   tab:sim-a-paired-autocorr     0.75      VAR + Glasso                         t           13.98           14.16
   tab:sim-a-paired-autocorr     0.75    VAR + PrecFact                   Delta_F           0.076           0.068
   tab:sim-a-paired-autocorr     0.75    VAR + PrecFact                         t            6.46            5.54
   tab:sim-a-paired-autocorr     0.75    VAR + PrecFact                      Wins           21/25           18/25
   tab:sim-a-paired-autocorr     0.75 Kastner (VAR+FSV)                   Delta_F           0.501           0.504
   tab:sim-a-paired-autocorr     0.89 PrecFact (no VAR)                   Delta_F           2.292           2.284
   tab:sim-a-paired-autocorr     0.89 PrecFact (no VAR)                        SE           0.121           0.117
   tab:sim-a-paired-autocorr     0.89 PrecFact (no VAR)                         t           18.95           19.59
   tab:sim-a-paired-autocorr     0.89   Glasso (no VAR)                   Delta_F           1.484           1.483
   tab:sim-a-paired-autocorr     0.89   Glasso (no VAR)                         t           16.93           16.89
   tab:sim-a-paired-autocorr     0.89      VAR + Glasso                   Delta_F           0.242           0.241
   tab:sim-a-paired-autocorr     0.89      VAR + Glasso                         t           13.34           13.19
   tab:sim-a-paired-autocorr     0.89    VAR + PrecFact                   Delta_F           0.071           0.064
   tab:sim-a-paired-autocorr     0.89    VAR + PrecFact                        SE           0.013           0.012
   tab:sim-a-paired-autocorr     0.89    VAR + PrecFact                         t            5.40            5.45
   tab:sim-a-paired-autocorr     0.89    VAR + PrecFact                      Wins           22/25           20/25
   tab:sim-a-paired-autocorr     0.89 Kastner (VAR+FSV)                   Delta_F           0.450           0.449
   tab:sim-a-paired-autocorr     0.89 Kastner (VAR+FSV)                         t            6.01            5.96
   tab:sim-a-paired-autocorr     0.89 Kastner (VAR+FSV)                      Wins           25/25           24/25
   tab:sim-a-paired-autocorr     0.95 PrecFact (no VAR)                   Delta_F           2.605           2.580
   tab:sim-a-paired-autocorr     0.95 PrecFact (no VAR)                        SE           0.133           0.126
   tab:sim-a-paired-autocorr     0.95 PrecFact (no VAR)                         t           19.59           20.44
   tab:sim-a-paired-autocorr     0.95   Glasso (no VAR)                   Delta_F           1.503           1.505
   tab:sim-a-paired-autocorr     0.95   Glasso (no VAR)                         t           16.51           16.58
   tab:sim-a-paired-autocorr     0.95      VAR + Glasso                   Delta_F           0.243           0.246
   tab:sim-a-paired-autocorr     0.95      VAR + Glasso                         t           12.87           12.93
   tab:sim-a-paired-autocorr     0.95    VAR + PrecFact                   Delta_F           0.081           0.073
   tab:sim-a-paired-autocorr     0.95    VAR + PrecFact                        SE           0.011           0.012
   tab:sim-a-paired-autocorr     0.95    VAR + PrecFact                         t            7.12            6.18
   tab:sim-a-paired-autocorr     0.95    VAR + PrecFact                      Wins           23/25           21/25
   tab:sim-a-paired-autocorr     0.95 Kastner (VAR+FSV)                   Delta_F           0.487           0.489
   tab:sim-a-paired-autocorr     0.95 Kastner (VAR+FSV)                         t            6.44            6.47
   tab:sim-a-paired-autocorr     0.95 Kastner (VAR+FSV)                      Wins           23/25           24/25
   tab:sim-a-paired-autocorr     0.99 PrecFact (no VAR)                   Delta_F           2.893           2.898
   tab:sim-a-paired-autocorr     0.99 PrecFact (no VAR)                        SE           0.149           0.145
   tab:sim-a-paired-autocorr     0.99 PrecFact (no VAR)                         t           19.42           20.01
   tab:sim-a-paired-autocorr     0.99   Glasso (no VAR)                        SE           0.092           0.091
   tab:sim-a-paired-autocorr     0.99   Glasso (no VAR)                         t           16.40           16.55
   tab:sim-a-paired-autocorr     0.99      VAR + Glasso                        SE           0.020           0.019
   tab:sim-a-paired-autocorr     0.99      VAR + Glasso                         t           12.25           12.54
   tab:sim-a-paired-autocorr     0.99    VAR + PrecFact                   Delta_F           0.084           0.067
   tab:sim-a-paired-autocorr     0.99    VAR + PrecFact                        SE           0.012           0.011
   tab:sim-a-paired-autocorr     0.99    VAR + PrecFact                         t            6.72            6.01
   tab:sim-a-paired-autocorr     0.99    VAR + PrecFact                      Wins           23/25           19/25
   tab:sim-a-paired-autocorr     0.99 Kastner (VAR+FSV)                   Delta_F           0.480           0.481
   tab:sim-a-paired-autocorr     0.99 Kastner (VAR+FSV)                        SE           0.073           0.072
   tab:sim-a-paired-autocorr     0.99 Kastner (VAR+FSV)                         t            6.55            6.65
   tab:sim-a-paired-autocorr     0.99 Kastner (VAR+FSV)                      Wins           25/25           24/25
 tab:sim-a-paired-samplesize       50 PrecFact (no VAR)                   Delta_F           1.581           1.547
 tab:sim-a-paired-samplesize       50 PrecFact (no VAR)                        SE           0.120           0.116
 tab:sim-a-paired-samplesize       50 PrecFact (no VAR)                         t           13.14           13.29
 tab:sim-a-paired-samplesize       50   Glasso (no VAR)                   Delta_F           0.547           0.544
 tab:sim-a-paired-samplesize       50   Glasso (no VAR)                        SE           0.039           0.036
 tab:sim-a-paired-samplesize       50   Glasso (no VAR)                         t           14.18           14.96
 tab:sim-a-paired-samplesize       50      VAR + Glasso                   Delta_F           0.062           0.059
 tab:sim-a-paired-samplesize       50      VAR + Glasso                        SE           0.025           0.024
 tab:sim-a-paired-samplesize       50      VAR + Glasso                         t            2.51            2.44
 tab:sim-a-paired-samplesize       50      VAR + Glasso                         p           0.019           0.022
 tab:sim-a-paired-samplesize       50      VAR + Glasso                      Wins           16/25           17/25
 tab:sim-a-paired-samplesize       50    VAR + PrecFact                   Delta_F           0.101           0.088
 tab:sim-a-paired-samplesize       50    VAR + PrecFact                        SE           0.065           0.064
 tab:sim-a-paired-samplesize       50    VAR + PrecFact                         t            1.56            1.37
 tab:sim-a-paired-samplesize       50    VAR + PrecFact                         p            0.13            0.18
 tab:sim-a-paired-samplesize       50    VAR + PrecFact                      Wins           19/25           18/25
 tab:sim-a-paired-samplesize       50 Kastner (VAR+FSV)                   Delta_F           0.147           0.144
 tab:sim-a-paired-samplesize       50 Kastner (VAR+FSV)                        SE           0.035           0.033
 tab:sim-a-paired-samplesize       50 Kastner (VAR+FSV)                         t            4.21            4.33
 tab:sim-a-paired-samplesize       50 Kastner (VAR+FSV)                         p         0.00031         0.00023
 tab:sim-a-paired-samplesize      100 PrecFact (no VAR)                   Delta_F           2.203           2.211
 tab:sim-a-paired-samplesize      100 PrecFact (no VAR)                        SE           0.128           0.130
 tab:sim-a-paired-samplesize      100 PrecFact (no VAR)                         t           17.21           17.02
 tab:sim-a-paired-samplesize      100   Glasso (no VAR)                   Delta_F           1.149           1.151
 tab:sim-a-paired-samplesize      100   Glasso (no VAR)                        SE           0.090           0.089
 tab:sim-a-paired-samplesize      100   Glasso (no VAR)                         t           12.78           13.00
 tab:sim-a-paired-samplesize      100      VAR + Glasso                   Delta_F           0.280           0.283
 tab:sim-a-paired-samplesize      100      VAR + Glasso                         t            9.19            9.58
 tab:sim-a-paired-samplesize      100    VAR + PrecFact                   Delta_F           0.103           0.105
 tab:sim-a-paired-samplesize      100    VAR + PrecFact                         t            4.89            5.06
 tab:sim-a-paired-samplesize      100    VAR + PrecFact                      Wins           19/25           22/25
 tab:sim-a-paired-samplesize      100 Kastner (VAR+FSV)                   Delta_F           0.345           0.348
 tab:sim-a-paired-samplesize      100 Kastner (VAR+FSV)                        SE           0.069           0.067
 tab:sim-a-paired-samplesize      100 Kastner (VAR+FSV)                         t            5.02            5.17
 tab:sim-a-paired-samplesize      200 PrecFact (no VAR)                   Delta_F           2.605           2.580
 tab:sim-a-paired-samplesize      200 PrecFact (no VAR)                        SE           0.133           0.126
 tab:sim-a-paired-samplesize      200 PrecFact (no VAR)                         t           19.59           20.44
 tab:sim-a-paired-samplesize      200   Glasso (no VAR)                   Delta_F           1.503           1.505
 tab:sim-a-paired-samplesize      200   Glasso (no VAR)                         t           16.51           16.58
 tab:sim-a-paired-samplesize      200      VAR + Glasso                   Delta_F           0.243           0.246
 tab:sim-a-paired-samplesize      200      VAR + Glasso                         t           12.87           12.93
 tab:sim-a-paired-samplesize      200    VAR + PrecFact                   Delta_F           0.081           0.073
 tab:sim-a-paired-samplesize      200    VAR + PrecFact                        SE           0.011           0.012
 tab:sim-a-paired-samplesize      200    VAR + PrecFact                         t            7.12            6.18
 tab:sim-a-paired-samplesize      200    VAR + PrecFact                      Wins           23/25           21/25
 tab:sim-a-paired-samplesize      200 Kastner (VAR+FSV)                   Delta_F           0.487           0.489
 tab:sim-a-paired-samplesize      200 Kastner (VAR+FSV)                         t            6.44            6.47
 tab:sim-a-paired-samplesize      200 Kastner (VAR+FSV)                      Wins           23/25           24/25
 tab:sim-a-paired-samplesize      500 PrecFact (no VAR)                   Delta_F           2.889           2.918
 tab:sim-a-paired-samplesize      500 PrecFact (no VAR)                        SE           0.143           0.150
 tab:sim-a-paired-samplesize      500 PrecFact (no VAR)                         t           20.22           19.50
 tab:sim-a-paired-samplesize      500   Glasso (no VAR)                   Delta_F           1.745           1.746
 tab:sim-a-paired-samplesize      500   Glasso (no VAR)                         t           18.94           18.89
 tab:sim-a-paired-samplesize      500      VAR + Glasso                   Delta_F           0.166           0.168
 tab:sim-a-paired-samplesize      500      VAR + Glasso                         t           17.47           16.49
 tab:sim-a-paired-samplesize      500    VAR + PrecFact                   Delta_F           0.044           0.033
 tab:sim-a-paired-samplesize      500    VAR + PrecFact                        SE           0.007           0.006
 tab:sim-a-paired-samplesize      500    VAR + PrecFact                         t            5.98            5.82
 tab:sim-a-paired-samplesize      500    VAR + PrecFact                      Wins           24/25           23/25
 tab:sim-a-paired-samplesize      500 Kastner (VAR+FSV)                   Delta_F           0.610           0.612
 tab:sim-a-paired-samplesize      500 Kastner (VAR+FSV)                        SE           0.088           0.089
 tab:sim-a-paired-samplesize     1000 PrecFact (no VAR)                   Delta_F           3.008           3.015
 tab:sim-a-paired-samplesize     1000 PrecFact (no VAR)                        SE           0.145           0.141
 tab:sim-a-paired-samplesize     1000 PrecFact (no VAR)                         t           20.77           21.33
 tab:sim-a-paired-samplesize     1000   Glasso (no VAR)                   Delta_F           1.871           1.870
 tab:sim-a-paired-samplesize     1000   Glasso (no VAR)                         t           20.31           20.29
 tab:sim-a-paired-samplesize     1000      VAR + Glasso                         t           19.98           20.23
 tab:sim-a-paired-samplesize     1000    VAR + PrecFact                   Delta_F           0.033           0.023
 tab:sim-a-paired-samplesize     1000    VAR + PrecFact                        SE           0.007           0.004
 tab:sim-a-paired-samplesize     1000    VAR + PrecFact                         t            4.52            5.14
 tab:sim-a-paired-samplesize     1000    VAR + PrecFact                         p         0.00014         <10^{-4
 tab:sim-a-paired-samplesize     1000 Kastner (VAR+FSV)                         t            7.64            7.62
          tab:sim-b-autocorr      0.1 PrecFact (no VAR)                      Sens    0.621(0.136)    0.633(0.120)
          tab:sim-b-autocorr      0.1 PrecFact (no VAR)                      Spec    0.982(0.012)    0.985(0.010)
          tab:sim-b-autocorr      0.1 PrecFact (no VAR)                      Frob    1.095(0.162)    1.069(0.159)
          tab:sim-b-autocorr      0.1 PrecFact (no VAR)              Delta_F (SE)    0.027(0.017)    0.042(0.012)
          tab:sim-b-autocorr      0.1   Glasso (no VAR)              Delta_F (SE)   -0.435(0.029)   -0.393(0.030)
          tab:sim-b-autocorr      0.1      VAR + Glasso              Delta_F (SE)   -0.435(0.028)   -0.393(0.030)
          tab:sim-b-autocorr      0.1    VAR + PrecFact                      Sens    0.635(0.132)    0.635(0.129)
          tab:sim-b-autocorr      0.1    VAR + PrecFact                      Spec    0.983(0.013)    0.985(0.012)
          tab:sim-b-autocorr      0.1    VAR + PrecFact                      Frob    1.076(0.171)    1.068(0.179)
          tab:sim-b-autocorr      0.1    VAR + PrecFact              Delta_F (SE)    0.008(0.018)    0.042(0.009)
          tab:sim-b-autocorr      0.1 Kastner (VAR+FSV)              Delta_F (SE)   -0.174(0.027)   -0.133(0.021)
          tab:sim-b-autocorr      0.1             VAPOR                      Sens    0.590(0.139)    0.613(0.119)
          tab:sim-b-autocorr      0.1             VAPOR                      Spec    0.988(0.011)    0.992(0.008)
          tab:sim-b-autocorr      0.1             VAPOR                      Frob    1.068(0.166)    1.027(0.180)
          tab:sim-b-autocorr      0.5 PrecFact (no VAR)                      Sens    0.461(0.121)    0.454(0.129)
          tab:sim-b-autocorr      0.5 PrecFact (no VAR)                      Spec    0.977(0.023)    0.978(0.016)
          tab:sim-b-autocorr      0.5 PrecFact (no VAR)                      Frob    1.343(0.249)    1.337(0.243)
          tab:sim-b-autocorr      0.5 PrecFact (no VAR)              Delta_F (SE)    0.276(0.044)    0.311(0.040)
          tab:sim-b-autocorr      0.5   Glasso (no VAR)              Delta_F (SE)   -0.042(0.044)   -0.002(0.042)
          tab:sim-b-autocorr      0.5      VAR + Glasso              Delta_F (SE)   -0.429(0.028)   -0.388(0.029)
          tab:sim-b-autocorr      0.5    VAR + PrecFact                      Sens    0.625(0.132)    0.636(0.120)
          tab:sim-b-autocorr      0.5    VAR + PrecFact                      Spec    0.986(0.010)    0.984(0.012)
          tab:sim-b-autocorr      0.5    VAR + PrecFact                      Frob    1.076(0.171)    1.082(0.183)
          tab:sim-b-autocorr      0.5    VAR + PrecFact              Delta_F (SE)    0.009(0.013)    0.055(0.009)
          tab:sim-b-autocorr      0.5 Kastner (VAR+FSV)              Delta_F (SE)   -0.185(0.022)   -0.145(0.022)
          tab:sim-b-autocorr      0.5             VAPOR                      Sens    0.595(0.127)    0.581(0.141)
          tab:sim-b-autocorr      0.5             VAPOR                      Spec    0.991(0.007)    0.991(0.009)
          tab:sim-b-autocorr      0.5             VAPOR                      Frob    1.067(0.160)    1.027(0.176)
          tab:sim-b-autocorr     0.75 PrecFact (no VAR)                      Sens    0.354(0.159)    0.352(0.129)
          tab:sim-b-autocorr     0.75 PrecFact (no VAR)                      Spec    0.948(0.025)    0.949(0.026)
          tab:sim-b-autocorr     0.75 PrecFact (no VAR)                      Frob    1.843(0.333)    1.853(0.337)
          tab:sim-b-autocorr     0.75 PrecFact (no VAR)              Delta_F (SE)    0.796(0.061)    0.832(0.063)
          tab:sim-b-autocorr     0.75   Glasso (no VAR)              Delta_F (SE)    0.423(0.048)    0.450(0.047)
          tab:sim-b-autocorr     0.75      VAR + Glasso              Delta_F (SE)   -0.409(0.031)   -0.382(0.029)
          tab:sim-b-autocorr     0.75    VAR + PrecFact                      Sens    0.619(0.123)    0.649(0.130)
          tab:sim-b-autocorr     0.75    VAR + PrecFact                      Spec    0.986(0.012)    0.984(0.012)
          tab:sim-b-autocorr     0.75    VAR + PrecFact                      Frob    1.111(0.179)    1.090(0.183)
          tab:sim-b-autocorr     0.75    VAR + PrecFact              Delta_F (SE)    0.064(0.018)    0.069(0.011)
          tab:sim-b-autocorr     0.75 Kastner (VAR+FSV)              Delta_F (SE)   -0.142(0.015)   -0.116(0.016)
          tab:sim-b-autocorr     0.75             VAPOR                      Sens    0.590(0.152)    0.592(0.139)
          tab:sim-b-autocorr     0.75             VAPOR                      Spec    0.991(0.009)    0.992(0.008)
          tab:sim-b-autocorr     0.75             VAPOR                      Frob    1.047(0.176)    1.021(0.168)
          tab:sim-b-autocorr     0.89 PrecFact (no VAR)                      Sens    0.302(0.139)    0.298(0.126)
          tab:sim-b-autocorr     0.89 PrecFact (no VAR)                      Spec    0.923(0.035)    0.925(0.033)
          tab:sim-b-autocorr     0.89 PrecFact (no VAR)                      Frob    2.305(0.432)    2.308(0.443)
          tab:sim-b-autocorr     0.89 PrecFact (no VAR)              Delta_F (SE)    1.233(0.082)    1.268(0.086)
          tab:sim-b-autocorr     0.89   Glasso (no VAR)              Delta_F (SE)    0.578(0.031)    0.609(0.030)
          tab:sim-b-autocorr     0.89      VAR + Glasso              Delta_F (SE)   -0.433(0.033)   -0.402(0.032)
          tab:sim-b-autocorr     0.89    VAR + PrecFact                      Sens    0.609(0.152)    0.652(0.132)
          tab:sim-b-autocorr     0.89    VAR + PrecFact                      Spec    0.982(0.013)    0.982(0.012)
          tab:sim-b-autocorr     0.89    VAR + PrecFact                      Frob    1.098(0.189)    1.082(0.172)
          tab:sim-b-autocorr     0.89    VAR + PrecFact              Delta_F (SE)    0.026(0.016)    0.041(0.011)
          tab:sim-b-autocorr     0.89 Kastner (VAR+FSV)              Delta_F (SE)   -0.174(0.021)   -0.143(0.022)
          tab:sim-b-autocorr     0.89             VAPOR                      Sens    0.592(0.120)    0.585(0.126)
          tab:sim-b-autocorr     0.89             VAPOR                      Spec    0.990(0.010)    0.991(0.009)
          tab:sim-b-autocorr     0.89             VAPOR                      Frob    1.072(0.194)    1.041(0.183)
          tab:sim-b-autocorr     0.95 PrecFact (no VAR)                      Sens    0.267(0.133)    0.263(0.132)
          tab:sim-b-autocorr     0.95 PrecFact (no VAR)                      Spec    0.901(0.045)    0.909(0.033)
          tab:sim-b-autocorr     0.95 PrecFact (no VAR)                      Frob    2.596(0.517)    2.610(0.520)
          tab:sim-b-autocorr     0.95 PrecFact (no VAR)              Delta_F (SE)    1.537(0.102)    1.581(0.099)
          tab:sim-b-autocorr     0.95   Glasso (no VAR)              Delta_F (SE)    0.624(0.027)    0.653(0.024)
          tab:sim-b-autocorr     0.95      VAR + Glasso              Delta_F (SE)   -0.421(0.030)   -0.392(0.029)
          tab:sim-b-autocorr     0.95    VAR + PrecFact                      Sens    0.636(0.132)    0.630(0.157)
          tab:sim-b-autocorr     0.95    VAR + PrecFact                      Frob    1.103(0.193)    1.085(0.169)
          tab:sim-b-autocorr     0.95    VAR + PrecFact              Delta_F (SE)    0.045(0.014)    0.055(0.008)
          tab:sim-b-autocorr     0.95 Kastner (VAR+FSV)              Delta_F (SE)   -0.160(0.020)   -0.131(0.022)
          tab:sim-b-autocorr     0.95             VAPOR                      Sens    0.575(0.145)    0.588(0.111)
          tab:sim-b-autocorr     0.95             VAPOR                      Spec    0.990(0.008)    0.992(0.008)
          tab:sim-b-autocorr     0.95             VAPOR                      Frob    1.059(0.175)    1.030(0.172)
          tab:sim-b-autocorr     0.99 PrecFact (no VAR)                      Sens    0.297(0.124)    0.305(0.120)
          tab:sim-b-autocorr     0.99 PrecFact (no VAR)                      Spec    0.869(0.068)    0.877(0.062)
          tab:sim-b-autocorr     0.99 PrecFact (no VAR)                      Frob    2.947(0.686)    2.909(0.679)
          tab:sim-b-autocorr     0.99 PrecFact (no VAR)              Delta_F (SE)    1.871(0.137)    1.872(0.137)
          tab:sim-b-autocorr     0.99   Glasso (no VAR)              Delta_F (SE)    0.617(0.025)    0.656(0.023)
          tab:sim-b-autocorr     0.99      VAR + Glasso              Delta_F (SE)   -0.440(0.033)   -0.400(0.027)
          tab:sim-b-autocorr     0.99    VAR + PrecFact                      Sens    0.633(0.170)    0.629(0.103)
          tab:sim-b-autocorr     0.99    VAR + PrecFact                      Spec    0.981(0.015)    0.983(0.012)
          tab:sim-b-autocorr     0.99    VAR + PrecFact                      Frob    1.101(0.174)    1.088(0.180)
          tab:sim-b-autocorr     0.99    VAR + PrecFact              Delta_F (SE)    0.025(0.015)    0.051(0.013)
          tab:sim-b-autocorr     0.99 Kastner (VAR+FSV)              Delta_F (SE)   -0.188(0.021)   -0.149(0.018)
          tab:sim-b-autocorr     0.99             VAPOR                      Sens    0.599(0.139)    0.575(0.150)
          tab:sim-b-autocorr     0.99             VAPOR                      Spec    0.988(0.010)    0.992(0.007)
          tab:sim-b-autocorr     0.99             VAPOR                      Frob    1.076(0.186)    1.037(0.163)
        tab:sim-b-samplesize       50 PrecFact (no VAR)                      Sens    0.078(0.078)    0.082(0.086)
        tab:sim-b-samplesize       50 PrecFact (no VAR)                      Spec    0.962(0.032)    0.965(0.030)
        tab:sim-b-samplesize       50 PrecFact (no VAR)                      Frob    2.531(0.548)    2.534(0.546)
        tab:sim-b-samplesize       50 PrecFact (no VAR)              Delta_F (SE)    0.987(0.103)    0.997(0.103)
        tab:sim-b-samplesize       50   Glasso (no VAR)              Delta_F (SE)    0.141(0.020)    0.148(0.021)
        tab:sim-b-samplesize       50      VAR + Glasso              Delta_F (SE)   -0.260(0.021)   -0.253(0.020)
        tab:sim-b-samplesize       50    VAR + PrecFact                      Sens    0.115(0.084)    0.124(0.101)
        tab:sim-b-samplesize       50    VAR + PrecFact                      Spec    0.994(0.009)    0.994(0.008)
        tab:sim-b-samplesize       50    VAR + PrecFact                      Frob    1.774(0.220)    1.770(0.216)
        tab:sim-b-samplesize       50    VAR + PrecFact              Delta_F (SE)    0.230(0.039)    0.232(0.039)
        tab:sim-b-samplesize       50 Kastner (VAR+FSV)              Delta_F (SE)    0.028(0.034)    0.035(0.034)
        tab:sim-b-samplesize       50             VAPOR                      Frob    1.544(0.098)    1.537(0.107)
        tab:sim-b-samplesize      100 PrecFact (no VAR)                      Sens    0.170(0.081)    0.166(0.085)
        tab:sim-b-samplesize      100 PrecFact (no VAR)                      Spec    0.924(0.043)    0.929(0.044)
        tab:sim-b-samplesize      100 PrecFact (no VAR)                      Frob    2.584(0.558)    2.571(0.525)
        tab:sim-b-samplesize      100 PrecFact (no VAR)              Delta_F (SE)    1.246(0.111)    1.281(0.105)
        tab:sim-b-samplesize      100   Glasso (no VAR)              Delta_F (SE)    0.350(0.024)    0.398(0.025)
        tab:sim-b-samplesize      100      VAR + Glasso              Delta_F (SE)   -0.414(0.025)   -0.366(0.028)
        tab:sim-b-samplesize      100    VAR + PrecFact                      Sens    0.316(0.162)    0.313(0.134)
        tab:sim-b-samplesize      100    VAR + PrecFact                      Spec    0.990(0.013)    0.992(0.014)
        tab:sim-b-samplesize      100    VAR + PrecFact                      Frob    1.389(0.205)    1.362(0.176)
        tab:sim-b-samplesize      100    VAR + PrecFact              Delta_F (SE)    0.050(0.024)    0.072(0.019)
        tab:sim-b-samplesize      100 Kastner (VAR+FSV)              Delta_F (SE)   -0.187(0.024)   -0.138(0.028)
        tab:sim-b-samplesize      100             VAPOR                      Sens    0.181(0.124)    0.162(0.106)
        tab:sim-b-samplesize      100             VAPOR                      Spec    0.997(0.004)    0.998(0.006)
        tab:sim-b-samplesize      100             VAPOR                      Frob    1.339(0.143)    1.290(0.156)
        tab:sim-b-samplesize      200 PrecFact (no VAR)                      Sens    0.267(0.133)    0.263(0.132)
        tab:sim-b-samplesize      200 PrecFact (no VAR)                      Spec    0.901(0.045)    0.909(0.033)
        tab:sim-b-samplesize      200 PrecFact (no VAR)                      Frob    2.596(0.517)    2.610(0.520)
        tab:sim-b-samplesize      200 PrecFact (no VAR)              Delta_F (SE)    1.537(0.102)    1.581(0.099)
        tab:sim-b-samplesize      200   Glasso (no VAR)              Delta_F (SE)    0.624(0.027)    0.653(0.024)
        tab:sim-b-samplesize      200      VAR + Glasso              Delta_F (SE)   -0.421(0.030)   -0.392(0.029)
        tab:sim-b-samplesize      200    VAR + PrecFact                      Sens    0.636(0.132)    0.630(0.157)
        tab:sim-b-samplesize      200    VAR + PrecFact                      Frob    1.103(0.193)    1.085(0.169)
        tab:sim-b-samplesize      200    VAR + PrecFact              Delta_F (SE)    0.045(0.014)    0.055(0.008)
        tab:sim-b-samplesize      200 Kastner (VAR+FSV)              Delta_F (SE)   -0.160(0.020)   -0.131(0.022)
        tab:sim-b-samplesize      200             VAPOR                      Sens    0.575(0.145)    0.588(0.111)
        tab:sim-b-samplesize      200             VAPOR                      Spec    0.990(0.008)    0.992(0.008)
        tab:sim-b-samplesize      200             VAPOR                      Frob    1.059(0.175)    1.030(0.172)
        tab:sim-b-samplesize      500 PrecFact (no VAR)                      Sens    0.386(0.130)    0.383(0.146)
        tab:sim-b-samplesize      500 PrecFact (no VAR)                      Spec    0.861(0.065)    0.857(0.073)
        tab:sim-b-samplesize      500 PrecFact (no VAR)                      Frob    2.604(0.575)    2.593(0.587)
        tab:sim-b-samplesize      500 PrecFact (no VAR)              Delta_F (SE)    1.695(0.121)    1.681(0.123)
        tab:sim-b-samplesize      500   Glasso (no VAR)              Delta_F (SE)    0.765(0.034)    0.762(0.031)
        tab:sim-b-samplesize      500      VAR + Glasso              Delta_F (SE)   -0.505(0.034)   -0.508(0.031)
        tab:sim-b-samplesize      500    VAR + PrecFact                      Sens    0.735(0.124)    0.746(0.114)
        tab:sim-b-samplesize      500    VAR + PrecFact                      Spec    0.980(0.012)    0.977(0.011)
        tab:sim-b-samplesize      500    VAR + PrecFact                      Frob    0.959(0.175)    0.940(0.170)
        tab:sim-b-samplesize      500    VAR + PrecFact              Delta_F (SE)    0.049(0.025)    0.028(0.014)
        tab:sim-b-samplesize      500 Kastner (VAR+FSV)              Delta_F (SE)   -0.178(0.024)   -0.181(0.022)
        tab:sim-b-samplesize      500             VAPOR                      Sens    0.751(0.117)    0.752(0.106)
        tab:sim-b-samplesize      500             VAPOR                      Spec    0.982(0.009)    0.981(0.013)
        tab:sim-b-samplesize      500             VAPOR                      Frob    0.910(0.196)    0.913(0.177)
        tab:sim-b-samplesize     1000 PrecFact (no VAR)                      Sens    0.453(0.144)    0.450(0.137)
        tab:sim-b-samplesize     1000 PrecFact (no VAR)                      Spec    0.828(0.068)    0.833(0.072)
        tab:sim-b-samplesize     1000 PrecFact (no VAR)                      Frob    2.621(0.542)    2.623(0.548)
        tab:sim-b-samplesize     1000 PrecFact (no VAR)              Delta_F (SE)    1.741(0.110)    1.766(0.116)
        tab:sim-b-samplesize     1000   Glasso (no VAR)              Delta_F (SE)    0.802(0.029)    0.825(0.031)
        tab:sim-b-samplesize     1000      VAR + Glasso              Delta_F (SE)   -0.590(0.041)   -0.566(0.040)
        tab:sim-b-samplesize     1000    VAR + PrecFact                      Sens    0.773(0.097)    0.777(0.097)
        tab:sim-b-samplesize     1000    VAR + PrecFact                      Spec    0.974(0.015)    0.972(0.014)
        tab:sim-b-samplesize     1000    VAR + PrecFact                      Frob    0.900(0.194)    0.898(0.207)
        tab:sim-b-samplesize     1000    VAR + PrecFact              Delta_F (SE)    0.020(0.018)    0.041(0.019)
        tab:sim-b-samplesize     1000 Kastner (VAR+FSV)              Delta_F (SE)   -0.193(0.028)   -0.170(0.023)
        tab:sim-b-samplesize     1000             VAPOR                      Sens    0.766(0.105)    0.789(0.111)
        tab:sim-b-samplesize     1000             VAPOR                      Frob    0.880(0.206)    0.857(0.200)
                   tab:sim-c        1 PrecFact (no VAR)                      MSPE  0.9087(0.1253)  0.9093(0.1250)
                   tab:sim-c        1 PrecFact (no VAR) Delta_{mathrm{MSPE}} (SE) -0.0008(0.0042) -0.0004(0.0037)
                   tab:sim-c        1    VAR + PrecFact                      MSPE  0.9132(0.1265)  0.9041(0.1295)
                   tab:sim-c        1    VAR + PrecFact Delta_{mathrm{MSPE}} (SE)  0.0037(0.0039) -0.0055(0.0025)
                   tab:sim-c        1             VAPOR                      MSPE  0.9095(0.1275)  0.9097(0.1264)
                   tab:sim-c        1 Kastner (VAR+FSV) Delta_{mathrm{MSPE}} (SE) -0.0575(0.0104) -0.0577(0.0104)
                   tab:sim-c        1 Sample covariance Delta_{mathrm{MSPE}} (SE)  0.0095(0.0141)  0.0093(0.0141)
                   tab:sim-c       10 PrecFact (no VAR)                      MSPE  0.9057(0.1416)  0.9007(0.1437)
                   tab:sim-c       10 PrecFact (no VAR) Delta_{mathrm{MSPE}} (SE) -0.0022(0.0046) -0.0078(0.0031)
                   tab:sim-c       10    VAR + PrecFact                      MSPE  0.9027(0.1331)  0.8992(0.1325)
                   tab:sim-c       10    VAR + PrecFact Delta_{mathrm{MSPE}} (SE) -0.0052(0.0061) -0.0094(0.0062)
                   tab:sim-c       10             VAPOR                      MSPE  0.9079(0.1408)  0.9085(0.1413)
                   tab:sim-c       10 Kastner (VAR+FSV) Delta_{mathrm{MSPE}} (SE) -0.0676(0.0164) -0.0683(0.0171)
                   tab:sim-c       10 Sample covariance Delta_{mathrm{MSPE}} (SE) -0.0094(0.0209) -0.0100(0.0212)
                   tab:sim-c       80 PrecFact (no VAR)                      MSPE  0.7017(0.2131)  0.7066(0.2062)
                   tab:sim-c       80 PrecFact (no VAR) Delta_{mathrm{MSPE}} (SE) -0.1506(0.0326) -0.1462(0.0330)
                   tab:sim-c       80    VAR + PrecFact                      MSPE  0.8520(0.3533)  0.8531(0.3504)
                   tab:sim-c       80    VAR + PrecFact Delta_{mathrm{MSPE}} (SE) -0.0002(0.0108)  0.0003(0.0071)
                   tab:sim-c       80             VAPOR                      MSPE  0.8522(0.3335)  0.8527(0.3393)
                   tab:sim-c       80 Kastner (VAR+FSV) Delta_{mathrm{MSPE}} (SE) -0.1096(0.0179) -0.1101(0.0176)
                   tab:sim-c       80 Sample covariance Delta_{mathrm{MSPE}} (SE) -0.3854(0.0684) -0.3859(0.0690)
                   tab:sim-d      0.5 PrecFact (no VAR)                      Sens    0.117(0.054)    0.123(0.062)
                   tab:sim-d      0.5 PrecFact (no VAR)                      Spec    0.872(0.052)    0.888(0.048)
                   tab:sim-d      0.5 PrecFact (no VAR)                      Frob    2.673(0.721)    2.651(0.708)
                   tab:sim-d      0.5 PrecFact (no VAR)              Delta_F (SE)    2.166(0.146)    2.146(0.145)
                   tab:sim-d      0.5   Glasso (no VAR)              Delta_F (SE)    0.000(0.015)    0.002(0.015)
                   tab:sim-d      0.5      VAR + Glasso              Delta_F (SE)   -0.036(0.013)   -0.035(0.013)
                   tab:sim-d      0.5    VAR + PrecFact                      Sens    0.011(0.029)    0.010(0.027)
                   tab:sim-d      0.5    VAR + PrecFact                      Spec    0.999(0.003)    1.000(0.002)
                   tab:sim-d      0.5    VAR + PrecFact                      Frob    0.513(0.124)    0.512(0.124)
                   tab:sim-d      0.5    VAR + PrecFact              Delta_F (SE)    0.005(0.009)    0.006(0.009)
                   tab:sim-d      0.5 Kastner (VAR+FSV)              Delta_F (SE)    0.061(0.016)    0.063(0.015)
                   tab:sim-d      0.5             VAPOR                      Sens    0.012(0.030)    0.015(0.034)
                   tab:sim-d      0.5             VAPOR                      Frob    0.507(0.106)    0.506(0.109)
                   tab:sim-d        1 PrecFact (no VAR)                      Sens    0.170(0.066)    0.166(0.073)
                   tab:sim-d        1 PrecFact (no VAR)                      Spec    0.878(0.052)    0.884(0.053)
                   tab:sim-d        1 PrecFact (no VAR)                      Frob    2.876(0.700)    2.873(0.689)
                   tab:sim-d        1 PrecFact (no VAR)              Delta_F (SE)    2.198(0.136)    2.192(0.134)
                   tab:sim-d        1   Glasso (no VAR)              Delta_F (SE)    0.738(0.061)    0.736(0.061)
                   tab:sim-d        1      VAR + Glasso              Delta_F (SE)    0.139(0.022)    0.136(0.022)
                   tab:sim-d        1    VAR + PrecFact                      Sens    0.325(0.186)    0.329(0.185)
                   tab:sim-d        1    VAR + PrecFact                      Spec    0.997(0.004)    0.998(0.004)
                   tab:sim-d        1    VAR + PrecFact                      Frob    0.772(0.104)    0.763(0.102)
                   tab:sim-d        1    VAR + PrecFact              Delta_F (SE)    0.093(0.012)    0.082(0.011)
                   tab:sim-d        1 Kastner (VAR+FSV)              Delta_F (SE)    0.209(0.031)    0.207(0.032)
                   tab:sim-d        1             VAPOR                      Sens    0.354(0.172)    0.351(0.177)
                   tab:sim-d        1             VAPOR                      Spec    0.999(0.003)    0.998(0.003)
                   tab:sim-d        1             VAPOR                      Frob    0.679(0.077)    0.681(0.077)
                   tab:sim-d      1.5 PrecFact (no VAR)                      Sens    0.244(0.110)    0.238(0.111)
                   tab:sim-d      1.5 PrecFact (no VAR)                      Spec    0.863(0.065)    0.864(0.061)
                   tab:sim-d      1.5 PrecFact (no VAR)                      Frob    3.296(0.689)    3.269(0.656)
                   tab:sim-d      1.5 PrecFact (no VAR)              Delta_F (SE)    2.605(0.133)    2.580(0.126)
                   tab:sim-d      1.5   Glasso (no VAR)              Delta_F (SE)    1.503(0.091)    1.505(0.091)
                   tab:sim-d      1.5      VAR + Glasso              Delta_F (SE)    0.243(0.019)    0.246(0.019)
                   tab:sim-d      1.5    VAR + PrecFact                      Sens    0.573(0.188)    0.555(0.181)
                   tab:sim-d      1.5    VAR + PrecFact                      Frob    0.772(0.082)    0.762(0.082)
                   tab:sim-d      1.5    VAR + PrecFact              Delta_F (SE)    0.081(0.011)    0.073(0.012)
                   tab:sim-d      1.5 Kastner (VAR+FSV)              Delta_F (SE)    0.487(0.076)    0.489(0.076)
                   tab:sim-d      1.5             VAPOR                      Sens    0.557(0.162)    0.550(0.169)
                   tab:sim-d      1.5             VAPOR                      Frob    0.691(0.069)    0.689(0.071)
                   tab:sim-d      2.5 PrecFact (no VAR)                      Sens    0.354(0.136)    0.333(0.147)
                   tab:sim-d      2.5 PrecFact (no VAR)                      Spec    0.849(0.059)    0.849(0.049)
                   tab:sim-d      2.5 PrecFact (no VAR)                      Frob    3.888(0.747)    3.882(0.777)
                   tab:sim-d      2.5 PrecFact (no VAR)              Delta_F (SE)    3.164(0.142)    3.158(0.146)
                   tab:sim-d      2.5      VAR + Glasso              Delta_F (SE)    0.244(0.014)    0.244(0.013)
                   tab:sim-d      2.5    VAR + PrecFact                      Sens    0.751(0.106)    0.751(0.101)
                   tab:sim-d      2.5    VAR + PrecFact                      Spec    0.998(0.005)    0.999(0.002)
                   tab:sim-d      2.5    VAR + PrecFact                      Frob    0.708(0.079)    0.702(0.066)
                   tab:sim-d      2.5    VAR + PrecFact              Delta_F (SE)   -0.017(0.019)   -0.022(0.019)
                   tab:sim-d      2.5 Kastner (VAR+FSV)              Delta_F (SE)    1.105(0.129)    1.105(0.130)
                   tab:sim-d      2.5             VAPOR                      Sens    0.725(0.110)    0.718(0.114)
                   tab:sim-d      2.5             VAPOR                      Spec    1.000(0.002)    1.000(0.001)
                   tab:sim-d        4 PrecFact (no VAR)                      Sens    0.394(0.142)    0.400(0.146)
                   tab:sim-d        4 PrecFact (no VAR)                      Spec    0.828(0.056)    0.835(0.060)
                   tab:sim-d        4 PrecFact (no VAR)                      Frob    4.391(0.806)    4.387(0.859)
                   tab:sim-d        4 PrecFact (no VAR)              Delta_F (SE)    3.510(0.157)    3.500(0.167)
                   tab:sim-d        4   Glasso (no VAR)              Delta_F (SE)    2.889(0.118)    2.882(0.117)
                   tab:sim-d        4      VAR + Glasso              Delta_F (SE)    0.090(0.023)    0.083(0.023)
                   tab:sim-d        4    VAR + PrecFact                      Sens    0.837(0.088)    0.829(0.090)
                   tab:sim-d        4    VAR + PrecFact                      Spec    0.996(0.007)    0.997(0.006)
                   tab:sim-d        4    VAR + PrecFact                      Frob    0.691(0.082)    0.680(0.080)
                   tab:sim-d        4    VAR + PrecFact              Delta_F (SE)   -0.189(0.032)   -0.207(0.032)
                   tab:sim-d        4 Kastner (VAR+FSV)              Delta_F (SE)    1.573(0.150)    1.566(0.148)
                   tab:sim-d        4             VAPOR                      Sens    0.783(0.092)    0.778(0.102)
                   tab:sim-d        4             VAPOR                      Frob    0.880(0.135)    0.887(0.137)
                   tab:sim-e    50,50        VAPOR (DL)             |bar R-R_0|_F      3.97(0.41)      3.95(0.41)
                   tab:sim-e    50,50        VAPOR (DL)                      Sens      0.29(0.07)      0.30(0.08)
                   tab:sim-e    50,50        VAPOR (DL)                      Spec      0.94(0.01)      0.95(0.01)
                   tab:sim-e    50,50     VAPOR (ridge)             |bar R-R_0|_F      3.96(0.41)      3.96(0.40)
                   tab:sim-e    50,50     VAPOR (ridge)     |hat B-B_0|_F/|B_0|_F      1.10(0.55)      1.09(0.54)
                   tab:sim-e    50,50     VAPOR (ridge)                      Sens      0.28(0.07)      0.29(0.07)
                   tab:sim-e    50,50     VAPOR (ridge)                      Spec      0.94(0.01)      0.95(0.01)
                   tab:sim-e    50,75        VAPOR (DL)             |bar R-R_0|_F      3.26(0.50)      3.25(0.50)
                   tab:sim-e    50,75     VAPOR (ridge)             |bar R-R_0|_F      3.24(0.48)      3.24(0.49)
                   tab:sim-e    50,75     VAPOR (ridge)     |hat B-B_0|_F/|B_0|_F      0.59(0.31)      0.58(0.30)
                   tab:sim-e    50,75     VAPOR (ridge)                      Sens      0.44(0.09)      0.44(0.10)
                   tab:sim-e   50,100        VAPOR (DL)             |bar R-R_0|_F      2.68(0.43)      2.69(0.43)
                   tab:sim-e   50,100     VAPOR (ridge)             |bar R-R_0|_F      2.62(0.42)      2.63(0.44)
                   tab:sim-e  100,100        VAPOR (DL)             |bar R-R_0|_F      6.88(1.31)      6.85(1.30)
                   tab:sim-e  100,100        VAPOR (DL)     |hat B-B_0|_F/|B_0|_F      0.64(0.07)      0.64(0.08)
                   tab:sim-e  100,100        VAPOR (DL)                      Sens      0.46(0.08)      0.47(0.09)
                   tab:sim-e  100,100     VAPOR (ridge)             |bar R-R_0|_F      6.98(1.26)      6.93(1.22)
                   tab:sim-e  100,100     VAPOR (ridge)     |hat B-B_0|_F/|B_0|_F      0.93(0.33)      0.94(0.33)
                   tab:sim-e  100,100     VAPOR (ridge)                      Sens      0.46(0.09)      0.47(0.08)
                   tab:sim-e  150,150        VAPOR (DL)             |bar R-R_0|_F      8.17(1.42)      8.10(1.30)
                   tab:sim-e  150,150        VAPOR (DL)                      Sens      0.57(0.06)      0.58(0.06)
                   tab:sim-e  150,150     VAPOR (ridge)             |bar R-R_0|_F      8.06(1.31)      8.12(1.37)
                   tab:sim-e  150,150     VAPOR (ridge)     |hat B-B_0|_F/|B_0|_F      0.79(0.26)      0.78(0.25)
                   tab:sim-e  150,150     VAPOR (ridge)                      Sens      0.58(0.05)      0.58(0.06)
                   tab:sim-f    50,35        VAPOR (DL)             |bar R-R_0|_F      4.37(0.38)      4.36(0.36)
                   tab:sim-f    50,35        VAPOR (DL)     |hat B-B_0|_F/|B_0|_F      0.83(0.08)      0.84(0.08)
                   tab:sim-f    50,35        VAPOR (DL)                      Sens      0.19(0.05)      0.20(0.06)
                   tab:sim-f    50,35     VAPOR (ridge)             |bar R-R_0|_F      4.32(0.37)      4.32(0.38)
                   tab:sim-f    50,35     VAPOR (ridge)     |hat B-B_0|_F/|B_0|_F      0.90(0.26)      0.90(0.27)
                   tab:sim-f    50,35     VAPOR (ridge)                      Sens      0.19(0.05)      0.20(0.06)
                   tab:sim-f   100,70        VAPOR (DL)             |bar R-R_0|_F      8.98(1.04)      8.95(1.05)
                   tab:sim-f   100,70        VAPOR (DL)     |hat B-B_0|_F/|B_0|_F      0.80(0.06)      0.79(0.05)
                   tab:sim-f   100,70        VAPOR (DL)                      Sens      0.23(0.05)      0.23(0.04)
                   tab:sim-f   100,70     VAPOR (ridge)             |bar R-R_0|_F      8.96(1.05)      8.98(1.06)
                   tab:sim-f  150,105        VAPOR (DL)             |bar R-R_0|_F     13.50(1.54)     13.47(1.55)
                   tab:sim-f  150,105        VAPOR (DL)     |hat B-B_0|_F/|B_0|_F      0.77(0.06)      0.77(0.05)
                   tab:sim-f  150,105        VAPOR (DL)                      Sens      0.23(0.05)      0.24(0.05)
                   tab:sim-f  150,105     VAPOR (ridge)             |bar R-R_0|_F     13.52(1.52)     13.56(1.54)
                   tab:sim-f  150,105     VAPOR (ridge)     |hat B-B_0|_F/|B_0|_F      0.72(0.10)      0.72(0.11)
```

## Simulation run provenance

sim_g's three settings were run as independent processes; sim_g.R re-seeds from the replicate index alone, so per-setting runs reproduce the sequential result exactly (verified: rerunning p = 50 gave results identical at tolerance 0). No floor-hit counts were recorded for this simulation run: the driver that logged them was stopped before writing simulations_floor_report.rds, and the per-setting runner does not wrap vapor().

*Source: run logs*

```
      simulation                          status         finished minutes
       A (sim_a)                        complete 2026-10-20 17:59      NA
       B (sim_b)                        complete 2026-10-20 20:48      NA
       C (sim_c)                        complete 2026-10-20 21:01      NA
       D (sim_d)                        complete 2026-10-21 14:31      NA
       E (sim_e)                        complete 2026-10-21 20:47      NA
  F (sim_f) p=50                        complete 2026-10-22 09:31    26.7
 F (sim_f) p=100                        complete 2026-10-21 16:48    74.3
 F (sim_f) p=150                        complete 2026-10-21 20:24   290.5
 F (sim_f) p=200 CANCELLED at replicate 19 of 25                -      NA
```
