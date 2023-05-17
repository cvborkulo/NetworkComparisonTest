# NCT works when directly inputting data

    Code
      summary(NCT_a)
    Output
       INDEPENDENT GROUPS BINARY NETWORK COMPARISON TEST 
      
       P-VALUE CORRECTION: none 
      
       NETWORK INVARIANCE TEST 
       Test statistic M: 0.5530656 
       p-value 0.5454545 
      
       GLOBAL STRENGTH INVARIANCE TEST 
       Global strength per group:  8.27092 7.016815 
       Test statistic S:  1.254105 
       p-value 0.2727273
      
       EDGE INVARIANCE TEST 
        Var1 Var2 p-value Test statistic E
      1   V1   V2       1                0
      2   V3   V6       1                0

# NCT works with output from estimateNetwork

    Code
      summary(NCT_b)
    Output
       INDEPENDENT GROUPS GAUSSIAN NETWORK COMPARISON TEST 
      
       P-VALUE CORRECTION: none 
      
       NETWORK INVARIANCE TEST 
       Test statistic M: 0.5530656 
       p-value 0.5454545 
      
       GLOBAL STRENGTH INVARIANCE TEST 
       Global strength per group:  8.27092 7.016815 
       Test statistic S:  1.254105 
       p-value 0.2727273
      
       EDGE INVARIANCE TEST 
        Var1 Var2 p-value Test statistic E
      1   V1   V2       1                0
      2   V3   V6       1                0

# NCT works when testing node strength

    Code
      summary(NCT_c)
    Output
       DEPENDENT GROUPS GAUSSIAN NETWORK COMPARISON TEST 
      
       P-VALUE CORRECTION: none 
      
       NETWORK INVARIANCE TEST 
       Test statistic M: 0.5530656 
       p-value 0.6363636 
      
       GLOBAL EXPECTED INFLUENCE INVARIANCE TEST 
       Global EI per group:  3.209379 1.651196 
       Test statistic S:  1.558182 
       p-value 0.09090909
      
       EDGE INVARIANCE TEST 
        Var1 Var2 p-value Test statistic E
      1   V1   V2       1                0
      2   V3   V6       1                0
      
      
       CENTRALITY INVARIANCE TEST 
       Nodes tested: V1 V2 V3 V4 V5 V6 
       Centralities tested: strength
       Test statistics C: 
           strength
      V1 -0.1520385
      V2  0.3358209
      V3  0.2329733
      V4  0.9527754
      V5  0.4287183
      V6  0.7099606
      
       p-values: 
          strength
      V1 0.8181818
      V2 0.7272727
      V3 0.4545455
      V4 0.1818182
      V5 0.3636364
      V6 0.3636364

# NCT works when testing expected influence

    Code
      summary(NCT_d)
    Output
       DEPENDENT GROUPS GAUSSIAN NETWORK COMPARISON TEST 
      
       P-VALUE CORRECTION: none 
      
       NETWORK INVARIANCE TEST 
       Test statistic M: 0.5530656 
       p-value 0.6363636 
      
       GLOBAL EXPECTED INFLUENCE INVARIANCE TEST 
       Global EI per group:  3.209379 1.651196 
       Test statistic S:  1.558182 
       p-value 0.09090909
      
       EDGE INVARIANCE TEST 
        Var1 Var2 p-value Test statistic E
      1   V1   V2       1                0
      2   V3   V6       1                0
      
      
       CENTRALITY INVARIANCE TEST 
       Nodes tested: V1 V2 V3 V4 V5 V6 
       Centralities tested: expectedInfluence
       Test statistics C: 
         expectedInfluence
      V1         0.1520385
      V2         0.3358209
      V3         0.2329733
      V4         1.3298653
      V5         0.4606148
      V6         0.6050513
      
       p-values: 
         expectedInfluence
      V1        0.81818182
      V2        0.81818182
      V3        0.45454545
      V4        0.09090909
      V5        0.72727273
      V6        0.27272727

