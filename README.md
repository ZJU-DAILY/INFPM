<h3>
	<center>Influence Persistence Maximization in Temporal Social Networks</center>
</h3>

### Introduction

1. This repository contains the full version of our paper.
2. This repository contains the codes and datasets used in our paper.

### Datasets

We use 7 publicly available real-world temporal networks, including CollegeMsg, EmailCore, MathOverflow, AskUbuntu, SuperUser, WikiTalk,  and StackOverflow.

All the datasets are available on SNAP [1].

*[1] Jure Leskovec and Andrej Krevl. 2014. SNAP Datasets: Stanford large network dataset collection. http://snap.stanford.edu/data.*


### Algorithms

All methods are implemented in C++ and executed on a server with Intel(R) Core(TM) 3.70GHz CPU and 128GB of RAM.

1. RevG, LRep: ReverseGreedy and LazyReplace algorithms that utilize the reverse influence sampling (RIS) method [2] to
estimate influence spread.

```shell
./IPMAX -dataset=<dataset root> -R=<number of RR sets> -k=<number of seeds (percentage of number of nodes)> -theta=<persistence threshold> -alpha=<influence threshold> -alg=<RG/LR>
```

*[2] Christian Borgs, Michael Brautbar, Jennifer Chayes, and Brendan Lucier. 2014. Maximizing social influence in nearly optimal time. In Proceedings of the twenty-fifth annual ACM-SIAM symposium on Discrete algorithms. SIAM, 946–957.*


2. RevG+, LRep+: ReverseGreedy and LazyReplace algorithms that employ the influence spread computation method proposed in this paper.

```shell
./IPMAX -dataset=<dataset root> -R=<number of SCC-RR sets> -RSCC=50 -k=<number of seeds (percentage of number of nodes)> -theta=<persistence threshold> -alpha=<influence threshold> -alg=<RGSCC/LRSCC>
```

3. WRG, WLR: ReverseGreedy and LazyReplace algorithms that employ PssW to compute window-based influence persistence.

```shell
./IPMAX -dataset=<dataset root> -R=<number of RR sets> -W=<window size> -k=<number of seeds (percentage of number of nodes)> -theta=<persistence threshold> -alpha=<influence threshold> -alg=<RG/LR>
```


### Running Environment

A 64-bit Linux-based OS. 

GCC 4.7.2 and later.
