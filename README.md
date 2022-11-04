# TreeTerminus

TreeTerminus is a program for grouping transcripts into distinct trees for an RNASeq experiment, taking as input [Salmon](https://github.com/COMBINE-lab/salmon) quantified RNA-Seq files. The leaves of the trees represent individual transcripts and internal nodes represent an aggregation of a transcript set. Across the samples, the uncertainty associated with the abundance estimate of the node decreases on an average, ascending the trees.

The transcript trees can be obtained either for each sample or across all the samples in an RNA-Seq experiment. To obtain the summarized trees over all the samples in the RNA-Seq experiment, we provide two modes - **Mean** and **Consensus**. The **Mean** trees are obtained by using mean of reduction in inferential variance across all the samples. To obtain the **Consensus** trees, first the trees are obtained for each sample and then after processing, are passed as an input to the consensus tree algorithm. [PHYLIP](https://evolution.genetics.washington.edu/phylip.html)'s implementation of the majority rule extended tree consensus algorithm has been used.


## Building TreeTerminus
TreeTerminus is implemented in [Rust](https://www.rust-lang.org/) and is built on top of [Terminus](https://github.com/COMBINE-lab/terminus) codebase. TreeTerminus uses the [cargo](https://github.com/rust-lang/cargo) build system and package manager.  You will need to have rust (ideally v1.54 or greater) installed. To build TreeTerminus, follow the steps below:

```
git clone git@github.com:COMBINE-lab/TreeTerminus.git
cd TreeTerminus
cargo build --release
```

To check if the build happens successfully, from the top-level directory of the repository run:
```
target/release/treeterminus
``` 
The above command should not produce any errors and provide the usage information.

## Using TreeTerminus
TreeTerminus provides two sub-commands - `group` and `consensus`.

The `group` step is used to generate transcript trees either for each sample or the **Mean** trees across samples for an RNASeq experiment. The `consensus` step outputs a set of **Consensus** trees for the RNASeq experiment using the majority rule extended consensus tree algorithm. It requires that the `group` step should be run individually on each sample in the experiment.

### Group
To obtain the transcript trees for a single sample in an RNA-Seq experiment, run `group` from the parent directory of `TreeTerminus` as:

```
target/release/treeterminus group -d <salmon_dir> -o <out_dir> --mean_inf false
```
The option `-d` denotes the path to `salmon` quantified sample directory, `-o` denotes the path where `TreeTerminus` output will be directed, `--mean_inf` is a boolean indicating whether **Mean** tree will be constructed. The final trees will be stored in the file `group_nwk.txt`, inside a subdirectory in the `out_dir` directory that belongs to a sample in the RNA-Seq experiment.

To obtain the **Mean** trees for an an RNA-Seq experiment, run `group` from the parent directory of `TreeTerminus` as:
```
target/release/treeterminus group -d <salmon_dir> -o <out_dir> --mean_inf true
```

Here to `-d` argument, provide the directory that contains all the `salmon` quantified samples of interest in the RNA-Seq experiment, rather than just a single sample. The final trees will be stored in the file `cluster_nwk.txt`, inside the `out_dir` directory.

The information about the other arguments that can be provided to `group`, can be obtained by running:
```
target/release/treeterminus group -h
```

### Consensus
To obtain consensus trees for the RNA-Seq Experiment, run `consensus` from the parent directory of `TreeTerminus` as: 
```
target/release/treeterminus consensus -d <salmon_dir> -o <out_dir> 
```
The inputs to the `-d` and `-o` flags are the same as those defined for the `group` step that was used to generate the **Mean** tree.

The information about the other arguments that can be provided to `consensus`, can be obtained by running:
```
target/release/treeterminus consensus -h
```

**Note** - A limitation of the `consensus` mode only one instance of it can be run at any given moment, aka do not run this on two or more experiments simultaneously.

### Example:
Let us assume the following directory structure of a parent directory:
- TreeTerminus
- SalmonQuant
    - SampleA
    - SampleB
    - SampleC
- TreeTermOut

`SalmonQuant` contains the salmon quantified files for an RNA-Seq experiment, with the individual sample folders inside it being `SampleA`, `SampleB` and `SampleC`. `TreeTermOut` is where we will keep the output of `TreeTerminus`.  

To run `TreeTerminus`, first move to its directory as:
```
cd TreeTerminus
```

To run `group` step on a single sample such as `SampleA`, run:
```
target/release/treeterminus group -d ../SalmonQuant/SampleA -o ../TreeTermOut --mean_inf false
```

To obtain the **Mean** trees for the experiment, run:     
```
target/release/treeterminus group -d ../SalmonQuant -o ../TreeTermOut --mean_inf true
```

To obtain the **Consensus** trees for the experiment, run:     
```
target/release/treeterminus group -d ../SalmonQuant -o ../TreeTermOut
```

### Cite


### Authors
Noor Pratap Singh, Michael Love, Rob Patro

### Downstream Analysis
We have also created an R Package [beaveR](https://github.com/NPSDC/beaveR) that parses the output of `TreeTerminus`, implements DP algorithms for finding an optimal cut, and provides helper functions to obtain
useful statistics for subtrees within the TreeTerminus-derived tree structures.

### Reproducing the experiments in the paper
The code to reproduce the experiments and figures in the paper are provided at https://github.com/NPSDC/treeterm-paper.