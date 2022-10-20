# TreeTerminus

TreeTerminus is a program for grouping transcripts into tree structures for an RNASeq experiment, taking as input [Salmon](https://github.com/COMBINE-lab/salmon) quantified RNA-Seq files. The leaves of the trees represent individual transcripts and internal nodes represent an aggregation of a transcript set. Across the samples, the uncertainty associated with the abundance estimate of the node decreases on average, ascending the tree.

### Building TreeTerminus
TreeTerminus is implemented in [Rust](https://www.rust-lang.org/) and is built on top of [Terminus](https://github.com/COMBINE-lab/terminus) codebase. Terminus uses the [cargo](https://github.com/rust-lang/cargo) build system and package manager.  To build terminus from source, you will need to have rust (ideally v1.40 or greater) installed.  To build TreeTerminus, follow the below steps:

```
$ cargo build --release
```



