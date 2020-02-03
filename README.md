# Diploidocus: Diploid genome assembly analysis tools

Diploidocus is a sequence analysis toolkit for a number of different analyses related to diploid genome assembly.
The main suite of analyses combines long read depth profiles, short read kmer analysis, assembly kmer analysis,
BUSCO gene prediction and contaminant screening for a number of assembly tasks including genome size prediction,
contamination identification, haplotig identification/removal and low quality contig/scaffold filtering. In addition,
Diploidocus has functions for removing redundancy, generating a non-redundant pseudo-diploid assembly with primary
and secondary scaffolds from 10x pseudohap output, and creating an _in silico_ diploid set of long reads from two
haploid parents (for testing phasing etc.).

Please note that Diploidocus is still in development and documentation is currently a bit sparse.

The different run modes are set using `runmode=X`:

* `diploidocus` default run mode will run `gensize`, `telomeres`, `vecscreen` and `purgehap` analysis
* `gensize` uses BUSCO results, a BAM file and read file(s) to predict the genome size of the organism
* `purgehap` filters scaffolds based on post-processing of purge_haplotigs
* `telomeres` performs a regex telomere search based on method of https://github.com/JanaSperschneider/FindTelomeres
* `vecscreen` searches for contaminants and flags/masks/removes identified scaffolds
* `sortnr` performs an all-by-all mapping with minimap2 and then removes redundancy
* `diphap` splits a pseudodiploid assembly into primary and alternative scaffolds
* `diphapnr` runs `sortnr` followed by `diphap`
* `insilico` generates balanced diploid combined reads from two sequenced haploid parents

See <https://slimsuite.github.io/diploidocus/> for details of each mode. General SLiMSuite run documentation can be
found at <https://github.com/slimsuite/SLiMSuite>.

**NOTE:** Diploidocus is under development and documentation might be a bit sparse. Please contact the author or
post an issue on GitHub if you have any questions.

Diploidocus is available as part of SLiMSuite, or via a standalone GitHub repo at
<https://github.com/slimsuite/diploidocus>.

---

## Running Diploidocus

Diploidocus is written in Python 2.x and can be run directly from the commandline:

    python $CODEPATH/diploidocus.py [OPTIONS]

If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
directory. If running from the standalone [Diploidocus git repo](https://github.com/slimsuite/diploidocus), `$CODEPATH`
will be the path the to `code/` directory. Please see details in the [Diploidocus git repo](https://github.com/slimsuite/diploidocus)
for running on example data.

For most modes, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
environment `$PATH` or given to Diploidocus with the `minimap2=PROG` setting.


