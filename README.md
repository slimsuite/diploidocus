# Diploidocus: Diploid genome assembly analysis tools

This module codes for a number of diploid genome assembly manipulation and analysis functions. The different run
modes are set using `runmode=X`:

* `sortnr` performs an all-by-all mapping with minimap2 and then removes redundancy
* `diphap` splits a pseudodiploid assembly into primary and alternative scaffolds
* `diphapnr` runs `sortnr` followed by `diphap`
* `purgehap` [coming soon!] filters scaffolds based on post-processing of purge_haplotigs
* `vecscreen` [coming soon!] searches for contaminants and flags/masks/removes identified scaffolds
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

For `sortnr` and `diphapnr` mode, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
environment `$PATH` or given to Diploidocus with the `minimap2=PROG` setting.


