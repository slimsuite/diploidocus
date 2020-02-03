#!/usr/bin/python

# See below for name and description
# Copyright (C) 2017 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       Diploidocus
Description:  Diploid genome assembly analysis toolkit
Version:      0.7.0
Last Edit:    02/02/20
Copyright (C) 2017  Richard J. Edwards - See source code for GNU License Notice

Function:
    Diploidocus is a sequence analysis toolkit for a number of different analyses related to diploid genome assembly.
    The main suite of analyses combines long read depth profiles, short read kmer analysis, assembly kmer analysis,
    BUSCO gene prediction and contaminant screening for a number of assembly tasks including genome size prediction,
    contamination identification, haplotig identification/removal and low quality contig/scaffold filtering. In addition,
    Diploidocus has functions for removing redundancy, generating a non-redundant pseudo-diploid assembly with primary
    and secondary scaffolds from 10x pseudohap output, and creating an in silico diploid set of long reads from two
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

    Diploidocus is available as part of SLiMSuite, or via a standalone GitHub repo at
    <https://github.com/slimsuite/diploidocus>.

Run Modes:
    ### ~ Main Diploidocus filtering [runmode=diploidocus] ~ ###
    Diploidocus builds on the PurgeHaplotigs classifications to use the depth bins, KAT assembly kmer frequencies and
    BUSCO results to reclassify scaffolds and partition them into:

    * `*.diploidocus.fasta` = the scaffolds kept for the next round of PurgeHap
    * `*.core.fasta` = the same set of scaffolds, minus repeats
    * `*.junk.fasta` = low coverage scaffolds, removed as junk
    * `*.quarantine.fasta` = putative haplotigs, removed from the assembly but retained for reference.

    **NOTE:** PurgeHaplotigs was not using zero coverage bases in its percentages. This is now fixed by Diploidocus.

    _Further details coming soon!_

    ---
    ### ~ Cycled Diploidocus filtering [runmode=dipcycle] ~ ###
    Diploidocus can be automated to cycle through repeated rounds of the main purge_haplotigs/filtering method until
    no further scaffolds are removed. Each cycle will run Diploidocus with a numerical suffix, e.g. `$BASEFILE.1.*`,
    using the `*.diploidocus.fasta` output from the previous cycle as input. Cycling will continue until no further
    scaffolds are filtered into either `*.quarantine.fasta` or `*.junk.fasta`.

    Final outputs from the final cycle will then be copied under the original `$BASEFILE` prefix:

    * `$BASEFILE.diploidocus.tdt` = Final ratings for the final set of scaffolds. (See earlier cycles for purged.)
    * `$BASEFILE.diploidocus.fasta` = the scaffolds kept from the final Diploidocus cycle
    * `$BASEFILE.core.fasta` = the same set of scaffolds, minus repeats
    * `$BASEFILE.quarantine.fasta` = concatenated purged scaffolds from all Diploidocus cycles.
    * `$BASEFILE.junk.fasta` = concatenated low coverage and low quality scaffolds, removed as junk, from all cycles.

    ---
    ### ~ Genome size prediction [runmode=gensize] ~ ###

    ---
    ### ~ Running Purge_haplotigs using BUSCO-guided cutoffs [runmode=purgehap] ~ ###
    _Coming soon!_

    ---
    ### ~ Telomere finding [runmode=telomere] ~ ###
    _Details coming soon!_

    ---
    ### ~ Vector/contamination screening [runmode=vecscreen] ~ ###
    _Details coming soon!_

    ---
    ### ~ Sorted non-redundant assembly cleanup [runmode=sortnr] ~ ###

    The sorted non-redundant assembly cleanup mode (`runmode=sortnr`) screens out any sequences that are 100% gap,
    then removes any sequences that are 100% redundant with other sequences in the input. This includes full and
    partial matches, i.e. if sequence X is wholly contained within sequence Y then X will be removed.

    First, sequences are loaded from the file given with `seqin=FILE` and any [rje_seqlist](http://rest.slimsuite.unsw.edu.au/seqlist)
    filters and sequence sorting are applied to the input. Sequences that are 100% Ns are removed and any gaps
    exceeding 10 nt are reduced to 10 `N`s (`NNNNNNNNNN`) to prevent minimap2 from splitting sequences on long gaps.
    These gap-reduced sequences are output to `$BASEFILE.tmp.fasta` and used for an all-by-all minimap2 search.

    By default, minimap2 is run with the options to generate a `$BASEFILE.tmp.paf` file:

        --cs -p 0.0001 -t 4 -x asm20 -N 250

    To modify minimap2 search settings, please see the [rje_paf](http://rest.slimsuite.unsw.edu.au/rje_paf)
    documentation.

    **NOTE:** These run options can probably be made more stringent to speed up minimap2 without loss of function.
    Future releases may alter defaults accordingly.

    Minimap2 output is parsed to identify scaffold-scaffold matches. Self-hits are ignored.
    The minimum (gap-reduced) sequence length is used as a rapid parsing filter: any minimap2 matches that are less
    than 95% of the query sequence (`Length`+`nn` fields) or less that 100% identity (`Identity`+`nn`)/(`Length`+`nn`)
    are filtered during parsing.

    **NOTE:** Future releases may feature an option to reduce the global percentage identity cut-off. Please contact
    the author if you wish to see this implemented.

    Minimap2 hits are then processed reverse-sorted by reference sequence size (e.g. scaffold length). Any hits
    where either sequence has already been filtered are skipped. Otherwise, if the match (as determined by the
    length of `:` regions in the CS string) matches the query length, the Query sequence will be flagged for
    remove as "identical to" or "contained within" the Hit. (Mutually partial overlapping exact matches are NOT
    filtered.) Filtered IDs and their matches are output to `$BASEFILE.redundant.txt`.

    Once all sequences have been filtered, the remaining sequences are output to: `$BASEFILE.nr.fasta`.

    **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
    add the `sortseq=invsize` command to the Diploidocus run command.

    Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to:
    `$BASEFILE.summarise.tdt`.

    Temporary gap-reduced and minimap2 PAF files are deleted unless running in `debug` or
    `dev` modes.

    ---
    ### ~ Pseudodiploid to primary and alternative haploptigs [runmode=diphap(nr)] ~ ###

    This protocol is based on 10x assemblies made for multiple organisms with supernova v2.0.0 and supernova v2.1.1.
    In each case, some redundancy in output was discovered (a) within pseudohap output, and (b) in terms of fully
    homozygous (identical) scaffolds shared by both haplotigs. It was also not entirely clear on what basis a
    particular haplotype was assigned to pseudohap1 or pseudohap2.

    The general workflow therefore sought to remove redundancy, generate a set of primary scaffolds based on scaffold
    length, and generate a non-redundant set of alternative scaffolds where heterozygosity exists. If `diphapnr` mode
    is used, the full workflow is implement by first running the `sortnr` workflow described above. In the reduced
    `diphap` mode, redundancy is not removed first.

    Sequences are loaded and matching haplotigs identified based on their names. Sequence names MUST end `HAP(\d+)`,
    where `(\d+)` indicates an integer that matches up haplotigs (as produced by supernova pseudohap2 output, for
    example). This is **not** a pipeline to identify haplotig pairs, it is purely for splitting identified
    haplotigs into primary and alternative assemblies.

    Processing itself is quite simple. Haplotig pairs are identified based on matching `HAP(\d+)` numbers. Where a
    single haplotig is found, it is assigned as `diploid`, under the assumption that the two haplotigs were identical
    and one was removed. (It is possible that only one parent had this scaffold, e.g. sex chromosomes, so some post-
    processing of descriptions may be required.) If two haplotigs with the same number are identified, the longest
    is assigned to `haploidA` and the shorter `haploidB`.

    The **Primary Assemmbly** is then compiled from all `haploidA` and `diploid` sequences. These are given `pri`
    prefixes and output to `$BASEFILE.pri.fasta`. The **Alternative** comprised of all `haploidB` sequences is output
    to `$BASEFILE.alt.fasta`. If redundancy has been removed, this will likely be a subset of the full assembly. The
    combined set of all primary and alternative sequences is output to `$BASEFILE.dipnr.fasta`.

    **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
    add the `sortseq=invsize` command to the Diploidocus run command.

    Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to
    `$BASEFILE.summarise.tdt`:

    * `$BASEFILE.dipnr.fasta` = Combined pseudodiploid with `haploidA`, `haploidB` and `diploid` annotation.
    * `$BASEFILE.pri.fasta` = Primary assembly with `haploidA` and `diploid` sequences.
    * `$BASEFILE.alt.fasta` = Alternative assembly with `haploidB` sequences.


    ---
    ### ~ In silico diploid generator [runmode=insilico] ~ ###

    This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
    parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
    parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
    identifier table.)

    A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
    unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
    selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
    two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
    This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
    no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
    subreads will be added from the other parent if they reduce the difference in cumulative output for each parent, or
    until `lenfilter=X` is reached.

    Final output will be a `*.LXXXRQXX.fasta` file in which each parent has a similar total sequence content and for
    which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
    assemblies, where one parent has higher quality data than the other.

    NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
    higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
    minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
    relaxed. Similarly, only sequences above `lenfilter=X` in length will be output. These are the figures given in the
    `LXXXRQXX` part of the output file, e.g. defaults of RQ>=0.84 and Len>=500 generates `*.L500RQ84.fas`.


Commandline:
    ### ~ Main Diploidocus run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    runmode=X       : Diploidocus run mode [sortnr/diphap/diphapnr/purgehap/telomere/vecscreen/insilico/gensize]
    basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
    summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    genomesize=INT  : Haploid genome size (bp) [0]
    dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
    tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
    ### ~ Genome size prediction & Purge haplotigs options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb file types matching reads for minimap2 mapping [ont]
    readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
    quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    kmerreads=FILELIST : File of high quality reads for KAT kmer analysis []
    10xtrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer Reads data for KAT analysis [False]
    scdepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
    minmedian=INT   : Minimum median depth coverage to avoid low coverage filter [3]
    phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
    phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
    phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
    zeroadjust=T/F  : Add zero coverage bases to purge_haplotigs LowPerc and adjust total [True]
    includegaps=T/F : Whether to include gaps in the zero coverage bases for adjustment (see docs) [False]
    mingap=INT      : Minimum length of a stretch of N bases to count as a gap for exclusion [10]
    purgemode=X     : Rules used for purgehap analysis (simple/complex/nala) [complex]
    ### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    telofwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
    telorev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]
    telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
    teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
    ### ~ VecScreen options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    screendb=FILE   : File of vectors/contaminants to screen out using blastn and VecScreen rules []
    screenmode=X    : Action to take following vecscreen searching (report/mask/trim/purge) [report]
    minvechit=INT   : Minimum length for a screendb match [50]
    efdr=NUM        : Expected FDR threshold for VecScreen queries (0 is no filter) [1.0]
    ### ~ SortNR filtering/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    checkcov=PERC   : Percentage coverage for double-checking partial exact matches [95]
    seqout=FILE     : Output sequence assembly [$BASEFILE.nr.fasta]
    ### ~ In silico diploid input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    rqfilter=X      : Minimum RQ for output subreads [0]
    lenfilter=X     : Min read length for filtered subreads [500]
    parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, subprocess, sys, time, shutil
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_forker, rje_obj, rje_rmd, rje_seqlist, rje_sequence, rje_paf #, rje_genomics
import rje_blast_V2 as rje_blast
import smrtscape
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Fixed bugs with parent basefile, genome size default and Sequel data parsing.
    # 0.2.0 - Added sortnr minimap2 run mode, and diphap primary/alternative assembly splitting.
    # 0.3.0 - Added summarising of sequences post-run. Improved documents. Moved to tools/.
    # 0.4.0 - Added VecScreen mode.
    # 0.5.0 - Added Telomere finding mode based on https://github.com/JanaSperschneider/FindTelomeres.
    # 0.6.0 - Added GenomeSize option.
    # 0.7.0 - Added PurgeHaplotigs, minimum vescreen hit length, eFDR and vecscreen coverage.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add download of NCBI Vector Database if vecdb=ncbi.
    # [ ] : Need to add a eukaryote mode and/or minlen for VecScreen - too many expected hits with NCBI rules.
    # [ ] : Add eFDR calculation and filtering to VecScreen.
    # [ ] : Implement long read BUSCO read depth and genome size prediction.
    # [ ] : Implement running of Purge haplotigs using BUSCO read depths to set cutoffs.
    # [Y] : Fix generation of warnings and error when Description is not available!
    # [ ] : Change purgehap to run purge_haplotigs only
    # [ ] : Add saving and reloading of total read count and depth?
    # [ ] : Add renaming of sequences to the dipnr analysis mode.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Diploidocus', '0.7.0', 'February 2020', '2017')
    description = 'Diploid genome assembly analysis tools.'
    author = 'Dr Richard J. Edwards.'
    comments = ['NOTE: telomere finding rules are based on https://github.com/JanaSperschneider/FindTelomeres',
                'This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show Minimap2 run (rje_paf) commandline options?',default='N'): out.verbose(-1,4,text=rje_paf.__doc__)
            if rje.yesNo('Show SeqList commandline options?',default='N'): out.verbose(-1,4,text=rje_seqlist.__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: print(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: print('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: print('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
paf_defaults = {'N':'250','p':'0.0001','x':'asm20'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Diploidocus Class                                                                                       #
#########################################################################################################################
class Diploidocus(rje_obj.RJE_Object):
    '''
    Diploidocus Class. Author: Rich Edwards (2019).

    Str:str
    - BAM=FILE        : BAM file of reads mapped onto assembly [$BASEFILE.bam]
    - BUSCO=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    - Parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    - Parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    - PurgeMode=X     : Rules used for purgehap analysis (simple/complex/nala) [complex]
    - RunMode=X       : Diploidocus run mode [insilico/sortnr/diphap/vecscreen]
    - ScreenDB=FILE   : File of vectors/contaminants to screen out using blastn and VecScreen rules []
    - ScreenMode=X    : Action to take following vecscreen searching (report/mask/trim/purge) [report]
    - SeqIn=FILE      : Input sequence assembly (sortnr/diphap modes) []
    - SeqOut=FILE     : Output sequence assembly [$BASEFILE.fasta]
    - TeloFwd=X      : Basic telomere sequence for search [C{2,4}T{1,2}A{1,3}]
    - TeloRev=X      : Basic telomere sequence for search [TTAGGG]
    - TmpDir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]

    Bool:boolean
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - IncludeGaps=T/F : Whether to include gaps in the zero coverage bases for adjustment (see docs) [False]
    - QuickDepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    - Summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    - ZeroAdjust=T/F  : Add zero coverage bases to purge_haplotigs LowPerc and adjust total [True]
    - 10xTrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer reads data [False]

    Int:integer
    - GenomeSize=INT  : Haploid genome size (bp) [0]
    - LenFilter=X     : Min read length for filtered subreads [500]
    - MinGap=INT      : Minimum length of a stretch of N bases to count as a gap for exclusion [10]
    - MinMedian=INT   : Minimum median depth coverage to avoid low coverage filter [3]
    - MinVecHit=INT   : Minimum length for a screendb match [50]
    - PHLow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
    - PHMid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
    - PHHigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
    - ReadBP=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
    - SCDepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
    - TeloSize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]

    Num:float
    - CheckCov=PERC   : Percentage coverage for double-checking partial exact matches [95]
    - eFDR=NUM        : Expected FDR threshold for VecScreen queries (0 is no filter) [1.0]
    - RQFilter=X      : Minimum RQ for output subreads [0]
    - TeloPerc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]

    File:file handles with matching str filenames
    
    List:list
    - KmerReads=FILELIST   : File of reads for KAT kmer analysis []
    - Reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType=LIST   : List of ont/pb file types matching reads for minimap2 mapping [ont]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database Object
    - Forker = Forker Object
    - SeqIn = SeqList Object for input sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BAM','BUSCO','GenomeSize','Parent1','Parent2','PurgeMode','RunMode','ScreenDB','ScreenMode','SeqIn','SeqOut','DebugStr','TeloFwd','TeloRev','TmpDir']
        self.boollist = ['DocHTML','IncludeGaps','QuickDepth','Summarise','ZeroAdjust','10xTrim']
        self.intlist = ['GenomeSize','LenFilter','MinGap','MinMedian','MinVecHit','PHLow','PHMid','PHHigh','SCDepth','ReadBP','TeloSize']
        self.numlist = ['CheckCov','eFDR','RQFilter','TeloPerc']
        self.filelist = []
        self.listlist = ['KmerReads','Reads','ReadType']
        self.dictlist = []
        self.objlist = ['Forker','SeqIn']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'PurgeMode':'complex','RunMode':'diploidocus','ScreenMode':'report','TeloFwd':'C{2,4}T{1,2}A{1,3}','TeloRev':'','TmpDir':'./tmpdir/'})
        self.setBool({'DocHTML':False,'IncludeGaps':False,'QuickDepth':False,'Summarise':True,'ZeroAdjust':True,'10xTrim':False})
        self.setInt({'LenFilter':500,'MinMedian':3,'MinVecHit':50,'GenomeSize':0,'ReadBP':0,'TeloSize':50,'MinGap':10})
        self.setNum({'CheckCov':95.0,'eFDR':1.0,'RQFilter':0,'TeloPerc':50.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.list['ReadType'] = ['ont']
        self.obj['SeqIn'] = None
        self.obj['Forker']  = rje_forker.Forker(self.log,['logfork=F']+self.cmd_list)
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['GenomeSize','DebugStr','PurgeMode','RunMode','ScreenMode','TeloFwd','TeloRev'])   # Normal strings
                self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['BAM','Parent1','Parent2','ScreenDB','SeqIn','SeqOut','BUSCO'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DocHTML','IncludeGaps','QuickDepth','Summarise','ZeroAdjust','10xTrim'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['LenFilter','MinGap','MinMedian','MinVecHit','PHLow','PHMid','PHHigh','ReadBP','SCDepth','TeloSize'])   # Integers
                self._cmdReadList(cmd,'float',['eFDR','RQFilter']) # Floats
                self._cmdReadList(cmd,'perc',['CheckCov','TeloPerc']) # Percentage
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ReadType'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['KmerReads','Reads']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStrLC('GenomeSize'):
         try: self.setInt({'GenomeSize':rje_seqlist.bpFromStr(self.getStrLC('GenomeSize'))})
         except:
            self.errorLog('Problem with GenomeSize. Setting genomesize=0')
            self.setInt({'GenomeSize':0})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main Diploidocus run method
        '''
        # Diploidocus: Diploid genome assembly analysis tools

        Diploidocus is a sequence analysis toolkit for a number of different analyses related to diploid genome assembly.
        The main suite of analyses combines long read depth profiles, short read kmer analysis, assembly kmer analysis,
        BUSCO gene prediction and contaminant screening for a number of assembly tasks including genome size prediction,
        contamination identification, haplotig identification/removal and low quality contig/scaffold filtering. In addition,
        Diploidocus has functions for removing redundancy, generating a non-redundant pseudo-diploid assembly with primary
        and secondary scaffolds from 10x pseudohap output, and creating an in silico diploid set of long reads from two
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

        Diploidocus is available as part of SLiMSuite, or via a standalone GitHub repo at
        <https://github.com/slimsuite/diploidocus>.

        ---

        # Running Diploidocus

        Diploidocus is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/diploidocus.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [Diploidocus git repo](https://github.com/slimsuite/diploidocus), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [Diploidocus git repo](https://github.com/slimsuite/diploidocus)
        for running on example data.

        For `sortnr` and `diphapnr` mode, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
        environment `$PATH` or given to Diploidocus with the `minimap2=PROG` setting.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Main Diploidocus run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        runmode=X       : Diploidocus run mode [sortnr/diphap/diphapnr/purgehap/telomere/vecscreen/insilico/gensize]
        basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
        summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
        genomesize=INT  : Haploid genome size (bp) [0]
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
        ### ~ Genome size prediction & Purge haplotigs options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
        bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
        reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype=LIST   : List of ont/pb file types matching reads for minimap2 mapping [ont]
        readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
        quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
        kmerreads=FILELIST : File of high quality reads for KAT kmer analysis []
        10xtrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer Reads data for KAT analysis [False]
        scdepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
        minmedian=INT   : Minimum median depth coverage to avoid low coverage filter [3]
        phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
        phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
        phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
        zeroadjust=T/F  : Add zero coverage bases to purge_haplotigs LowPerc and adjust total [True]
        includegaps=T/F : Whether to include gaps in the zero coverage bases for adjustment (see docs) [False]
        mingap=INT      : Minimum length of a stretch of N bases to count as a gap for exclusion [10]
        purgemode=X     : Rules used for purgehap analysis (simple/complex/nala) [complex]
        ### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        telofwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
        telorev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]
        telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
        teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
        ### ~ VecScreen options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        screendb=FILE   : File of vectors/contaminants to screen out using blastn and VecScreen rules []
        screenmode=X    : Action to take following vecscreen searching (report/mask/trim/purge) [report]
        minvechit=INT   : Minimum length for a screendb match [50]
        efdr=NUM        : Expected FDR threshold for VecScreen queries (0 is no filter) [1.0]
        ### ~ SortNR filtering/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        checkcov=PERC   : Percentage coverage for double-checking partial exact matches [95]
        seqout=FILE     : Output sequence assembly [$BASEFILE.nr.fasta]
        ### ~ In silico diploid input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        rqfilter=X      : Minimum RQ for output subreads [0]
        lenfilter=X     : Min read length for filtered subreads [500]
        parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
        parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
        See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        # Diploidocus run modes

        Details for the main Diploidocus run modes are given below.

        **NOTE:** Diploidocus is under development and documentation might be a bit sparse. Please contact the author or
        post an issue on GitHub if you have any questions.

        ---

        ## Main Diploidocus filtering [runmode=diploidocus]
        Diploidocus builds on the PurgeHaplotigs classifications to use the depth bins, KAT assembly kmer frequencies and
        BUSCO results to reclassify scaffolds and partition them into:

        * `*.diploidocus.fasta` = the scaffolds kept for the next round of PurgeHap
        * `*.core.fasta` = the same set of scaffolds, minus repeats
        * `*.junk.fasta` = low coverage scaffolds, removed as junk
        * `*.quarantine.fasta` = putative haplotigs, removed from the assembly but retained for reference.

        **NOTE:** PurgeHaplotigs was not using zero coverage bases in its percentages. This is now fixed by Diploidocus.

        _Further details coming soon!_

        ---

        ## Cycled Diploidocus filtering [runmode=dipcycle]
        Diploidocus can be automated to cycle through repeated rounds of the main purge_haplotigs/filtering method until
        no further scaffolds are removed. Each cycle will run Diploidocus with a numerical suffix, e.g. `$BASEFILE.1.*`,
        using the `*.diploidocus.fasta` output from the previous cycle as input. Cycling will continue until no further
        scaffolds are filtered into either `*.quarantine.fasta` or `*.junk.fasta`.

        Final outputs from the final cycle will then be copied under the original `$BASEFILE` prefix:

        * `$BASEFILE.diploidocus.tdt` = Final ratings for the final set of scaffolds. (See earlier cycles for purged.)
        * `$BASEFILE.diploidocus.fasta` = the scaffolds kept from the final Diploidocus cycle
        * `$BASEFILE.core.fasta` = the same set of scaffolds, minus repeats
        * `$BASEFILE.quarantine.fasta` = concatenated purged scaffolds from all Diploidocus cycles.
        * `$BASEFILE.junk.fasta` = concatenated low coverage and low quality scaffolds, removed as junk, from all cycles.

        ---

        ## Genome size prediction [runmode=gensize]
        _Details coming soon!_

        ---

        ## Running Purge_haplotigs using BUSCO-guided cutoffs [runmode=purgehap]
        _Coming soon! (Currently will run full diploidocus mode)_

        ---

        ## Telomere finding [runmode=telomere]
        _Details coming soon!_

        ---

        ## Vector/contamination screening [runmode=vecscreen]
        _Details coming soon!_

        ---

        ## Sorted non-redundant assembly cleanup [runmode=sortnr]

        The sorted non-redundant assembly cleanup mode (`runmode=sortnr`) screens out any sequences that are 100% gap,
        then removes any sequences that are 100% redundant with other sequences in the input. This includes full and
        partial matches, i.e. if sequence X is wholly contained within sequence Y then X will be removed.

        First, sequences are loaded from the file given with `seqin=FILE` and any [rje_seqlist](http://rest.slimsuite.unsw.edu.au/seqlist)
        filters and sequence sorting are applied to the input. Sequences that are 100% Ns are removed and any gaps
        exceeding 10 nt are reduced to 10 `N`s (`NNNNNNNNNN`) to prevent minimap2 from splitting sequences on long gaps.
        These gap-reduced sequences are output to `$BASEFILE.tmp.fasta` and used for an all-by-all minimap2 search.

        By default, minimap2 is run with the options to generate a `$BASEFILE.tmp.paf` file:

            --cs -p 0.0001 -t 4 -x asm20 -N 250

        To modify minimap2 search settings, please see the [rje_paf](http://rest.slimsuite.unsw.edu.au/rje_paf)
        documentation.

        **NOTE:** These run options can probably be made more stringent to speed up minimap2 without loss of function.
        Future releases may alter defaults accordingly.

        Minimap2 output is parsed to identify scaffold-scaffold matches. Self-hits are ignored.
        The minimum (gap-reduced) sequence length is used as a rapid parsing filter: any minimap2 matches that are less
        than 95% of the query sequence (`Length`+`nn` fields) or less that 100% identity (`Identity`+`nn`)/(`Length`+`nn`)
        are filtered during parsing.

        **NOTE:** Future releases may feature an option to reduce the global percentage identity cut-off. Please contact
        the author if you wish to see this implemented.

        Minimap2 hits are then processed reverse-sorted by reference sequence size (e.g. scaffold length). Any hits
        where either sequence has already been filtered are skipped. Otherwise, if the match (as determined by the
        length of `:` regions in the CS string) matches the query length, the Query sequence will be flagged for
        remove as "identical to" or "contained within" the Hit. (Mutually partial overlapping exact matches are NOT
        filtered.) Filtered IDs and their matches are output to `$BASEFILE.redundant.txt`.

        Once all sequences have been filtered, the remaining sequences are output to: `$BASEFILE.nr.fasta`.

        **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
        add the `sortseq=invsize` command to the Diploidocus run command.

        Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to:
        `$BASEFILE.summarise.tdt`.

        Temporary gap-reduced and minimap2 PAF files are deleted unless running in `debug` or
        `dev` modes.

        ---

        ## Pseudodiploid to primary and alternative haploptigs [runmode=diphap(nr)]

        This protocol is based on 10x assemblies made for multiple organisms with supernova v2.0.0 and supernova v2.1.1.
        In each case, some redundancy in output was discovered (a) within pseudohap output, and (b) in terms of fully
        homozygous (identical) scaffolds shared by both haplotigs. It was also not entirely clear on what basis a
        particular haplotype was assigned to pseudohap1 or pseudohap2.

        The general workflow therefore sought to remove redundancy, generate a set of primary scaffolds based on scaffold
        length, and generate a non-redundant set of alternative scaffolds where heterozygosity exists. If `diphapnr` mode
        is used, the full workflow is implement by first running the `sortnr` workflow described above. In the reduced
        `diphap` mode, redundancy is not removed first.

        Sequences are loaded and matching haplotigs identified based on their names. Sequence names MUST end `HAP(\d+)`,
        where `(\d+)` indicates an integer that matches up haplotigs (as produced by supernova pseudohap2 output, for
        example). This is **not** a pipeline to identify haplotig pairs, it is purely for splitting identified
        haplotigs into primary and alternative assemblies.

        Processing itself is quite simple. Haplotig pairs are identified based on matching `HAP(\d+)` numbers. Where a
        single haplotig is found, it is assigned as `diploid`, under the assumption that the two haplotigs were identical
        and one was removed. (It is possible that only one parent had this scaffold, e.g. sex chromosomes, so some post-
        processing of descriptions may be required.) If two haplotigs with the same number are identified, the longest
        is assigned to `haploidA` and the shorter `haploidB`.

        The **Primary Assemmbly** is then compiled from all `haploidA` and `diploid` sequences. These are given `pri`
        prefixes and output to `$BASEFILE.pri.fasta`. The **Alternative** comprised of all `haploidB` sequences is output
        to `$BASEFILE.alt.fasta`. If redundancy has been removed, this will likely be a subset of the full assembly. The
        combined set of all primary and alternative sequences is output to `$BASEFILE.dipnr.fasta`.

        **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
        add the `sortseq=invsize` command to the Diploidocus run command.

        Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to
        `$BASEFILE.summarise.tdt`:

        * `$BASEFILE.dipnr.fasta` = Combined pseudodiploid with `haploidA`, `haploidB` and `diploid` annotation.
        * `$BASEFILE.pri.fasta` = Primary assembly with `haploidA` and `diploid` sequences.
        * `$BASEFILE.alt.fasta` = Alternative assembly with `haploidB` sequences.


        ---

        ## In silico diploid generator [runmode=insilico]

        This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
        parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
        parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
        identifier table.)

        A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
        unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
        selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
        two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
        This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
        no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
        subreads will be added from the other parent if they reduce the difference in cumulative output for each parent, or
        until `lenfilter=X` is reached.

        Final output will be a `*.LXXXRQXX.fasta` file in which each parent has a similar total sequence content and for
        which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
        assemblies, where one parent has higher quality data than the other.

        NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
        higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
        minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
        relaxed. Similarly, only sequences above `lenfilter=X` in length will be output. These are the figures given in the
        `LXXXRQXX` part of the output file, e.g. defaults of RQ>=0.84 and Len>=500 generates `*.L500RQ84.fas`.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.getBool('DocHTML'): return self.docHTML()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#MODE',self.getStrLC('RunMode'))
            if self.getStrLC('RunMode') in ['sortnr','seqnr','nrseq','nr']: return self.sortNR()
            elif self.getStrLC('RunMode') in ['diphap','diphapnr']: return self.dipHap()
            elif self.getStrLC('RunMode') == 'insilico':
                return self.inSilicoHybrid()
            elif self.getStrLC('RunMode') == 'vecscreen': return self.vecScreen()
            elif self.getStrLC('RunMode').startswith('telomere'): return self.findTelomeres()
            elif self.getStrLC('RunMode') == 'diploidocus': return self.diploidocus()
            elif self.getStrLC('RunMode').startswith('purgehap'): return self.purgeHaplotigs()
            elif self.getStrLC('RunMode') in ['gensize','genomesize']: return self.genomeSize()
            elif self.getStrLC('RunMode') in ['dipcycle','purgecycle']: return self.purgeCycle()
            else: raise ValueError('RunMode="%s" not recognised!' % self.getStrLC('RunMode'))
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def purgeCycle(self):  ### Repeat Diploidocus purge cycles to convergence (none removed).                    # v0.7.0
        '''
        ### ~ Cycled Diploidocus filtering [runmode=dipcycle] ~ ###
        Diploidocus can be automated to cycle through repeated rounds of the main purge_haplotigs/filtering method until
        no further scaffolds are removed. Each cycle will run Diploidocus with a numerical suffix, e.g. `$BASEFILE.1.*`,
        using the `*.diploidocus.fasta` output from the previous cycle as input. Cycling will continue until no further
        scaffolds are filtered into either `*.quarantine.fasta` or `*.junk.fasta`.

        Final outputs from the final cycle will then be copied under the original `$BASEFILE` prefix:

        * `$BASEFILE.diploidocus.tdt` = Final ratings for the final set of scaffolds. (See earlier cycles for purged.)
        * `$BASEFILE.diploidocus.fasta` = the scaffolds kept from the final Diploidocus cycle
        * `$BASEFILE.core.fasta` = the same set of scaffolds, minus repeats
        * `$BASEFILE.quarantine.fasta` = concatenated purged scaffolds from all Diploidocus cycles.
        * `$BASEFILE.junk.fasta` = concatenated low coverage and low quality scaffolds, removed as junk, from all cycles.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newbase = basefile = self.baseFile(strip_path=True)
            cycle = 0
            prevseqx = 0
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])

            ### ~ [2] ~ Cycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while seqlist.seqNum() != prevseqx:
                # Increment cycle and make new basefile
                prevseqx = seqlist.seqNum()
                cycle += 1
                oldbase = newbase
                newbase = '{}.{}'.format(basefile,cycle)
                seqout = '{}.diploidocus.fasta'.format(newbase)
                self.printLog('#SEQIN','Using {} as input for Diploidocus cycle {}.'.format(seqin,cycle))
                self.printLog('#CYCLE','Cycle {}: running as {}.*'.format(cycle,newbase))
                # Run Diploidocus
                info = makeInfo()
                cyccmd = self.cmd_list+['basefile={}'.format(newbase),'i=-1','runmode=diploidocus','seqin=%s' % seqin]
                cmd_list = rje.getCmdList(cyccmd,info=info)   # Reads arguments and load defaults from program.ini
                out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
                out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
                out.printIntro(info)                                # Prints intro text using details from Info object
                log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
                if self.debugging():
                    Diploidocus(log,['dna=T']+cmd_list+['i=1']).run()
                else:
                    Diploidocus(log,['dna=T']+cmd_list+['i=-1']).run()
                # Check and process output
                if not rje.exists(seqout):
                    raise IOError('Expected %s output for Cycle %s not found! Check %s.log. Aborting run.' % (seqout,cycle,newbase))
                self.printLog('#CYCLE','Cycle {} complete. See {}.log for details.'.format(cycle,newbase))
                seqin = seqout
                seqlist = self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file','seqin=%s' % seqin])

            ### ~ [3] Tidy up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#CYCLE','Convergence achieved after cycle {}'.format(cycle))
            for ext in ['diploidocus.tdt','diploidocus.fasta','core.fasta']:
                if rje.exists('{}.{}'.format(newbase,ext)):
                    rje.backup(self,'{}.{}'.format(basefile,ext))
                    shutil.copy('{}.{}'.format(newbase,ext),'{}.{}'.format(basefile,ext))
                    self.printLog('#COPY','Output copied: {}.{} -> {}.{}'.format(newbase,ext,basefile,ext))
            for ext in ['quarantine.fasta','junk.fasta']:
                cycfiles = []
                for i in range(1,cycle+1):
                    cycbase = '{}.{}'.format(basefile,i)
                    if rje.exists('{}.{}'.format(cycbase,ext)): cycfiles.append('{}.{}'.format(cycbase,ext))
                if cycfiles:
                    rje.backup(self,'{}.{}'.format(basefile,ext))
                    os.system('cat {} > {}.{}'.format(' '.join(cycfiles),basefile,ext))
                    self.printLog('#COPY','Output copied: {} *.{} files -> {}.{}'.format(len(cycfiles),ext,basefile,ext))

            ### ~ [4] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist.getBool('Summarise'):
                for seqset in ['diploidocus','core','quarantine','junk']:
                    setfas = '%s.%s.fasta' % (basefile,seqset)
                    if not rje.baseFile(setfas): continue
                    seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % setfas,'autofilter=F']
                    rje_seqlist.SeqList(self.log,seqcmd)

        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if not self.baseFile(return_none=''):
                if self.getStrLC('SeqIn'): self.baseFile(rje.baseFile(self.getStr('SeqIn'),strip_path=True))
                else: self.baseFile('diploidocus')
            self.printLog('#BASE','Output file basename: %s' % self.baseFile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def seqinObj(self): ### Returns the a SeqList object for the SeqIn file
        '''
        Returns the a SeqList object for the SeqIn file.
        :return: self.obj['SeqIn']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['SeqIn']:
                self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
        except:
            self.errorLog('Diploidocus.seqinObj() error')
        return self.obj['SeqIn']
#########################################################################################################################
    def docHTML(self):  ### Generate the PAFScaff Rmd and HTML documents.                                        # v0.1.0
        '''Generate the PAFScaff Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            info = self.log.obj['Info']
            prog = '%s V%s' % (info.program,info.version)
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
            #!# Replace this with documentation text?
            rtxt += string.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2019 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown Diploidocus documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def loggedSysCall(self,cmd,syslog=None,stderr=True,append=True,verbosity=1,nologline=None):    ### Makes a system call, catching output in log file
        '''
        Makes a system call, catching output in log file.
        :param cmd:str = System call command to catch
        :param syslog:str [None] = Filename for system log in which to capture output. Will use $BASEFILE.sys.log if None
        :param stderr:bool [True] = Whether to also capture the STDERR
        :param append:bool [False] = Whether to append the log if it exists
        :param verbosity:int [1] = Verbosity level at which output also goes to screen (tee, not redirect)
        :param nologline:str [None] = Default logline returned if nothing is in the log
        :return: last line of syslog output
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# System log filename
            if not syslog: syslog = '{}.sys.log'.format(self.baseFile())
            #i# Setup log with command and identify position to capture extra content
            headline = '[{}] {}\n'.format(rje.dateTime(),cmd)
            if append: open(syslog,'a').write(headline)
            else:
                rje.backup(self,syslog,appendable=False)
                open(syslog,'w').write(headline)
            #i# Identify position at end of file
            fend = rje.endPos(filename=syslog)
            #i# Generate system command
            if stderr:
                if ' > ' in cmd: self.warnLog('Cannot capture stderr or stdout for command: {}'.format(cmd))
                else: cmd = '{} 2>&1'.format(cmd)
            if append:
                if self.v() >= verbosity: cmd = '{} | tee -a {}'.format(cmd,syslog)
                else: cmd = '{} >> {}'.format(cmd,syslog)
            else:
                if self.v() >= verbosity: cmd = '{} | tee {}'.format(cmd,syslog)
                else: cmd = '{} > {}'.format(cmd,syslog)
            ### ~ [2] ~ Process System Call ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SYS',cmd)
            os.system(cmd)
            logline = nologline
            if not logline:
                logline = 'WARNING: No run log output!'
            SYSLOG = open(syslog,'r')
            SYSLOG.seek(fend)
            for readline in SYSLOG.readlines():
                readline = rje.chomp(readline)
                if readline: logline = readline
            self.printLog('#SYSEND',logline)
            return logline
        except:
            self.errorLog('Diploidocus.loggedSysCall() error')
            return None
#########################################################################################################################
    ### <3> ### In silico hybrid run method                                                                             #
#########################################################################################################################
    def inSilicoHybrid(self):  ### Filter and combine subreads from parent and output to fasta file.
        '''
        Filter and combine subreads from parent and output to fasta file.

        This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
        parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
        parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
        identifier table.)

        A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
        unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
        selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
        two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
        This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
        no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
        subreads will be added from the other parent if they reduce the difference in cumulative output for each parent.

        Final output will be a `*.subreads.fasta` file in which each parent has a similar total sequence content and for
        which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
        assemblies, where one parent has higher quality data than the other.

        NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
        higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
        minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
        relaxed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Parent 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 1 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent1: %s' % self.getStr('Parent1'))
            base1 = rje.baseFile(self.getStr('Parent1'))
            parent1 = smrtscape.SMRTSCAPE(self.log,self.cmd_list+['batch=%s' % self.getStr('Parent1'),'basefile=%s' % base1])
            parent1.setup()
            udb1 = parent1.udb()
            cdb = parent1.db('smrt',add=True,mainkeys=['Name'])
            cdb.dataFormat({'SMRT':'int'})
            cx = cdb.entryNum()
            ## ~ [0a] Parent 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 2 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent2: %s' % self.getStr('Parent2'))
            base2 = rje.baseFile(self.getStr('Parent2'))
            parent2 = smrtscape.SMRTSCAPE(self.log,self.cmd_list+['batch=%s' % self.getStr('Parent2'),'basefile=%s' % base2])
            parent2.setup()
            udb2 = parent2.udb()
            cdb2 = parent2.db('smrt',add=True,mainkeys=['Name'])
            cdb2.dataFormat({'SMRT':'int'})
            # Shift all of the Parent2 SMRT IDs to avoid conflict with Parent1
            for entry in cdb2.entries() + udb2.entries(): entry['SMRT'] = entry['SMRT'] + cx
            cdb = parent1.db().mergeTables(cdb,cdb2)
            ## ~ [0c] Output Sequence File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ DIPLOIDOCUS SUBREADS ~~~~~~~~~~~~~~~~~~~~ #')
            minlen = self.getInt('LenFilter')
            minrq = self.getNum('RQFilter')
            rqstr = '%s' % minrq
            filtfile = '%s.L%sRQ%s.fasta' % (self.baseFile(),minlen,rqstr[2:])
            ## ~ [0d] Input Sequence Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqbatch = []   # List of SeqList objects
            self.printLog('#BATCH','%s sequence files to process.' % rje.iLen(parent1.list['Batch']+parent2.list['Batch']))
            for seqfile in parent1.list['Batch']+parent2.list['Batch']:
                seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=F','seqin=%s' % seqfile,'autofilter=F']
                seqbatch.append(rje_seqlist.SeqList(self.log,seqcmd))
            self.printLog('#BATCH','%s sequence files to summarise.' % rje.iLen(seqbatch))
            if not seqbatch: raise IOError('No batch input fasta files found! Make sure parentN=FILE settings given *.fofn.')
            ## ~ [0e] Setup subread lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elists = [udb1.sortedEntries('Len',reverse=True),udb2.sortedEntries('Len',reverse=True)]
            plen = [0,0]    # Summed lengths for each parent
            pseq = [0,0]    # Total sequence number for each parent
            prq = [0,0]     # Total sequence RQ for each parent (convert to mean)
            if not elists[0] or not elists[1]: raise ValueError('No Unique ZMW subreads for one or both parents!')
            lastlen = max(elists[0][0]['Len'],elists[1][0]['Len'])    # Length of last selected read
            for elist in elists:
                while elist and elist[0]['RQ'] < minrq: elist.pop(0)
            if not elists[0] or not elists[1]: raise ValueError('No Unique ZMW subreads for one or both parents!')
            nextp = 0       # Index of next parent to use
            if elists[0][0]['Len'] < elists[1][0]['Len']: nextp = 1

            ### ~ [1] Filter and Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Filter Unique Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            zmwlist = []    # List of (smrt,zmw) meeting filtering criteria
            ux = 0.0; utot = len(elists[0])+len(elists[1])
            while lastlen:
                self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                elist = elists[nextp]
                while elist and elist[0]['RQ'] < minrq: elist.pop(0); ux += 100.0
                if elist and elist[0]['Len'] < minlen: ux += 100.0 * len(elist); elist = []
                if not elist: nextp = 1 - nextp; break  # Finish
                entry = elist.pop(0); ux += 100.0
                zmwlist.append((entry['SMRT'],entry['ZMW'],entry['Pos']))
                plen[nextp] += entry['Len']
                prq[nextp] += entry['RQ']
                pseq[nextp] += 1
                if plen[1-nextp] <= plen[nextp]: nextp = 1 - nextp
                lastlen = entry['Len']
            ## ~ [1b] Final processing of last reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while elists[nextp]:
                elist = elists[nextp]
                while elist and elist[0]['RQ'] < minrq:
                    self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                    elist.pop(0); ux += 100.0
                while elist and elist[0]['Len'] >= minlen:
                    self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                    entry = elist.pop(0); ux += 100.0
                    pdiff = rje.modulus(plen[0]-plen[1])
                    ediff = rje.modulus(plen[nextp]+entry['Len']-plen[1-nextp])
                    if ediff >= pdiff: elists[nextp] = []; break    #Finish!
                    zmwlist.append((entry['SMRT'],entry['ZMW'],entry['Pos']))
                    plen[nextp] += entry['Len']
                    prq[nextp] += entry['RQ']
                    pseq[nextp] += 1
            self.printLog('\r#DIP','Diploidising subreads complete: %s subreads to output.' % rje.iLen(zmwlist))
            self.printLog('\r#DIP','%s: %s seq; %s bp (%.1fX); %.3f mean RQ.' % (self.getStr('Parent1'),rje.iStr(pseq[0]),rje.iStr(plen[0]),1.0*plen[0]/self.getInt('GenomeSize'),prq[0]/pseq[0]))
            self.printLog('\r#DIP','%s: %s seq; %s bp (%.1fX); %.3f mean RQ.' % (self.getStr('Parent2'),rje.iStr(pseq[1]),rje.iStr(plen[1]),1.0*plen[1]/self.getInt('GenomeSize'),prq[1]/pseq[1]))
            ## ~ [1b] Extract Filtered Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,filtfile)
            SEQOUT = open(filtfile,'w')
            sx = 0.0; stot = 0; sn = len(seqbatch); fx = 0
            for seqlist in seqbatch:
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/9/0_3967 RQ=0.784
                si = 100.0/seqlist.seqNum(); stot += seqlist.seqNum()
                for seq in seqlist.seqs():
                    self.progLog('\r#OUT','Extracting subreads: %.2f%%' % (sx/sn)); sx += si
                    (name,sequence) = seqlist.getSeq(seq)
                    try: [smrt,zmw,pos,rq] = string.split(string.replace(name,'/',' '))
                    except:
                        [smrt,zmw,pos] = string.split(string.replace(name,'/',' '))
                        rq = minrq
                    if (cdb.data(smrt)['SMRT'],int(zmw),pos) not in zmwlist: continue
                    SEQOUT.write('>%s\n%s\n' % (name,sequence)); fx += 1
            self.printLog('\r#OUT','Saved %s filtered subreads to %s.' % (rje.iStr(fx),filtfile))

            ### ~ [2] Summarise Filtered File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % filtfile,'autofilter=F']
            rje_seqlist.SeqList(self.log,seqcmd)

            return True
        except: self.errorLog('%s.run error' % self.prog()); return False
#########################################################################################################################
    ### <4> ### SortNR & DipHap run methods                                                                             #
#########################################################################################################################
    def sortNR(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            debugstr = 'NOTDEBUGGING'
            if self.getStrLC('DebugStr'): debugstr = self.getStr('DebugStr')
            ## ~ [1a] Input Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F'])
            tmpseq = '%s.tmp.fasta' % self.baseFile()
            TMPSEQ = open(tmpseq,'w')
            #i# seqdict will store all the sequences, mapped by shortname onto PAF results
            #i# As sequences are filtered, they will be removed from seqdict and added to badseq. Remaining sequences will be saved.
            seqdict = seqin.seqNameDic()
            seqin.dict['Filter']['NullSeq'] = 0
            seqin.nextSeq()
            minlen = seqin.seqNonX()
            sx = 0.0; stot = seqin.seqNum()
            while seqin.currSeq():
                self.progLog('\r#LEN','Scanning sequences: %.1f%%' % (sx/stot)); sx += 100
                if seqin.seqNonX():
                    (name,sequence) = seqin.currSeq()
                    sequence = string.join( re.split('[Nn]{10}[Nn]+',sequence), 'NNNNNNNNNN')
                    TMPSEQ.write('>%s\n%s\n' % (name,sequence))
                    minlen = min(minlen,len(sequence))
                else:
                    seqin.dict['Filter']['NullSeq'] += 1
                    seqdict.pop(seqin.shortName())
                if not seqin.nextSeq(): break
            TMPSEQ.close()
            self.printLog('\r#NULL','%s 100%% N sequences filtered.' % (rje.iStr(seqin.dict['Filter']['NullSeq'])))
            self.printLog('#TMP','%s gap-reduced sequences output to temporary file %s for minimap2 search' % (rje.iLen(seqdict),tmpseq))
            self.debug(minlen)

            #!# There is an issue of long gaps causing problems, so need to replace all gaps with 10 Ns prior to PAF
            #!# generation, and make sure that files are subsequently cleaned up (unless debug=T or dev=T).

            ## ~ [1b] MiniMap PAF generation and parsing object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pafin = '%s.tmp.paf' % self.baseFile()
            self.obj['PAF'] = paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'basefile=%s' % self.baseFile(),'seqin=%s' % tmpseq,'reference=%s' % tmpseq,'sortseq=None'])
            paf.obj['DB'] = self.obj['DB']
            paf.setNum({'MinLocID':100})
            paf.setInt({'MinLocLen':int(minlen*self.getPerc('CheckCov'))})
            paf.dict['MapOpt'] = rje.combineDict(paf_defaults,paf.dict['MapOpt'],overwrite=True)
            paf.setup()
            if not rje.exists(pafin) or self.force():
                rje.backup(self,pafin)
                paf.minimap2()
            # Parse PAF with initial minloclen and minlocid filtering
            pafdb = paf.parsePAF(debugstr=debugstr)
            pafdb.newKey(['Hit','Qry','#'])

            ### ~ [2] NR Sequence filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nrfile = '%s.redundant.txt' % self.baseFile()
            rje.backup(self,nrfile)
            NRFILE = open(nrfile,'w')
            #i# pafhead = Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality .. cs
            seqin.dict['Filter']['NR'] = 0
            px = 0.0; ptot = pafdb.entryNum()
            #i# Work through PAF table starting with longest hits and investigate status of queries
            for pentry in pafdb.sortedEntries('SbjLen',reverse=True):
                self.progLog('\r#NR','Filtering redundant sequences from PAF CS strings: %.2f%%' % (px/ptot)); px += 100
                qry = pentry['Qry']
                hit = pentry['Hit']
                if qry.endswith(debugstr): self.bugPrint(pentry)
                #i# Ignore self-hits
                if qry == hit: continue
                #i# Skip any hits already filtered
                if qry not in seqdict: continue
                if hit not in seqdict: continue
                cstats = paf.statsFromCS(pentry['cs'])
                if qry.endswith(debugstr): self.bugPrint(cstats)
                #!# Might need to use an actual value rather than a percentage!
                if self.getPerc('CheckCov') * pentry['QryLen'] <= cstats[':'] < pentry['QryLen']:  # Hit might be redundant with query
                    qryseq = seqin.getSeq(seqdict[qry])[1]
                    hitseq = seqin.getSeq(seqdict[hit])[1]
                    if pentry['Strand'] == '-': qryseq = rje_sequence.reverseComplement(qryseq)
                    if qryseq in hitseq: cstats[':'] = pentry['QryLen']
                    if qry.endswith(debugstr): self.bugPrint('\n%s\n vs \n%s\n = %s' % (qryseq,hitseq,qryseq in hitseq))
                if cstats[':'] == pentry['QryLen']:  # Hit is redundant with query
                    rtxt = qry
                    if pentry['Strand'] == '-': rtxt = rtxt + ' (RevComp)'
                    if cstats[':'] == pentry['SbjLen']: rtxt += ' = identical to '
                    else: rtxt += ' = contained within '
                    rtxt += hit
                    NRFILE.write('%s\n' % rtxt)
                    seqin.dict['Filter']['NR'] += 1
                    seqdict.pop(qry)
                if qry.endswith(debugstr): self.debug('%s in seqdict?: %s' % (qry,qry in seqdict))
            NRFILE.close()
            self.printLog('\r#NR','%s redundant sequences filtered from PAF CS strings: see %s.' % (rje.iStr(seqin.dict['Filter']['NR']),nrfile))

            ### ~ [3] Save filtered sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodseq = seqdict.values()
            goodseq.sort()
            seqin.list['Seq'] = goodseq
            if not seqin.getStrLC('SeqOut'):
                seqin.setStr({'SeqOut':'%s.nr.fasta' % self.baseFile()})
            seqfiles = [seqin.getStr('SeqIn'),seqin.getStr('SeqOut')]
            seqin.saveSeq()
            seqin.setStr({'SeqIn':seqin.getStr('SeqOut')})
            ## ~ [3a] Option to delete PAF files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for tmp in ['paf','paf.cmd','fasta','fasta.index']:
                tmpfile = '%s.tmp.%s' % (self.baseFile(),tmp)
                if rje.exists(tmpfile) and (not (self.dev() or self.debugging()) or (self.i() > 0 and rje.yesNo('Delete %s?' % tmpfile))):
                    os.unlink(tmpfile)
                    self.printLog('#CLEAN','%s deleted' % tmpfile)
            ## ~ [3b] Summarise sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Summarise'): rje_seqlist.batchSummarise(self,seqfiles,save=True,overwrite=self.force())
            return seqin     # SortNR successful
        except: self.errorLog('Problem during %s.sortNR.' % self.prog()); return False
#########################################################################################################################
    def dipHap(self):   ### Update sequence names with haploid[AB] and diploid, then output along with pri and alt files.
        '''Update sequence names.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Setup input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            overwrite = self.force()
            seqfiles = []
            if self.getStrLC('RunMode') == 'diphapnr':
                seqlist = self.sortNR()
                overwrite = False
            else:
                seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F'])
                seqfiles = [seqlist.getStr('SeqIn')]
            ## ~ [0b] Setup output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            SEQOUT = open(self.baseFile() + '.dipnr.fasta','w')
            PRIOUT = open(self.baseFile() + '.pri.fasta','w')
            ALTOUT = open(self.baseFile() + '.alt.fasta','w')
            seqfiles += [self.baseFile() + '.dipnr.fasta',self.baseFile() + '.pri.fasta',self.baseFile() + '.alt.fasta']
            ### ~ [1] Allocate primary and alternative scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            diplist = []
            haplist = []
            dx = px = ax = sx = 0.0; stot = seqlist.seqNum() * 2
            for seq in seqlist.seqs():
                self.progLog('\r#DIPHAP','Pseudodiploid haplotig assignment: %.1f%%' % (sx/stot)); sx += 100
                hap = rje.matchExp('HAP(\d+)',seqlist.shortName(seq))
                if hap in diplist: haplist.append(hap)
                else: diplist.append(hap)
            for seq in seqlist.seqs():
                self.progLog('\r#DIPHAP','Pseudodiploid haplotig assignment: %.1f%%' % (sx/stot)); sx += 100
                sname = string.split(seqlist.shortName(seq),'_')
                hap = rje.matchExp('HAP(\d+)',seqlist.shortName(seq))
                if hap in haplist:
                    if hap in diplist: haptxt = 'haploidA'; diplist.remove(hap); sname[0] = 'pri%s' % hap; px += 1
                    else: haptxt = 'haploidB'; sname[0] = 'alt%s' % hap; ax += 1
                else: haptxt = 'diploid'; sname[0] = 'pri%s' % hap; dx += 1; px += 1
                sname = string.join(sname,'_')
                SEQOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                if sname[:3] == 'pri':
                    PRIOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                else:
                    ALTOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                #self.printLog('#DIP','%s: %s\n' % (sname,haptxt))
            self.printLog('\r#DIPHAP','Pseudodiploid haplotig assignment complete!')
            SEQOUT.close()
            PRIOUT.close()
            ALTOUT.close()
            self.printLog('#OUT','%s Primary assembly haplotigs (inc. %s Diploid) output to %s.pri.fasta' % (rje.iStr(px),rje.iStr(dx),self.baseFile()))
            self.printLog('#OUT','%s Alternative assembly haplotigs output to %s.alt.fasta' % (rje.iStr(ax),self.baseFile()))
            self.printLog('#OUT','%s Combined Primary+Alternative assembly haplotigs output to %s.dipnr.fasta' % (rje.iStr(ax+px),self.baseFile()))
            ### ~ [2] Summarise sequence files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Summarise'): rje_seqlist.batchSummarise(self,seqfiles,save=True,overwrite=overwrite)
            return True
        except: self.errorLog('Problem during %s.dipHap.' % self.prog()); return False
#########################################################################################################################
    ### <5> ### VecScreen methods                                                                                       #
#########################################################################################################################
    def vecScreen(self):    ### Screen Scaffolds for possible contaminants
        '''
        ## Vector contamination screening [runmode=vecscreen]

        Screen Scaffolds for possible contaminants.

        First, a `blastn` search of the `screendb=FILE` sequences is performed, using the NCBI VecScreen search and
        match strategy (see below). These parameters can be modified using `blastopt=X` and/or `optionfile=FILE`,
        which will be appended to the run command. Alternatively, the `$BASEFILE.vecscreen.blast` file produced can be
        pre-made with different parameters (if `force=F`).

        Results are then parsed into a local hits table and rated according to the strength of a match. This is performed
        iteratively, re-assigning internal matches as "proximal" if they are withing 25 bases of another match (of any
        strength) and re-rating the strength of the match, until no ratings changes occur. Two additional fields are
        added to the local hits table during this process: `MatchPos` and `MatchStr`. Once all assignments have been
        made, segments of the assembly between two matches and/or sequence ends are added to the table as `Suspect`.

        `MatchPos` will have a value of:

        * `Terminal` = within 25 bp of either end of the sequence.
        * `Proximal` = within 25 bp of a vecsreen match (`Weak`, `Moderate` or `Strong`).
        * `Internal` = over 25 bp of a sequence end or vecsreen match.
        * `Suspect` = Segments added as `Suspect`.

        `MatchStr` will have a value of:

        * `Strong` = `Terminal`/`Proximal` match with Score >= 24, or `Internal` match with Score >= 30.
        * `Moderate` = `Terminal`/`Proximal` match with Score 19 to 23, or `Internal` match with Score 25 to 29.
        * `Weak` = `Terminal`/`Proximal` match with Score 16 to 18, or `Internal` match with Score 23 to 24.
        * `Suspect` = Any segment of fewer than 50 bases between two vector matches or between a match and an end.

        ### VecScreen parameters (from NCBI website)

        The VecScreen parameters are pre-set using blastn options: -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000

        VecScreen Match Categories
        Vector contamination usually occurs at the beginning or end of a sequence; therefore, different criteria are
        applied for terminal and internal matches. VecScreen considers a match to be terminal if it starts within 25
        bases of the beginning of the query sequence or stops within 25 bases of the end of the sequence. Matches that
        start or stop within 25 bases of another match are also treated like terminal matches. Matches are categorized
        according to the expected frequency of an alignment with the same score occurring between random sequences.

        Strong Match to Vector
        (Expect 1 random match in 1,000,000 queries of length 350 kb.)
        Terminal match with Score >= 24.
        Internal match with Score >= 30.

        Moderate Match to Vector
        (Expect 1 random match in 1,000 queries of length 350 kb.)
        Terminal match with Score 19 to 23.
        Internal match with Score 25 to 29.

        Weak Match to Vector
        (Expect 1 random match in 40 queries of length 350 kb.)
        Terminal match with Score 16 to 18.
        Internal match with Score 23 to 24.

        Segment of Suspect Origin
        Any segment of fewer than 50 bases between two vector matches or between a match and an end.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            #i#covdb = self.db().addEmptyTable('screencov',['Query','Hit','HitLen','Coverage','CovPC'],['Query','Hit'],log=self.debugging())
            screencov = '{}.screencov.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and rje.checkForFiles(filelist=[screencov],basename='',log=self.log):
                return screencov
            forks = self.getInt('Forks')
            bfile = '%s.vecscreen.blast' % self.baseFile()
            locfile = '%s.Local.tdt' % self.baseFile()
            if not rje.exists(self.getStr('SeqIn')):
                raise IOError('Diploidocus VecScreen mode needs input assembly (seqin=FILE)')
            if not rje.exists(bfile) and not rje.exists(self.getStr('ScreenDB')):
                raise IOError('Diploidocus VecScreen mode needs screendb=FILE set or exisiting %s file' % bfile)
            blast = rje_blast.BLASTRun(log=self.log,cmd_list=['blaste=700','keepblast=T','blastgz=T']+self.cmd_list+['bitscore=F'])
            self.obj['BLAST'] = blast
            mytables = self.db().tables()
            blast.obj['DB'] = self.db()
            blast.setStr({'Type':'blastn','InFile':self.getStr('ScreenDB'),'BlastRes':bfile,'Name':bfile,'DBase':self.getStr('SeqIn')})
            blast.setup(load=False)
            ## ~ [1a] ~ Check existing results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            runblast = self.force() or not rje.exists(locfile)
            if rje.exists(locfile):
                if self.force(): self.printLog('#FORCE','%s found but will be regenerated. (force=T)' % locfile)
                else: self.printLog('#LOCAL','%s found - BLAST will be skipped. (force=F)' % locfile)
            if runblast and not self.force() and blast.checkBLAST():                ### BLAST Results exist
                if blast.getBool('IgnoreDate'): runblast = False       ### Don't check age!
                elif rje.isYounger(blast.getStr('DBase'),blast.getStr('Name')) == blast.getStr('Name') and rje.isYounger(blast.getStr('InFile'),blast.getStr('Name')) == blast.getStr('Name'): runblast = False
            ## ~ [1b] ~ Setup search database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if runblast:
                blast.formatDB(protein=False,force=self.force(),details=self.debugging())
            ## ~ [1c] Input Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F'])
            #i# seqdict will store all the sequences, mapped by shortname
            seqdict = seqin.seqNameDic()

            ### ~ [2] ~ Perform BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # First, a `blastn` search of the `screendb=FILE` sequences is performed, using the NCBI VecScreen search and
            # match strategy (see below). These parameters can be modified using `blastopt=X` and/or `optionfile=FILE`,
            # which will be appended to the run command. Alternatively, the `$BASEFILE.vecscreen.blast` file produced can be
            # pre-made with different parameters (if `force=F`).
            if runblast:
                veccmd = ' -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -searchsp 1750000000000'
                command = blast.blastPath() + blast.getStr('Type') + veccmd
                command += ' -evalue %e' % blast.getNum('E-Value')
                command += ' -query %s' % blast.getStr('InFile')
                command += ' -db %s' % blast.getStr('DBase')
                command += ' -out %s' % blast.getStr('Name')
                if forks > 1: command += ' -num_threads %d' % forks
                if blast.getStrLC('BLASTOpt'): command = '%s %s' % (command,blast.getStr('BLASTOpt'))
                if rje.exists(blast.getStr('OptionFile')):
                    for line in open(blast.getStr('OptionFile'),'r').readlines(): command = '%s %s' % (command,rje.chomp(line))
                blast.str['BLASTCmd'] = command
                self.printLog('\r#SYS',command)
                os.system(command)
                if not blast.checkBLAST(): raise IOError('Problem with VecSreen BLAST results file "%s"' % bfile)
            ## ~ [2a] Read Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vecdb = None
            if rje.exists(locfile) and not self.force():
                blast.db().list['Tables'] = mytables
                vecdb = self.db().addTable(locfile,['Query','Hit','AlnID'],name='vecscreen')
                blast.formatTables()
                vecdb.dataFormat({'QRank':'int'})
            else:
                keepaln = False     #!# Might want to add alignment output later
                rmblast = not blast.getBool('KeepBLAST')
                #self.debug('Delete BLAST: %s' % rmblast)
                blast.readBLAST(resfile=blast.getStr('BlastRes'),clear=True,gablam=False,local=True,keepaln=keepaln,unlink=rmblast)
                blast.formatTables()
            ## ~ [2b] Add eFDR calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                vecdb = self.db('Local')    #?# Should this be a copy?
                #i# Keys: ['Query','Hit','AlnID']
                #i# Fields: ['Query','Hit','AlnID','Score','Expect','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd']
                vecdb.keepFields(['Query','Hit','AlnID','Score','Expect','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd'])
                self.printLog('#EFDR','Calculating expected FDR (eFDR) based on Expect scores and observed hits')
                vecdb.rankFieldByIndex('Query','Expect',newfield='QRank',rev=True,absolute=True,lowest=False,warn=True,highest=True)
                vecdb.makeField('Expect/QRank','eFDR')
                vecdb.saveToFile()
                vecdb.setStr({'Name':'vecscreen'})
            ## ~ [2c] Filter Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vecdb.dataFormat({'eFDR':'float'})
            if self.getNum('eFDR') > 0:
                #self.printLog('#EFDR','Filtering entries with eFDR<{}'.format(self.getNum('eFDR')))
                #vecdb.dropEntries('eFDR<{}'.format(self.getNum('eFDR')))
                ex = 0.0; etot = vecdb.entryNum()
                for ekey, data in vecdb.data().items():
                    self.progLog('\r#EFDR','Filtering entries with eFDR>{}: {:.2f}%'.format(self.getNum('eFDR'),ex/etot)); ex += 100.0
                    if data['eFDR'] > self.getNum('eFDR'): vecdb.data().pop(ekey)
                    #else: self.debug(data)
                self.printLog('\r#EFDR','Filtered entries with eFDR<{}: {} -> {} entries'.format(self.getNum('eFDR'),rje.iStr(etot),rje.iStr(vecdb.entryNum())))
            if self.getNum('MinVecHit') > 0:
                self.printLog('#HITLEN','Filtering entries with Length<{}'.format(self.getInt('MinVecHit')))
                vecdb.dropEntries('Length<{}'.format(self.getInt('MinVecHit')))

            ### ~ [3] Process Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            proxbp = 25
            suspectbp = 50
            # Results are then parsed into a local hits table and rated according to the strength of a match. This is performed
            # iteratively, re-assigning internal matches as "proximal" if they are withing 25 bases of another match (of any
            # strength) and re-rating the strength of the match, until no ratings changes occur. Two additional fields are
            # added to the local hits table during this process: `MatchPos` and `MatchStr`. Once all assignments have been
            # made, segments of the assembly between two matches and/or sequence ends are added to the table as `Suspect`.
            #
            # `MatchPos` will have a value of:
            #
            # * `Terminal` = within 25 bp of either end of the sequence.
            # * `Proximal` = within 25 bp of a vecsreen match (`Weak`, `Moderate` or `Strong`).
            # * `Internal` = over 25 bp of a sequence end or vecsreen match.
            # * `Suspect` = Segments added as `Suspect`.
            #
            # `MatchStr` will have a value of:
            #
            # * `Strong` = `Terminal`/`Proximal` match with Score >= 24, or `Internal` match with Score >= 30.
            # * `Moderate` = `Terminal`/`Proximal` match with Score 19 to 23, or `Internal` match with Score 25 to 29.
            # * `Weak` = `Terminal`/`Proximal` match with Score 16 to 18, or `Internal` match with Score 23 to 24.
            # * `Suspect` = Any segment of fewer than 50 bases between two vector matches or between a match and an end.
            ## ~ [3a] Setup table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Add Strand and make start/end consistent
            vecdb.addField('Strand',evalue='+')
            for ventry in vecdb.entries():
                if ventry['SbjStart'] > ventry['SbjEnd']: (ventry['SbjStart'],ventry['SbjEnd'],ventry['Strand']) = (ventry['SbjEnd'],ventry['SbjStart'],'-')
            ## Add MatchPos and MatchType
            vecdb.addField('MatchPos',evalue='Internal')
            vecdb.addField('MatchStr',evalue='None')
            #i# MatchStart and MatchEnd have the Subject Start/End if a Match
            vecdb.addField('MatchStart',evalue=0)
            vecdb.addField('MatchEnd',evalue=0)
            #i# Terminal and Internal give the match ratings depending on position
            vecdb.addField('Terminal',evalue='None')
            vecdb.addField('Internal',evalue='None')
            ## ~ [3b] Filter < Weak matches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for ventry in vecdb.entries():
                if ventry['Score'] >= 30: ventry['Internal'] = 'Strong'
                elif ventry['Score'] >= 25: ventry['Internal'] = 'Moderate'
                elif ventry['Score'] >= 23: ventry['Internal'] = 'Weak'
                if ventry['Score'] >= 24: ventry['Terminal'] = 'Strong'
                elif ventry['Score'] >= 19: ventry['Terminal'] = 'Moderate'
                elif ventry['Score'] >= 16: ventry['Terminal'] = 'Weak'
            vecdb.dropEntriesDirect('Terminal',['None'])
            ## ~ [3c] Re-order based on Hit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# New key is for determining proximity to other matches.
            vecdb.newKey(['Hit','SbjStart','SbjEnd','Query','AlnID'])
            ## ~ [3d] Add Terminal ratings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x# vecdb.indexReport('Hit')
            for sname in vecdb.index('Hit'):
                hitlen = seqin.seqLen(seqdict[sname])
                for ventry in vecdb.indexEntries('Hit',sname):
                    if ventry['SbjStart'] <= proxbp: ventry['MatchPos'] = 'Terminal'
                    elif hitlen - ventry['SbjEnd'] < proxbp: ventry['MatchPos'] = 'Terminal'
                    if ventry['MatchPos'] == 'Terminal': ventry['MatchStr'] = ventry['Terminal']
                    else: ventry['MatchStr'] = ventry['Internal']
                    if ventry['MatchStr'] != 'None': (ventry['MatchStart'],ventry['MatchEnd']) = (ventry['SbjStart'],ventry['SbjEnd'])
            ## ~ [3e] Cycle and rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cycling = True
            while cycling:
                cycling = False     #i# Set to True if ratings updated.
                sx = 0.0; stot = len(vecdb.index('Hit'))
                for sname in vecdb.index('Hit'):
                    vhits = vecdb.indexEntries('Hit',sname)[0:]
                    sstep = 100.0 / len(vhits)
                    for ventry in vhits[0:]:
                        self.progLog('\r#MATCH','Updating proximal matches: %.1f%%' % (sx/stot)); sx += sstep
                        if ventry['MatchPos'] in ['Terminal','Proximal']: continue
                        #self.bugPrint(vecdb.entrySummary(ventry))
                        #i# Identify a match within the required distance #i#
                        proximal = False
                        proxstart = max(1,ventry['SbjStart'] - proxbp + 1)
                        proxend = ventry['SbjEnd'] + proxbp - 1
                        for vhit in vhits[0:]:
                            if vhit['MatchEnd'] < proxstart: vhits.remove(vhit)
                            if vhit['MatchStart'] <= proxend:
                                #self.debug(vecdb.entrySummary(vhit,collapse=True))
                                proximal = True
                                break
                        if proximal:
                            ventry['MatchPos'] = 'Proximal'
                            ventry['MatchStr'] = ventry['Terminal']
                            (ventry['MatchStart'],ventry['MatchEnd']) = (ventry['SbjStart'],ventry['SbjEnd'])
                            cycling = True
                if cycling: self.printLog('\r#MATCH','Updating proximal matches: still cycling.')
                else: self.printLog('\r#MATCH','Updating proximal matches complete.')
            vecdb.dropEntriesDirect('MatchStr',['None'])
            ## ~ [3g] Add suspect regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sx = 0.0; stot = len(vecdb.index('Hit')); suspectx = 0
            for sname in vecdb.index('Hit'):
                self.progLog('\r#REGION','Identifying suspect regions: %.1f%%' % (sx/stot)); sx += 100.0
                hitlen = seqin.seqLen(seqdict[sname])
                matches = []    # (start,end) tuple list
                vhits = vecdb.indexEntries('Hit',sname)[0:]
                for ventry in vhits[0:]:
                    matches.append((ventry['MatchStart'],ventry['MatchEnd']))
                matches.sort()
                matches = rje.collapseTupleList(matches)
                suspects = rje.invertTupleList(matches,minx=1,maxx=hitlen)
                for (starti,endi) in suspects:
                    if (endi - starti) < suspectbp:
                        vecdb.addEntry({'Query':'Suspect','Hit':sname,'AlnID':0,'SbjStart':starti,'SbjEnd':endi,
                                        'Score':0,'Identity':0,'Length':endi - starti + 1,
                                        'MatchPos':'Suspect','MatchStr':'Suspect'})
                        suspectx += 1
            self.printLog('\r#REGION','Identifying suspect regions complete: %s regions' % rje.iStr(suspectx))

            ### ~ [4] Save Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vecdb.dropFields(['MatchStart','MatchEnd','Terminal','Internal'])
            vecdb.saveToFile()

            ### ~ [5] Final Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Need to look at numbers per contaminant to help detect false positives.
            vecdb.indexReport('Query')
            ## ~ [5a] Calculate percentage coverage per query/hit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            covdb = self.db().addEmptyTable('screencov',['Query','Hit','HitLen','Coverage','CovPC'],['Query','Hit'],log=self.debugging())
            qhcov = {}  # Dictionary of { (Qry,Hit) : [(hitstart,hitend) , (...)]
            vx = 0.0; vtot = vecdb.entryNum()
            for ventry in vecdb.entries():
                self.progLog('\r#COV','Calculating percentage coverage per query/hit pair: {:.2f}%'.format(vx/vtot)); vx += 50.0
                qry = ventry['Query']
                hit = ventry['Hit']
                qh = (qry,hit)
                if qh not in qhcov: qhcov[qh] = []
                qhcov[qh].append((ventry['SbjStart'],ventry['SbjEnd']))
            qx = 0.0; qtot = len(qhcov)
            for qh, covlist in qhcov.items():
                self.progLog('\r#COV','Calculating percentage coverage per query/hit pair: {:.2f}%'.format(50.0+(qx/qtot))); qx += 50.0
                covsum = 0
                for (ci,cj) in rje.collapseTupleList(covlist,joindistance=1,overlaps=False):
                    covsum += cj - ci + 1
                hlen = seqin.seqLen(seqdict[qh[1]])
                covdb.addEntry({'Query':qh[0], 'Hit':qh[1], 'HitLen':hlen, 'Coverage':covsum, 'CovPC':rje.dp(100.0 * covsum /  hlen,2)})
            self.printLog('\r#COV','Calculation of percentage coverage per query/hit pair complete')
            ## ~ [5b] Calculate total contamination coverage for each Hit ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hcov = {}
            qx = 0.0; qtot = len(qhcov)
            for qh, covlist in qhcov.items():
                self.progLog('\r#COV','Calculating percentage coverage per scaffold: {:.2f}%'.format((qx/qtot))); qx += 50.0
                if qh[1] not in hcov: hcov[qh[1]] = []
                hcov[qh[1]] += covlist
            qx = 0.0; qtot = len(hcov)
            for hit, covlist in hcov.items():
                self.progLog('\r#COV','Calculating percentage coverage per scaffold: {:.2f}%'.format(50.0+(qx/qtot))); qx += 50.0
                covsum = 0
                for (ci,cj) in rje.collapseTupleList(covlist,joindistance=1,overlaps=False):
                    covsum += cj - ci + 1
                hlen = seqin.seqLen(seqdict[hit])
                covdb.addEntry({'Query':'TOTAL', 'Hit':hit, 'HitLen':hlen, 'Coverage':covsum, 'CovPC':rje.dp(100.0 * covsum /  seqin.seqLen(seqdict[hit]),2)})
            self.printLog('\r#COV','Calculation of percentage coverage per scaffold complete')
            covdb.saveToFile()
            return screencov

        except IOError:
            self.errorLog('Diploidocus.vecSreen() error')
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    ### <6> ### Genome Size and Purge Haplotigs methods                                                                 #
#########################################################################################################################
# echo EOG MPileup Depth > busco.depth.txt
# for REGION in $( grep Complete $BUSCO | cut -f 1,3-5 | sed 's/\t/:/g'); do
#   EOG=$(echo $REGION | awk -F ':' '{print $1;}')
#   LOC=$(echo $REGION | awk -F ':' '{print $2;}')
#   START=$(echo $REGION | awk -F ':' '{print $3;}')
#   END=$(echo $REGION | awk -F ':' '{print $4;}')
#   echo $REGION
#   DEP=$(samtools view -h -F 0x100 $BAM ${LOC}:${START}-${END} | samtools depth -  | awk -v pos="$START" '$2 >= pos' | awk -v pos="$END" '$2 <= pos' | awk '{print $3;}' | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')
#   PILE=$(samtools view -b -h -F 0x100 $BAM ${LOC}:${START}-${END} | samtools mpileup -BQ0 - | awk -v pos="$START" '$2 >= pos' | awk -v pos="$END" '$2 <= pos' |  awk '{print $4;}' | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')
#   echo $EOG $PILE $DEP | tee -a busco.depth.txt
# done
# echo 'Mpileup:'
# awk '{print $2;}' busco.depth.txt | sort | uniq -c | sort -nr | head
# echo 'Depth:'
# awk '{print $3;}' busco.depth.txt | sort | uniq -c | sort -nr | head
#########################################################################################################################
    def genomeSizeFromModeFile(self,dephist,scdepth=False):   ### Uses read depth from BUSCO single copy genes to predict genome size
        '''
        Uses read depth from BUSCO single copy genes to predict genome size.
        >> dephist:str = *.dephist.tdt file previously generated from BUSCO analysis
        >> scdepth:bool [False] = Whether to return single copy read depth only (w/o Genome Size prediction)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            #i# NOTE: readbp must have been pre-calculated for this method. Use genomeSize() if not.
            readbp = self.getInt('ReadBP')
            if scdepth and not self.getInt('GenomeSize') and readbp: scdepth = False
            ## ~ [1a] ~ Load table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# headers=['Method','X','n','Modes'],
            mdb = db.addTable(dephist,mainkeys=['Method','X'],expect=True,name='dephist')
            mdb.dataFormat({'X':'int','n':'int','Modes':'int'})
            if depmethod not in mdb.index('Method'):
                self.printLog('#METHOD','Previous mode calculation did not use {}'.format(depmethod))
                return False
            ### ~ [2] Calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Generate mode and mode of modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            modelist = []
            for entry in mdb.indexEntries('Method',depmethod):
                n = entry['Modes']
                X = entry['X']
                modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'ModeOfModes':modelist[0][1]})
            self.printLog('#MODE','BUSCO SC mode of samtools {} modal read depth: {}X'.format(depmethod,self.getInt('ModeOfModes')))
            for entry in mdb.indexEntries('Method',depmethod):
                n = entry['n']
                X = entry['X']
                modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'BUSCOMode':modelist[0][1]})
            self.printLog('#MODE','BUSCO SC combined samtools {} modal read depth: {}X'.format(depmethod,self.getInt('BUSCOMode')))
            if scdepth: return self.getInt('BUSCOMode')
            ### ~ [3] Calculate genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            estgensize = int(0.5+(float(readbp) / self.getInt('BUSCOMode')))
            self.setInt({'EstGenomeSize':estgensize})
            if not self.getInt('GenomeSize'): self.setInt({'GenomeSize':estgensize})
            self.printLog('#GSIZE','Estimated genome size ({} at {}X): {}'.format(rje_seqlist.dnaLen(readbp,dp=0,sf=3),self.getInt('BUSCOMode'),rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))
            return self.getInt('BUSCOMode')
        except:
            self.errorLog('Diploidocus.genomeSizeFromModeFile() error')
            return False
#########################################################################################################################
    def genomeSize(self,scdepth=False):   ### Uses read depth from BUSCO single copy genes to predict genome size
        '''
        Uses read depth from BUSCO single copy genes to predict genome size.
        >> scdepth:bool [False] = Whether to return single copy read depth only (w/o Genome Size prediction)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            forker = self.obj['Forker']
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            readbp = self.getInt('ReadBP')
            if scdepth and not self.getInt('GenomeSize') and (readbp or self.list['Reads']): scdepth = False
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj()
            seqdict = seqlist.seqNameDic()
            ## ~ [1a] Check Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bamfile = self.getStr('BAM')
            if not self.getStrLC('BAM'): bamfile = self.baseFile(strip_path=True) + '.bam'
            if not rje.exists(bamfile): raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            self.printLog('#BAM',bamfile)
            busco = self.getStr('BUSCO')
            if not self.getStrLC('BUSCO'): busco = 'full_table_{}.busco.tsv'.format(self.baseFile(strip_path=True))
            if not rje.exists(busco): raise IOError('Cannot find BUSCO full table "{}" (busco=FILE)'.format(busco))
            self.printLog('#BUSCO',busco)
            if not readbp and not self.list['Reads'] and not scdepth: raise IOError('Cannot find any read files (reads=FILELIST) and no readbp=X set.')
            ## ~ [1b] Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not readbp and self.list['Reads'] and not scdepth:
                for rfile in self.list['Reads']:
                    countfile = '{}.basecount.txt'.format(rfile)
                    if rje.exists(countfile) and not self.force():
                        filebp = int(rje.chomp(open(countfile,'r').readline()))
                        self.printLog('#READBP','%s: %s bases in %s' % (countfile,rje.iStr(filebp),rfile))
                        readbp += filebp
                        continue
                    self.progLog('\r#READBP','Counting bases in {}...'.format(rfile))
                    gzip = rfile.endswith('.gz')
                    fastq = rfile.endswith('q.gz') or rfile.endswith('q')
                    if gzip and fastq:
                        filebp = int(os.popen("zcat %s | grep '^+$' -B1 | grep -v '^+$' | grep -v '^--' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    elif fastq:
                        filebp = int(os.popen("grep '^+$' -B1 %s | grep -v '^+$' | grep -v '^--' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    elif gzip:
                        filebp = int(os.popen("zcat %s | grep -v '^>' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    else:
                        filebp = int(os.popen("grep -v '^>' %s | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    self.printLog('#READBP','%s bases counted in %s' % (rje.iStr(filebp),rfile))
                    open(countfile,'w').write('%d\n' % filebp)
                    readbp += filebp
                self.setInt({'ReadBP':readbp})
            if not scdepth: self.printLog('#READBP','Total base count: %s' % (rje.iStr(readbp)))
            ## ~ [1c] Load BUSCO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ext = rje.delimitExt(db.getStr('Delimit'))
            scfiles = ['.busco{}.{}'.format(depmethod,ext), '.buscodep.{}'.format(ext), '.dephist.{}'.format(ext)]
            if rje.checkForFiles(filelist=scfiles,basename=db.baseFile(),log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.'):
                runok = self.genomeSizeFromModeFile('{}.dephist.{}'.format(db.baseFile(),ext),scdepth)
                if runok: return runok
                else: self.printLog('#SCDEP','Unable to establish modal SC depth from {}'.format(scfiles[-1]))
            fullhead = ['BuscoID','Status','Contig','Start','End','Score','Length']
            bdb = db.addTable(busco,mainkeys='auto',headers=fullhead,expect=True,name='busco{}'.format(depmethod))
            bdb.dropEntriesDirect('Status',['Complete'],inverse=True)
            dephead = ['Method','BuscoID','X','n']
            depdb = db.addEmptyTable('buscodep',dephead,['Method','BuscoID','X'],log=self.debugging())
            ## ~ [1d] Temp directory for forked depths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpdir = rje.makePath(self.getStr('TmpDir'),wholepath=False)
            if not rje.exists(tmpdir): rje.mkDir(self,tmpdir)
            ## ~ [1e] Check programs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.popen('samtools --version').read():
                self.printLog('#SYS',' '.join(os.popen('samtools --version').read().split()))
            else:
                raise IOError('Cannot open samtools: check installation and/or module loading')

            ### ~ [2] Cycle through and fork out the depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []
            for seqname in bdb.index('Contig'):
                if seqname not in seqdict:
                    bx = rje.iLen(bdb.index('Contig')[seqname])
                    self.warnLog('"%s" (%s BUSCO-complete sequences) not found in "%s". May have been removed in previous purge cycle.' % (seqname,bx,seqin))
                    missing.append(seqname)
            if missing:
                bdb.dropEntriesDirect('Contig',missing,log=False)
            ## ~ [2a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cleanup = 0; skipped = 0
            forker.list['ToFork'] = []
            for bentry in bdb.entries():
                # if bentry['Contig'] not in seqdict:
                #     self.warnLog('BUSCO-complete sequence "%s" not found in "%s". May have been removed in previous purge cycle.' % (bentry['Contig'],seqin))
                #     bdb.dropEntry(bentry)
                #     continue
                tmpfile = '{}{}.{}.{}.tmp'.format(tmpdir,basefile,bentry['BuscoID'],depmethod)
                if rje.exists(tmpfile):
                    #?# Add checking for completeness #?#
                    if not self.force() and len(open(tmpfile,'r').readline().split()) == 2: skipped += 1; continue
                    else: os.unlink(tmpfile); cleanup += 1
                if depmethod == 'mpileup':
                    forker.list['ToFork'].append("samtools view -b -h -F 0x100 %s %s:%s-%s | samtools mpileup -BQ0 - 2> /dev/null | awk -v pos=\"%s\" '$2 >= pos' | awk -v pos=\"%s\" '$2 <= pos' | awk '{print $4;}' | sort | uniq -c | sort -nr > %s" % (bamfile,bentry['Contig'],bentry['Start'],bentry['End'],bentry['Start'],bentry['End'],tmpfile))
                else:
                    forker.list['ToFork'].append("samtools view -h -F 0x100 %s %s:%s-%s | samtools depth - | awk -v pos=\"%s\" '$2 >= pos' | awk -v pos=\"%s\" '$2 <= pos' | awk '{print $3;}' | sort | uniq -c | sort -nr > %s" % (bamfile,bentry['Contig'],bentry['Start'],bentry['End'],bentry['Start'],bentry['End'],tmpfile))
            self.printLog('#BUSCO','{} Complete BUSCO gene regions queued for forking ({} existing files deleted); {} existing results skipped'.format(rje.iLen(forker.list['ToFork']),rje.iStr(cleanup),rje.iStr(skipped)))
            ## ~ [2b] Fork out depth analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for forkcmd in forker.list['ToFork']:
                        self.printLog('#SYS',forkcmd)
                        os.system(forkcmd)
                elif forker.run():
                    self.printLog('#FORK','Forking of BUSCO depth analysis completed.')
                else:
                    try:
                        self.errorLog('Samtools forking did not complete',printerror=False,quitchoice=True)
                    except:
                        raise RuntimeError('Samtools forking did not complete')

            ### ~ [3] Read in and process depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            depmodes = {}   # Dictionary of depth modes
            depcounts = {}  # Dictionary of depth counts
            ## ~ [3a] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb.addField('Mode')
            for bentry in bdb.entries():
                tmpfile = '{}{}.{}.{}.tmp'.format(tmpdir,basefile,bentry['BuscoID'],depmethod)
                try:
                    if not rje.exists(tmpfile):
                        raise IOError('Cannot find {}'.format(tmpfile))
                    #?# Add full processing and calculation of mean coverage?
                    try:
                        bentry['Mode'] = int(open(tmpfile,'r').readline().split()[1])
                    except:
                        self.errorLog('Something has gone wrong with "%s"' % tmpfile)
                        raise
                    if bentry['Mode'] not in depmodes: depmodes[bentry['Mode']] = 0
                    depmodes[bentry['Mode']] += 1
                    for dline in open(tmpfile,'r').readlines():
                        data = rje.chomp(dline).split()
                        n = int(data[0])
                        X = int(data[1])
                        depdb.addEntry({'Method':depmethod,'BuscoID':bentry['BuscoID'],'n':n,'X':X})
                        if X not in depcounts: depcounts[X] = 0
                        depcounts[X] += n
                except:
                    self.errorLog('Samtools depth result processing error',quitchoice=True)
                    continue
            ## ~ [3b] Generate mode and mode of modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            modelist = []
            for X, n in depmodes.items(): modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'ModeOfModes':modelist[0][1]})
            self.printLog('#MODE','BUSCO SC mode of samtools {} modal read depth: {}X'.format(depmethod,self.getInt('ModeOfModes')))
            modelist = []
            for X, n in depcounts.items(): modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'BUSCOMode':modelist[0][1]})
            self.printLog('#MODE','BUSCO SC combined samtools {} modal read depth: {}X'.format(depmethod,self.getInt('BUSCOMode')))
            ## ~ [3c] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb.saveToFile()
            #!# Use the depdb table to make a box plot of the normalised depth counts?
            depdb.saveToFile()
            ## ~ [3d] Table of total and mode counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            depdb.compress(['Method','X'],rules={'n':'sum'})
            depdb.dropField('BuscoID')
            depdb.addField('Modes',evalue=0)
            for X, n in depmodes.items(): depdb.data((depmethod,X))['Modes'] = n
            depdb.rename('dephist')
            depdb.saveToFile()
            if scdepth: return self.getInt('BUSCOMode')

            ### ~ [4] Calculate genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            estgensize = int(0.5+(float(readbp) / self.getInt('BUSCOMode')))
            self.setInt({'EstGenomeSize':estgensize})
            if not self.getInt('GenomeSize'): self.setInt({'GenomeSize':estgensize})
            self.printLog('#GSIZE','Estimated genome size ({} at {}X): {}'.format(rje_seqlist.dnaLen(readbp,dp=0,sf=3),self.getInt('BUSCOMode'),rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))

            return self.getInt('BUSCOMode')
        except SystemExit: raise    # Child
        except:
            self.errorLog('Diploidocus.genomeSize() error')
            return False
#########################################################################################################################
    def assemblyMinimap(self):  ### Performs assembly versus self minimap2 search
        '''
        Performs assembly versus self minimap2 search
        :return: True/False
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

            #!# Modify the code below #!#

            #!# NOTE: This isn't currently used for anything. #!# Purge_haplotigs hits used #!#

            #mapopt=CDICT    : Dictionary of updated minimap2 options [N:250,p:0.0001,x:asm20]
            paf_defaults = {'N':'250','p':'0.0001','x':'asm20'}
            if self.getStr('Mapper') in ['minimap','minimap2']:
                pafin = '%s.paf' % self.baseFile()
                paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'basefile=%s' % self.baseFile(),'seqin=%s' % self.getStr('QueryDB'),'reference=%s' % self.getStr('SearchDB')])
                paf.setInt({'MinLocLen':self.getInt('LocalMin')})
                paf.setNum({'MinLocID':self.getNum('LocalIDMin')})
                defaults = paf_defaults
                paf.dict['MapOpt'] = rje.combineDict(defaults,paf.dict['MapOpt'],overwrite=True)
                paf.setup()
                if not rje.exists(pafin) or self.force():
                    rje.backup(self,pafin)
                    paf.minimap2()
                    #!# Add localaln=T to PAF #!#
                paf.run()
                paf.db('hitunique').rename('unique')
                #!# Change code to handle Qry and AlnNum #!#
                paf.db('unique').renameField('Qry','Query')
                paf.db('unique').renameField('AlnNum','AlnID')
                paf.db('unique').saveToFile()
                self.obj['BLAST'] = rje_blast.BLASTRun(log=self.log,cmd_list=['blastf=F']+self.cmd_list+['checkblast=F'])
                self.obj['BLAST'].obj['DB'] = self.obj['DB'] = paf.db()
                self.debug('PAF to GABLAM conversion complete')
                if self.getStr('NRStat') == 'OrderedAlnID':
                    self.printLog('#NRSTAT','NRStat updated to AlnID for minimap2 run.')
                    self.setStr({'NRStat':'AlnID'})
                return True
        except:
            self.errorLog('Diploidocus.assemblyMinimap() error')
            return False
#########################################################################################################################
    def longreadMinimap(self):  ### Performs long read versus assembly minimap2 and converts to BAM file
        '''
        Performs long read versus assembly minimap2 and converts to BAM file
        :return: bamfile/None
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['PAF'] = paf = rje_paf.PAF(self.log, self.cmd_list)
            bamfile = self.baseFile() + '.bam'
            mmv = os.popen('%s --version' % paf.getStr('Minimap2')).read()
            if not mmv: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting (%s)' % paf.getStr('Minimap2'))
            rje.backup(self,bamfile)

            ### ~ [2] ~ Generate individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Check these BAM files have headers! #!#
            bamlist = []; rx = 0
            if not self.list['ReadType']:
                self.warnLog('Read Type not given (pb/ont): check readtype=LIST. Will use "ont".')
            elif len(self.list['ReadType']) == 1 and len(self.list['Reads']) != 1:
                self.printLog('#READS','Using "%s" as read type for all long reads' % self.list['ReadType'][0])
            elif len(self.list['ReadType']) != len(self.list['Reads']):
                self.warnLog('reads=FILELIST vs readtype=LIST length mismatch: check readtype=LIST. Will cycle if needed.')
            for readfile in self.list['Reads']:
                if not rje.exists(readfile): raise IOError('Read file "{}" not found!'.format(readfile))
                if self.list['ReadType']:
                    try: rtype = self.list['ReadType'][rx]; rx +=1
                    except: rtype = self.list['ReadType'][0]; rx = 1
                    if rtype not in ['ont','pb']:
                        self.warnLog('Read Type "%s" not recognised (pb/ont): check readtype=LIST. Will use "ont".' % rtype)
                        rtype = 'ont'
                else: rtype = 'ont'
                prefix = '{}.{}'.format(rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(readfile,strip_path=True))
                maplog = '{}.log'.format(prefix)
                ## ~ [2a] Make SAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                maprun = '{} -t {} --secondary=no -o {}.sam -L -ax map-{} {} {}'.format(paf.getStr('Minimap2'),self.threads(),prefix,rtype,self.getStr('SeqIn'),readfile)
                # self.printLog('#SYS',maprun)
                # if self.v() < 1:
                #     mapsys = subprocess.Popen(maprun, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #     stdout, stderr = mapsys.communicate()
                #     if stderr:
                #         logline = rje.chomp(stderr).split('\n')[-1]
                #         self.printLog('#SYSEND',logline)
                #     else: self.printLog('#SYSEND','No STDERR output produced during minimap2 mapping!')
                # else:
                #     #i# Generally useful to see minimap2 progress
                #     os.system(maprun)
                logline = self.loggedSysCall(maprun,maplog,append=False)
                #!# Add check that run has finished #!#
                ## ~ [2b] Converting SAM to BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #sam2bam = 'samtools view -bo {}.tmp.bam -@ {} -S {}.sam'.format(prefix,self.threads()-1,prefix)
                self.printLog('#BAM','Converting SAM to BAM. Using a single thread due to past issues of missing data.')
                sam2bam = 'samtools view -bo {}.tmp.bam -S {}.sam'.format(prefix,prefix)
                logline = self.loggedSysCall(sam2bam,maplog,append=True,nologline='No stdout from sam2bam')
                #!# Add check that run has finished #!#
                ## ~ [2c] Sorting BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.printLog('#BAM','Sorting BAM file.')
                sortbam = '{}.bam'.format(prefix)
                bamsort = 'samtools sort -@ {} -o {}.bam -m 6G {}.tmp.bam'.format(self.threads()-1,prefix,prefix)
                logline = self.loggedSysCall(bamsort,maplog,append=True)
                #!# Add check that run has finished #!#
                if not rje.exists(sortbam): raise IOError('Sorted BAM file "%s" not generated' % sortbam)
                os.unlink('{}.sam'.format(prefix))
                os.unlink('{}.tmp.bam'.format(prefix))
                bamlist.append(sortbam)

            ### ~ [3] ~ Merge individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(bamlist) > 1:
                # samtools merge - merges multiple sorted input files into a single output.
                bammerge = 'samtools merge -@ {} {} {}'.format(self.threads()-1,bamfile,' '.join(bamlist))
                logline = self.loggedSysCall(bammerge,append=True)
                if not rje.exists(bamfile): raise IOError('Merged BAM file "%s" not generated' % bamfile)
                for sortbam in bamlist: os.unlink(sortbam)
            else: os.rename(bamlist[0],bamfile)
            ## ~ [3a] Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            return bamfile
        except:
            self.errorLog('Diploidocus.longreadMinimap() error')
            return None
#########################################################################################################################
    def diploidocus(self):   ### Combines purge_haplotigs data with other stats to perform assembly filtering
        '''
        Diploidocus builds on the PurgeHaplotigs classifications to use the depth bins, KAT assembly kmer frequencies and
        BUSCO results to reclassify scaffolds and partition them into:

        * `*.diploidocus.fasta` = the scaffolds kept for the next round of PurgeHap
        * `*.core.fasta` = the same set of scaffolds, minus repeats
        * `*.junk.fasta` = low coverage scaffolds, removed as junk
        * `*.quarantine.fasta` = putative haplotigs, removed from the assembly but retained for reference.

        **NOTE:** PurgeHaplotigs was not using zero coverage bases in its percentages. This is now fixed by Diploidocus.

        The first step in the classification decision process is to load the data and reduce it to the core data frame
        for analysis. `SeqName` will be used for rownames to help xref the full data.

        :return:
        '''
        #!# Will need to split this up at some point
        return self.purgeHaplotigs()
#########################################################################################################################
    def purgeHaplotigs(self):   ### Combines purge_haplotigs data with other stats to perform assembly filtering
        '''
        The Pilon-polished genome underwent a final scaffold cleanup to generate a high-quality core assembly, remove
        low-coverage artefacts and haplotig sequences, and annotate remaining scaffolds with potential issues.



        NOTE: The current method should be divided into submethods and renamed as the primary DiploidocusHocusPocus
        method that runs everything!

        NOTE: Purge_Haplotigs only calculates coverage based on non-zero coverage bases!
        Add includegaps=T/F and mingap=INT settings to include gaps in percentages (including lowcov) or not

        NOTE: Gaps are still not considered for the low coverage (Median depth) filter - set this to 0 if you are
        worried about retention of very gap-dense regions.

        NOTE: After adjustment, this additional filter may not be required!

        :return:
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## ~ [1a] IO names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            basefile = self.getStr('Basefile')
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
            seqdict = seqlist.seqNameDic()
            bamfile = self.getStr('BAM')
            purgemode = self.getStrLC('PurgeMode')
            ## ~ [1b] Programs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Make more nuanced. Add a partial=T/F mode, which can run with some of the programs? Or useX settings?
            if not os.popen('purge_haplotigs hist 2>&1').read():
                self.warnLog('Cannot run "purge_haplotigs hist": check installation or pre-generation of files')
            for program in ['kat','samtools','pileup.sh']:
                if not os.popen('{} --version 2>&1'.format(program)).read():
                    self.warnLog('Cannot run "{} --version": check installation or pre-generation of files'.format(program))
            ## ~ [1c] BAM file check/generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('BAM'):
                bamfile = self.baseFile(strip_path=True) + '.bam'
                if not rje.exists(bamfile):
                    self.printLog('#BAM','Cannot find BAM file "{}": generating'.format(bamfile))
                    bamfile = self.longreadMinimap()
            if not rje.exists(bamfile): raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            self.setStr({'BAM':bamfile})
            baifile = '{}.bai'.format(bamfile)
            #i# NOTE: If Diploidocus keeps remaking files, switch ignoredate=T
            rje.checkForFiles(filelist=[bamfile,baifile],basename='',log=self.log,cutshort=False,ioerror=False)
            if self.needToRemake(baifile,bamfile):
                makebai = 'samtools index -b {} {}.bai'.format(bamfile,bamfile)
                logline = self.loggedSysCall(makebai,append=True)
                #os.system('samtools index -b {} {}.bai'.format(bamfile,bamfile))

            ### ~ [2] ~ Run PurgeHapolotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Establish SC read depth using samtools ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # scdepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
            scdepth = self.getInt('SCDepth')
            if self.getInt('SCDepth'):
                self.printLog('#SCDEP','Using loaded single copy read depth = {}X'.format(scdepth))
            else:
                scdepth = self.genomeSize(scdepth=True)
                self.printLog('#SCDEP','Using BUSCO-derived single copy read depth = {}X'.format(scdepth))
                if not scdepth: raise ValueError('Failed to establish SC read depth')
            ## ~ [2b] ~ Setup purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
            if self.getInt('PHLow') <= 0: self.setInt({'PHLow': int(float(scdepth)/4.0) })
            phlow = self.getInt('PHLow')
            # phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
            if self.getInt('PHMid') <= 0:
                dupdepth = scdepth/2.0
                self.setInt({'PHMid': int(1.5 * dupdepth) })
            phmid = self.getInt('PHMid')
            # phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
            if self.getInt('PHHigh') <= 0:
                self.setInt({'PHHigh': scdepth * 2 })
            phhigh = self.getInt('PHHigh')
            #?# Add checks and warnings of cutoff conflicts
            ## ~ [2c] ~ Run purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gencov = '{}.gencov'.format(bamfile)
            covstats = '{}.purge.coverage_stats.csv'.format(basefile)
            purge = '{}.purge.reassignments.tsv'.format(basefile)
            rje.checkForFiles(filelist=[gencov,covstats,purge],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            #i# The -depth setting will be increased from 200 to 2xphhigh if >100
            phcmd1 = 'purge_haplotigs hist -b {} -g {} -t {} -d {}'.format(bamfile,seqin,self.threads(),max(200,2*phhigh))
            if self.needToRemake(gencov,bamfile):
                logline = self.loggedSysCall(phcmd1,append=True)
            #!# Option to update the automatically set cutoffs
            self.printLog('#PHDEP','Low=%dX; Mid=%dX; High=%dX. (SC=%dX)' % (phlow,phmid,phhigh,scdepth))
            phcmd2 = 'purge_haplotigs cov -i {}.gencov -l {} -m {} -h {} -o {}.purge.coverage_stats.csv -j 80 -s 80'.format(bamfile,phlow,phmid,phhigh,basefile)
            if self.needToRemake(covstats,gencov):
                logline = self.loggedSysCall(phcmd2,append=True)
            else: self.printLog('#NOTE','Reusing existing %s on assumption that cutoffs have not changed' % covstats)
            phcmd3 = 'purge_haplotigs purge -g {} -c {}.purge.coverage_stats.csv -t {} -o {}.purge -a 95'.format(seqin,basefile,self.threads(),basefile)
            if self.needToRemake(purge,covstats):
                logline = self.loggedSysCall(phcmd3,append=True)
            if self.getStrLC('RunMode').startswith('purgehap'):
                self.printLog('#PURGE','purge_haplotigs run complete. Use runmode=diploidocus for additional filtering')
                return True

            ### ~ [3] ~ Run KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] KmerReads data versus Assembly KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Add --5ptrim 16,0 setting, or additional kat options on commandline?
            katfile = '{}.kat-stats.tsv'.format(basefile)
            # seq_name        median  mean    gc%     seq_length      kmers_in_seq    invalid_kmers   %_invalid       non_zero_kmers  %_non_zero      %_non_zero_corrected
            # NALACHR1.01 RevComp HiC_scaffold_38_pilon Len=123Mb; PBdep=39.8X; Nala German Shepherd Dog Chromosome 1 28      111164.29792    0.41661 122971501       122971475       54368   0.04421 122619711       99.71395        99.75805
            katcvg =  '{}.kat-counts.cvg'.format(basefile)
            rje.checkForFiles(filelist=[katfile,katcvg],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            # Establish SC kmer count using BUSCO too
            # ${PREFIX}-counts.cvg
            # >NALACHR1.01 RevComp HiC_scaffold_38_pilon Len=123Mb; PBdep=39.8X; Nala German Shepherd Dog Chromosome 1
            # 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 11 528 585 553 470 111 ...
            katcall = 'kat sect -t {} -o {}.kat {} {}'.format(self.threads(),basefile,seqin,' '.join(self.list['KmerReads']))
            if self.getBool('10xTrim'):
                trim5 = ['16'] + ['0'] * (len(self.list['KmerReads']) - 1)
                trim5 = ','.join(trim5)
                katcall = 'kat sect -t {} --5ptrim {} -o {}.kat {} {}'.format(self.threads(),trim5,basefile,seqin,' '.join(self.list['KmerReads']))
            if not self.list['KmerReads']:
                self.printLog('#KAT','Cannot use KAT kmer analysis without KmerReads data')
            if self.list['KmerReads'] and (self.force() or not (rje.exists(katfile) and rje.exists(katcvg))):
                self.printLog('#SYS',katcall)
                #i# Catching completion in case KAT hangs after running
                KAT = os.popen(katcall)
                while not KAT.readline().startswith('Total runtime'): continue
                KAT.close()
            ## ~ [3b] Assembly versus self KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            selfkat = '{}.selfkat-stats.tsv'.format(basefile)
            selfcvg =  '{}.selfkat-counts.cvg'.format(basefile)
            rje.checkForFiles(filelist=[selfkat,selfcvg],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            katcall = 'kat sect -t {} -o {}.selfkat {} {}'.format(self.threads(),basefile,seqin,seqin)
            if (self.force() or not (rje.exists(selfkat) and rje.exists(selfcvg))):
                self.printLog('#SYS',katcall)
                #i# Catching completion in case KAT hangs after running
                KAT = os.popen(katcall)
                while not KAT.readline().startswith('Total runtime'): continue
                KAT.close()

            ### ~ [4] ~ Run BBMap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # module add java/8u201-jre bbmap/38.51
            depfile = '{}.depth.tdt'.format(basefile)
            depstat = '{}.depth.stats'.format(basefile)
            rje.checkForFiles(filelist=[depfile,depstat],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            bbmap = 'samtools view -h {} | pileup.sh in=stdin out={}.depth.tdt 2>&1 | tee {}.depth.stats'.format(bamfile,basefile,basefile)
            if self.force() or not (rje.exists(depfile) and rje.exists(depstat)):
                logline = self.loggedSysCall(bbmap,append=True)
            #bblines = os.popen('samtools view -h {} | pileup.sh in=stdin out={}.depth.tdt 2>&1 | tee {}.depth.stats'.format(bamfile,basefile,basefile)).readlines()

            ### ~ [5] ~ Add VecScreen, Telomere and BUSCO Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [5a] Vecscreen coverage file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            screencov = '{}.screencov.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if self.getStrLC('ScreenDB'):
                self.printLog('#VEC','VecScreen screendb=FILE given -> running VecScreen')
                screencov = self.vecScreen()
            elif not rje.checkForFiles(filelist=[screencov],basename='',log=self.log):
                screencov = None
            vecdb = self.db('screencov')
            if screencov:
                if not vecdb: vecdb = db.addTable(screencov,mainkeys=['Query','Hit'],expect=True,name='screencov')
                vecdb.dropEntriesDirect('Query',['TOTAL'],inverse=True)
                vecdb.newKey(['Hit'])
                vecdb.setFields(['Hit','Coverage','CovPC'])
                vecdb.renameField('Coverage','ScreenCov')
                vecdb.renameField('CovPC','ScreenPerc')
            else:
                screencov = '{}.screencov.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            ## ~ [5b] Find Telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            teldb = self.findTelomeres()
            ## ~ [5c] BUSCO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            busdb = None
            busco = self.getStr('BUSCO')
            if not self.getStrLC('BUSCO'): busco = 'full_table_{}.busco.tsv'.format(self.baseFile(strip_path=True))
            if rje.checkForFiles(filelist=[busco],basename='',log=self.log,missingtext='Not found: no BUSCO counts incorporated'):
                self.printLog('#BUSCO',busco)
                fullhead = ['BuscoID','Status','Contig','Start','End','Score','Length']
                busdb = db.addTable(busco,mainkeys='auto',headers=fullhead,expect=True,name='busco')
                for rating in ['Complete','Duplicated','Fragmented','Missing']: busdb.addField(rating,evalue=0)
                for bentry in busdb.entries(): bentry[bentry['Status']] = 1
                busdb.keepFields(['#','Contig','Complete','Duplicated','Fragmented'])
                busdb.compress(['Contig'],default='sum')
                busdb.dropField('#')
                missing = 0
                for bentry in busdb.entries():
                    if bentry['Contig'] not in seqdict: busdb.dropEntry(bentry); missing += 1
                self.printLog('#BUSCO','{} BUSCO contigs dropped: not found in {}'.format(rje.iStr(missing),seqin))

            ### ~ [6] ~ Compile Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            joinlist = []
            # >> join:list of (Table,Field[,Fieldlist]) tuples, where Table is a table name and Field is a Field name or
            #     formula to be used for the join. Fieldlist is an optional list of Fields from that Table to include in the
            #     new table. If Field does not exist, it will be added. (Field may be a Formula.)
            # # ~ [6a] Summarise Check Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # depfile = '{}.depth.tdt'.format(basefile)
            # covstats = '{}.purge.coverage_stats.csv'.format(basefile)
            # purge = '{}.purge.reassignments.tsv'.format(basefile)
            # katfile = '{}.kat-stats.tsv'.format(basefile)
            # selfkat = '{}.selfkat-stats.tsv'.format(basefile)
            rje.checkForFiles(filelist=[depfile,covstats,purge,katfile,selfkat,busco,screencov],basename='',log=self.log,missingtext='Not found: excluded from classification')

            ## ~ [6b] BBMap Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.depth.tdt <==
            # #ID     Avg_fold        Length  Ref_GC  Covered_percent Covered_bases   Plus_reads      Minus_reads     Read_GC Median_fold     Std_Dev
            # Scaff10x_0  29.6599 14742825        0.0000  99.9948 14742058        28191   27814   0.4111  29      10.10
            depdb = db.addTable(depfile,mainkeys=['#ID'],expect=True,name='depth')
            depdb.renameField('#ID','SeqName')
            depdb.renameField('Length','SeqLen')
            depdb.setFields(['SeqName','SeqLen','Median_fold','Avg_fold','Covered_percent','Covered_bases','Plus_reads','Minus_reads','Read_GC'])
            joinlist.append((depdb,'SeqName'))

            ## ~ [6c] PurgeHaplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# NOTE: Had a problem with a sequence missing from *.coverage_stats.csv
            #  - will need some default values and warnings
            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.purge.coverage_stats.csv <==
            # #contig,contig_reassign,bases_hap_dip,bases_low_high,bases_all,perc_low_coverage,perc_hap_coverage,perc_dip_coverage,perc_high_coverage
            # Scaff10x_0,,14719523,37459,14756982,0.004,9.823,89.923,0.250
            phcovdb = db.addTable(covstats,mainkeys=['#contig'],expect=True,name='phcov')
            phcovdb.renameField('#contig','SeqName')
            for cov in ['Low','Hap','Dip','High']:
                phcovdb.renameField('perc_{}_coverage'.format(cov.lower()), '{}Perc'.format(cov))
            phcovdb.setFields(['SeqName','LowPerc','HapPerc','DipPerc','HighPerc'])
            joinlist.append((phcovdb,'SeqName'))

            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.purge.reassignments.tsv <==
            # #reassigned_contig      top_hit_contig  second_hit_contig       best_match_coverage     max_match_coverage      reassignment
            # Scaff10x_1000   Scaff10x_562    Scaff10x_676    98.37   918.99  REPEAT
            purgedb = db.addTable(purge,mainkeys=['#reassigned_contig'],expect=True,name='purge')
            purgedb.renameField('#reassigned_contig','SeqName')
            purgedb.renameField('top_hit_contig','TopHit')
            purgedb.renameField('second_hit_contig','SecHit')
            purgedb.renameField('best_match_coverage','TopHitCov')
            purgedb.renameField('max_match_coverage','MaxHitCov')
            purgedb.renameField('reassignment','PurgeHap')
            purgedb.index('TopHit')
            purgedb.index('SecHit')
            purgedb.addFields(['TopNum','SecNum'],evalue=0)
            for entry in purgedb.entries():
                if entry['SeqName'] in purgedb.index('TopHit'): entry['TopNum'] = len(purgedb.index('TopHit')[entry['SeqName']])
                if entry['SeqName'] in purgedb.index('SecHit'): entry['SecNum'] = len(purgedb.index('SecHit')[entry['SeqName']])
            joinlist.append((purgedb,'SeqName'))

            ## ~ [6d] KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.kat-stats.tsv <==
            # seq_name        median  mean    gc%     seq_length      kmers_in_seq    invalid_kmers   %_invalid       non_zero_kmers  %_non_zero      %_non_zero_corrected
            # Scaff10x_0      1       1854.52173      0.41236 14757069        14757043        0       0.00000 14757043        100.00000       100.00000
            selfdb = db.addTable(selfkat,mainkeys=['seq_name'],expect=True,name='selfkat')
            selfdb.renameField('seq_name','SeqName')
            selfdb.renameField('median','SelfMedK')
            selfdb.renameField('mean','SelfAvgK')
            selfdb.setFields(['SeqName','SelfMedK','SelfAvgK'])
            joinlist.append((selfdb,'SeqName'))
            for entry in selfdb.entries(): entry['SeqName'] = entry['SeqName'].split()[0]
            selfdb.remakeKeys()

            #!# Add capacity for running without KAT -> MedK=-1
            katdb = db.addTable(katfile,mainkeys=['seq_name'],expect=False,name='kat')
            if katdb:
                katdb.renameField('seq_name','SeqName')
                katdb.renameField('median','MedK')
                katdb.renameField('mean','AvgK')
                katdb.renameField('gc%','SeqGC')
                katdb.renameField('%_non_zero_corrected','KPerc')
                katdb.setFields(['SeqName','MedK','AvgK','SeqGC','KPerc'])
                for entry in katdb.entries(): entry['SeqName'] = entry['SeqName'].split()[0]
                katdb.remakeKeys()
            else:
                katdb = db.addEmptyTable('kat',['SeqName','MedK','AvgK','SeqGC','KPerc'],['SeqName'])
                for entry in selfdb.entries(): katdb.addEntry({'SeqName':entry['SeqName'],'MedK':-1,'AvgK':-1,'SeqGC':-1,'KPerc':-1})
            joinlist.append((katdb,'SeqName'))

            ## ~ [6e] Other ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if vecdb:
                vecdb.renameField('Hit','SeqName')
                joinlist.append((vecdb,'SeqName'))
            if teldb:
                teldb.renameField('Name','SeqName')
                teldb.dropField('SeqLen')
                joinlist.append((teldb,'SeqName'))
            if busdb:
                busdb.renameField('Contig','SeqName')
                joinlist.append((busdb,'SeqName'))

            ## ~ [6x] Join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dipdb = db.joinTables(name='diploidocus',join=joinlist,newkey=['SeqName'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True,warnings=True)
            dipdb.remakeKeys()
            mingap = max(1,self.getInt('MinGap'))
            self.printLog('#GAPS','Identifying all runs of %d+ Ns as gaps' % mingap)
            if self.getBool('ZeroAdjust'):
                if self.getBool('IncludeGaps'):
                    self.printLog('#GAPS','Gap regions will be included in zero-coverage purge_haplotigs percentage adjustments (includegaps=T)')
                else:
                    self.printLog('#GAPS','Gap regions will be excluded from zero-coverage purge_haplotigs percentage adjustments (includegaps=F)')
            dipdb.addFields(['N_bases','Gap_bases'],evalue=-1)
            dipdb.addField('SeqDesc')
            ex = 0.0; etot = dipdb.entryNum(); gapwarn = 0; gaptot = 0
            for entry in dipdb.entries():
                self.progLog('\r#GAPS','Adding sequence gap counts and descriptions: %.1f%%' % (ex/etot)); ex += 100.0
                gapn = 0
                try:
                    seqname, seq = seqlist.getSeq(seqdict[entry['SeqName']])
                    seq = seq.upper()
                    entry['N_bases'] = seq.count('N')
                    seqdata = string.split(seqname,maxsplit=1)
                    if len(seqdata) > 1: entry['SeqDesc'] = seqdata[1]
                    gapn = len(''.join(re.findall('N{%d,}' % mingap,seq)))
                    gaptot += gapn
                    if gapn > (len(seq)/2.0): gapwarn += 1
                except:
                    self.warnLog('Could not extract description for %s!' % entry['SeqName'])
                    entry['SeqDesc'] = 'ERROR!'
                    open('seqdict.tmp','w').write('{}\n'.format(seqdict).replace(',','\n'))
                entry['Gap_bases'] = gapn
            self.printLog('\r#GAPS','Added %s bp sequence gaps (%d+ Ns) and descriptions.' % (rje.iStr(gaptot),mingap))
            if gapwarn: self.warnLog('{} sequences have >50% gaps. Check use of minmedian=X'.format(rje.iStr(gapwarn)))
            self.printLog('\r#FIELDS',', '.join(dipdb.fields()))
            #i# Tidy up join
            dipdb.fillBlanks(blank='False',fields=['Tel5','Tel3'],fillempty=True,prog=True,log=True)
            dipdb.fillBlanks(blank=-1,fields=['Trim5','Trim3'],fillempty=True,prog=True,log=True)
            dipdb.fillBlanks(blank=0,fields=['ScreenCov','Complete','Duplicated','Fragmented'],fillempty=True,prog=True,log=True)
            dipdb.fillBlanks(blank=0.0,fields=['TelPerc','ScreenPerc'],fillempty=True,prog=True,log=True)
            for entry in dipdb.entries():
                if entry['TopHitCov'] == '-': entry['TopHitCov'] = 0.0
                if entry['MaxHitCov'] == '-': entry['MaxHitCov'] = 0.0
            #!# Reorder dipdb fields
            fields = {'str':['SeqName', 'TopHit', 'SecHit', 'PurgeHap'],
                    'int': ['SeqLen', 'Median_fold', 'Covered_bases', 'Plus_reads', 'Minus_reads','TopNum','SecNum', 'SelfMedK', 'MedK', 'ScreenCov', 'Trim5', 'Trim3', 'Complete', 'Duplicated', 'Fragmented','Gap_bases','N_bases'],
                    'num': ['Avg_fold', 'Covered_percent', 'Read_GC', 'LowPerc', 'HapPerc', 'DipPerc', 'HighPerc', 'TopHitCov', 'MaxHitCov', 'SelfAvgK', 'AvgK', 'SeqGC', 'KPerc', 'ScreenPerc', 'TelPerc'],
                    'bool': ['Tel5', 'Tel3']}
            reformat = {}
            for ftype, tfields  in fields.items():
                for field in tfields: reformat[field] = ftype
            dipdb.dataFormat(reformat)

            ## ~ [6z] ~ Zero Adjustment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('ZeroAdjust'):
                ex = 0.0; etot = dipdb.entryNum(); adjx = 0
                for d in dipdb.entries():
                    self.progLog('\r#ZERO','Adjusting purge_hapotigs percentages for zero coverage: %.1f%%' % (ex/etot)); ex += 100.0
                    zerobp = d['SeqLen'] - d['Covered_bases']
                    if not self.getBool('IncludeGaps'): zerobp = max(0, zerobp - d['Gap_bases'])
                    if zerobp < 0: self.warnLog('Covered_bases < SeqLen for {}!'.format(d['SeqName'])); continue
                    if not zerobp: continue
                    if d['Covered_bases'] < 1:
                        d['LowPerc'] = 100.0
                        d['HapPerc'] = d['DipPerc'] = d['HighPerc'] = 0.0
                        adjx += 1
                        continue
                    norm = d['Covered_bases'] / 100.0
                    bpbins = [d['LowPerc'], d['HapPerc'], d['DipPerc'], d['HighPerc']]
                    try:
                        bpbins = [zerobp + (norm * d['LowPerc']), norm * d['HapPerc'], norm * d['DipPerc'], norm * d['HighPerc']]
                    except:
                        self.debug(d)
                        self.errorLog(self.wisdom())
                    newsum = sum(bpbins)
                    norm = 100.0 / newsum
                    d['LowPerc'] = norm * bpbins[0]
                    d['HapPerc'] = norm * bpbins[1]
                    d['DipPerc'] = norm * bpbins[2]
                    d['HighPerc'] = norm * bpbins[3]
                    adjx += 1
                self.printLog('\r#ZERO','Adjusted purge_hapotigs percentages for zero coverage in %s of %s sequences' % (rje.iStr(adjx),rje.iStr(etot)))

            ### ~ [7] ~ Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DIPLOIDOCUS SCAFFOLD RATING',line='=')
            dipdb.addField('Class', evalue='-')
            dipdb.addField('TopClass', evalue='-')
            dipdb.addField('SecClass', evalue='-')
            dipdb.addField('Rating', evalue='UNRATED')
            dipdb.addField('TopRating', evalue='-')
            dipdb.addField('SecRating', evalue='-')
            ## ~ [7a] ~ Cutoffs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            minmedian = self.getInt('MinMedian')
            self.printLog('#MINDEP','Low coverage min. median depth filter: {:d}X'.format(minmedian))
            mindipcov = 20  # Minimim DipPerc coverage to assign a "KEEP" rating
            self.printLog('#DIPCOV','Min. DipPerc coverage to assign as PRIMARY: {:d}%'.format(mindipcov))
            artefactcov = 80   # Min High/Low coverage to assign as artefacts
            self.printLog('#LOWCOV','Min. High/Low coverage to assign as artefacts: {:d}%'.format(artefactcov))
            hpurgecov = 80  # Min TopHitCov to assign HPURGE
            self.printLog('#HPURGE','Min. TopHitCov for HPURGE rating: {:d}%'.format(hpurgecov))
            rephitcov = 250  # Min MaxHitCov for REPEAT
            self.printLog('#REPEAT','Min. MaxHitCov for REPEAT rating: {:d}%'.format(rephitcov))
            pureperc = 80  # Depth class percentage to classify as pure
            self.printLog('#PURE','Depth class percentage to classify as pure: {:d}%'.format(pureperc))
            covwarning = 95
            self.printLog('#COV','Read coverage threshold for warnings: {:d}%'.format(covwarning))
            purge = 'normal'
            ## ~ [7b] ~ Pre-Rating Classifcation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in dipdb.entries():
                #i# Set up classification list to join
                eclass = ['WEAK','Dip']
                #i# Establish primary depth class
                for dep in ['Dip','Hap','High','Low']:
                    if entry['{}Perc'.format(dep)] >= pureperc: eclass[0] = 'PURE'
                    elif entry['{}Perc'.format(dep)] >= 50: eclass[0] = 'GOOD'
                if entry['Median_fold'] < minmedian: eclass[0] = 'LOWX'
                for dep in ['Hap','High','Low']:
                    if entry['{}Perc'.format(dep)] > entry['{}Perc'.format(eclass[1])]: eclass[1] = dep
                if eclass[1] == 'High': eclass[1] = 'EXS'
                else: eclass[1] = eclass[1].upper()
                #i# TopHit, SecHit and HitCov - Homology/Repeat status
                if not entry['TopHit'] or entry['TopHit'] == '-': eclass.append('UNIQ')
                elif entry['TopHitCov'] < 50: eclass.append('PART')
                elif not entry['SecHit'] or entry['SecHit'] == '-': eclass.append('HAPL')
                elif entry['MaxHitCov'] < rephitcov: eclass.append('HOMO')
                else: eclass.append('REPT')
                #i# TopHitNum and SecHitNum (v bad sequences will not be top hits)
                if entry['TopNum'] > 0: eclass.append('TOP')
                elif entry['SecNum'] > 0: eclass.append('SEC')
                else: eclass.append('NON')
                #i# KAT based self-assessment
                if entry['SelfMedK'] == 1: eclass.append('PRI')
                elif entry['SelfMedK'] == 2: eclass.append('ALT')
                else: eclass.append('REP')
                #i# BUSCO Genes (can help decide about risk)
                if entry['Complete'] > 0 and entry['Complete'] > entry['Duplicated']: eclass.append('COMP')
                elif entry['Duplicated'] > 0: eclass.append('DUPL')
                elif entry['Fragmented'] > 0: eclass.append('FRAG')
                else: eclass.append('NONE')
                #i# Extras
                if entry['TelPerc'] > 0: eclass.append('+TEL')
                if entry['ScreenCov'] > 0: eclass.append('+VEC')
                #i# Join
                entry['Class'] = '|'.join(eclass)
            dipdb.indexReport('Class')

            ## ~ [7b] ~ Default Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in dipdb.entries():
                tentry = dipdb.data(entry['TopHit'],expect=False)
                sentry = dipdb.data(entry['SecHit'],expect=False)
                if tentry: entry['TopClass'] = tentry['Class']
                if sentry: entry['SecClass'] = sentry['Class']
                # PURITY | DEPTH | HOM | TOP | MEDK | BUSCO
                cdata = entry['Class'].split('|')
                rep = cdata[4] == "REP"

                ### --- INITIAL LOW QUALITY FILTER --- ###
                # * `CONTAMINATION` = Scaffolds with 50%+ identified contamination
                if entry['ScreenPerc'] >= 50: entry['Rating'] = 'CONTAMINATION'
                #  - `Median_fold` < 3
                elif entry['Class'].startswith('LOWX'): entry['Rating'] = 'LOWCOV'
                #  - `[GOOD/PURE]|LOW` and `MedK`<1
                elif entry['Class'].startswith('PURE|LOW') and entry['MedK'] == 0: entry['Rating'] = 'LOWCOV'
                elif entry['Class'].startswith('GOOD|LOW') and entry['MedK'] == 0: entry['Rating'] = 'LOWCOV'
                # - `LOW` & `REP`
                elif cdata[1] == 'LOW' and rep: entry['Rating'] = 'LOWCOV'

                ### --- KEEP --- ###
                # * `QUALITY` = Highest quality scaffolds: pure duploid, complete BUSCOs, no Duplicated BUSCOs
                #  - `PURE|DIP|UNIQ` & `PRI|COMP` & `Duplicated` = 0
                elif entry['Class'].startswith('PURE|DIP|UNIQ') and "PRI|COMP" in entry['Class'] and entry['Duplicated'] == 0: entry['Rating'] = 'QUALITY'
                #* `FINAL` = As Quality but `PRI|FRAG` and `PRI|NONE` also allowed
                elif entry['Class'].startswith('PURE|DIP|UNIQ') and cdata[4] == 'PRI' and entry['Duplicated'] == 0: entry['Rating'] = 'FINAL'
                #* `CORE` = Predominantly Diploid with <50% covered by TopHit and SelfMedK=1
                elif cdata[0] in ("PURE","GOOD") and cdata[1] == "DIP" and cdata[2] in ("UNIQ","PART") and cdata[4] == "PRI": entry['Rating'] = 'CORE'
                #* `COREHAP` = Predominantly haploid but <50% covered by TopHit and 1+ Complete BUSCOs
                #  - `[PURE/GOOD]|HAP` & `UNIQ/PART` & `Complete` > 0
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and entry['Complete'] > 0: entry['Rating'] = 'COREHAP'
                #* `PRIMARY` = Putative primary scaffold but with possible alternative scaffolds still in assembly and/or low quality regions
                elif cdata[1] == "DIP" and cdata[2] in ("UNIQ","PART","HAPL","HOMO") and not rep: entry['Rating'] = 'PRIMARY'

                #* `PRIRPT` = Putative primary scaffold but >50% repeated
                #  - DIP with UNIQ/PART and REP
                elif cdata[1] == "DIP" and cdata[2] in ("UNIQ","PART") and rep: entry['Rating'] = 'PRIRPT'
                #   - `[PURE/GOOD]|DIP|[HAPL/HOMO]` & `REP`
                elif cdata[0] in ("PURE","GOOD") and cdata[1] == "DIP" and rep: entry['Rating'] = 'PRIRPT'

                #* `COLLAPSED` = High coverage scaffolds representing putative collapsed repeats. Check for other contamination.
                #  - `EXS` & `UNIQ`
                elif cdata[1] == "EXS" and cdata[2] in ("UNIQ","PART"): entry['Rating'] = 'COLLAPSED'
                # - `[GOOD/PURE]|EXS` & `TOP`
                elif cdata[0] in ("PURE","GOOD") and cdata[1] == "EXS" and cdata[3] == "TOP": entry['Rating'] = 'COLLAPSED'
                #  - `WEAK|EXS` & (`HapPerc`+`LowPerc`)<=`DipPerc`
                elif cdata[1]  == "EXS" and cdata[0] == "WEAK" and (entry['HapPerc'] + entry['LowPerc']) <= entry['DipPerc']: entry['Rating'] = 'COLLAPSED'

                # * `REPEAT` = Predominantly Diploid scaffolds that have major signs of redundancy, probably due to presence of alternative contigs
                #   - `DIP|REPT`
                elif "|DIP|REPT|" in entry['Class']: entry['Rating'] = 'REPEAT'
                #- `WEAK|DIP|[HAPL/HOMO]` & `REP`
                elif cdata[0] == 'WEAK' and cdata[1] == 'DIP' and cdata[2] in ("HAPL","HOMO") and rep: entry['Rating'] = 'REPEAT'

                #* `HAPLOID` = Predominantly haploid coverage but enough unique sequence to keep for now? (Option to quarantine?) Might be very heterozygous Alternatice haplotigs
                #  - `HAP` `UNIQ/PART` & `PRI`
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and cdata[4] == "PRI": entry['Rating'] = 'HAPLOID'
                #* `HAPLOTIG` = Predominantly haploid coverage but enough unique sequence to keep for now and/or a TOP hit? Probable Alternative haplotig
                #  - `HAP` `UNIQ/PART` & `ALT`
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and cdata[4] == "ALT": entry['Rating'] = 'HAPLOTIG'
                elif cdata[1] == "HAP" and cdata[3] == "TOP" and not rep: entry['Rating'] = 'HAPLOTIG'
                #* `HAPRPT` = Low quality scaffold that is probably most repeats, but not bad enough to dump outright
                #  - `HAP` & `UNIQ/PART` & `REP`
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and rep: entry['Rating'] = 'HAPRPT'
                elif cdata[1] == "HAP" and cdata[3] == "TOP" and rep: entry['Rating'] = 'HAPRPT'
                ## Other HAPLOTIGs
                #  - `HAP/LOW` & `HAPL/HOMO/REPT` & `TopClass`-`DIP`
                elif cdata[1] in ("LOW","HAP") and cdata[2] not in ("UNIQ","PART") and '|DIP|' in entry['TopClass'] and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPLOTIG'
                elif cdata[0] in ('PURE','GOOD') and '|HAP|HOMO|' in entry['Class'] and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPLOTIG'
                # Assess remaining HAP scaffolds
                # - `HAP` and dipdb$TopHitCov < hpurgecov = HAPREPEAT
                elif cdata[1] == "HAP" and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPLOTIG'

                #  - `WEAK|EXS` & (`HapPerc`+`LowPerc`)>`DipPerc`
                elif cdata[1] == "EXS" and cdata[0] == "WEAK" and (entry['HapPerc'] + entry['LowPerc']) > entry['DipPerc'] and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPRPT'
                # - `[GOOD/PURE]|EXS[HAPL/HOMO/REPT]` & not `TOP`
                elif cdata[1] == "EXS" and cdata[0] != "WEAK" and cdata[2] not in ("UNIQ","PART") and cdata[3] != "TOP" and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'REPEAT'

                ### --- DUMP --- ###
                #* `RPURGE` = Messy scaffolds that are largely repeats and are sufficiently redundant/low quality to purge
                #  - `WEAK|EXS` & (`HapPerc`+`LowPerc`)>`DipPerc`
                elif cdata[1] == "EXS" and cdata[0] == "WEAK" and (entry['HapPerc'] + entry['LowPerc']) > entry['DipPerc'] and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'RPURGE'
                # - `[GOOD/PURE]|EXS[HAPL/HOMO/REPT]` & not `TOP`
                elif cdata[1] == "EXS" and cdata[0] != "WEAK" and cdata[2] not in ("UNIQ","PART") and cdata[3] != "TOP" and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'RPURGE'

                # - `LOW` & `HAPL/HOMO/REPT` & NOT `TOP`
                elif cdata[1] == 'LOW' and cdata[2] not in ['UNIQ','PART'] and cdata[3] != 'TOP': entry['Rating'] = 'LOWCOV'

                #* `HPURGE` = Clear candidate haplotig to purge
                #  - `HAP/LOW` & `HAPL/HOMO/REPT` & `TopClass`-`DIP`
                elif cdata[1] in ("LOW","HAP") and cdata[2] not in ("UNIQ","PART") and '|DIP|' in entry['TopClass'] and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                elif cdata[0] in ('PURE','GOOD') and '|HAP|HOMO|' in entry['Class'] and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                # - `HAP` and dipdb$TopHitCov >= hpurgecov = HPURGE
                elif cdata[1] == "HAP" and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'

                ### ---- LEFTOVERS ---- ###
                # * `LOWQUAL` = Unconvincing scaffolds that do not fall into a clear class but are not bad enough to dump outright
                #   - `LOW` & otherwise escaping `LOWCOV` filtering
                elif cdata[1] == "LOW": entry['Rating'] = 'LOWQUAL'

                ### Change HAPLOTIG to HAPRPT
                if entry['Rating'] == "HAPLOTIG" and rep: entry['Rating'] = 'HAPRPT'

            ## ~ [7b] ~ Basic Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if purgemode == 'simple':
                for entry in dipdb.entries():
                    # PURITY | DEPTH | HOM | TOP | MEDK | BUSCO
                    cdata = entry['Class'].split('|')
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    #i# Low quality artefacts (e.g. less than 50% of the scaffold covered by short-read kmers)
                    if entry['MedK'] == 0: entry['Rating'] = 'LOWQUAL'; continue
                    #i# KEEP set are those with DIP, EXS, UNIQ, PART, TOP or COMP = REPT/REP = REPEAT
                    if cdata[1] in ['DIP','EXS'] or cdata[2] in ['UNIQ','PART'] or cdata[3] == 'TOP' or entry['Complete'] > 0:
                        # PRIMARY
                        if cdata[1] in ['DIP'] or cdata[2] in ['UNIQ'] or cdata[5] == 'COMP':  #?# entry['Complete'] > 0:
                            if cdata[2] == 'REPT' or cdata[4] == 'REP': entry['Rating'] = 'PRIRPT'
                            else: entry['Rating'] = 'PRIMARY'
                        # COLLAPSED
                        elif cdata[1] == 'EXS': entry['Rating'] = 'COLLAPSED'
                        else:
                            if cdata[2] == 'REPT' or cdata[4] == 'REP': entry['Rating'] = 'REPEAT'
                            else: entry['Rating'] = 'HAPLOID'
                    #i# QUARANTINE any HAP
                    elif cdata[1] == 'HAP':
                        if cdata[2] == 'REPT' or cdata[4] == 'REP': entry['Rating'] = 'RPURGE'
                        else: entry['Rating'] = 'HPURGE'
                    #i# JUNK anything else LOW
                    elif cdata[1] == 'LOW': entry['Rating'] = 'LOWCOV'
                    else: raise ValueError('No rating: %s' % entry)
                    if entry['SeqName'] == 'cam006_MBG479__MBG479C13.006':
                        self.deBug('Rating: %s' % entry)

            ## ~ [7b] ~ Basic Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode == 'dev':
                for entry in dipdb.entries():
                    tentry = dipdb.data(entry['TopHit'],expect=False)
                    sentry = dipdb.data(entry['SecHit'],expect=False)
                    if tentry: entry['TopClass'] = tentry['Class']
                    if sentry: entry['SecClass'] = sentry['Class']
                    #i# INITIAL FILTERS #i#
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    #i# Any scaffolds with <80% at diploid read depth were identified by PurgeHaplotigs for reassignment.
                    # Scaffolds with 80%+ bases in the low/haploid coverage bins and 95%+ of their length mapped by
                    # PurgeHaplotigs onto another scaffold were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                    elif (entry['HapPerc']) >= 50 and entry['TopHitCov'] >= hpurgecov and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    # Any other scaffolds exceeding the mindipcov DipPerc coverage assigned to KEEP for further assignment
                    elif (entry['DipPerc']) >= mindipcov and entry['SelfMedK'] == 1: entry['Rating'] = 'PRIMARY'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'REPEAT'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov: entry['Rating'] = 'KEEP'  # 'PRIMARY' ?
                    # Any other scaffolds with 80%+ low coverage bases were filtered as Low Coverage.
                    elif (entry['LowPerc']) >= artefactcov: entry['Rating'] = 'LOWCOV'
                    # Any other scaffolds with 80%+ high coverage bases are flagged as COLLAPSED.
                    elif (entry['HighPerc']) >= artefactcov: entry['Rating'] = 'COLLAPSED'
                    # Scaffolds that are mostly Haploid Depth are marked as HAPLOTIGS if mostly present 2+ times, or HAPLOID if present once
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 1: entry['Rating'] = 'HAPLOID'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif entry['HapPerc'] < 50 and (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'
                    # Scaffolds with at least 80% hitting other scaffolds are HAPLOTIG or REPEAT
                    elif entry['TopHitCov'] >= hpurgecov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'HAPREPEAT'
                    elif entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HAPLOTIG'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'

            ## ~ [7x] ~ Basic Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode == 'crude':
                for entry in dipdb.entries():
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    #i# Any scaffolds with <80% at diploid read depth were identified by PurgeHaplotigs for reassignment.
                    # Scaffolds with 80%+ bases in the low/haploid coverage bins and 95%+ of their length mapped by
                    # PurgeHaplotigs onto another scaffold were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                    elif (entry['HapPerc']) >= 50 and entry['TopHitCov'] >= hpurgecov and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    # Any other scaffolds exceeding the mindipcov DipPerc coverage assigned to KEEP for further assignment
                    elif (entry['DipPerc']) >= mindipcov and entry['SelfMedK'] == 1: entry['Rating'] = 'PRIMARY'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'REPEAT'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov: entry['Rating'] = 'KEEP'  # 'PRIMARY' ?
                    # Any other scaffolds with 80%+ low coverage bases were filtered as Low Coverage.
                    elif (entry['LowPerc']) >= artefactcov: entry['Rating'] = 'LOWCOV'
                    # Any other scaffolds with 80%+ high coverage bases are flagged as COLLAPSED.
                    elif (entry['HighPerc']) >= artefactcov: entry['Rating'] = 'COLLAPSED'
                    # Scaffolds that are mostly Haploid Depth are marked as HAPLOTIGS if mostly present 2+ times, or HAPLOID if present once
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 1: entry['Rating'] = 'HAPLOID'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif entry['HapPerc'] < 50 and (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'
                    # Scaffolds with at least 80% hitting other scaffolds are HAPLOTIG or REPEAT
                    elif entry['TopHitCov'] >= hpurgecov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'HAPREPEAT'
                    elif entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HAPLOTIG'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'

            ## ~ [7x] ~ Add Nala Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode == 'nala':
                self.errorLog('#PURGE','purgemode=nala not yet implemented!')
                # Low-coverage filter
                #
                # The TOW5157A1 library PacBio subreads (12.5M subreads; 108Gb) were mapped onto the assembly using Minimap2 v2.16
                # (-ax map-pb --secondary=no) [REF]. Initial read depth analysis was performed with BBMap v38.51 [REF] pileup.sh.
                # Any scaffolds with median coverage less than 3 (e.g. less than 50% of the scaffold covered by at least 3 reads)
                # were filtered as Low Coverage.
                #
                # Purge Haplotigs analysis - round 1
                #
                # Subreads were re-mapped on to the remaining 837 scaffolds and processed with PurgeHaplotigs v20190612 [REF]
                # (implementing Perl v5.28.0, BEDTools v2.27.1 [REF], R v3.5.3, and SAMTools v1.9 [REF]). Based on the
                # PurgeHaplotigs depth histogram, low-, mid- and high-depth thresholds were set to 5X, 30X and 80X. Any scaffolds
                # with <80% at diploid read depth were identified by PurgeHaplotigs for reassignment. Scaffolds with 80%+ bases in
                # the low/haploid coverage bins and 95%+ of their length mapped by PurgeHaplotigs onto another scaffold were
                # filtered as haplotigs or assembly artefacts. Any other scaffolds with 80%+ low coverage bases were filtered as
                # Low Coverage.
                #
                # Purge Haplotigs analysis - round 2
                #
                # Subreads were re-mapped on to the remaining 558 scaffolds for a second round of slightly more stringent
                # PurgeHaplotigs analysis. No more scaffolds with 80%+ low coverage bases were identified. Any scaffold with 80%+
                # bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts. Scaffolds with 20%+
                # diploid coverage were marked as retention as probable diploids. Scaffolds with <20% diploid coverage and 50%+
                # high coverage were marked as probable collapsed repeats. A single remaining Scaffold marked as JUNK by
                # PurgeHaplotigs (over 80% low/high coverage) was also filtered as a probable artefact.
                #
                # Purge Haplotigs analysis - round 3
                #
                # Subreads were re-mapped on to the remaining 430 scaffolds for a third round of PurgeHaplotigs analysis. No
                # further scaffolds were identified for filtering.
                #
                # Results
                # Of the 1,057 Pilon-polished scaffolds, 627 (59.3%) totalling 20.1 Mb (0.8%) of the assembly were removed by the
                # final scaffold cleanup. 220 scaffolds were removed in the initial Low Coverage filter. Based on the
                # PurgeHaplotigs depth histogram, a diploid ("1X") read depth of approx. 40X was identified and used to infer a
                # ("0.5X") haploid depth of 20X and a mid-depth point for PurgeHaplotigs of 30X. Based on the histogram, the
                # low- and high-depth PurgeHaplotigs were set to 5X and 80X, respectively. Based on these data, a further 11
                # Scaffolds were filtered for low coverage and 268 were filtered as haplotigs or assembly artefacts. Following the
                # second round of PurgeHaplotigs, a further 128 Scaffolds were filtered as haplotigs or assembly artefacts. None of
                # the remaining 430 scaffolds were identified for filtering following additional PurgeHaplotigs analysis.

            ## ~ [7x] ~ Add Nala Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode not in ['','none','default','diploidocus','complex']:
                self.errorLog('#PURGE','purgemode={} not recognised! Using purgemode=complex'.format(purgemode))

            ## ~ [7c] ~ Add Top and SecHit Ratings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in dipdb.entries():
                tentry = dipdb.data(entry['TopHit'],expect=False)
                sentry = dipdb.data(entry['SecHit'],expect=False)
                if tentry: entry['TopRating'] = tentry['Rating']
                if sentry: entry['SecRating'] = sentry['Rating']

            ## ~ [7d] ~ Add Set for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dipdb.addField('Set',evalue='-')
            keep = ['QUALITY','FINAL','CORE','COREHAP','PRIMARY','PRIRPT','HAPLOID','HAPLOTIG']
            repeat = ['COLLAPSED', 'REPEAT', 'HAPRPT']
            quarantine = ['HPURGE', 'RPURGE']
            junk = ['LOWCOV', 'LOWQUAL', 'CONTAMINATION']

            for entry in dipdb.entries():
                if entry['Rating'] in keep: entry['Set'] = 'keep'
                elif entry['Rating'] in repeat: entry['Set'] = 'repeat'
                elif entry['Rating'] in quarantine: entry['Set'] = 'quarantine'
                elif entry['Rating'] in junk: entry['Set'] = 'junk'
                else: self.warnLog('Unknown Rating for %s: %s' % (entry['SeqName'],entry['Rating']))

            ## ~ [7d] ~ Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Will need to rationalise these warnings: so many!
            dipdb.index('TopHit')
            # * All `CONTAMINATION` scaffolds
            if 'CONTAMINATION' in dipdb.index('Rating'):
                self.warnLog('{} scaffolds had 50%+ contamination: {}'.format(rje.iLen(dipdb.index('Rating')['CONTAMINATION']),', '.join(dipdb.index('Rating')['CONTAMINATION'])))
            for entry in dipdb.entries():
                if entry['Set'] != 'junk':
                    # * `MedK` < 1 = Scaffolds where most of the scaffold is not in Illumina data.
                    if entry['MedK'] == 0:
                        self.warnLog('MedK warning: less than half of {} {} has kmerread support.'.format(entry['Rating'],entry['SeqName']))
                    # * `ScreenCov > 0` = Scaffolds with possible vector contamination
                    if entry['ScreenPerc'] > 0:
                        self.warnLog('{} {} has {}% contamination detected'.format(entry['Rating'],entry['SeqName'],entry['ScreenPerc']))
                    # * Covered_percent < 95
                    if entry['Covered_percent'] < covwarning:
                        self.warnLog('{} {} only has {}% read coverage'.format(entry['Rating'],entry['SeqName'],entry['Covered_percent']))
                # * Any `HAPLOID` scaffolds with `Duplicated`
                if entry['Rating'] in ['COREHAP','HAPLOID','HAPLOTIG'] and entry['Duplicated'] > 0:
                    self.warnLog('{} {} has {} Duplicated BUSCO: evidence for redundant haplotig?'.format(entry['Rating'],entry['SeqName'],entry['Duplicated']))
                if entry['Set'] in ['junk','quarantine']:
                    # * Any purged scaffolds with `TOP`
                    if entry['TopNum'] > 0:
                        self.warnLog('{} ({}) {} is the Best PurgeHap hit for {} Scaffolds: {}'.format(entry['Rating'],entry['Set'],entry['SeqName'],entry['TopNum'],', '.join(dipdb.index('TopHit')[entry['SeqName']])))
                #?# Think about adding warning for:
                # * Haploid scaffolds that could be Alternative haplotigs, or they could be Haploid alternative regions,
                # e.g. Sex chromosomes. QUESTION: Should be extract one sex chromosome as an alternative contig?!

            ## ~ [7e] ~ Save Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.headLog('Rating Summary',line='-')
            dipdb.indexReport('Rating')
            self.headLog('Sequence Sets',line='-')
            dipdb.indexReport('Set')
            dipdb.saveToFile(sfdict={'LowPerc':4, 'HapPerc':4, 'DipPerc':4, 'HighPerc':4})

            for seqset in ['diploidocus','core','quarantine','junk']:
                rje.backup(self,'%s.%s.fasta' % (basefile,seqset),appendable=False)
                open('%s.%s.fasta' % (basefile,seqset),'w')
            seqx = {'diploidocus':0,'core':0,'quarantine':0,'junk':0}
            for rating in keep:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.diploidocus.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    open('{}.core.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['diploidocus'] += 1; seqx['core'] += 1
            for rating in repeat:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.diploidocus.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['diploidocus'] += 1
            for rating in quarantine:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.quarantine.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['quarantine'] += 1
            for rating in junk:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.junk.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['junk'] += 1
            for seqset in ['diploidocus','core','quarantine','junk']:
                self.printLog('#SEQ','%s sequences saved to %s.%s.fasta' % (rje.iStr(seqx[seqset]),basefile,seqset))

            if seqlist.getBool('Summarise'):
                for seqset in ['diploidocus','core','quarantine','junk']:
                    if not seqx[seqset]: continue
                    setfas = '%s.%s.fasta' % (basefile,seqset)
                    seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % setfas,'autofilter=F']
                    rje_seqlist.SeqList(self.log,seqcmd)
        except:
            self.errorLog('Diploidocus.purgeHaplotigs() error')
            return False
#########################################################################################################################
    # def joinTables(self,name='',join=[],newkey=[],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True,warnings=True):   ### Makes a new table using join of [(Table,Field[,Fieldlist])]
    #     '''
    #     Makes a new table by joining existing tables, using the given fields, in priority order listed. If a Field from a
    #     latter table already exists, it will be added as "Table_Field" instead. If multiple combinations arise for the
    #     same new key, only the first will be kept.
    #     >> name:str [''] = Name for new table. If not given will become "TableX"
    #     >> join:list of (Table,Field[,Fieldlist]) tuples, where Table is a table name and Field is a Field name or
    #         formula to be used for the join. Fieldlist is an optional list of Fields from that Table to include in the
    #         new table. If Field does not exist, it will be added. (Field may be a Formula.)
    #     >> newkey:list [] = If None, will make a new "AutoID" key field. Else, will use the given Fields.
    #     >> cleanup:bool [True] = If True will delete any Fields generated just to make the join
    #     >> delimit:str ['\t'] = Delimiter to be used to join the key fields
    #     >> empties:bool [True] = Whether to keep entries that do not link to 1+ tables with empty values or delete them.
    #     >> check:bool [False] = Whether to check for entries that are not being joined.
    #     >> keeptable:bool [True] = Whether to add new table to self.list['Tables']
    #     '''
#########################################################################################################################
    ### <7> ### Telomere methods                                                                                        #
#########################################################################################################################
### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#     telomere=X      : Basic telomere sequence for search [TTAGGG]
#     telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
#     teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
# parser.add_argument("-w", "--window", type=int, help="This defines the number of first and last nucleotides that will get scanned for telomeric repeats (default: 50).")
# parser.add_argument("-c", "--cutoff", type=float, help='''A telomere is detected if >= c%% of the first (last) nucleotides are telomeric repeats (default: 50%%).''')
#########################################################################################################################
    def findTelomere(self,sequence):    ### Looks for telomeres in nucleotide sequence using telomere regex
        '''
        Looks for telomeres in nucleotide sequence using telomere regex search of sequence ends. Returns a dictionary of
        whether the ends have telomeres and how much was trimmed off as being Ns (5' and 3').
        Based on https://github.com/JanaSperschneider/FindTelomeres.
        >> sequence:str = DNA sequence to search
        << returns dictionary of {'tel5':T/F,'tel3':T/F,'trim5':INT,'trim3':INT}
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            trim5 = 0
            trim3 = 0
            telomere_at_start, telomere_at_end = False, False
            tel_forward, tel_reverse = self.getStrUC('TeloFwd'), self.getStrUC('TeloRev')
            sequence = sequence.upper()
            WINDOW = self.getInt('TeloSize')
            REPEAT_CUTOFF = self.getNum('TeloPerc')

            for index, position in enumerate(sequence):
                if position != 'N':
                    trim5 = index
                    break
            start_of_sequence_withoutNs = trim5

            for index, position in enumerate(reversed(sequence)):
                if position != 'N':
                    trim3 = index
                    break
            end_of_sequence_withoutNs = len(sequence) - trim3

            # Look for telomeric repeats at the start of the sequence
            telomeric_repeats = re.findall(tel_forward, sequence[start_of_sequence_withoutNs:start_of_sequence_withoutNs+WINDOW])
            # Calculate the % of nucleotides that are part of telomeric repeats
            percent_telomeric_repeats_start = 100.0*sum([len(repeat) for repeat in telomeric_repeats])/float(WINDOW)

            # Look for telomeric repeats at the end of the sequence
            telomeric_repeats = re.findall(tel_reverse, sequence[(end_of_sequence_withoutNs-WINDOW):end_of_sequence_withoutNs])
            # Calculate the % of nucleotides that are part of telomeric repeats
            percent_telomeric_repeats_end = 100.0*sum([len(repeat) for repeat in telomeric_repeats])/float(WINDOW)

            # If more than half of nucleotides at the start/end are telomeric repeats
            if percent_telomeric_repeats_start >= REPEAT_CUTOFF:
                telomere_at_start = True
            if percent_telomeric_repeats_end >= REPEAT_CUTOFF:
                telomere_at_end = True

            # Calculate total percentage telomeres (does not enforce terminal sequences)
            telperc = 0.0
            if telomere_at_start or telomere_at_end:
                telomeric_repeats = re.findall(tel_forward, sequence) + re.findall(tel_reverse, sequence)
                telperc = 100.0 * sum([len(repeat) for repeat in telomeric_repeats]) / float(len(sequence) - sequence.count('N'))

            return {'Tel5':telomere_at_start, 'Tel3':telomere_at_end, 'Trim5':trim5, 'Trim3':trim3, 'TelPerc':telperc}
        except:
            self.errorLog('Diploidocus.findTelomere() error'); raise
#########################################################################################################################
    def findTelomeres(self,save=True,keepnull=False):
        '''
        ## Telomere finding [runmode=telomere]

        Canonical motif is TTAGGG/CCCTAA, but one might see variation, so this searches for regex based on the code at
        https://github.com/JanaSperschneider/FindTelomeres.


        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            #i#teldb = self.db().addEmptyTable('telomeres',['Name','SeqLen','Tel5','Tel3','Trim5','Trim3','TelPerc'],['Name'],log=self.debugging())
            telfile = '{}.telomeres.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and rje.checkForFiles(filelist=[telfile],basename='',log=self.log):
                teldb = db.addTable(telfile,name='telomeres',mainkeys=['Name'])
                teldb.dataFormat({'SeqLen':'int','Trim5':'int','Trim3':'int','TelPerc':'num'})
                return teldb
            forks = self.getInt('Forks')
            ## ~ [1a] ~ Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.exists(self.getStr('SeqIn')):
                raise IOError('Diploidocus Telomere mode needs input assembly (seqin=FILE)')
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F'])
            ## ~ [1b] ~ Results table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            teldb = self.db().addEmptyTable('telomeres',['Name','SeqLen','Tel5','Tel3','Trim5','Trim3','TelPerc'],['Name'],log=self.debugging())
            telomeres = []  # List of sequences with telomeres
            tel5 = tel3 = telboth = 0

            ### ~ [2] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = seqin.seqNum()
            while seqin.nextSeq():
                self.progLog('\r#TELO','Analysing {} sequences for telomeric repeats: {:.2f}%'.format(rje.iStr(stot),sx/stot)); sx += 100.0
                sname = seqin.shortName()
                sequence = seqin.seqSequence()
                tentry = teldb.addEntry(rje.combineDict({'Name':sname,'SeqLen':len(sequence)},self.findTelomere(sequence)))
                # Add reporting in verbose mode?
                if tentry['Tel5'] or tentry['Tel3']:
                    telomeres.append(sname)
                    if tentry['Tel5']: tel5 += 1
                    if tentry['Tel3']: tel3 += 1
                    if tentry['Tel5'] and tentry['Tel3']: telboth += 1
                    if self.v() > 0 and tentry['Tel5']:
                        self.printLog('\r#TELO','5\' Telomeric repeat found in {} ({} terminal Ns)'.format(sname,tentry['Trim5']))
                    if self.v() > 0 and tentry['Tel3']:
                        self.printLog('\r#TELO','3\' Telomeric repeat found in {} ({} terminal Ns)'.format(sname,tentry['Trim3']))
            if telomeres: self.printLog('\r#TELO','Telomeric repeats found in {} of {} sequences: {} at 5\', {} at 3\' ({} both)'.format(rje.iLen(telomeres),rje.iStr(stot),tel5,tel3,telboth))
            else: self.printLog('\r#TELO','Telomeric repeats found in {} of {} sequences.'.format(rje.iLen(telomeres),rje.iStr(stot)))

            ### ~ [3] Save Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not keepnull:
                for ekey in list(teldb.dict['Data'].keys()):
                    if ekey not in telomeres: teldb.dict['Data'].pop(ekey)
            if save: teldb.saveToFile(sfdict={'TelPerc':4})

            return teldb
        except:
            self.errorLog('Diploidocus.findTelomeres() error')
            return False
#########################################################################################################################
### End of SECTION II: Diploidocus Class                                                                                #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: Diploidocus(mainlog,['dna=T']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
