#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Description:  Diploid genome assembly analysis tools
Version:      0.3.0
Last Edit:    06/09/19
Copyright (C) 2017  Richard J. Edwards - See source code for GNU License Notice

Function:
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

    Diploidocus is available as part of SLiMSuite, or via a standalone GitHub repo at
    <https://github.com/slimsuite/diploidocus>.

Commandline:
    ### ~ Main Diploidocus run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    runmode=X       : Diploidocus run mode [sortnr/diphap/diphapnr/purgehap/vecscreen/insilico]
    basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
    genomesize=INT  : Haploid genome size (bp) [13.1e6]
    summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
    ### ~ Assembly processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly (sortnr/diphap modes) []
    ### ~ SortNR filtering/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    checkcov=PERC   : Percentage coverage for double-checking partial exact matches [95]
    seqout=FILE     : Output sequence assembly [$BASEFILE.nr.fasta]
    ### ~ In silico input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
    ### ~ Read filtering options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    rqfilter=X      : Minimum RQ for output subreads [0]
    lenfilter=X     : Min read length for filtered subreads [500]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_rmd, rje_seqlist, rje_sequence, rje_paf #, rje_genomics
import smrtscape
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Fixed bugs with parent basefile, genome size default and Sequel data parsing.
    # 0.2.0 - Added sortnr minimap2 run mode, and diphap primary/alternative assembly splitting.
    # 0.3.0 - Added summarising of sequences post-run. Improved documents. Moved to tools/.
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
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Diploidocus', '0.3.0', 'September 2019', '2017')
    description = 'Diploid genome assembly analysis tools.'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
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
    Diploidocus Class. Author: Rich Edwards (2015).

    Str:str
    - Parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    - Parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    - RunMode=X       : Diploidocus run mode [insilico/sortnr/diphap/vecscreen]
    - SeqIn=FILE      : Input sequence assembly (sortnr/diphap modes) []
    - SeqOut=FILE     : Output sequence assembly [$BASEFILE.fasta]

    Bool:boolean
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - Summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]

    Int:integer
    - GenomeSize=INT  : Haploid genome size (bp) [13.1e6]
    - LenFilter=X     : Min read length for filtered subreads [500]

    Num:float
    - CheckCov=PERC   : Percentage coverage for double-checking partial exact matches [95]
    - RQFilter=X      : Minimum RQ for output subreads [0]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Parent1','Parent2','RunMode','SeqIn','SeqOut','DebugStr']
        self.boollist = ['DocHTML','Summarise']
        self.intlist = ['GenomeSize','LenFilter']
        self.numlist = ['CheckCov','RQFilter']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'DocHTML':False,'Summarise':True})
        self.setInt({'LenFilter':500,'GenomeSize':13.1e6})
        self.setNum({'CheckCov':95.0,'RQFilter':0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'str',['RunMode','DebugStr'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Parent1','Parent2','SeqIn','SeqOut'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DocHTML','Summarise'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['GenomeSize','LenFilter'])   # Integers
                self._cmdReadList(cmd,'float',['RQFilter']) # Floats
                self._cmdReadList(cmd,'perc',['CheckCov']) # Percentage
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main Diploidocus run method
        '''
        # Diploidocus: Diploid genome assembly analysis tools

        This module codes for a number of diploid genome assembly manipulation and analysis functions. The different run
        modes are set using `runmode=X`:

        * `sortnr` performs an all-by-all mapping with minimap2 and then removes redundancy
        * `diphap` splits a pseudodiploid assembly into primary and alternative scaffolds
        * `diphapnr` runs `sortnr` followed by `diphap`
        * `purgehap` [coming soon!] filters scaffolds based on post-processing of purge_haplotigs
        * `vecscreen` [coming soon!] searches for contaminants and flags/masks/removes identified scaffolds
        * `insilico` generates balanced diploid combined reads from two sequenced haploid parents

        See sections below for details of each mode. General SLiMSuite run documentation can be found at
        <https://github.com/slimsuite/SLiMSuite>.

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
        runmode=X       : Diploidocus run mode [sortnr/diphap/diphapnr/purgehap/vecscreen/insilico]
        basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
        genomesize=INT  : Haploid genome size (bp) [13.1e6]
        summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        ### ~ Assembly processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly (sortnr/diphap modes) []
        ### ~ SortNR filtering/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        checkcov=PERC   : Percentage coverage for double-checking partial exact matches [95]
        seqout=FILE     : Output sequence assembly [$BASEFILE.nr.fasta]
        ### ~ In silico input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
        parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
        See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
        ### ~ Read filtering options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        rqfilter=X      : Minimum RQ for output subreads [0]
        lenfilter=X     : Min read length for filtered subreads [500]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        # Diploidocus run modes

        Details for the main Diploidocus run modes are given below.

        **NOTE:** Diploidocus is under development and documentation might be a bit sparse. Please contact the author or
        post an issue on GitHub if you have any questions.

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
        ## Vector contamination screening [runmode=vecscreen]
        _Coming soon!_

        ---
        ## Purge haplotigs post-processing [runmode=purgehap]
        _Coming soon!_

        ---
        ## In silico diploid [runmode=insilico]

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
            if self.getStrLC('RunMode') in ['sortnr','seqnr','nrseq','nr']: return self.sortNR()
            elif self.getStrLC('RunMode') in ['diphap','diphapnr']: return self.dipHap()
            elif self.getStrLC('RunMode') == 'insilico':
                return self.inSilicoHybrid()
            else: raise ValueError('RunMode="%s" not recognised!' % self.getStrLC('RunMode'))
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if not self.baseFile(return_none=''):
                if self.getStrLC('SeqIn'): self.baseFile(rje.baseFile(self.getStr('SeqIn')))
                else: self.baseFile('diploidocus')
            self.printLog('#BASE','Output file basename: %s' % self.baseFile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
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
            parent1 = smrtscape.SMRTSCAPE(self.log,['genomesize=13.1e6']+self.cmd_list+['batch=%s' % self.getStr('Parent1'),'basefile=%s' % base1])
            parent1.setup()
            udb1 = parent1.udb()
            cdb = parent1.db('smrt',add=True,mainkeys=['Name'])
            cdb.dataFormat({'SMRT':'int'})
            cx = cdb.entryNum()
            ## ~ [0a] Parent 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 2 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent2: %s' % self.getStr('Parent2'))
            base2 = rje.baseFile(self.getStr('Parent2'))
            parent2 = smrtscape.SMRTSCAPE(self.log,['genomesize=13.1e6']+self.cmd_list+['batch=%s' % self.getStr('Parent2'),'basefile=%s' % base2])
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
