"""
 fjoin.py is the reference implementation of the feature join 
 algorithm described in the paper:
 "fjoin: Simple and Efficient Computation of Feature Overlap"
 J. Comp. Bio., 13(8), Oct. 2006, pp 1457-1464.

 FJoin finds overlaps (or, in general, proximity-based) pairs
 of features, given two feature sets. The basic outline of the
 algorithm is:
	1. If needed, sort feature set 1
	2. If needed, sort feature set 2
	3. Move a sliding window along the sorted feature sets, 
	detecting and outputting qualifying feature pairs as 
	we go.  For details, see the paper.

 This implementation is a command-line filter. A typical
 invocation might look like this:

    % python fjoin.py -1 besthits.gff -2 genes.gff -s both > output.txt

 For help, type:

    % python fjoin.py -h

 To call FJoin from a Python program:

    from fjoin import FJoin
    ...
    argv = ["-1", "besthits.gff", 
            "-2", "genes.gff", 
            "-s", "both", 
	    "-o", "output.txt"]
    FJoin( argv ).go()

 INPUTS:
   The two input files are specified with the "-1" and "-2"
   options. One of these files may be standard in; use
   "-" as the file name.

        - Inputs must be in delimited, tabular ASCII format.
	- The default separator character is TAB. You can
	  change this with --separator1 and --separator2 options
	- Comment lines begin with a HASH ('#') character. You can
	  change this with the --comment1 and --comment2 options.
	- Blank lines are ignored. 
	- All rows must have the same number of columns, as
	  determined by the first row of data. Rows having a
	  different number of columns are reported and skipped.
	- Column numbers start at 1.
	- By default, inputs are assumed to be in GFF
	  format.  wherein: column 1 contains the chromosome,
	  columns 4 and 5 contain the start and end coordinates,
	  and column 7 contains the strand.
	- You can change these column assignments with
	  the --columns1 and --columns2 options.
	- The --columns1 and --columns2 options also allow
	  fjoin to ignore strand or both strand and chromosome
	  in detecting overlap.

 OUTPUTS:
   By default, output is written to standard out. You may
   supply an output file name using the "-o" option.
   Output is a TAB-delimited ASCII file.
   One output line is generated for each pair of overlapping
   features in the input. For a given overlapping pair, (a,b),
   the output line contains the following columns:
	col 1		V(a,b), i.e., the amount of overlap 
				between a and b
	cols 2..N+1	the N columns from input 1
	cols N+2..N+M+1	the M columns from input 2

 MESSAGES:
   Status messages, run statistics, and errors are written to 
   standard error by default. You may supply a log file name
   using the "-l" option.
 
 REQUIREMENTS:
   Python 2.3 or later
   GNU sort (not required, but will be used if available)

 AUTHOR:
   Joel Richardson, Ph.D.
   Mouse Genome Informatics
   The Jackson Laboratory
   Bar Harbor, Maine 04609
   jer@informatics.jax.org

 DATE:
   May 2006

 VERSIONS:
   1.0 May 2006	Initial release
   1.1 June 2006	Revised based on suggestion by Jim Kadin. Makes
	the algorithm easier to prove correct. Also seems to improve 
	performance a bit though I don't know why. Main changes:
	1. Streams return sentinels after end of file. Sentinals
	have "infinite" start positions.
	2. Simplified loop structure / end conditions.
	Other changes:
	3. the go() method now returns the number of overlaps found
   1.2 March 2007	Implemented several improvements suggested by 
   	Malcolm Cook and Mike Coleman at the Stowers Institute.
	1. Expand support for delimited formats beyond GFF. User
	can specify which columns contain the needed data (chromosome,
	start, end, and strand) as well as the delimiters (e.g.,
	TAB, COMMA, PIPE, etc.) NOTE that options -i and -I have
	been removed; use options --columns1 and --columns2 instead.
	2. Provide an option to sort internally, if the user
	doesn't have GNU sort.
	3. Move comments into doc strings.
"""

#------------------------------------------------------------
#						VERSION
#------------------------------------------------------------
FJOIN="fjoin"
VERSION="1.2"

#------------------------------------------------------------
#						IMPORTS
#------------------------------------------------------------
import sys
import os
import string
import re
import time
from optparse import OptionParser

#------------------------------------------------------------
#						CONSTANTS
#------------------------------------------------------------
NL	= "\n"
TAB	= "\t"
HASH	= "#"
COMMA	= ","
PIPE	= "|"

# Column assignments for default format.
GFFCOLMAP = { 		
    	'chromosome' : 0,
    	'start'      : 3,
    	'end'        : 4,
    	'strand'     : 6
	}

# Error raised if input not sorted on start coordinate
SORT_ERROR = "Input not sorted"

#------------------------------------------------------------
#						CLASS FJOIN
#------------------------------------------------------------
class FJoin:
    """Given two sets of features, X and Y, finds all feature
    pairs, x and y, where x element of X and y element of Y,
    and x and y overlap by a specified amount.
    """

    #-------------------------------------
    def __init__(self, argv):
	"""Creates an FJoin instance from command line args.
	"""
	self.parser = None	# an OptionParser
	self.options = None	# the parsed options
	self.args = None	# the parsed positional args
	self.ofd = sys.stdout	# the output file desc
	self.lfd = sys.stderr	# the log file desc

	self.k = 1		# user's min overlap amount
	self.isPercent = False	# is k a percentage (True),
				# or absolute (False)?

	self.xStream = None	# the X stream
	self.yStream = None	# the Y stream
	self.streams = []	# X and Y as a list

	self.columnMapX = None	# maps names to column indexes for xStream
	self.columnMapY = None	# maps names to column indexes for yStream
	self.ignoreStrand = False # if True, only looks at chr,start,end
	self.ignoreBoth = False # if True, only looks at start,end

	self.key2Wx = {} 	# maps of chr/strand
	self.key2Wy = {}	# to window. per stream.

	self.lastPick = 0	# most recent pick
				# 0==x, 1==y

	self.CADJUST = 1	# coord.system adjustment
	    # For working in discrete coordinate systems,
	    # this should be set to integer 1. For continuous
	    # coordinate systems, this should be set to float 0.0.

	# Statistics
	self.nfeatures=[0,0]		# number of features
	self.wsizes=[0,0]		# sum of window sizes
	self.avgwsizes = [0.0,0.0]	# average window sizes
	self.maxwsizes = [0,0]		# maximum window sizes
	self.nconsidered = 0		# total number pairs considered
	self.nOutputRows = 0		# total number pairs found

	#
	self.parseCmdLine(argv)

    #-------------------------------------
    def initArgParser(self):
	""" Creates the OptionParser, and initializes the options.
	"""
	self.parser = OptionParser(self.usageString())

	self.parser.add_option("--version", 
	    action="store_true", dest="pVersion", default=False,
	    help="Print program name and version, and exit.")

	self.parser.add_option("-k", metavar="[+-]N[%]",
	    action="store", dest="k", default=None,
	    help="Minimum overlap amount. (Default: 1) " + \
	    "Examples: -k100 : overlap of at least 100 bases. " + \
	    "-k-1000 : separated by no more than 1 kb. " 
	       "-k50% : overlap of >= 50% of longer of pair. " + \
	       "-k-85% : overlap of >= 85% of shorter of pair. " + \
	       "")

	self.parser.add_option("-c", 
	    action="store_true", dest="continuous", default=False,
	    help="Use continuous coordinate system. " + \
	    "(Default: Uses discrete coordinate system.)")

	self.parser.add_option("-1", metavar="FILE",
	    action="store", dest="file1", default=None,
	    help="First input file. (Required.) Use \"-\" for standard in.")

	self.parser.add_option("-2", metavar="FILE",
	    action="store", dest="file2", default=None,
	    help="Second input file. (Required.) Use \"-\" for standard in.")

	self.parser.add_option("-s", metavar="1|2|both",
	    action="store", dest="sort", default=None,
	    help="Sort the input(s). Specify 1, 2, or both. " + \
	    "(Default: no sorting done.)")

	self.parser.add_option("--sortInternal", 
	    action="store_true", dest="sortInternal", default=False,
	    help="Forces any sorting to be done internally. " + \
	    "(Default: uses GNU sort, if available)")

	self.parser.add_option("-o", metavar="FILE",
	    action="store", dest="ofile", default="-",
	    help="Output file. (Default: standard out)")

	self.parser.add_option("-l", metavar="FILE",
	    action="store", dest="lfile", default="-",
	    help="Log file. This is where messages, stats, and errors go. " + \
	    "(Default: standard error)")

	self.parser.add_option("--columns1", metavar="COLS",
	    action="store", dest="columns1", default="",
	    help="Tells fjoin where to find chromosome,start,end,strand in file1. " + \
	    "COLS is a comma-separated list of column indices; column numbers " + \
	    "start at 1. The list may contain 2, 3, or 4 integers, " + \
	    "as follows: " + \
	    "If 2 integers (e.g., --columns1=4,5), they specify the start and end " + \
	    "columns, respectively; chromosome and strand are ignored. " + \
	    "If 3 integers (e.g., --columns1=1,4,5), they specify the chromosome, " + \
	    "start, and end columns, respectively; strand is ignored. " + \
	    "If 4 integers (e.g., --columns1=1,4,5,6), they specify the chromosome, " + \
	    "start, end, and strand columns, respectively. " + \
	    "(Default: assumes GFF format, which is equivalent to: --columns1=1,4,5,7)")

	self.parser.add_option("--columns2", metavar="COLS",
	    action="store", dest="columns2", default="",
	    help="See help for option --columns1. " + \
	    "(Default: assumes GFF format, which is equivalent to: --columns2=1,4,5,7)")

	self.parser.add_option("--separator1", metavar="CHAR",
	    action="store", dest="separator1", default=TAB,
	    help="The column separator for input 1. (Default: TAB character)")
	    
	self.parser.add_option("--separator2", metavar="CHAR",
	    action="store", dest="separator2", default=TAB,
	    help="The column separator for input 2. (Default: TAB character)")

	self.parser.add_option("--comment1", metavar="STRING",
	    action="store", dest="comment1", default=HASH,
	    help="The character or string signifying a comment line in input 1. " + \
	    "A comment line begins with STRING and ends at the next newline. " + \
	    "Comment lines are ignored and are not preserved in the output. " + \
	    "(Default: the HASH character, '#')")
	    
	self.parser.add_option("--comment2", metavar="STRING",
	    action="store", dest="comment2", default=HASH,
	    help="The character or string signifying a comment line in input 2. " + \
	    "(Default: the HASH character, '#')")
	    
    #-------------------------------------
    def validate(self):
	"""Checks that args are valid, opens files, and generally gets set up.
	"""
	if self.options.lfile == "-":
	    self.lfd = sys.stderr
	else:
	    self.lfd = open(self.options.lfile, 'w')

	if self.options.pVersion:
	    self.printVersion()

	self.log("Starting fjoin...\n", True)

	if self.options.file1 is None or self.options.file2 is None:
	    self.parser.error("Must specify both input files.")
	if self.options.file1 == "-" and self.options.file2 == "-":
	    self.parser.error("Cannot use stdin for both inputs.")

	self.columnMapX = self.mkColumnMap(self.options.columns1)
	self.columnMapY = self.mkColumnMap(self.options.columns2)

	if self.options.sort not in [None,"1","2","both"]:
	    self.parser.error("-s argument must be '1', '2', or 'both'.")

	sort1 = self.options.sort in ["1","both"]
	if sort1:
	    if self.options.sortInternal:
		sort1 = "internal"
		self.log("Internally sorting input: %s\n"%self.options.file1)
	    else:
		self.log("Sorting input: %s\n" % self.options.file1)

	sort2 = self.options.sort in ["2","both"]
	if sort2:
	    if self.options.sortInternal:
		sort2 = "internal"
		self.log("Internally sorting input: %s\n"% self.options.file2)
	    else:
		self.log("Sorting input: %s\n" % self.options.file2 )


	self.xStream = \
	    FJoinStream(self.options.file1, sort1, self.options.continuous, \
	    	self.columnMapX, self.options.separator1, self.options.comment1, self )
	self.yStream = \
	    FJoinStream(self.options.file2, sort2, self.options.continuous, \
	    	self.columnMapY, self.options.separator2, self.options.comment2, self )
	self.streams = [self.xStream, self.yStream]

	if self.options.ofile == "-":
	    self.ofd = sys.stdout
	else:
	    self.ofd = open(self.options.ofile, 'w')


	self.isPercent = False
	if self.options.k is None:
	    # default is 1 base overlap
	    self.k = 1
	elif self.options.k[-1] == "%":
	    # percent overlap
	    self.k = int(self.options.k[:-1])
	    if self.k < -100 or self.k > 100:
		self.parser.error( \
			"Percentage must be between 0 and 100, inclusive.")
	    self.isPercent = True
	else:
	    self.k = int(self.options.k)

	if self.options.continuous:
	    self.CADJUST = float(0)
	else:
	    self.CADJUST = int(1)

    #-------------------------------------
    def go(self):
    	"""fjoin's main loop. On each iteration, advance
	the lagging stream, and scan the opposite window.
    	Returns number of overlapping pairs output.
	"""
	x = self.xStream.next()
	y = self.yStream.next()
	while x is not self.xStream.sentinel \
	or  y is not self.yStream.sentinel:
	    if x.start < y.start:
		self.lastPick = 0
		self.scan(x, self.getWindow(x, self.key2Wx), 
			  y, self.getWindow(x, self.key2Wy))
		x = self.xStream.next()
	    else:
		self.lastPick = 1
		self.scan(y, self.getWindow(y, self.key2Wy),
			  x, self.getWindow(y, self.key2Wx)) 
		y = self.yStream.next()

	self.finishStats()
	return self.nOutputRows

    #-------------------------------------
    def scan(self, f, Wf, g, Wg):
	"""Scans window Wg for features overlapping f.
	May remove features from Wg. May add f to Wf.
	"""
	self.accumStats(Wg)
	ig = Wg.iterator()
	for g2 in ig:
	    self.nconsidered += 1
	    if self.leftOf( g2, f.start ):
		ig.remove()
	    else:
	    	v = self.overlaps(f, g2)
		if v is not None:
		    self.output(v, f, g2)
	#
	if not self.leftOf(f, g.start):
	    Wf.append(f)

    #-------------------------------------
    def output(self, v, f1, f2):
	"""Outputs an overlapping pair. v is the
	amount of overlap. f1 and f2 are the 
	features. Reverses order of f1/f2 if
	necessary to maintain consistent ordering
	of x's and y's.
	"""
	if self.lastPick == 1:
	    f1,f2 = f2,f1
	outrow = [str(v)] + f1.row + f2.row
	s = TAB.join(outrow) + NL
	self.ofd.write(s)
	self.nOutputRows += 1

    #-------------------------------------
    def overlaps(self, a, b):
	"""If a and b overlap by at least k (specified on cmd line),
	returns the amount of overlap. Otherwise, returns None. 
	"""
	k = self.k
	if self.isPercent:
	    la = self.length(a)
	    lb = self.length(b)
	    if k < 0:
		k = (-k * min(la,lb))/100
	    else:
		k = (k * max(la,lb))/100
	v = min(a.end,b.end) - max(a.start,b.start) + self.CADJUST
	if v >= k:
	    return v
	else:
	    return None

    #-------------------------------------
    def length(self, a):
	"""Returns the length of feature a.
	"""
	return a.end - a.start + self.CADJUST

    #-------------------------------------
    def leftOf(self, a, p):
	"""Returns true iff a is left of position p, adjusted for k.
	"""
	k = self.k
	if self.isPercent:
	    k = 0
	return (a.end < (p + k - self.CADJUST))

    #-------------------------------------
    def timestamp(self):
	"""Returns the current time, formatted as a string.
	"""
	return time.asctime(time.localtime(time.time()))

    #-------------------------------------
    def log(self, s, ts=False):
	"""Writes a message to the log. If ts is True, prepends a timestamp.
	"""
	if ts:
	    self.lfd.write(NL+"Log message: "+self.timestamp()+NL)
	self.lfd.write(s)

    #-------------------------------------
    def parseCmdLine(self, argv):
	"""Parses and validates command line arguments.
    	"""
        self.initArgParser()
        (self.options, self.args) = self.parser.parse_args(argv)
        self.validate()

    #-------------------------------------
    def usageString(self):
	"""Returns a string describing how to invoke this program.
	"""
	return "\n\t%prog -h\n\t%prog -v\n\t%prog [options] -1 FILE -2 FILE"

    #-------------------------------------
    def getVersion(self):
	"""Returns name and version string for this implementation of fjoin.
	"""
	return "%s-%s" % (FJOIN,VERSION)
    #-------------------------------------
    def printVersion(self):
	"""Prints the version number and exits.
	"""
	vers = self.getVersion()
	print vers
	sys.exit(0)

    #-------------------------------------
    def mkColumnMap(self, colsArg):
	if colsArg == "":
	    return GFFCOLMAP

	colMap = {}
	cols = map(int, colsArg.split(COMMA))
	if len(cols) == 2:
	    colMap['start']      = cols[0] - 1
	    colMap['end']        = cols[1] - 1
	    self.ignoreBoth = True
	elif len(cols) == 3:
	    colMap['chromosome'] = cols[0] - 1
	    colMap['start']      = cols[1] - 1
	    colMap['end']        = cols[2] - 1
	    self.ignoreStrand = True
	elif len(cols) == 4:
	    colMap['chromosome'] = cols[0] - 1
	    colMap['start']      = cols[1] - 1
	    colMap['end']        = cols[2] - 1
	    colMap['strand']     = cols[3] - 1
	else:
	    self.parser.error("Must specify 2, 3, or 4 columns.")

	return colMap

    #-------------------------------------
    def getWindow(self, f, dict):
	"""Returns the window for this key within dict.
	If not found, creates/enters a new window
	under this key.
	"""
	key = self.makeFKey(f)
	if not dict.has_key(key):
	    w = DLL()
	    dict[key] = w
	else:
	    w = dict[key]
	return w

    #-------------------------------------
    def makeFKey(self, f):
	"""Returns key constructed from feature f's 
	chromosome and strand. (Used for finding the right
	window for this feature.)
	"""
	c = f.chromosome
	s = f.strand
	if self.ignoreBoth:
	    c = "*"
	if self.ignoreBoth or self.ignoreStrand:
	    s = "*"
	return (c,s)

    #-------------------------------------
    def accumStats(self, W):
	"""Accumulates statistics for window W.
	"""
	ws = len(W)
	i = 1-self.lastPick
	self.nfeatures[ self.lastPick ] += 1
	self.wsizes[ i ] += ws
	self.maxwsizes[ i ] = max( self.maxwsizes[i], ws )

    #-------------------------------------
    def finishStats(self):
	"""Finalizes and reports statistics.
	"""
	if self.nfeatures[1] > 0:
	    self.avgwsizes[0] = float(self.wsizes[0]) / self.nfeatures[1]

	if self.nfeatures[0] > 0:
	    self.avgwsizes[1] = float(self.wsizes[1]) / self.nfeatures[0]

	self.log("Finished merge...\n", True)
	self.log("Total # features (X,Y): (%d,%d)\n" % \
		(self.nfeatures[0],self.nfeatures[1]))
	self.log("Pairs considered: %d\n" % self.nconsidered)
	self.log("Pairs output: %d\n" % self.nOutputRows)
	if self.nconsidered > 0:
	    self.log("Efficiency: %1.2f\n" % (float(self.nOutputRows)/self.nconsidered))
	self.log("Average window size (X,Y): (%1.2f,%1.2f)\n" % \
		(self.avgwsizes[0], self.avgwsizes[1]))
	self.log("Maximum window size (X,Y): (%d,%d)\n" % \
		(self.maxwsizes[0], self.maxwsizes[1]))

    #-------------------------------------
    def closeFiles(self):
    	"""Closes the input streams. 
	"""
	self.xStream.close()
	self.yStream.close()


#------------------------------------------------------------
#					CLASS FJOINSTREAM
#------------------------------------------------------------
class FJoinStream:
    """Encapsulates a sorted stream of features.
    Option to sort the stream, using the assumed sort
    utilities (see globals: SORT_PATH and SORT_CMD)
    Maintains current position and checks sortedness.
    
    Modifications:
      June 2006 - next() returns sentinel feature after 
         last real feature has been returned.
    """

    #-------------------------------------
    def __init__(self, fname, sortit, contin, columnMap, separator, commentStr, logger):
	"""Initialize a stream from the given file name. If sortit
	is true, sort the stream on start coordinate. If contin is
	true, read the coordinates as floats; otherwise, as ints.
	"""
	self.isClosed = False
	self.logger = logger
	self.separator = separator
	self.commentStr = commentStr
	self.currentPos = None
	self.currentLineNum = 0
	self.currentRowNum = 0
	self.continuous = contin
	self.finished = False
	self.columnMap = columnMap
	self.sentinel = Feature( ["9999999999"] * 9, False, GFFCOLMAP )
	self.sortIt = sortit
	self.sortInternal = False
	self.sortInternalCounter = 0
	self.sortCommand = None

	if fname == "-":
	    self.fname = "<stdin>"
	else:
	    self.fname = fname

	if self.sortIt:
	    self.initSortCmd()
	    if self.sortInternal:
		# Try to load the file and
	        # sort internally.
		self.loadNsort()
	    else:
	        # External sort. Invoke the sort command on the input,
		# and read its output from a pipe.
		if fname == "-":
		    self.fd = os.popen(self.sortCommand)
		else:
		    self.fd = os.popen(self.sortCommand + " " + fname)
	else:
	    if fname == "-":
		self.fd = sys.stdin
	    else:
		self.fd = open(fname, 'r')

    #-------------------------------------
    def close(self):
    	if not self.isClosed:
	    self.isClosed = True
	    self.fd.close()

    #-------------------------------------
    def initSortCmd(self):
	"""Checks to see if we can sort using GNU sort, or if we have
	to sort internally. If we can use GNU sort, initializes self.sortCommand.
	Otherwise, leaves it as None, and sets self.sortInternal to True.
	"""
	if self.sortIt == "internal":
	    self.sortInternal = True
	    return

	p = os.popen('sort --version')
	response = p.read()
	p.close()
	if "not found" in response or "GNU" not in response:
	    self.sortInternal = True
	    self.logger.log(
	      "Did not find GNU sort; attempting internal sort for file: %s.\n"%self.fname)
	else:
	    istart = self.columnMap['start']
	    self.sortCommand = "sort -k%dn -t \"%s\" " % (istart+1,self.separator[0])
	    self.logger.log("SORT CMD: " + self.sortCommand + NL)

	
    #-------------------------------------
    def loadNsort(self):
	"""Loads the input file into a list of features, and sorts them internally."""
	if self.fname == "<stdin>":
	    self.fd = sys.stdin
	else:
	    self.fd = open(self.fname, 'r')
	self.allFeatures = []
	f = self.makeFeature(self.__nextRow())
	while(f is not self.sentinel):
	    self.allFeatures.append(f)
	    f = self.makeFeature(self.__nextRow())
	self.fd.close()

	self.allFeatures.sort( key = lambda x : x.start )

    #-------------------------------------
    def getCurrentPosition(self):
	"""Returns the start coordinate of the most recently returned feature."""
	return self.currentPos

    #-------------------------------------
    def __nextRow(self):
	"""Reads next row of data, and returns it as a list
	of tokens. Skips blank lines and comment lines.
	"""
	if self.finished:
	    return None
	line = self.fd.readline()
	while line:
	    self.currentLineNum += 1
	    if not (line==NL or line.startswith(self.commentStr)):
		tokens = string.split(line, self.separator)
		tokens[-1] = tokens[-1][:-1]
		self.currentLineNum += 1
		return tokens
	    line = self.fd.readline()
	#
	self.finished = True
	return None
	
    #-------------------------------------
    def makeFeature(self, row):
	"""Wraps row (a list of tokens) in a Feature object. If row
	is None, returns the sentinel Feature.
	"""
	if row is None:
	    return self.sentinel
	else:
	    return Feature(row, self.continuous, self.columnMap)

    #-------------------------------------
    def next(self):
	"""Advances the stream to the next row and returns it
	as a feature object. If the associated file has reached
	EOF, return the sentinel object.
	"""
	if self.sortInternal:
	    if self.sortInternalCounter < len(self.allFeatures):
		f = self.allFeatures[self.sortInternalCounter]
		self.sortInternalCounter += 1
	    else:
	        f = self.sentinel
	else:
	    f = self.makeFeature(self.__nextRow())

	if self.currentPos is None or self.currentPos <= f.start:
	    self.currentPos = f.start
	else:
	    raise SORT_ERROR, "Inversion (" + str(self.currentPos) + " > " + str(f.start) + ") detected at line " \
	    	+ str(self.currentLineNum) \
		+ " of input " + self.fname
	return f

#------------------------------------------------------------
#					CLASS FEATURE
#------------------------------------------------------------
class Feature:
    """A Feature is a wrapper around the list of column values
    in a line of input. It defines the field names: chromosome,
    start, end, and strand. 
    """ 

    def __init__(self, row, contin, columnMap):
	"""Creates an instance given row, a list of tokens representing
	the data columns, and contin, a flag indicating whether coordinates
	are floating point (True) or int (False). The column map maps from
	predefined names ('chromosome','start','end','strand') to column
	numbers. Keys 'start' and 'end' are required. Keys 'chromosome'
	and 'strand' are optional.
	"""
	self.row = row
	if columnMap.has_key('chromosome'):
	    self.chromosome = row[columnMap['chromosome']]
	else:
	    self.chromosome = '*'

	if columnMap.has_key('strand'):
	    self.strand     = row[columnMap['strand']]
	else:
	    self.strand     = '*'

	if contin:
	    self.start = float(row[columnMap['start']])
	    self.end   = float(row[columnMap['end']])
	else:
	    self.start = int(row[columnMap['start']])
	    self.end   = int(row[columnMap['end']])

#------------------------------------------------------------
#					CLASS DLL
#------------------------------------------------------------
class DLL:
    """A quick doubly-linked-list implementation. Why do this
    when Python has builtin lists?? Because we need a list
    structure that supports O(1) deletion during an iteration.
    This is needed by fjoin's scan() proceedure to remove items
    from the window.  For Python lists, deletion can be O(n).
    
    A DLL object points to first and last list elements,
    and keeps the count of elements.
    Each element is a DLLElement, which points to the
    next and previous elements, and holds the user's data item.
    Class DLLIter allows forward iteration and removal of
    items during iteration.
    """

    def __init__(self):
	"""Creates an empty list."""
	self.first = None
	self.last = None
	self.n = 0

    def __len__(self):
	"""Returns the number of list items."""
	return self.n

    def __iter__(self):
	"""Allows iteration over DLL instances."""
	return DLLIter(self)

    def iterator(self):
	"""Allows iteration over DLL instances."""
	return self.__iter__()

    def append(self, x):
	"""Appends item x to the end of the DLL."""
	x = DLLElement(x)
	x.next = None
	if self.last is None:
	    x.prev = None
	    self.last = x
	    self.first = x
	    self.n = 1
	else:
	    x.prev = self.last
	    x.prev.next = x
	    self.last = x
	    self.n += 1

#------------------------------------------------------------
#					CLASS DLLELEMENT
#------------------------------------------------------------
class DLLElement:
    """Linked list element. Points to next and previous elements,
    and holds users data object.
    """
    def __init__(self, x):
	"""Creates a list elment that refers to x. next/prev
	pointers are None.
	"""
	self.next = None
	self.data = x
	self.prev = None

#------------------------------------------------------------
#					CLASS DLLITER
#------------------------------------------------------------
class DLLIter:
    """Iterator over a DLL instance. Calls to next() return the
    user data object from the next list element. During
    an iteration, a call to remove() will remove the most
    recently returned element.
    """

    def __init__(self, dll):
	"""Initializes iterator over list instance, dll."""
	self.dll = dll		# the doubly-linked list
	self.p = dll.first	# current position cursor
	self.lastp = None	# most recently returned item

    def __iter__(self):
	return self

    def next(self):
	"""Returns list item value under the cursor and advances the cursor."""
	if self.p is None:
	    raise StopIteration
	self.lastp = self.p
	self.p = self.p.next
	return self.lastp.data

    def remove(self):
	"""Removes from list the last item returned."""
	p = self.lastp
	if p is None:
	    raise "Error", "Nothing to remove."
	self.lastp = None # prevent repeated removals

	if p is self.dll.first:
	    self.dll.first = p.next
	if p is self.dll.last:
	    self.dll.last = p.prev
	if p.prev is not None:
	    p.prev.next = p.next
	if p.next is not None:
	    p.next.prev = p.prev
	self.dll.n -= 1

#------------------------------------------------------------
#						 MAIN
#------------------------------------------------------------
if __name__ == "__main__":
    FJoin(sys.argv).go()
