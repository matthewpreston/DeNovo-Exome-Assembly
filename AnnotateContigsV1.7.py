# De novo consensus assembler. Takes contigs from a blasted .xml file, queries
# internet for reference and splicing instructions, & attempts to create exome.
# Five modes (and combinations of the five) are available:
#    -m a OR -m A    Outputs reference and aligned contigs
#    -m c OR -m C    Outputs reference, aligned contigs, and assembled gene
#    -m d OR -m D    Outputs statistical data (breadth of coverage)
#    -m f OR -m F    Outputs reference and assembled gene from contigs
#    -m s OR -m S    Outputs spliced reference and spliced assembled gene
# 
# Author: Matt Preston and Xiaotong Yang
# Created on: May 16, 2016        V1.0
# Revised on: May 17, 2016        V1.1 Added splicing and mapping locations
#             May 22, 2016        V1.2 Added protein records
#             May 26, 2016        V1.3 Added modes, gene assembly, and splicing
#             May 30, 2016        V1.4 Split into modules, added offline mode
#             Jun 22, 2016        V1.5 Added stats for breadth and depth of cov.
#             Jul 26, 2016        V1.6 Upgraded consensus building
#             Sep 21, 2016        V1.7 Added identity and similarity thresholds

help="""\nUsage: [options] -x <blast.xml> -m <MODE>

WORD OF WARNING: When using this program on large jobs (100+ entries in XML 
file), please do so on weekends or between 9:00PM and 5:00AM Eastern time
during weekdays. Failure to comply may result in an IP address being blocked
from accessing NCBI. More info can be found on:
http://www.ncbi.nlm.nih.gov/books/NBK25497/

If using offline mode, go nuts at any time.

Parameters:
    -x FILE, --xml=FILE         Reads annotations from XML file

    -m MODE, --mode=MODE        Use 'a'/'A' for outputting alignment, 'c'/'C' 
                                for seeing contigs assemble into gene, 'f'/'F'
                                for assembled gene, 's'/'S' for spliced gene,
                                and 'd'/'D' for statistical data

Options:
    -a DIR, --alignment_dir=DIR Will output alignments to DIR if 
                                MODE contains "a" or "A" [Alignments/]

    -c FLOAT, --coverage=FLOAT  Minimum coverage threshold for alignment [0.7]
    
    -C DIR, --contigs_dir=DIR   Will output contigs to DIR if MODE contains "c"
                                or "C" [Contigs/]

    -d, --debug                 Sets to debug mode [False]

    -e EMAIL, --email=EMAIL     For NCBI to contact if problem [Mine dammit]

    -f DIR, --fixed_dir=DIR     Will output ref and 'fixed'/assembled gene to
                                DIR if MODE contains "f" or "F" [Fixed/]

    -i FLOAT, --identity=FLOAT  Lower threshold to keep contiguous sequences
                                with sufficient sequence identity (exact
                                matches) [0.7]
                                
    -p NUM, --protein_num=NUM   Max proteins per alignment to be written [all]
                                Note: only largest proteins will be output
    
    -r, --reuse                 Will reuse output directory [False]
    
    -R REF, --reference=REF     Instead of querying NCBI for references, input
                                REF for offline (and faster) use [None]
    
    -s DIR, --spliced_dir=DIR   Will output spliced ref and assembled gene to 
                                DIR if MODE contains "s" or "S" [Spliced/]
    
    -S FLOAT,--similarity=FLOAT Lower threshold to keep contiguous sequences
                                with sufficient sequence similarity (similar
                                matches) [0.7]
    
    -t FLOAT, --threshold=FLOAT Threshold needed to exceed for creating a
                                consensus nucleotide [0.7]
"""

import getopt
import sys
from AnnotateContigs.TerminalCommands import *

try:
    CheckModule("Bio")
except:
    sys.path.append("C:/Python27/Lib/site-packages") # for me

from Bio import SeqIO
from Bio.Blast import NCBIXML

def usage():
    global help
    print(help)

def main(argv):
    opt_str = "ha:c:C:de:f:i:m:p:rR:s:S:t:x:"
    opt_list = [
                "help",
                "alignment_dir=",
                "contigs_dir=",
                "coverage=",
                "debug",
                "email=",
                "fixed_dir=",
                "identity=",
                "mode=",
                "protein_num=",
                "reuse",
                "reference=",
                "threshold=",
                "spliced_dir=",
                "similarity=",
                "xml="
    ]
    try:
        opts, args = getopt.getopt(argv, opt_str, opt_list)
    except getopt.GetoptError:
        FPrint(sys.stderr, "GetoptError occured")
        usage()
        sys.exit(2)
    alignment_dir = "Alignments"
    blast_file = None
    coverage = 0.7
    contigs_dir = "Contigs"
    email = "matthew.preston@mail.utoronto.ca"
    fixed_dir = "Fixed"
    debug = False
    identity = 0.7
    mode = None
    protein_num = float("inf")
    reference = None
    reuse = False
    spliced_dir = "Spliced"
    similarity = 0.7
    threshold = 0.7
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-a", "--alignment_dir"):
            alignment_dir = arg
        elif opt in ("-c", "--coverage"):
            coverage = float(arg)
        elif opt in ("-C", "--contigs_dir"):
            contigs_dir = arg
        elif opt in ("-d", "--debug"):
            debug = True
        elif opt in ("-e", "--email"):
            email = arg
        elif opt in ("-f", "--fixed_dir"):
            fixed_dir = arg
        elif opt in ("-i", "--identity"):
            identity = float(arg)
        elif opt in ("-m", "--mode"):
            mode = list(arg)
        elif opt in ("-p", "--protein_num"):
            protein_num = int(arg)
        elif opt in ("-r", "--reuse"):
            reuse = True
        elif opt in ("-R", "--reference"):
            reference = arg
        elif opt in ("-s", "--spliced_dir"):
            spliced_dir = arg
        elif opt in ("-S", "--similarity"):
            similarity = float(arg)
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)
        elif opt in ("-x", "--xml"):
            blast_file = arg
        else:
            FPrint(sys.stderr, "Unhandled opt: %s, exiting" % opt)
            sys.exit(2)

    # Check parameters
    CheckParam(blast_file, test=CheckFiles, 
               fail_str="Bad blast file: %s" % blast_file,
               exit=True)
    CheckParam(mode, fail_str="Please input mode:(a/A)/(c/C)/(d/D)/(f/F)/(s/S)",
               exit=True)
    CheckParam(coverage, test=(lambda x: 0<=x and x<=1),
               fail_str="Invalid coverage value: %f (0 <= cov <= 1)" % coverage,
               exit=True)
    CheckParam(identity, test=(lambda x: 0<=x and x<=1),
               fail_str="Invalid identity value: %f (0 <= iden <= 1)" 
                   % identity,
               exit=True)
    CheckParam(protein_num, test=(lambda x: x>=0),
               fail_str="Invalid protein number: %f (num => 0)" % protein_num, 
               exit=True)
    CheckParam(similarity, test=(lambda x: 0<=x and x<=1),
               fail_str="Invalid similarity value: %f (0 <= sim <= 1)" 
                   % similarity,
               exit=True)
    CheckParam(threshold, test=(lambda x: x>=0.5),
               fail_str="Invalid threshold: %f (thresh => 0.5)" % threshold, 
               exit=True)
    if reference:
        CheckParam(reference, test=CheckFiles,
                   fail_str="Bad reference file: %s" % reference,
                   exit=True)

    # Set up proper alignment
    if reference:
        from AnnotateContigs import AlignToRef
        ref_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
        alignments = AlignToRef.Alignments(identity, similarity, threshold, debug)
    else: # TODO: update AlignInternet one day...
        from AnnotateContigs import AlignInternet
        AlignInternet.Entrez.email = email
        alignments = AlignInternet.Alignments(identity, similarity, threshold,
											  debug)

    try:
        # Let's go boys
        blast_handle = open(blast_file, 'r')
        blast_records = NCBIXML.parse(blast_handle)

        # Add alignments to dict
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                if reference:
                    alignments.AddAlignment(alignment, ref_dict)
                else:
                    alignments.AddAlignment(alignment)
        
        # Create output folders
        if 'a' in mode or 'A' in mode:
            alignment_dir = CreateDir(alignment_dir, reuse)
        if 'c' in mode or 'C' in mode:
            contigs_dir = CreateDir(contigs_dir, reuse)
        if 'f' in mode or 'F' in mode:
            fixed_dir = CreateDir(fixed_dir, reuse)
        if 's' in mode or 'S' in mode:
            spliced_dir = CreateDir(spliced_dir, reuse)
        
        # Output alignments/assembled gene to respective folder
        for key, record in alignments.iteritems():
            if reference:
                output_title = "%s" % key
            else:
                output_title = "%s_Accession_%s" % (record.gene, key)
            if 'a' in mode or 'A' in mode:
                output_file = alignment_dir + "/" + output_title
                alignments.OutputAlignment(key, output_file, "fasta",
                                           coverage, protein_num)
            if 'c' in mode or 'C' in mode:
                output_file = contigs_dir + "/" + output_title
                alignments.OutputContigs(key, output_file, "fasta",
                                         coverage, protein_num)
            if 'f' in mode or 'F' in mode:
                output_file = fixed_dir + "/" + output_title
                alignments.OutputFixed(key, output_file, "fasta",
                                       coverage, protein_num)
            if 's' in mode or 'S' in mode:
                output_file = spliced_dir + "/" + output_title
                alignments.OutputSpliced(key, output_file, "fasta",
                                         coverage, protein_num)
        
        # Output statistical data
        if 'd' in mode or 'D' in mode:
            alignments.OutputStats("Summary_Results", coverage)
        
        # Finished
        #alignments.PrintData(verbosity=1) #Manually toggle this
        #print("All done")
    except Exception:
        FPrint(sys.stderr, "An exception occured:\n")
        raise
        sys.exit(2)
    finally:
        blast_handle.close()

if __name__ == "__main__":
    main(sys.argv[1:])