import sys
from Consensus import ConsensusSeq
from MergeLocations import MergeLocations
from TerminalCommands import CheckFileName, CheckModule, FPrint

try:
    CheckModule("Bio")
except:
    sys.path.append("C:/Python27/Lib/site-packages") # for me    
    
from Bio import SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

class ContigRecord:
    """
    Holds info on a contig's sequence and it's mapping location on the reference
    gene.
    Initialization:
        seq(string)        Sequence
        start(int)         Location where mapping starts
        end(int)           Location where mapping ends
    """
    def __init__(self, seq, start, end):
        self.seq = seq
        self.start = start
        self.end = end

class AlignmentParser:
    """
    Contains methods to parse alignment into gene, accession, version, id,
    reference sequence, and contig.
    Initialization:
        debug(bool)   If True, then print debugging info
    Members:
        debug         If True, then print debugging info
    Methods:
        CheckIndel(ref, contig)
            Returns a fixed contig with indels removed
        GetGene(alignment)
            Returns gene abbreviation from alignment
        GetReference(alignment, ref_dict)
            Returns reference sequence retrieved from ref_dict
        GetContigs(alignment, identity, similarity)
            Returns a list of contigs that align to the reference
        GetMappingandContigs(alignment, identity, similarity)
            Returns a list of contigs and the extent of mapping
    """
    def __init__(self, debug=False):
        self.debug = debug
    
    def CheckIndel(self, ref, contig):
        """Rids of indels due to sequencing error in contig"""
        # Remove insertion (when aligned, ref will contain '-')
        ref = list(ref)
        contig = list(contig)
        assert len(ref) == len(contig), "CheckIndel: Ref & contig not same size"
        for i in range(len(ref)):
            if ref[i] == '-':
                contig[i] = '!'
        ref = [n for n in ref if n != '-']
        contig = [n for n in contig if n != '!']
        # Remove deletion (fill with 'N')
        for i in range(len(contig)):
            if contig[i] == '-':
                contig[i] = 'N'
        return "".join(contig)
    
    def GetGene(self, alignment):
        """Gets gene"""
        try:
            # Usually syntax is Species_Gene_Stuff...
            # temp = alignment.title.split(" ", 1)[-1].split("_", 2)
            # gene = temp[0] + "_" + temp[1]
            gene = alignment.title.split(" ", 1)[-1].split("_", 2)[1]
            
            # Reject genes with chars that may confuse file creation
            CheckFileName(gene)
            if gene.find(" ") != -1:
                raise Exception
        except:
            gene = "NA"

        if self.debug:
            print("Gene: %s" % gene)
        return gene
    
    def GetReference(self, alignment, ref_dict):
        """Gets reference sequence (from FASTA)"""
        title = alignment.title.split(" ")[-1]
        return str(ref_dict[title].seq)
    
    def GetContigs(self, alignment, identity, similarity):
        """Returns a list of contigs that align to the reference"""
        contigs = []
        for hsp in alignment.hsps:
            contig = self.CheckIndel(hsp.sbjct, hsp.query)
            # Ensure contig has sufficient sequence identity (exact matches)
            if (float(hsp.identities) / hsp.align_length) < identity:
                continue
            # Ensure contig has sufficient sequence similarity (similar matches)
            if (float(hsp.positives) / hsp.align_length) < similarity:
                continue
            # Collect start and end points of the forward contig
            if hsp.sbjct_start < hsp.sbjct_end:
                contigs.append(ContigRecord(contig, hsp.sbjct_start,
                                            hsp.sbjct_end))
            if hsp.sbjct_start > hsp.sbjct_end:
                # Find reverse complement of contig for it to be forward
                temp = Seq(contig, IUPAC.ambiguous_dna)
                temp = temp.reverse_complement()
                contigs.append(ContigRecord(str(temp),
                                            hsp.sbjct_end ,hsp.sbjct_start))
        return contigs
    
    def GetMappingandContigs(self, alignment, identity, similarity):
        """Returns a list of contigs and the extent of mapping"""
        contigs = self.GetContigs(alignment, identity, similarity)
        temp = []
        for contig in contigs:
            temp.append((contig.start, contig.end))
        mapped = MergeLocations(temp)
        return mapped, contigs
    
class AlignmentRecord:
    """
    Stores all alignments for a certain gene
    Initialization:
        alignment(Alignment)    Contains info about gene and contigs aligned
        ref_dict(dict)          Dictionary created from reference FASTA
        identity(float)         Minimum identity needed to keep contig [0.7]
        similarity(float)       Minimum similarity needed to keep contig [0.7]
        alignmentParser
            (AlignmentParser)   Existing AlignmentParser object to use [None]
        debug(bool)             If True, then print debugging info [False]
    Members:
        gene         Gene name
        ref_seq      Reference sequence
        mapped       Locations of where the reference was mapped by contigs
        contigs      List of contigs that map to this gene
        identity     The minimum identity needed to keep contig
        similarity   The minimum similarity needed to keep contig
        ap           An AlignmentParser object to retrieve record attributes
        debug        If True, then print debugging info
    Methods:
        AddContigs(alignment)       Adds contig to list of contigs
        SetDebug(debug)             Sets debugging mode
    """
    def __init__(self, alignment, ref_dict, identity=0.7, similarity=0.7,
                 alignmentParser=None, debug=False):
        self.debug = debug
        self.gene = None
        self.ref_seq = ""
        self.mapped = []
        self.contigs = []
        self.identity = identity
        self.similarity = similarity
        if alignmentParser is not None:
            self.ap = alignmentParser
        else:
            self.ap = AlignmentParser(self.debug)
        try:
            self.gene = self.ap.GetGene(alignment)
            self.ref_seq = self.ap.GetReference(alignment, ref_dict)
            self.mapped, self.contigs = self.ap.GetMappingandContigs(alignment,
                                                                self.identity,
                                                                self.similarity)
        except Exception as err:
            # Report error but don't raise Exception (probably will change this)
            template = "An exception of type {0} occured. Arguments:\n{1!r}"
            message = template.format(type(err).__name__, err.args)
            FPrint(sys.stderr, message)
            FPrint(sys.stderr, "Unable to parse gene: %s " % self.gene)

    def AddContigs(self, alignment):
        """Adds contig to internal list of contigs as well as map them"""
        temp_mapped, temp_contigs = self.ap.GetMappingandContigs(alignment,
                                                                 self.identity,
                                                                 self.similarity)
        self.mapped = MergeLocations(self.mapped + temp_mapped)
        self.contigs += temp_contigs
        
    def SetDebug(self, debug):
        """Set debugging mode"""
        self.debug = debug
    
class Alignments:
    """
    Contains dictionary to store alignments common to one gene. Can access
    individual alignments from accession number (gene abbreviations are
    difficult to parse).
    Initialization:
        identity(float)      The minimum identity needed to keep contig
        similarity(float)    The minimum similarity needed to keep contig
        threshold(float)     The threshold needed to create consensus
                             See Consensus.py in this dir
        debug(bool)          If True, then print debugging info
    Members:
        alignments_dict      Stores accession number and alignment record
        identity             The minimum identity needed to keep contig
        similarity           The minimum similarity needed to keep contig
        threshold            The threshold needed to create consensus
        ap                   An AlignmentParser object to parse gene names
        debug                If True, then print debugging info
    Methods:
        AddAlignment(alignment, ref_dict)
            Adds alignment to dictionary
        Breadth(gene)
            Returns bases mapped and breadth of coverage for a gene
        OutputAlignment(key, handle_title, format, min_coverage=0,
                        max_proteins=float("inf")):
            Writes reference and contigs to handle
        OutputContigs(key, handle_title, format, min_coverage=0,
                      max_proteins=float("inf")):
            Writes reference, aligned contigs, and assembled gene
        OutputFixed(key, handle_title, format, min_coverage=0,
                    max_proteins=float("inf")):
            Writes reference and 'fixed' reference created by merging contigs
        OutputSpliced(key, handle_title, format, min_coverage=0,
                      max_proteins=float("inf")):
            Writes spliced reference and assembled gene
        OutputStats(handle_title):
            Writes breadth and depth of coverage per gene in .csv
        PrintData(stream=sys.stdout)
            Prints some data to stream
    """
    def __init__(self, identity=0.7, similarity=0.7, threshold=0.7,
                 debug=False):
        self.alignments_dict = {}
        self.identity = identity
        self.similarity = similarity
        self.threshold = threshold
        self.ap = AlignmentParser(debug)
        self.debug = debug
    
    def iteritems(self):
        for k, v in self.alignments_dict.iteritems():
            yield (k, v)
    
    def iterkeys(self):
        for k in self.alignments_dict.iterkeys():
            yield k

    def itervalues(self):
        for v in self.alignments_dict.itervalues():
            yield v

    def _CheckCriteria(self, record, handle, min_coverage):
        """Sees if contigs exist and coverage is within limits"""
        # Check empty contigs list
        if not record.contigs:
            if self.debug:
                print("Contigs empty, not writing to %s" % handle)
            return False
        # Check coverage
        m_cov = float(sum([p[1]-p[0]+1 for p in record.mapped]))
        l_cov = float(len(record.ref_seq))
        coverage = m_cov / l_cov
        if coverage < min_coverage:
            if self.debug:
                print("Coverage(%.3f) < threshold(%.3f), not writing to %s" 
                       % (coverage, min_coverage, handle))
            return False
        # OK to output
        if self.debug:
            print("Coverage(%.3f) => threshold(%.3f), writing to file %s" 
                  % (coverage, min_coverage, handle))
        return True
    
    def _CreateGene(self, ref_seq, contigs):
        """Returns gene created by assembling contigs"""
        records = ["-"*(c.start-1) + c.seq + "-"*(len(ref_seq)-c.end)
                   for i,c in enumerate(contigs)]
        return "".join(ConsensusSeq(records, self.threshold).GetConsensus())
        
    def _CreateFL(self, locations):
        """Creates a CompoundLocation object, which is made of a list
        of FeatureLocation objects for ease of output"""
        result = []
        for l in locations:
            result.append(FeatureLocation(l[0], l[1]))
        return sum(result)
    
    def AddAlignment(self, alignment, ref_dict):
        """Adds alignment to dictionary"""
        if self.debug:
            print("-"*80)
        key = self.ap.GetGene(alignment)
        if not key in self.alignments_dict.keys():
            if self.debug:
                print("Adding new record %s" % key)
            self.alignments_dict[key] = AlignmentRecord(alignment, ref_dict,
                                                        self.identity,
                                                        self.similarity,
                                                        self.ap, self.debug)
        else:
            if self.debug:
                print("Updating %s" % key)
            self.alignments_dict[key].AddContigs(alignment)
            
    def Breadth(self, record):
        """Returns bases mapped for a gene"""
        result = 0
        for m in record.mapped:
            result += m[1] - m[0] + 1
        return result
        
    def Depth(self, gene):# DEPRECATED
        """Returns sum of depths and depth of coverage for a gene"""
        result = 0
        for c in self.alignments_dict[gene].contigs:
            result += c.end - c.start + 1
        length = len(self.alignments_dict[gene].ref_seq)
        return (result, float(result) * 100 / length)

    def OutputAlignment(self, key, handle_title, format, min_coverage=0,
                        max_proteins=float("inf")):
        """Writes reference and contigs to handle"""
        handle = handle_title + "." + format
        record = self.alignments_dict[key]
        # Check criteria for output
        if not self._CheckCriteria(record, handle, min_coverage):
            return
        # OK to output
        records = []
        id = record.gene
        dscrpt = str(FeatureLocation(1, len(record.ref_seq)))
        temp = SeqRecord(Seq(record.ref_seq, generic_dna),
                         id=id, description=dscrpt)
        records.append(temp)
        for i, contig in enumerate(record.contigs):
            dscrpt = str(FeatureLocation(contig.start, contig.end))
            temp = SeqRecord(Seq(contig.seq, generic_dna),
                             id=id+"_contig_%i"%(i+1), description=dscrpt)
            records.append(temp)
        SeqIO.write(records, handle, format)

    def OutputContigs(self, key, handle_title, format, min_coverage=0,
                      max_proteins=float("inf")):
        """Writes reference, aligned contigs, and assembled gene"""
        handle = handle_title + "." + format
        record = self.alignments_dict[key]
        # Check criteria for output
        if not self._CheckCriteria(record, handle, min_coverage):
            return
        # Create fixed gene
        gene = self._CreateGene(record.ref_seq, record.contigs)
        # OK to output
        id = record.gene
        dscrpt = str(FeatureLocation(1, len(record.ref_seq)))
        # Hope you like well-structured code
        # Essentially: [ref, contigs (padded with '-'), assembled gene]
        records = [SeqRecord(Seq(record.ref_seq, generic_dna),
                             id=id, description=dscrpt)] + \
                  [SeqRecord(Seq("-"*(c.start-1) + \
                                     c.seq + \
                                     "-"*(len(record.ref_seq)-c.end),
                                 generic_dna),
                             id="contig_%d"%i,
                             description="%d:%d"%(c.start,c.end))
                       for i,c in enumerate(record.contigs)] + \
                  [SeqRecord(Seq(gene, generic_dna), id=id+"_spliced",
                             description=str(self._CreateFL(record.mapped)))]
        SeqIO.write(records, handle, format)

    def OutputFixed(self, key, handle_title, format, min_coverage=0,
                    max_proteins=float("inf")):
        """Writes reference and 'fixed' reference created by merging contigs"""
        handle = handle_title + "." + format
        record = self.alignments_dict[key]
        # Check criteria for output
        if not self._CheckCriteria(record, handle, min_coverage):
            return
        # Create fixed gene
        gene = self._CreateGene(record.ref_seq, record.contigs)
        # OK to output
        id = record.gene
        dscrpt = str(FeatureLocation(1, len(record.ref_seq)))
        records = [SeqRecord(Seq(record.ref_seq, generic_dna),
                             id=id, description=dscrpt),
                   SeqRecord(Seq(gene, generic_dna), id=id+"_fix",
                             description=str(self._CreateFL(record.mapped)))]
        SeqIO.write(records, handle, format)
            
    def OutputSpliced(self, key, handle_title, format, min_coverage=0,
                      max_proteins=float("inf")):
        """Writes spliced reference and assembled gene"""
        handle = handle_title + "." + format
        record = self.alignments_dict[key]
        # Check criteria for output
        if not self._CheckCriteria(record, handle, min_coverage):
            return
        # Create fixed gene
        gene = self._CreateGene(record.ref_seq, record.contigs)
        # OK to output
        id = record.gene
        dscrpt = str(FeatureLocation(1, len(record.ref_seq)))
        records = [SeqRecord(Seq(record.ref_seq, generic_dna),
                             id=id, description=dscrpt),
                   SeqRecord(Seq(gene, generic_dna), id=id+"_spliced",
                             description=str(self._CreateFL(record.mapped)))]
        SeqIO.write(records, handle, format)
    
    def OutputStats(self, handle_title, min_coverage=0):
        """Writes breadth and depth of coverage per gene in .csv"""
        handle = handle_title + ".csv"
        sorted_keys = self.alignments_dict.keys()
        sorted_keys.sort(key=(lambda x: x.lower()))
        lengths = []
        total_bases_mapped = []
        total_breadth_of_coverage = []

        with open(handle, 'w') as o:
            #o.write("Gene,Length,Bases Mapped,Breadth of Coverage,"
            #        + "Sum of Depths,Depth of Coverage\n")
            o.write("Gene,Length,Bases Mapped,Breadth of Coverage\n")
            if not sorted_keys:
                o.write("No genes found\n")
                return
            for gene in sorted_keys:
                record = self.alignments_dict[gene]

                # Check criteria for output
                if not self._CheckCriteria(record, handle, min_coverage):
                    return

                # Get gene name
                o.write("%s," % gene)

                # Get gene length and store
                length = len(record.ref_seq)
                lengths.append(length)
                o.write("%i," % length)

                # Get the number of bases mapped to this gene
                bm = self.Breadth(record)
                total_bases_mapped.append(bm)
                o.write("%i," % bm)

                # Get the breadth of coverage of this gene
                boc = float(bm) / length
                total_breadth_of_coverage.append(boc)
                o.write("%f\n" % boc)

            # Calculate averages
            num = len(sorted_keys)
            if num:
                o.write("\nAverage,")
                o.write("%f," % (float(sum(lengths)) / num))
                o.write("%f," % (float(sum(total_bases_mapped)) / num))
                o.write("%f," % (float(sum(total_breadth_of_coverage)) / num))
            
    def PrintData(self, stream=sys.stdout, verbosity=0):
        """Prints some data to stream"""
        for key in self.alignments_dict.keys():
            FPrint(stream, "-"*80)
            FPrint(stream, "Gene: %s contains %i contigs" % 
                   (key, len(self.contigs)))
            if verbosity == 1:
                for i, contig in enumerate(self.contigs):
                    if len(contig.seq) > 50:
                        FPrint(stream, "Contig %i: %s... %i to %i" % 
                              (i+1,contig.seq[:51],contig.start,contig.end))
                    else:
                        FPrint(stream, "Contig %i: %s %i to %i" %
                               (i+1, contig.seq, contig.start, contig.end))
            FPrint(stream, "Mapped:", self.mapped)
            cov  = float(sum([p[1] - p[0] for p in self.mapped]))\
                / float(len(self.ref_seq))
            FPrint(stream, "Coverage: %i" % cov)