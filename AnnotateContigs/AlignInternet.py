import sys
import time
from Consensus import ConsensusSeq
from MergeLocations import MergeLocations
from TerminalCommands import CheckFileName, CheckModule, FPrint

try:
    from urllib.error import HTTPError   # for Python 3
except ImportError as err:
    from urllib2 import HTTPError        # for Python 2

try:
    CheckModule("Bio")
except:
    sys.path.append("C:/Python27/Lib/site-packages") # for me    
    
from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.GenBank import _FeatureConsumer
from Bio.GenBank.utils import FeatureValueCleaner
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

def ParseRecord(record, description):
    """
    Splices record sequence according to description splicing instructions
    Input:
        reference(SeqRecord)        SeqRecord object for reference
        description(string)         String with splicing instructions
    """
    print(description)
    # Some class which contains a method to splice
    consumer = _FeatureConsumer(use_fuzziness=1,
                    feature_cleaner=FeatureValueCleaner())
    
    # Must format location into a string which the consumer can read and splice
    location = description.split(" ")[-1]
    location = location.replace("{","(").replace("}",")")
    location = location.replace("[","").replace("]","")
    location = location.replace("(+)","").replace("(-)","").replace(":","..")
    
    # Doesn't have to be "CDS", change it to whatever may seem suitable
    # However, you must input some string for the .location method to work
    consumer.feature_key("CDS")
    # Input splicing instructions
    consumer.location(location)
    
    # Returns a SeqRecord object (.extract method is the splicing function)
    return consumer.data.features[0].location.extract(record)
    
class ProteinRecord:
    """
    Holds NCBI protein ID, splicing instructions (location), locations of where
    the protein was mapped, and the contigs which map to it.
    Members:
        protein_id          Protein ID
        location            Splicing instructions
        mapped              Locations of where the protein was mapped by contigs
        contigs             List of contigs that map to this protein
    """
    def __init__(self):
        self.protein_id = None
        self.location = None
        self.mapped = None
        self.contigs = []
    
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

class CodingLocationRecord:
    """
    Holds info on protein id and coding locations on a gene
    Initialization:
        protein_id(string)         Protein ID
        locations(list)            Location where mapping starts
            Note: list contains either a SeqFeature or CompoundLocation objects.
            To get at each individual range of contiguous coding sequences, use:
            #Ex For Accession: BK006917, 1st CDS is:
            #join(4332..6749,339984..340042,351193..351281,363963..364216)
            >>> record = CodingLocationRecord(id, location)
            >>> for part in record.location:
            ...        print part,
            ...
            [4331:6749](+) [339983:340042](+) [351192:351281](+) 
            [363962:364216](+)
    """
    def __init__(self, protein_id, location):
        self.protein_id = protein_id
        self.location = location
    
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
        GetAccessionAndVersion(alignment)
            Returns accession and version numbers from alignment
        GetIDAndReference(accession)
            Returns id and reference sequence from accession number
        GetLocations(id)
            Returns a list of coding locations
        GetContigs(alignment)
            Returns a list of contigs that align to the reference
        GetMapping(locations, contigs) #May be depricated now
            Returs a list regarding the extent of mapping
        GetMappingandContigs(alignment, location)
            Returns a list of contigs and the extent of mapping
        GetProteins(alignment, id=None):
            Returns a list of proteins/isoforms transcribed by the gene
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
            temp = alignment.title
            # Usually, gene is surrounded by parentheses
            if temp.find("(") == -1 or temp.find(")") == -1:
                raise Exception
            
            # Usually, gene is in the last of these sets of parentheses
            gene = temp.split("(")[-1].split(")", 1)[0]
            
            # Reject genes with chars that may confuse file creation
            CheckFileName(gene)
            if gene.find(" ") != -1:
                raise Exception
        except:
            gene = "NA"

        if self.debug:
            print("Gene: %s" % gene)
        return gene
    
    def GetAccessionAndVersion(self, alignment):
        """Gets accession and version"""
        # NCBI gets rid of GI numbers in Sept 2016, so try statement
        # treats it as if GI numbers are still in effect, while except
        # statement believes they are phased out
        try:
            temp = alignment.title.split("|", 4)[3].split(".", 1)
        except IndexError:
            temp = alignment.title.split(" ", 1)[0].split(".", 1)
        accession = temp[0]
        version = temp[1]
        if self.debug:
            print("Accession.version: %s.%s" % (accession, version))
        return (accession, version)
    
    def GetIDAndReference(self, accession):
        """Gets reference sequence (from internet) and id"""
        # Get ID
        if self.debug:
            print("Attempting to get id")
        attempt = 1
        while attempt <= 3:
            try:
                net_handle = Entrez.esearch(db="nucleotide",
                                            term="%s[ACCN]" % accession)
                if self.debug:
                    print("Got response, parsing id")
                record = Entrez.read(net_handle)
                if len(record["IdList"]):
                    id = record["IdList"][0]
                    if self.debug:
                        print("ID: %s" % id)
                else:
                    if self.debug:
                        print("Unable to find %s (may be removed from NCBI" %
                              accession)
                    return (None, "")
                break
            except HTTPError as err:
                # NCBI/Internet problem
                if 500 <= err.code <= 599:
                    if self.debug:
                        print("Received error from server %s" % err)
                        print("Attempt %i of 3" % attempt)
                        attempt += 1
                        time.sleep(15)
                else:
                    FPrint(sys.stderr, "NCBI/Internet problem unresolvable,"
                            + "exiting\n")
                    raise
            except IndexError:
                # Improper accession number used
                FPrint(sys.stderr, "Accession number: %s incorrect in XML\n"
                       % accession)
                raise
            except Exception as err:
                # BioPython throws a random exception, try again
                if self.debug:
                    print("Received error from BioPython %s" % err)
                    print("Attempt %i of 3" % attempt)
                    attempt += 1
                    time.sleep(15)
    
        # Get reference sequence
        if self.debug:
            print("Attempting to get ref_seq")
        attempt = 1
        while attempt <= 3:
            try:
                net_handle = Entrez.efetch(db="nucleotide", id=id,
                                           rettype = "fasta", retmode="text")
                if self.debug:
                    print("Got response, parsing ref_seq")
                record = SeqIO.read(net_handle, "fasta")
                ref = str(record.seq)
                break
            except HTTPError as err:
                # NCBI/Internet problem
                if 500 <= err.code <= 599:
                    if self.debug:
                        print("Received error from server %s" % err)
                        print("Attempt %i of 3" % attempt)
                        attempt += 1
                        time.sleep(15)
                else:
                    FPrint(sys.stderr, "NCBI/Internet problem unresolvable,"
                            + "exiting\n")
                    raise
            except Exception as err:
                # BioPython throws a random exception, try again
                if self.debug:
                    print("Received error from BioPython %s" % err)
                    print("Attempt %i of 3" % attempt)
                    attempt += 1
                    time.sleep(15)
        net_handle.close()
        return (id, ref)
    
    def GetLocations(self, id):
        """Returns a list of coding locations"""
        # Don't try this at home kids! $ cd ~; rm -Rf *
        if self.debug:
            print("Attempting to get coding locations")
        attempt = 1
        while attempt <= 3:
            try:
                net_handle = Entrez.efetch(db="nucleotide", id=id,
                                           rettype="gb", retmode="text")
                if self.debug:
                    print("Got response, parsing coding locations")
                record = SeqIO.read(net_handle, "genbank")
                break
            except HTTPError as err:
                # NCBI/Internet problem
                if 500 <= err.code <= 599:
                    if self.debug:
                        print("Received error from server %s" % err)
                        print("Attempt %i of 3" % attempt)
                        attempt += 1
                        time.sleep(15)
                else:
                    FPrint(sys.stderr, "NCBI/Internet problem unresolvable,"
                            + "exiting\n")
                    raise
            except Exception as err:
                # BioPython throws a random exception, try again
                if self.debug:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    attempt += 1
                    time.sleep(15)
        net_handle.close()
        # Parse these locations kthxbye
        locations = []
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    if "protein_id" in feature.qualifiers.keys():
                        locations.append(CodingLocationRecord(
                                         feature.qualifiers["protein_id"][0],
                                         feature.location))
                    else:
                        locations.append(CodingLocationRecord("NA",
                                         feature.location))

        else:
            raise Exception("No features were found, therefore no CDS found")
        if self.debug:
            if not locations:
                print("No CDS's were found")
        return locations
        
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
    
    def GetMapping(self, locations, contigs): # DEPRECATED
        """Returs a list regarding the extent of mapping"""
        # Fix BioPython's broken crap with the -1 and +1 terms
        mapping = [None] * len(locations)
        for i in range(len(locations)):
            for c_rec in contigs:
                for exon in locations[i].location.parts:
                    if c_rec.start - 1 in exon:
                        if c_rec.end - 1 in exon:
                            if mapping[i]:
                                mapping[i] += FeatureLocation(
                                                c_rec.start, c_rec.end, 1)
                            else:
                                mapping[i] = FeatureLocation(
                                                c_rec.start, c_rec.end, 1)
                        else:
                            if mapping[i]:
                                mapping[i] += FeatureLocation(
                                                c_rec.start, exon.end, 1)
                            else:
                                mapping[i] = FeatureLocation(
                                                c_rec.start, exon.end, 1)                                
                    else:
                        if c_rec.end - 1 in part:
                            if mapping[i]:
                                mapping[i] += FeatureLocation(
                                                exon.start + 1, c_rec.end, 1)
                            else:
                                mapping[i] = FeatureLocation(
                                                exon.start + 1, c_rec.end, 1)
        if self.debug:
            for l, m in zip(locations, mapping):
                print("Protein %s: %s" % (l.protein_id, m))
        return mapping
        
    def GetMappingandContigs(self, alignment, location, identity, similarity):
        """Returns a list of contigs and the extent of mapping"""
        final_contigs = []
        temp = []
        for contig in self.GetContigs(alignment, identity, similarity):
            for exon in location.parts:
                # BioPython's FeatureLocation class is stupid. When an exon
                # covers from [183, 295] (including nts at position 183 & 295),
                # the FeatureLocation class created from calling Entrez.efetch()
                # and SeqIO.parse() yields: [182:295](+). Testing whether a 
                # contig falls within this class yields:
                #
                # >>> 181 in exon
                # False
                # >>> 182 in exon
                # True (Should return False for our purposes)
                # >>> 183 in exon
                # True
                # >>> 294 in exon
                # True
                # >>> 295 in exon
                # False (Should return True for our purposes)
                #
                # As a solution, subtract 1 from end points for inclusion tests
                # Took some fiddling really to get it to work as intended
                if contig.start-1 in exon:
                    final_contigs.append(contig)
                    if contig.end-1 in exon:
                        temp.append((contig.start, contig.end))
                    else:
                        temp.append((contig.start, exon.end))
                elif contig.end-1 in exon:
                    final_contigs.append(contig)
                    temp.append((exon.start+1, contig.end))
        mapped = MergeLocations(temp)
        return mapped, final_contigs
    
    def GetProteins(self, alignment, id, identity, similarity):
        """Dangnabbit, get me them darned list o' proteins via them fancy
        alignment records those new young folk be talking about"""
        #contigs = self.GetContigs(alignment)
        #mapped = []
        proteins = []
        # Initialize ProteinRecord
        for record in self.GetLocations(id):
            protein = ProteinRecord()
            protein.protein_id = record.protein_id
            protein.location = record.location
            protein.mapped, protein.contigs = self.GetMappingandContigs(
                                                alignment, protein.location,
                                                identity, similarity)
            proteins.append(protein)
            if self.debug:
                print("Protein: %s; location: %s" % (protein.protein_id,
                                                     protein.location))
                print("Map: %s" % protein.mapped)
                print("Number of contigs: %i" % len(protein.contigs))
                for contig in protein.contigs:
                    if len(contig.seq) > 50:
                        print("Contig: %s..." % contig.seq[:51]),
                    else:
                        print("Contig: %s" % contig.seq),
                    print("%ito%i" % (contig.start, contig.end))
        return proteins

class AlignmentRecord:
    """
    Stores all alignments for a certain gene
    Initialization:
        alignment(Alignment)    Contains info about gene and contigs aligned
        identity(float)         Minimum identity needed to keep contig [0.7]
        similarity(float)       Minimum similarity needed to keep contig [0.7]
        alignmentParser
            (AlignmentParser)   Existing AlignmentParser object to use [None]
        debug(bool)             If True, then print debugging info [False]
    Members:
        gene         Gene name
        accession    Accession number (GI numbers get phased out by Sept 2016)
        version      Version number
        id           ID used to look up NCBI record
        ref_seq      Reference sequence
        proteins     List of proteins/isoforms expressed within the gene
        debug        If True, then print debugging info
    Methods:
        AddContigs(alignment)       Adds contig to list of contigs
        ContigNumber(protein_id)    Returns number of contigs given protein ID
        ProteinNumber()             Returns number of proteins
        ExomeLengths()              Returns dictionary of {entry:protein length}
        LongestProteins(num)        Returns a list of the $num longest proteins
        SetDebug(debug)             Sets debugging mode
    """
    def __init__(self, alignment, identity=0.7, similarity=0.7,
                 alignmentParser=None, debug=False):
        self.debug = debug
        self.gene = None
        self.accession = None
        self.version = None
        self.id = None
        self.ref_seq = ""
        self.proteins = []
        self.identity = identity
        self.similarity = similarity
        if alignmentParser is not None:
            self.ap = alignmentParser
        else:
            self.ap = AlignmentParser(self.debug)
        try:
            self.gene = self.ap.GetGene(alignment)
            self.accession, self.version = self.ap.GetAccessionAndVersion(alignment)
            self.id, self.ref_seq = self.ap.GetIDAndReference(self.accession)
            if self.id is not None:
                self.proteins = self.ap.GetProteins(alignment, self.id,
                                                    self.identity, self.similarity)
        except Exception as err:
            # Report error but don't raise Exception (probably will change this)
            sys.stdout.flush()
            template = "An exception of type {0} occured. Arguments:\n{1!r}"
            message = template.format(type(err).__name__, err.args)
            FPrint(sys.stderr, message)
            FPrint(sys.stderr, "\nUnable to parse gene: %s\n" % self.gene +
                   "Accession: %s\n" % self.accession +
                   "Version: %s\n" % self.version + 
                   "ID: %s\n" % self.id)
            raise

    def AddContigs(self, alignment):
        """Adds contig to internal list of contigs as well as map them"""
        for protein in self.proteins:
            temp_mapped, temp_contigs = self.ap.GetMappingandContigs(alignment,
                                            protein.location,
                                            self.identity,
                                            self.similarity)
            protein.mapped = MergeLocations(protein.mapped + temp_mapped)
            protein.contigs += temp_contigs
        
    def ContigNumber(self, protein_id):
        """Returns number of contigs given protein ID"""
        for protein in self.proteins:
            if protein.protein_id == protein_id:
                return len(protein.contigs)

    def ExomeLengths(self):
        """Returns a dictionary of {entry: protein length}"""
        result = {}
        for i in range(len(self.proteins)):
            length = 0
            for exon in self.proteins[i].location.parts:
                length += len(exon)
            result[i] = length
        return result

    def LongestProteins(self, num):
        """Returns a list of the $num longest proteins"""
        entries = [(k,v) for k,v in self.ExomeLengths().iteritems()]
        entries.sort(key=(lambda (k,v): v), reverse=True)
        result = []
        for entry, length in entries[:num]:
            result.append(self.proteins[entry])
        return result
        
    def ProteinNumber(self):
        """Returns number of proteins"""
        return len(self.proteins)
    
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
        debug                If True, then print debugging info
    Methods:
        AddAlignment(alignment)
            Adds alignment to dictionary
        OutputAlignment(accession, handle, format, min_coverage=0,
                        max_proteins=float("inf"))
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

    def _CheckCriteria(self, protein, handle, min_coverage):
        """Sees if contigs exist and coverage is within limits"""
        # Check empty contigs list
        if not protein.contigs:
            if self.debug:
                print("Contigs empty, not writing to %s" % handle)
            return False
        # Check coverage
        if not protein.mapped:
            if self.debug:
                print("Empty map, not writing to %s" % handle)
            return False
        m_cov = float(sum([p[1]-p[0] for p in protein.mapped]))
        l_cov = float(sum([p.end-p.start for p in protein.location.parts]))
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

    def _CreateFL(self, locations):
        """Creates a CompoundLocation object, which is made of a list
        of FeatureLocation objects for ease of output"""
        result = []
        for l in locations:
            result.append(FeatureLocation(l[0], l[1]))
        return sum(result)
    
    def _CreateGene(self, ref_seq, contigs):
        """Returns gene created by assembling contigs"""
        records = ["-"*(c.start-1) + c.seq + "-"*(len(ref_seq)-c.end)
                   for i,c in enumerate(contigs)]
        return "".join(ConsensusSeq(records, self.threshold).GetConsensus())
            
    def _FixLocation(self, location):
        """Increments start by one since parsing it subtracts one somehow"""
        fixed = []
        for part in location.parts:
            temp = FeatureLocation(part.start+1, part.end, part.strand)
            fixed.append(temp)
        return sum(fixed)
    
    def AddAlignment(self, alignment):
        """Adds alignment to dictionary"""
        if self.debug:
            print("-"*80)
        accession, version = self.ap.GetAccessionAndVersion(alignment)
        if not accession in self.alignments_dict.keys():
            if self.debug:
                print("Adding new record")
            record = AlignmentRecord(alignment, self.identity, self.similarity, 
                                     self.ap, self.debug)
            self.alignments_dict[accession] = record
        else:
            if self.debug:
                print("Updating %s" % accession)
            self.alignments_dict[accession].AddContigs(alignment)

    def Breadth(self, record):
        """Returns bases mapped for a gene"""
        result = 0
        for m in record.mapped:
            result += m[1] - m[0] + 1
        return result

    def OutputAlignment(self, accession, handle_title, format, min_coverage=0,
                        max_proteins=float("inf")):
        """Writes reference and contigs to handle"""
        record = self.alignments_dict[accession]
        num = min(max_proteins, len(record.proteins))
        for protein in record.LongestProteins(num):
            handle = handle_title+"_Protein_"+protein.protein_id+"."+format
            # Check criteria for output
            if not self._CheckCriteria(protein, handle, min_coverage):
                continue
            # OK to output
            records = []
            id = record.gene + "_Accession_" + accession + "_Protein_" + \
                 protein.protein_id
            dscrpt = str(self._FixLocation(protein.location)).replace(" ","")
            temp = SeqRecord(Seq(record.ref_seq, generic_dna), id=id,
                             description=dscrpt)
            records.append(temp)
            for i, contig in enumerate(protein.contigs):
                dscrpt = str(FeatureLocation(contig.start+1, contig.end, 1))
                temp = SeqRecord(Seq(contig.seq, generic_dna),
                                 id=id+"_contig_%i"%(i+1), description=dscrpt)
                records.append(temp)
            SeqIO.write(records, handle, format)
    
    def OutputContigs(self, accession, handle_title, format, min_coverage=0,
                      max_proteins=float("inf")):
        """Writes reference, aligned contigs, and assembled gene"""
        record = self.alignments_dict[accession]
        num = min(max_proteins, len(record.proteins))
        for protein in record.LongestProteins(num):
            handle = handle_title+"_Protein_"+protein.protein_id+"."+format
            # Check criteria for output
            if not self._CheckCriteria(protein, handle, min_coverage):
                return
            # Create fixed gene
            gene = self._CreateGene(record.ref_seq, protein.contigs)
            # OK to output
            id = record.gene + "_Accession_" + accession + "_Protein_" + \
                 protein.protein_id
            dscrpt = str(self._FixLocation(protein.location)).replace(" ","")
            # Hope you like well-structured code
            # Essentially: [ref, contigs (padded with '-'), assembled gene]
            records = [SeqRecord(Seq(record.ref_seq, generic_dna),
                                 id=id,
                                 description=dscrpt)] + \
                      [SeqRecord(Seq("-" * (c.start-1) + \
                                     c.seq + \
                                     "-" * (len(record.ref_seq)-c.end),
                                 generic_dna),
                                 id="contig_%d" % i,
                                 description="%d:%d" % (c.start,c.end))
                           for i,c in enumerate(protein.contigs)] + \
                      [SeqRecord(Seq(gene, generic_dna),
                                 id=id + "_spliced",
                                 description = \
                                     str(self._CreateFL(protein.mapped)))]
            SeqIO.write(records, handle, format)
    
    def OutputFixed(self, accession, handle_title, format, min_coverage=0,
                    max_proteins=float("inf")):
        """Writes reference and 'fixed' reference created by merging contigs"""
        # Only output "$max_proteins" largest proteins
        record = self.alignments_dict[accession]
        num = min(max_proteins, len(record.proteins))
        for protein in record.LongestProteins(num):
            handle = handle_title+"_Protein_"+protein.protein_id+"."+format
            # Check criteria for output
            if not self._CheckCriteria(protein, handle, min_coverage):
                continue
            # Create fixed gene
            gene = self._CreateGene(record.ref_seq, protein.contigs)
            # OK to output
            id = record.gene + "_Accession_" + accession + "_Protein_" + \
                 protein.protein_id
            dscrpt = str(self._FixLocation(protein.location)).replace(" ","")
            records = [SeqRecord(Seq(record.ref_seq, generic_dna), id=id,
                                 description=dscrpt),
                       SeqRecord(Seq(gene, generic_dna), id=id+"_fix",
                                 description = \
                                     str(self._CreateFL(protein.mapped)))]
            SeqIO.write(records, handle, format)
            
    def OutputSpliced(self, accession, handle_title, format, min_coverage=0,
                      max_proteins=float("inf")):
        """Writes spliced reference and assembled gene"""
        # Only output "$max_proteins" largest proteins
        record = self.alignments_dict[accession]
        num = min(max_proteins, len(record.proteins))
        for protein in record.LongestProteins(num):
            handle = handle_title+"_Protein_"+protein.protein_id+"."+format
            # Check criteria for output
            if not self._CheckCriteria(protein, handle, min_coverage):
                continue
            # Create fixed gene
            gene = self._CreateGene(record.ref_seq, protein.contigs)
            # OK to output
            id = record.gene + "_Accession_" + accession + "_Protein_" + \
                 protein.protein_id
            dscrpt = str(self._FixLocation(protein.location)).replace(" ","")
            ref_record = SeqRecord(Seq(record.ref_seq, generic_dna), id=id,
                                   description=dscrpt)
            gene_record = SeqRecord(Seq(gene, generic_dna), id=id+"_spliced",
                                    description = \
                                        str(self._CreateFL(protein.mapped)))
            records = [ParseRecord(ref_record, dscrpt),
                       ParseRecord(gene_record, dscrpt)]
            SeqIO.write(records, handle, format)
    
    def OutputStats(self, handle_title, min_coverage):
        """Writes breadth and depth of coverage per gene in .csv"""
        handle = handle_title + ".csv"
        sorted_keys = self.alignments_dict.keys()
        sorted_keys.sort(key=(lambda x: x.lower()))
        gene_lengths = []
        protein_lengths = []
        total_bases_mapped = []
        total_breadth_of_coverage = []
        proteins = 0

        with open(handle, 'w') as o:
            #o.write("Gene,Length,Bases Mapped,Breadth of Coverage,"
            #        + "Sum of Depths,Depth of Coverage\n")
            o.write("Gene,Accession,Gene Length,Protein ID,Protein Length," + \
                    "Bases Mapped,Breadth of Coverage\n")
            if not sorted_keys:
                o.write("No genes found\n")
                return
            for key in sorted_keys:
                record = self.alignments_dict[key]
                gene_length = len(record.ref_seq)
                proteins_processed = 0
                for protein in record.proteins:
                    # Check criteria for output
                    if not self._CheckCriteria(protein, handle, min_coverage):
                        continue

                    # Write gene name and its length
                    o.write("%s,%s,%i," % 
                            (record.gene, record.accession, gene_length))

                    # Write protein ID and its length
                    protein_length = len(protein.location)
                    protein_lengths.append(protein_length)
                    o.write("%s,%i," % (protein.protein_id, protein_length))

                    # Get the number of bases mapped to this gene
                    bm = self.Breadth(protein)
                    total_bases_mapped.append(bm)
                    o.write("%i," % bm)

                    # Get the breadth of coverage of this gene
                    boc = float(bm) / protein_length
                    total_breadth_of_coverage.append(boc)
                    o.write("%f\n" % boc)
                    proteins_processed += 1
                if proteins_processed:
                    gene_lengths.append(gene_length)
                    proteins += proteins_processed

            # Calculate averages
            if len(gene_lengths):
                o.write("\nAverage,,")
                o.write("%f,," % (float(sum(gene_lengths)) / len(sorted_keys)))
                o.write("%f," % (float(sum(protein_lengths)) / proteins))
                o.write("%f," % (float(sum(total_bases_mapped)) / proteins))
                o.write("%f,"%(float(sum(total_breadth_of_coverage))/proteins))
    
    def PrintData(self, stream=sys.stdout, verbosity=0):
        """Prints some data to stream"""
        for accession in self.alignments_dict.keys():
            FPrint(stream, "-"*80)
            FPrint(stream, "Accession: %s" % accession)
            for protein in self.alignments_dict[accession].proteins:
                
                FPrint(stream, "Protein: %s contains %i contigs" % 
                       (protein.protein_id, len(protein.contigs)))
                if verbosity == 1:
                    for i, contig in enumerate(protein.contigs):
                        if len(contig.seq) > 50:
                            FPrint(stream, "Contig %i: %s... %i to %i" % 
                                  (i+1,contig.seq[:51],contig.start,contig.end))
                        else:
                            FPrint(stream, "Contig %i: %s %i to %i" %
                                   (i+1, contig.seq, contig.start, contig.end))
                FPrint(stream, "Splicing instructions:", protein.location)
                FPrint(stream, "Mapped:", protein.mapped)
                if protein.mapped:
                    cov = sum([p.start-p.end for p in protein.mapped.parts]) \
                        / sum([p.start-p.end for p in protein.location.parts])
                else:
                    cov = 0
                FPrint(stream, "Coverage: %i" % cov)