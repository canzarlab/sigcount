// ==========================================================================
//                                 transemble
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Sandro Andreotti <sandro.andreotti@fu-berlin.de>
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Read a Sam file and generate input for RNASeqQuantification.
// ==========================================================================

#include <iostream>
#include <cstdlib>
#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.

#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/store.h>
#include <seqan/arg_parse.h>

#include <seqan/bam_io/bam_scanner_cache.h>

//#define USE_FRAGSTORE
//#define GENERATEGAPS

using namespace seqan;

size_t INVALID_ID = static_cast<size_t>(-1);

enum ReadMode
{
    SINGLE,
    PAIRED,
    MIXED
};

enum SpliceType
{
    NONE,
    LEFT,
    RIGHT,
    FULL
};

enum GroupMode
{
    GROUPANNOT,
    GROUPMAPPINGS,
    GROUPBOTH
};

struct Options
{
    bool showHelp;
    bool showVersion;
    unsigned verbosity;
    
    CharString samInFile;
    CharString bamIndexFile;
    CharString gtfInFile;
    CharString gtfOutFile;
    CharString countFile;
    CharString statFile;
    CharString matFile;
    
    //intervals with a gap smaller than this are merged
    int minIntronLength;
    int minMatchLength;
    int maxGapLength;
    ReadMode rMode;
    GroupMode gMode;
    bool skipMulti;
    unsigned boundaryTolerance;
    bool useStrandInfo;
    
    Options()
    {
        showHelp = false;
        showVersion = false;
        verbosity = 0;
        minIntronLength = 50;
        minMatchLength = 1;
        maxGapLength = std::numeric_limits<int>::max();
        skipMulti = false;
        rMode = SINGLE;
        gMode = GROUPBOTH;
        boundaryTolerance = 0;
        useStrandInfo = false;
    }
};

struct LibInfo
{
    //mean read length
    double rLenMean;
    //number of reads used for estimating
    unsigned numR;
    //mean fragment length
    double fLenMean;
    //fragment length sdev
    double fLenSdev;
    //number of fragments used for estimating
    unsigned numFLen;
    //total number of mapped fragments
    unsigned totalFrags;
};

void setupCommandLineParser(ArgumentParser & parser, Options const &)
{
    setAppName(parser, "sigCount");
    std::string version = "0.2";
#ifdef SEQAN_REVISION
    version += std::string(" [") + std::string(SEQAN_REVISION) + "]";
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
	setVersion(parser, version);
    
    
    //////////////////////////////////////////////////////////////////////////////
    // Define options
#ifdef SEQAN_HAS_ZLIB
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIalignments.[bam]\\fP> <\\fIannotation.[gtf|gff]\\fP> <\\fIoutfile-prefix\\fP>");
    addDescription(parser, "This program parses short read alignments and generates input for CIDANE");
#else
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIalignments.[sam]\\fP> <\\fIannotation.[gtf|gff]\\fP> <\\fIoutfile-prefix\\fP>");
    addDescription(parser, "This program parses short read alignments and generates input for CIDANE. Install ZLIB and recompile to read bam files");
#endif
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "alignments"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "annotation"));
    setValidValues(parser, 1, seqan::GffFileIn::getFileExtensions());
    
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_PREFIX, "outfile-prefix"));
    setValidValues(parser, 0, "bam");
    //setValidValues(parser, 0, seqan::BamFileIn::getFileExtensions());
    
    /*
     addOption(parser, seqan::ArgParseOption("a", "annotation", "GFF/GTF file containing subexon information",
     ArgParseArgument::INPUT_FILE, "IN"));
     setValidValues(parser, "annotation", seqan::GffFileIn::getFileExtensions());
     
     
     addOption(parser, seqan::ArgParseOption("i", "intron", "minimum required intron size - default: 50",
     ArgParseArgument::INTEGER, "INT"));
     
     
     addOption(parser, seqan::ArgParseOption("m", "mode", "read mode",
     ArgParseArgument::STRING, "STRING"));
     setValidValues(parser, "mode", "SINGLE PAIRED MIXED");
     setDefaultValue(parser, "mode", "MIXED");
     
     addOption(parser, seqan::ArgParseOption("g", "group", "group mode to group subexons to loci",
     ArgParseArgument::STRING, "STRING"));
     setValidValues(parser, "group", "ANNOTATION MAPPINGS BOTH");
     setDefaultValue(parser, "group", "BOTH");
     */
    addOption(parser, seqan::ArgParseOption("s", "strand-specific", "use strand information in provided bam if contained. (Searches for TopHat's XS:A tags)"));
    
    addOption(parser, seqan::ArgParseOption("l", "min-match", "minimum length of every mapped fragment of a read",
                                            ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "min-match", 1);
    
    addOption(parser, seqan::ArgParseOption("b", "boundary-tolerance", "Tolerated deviation of a reads implied exon-intron boundary from an annotated one",
                                            ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "boundary-tolerance", 0);
    
    addOption(parser, seqan::ArgParseOption("u", "max-gap", "maximum gap between split parts of read or between read pairs to be discarded. (default: inf)",
                                            ArgParseArgument::INTEGER, "INT"));
    
    //    addSection(parser, "verbosity");
    //    addOption(parser, seqan::ArgParseOption("v",  "verbose", "verbose mode"));
    //    addOption(parser, seqan::ArgParseOption("vv", "vverbose", "very verbose mode"));
    
    //requiredArguments(parser, 3);
    
}

int parseCommandLineAndCheck(Options & options,
                             ArgumentParser & parser,
                             int argc,
                             char const ** argv)
{
    ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    
    //    if (isSet(parser, "verbose"))
    //        options.verbosity = 2;
    //    if (isSet(parser, "vverbose"))
    //        options.verbosity = 3;
    
    //    if (isSet(parser, "annotation"))
    //        getOptionValue(options.gtfInFile, parser, "annotation");
    
    getArgumentValue(options.gtfInFile, parser, 1);
    
    if (isSet(parser, "strand-specific"))
    {
        options.useStrandInfo = true;
    }
    
    if (isSet(parser, "min-match"))
    {
        getOptionValue(options.minMatchLength, parser, "min-match");
    }
    if (isSet(parser, "max-gap"))
    {
        getOptionValue(options.maxGapLength, parser, "max-gap");
    }
    
    /*
     if (isSet(parser, "mode"))
     {
     std::string modetmp;
     getOptionValue(modetmp, parser, "mode");
     if (modetmp == "SINGLE")
     options.rMode = SINGLE;
     else if (modetmp == "PAIRED")
     options.rMode = PAIRED;
     else if (modetmp == "MIXED")
     options.rMode = MIXED;
     }
     */
    /*
     if (isSet(parser, "group"))
     {
     std::string modetmp;
     getOptionValue(modetmp, parser, "group");
     if (modetmp == "ANNOTATION")
     options.gMode = GROUPANNOT;
     else if (modetmp == "MAPPINGS")
     options.gMode = GROUPMAPPINGS;
     else if (modetmp == "BOTH")
     options.gMode = GROUPBOTH;
     }
     */
    /*
     if (isSet(parser, "intron"))
     getOptionValue(options.minIntronLength, parser, "intron");
     */
    
    if (isSet(parser, "boundary-tolerance"))
        getOptionValue(options.boundaryTolerance, parser, "boundary-tolerance");
    
    getArgumentValue(options.samInFile, parser, 0);
    
    options.bamIndexFile = options.samInFile;
    options.bamIndexFile += ".bai";
    
    
    
    CharString prefix;
    getArgumentValue(prefix, parser, 2);
    
    options.countFile = prefix;
    append(options.countFile, ".cnt");
    options.matFile = prefix;
    append(options.matFile, ".mat");
    options.statFile = prefix;
    append(options.statFile, ".stat");
    options.gtfOutFile = prefix;
    append(options.gtfOutFile, "_refined.gtf");
    
    return 0;
}

typedef long int TIntervalPosType;
struct MyCargo
{
    MyCargo() : id(0), type(NONE){}
    MyCargo(size_t id, SpliceType type) : id(id), type(type){}
    size_t id;
    SpliceType type;
};

typedef unsigned TIntervalCargoType;
//typedef MyCargo TIntervalCargoType;

typedef IntervalAndCargo<TIntervalPosType, TIntervalCargoType> TInterval;
namespace seqan
{
    bool operator< (const TInterval &left, const TInterval &right)
    {
        if (left.i1 == right.i1)
            return left.i2 < right.i2;
        return left.i1 < right.i1;
    }
    bool operator== (const TInterval &left, const TInterval &right)
    {
        return (left.i1 == right.i1 && left.i2 == right.i2);
    }
}

typedef String<TInterval, Alloc<Exact> > TIntervalList;
typedef Iterator<TIntervalList>::Type TIntervalListIterator;
typedef Iterator<TIntervalList const>::Type TIntervalListConstIterator;
typedef std::unordered_map<int, int> TLenDistMap;

struct TSplitReadIntervals
{
    TSplitReadIntervals(){}
    
    TSplitReadIntervals(const TSplitReadIntervals &rhs):
    readId(rhs.readId),
    strand (rhs.strand),
    intervalsLeft(rhs.intervalsLeft, Exact()),
    intervalsRight(rhs.intervalsRight, Exact())
    {
        shrinkToFit(intervalsLeft);
        shrinkToFit(intervalsRight);
    }
    
    unsigned readId;
    char strand;
    TIntervalList intervalsLeft;
    TIntervalList intervalsRight;
};

typedef String<unsigned> TIndexList;
typedef Iterator<TIndexList>::Type TIndexListIterator;
typedef Iterator<TIndexList const>::Type TIndexListConstIterator;

struct TSplitReadIndices
{
    TSplitReadIndices() {}
    TSplitReadIndices(const TSplitReadIndices &rhs):
    readId(rhs.readId),
    indicesLeft(rhs.indicesLeft, Exact()),
    indicesRight(rhs.indicesRight, Exact()),
    minGap(rhs.minGap)
    {
        shrinkToFit(indicesLeft);
        shrinkToFit(indicesRight);
    }
    
    unsigned readId;
    TIndexList indicesLeft;
    TIndexList indicesRight;
    int minGap;
};

typedef String<TSplitReadIntervals> TSplitReadIntervalsList;
typedef Iterator<TSplitReadIntervalsList const>::Type TSplitReadIntervalsListConstIterator;
typedef std::map<CharString, TSplitReadIntervalsList> TSplitReadIntervalsMap;

typedef String<TSplitReadIndices> TSplitReadIndicesList;
typedef Iterator<TSplitReadIndicesList const>::Type TSplitReadIndicesListConstIterator;

typedef std::set<unsigned> TSplitPositionList;
typedef TSplitPositionList::iterator TSplitPositionListIter;
typedef TSplitPositionList::const_iterator TSplitPositionListConstIter;

typedef String<std::pair<TIntervalPosType, TIntervalPosType> > TJunctionList;
typedef Iterator<TJunctionList>::Type  TJunctionListIterator;
typedef Iterator<TJunctionList const>::Type TJunctionListIteratorConst;

typedef std::pair<TIndexList, TIndexList> TCountTuple;
struct TCount
{
    TCount() : totalCount(0.), uniqueCount(0){}
    
    double totalCount;
    unsigned uniqueCount;
};

typedef std::map<TCountTuple, TCount> TTupleCountList;
typedef TTupleCountList::const_iterator TTupleCountListConstIter;

typedef std::map<TCountTuple, std::map<int, double> > TTupleCountGapList;
typedef std::map<int, double>::const_iterator TGapCountMapConstIter;


typedef String<bool> TBoolString;
typedef std::map<unsigned, short> TMultiReadMap;

typedef std::map<CharString, TIntervalList> TSubexonIntervalsMap;
typedef TSubexonIntervalsMap::const_iterator TSubexonIntervalsMapConstIter;
typedef std::map<CharString, unsigned>TReadIdMap;

TReadIdMap readNameStore;
unsigned genericId = 0;

template <
typename TBlockBoundaries,
typename TCigarString,
typename TContigPos
>
inline void
getBlockBoundaries(
                    TBlockBoundaries &blockBoundaries,
                    TCigarString const & cigar,
                    TContigPos startPos)
{
    TContigPos lp = startPos;
    bool gap = false;
    
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        switch (cigar[i].operation)
        {
            case 'D':
            case 'M':
                if (gap || empty(blockBoundaries))
                {
                    appendValue(blockBoundaries, lp);
                    lp += cigar[i].count;
                    appendValue(blockBoundaries, lp);
                    gap = false;
                }
                else
                {
                    lp += cigar[i].count;
                    back(blockBoundaries) += cigar[i].count;
                }
                break;
                
            case 'N':
                gap = true;
                lp += cigar[i].count;
                break;
            case 'I':
                break;
            case 'S':
                break;
            case 'H':
                break;
            case 'P':
                break;
                //lp += cigar[i].count;
        }
    }
}

template <typename TBlockBoundaries>
inline void
removeSmallGaps(TBlockBoundaries &blockBoundaries, const Options & options)
{
    unsigned j = 2;
    for (unsigned i = 3; i < length(blockBoundaries); i += 2)
    {
        if (blockBoundaries[i - 1] - blockBoundaries[i - 2] < options.minIntronLength)
            blockBoundaries[j - 1] = blockBoundaries[i];
        else
        {
            blockBoundaries[j++] = blockBoundaries[i - 1];
            blockBoundaries[j++] = blockBoundaries[i];
        }
    }
    resize(blockBoundaries, j);
}


int groupMultiMapped(seqan::StringSet<seqan::CharString> &/*contigNameStore*/,
                     const Options &options,
                     String<size_t> &recordToId)
{
    typedef Iterator<String<size_t> const, Rooted>::Type TPosListConstIter;
    typedef std::map<CharString, String<size_t> > TSecondaryAlignedMap;
    typedef TSecondaryAlignedMap::iterator TSecondaryAlignedIter;
    
    TSecondaryAlignedMap secAligns;
    
    size_t numPri = 0; // number of primary reads
    size_t numSingleSecs = 0; //number of lone reads marked as secondary
    
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(options.samInFile)))
    {
        std::cerr << "Could not open " << options.samInFile << std::endl;
        return 1;
    }
    
    BamHeader header;
    readHeader(header, bamFileIn);
    
    size_t maxContigId = length(contigNames(context(bamFileIn)));
    if (maxContigId == 0)
    {
        std::cerr << "No contigs listed in header or header completely missing! Terminate" << std::endl;
        exit(1);
    }
    
    
    BamScannerCache cache;
    String<BamAlignmentRecord> records;
    size_t nextId = 0;
    
    for (size_t r = 0; !atEnd(bamFileIn); ++r)
    {
        appendValue(recordToId, INVALID_ID);
        readMultiRecords(records, bamFileIn, cache);
        if (empty(records))
            break;
        
        if (options.verbosity >= 1 && (r % 1000000) == 0)
        {
            std::cerr << '.' << std::flush;
        }
        
        if (!hasFlagSecondary(records[0]) && (length(records) == 1 || !hasFlagSecondary(records[1])))
        {
            recordToId[r] = nextId++;
        }
        else
        {
            appendValue(secAligns[records[0].qName], r);
            //std::cerr << secAligns.size() << std::endl;
        }
    }
    
    if (options.verbosity >= 1)
        std::cerr << "Done first Parse. Num Primary reads: " << nextId << std::endl;
    
    
    //set ids for secondary alignments
    setPosition(bamFileIn, 0);
    readHeader(header, bamFileIn);
    
    
    for (size_t i = 0; !atEnd(bamFileIn); ++i)
    {
        readMultiRecords(records, bamFileIn, cache);
        if (empty(records))
            break;
        
        if (options.verbosity >= 1 && (i % 1000000) == 0)
        {
            std::cerr << '.' << std::flush;
        }
        
        if (recordToId[i] != INVALID_ID)
        {
            TSecondaryAlignedIter it = secAligns.find(records[0].qName);
            if (it != secAligns.end())
            {
                ++numPri;
                TPosListConstIter lIt = begin(it->second);
                for (; !atEnd(lIt); ++lIt)
                {
                    recordToId[*lIt] = recordToId[i];
                }
                clear(it->second);
            }
        }
    }
    
    if (options.verbosity >= 1)
        std::cerr << "Done second Parse" << std::endl;
    
    if (options.verbosity >= 1)
    {
        size_t numNoId = 0;
        for (size_t i = 0; i < length(recordToId); ++i)
        {
            if (recordToId[i] == INVALID_ID)
                ++numNoId;
        }
        
        std::cerr << "Total missing Ids: " << numNoId << std::endl;
    }
    
    //in case that no alignment is flagged as primary assign ids to leftover secondary reads
    for (TSecondaryAlignedIter it = secAligns.begin(); it !=  secAligns.end(); ++it)
    {
        if (!empty(it->second))
        {
            if (length(it->second)==1)
            {
                ++numSingleSecs;
            }
            TPosListConstIter lIt = begin(it->second);
            for (; !atEnd(lIt); ++lIt)
            {
                recordToId[*lIt] = nextId;
            }
            ++nextId;
        }
    }
    
    if (options.verbosity >= 1)
        std::cerr << "Done third Parse. Num mapped reads: " << nextId << std::endl;
    
    if (options.verbosity >= 1)
    {
        size_t numNoId = 0;
        for (size_t i = 0; i < length(recordToId); ++i)
        {
            if (recordToId[i] == INVALID_ID)
                ++numNoId;
        }
        std::cerr << "Total missing Ids: " << numNoId << "  " << numSingleSecs << "  " << numPri << std::endl;
    }

    return nextId;
}

//returns the last parsed recordId
size_t getIntervalsForContig(TSplitReadIntervalsList &splitReadIntervals,
                             TIntervalList &intervals,
                             TJunctionList &junctions,
                             seqan::String<size_t> &recordToId,
                             TLenDistMap &rLenMap,
                             size_t firstRecordId,
                             unsigned contigId,
                             const Options & options)
{
    typedef FragmentStore<>							              				TFragmentStore;
    typedef TFragmentStore::TAlignedReadStore				TAlignedReadStore;
    typedef TFragmentStore::TContigStore             TContigStoreStore;
    typedef Value<TAlignedReadStore>::Type		   			TAlignedRead;
    typedef Value<TContigStoreStore>::Type				  	TContig;
    typedef TAlignedRead::TGapAnchors						    TReadAnchors;
    typedef TFragmentStore::TContigPos						    TContigPos;
    typedef TContig::TGapAnchors                     TContigAnchors;
    typedef TFragmentStore::TContigSeq               TContigSeq;
    typedef TFragmentStore::TReadSeq                 TReadSeq;
    typedef Gaps<TContigSeq, AnchorGaps<TContigAnchors> >     TContigGaps;
    typedef Gaps<TReadSeq, AnchorGaps<TReadAnchors> >         TReadGaps;
    //typedef TContig::TId                             TContigId;
    //typedef Pair<TContigId, size_t> TPair;
    //typedef std::map<CharString, String<TPair> > TSecondaryAlignedMap;
    //typedef TSecondaryAlignedMap::iterator TSecondaryAlignedIter;
    //typedef Iterator<String<TPair> const, Rooted>::Type TPairListConstIter;
    
    
    //TContigId nextId = storeIter->contigId;
    String<TContigPos> blockBoundariesLeft;
    String<TContigPos> blockBoundariesRight;
    
    TReadSeq readSeqLeft;
    TReadSeq readSeqRight;
    TReadAnchors readGapAnchorsLeft, readGapAnchorsRight;
    TContigSeq contig;
    TContigAnchors contigGapAnchorsLeft, contigGapAnchorsRight;
    
    //TSecondaryAlignedMap secAligns;
    
    unsigned unpaired = 0;
    unsigned skippedShort = 0;
    unsigned skippedGap = 0;
    
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(options.samInFile)))
    {
        std::cerr << "Could not open " << options.samInFile << std::endl;
        return firstRecordId;
    }
    
    BamHeader header;
    readHeader(header, bamFileIn);
    
    BamIOContext<> const & bamContext = context(bamFileIn);
    
    // Read BAI index.
    BamIndex<seqan::Bai> baiIndex;
    if (!open(baiIndex, toCString(options.bamIndexFile)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << options.bamIndexFile << "\n";
        return firstRecordId;
    }
    
    bool hasAlignments = false;
    if (!jumpToRegion(bamFileIn, hasAlignments, contigId, 0, contigLengths(bamContext)[contigId], baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << 0 << ":" << contigLengths(bamContext)[contigId] << "\n";
        return firstRecordId;
    }
    if (!hasAlignments) // No alignments here.
    {
        std::cerr << "No alignments for contig: " << contigNames(bamContext)[contigId] << std::endl;
        return firstRecordId;
    }
    
    
    BamScannerCache cache;
    String<BamAlignmentRecord> records;
    
    size_t maxContigId = length(contigNames(context(bamFileIn)));
    if (maxContigId == 0)
    {
        std::cerr << "No contigs listed in header or header completely missing! Terminate" << std::endl;
        exit(1);
    }
    
    for (size_t &r = firstRecordId; !atEnd(bamFileIn); ++r)
    {
        appendValue(recordToId, INVALID_ID);
        readMultiRecords(records, bamFileIn, cache);
        
        if (empty(records))
            break;
        
        if (records[0].rID != static_cast<int>(contigId)) //reached next contig
            break;
        
        if (options.verbosity >= 1 && (r % 1000000) == 0)
        {
            std::cerr << '.' << std::flush;
        }
        
        TContigPos beginPos, endPos;
        
        clear(readGapAnchorsLeft);
        unsigned queryLengthL =  _getQueryLength(records[0].cigar);
        ++rLenMap[queryLengthL]; // update for read length statistics
        readSeqLeft = records[0].seq;        
        if (empty(readSeqLeft))
            resize(readSeqLeft, queryLengthL, 'N');
        
        TReadGaps readGapsLeft(readSeqLeft, readGapAnchorsLeft);
        cigarToGapAnchorRead(readGapsLeft, records[0].cigar);
        
        beginPos = records[0].beginPos;
        endPos = beginPos + length(readGapsLeft);
        
        if (beginPos > endPos)
        {
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }
        
        clear(contigGapAnchorsLeft);
        TContigGaps	contigGapsLeft(contig, contigGapAnchorsLeft);
        cigarToGapAnchorContig(contigGapsLeft, records[0].cigar);
        
        clear(blockBoundariesLeft);
        clear(blockBoundariesRight);
        
        // we only need the correct lengths not the sequence itself
        if ((TContigPos)length(contig) < endPos)
            resize(contig, endPos);
        
        setBeginPosition(contigGapsLeft, beginPos);
        setEndPosition(contigGapsLeft, endPos);
        
        // get gap boundaries
        //getBlockBoundaries(blockBoundariesLeft, contigGapsLeft, readGapsLeft);
        getBlockBoundaries(blockBoundariesLeft, records[0].cigar, beginPos);
        
        bool matchTooShort = false;
        bool gapTooLong = false;
        for (unsigned i = 0; i + 1 < length(blockBoundariesLeft); i+=2)
        {
            if (blockBoundariesLeft[i+1] - blockBoundariesLeft[i] < options.minMatchLength)
            {
                matchTooShort = true;
                break;
            }
            
            if (i != 0 && blockBoundariesLeft[i] - blockBoundariesLeft[i-1] > options.maxGapLength)
            {
                gapTooLong = true;
                break;
            }
        }
        if (matchTooShort)
        {
            ++skippedShort;
            continue;
        }
        
        if (gapTooLong)
        {
            ++skippedGap;
            continue;
        }

        
        // remove small gaps and result in exon blocks
        //removeSmallGaps(blockBoundariesLeft, options);
        
        TContigPos beginPosRight = beginPos;
        TContigPos endPosRight = endPos;
        
        if (length(records) > 1)
        {
            {
                clear(readGapAnchorsRight);
                readSeqRight = records[1].seq;
                unsigned queryLengthR =  _getQueryLength(records[1].cigar);
                ++rLenMap[queryLengthR]; // update for read length statistics
                if (empty(readSeqRight))
                    resize(readSeqRight, queryLengthR, 'N');
                
                TReadGaps readGapsRight(readSeqRight, readGapAnchorsRight);
                cigarToGapAnchorRead(readGapsRight, records[1].cigar);
                
                beginPosRight = records[1].beginPos;
                endPosRight = beginPosRight + length(readGapsRight);
                
                if (beginPosRight > endPosRight)
                {
                    TContigPos temp = beginPosRight;
                    beginPosRight = endPosRight;
                    endPosRight = temp;
                }
                
                if (beginPosRight - endPos > options.maxGapLength)
                {
                    ++skippedGap;
                    continue;
                }

                clear(contigGapAnchorsRight);
                TContigGaps	contigGapsRight(contig, contigGapAnchorsRight);
                cigarToGapAnchorContig(contigGapsRight, records[1].cigar);
                
                if ((TContigPos)length(contig) < endPosRight)
                    resize(contig, endPosRight);
                
                setBeginPosition(contigGapsRight, beginPosRight);
                setEndPosition(contigGapsRight, endPosRight);
                
                // get gap boundaries
                //getBlockBoundaries(blockBoundariesRight, contigGapsRight, readGapsRight);
                getBlockBoundaries(blockBoundariesRight, records[1].cigar, beginPosRight);
                
                for (unsigned i = 0; i + 1 < length(blockBoundariesRight); i+=2)
                {
                    if (blockBoundariesRight[i+1] - blockBoundariesRight[i] < options.minMatchLength)
                    {
                        matchTooShort = true;
                        break;
                    }
                    
                    if (i != 0 && blockBoundariesRight[i] - blockBoundariesRight[i-1] > options.maxGapLength)
                    {
                        gapTooLong = true;
                        break;
                    }
                }
                
                if (matchTooShort)
                {
                    ++skippedShort;
                    continue;
                }
                
                if (gapTooLong)
                {
                    ++skippedGap;
                    continue;
                }
                // remove small gaps and result in exon blocks
                //removeSmallGaps(blockBoundariesRight, options);
            }
        }
        
        if (options.rMode == PAIRED && empty(blockBoundariesRight))
        {
            if (options.verbosity >= 2)
                std::cerr << "skipping unpaired read\t" << beginPos<<'\t'<<endPos<<'\t'<<beginPosRight<<'\t'<<endPosRight<<std::endl;
            ++unpaired;
            continue;
        }
        
        resize(splitReadIntervals, length(splitReadIntervals) + 1 /*, seqan::Generous()*/);
        TSplitReadIntervals & tmp = back(splitReadIntervals);
        tmp.readId = recordToId[r];
        
        // check for tophat strand info
        if (options.useStrandInfo)
        {
            tmp.strand = '0'; // default: unknown
            seqan::BamTagsDict tagsDict(records[0].tags);
            unsigned tagIdx = 0;

            if (findTagKey(tagIdx, tagsDict, "XS"))
            {
                char strandTag;
                if (extractTagValue(strandTag, tagsDict, tagIdx))
                {
                    tmp.strand = strandTag;
                    //std::cout << records[0].qName << " " << tmp.strand << std::endl;
                }
            }
        }
        
        reserve(tmp.intervalsLeft, length(blockBoundariesLeft)/2);
        for (unsigned i = 0; i + 1 < length(blockBoundariesLeft); i+=2)
        {
            //appendValue(intervals[contigId], TInterval(blockBoundariesLeft[i], blockBoundariesLeft[i+1], contigId));
            appendValue(intervals, TInterval(blockBoundariesLeft[i], blockBoundariesLeft[i+1], -1));
            appendValue(tmp.intervalsLeft, back(intervals));
            
            if (i > 0)
            {
                appendValue(junctions, std::make_pair(blockBoundariesLeft[i-1], blockBoundariesLeft[i]));
            }
        }
        shrinkToFit(tmp.intervalsLeft);
        //std::sort(begin(tmp.intervalsLeft), end(tmp.intervalsLeft));
        
        reserve(tmp.intervalsRight, length(blockBoundariesRight)/2);
        for (unsigned i = 0; i + 1 < length(blockBoundariesRight); i+=2)
        {
            //appendValue(intervals[contigId], TInterval(blockBoundariesRight[i], blockBoundariesRight[i+1], contigId));
            appendValue(intervals, TInterval(blockBoundariesRight[i], blockBoundariesRight[i+1], -1));
            appendValue(tmp.intervalsRight, back(intervals));
            if (i > 0)
            {
                appendValue(junctions, std::make_pair(blockBoundariesRight[i-1], blockBoundariesRight[i]));
            }
        }
        shrinkToFit(tmp.intervalsRight);
        //std::sort(begin(tmp.intervalsRight), end(tmp.intervalsRight));
    }
    
    //make junctions and interval a set
    std::sort(begin(junctions), end(junctions));
    TJunctionListIterator newEndIt = std::unique(begin(junctions), end(junctions));
    resize(junctions, newEndIt - begin(junctions));
    
    std::sort(begin(intervals), end(intervals));
    TIntervalListIterator newEndIt2 = std::unique(begin(intervals), end(intervals));
    resize(intervals, newEndIt2 - begin(intervals));
    
    
    if (options.verbosity >= 1)
        std::cerr << "Done set creation" << std::endl;
    
    if (unpaired && options.verbosity >= 1)
        std::cerr << "skipped " << unpaired << " unpaired reads" << std::endl;
    
    if (skippedShort && options.verbosity >= 1)
        std::cerr << "skipped " << skippedShort << " reads with too short matching segment" << std::endl;
    
    if (skippedGap && options.verbosity >= 1)
        std::cerr << "skipped " << skippedGap << " reads with too long gaps between split parts or read pairs" << std::endl;
    
    
    return firstRecordId;
}


//reads the intervals (corresponding to subexons) from refined.gff file
template <typename TIntervalsMap, typename TannoNameStoreInv>
void readIntervalsFromAnnotation(GffFileIn &inFile, TIntervalsMap &intervals, TannoNameStoreInv& annoNameStoreInv, std::map<size_t,std::vector<size_t> > & groupGenes, const Options & options)
{
    typedef seqan::FragmentStore<> TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename seqan::Iterator<TFragmentStore, seqan::AnnotationTree<> >::Type TAnnoTreeIter;
    typedef typename seqan::Value<TAnnotationStore>::Type TAnnotation;
    
    TFragmentStore store;
    unsigned SubExonType = 0;
    _storeAppendType(store, SubExonType, "subexon");
    unsigned geneType = 0;
    _storeAppendType(store, geneType, "gene");
    
    if (options.verbosity >= 1)
        std::cerr << "Loading GFF" << std::endl;
    
    readRecords(store, inFile);
    
    if (options.verbosity >= 1)
        std::cerr << "Loading Done" << std::endl;
    
    TAnnoTreeIter dfsIt = seqan::begin(store, seqan::AnnotationTree<>());
    
    if (options.verbosity >= 1)
        std::cerr << "storing subexons in InterValList" << std::endl;
    
    String<bool> abundantIds;
    
    size_t lastGene;
    while (!atEnd(dfsIt))
    {
        TAnnotation &anno = getAnnotation(dfsIt);
        if (anno.typeId == geneType)
            lastGene = dfsIt._id;
        
        if (anno.typeId == SubExonType)
        {
            seqan::CharString nodeIdStr;
            annotationGetValueByKey(nodeIdStr, store, anno, "NodeId");
            std::stringstream strstr(toCString(nodeIdStr));
            unsigned subexonId;
            strstr >> subexonId;
            groupGenes[lastGene].push_back(subexonId);
            if (subexonId >= length(abundantIds))
            {
                resize(abundantIds, subexonId + 1, false);
            }
            if (!abundantIds[subexonId])
            {
                abundantIds[subexonId] = true;
                if (anno.contigId >= length(intervals))
                    resize(intervals, anno.contigId + 1);
                
                appendValue(intervals[anno.contigId], TInterval(anno.beginPos, anno.endPos, subexonId));
            }
        }
        ++dfsIt;
    }
    
    for (size_t i = 0; i < length(intervals); ++i)//mapIt = intervals.begin(); mapIt != intervals.end(); ++mapIt)
    {
        std::sort(begin(intervals[i]), end(intervals[i]));
        if (options.verbosity >= 1)
            std::cerr << "For Chromosome " << store.contigNameStore[i] << "  found "  << length(intervals[i]) << " subexons" << std::endl;
    }
    
    for(std::map<size_t,std::vector<size_t> >::iterator it = groupGenes.begin(); it != groupGenes.end(); ++it)
    {
        std::sort(it->second.begin(), it->second.end());
        it->second.erase(std::unique(it->second.begin(), it->second.end()), it->second.end());
    }
    
    for (size_t i = 0; i < length(store.contigNameStore); ++i)
        annoNameStoreInv[store.contigNameStore[i]] = i;
}


struct IndexComparator
{
    IndexComparator(const TIntervalList & ptr, const String<unsigned> & inv): intervals(ptr), inverse(inv){}
    
    bool operator() (TIntervalPosType i, TIntervalPosType j)
    {
        return (intervals[inverse[i]].i1 < intervals[inverse[j]].i1);
    }
    const TIntervalList& intervals;
    const String<unsigned>& inverse;
};


//reads are transformed into our count format [exons_left_mate]^[exons_right_mate]
void intervalIndexTransform(const TIntervalList &subexons,
                            const TSplitReadIntervalsList &readsInter,
                            TSplitReadIndicesList &readsIndex,
                            TLenDistMap &innerDistMap,
                            const Options &options)
{
    //create one interval tree for each strand
    String<IntervalTree<TIntervalPosType, TIntervalCargoType> > intervalTrees;
    TIntervalList subexonsTmp;
    TIntervalListConstIterator intervalIt = begin(subexons);
    TIntervalListConstIterator intervalEnd = end(subexons);
    TIntervalList intervalsForward, intervalsReverse;
    
    unsigned boundaryTolerance = 0;
    
    enum ErrorCode
    {
        NOERROR                  = 0,
        NOTCONSECUTIVE1          = 1,
        OVERHANGING1             = 2,
        UNSUPPORTEDEXONBOUNDARY1 = 3,
        UNSUPPORTEDJUNCTION1     = 4,
        NOINTERVALSFOUND1        = 5,
        NOTCONSECUTIVE2          = 6,
        OVERHANGING2             = 7,
        UNSUPPORTEDEXONBOUNDARY2 = 8,
        UNSUPPORTEDJUNCTION2     = 9,
        NOINTERVALSFOUND2        = 10,
        INCONSISTENTPAIRS        = 11,
        WRONGSTRAND              = 12
    };
    
    unsigned doubleAligned = 0, unaligned = 0, invalidSplit = 0;
    
    while (intervalIt != intervalEnd)
    {
        if (intervalIt->i1 >= intervalIt->i2)
        {
            appendValue(intervalsReverse, TInterval(intervalIt->i2, intervalIt->i1, intervalIt->cargo));
            appendValue(subexonsTmp, TInterval(intervalIt->i2, intervalIt->i1, intervalIt->cargo));
        }
        else
        {
            appendValue(intervalsForward, *intervalIt);
            appendValue(subexonsTmp, *intervalIt);
        }
        ++intervalIt;
    }
    
    if (options.verbosity >= 2)
        std::cerr << "size of Reverse Interval Tree: " << length(intervalsReverse) << std:: endl;
    
    appendValue(intervalTrees, IntervalTree<TIntervalPosType, TIntervalCargoType> (intervalsForward));
    appendValue(intervalTrees, IntervalTree<TIntervalPosType, TIntervalCargoType> (intervalsReverse));
    reserve(readsIndex, length(readsIndex) + length(readsInter));
    
    //compute inverse intervals to access the right intervals with the cargos returned from findIntervals
    String<unsigned> intervalsInv;
    for (unsigned i = 0; i < length(subexons); ++i)
    {
        if (subexons[i].cargo >= length(intervalsInv))
        {
            resize(intervalsInv, subexons[i].cargo + 1, -1);
        }
        intervalsInv[subexons[i].cargo] = i;
    }
    
    TSplitReadIntervalsListConstIterator readsIt = begin(readsInter);
    String<TIntervalCargoType> results;
    
    IndexComparator comp(subexons, intervalsInv);
    
    while (readsIt != end(readsInter))
    {
        //int validRead = 0;
        ErrorCode error[2] = {NOERROR, NOERROR};
        for (unsigned strand = 0; strand < 2; ++strand)
        {
            if (options.useStrandInfo)
            {
                if ((strand == 0 && readsIt->strand == '-') ||
                    (strand == 1 && readsIt->strand == '+'))
                {
                    error[strand] = WRONGSTRAND;
                    continue;
                }
            }
            
            TIntervalListConstIterator interItL = begin(readsIt->intervalsLeft);
            TSplitReadIndices readTmp;
            readTmp.readId = readsIt->readId;
            //bool validForStrand = true;
            if(options.rMode == PAIRED && options.verbosity >= 2)
            {
                if(empty(readsIt->intervalsRight))
                    std::cerr << "No right intervals:  " << readsIt->readId << std::endl;
                
                if(empty(readsIt->intervalsLeft))
                    std::cerr << "No left intervals:  " << readsIt->readId << std::endl;
            }
            
            size_t numInter = 0;
            while (interItL != end(readsIt->intervalsLeft) )
            {
                clear(results);
                findIntervals(results, intervalTrees[strand], interItL->i1, interItL->i2);
                
                if (empty(results))
                {
                    error[strand] = NOINTERVALSFOUND1;
                    break;
                }
                std::sort(begin(results), end(results), comp);
                
                for (unsigned i = 0; i + 1 < length(results); ++i)
                {
                    if ( (subexonsTmp[intervalsInv[results[i]]].i2 != subexonsTmp[intervalsInv[results[i+1]]].i1))
                    {
                        error[strand] = NOTCONSECUTIVE1;
                        break;
                    }
                }
                if (error[strand] != NOERROR)
                {
                    break;
                }
                
                if (subexonsTmp[intervalsInv[results[0]]].i1 > interItL->i1 + boundaryTolerance || subexonsTmp[intervalsInv[back(results)]].i2 + boundaryTolerance < interItL->i2)
                {
                    error[strand] = OVERHANGING1;
                    break;
                }
                if (numInter != 0 && (_max(subexonsTmp[intervalsInv[results[0]]].i1, interItL->i1) - _min(subexonsTmp[intervalsInv[results[0]]].i1, interItL->i1)) > boundaryTolerance)
                {
                    //check whether left split position of the read matches the begin position of the first subexon
                    //std::cerr << subexons[intervalsInv[results[0]]].i1 << " " << subexons[intervalsInv[results[0]]].i2 << " L " << interItL->i1 << " " << interItL->i2 << std::endl;
                    error[strand] = UNSUPPORTEDEXONBOUNDARY1;
                    break;
                }
                if (numInter + 1 != length(readsIt->intervalsLeft) && (_max(subexonsTmp[intervalsInv[back(results)]].i2,interItL->i2)-_min(subexonsTmp[intervalsInv[back(results)]].i2,interItL->i2)) > boundaryTolerance)
                {
                    //check whether right split position of the read matches the end position of the last subexon
                    //std::cerr << subexons[intervalsInv[back(results)]].i1 << " " << subexons[intervalsInv[back(results)]].i2 << " R " << interItL->i1 << " " << interItL->i2 << std::endl;
                    error[strand] = UNSUPPORTEDEXONBOUNDARY1;
                    break;
                }
                append(readTmp.indicesLeft, results);
                ++interItL;
                ++numInter;
            }
            
            if (error[strand] == NOERROR)
            {
                numInter = 0;
                TIntervalListConstIterator interItR = begin(readsIt->intervalsRight);
                while (interItR != end(readsIt->intervalsRight) )
                {
                    clear(results);
                    findIntervals(results, intervalTrees[strand], interItR->i1, interItR->i2);
                    if (empty(results))
                    {
                        error[strand] = NOINTERVALSFOUND2;
                        break;
                    }
                    
                    std::sort(begin(results), end(results), comp);
                    
                    for (unsigned i = 0; i + 1 < length(results); ++i)
                    {
                        if ( (subexonsTmp[intervalsInv[results[i]]].i2 != subexonsTmp[intervalsInv[results[i+1]]].i1))
                        {
                            error[strand] = NOTCONSECUTIVE2;
                            break;
                        }
                    }
                    if (error[strand] != NOERROR)
                    {
                        break;
                    }
                    if (subexonsTmp[intervalsInv[results[0]]].i1 > interItR->i1 + boundaryTolerance || subexonsTmp[intervalsInv[back(results)]].i2 + boundaryTolerance < interItR->i2)
                    {
                        error[strand] = OVERHANGING2;
                        break;
                    }
                    if (numInter != 0 && (_max(subexonsTmp[intervalsInv[results[0]]].i1, interItR->i1)-_min(subexonsTmp[intervalsInv[results[0]]].i1, interItR->i1)) > boundaryTolerance)
                    {
                        //check whether left split position of the read matches the begin position of the subexon
                        error[strand] = UNSUPPORTEDEXONBOUNDARY2;
                        break;
                    }
                    if (numInter +1 != length(readsIt->intervalsRight) && (_max(subexonsTmp[intervalsInv[back(results)]].i2,interItR->i2) - _min(subexonsTmp[intervalsInv[back(results)]].i2,interItR->i2)) > boundaryTolerance)
                    {
                        //check whether right split position of the read matches the end position of the subexon
                        error[strand] = UNSUPPORTEDEXONBOUNDARY2;
                        break;
                    }
                    
                    append(readTmp.indicesRight, results);
                    ++interItR;
                    ++numInter;
                }
            }
            
            if (error[strand] == NOERROR)
            {
                //                for (size_t kk = 0; kk < length(readTmp.indicesLeft); ++kk)
                //                    std::cout << readTmp.indicesLeft[kk] << '\t';
                //                std::cout << " ^ ";
                //                for (size_t kk = 0; kk < length(readTmp.indicesRight); ++kk)
                //                    std::cout << readTmp.indicesRight[kk] << '\t';
                //                std::cout << std::endl;
                
#ifdef GENERATEGAPS
                if (!empty(readTmp.indicesLeft) && !empty(readTmp.indicesRight))
                {
                    //compute the lower bound on the insert size for the given mapping
                    int minIs = 0;
                    
                    if (back(readTmp.indicesLeft) == front(readTmp.indicesRight) ||
                        front(readsIt->intervalsRight).i1 <= back(readsIt->intervalsLeft).i2)
                    {
                        //pairs that cover the same subexon or even overlapping pairs
                        minIs = static_cast<int>(front(readsIt->intervalsRight).i1) - static_cast<int>(back(readsIt->intervalsLeft).i2);
                    }
                    else
                    {
                        //pairs without shared subexons
                        minIs = (subexonsTmp[intervalsInv[back(readTmp.indicesLeft)]].i2  - back(readsIt->intervalsLeft).i2) +
                        (front(readsIt->intervalsRight).i1 - subexonsTmp[intervalsInv[front(readTmp.indicesRight)]].i1);
                        //std::cout << "Universal: " << minIs << std::endl;
                    }
                    readTmp.minGap = minIs;
                }
#endif
                //infer fragment length if rightmost segment of left read equals leftmost segment of right read or they are neighbors
                if (!empty(readTmp.indicesRight))
                {
                    if (back(readTmp.indicesLeft) == front(readTmp.indicesRight)) //same segment
                    {
                        ++innerDistMap[static_cast<int>(front(readsIt->intervalsRight).i1) - static_cast<int>(back(readsIt->intervalsLeft).i2)];
                    }
                    else if (1 + back(readTmp.indicesLeft) == front(readTmp.indicesRight)) //neighboring segments
                    {
                        int gapLen = (static_cast<int>(subexonsTmp[intervalsInv[back(readTmp.indicesLeft)]].i2)  - static_cast<int>(back(readsIt->intervalsLeft).i2)) +
                                (static_cast<int>(front(readsIt->intervalsRight).i1) - static_cast<int>(subexonsTmp[intervalsInv[front(readTmp.indicesRight)]].i1));
                        ++innerDistMap[gapLen];
                    }
                }

                if (strand == 1)
                {
                    //in transcript orientation
                    std::reverse(begin(readTmp.indicesLeft),end(readTmp.indicesLeft));
                    
                    if (!empty(readTmp.indicesRight))
                    {
                        std::reverse(begin(readTmp.indicesRight),end(readTmp.indicesRight));
                        swap(readTmp.indicesLeft, readTmp.indicesRight);
                    }
                }
                
                if (std::adjacent_find(begin(readTmp.indicesLeft), end(readTmp.indicesLeft)) != end(readTmp.indicesLeft))
                {
                    //if this holds this read contains a junction not supported by our subexons
                    error[strand] = UNSUPPORTEDJUNCTION1;
                    ++invalidSplit;
                }
                else if (std::adjacent_find(begin(readTmp.indicesRight), end(readTmp.indicesRight)) != end(readTmp.indicesRight))
                {
                    //if this holds this read contains a junction not supported by our subexons
                    error[strand] = UNSUPPORTEDJUNCTION2;
                    ++invalidSplit;
                }
                //check for inconsistent junction in mates
                else if (!empty(readTmp.indicesRight))
                {
                    size_t posL(0), posR(1);
                    while (posL < length(readTmp.indicesLeft) && readTmp.indicesLeft[posL] < readTmp.indicesRight[0])
                        ++posL;
                    
                    if (posL < length(readTmp.indicesLeft) && readTmp.indicesLeft[posL] == readTmp.indicesRight[0])
                    {
                        ++posL;
                        bool valid = true;
                        while (valid && posL < length(readTmp.indicesLeft) && posR < length(readTmp.indicesRight))
                        {
                            valid = readTmp.indicesLeft[posL++] == readTmp.indicesRight[posR++];
                        }
                        
                        if (!valid)
                        {
                            error[strand] = INCONSISTENTPAIRS;
                            if (options.verbosity >= 2)
                            {
                                std::cerr << "inconsistent Pair! " << readsIt->readId << std::endl;
                                for (size_t p = 0; p < length(readTmp.indicesLeft); ++ p)
                                    std::cerr << readTmp.indicesLeft[p] << " - ";
                                std::cerr << " ^ ";
                                for (size_t p = 0; p < length(readTmp.indicesRight); ++ p)
                                    std::cerr << readTmp.indicesRight[p] << " - ";
                                std::cerr << std::endl;
                                for (size_t p = 0; p < length(readsIt->intervalsLeft); ++ p)
                                    std::cerr << "[" << readsIt->intervalsLeft[p].i1 << "," << readsIt->intervalsLeft[p].i2 << "] ";
                                std::cerr << " ^ ";
                                for (size_t p = 0; p < length(readsIt->intervalsRight); ++ p)
                                    std::cerr << "[" << readsIt->intervalsRight[p].i1 << "," << readsIt->intervalsRight[p].i2 << "] ";
                                std::cerr << std::endl;
                            }
                        }
                    }
                }
            }
            
            if (error[strand] == NOERROR)
            {
                appendValue(readsIndex, readTmp);
            }
        }//strand loop
        
        if (error[0] != NOERROR && error[1] != NOERROR)
        {
            ++unaligned;
            if (options.verbosity >= 2)
            {
                std::cerr << error[0] << "  ::  " << error[1] << '\t';
                
                for (size_t p = 0; p < length(readsIt->intervalsLeft); ++ p)
                    std::cerr << "[" << readsIt->intervalsLeft[p].i1 << "," << readsIt->intervalsLeft[p].i2 << "] ";
                std::cerr << " ^ ";
                for (size_t p = 0; p < length(readsIt->intervalsRight); ++ p)
                    std::cerr << "[" << readsIt->intervalsRight[p].i1 << "," << readsIt->intervalsRight[p].i2 << "] ";
                
                std::cerr << std::endl;
            }
        }
        if (error[0] == NOERROR && error[1] == NOERROR)
        {
            ++doubleAligned;
        }
        ++readsIt;
    }
    
    if (options.verbosity >= 1)
    {
        std::cerr << "Not Aligned: " << unaligned << " --  Two strand aligned: ";
        std::cerr << doubleAligned << " -- unsupportedSplit: " << invalidSplit << std::endl;
    }
}



void processIntervals(const TIntervalList &intervals,
                      const TJunctionList &junctions,
                      const Options & options,
                      TIntervalList &subexonIntervals,
                      size_t &nextIndex)
{
    if (empty(intervals))
    {
        return;
    }
    
    TIntervalListConstIterator itR = begin(intervals), itEnd = end(intervals);
    TIntervalListConstIterator itL = itR++;
    TIntervalList mergedIntervals;
    TSplitPositionList donor, acceptor, splits;
    
    
    // from the junctions get a list of sorted split positions
    TJunctionListIteratorConst juncsIt = begin(junctions);
    while (juncsIt != end(junctions))
    {
        splits.insert(juncsIt->first);
        splits.insert(juncsIt->second);
        donor.insert(juncsIt->first);
        acceptor.insert(juncsIt->second);
        ++juncsIt;
    }
    
    
    appendValue(mergedIntervals, *itL);
    ++itL;
    ++itR;
    // MERGING OVERLAPPING INTERVALS
    while (itL != itEnd)
    {
        bool newInter = (itL->i1 > back(mergedIntervals).i2 + options.minIntronLength);
        
        newInter = newInter || (itL->i1 > back(mergedIntervals).i2 &&
                                donor.count(back(mergedIntervals).i2) &&
                                acceptor.count(itL->i1));
        
        //if ( (itL->i1 <= back(mergedIntervals).i2 + options.minIntronLength)/* &&
        //    (itR == itEnd || ((itL->cargo != BOTH && itL->cargo!= RIGHT) && (itR->cargo != BOTH && itR->cargo != LEFT)))*/)
        if (!newInter) // merge with last interval
        {
            if (itL->i2 > back(mergedIntervals).i2)
            {
                // extend stored interval
                back(mergedIntervals).i2 = itL->i2;
            }
        }
        else //open new interval
        {
            appendValue(mergedIntervals, *itL);
        }
        ++itL;
        ++itR;
    }
    
    // SPLIT INTERVALS ACCORDING TO SPLIT READ ALIGNMENTS
    TSplitPositionListConstIter splitIter = splits.begin();
    TIntervalListIterator mergedIt = begin(mergedIntervals);
    while (splitIter != splits.end() && mergedIt != end(mergedIntervals) )
    {
        while (mergedIt != end(mergedIntervals) && *splitIter >= mergedIt->i2)
        {
            appendValue(subexonIntervals, *mergedIt);
            ++mergedIt;
        }
        
        if (*splitIter > mergedIt->i1)
        {
            appendValue(subexonIntervals, TInterval(mergedIt->i1, *splitIter, -1));
            mergedIt->i1 = *splitIter;
        }
        ++splitIter;
    }
    
    // add possible remaining intervals at the end and assign ids
    append(subexonIntervals, suffix(mergedIntervals, mergedIt) );
    
    TIntervalListIterator subIt = begin(subexonIntervals);
    while (subIt != end(subexonIntervals))
    {
        subIt->cargo = nextIndex++;
        //        bool isDonor = donor.count(subIt->i2);
        //        bool isAcceptor = donor.count(subIt->i2);
        //
        //        if (isDonor && isAcceptor)
        //            subIt->cargo = FULL;
        //        else if (isDonor)
        //            subIt->cargo = RIGHT;
        //        else if (isAcceptor)
        //            subIt->cargo = LEFT;
        //        subIt->cargo = FULL;
        ++subIt;
    }
}


void generateCounts(const String<TSplitReadIndicesList>& splitReads,
                    String<unsigned> & readCounts,
                    TTupleCountList & counts,
                    const Options & options)
{
    size_t unaligned = 0;
    for (size_t i = 0; i < length(splitReads); ++i)
    {
        TSplitReadIndicesListConstIterator readsIt = begin(splitReads[i]);
        while (readsIt != end(splitReads[i]))
        {
            if(readCounts[readsIt->readId] > 0)
            {
                if(!options.skipMulti || readCounts[readsIt->readId] == 1)
                {
                    counts[TCountTuple(readsIt->indicesLeft, readsIt->indicesRight)].totalCount += 1. / readCounts[readsIt->readId];
                }
                if (readCounts[readsIt->readId] == 1)
                {
                    counts[TCountTuple(readsIt->indicesLeft, readsIt->indicesRight)].uniqueCount += 1;
                }
            }
            else
            {
                ++unaligned;
            }
            ++readsIt;
        }
    }
    if (options.verbosity >= 1)
        std::cerr<<"Total number unaligned reads: " << unaligned <<std::endl;
}

void generateCounts(const String<TSplitReadIndicesList>& splitReads,
                    String<unsigned> & readCounts,
                    TTupleCountList & counts,
                    TTupleCountGapList & countsGaps,
                    const Options & options)
{
    size_t unaligned = 0;
    for (size_t i = 0; i < length(splitReads); ++i)
    {
        TSplitReadIndicesListConstIterator readsIt = begin(splitReads[i]);
        while (readsIt != end(splitReads[i]))
        {
            if(readCounts[readsIt->readId] > 0)
            {
                if(!options.skipMulti || readCounts[readsIt->readId] == 1)
                {
                    counts[TCountTuple(readsIt->indicesLeft, readsIt->indicesRight)].totalCount += 1. / readCounts[readsIt->readId];
                    countsGaps[TCountTuple(readsIt->indicesLeft, readsIt->indicesRight)][readsIt->minGap] += 1. / readCounts[readsIt->readId];
                    //std::cout << readsIt->minGap << std::endl;
                }
                if (readCounts[readsIt->readId] == 1)
                {
                    counts[TCountTuple(readsIt->indicesLeft, readsIt->indicesRight)].uniqueCount += 1;
                }
            }
            else
            {
                ++unaligned;
            }
            
            ++readsIt;
        }
    }
    if (options.verbosity >= 1)
        std::cerr<<"Total number unaligned reads: " << unaligned <<std::endl;
}

template<typename TCountFile>
void createCountFileHeader(TCountFile &countFile, const LibInfo &libInfo)
{
    countFile << "#ReadLen:\t" << floor(0.5 + libInfo.rLenMean) << "\n";
    countFile << "#FLen:\t" << floor(0.5 + libInfo.fLenMean) << "\n";
    countFile << "#FLenSdev:\t" << libInfo.fLenSdev << "\n";
    countFile << "#M:\t" << libInfo.totalFrags/1000000. << "\n";
}

template<typename TCountFile>
void writeCountFile(TCountFile &countFile,
                    const TTupleCountList &counts,
                    const Options &options,
                    const LibInfo &libInfo)
{
    //add header to file
    createCountFileHeader(countFile, libInfo);

    TTupleCountListConstIter countIter = counts.begin();
    size_t tmp = 0;
    while (countIter != counts.end())
    {
        TIndexListConstIterator tupleIt = begin(countIter->first.first);
        while (tupleIt != end(countIter->first.first))
        {
            countFile<< *tupleIt;
            ++tupleIt;
            
            if (tupleIt != end(countIter->first.first))
                countFile<<'-';
        }
        
        if (!empty(countIter->first.second))
        {
            countFile<<'^';
            tupleIt = begin(countIter->first.second);
            while (tupleIt != end(countIter->first.second))
            {
                countFile<< *tupleIt;
                ++tupleIt;
                
                if (tupleIt != end(countIter->first.second))
                    countFile<<'-';
            }
        }
        countFile << '\t' << countIter->second.totalCount << '\t' << countIter->second.uniqueCount << "\n";
        ++countIter;
        ++tmp;
    }
    if (options.verbosity >= 1)
        std::cerr << "written " << tmp << " cnt entries" << std::endl;
}

template<typename TCountFile>
void writeCountFile(TCountFile &countFile,
                    const TTupleCountList &counts,
                    TTupleCountGapList &countsGaps,
                    const Options &options,
                    const LibInfo &libInfo)
{        
    //add header to file
    createCountFileHeader(countFile, libInfo);

    TTupleCountListConstIter countIter = counts.begin();
    size_t tmp = 0;
    while (countIter != counts.end())
    {
        TIndexListConstIterator tupleIt = begin(countIter->first.first);
        while (tupleIt != end(countIter->first.first))
        {
            countFile<< *tupleIt;
            ++tupleIt;
            
            if (tupleIt != end(countIter->first.first))
                countFile<<'-';
        }
        
        if (!empty(countIter->first.second))
        {
            countFile<<'^';
            tupleIt = begin(countIter->first.second);
            while (tupleIt != end(countIter->first.second))
            {
                countFile<< *tupleIt;
                ++tupleIt;
                
                if (tupleIt != end(countIter->first.second))
                    countFile<<'-';
            }
        }
        countFile << '\t' << countIter->second.totalCount << '\t' << countIter->second.uniqueCount;
        
        if (countsGaps.count(countIter->first) == 0)
        {
            std::cerr << "Missing gap counts -- aborting! " << std::endl;
            throw;
        }
        
        TGapCountMapConstIter gapCountIt = countsGaps.find(countIter->first)->second.begin();
        TGapCountMapConstIter gapCountEnd = countsGaps.find(countIter->first)->second.end();
        while (gapCountIt != gapCountEnd)
        {
            countFile << '\t' << gapCountIt->first << '\t' << gapCountIt->second;
            ++gapCountIt;
        }
        countFile << "\n";
        
        ++countIter;
        ++tmp;
    }
    if (options.verbosity >= 1)
        std::cerr << "written " << tmp << " cnt entries" << std::endl;
}





template<typename TMatFile, typename TStatFile, typename TFragmentStore>
unsigned writeMatStat(TMatFile & matFile,
                      TStatFile & statFile,
                      TFragmentStore & store,
                      const String<TIntervalList>& subExonsMap,
                      const String<TSplitReadIndicesList>& splitReadsIndices,
                      const Options & options,
                      const std::map<size_t, std::vector<size_t> > & byGene = std::map<size_t, std::vector<size_t> >(),
                      GroupMode group = GROUPBOTH)
{
    typedef Graph<Undirected<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename seqan::Iterator<TFragmentStore, seqan::AnnotationTree<> >::Type TAnnoTreeIter;
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename seqan::Value<TAnnotationStore>::Type TAnnotation;
    
    //----------------------------------------------------------------
    //  BUILD GRAPH AND COMPUTE CONNECTED COMPONENTS
    //----------------------------------------------------------------
    
    //annotation store for output gtf file
    unsigned SubTypeId = 0;
    seqan::_storeAppendType(store, SubTypeId, "subexon");
    TAnnoTreeIter root = seqan::begin(store, seqan::AnnotationTree<>());
    
    unsigned nextGeneId = 1;
    
    for (size_t contigId = 0; contigId < length(subExonsMap); ++ contigId)
    {
        if (empty(subExonsMap[contigId]))
            continue;
        
        TGraph g;
        
        //inverse intervals to access the right intervals with the cargos returned from findIntervals
        String<unsigned> intervalsInv;
        
        const TIntervalList& subExons = subExonsMap[contigId];
        
        for (size_t i = 0; i < length(subExons); ++i)
        {
            addVertex(g);
            
            if (subExons[i].cargo >= length(intervalsInv))
            {
                resize(intervalsInv, subExons[i].cargo + 1, -1);
            }
            intervalsInv[subExons[i].cargo] = i;
        }
        
        if (group != GROUPANNOT)
        {
            size_t matchEdges = 0;
            TSplitReadIndicesListConstIterator indReadIt = begin(splitReadsIndices[contigId]);
            while (indReadIt != end(splitReadsIndices[contigId]))
            {
                TIndexListConstIterator segmentItL, segmentItR, segmentItEnd;
                segmentItL = segmentItR = begin(indReadIt->indicesLeft);
                segmentItEnd = end(indReadIt->indicesLeft);
                while (++segmentItR != segmentItEnd)
                {
                    size_t subExL = intervalsInv[*segmentItL];
                    size_t subExR = intervalsInv[*segmentItR];
                    
                    if(subExL == subExR)
                    {
                        if (options.verbosity >= 2)
                            std::cerr << "ERROR 1: " << subExL << "  " << length(indReadIt->indicesLeft) << "  " << indReadIt->readId << std::endl;
                        for (unsigned kk = 0; kk < length(indReadIt->indicesLeft); ++kk)
                            std::cerr << indReadIt->indicesLeft[kk] << std::endl;
                        exit(1);
                    }
                    
                    if (subExL >= numVertices(g) || subExR >= numVertices(g))
                    {
                        if (options.verbosity >= 2)
                            std::cerr << "ERROR: " << subExL << "  " << subExR << "  " << numVertices(g) << std::endl;
                    }
                    
                    if (!findEdge(g, subExL, subExR))
                    {
                        addEdge(g, subExL, subExR);
                        ++matchEdges;
                    }
                    
                    ++segmentItL;
                }
                
                if(!empty(indReadIt->indicesRight))
                {
                    segmentItL = segmentItR = begin(indReadIt->indicesRight);
                    segmentItEnd = end(indReadIt->indicesRight);
                    while (++segmentItR != segmentItEnd)
                    {
                        size_t subExL = intervalsInv[*segmentItL];
                        size_t subExR = intervalsInv[*segmentItR];
                        
                        if(subExL == subExR)
                        {
                            std::cerr << "ERROR2: " << subExL << "  " << length(indReadIt->indicesRight) << "  " << indReadIt->readId << std::endl;
                            exit(1);
                        }
                        
                        if (!findEdge(g, subExL, subExR))
                        {
                            addEdge(g, subExL, subExR);
                            ++matchEdges;
                        }
                        
                        ++segmentItL;
                    }
                    const TVertexDescriptor firstSegmentLeft = intervalsInv[front(indReadIt->indicesLeft)];
                    
                    if (firstSegmentLeft != intervalsInv[*segmentItL] && !findEdge(g, firstSegmentLeft, intervalsInv[*segmentItL]))
                        addEdge(g, firstSegmentLeft, intervalsInv[*segmentItL]);
                }
                ++indReadIt;
            }
            if (options.verbosity >= 2)
                std::cerr << "added " << matchEdges << " edges based on read mappings" << std::endl;
        }
        
        //connect subexons according to annotation
        if (group != GROUPMAPPINGS)
        {
            size_t annoEdges = 0;
            if (!byGene.empty())
            {
                for (std::map<size_t, std::vector<size_t> >::const_iterator it = byGene.begin(); it != byGene.end(); ++it)
                {
                    for (size_t sub = 1; sub < it->second.size(); ++sub)
                    {
                        size_t idl = it->second[sub - 1];
                        size_t idr = it->second[sub];
                        
                        //std::cerr << idl << "vs: " << intervalsInv[idl] << "   " << idr << "vs: " << intervalsInv[idr] << std::endl;
                        if (idl >= length(intervalsInv) || idr >= length(intervalsInv) || intervalsInv[idl] == INVALID_ID || intervalsInv[idr] == INVALID_ID)
                            continue;
                        //std::cerr << idl << "vs: " << intervalsInv[idl] << std::endl;
                        if (!findEdge(g, intervalsInv[idl], intervalsInv[idr]))
                        {
                            addEdge(g, intervalsInv[idl], intervalsInv[idr]);
                            //std::cerr << "adding: " << idl << " -- " << idr << std::endl;
                            ++annoEdges;
                        }
                    }
                }
            }
            if (options.verbosity >= 2)
                std::cerr << "added " << annoEdges << " edges based on annotation" << std::endl;
        }
        
        //compute connected components (our Loci afterwards)
        String<TVertexDescriptor> components;
        resizeVertexMap(components, g);
        connectedComponents(components, g);
        
        StringSet< String<unsigned> > connComponents;
        seqan::resize(connComponents, seqan::length(components));
        
        for (unsigned i = 0; i < length(components); ++i)
            appendValue(connComponents[components[i]], i);
        
        //generate MatFile and annotation file
        for (unsigned i = 0; i < length(connComponents); ++i)
        {
            if (empty(connComponents[i]))
                continue;
            
            std::sort(begin(connComponents[i]), end(connComponents[i]));
            
            if (subExons[connComponents[i][0]].i1 > subExons[connComponents[i][0]].i2)
                reverse(connComponents[i]);
            
            //add Locus to annotation File
            TAnnoTreeIter gene = seqan::createRightChild(root);
            TAnnotation &anno_gene = getAnnotation(gene);
            anno_gene.typeId = store.ANNO_GENE;
            std::stringstream strstr;
            strstr << "Locus_" << nextGeneId++;
            setName(gene, CharString(strstr.str()));
            
            //add Subexons to annotation File and MatFile
            //unsigned contigId = subExonsToContig[connComponents[i][0]];
            for (unsigned j = 0; j < length(connComponents[i]); ++j)
            {
                TAnnoTreeIter mrna = seqan::createRightChild(gene);
                TAnnotation &anno_mrna = getAnnotation(mrna);
                anno_mrna.typeId = store.ANNO_MRNA;
                std::stringstream jstrs;
                jstrs << j+1;
                setName(mrna, CharString(strstr.str() + "__" + jstrs.str()));
                
                TAnnoTreeIter exon = seqan::createRightChild(mrna);
                
                TAnnotation &anno_exon = getAnnotation(exon);
                anno_exon.typeId = store.ANNO_EXON;
                
                //Create a subexons having the same positions (required by our quantifier that searches only for subexons)
                TAnnoTreeIter subexon = seqan::createRightChild(exon);
                TAnnotation &annoSub = getAnnotation(subexon);
                annoSub.typeId = SubTypeId;
                annoSub.contigId = contigId;
                annoSub.beginPos = subExons[connComponents[i][j]].i1;
                annoSub.endPos = subExons[connComponents[i][j]].i2;
                std::stringstream nstrs;
                nstrs << subExons[connComponents[i][j]].cargo;
                assignValueByKey(subexon, "NodeId",CharString(nstrs.str()));
                assignValueByKey(subexon, "SpliceEnd","B");
            }
        }
        
        //generate MatFile
        for (unsigned i = 0; i < length(connComponents); ++i)
        {
            if (empty(connComponents[i]))
                continue;
            
            //add Locus
            matFile << '>' << i << std::endl;
            matFile << 'T';
            
            //add Subexons to MatFile
            for (unsigned j = 0; j < length(connComponents[i]); ++j)
                matFile << '\t' << subExons[connComponents[i][j]].cargo;
            
            matFile << std::endl;
        }
        
        //generate StatFile
        for (unsigned i = 0; i < length(subExons); ++i)
        {
            statFile << subExons[i].cargo << '\t' << -1 << '\t' << _abs((int)subExons[i].i2 - (int)subExons[i].i1) << std::endl;
        }
    }
    return 0;
}




template <typename TCountFile, typename TStatFile, typename TMatFile>
int doWork(GffFileIn & gffInFile,
           GffFileOut & /*gffOutFile*/,
           TCountFile & countFile,
           TStatFile & /*matFile*/,
           TMatFile & /*statFile*/,
           const Options & options,
           bool withAnno)
{
    
    typedef FragmentStore<>                                   TFragmentStore;
    typedef typename TFragmentStore::TContigStore             TContigStoreStore;
    typedef typename Value<TContigStoreStore>::Type           TContig;
    typedef typename TContig::TId                             TContigId;
    typedef std::map<CharString, unsigned>                   TAnnoNameStoreInv;
    
    String<TSplitReadIndicesList> splitReadsIndices;
    TIntervalList subExonIntervalsGlobal;
    String<TContigId> subExonsToContig;
    
    size_t nextSubexonsIndex = 0;
    
    String<TIntervalList> subexonIntervalsMap;
    std::map<size_t, std::vector<size_t> > byGene;
    TAnnoNameStoreInv annoNameStoreInv;

    //maps for estimating fragment length and read length distribution
    TLenDistMap mateInnerDistMap, readLenMap;


    
    
    if (withAnno)
    {
        //parse input annotation if available
        readIntervalsFromAnnotation(gffInFile, subexonIntervalsMap, annoNameStoreInv, byGene, options);
        std::cerr<<"Work with given annotation"<<std::endl;
    }
    
    StringSet<CharString> contigNameStore;
    //String<TSplitReadIntervalsList> splitReadsIntervalsMap;
    //String<TIntervalList> intervalsMap;
    //String<TJunctionList> junctionsMap;
    //double rt1 = seqan::sysTime();
    
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(options.samInFile)))
    {
        std::cerr << "Could not open " << options.samInFile << std::endl;
        return 1;
    }
    
    BamHeader header;
    readHeader(header, bamFileIn);
    
    contigNameStore = contigNames(context(bamFileIn));
    resize(splitReadsIndices, length(contigNameStore));
    
    //map bam records to read ids (to group multiple mappings)
    String<size_t> recordToId;
    unsigned totalMappedFrags = groupMultiMapped(contigNameStore, options, recordToId);
    
    size_t startRecordId = 0;
    
    for (size_t contigId = 0; contigId < length(contigNameStore); ++contigId)
    {
        double rt  = seqan::sysTime();
        
        TSplitReadIntervalsList splitReadsIntervals;
        TIntervalList intervals;
        TJunctionList junctions;
        
        CharString contigName = contigNameStore[contigId];
        
        size_t contigIdSubexons = contigId; // if using annotation the ids might be different
        
        if (withAnno)
        {
            TAnnoNameStoreInv::const_iterator aIt = annoNameStoreInv.find(contigName);
            if (aIt != annoNameStoreInv.end())
            {
                contigIdSubexons = aIt->second;
            }
            else
            {
                if (options.verbosity >= 1)
                {
                    std::cerr << "Annotation contains no entries for contig " << contigName << " skipping!" << std::endl;
                }
                continue;
            }
        }
        
        startRecordId = getIntervalsForContig(splitReadsIntervals, intervals, junctions, recordToId, readLenMap, startRecordId, contigId, options);
        
        if (empty(intervals))
        {
            if (options.verbosity >= 1)
                std::cerr << "no intervals parsed! Aborting." << std::endl;
            
            continue;
        }
        
        if (!withAnno)
        {
            if (options.verbosity >= 1)
                std::cerr<<"Working without given annotation"<<std::endl;
            
            processIntervals(intervals, junctions, options, subexonIntervalsMap[contigIdSubexons], nextSubexonsIndex);
        }
        
        const TIntervalList& subexonIntervals = subexonIntervalsMap[contigIdSubexons];
        
        if (options.verbosity >= 1)
            std::cerr << "Subexons for chromosome: " << length(subexonIntervals) << std::endl;
        
        double rtTrans = seqan::sysTime();
        intervalIndexTransform(subexonIntervals, splitReadsIntervals, splitReadsIndices[contigId], mateInnerDistMap, options);
        if (options.verbosity >= 1)
            std::cerr << "Transform took " << seqan::sysTime() - rtTrans << " seconds" << std::endl;
        
        if (options.verbosity >= 1)
            std::cerr << "Done interval transform. Transformed " << length(splitReadsIndices[contigId]) << " reads" << std::endl;
        
        if (options.verbosity >= 1)
            std::cerr << "Chromosome " << contigName << " took " << seqan::sysTime() - rt << " sec." << std::endl;
    }
    
    //for every read count how often it appears in split_reads_interval (due to mutliple mapping or overlapping subexons on reverse strands)
    String<unsigned> mappingCounts;
    
    double rtAdjCounts = seqan::sysTime();
    
    unsigned maxId = 0;
    //get maximum read id
    for (size_t i = 0; i < length(splitReadsIndices); ++i)
    {
        for (size_t j = 0; j < length(splitReadsIndices[i]); ++j)
        {
            if (maxId < splitReadsIndices[i][j].readId)
                maxId = splitReadsIndices[i][j].readId;
        }
    }
    
    resize(mappingCounts, maxId + 1, 0);
    
    for (size_t i = 0; i < length(splitReadsIndices); ++i)
    {
        for (size_t j = 0; j < length(splitReadsIndices[i]); ++j)
        {
            ++mappingCounts[splitReadsIndices[i][j].readId];
        }
    }
    if (options.verbosity >= 1)
        std::cerr << "Adjusting mapping counts took " << seqan::sysTime() - rtAdjCounts << " sec." << std::endl;
    
    double rtgenCounts = seqan::sysTime();
    TTupleCountList counts;
#ifdef GENERATEGAPS
    TTupleCountGapList countsGaps;
    generateCounts(splitReadsIndices, mappingCounts, counts, countsGaps, options);
#else
    generateCounts(splitReadsIndices, mappingCounts, counts, options);
#endif
    
    if (options.verbosity >= 1)
        std::cerr << "Count generation took " << seqan::sysTime() - rtgenCounts << " sec." << std::endl;

    //estimate fragment length distribution and avg. read length

    //first read length
    double meanRl = 0;
    unsigned numR = 0;
    for (TLenDistMap::const_iterator it = readLenMap.begin(); it != readLenMap.end(); ++it)
    {
      meanRl += it->first * it->second;
      numR += it->second;
    }
    meanRl /= numR;

    //Now fragment length
    double meanInner= 0;
    double sdevInner = 0;
    unsigned numInner = 0;

    for (TLenDistMap::const_iterator it = mateInnerDistMap.begin(); it != mateInnerDistMap.end(); ++it)
    {
      meanInner += it->first * it->second;
      numInner += it->second;
    }
    meanInner /= numInner;

    for (TLenDistMap::const_iterator it = mateInnerDistMap.begin(); it != mateInnerDistMap.end(); ++it)
    {
      double dist = fabs(it->first - meanInner);
      sdevInner += (it->second * dist * dist) / numInner;
    }
    sdevInner =  sqrt(sdevInner);

    LibInfo libInfo;
    libInfo.rLenMean = meanRl;
    libInfo.numR = numR;
    libInfo.fLenMean = meanInner + 2 * meanRl;
    libInfo.fLenSdev = sdevInner;
    libInfo.numFLen = numInner;
    libInfo.totalFrags = totalMappedFrags;

    std::cout << "Estimated paired-end reads inner distance distribution (i.e. fragment length - 2* read length) from " << numInner << " paired-end alignments:" << std::endl;
    std::cout << "Mean: " << meanInner << "\t" << "S_dev: " << sdevInner << std::endl;

    
#ifdef GENERATEGAPS
    writeCountFile(countFile, counts, countsGaps, options, libInfo);
#else
    writeCountFile(countFile, counts, options, libInfo);
#endif    

    //if (options.verbosity >= 1)
    //    std::cerr << "Done count outfile for contig" << std::endl;
    //writeMatStat(matFile, statFile, annoStore, subexonIntervalsMap, splitReadsIndices, options, byGene);
    //writeRecords(gffOutFile, annoStore);
    
    return 0;
}


int main(int argc, char const ** argv)
{
    using namespace seqan;
    
    // -----------------------------------------------------------------------
    // Handle Command Line
    // -----------------------------------------------------------------------
    
    // Setup command line parser.
    ArgumentParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int res = parseCommandLineAndCheck(options, parser, argc, argv);
    if (res != 0)
        return 1;
    if (options.showHelp || options.showVersion)
        return 0;
    
    // -----------------------------------------------------------------------
    // Do Work.
    // -----------------------------------------------------------------------
    seqan::GffFileOut gtfOutFile;
    /*
     if (!open(gtfOutFile, toCString(options.gtfOutFile)))
     {
     std::cerr << "Could not open " << options.gtfOutFile << std::endl;
     return 1;
     }
     
     std::ofstream matFile(toCString(options.matFile), std::ios::binary);
     if (!matFile.is_open())
     {
     std::cerr << "Could not open " << options.matFile << std::endl;
     return 1;
     }
     
     std::ofstream statFile(toCString(options.statFile), std::ios::binary);
     if (!statFile.is_open())
     {
     std::cerr << "Could not open " << options.statFile << std::endl;
     return 1;
     }
     */
    std::ofstream matFile, statFile; // dummy, not used
    
    std::ofstream countFile(toCString(options.countFile));
    if (!countFile.is_open())
    {
        std::cerr << "Could not open " << options.countFile << std::endl;
        return 1;
    }
    
    // does bam index exist?
    std::ifstream bamIndexTest(toCString(options.bamIndexFile));
    if (!bamIndexTest)
    {
        std::cerr << "Could not open required bam index file" << options.bamIndexFile << std::endl;
        std::cerr << "Please create index using samtools" << std::endl;
        return 1;
    }
    else
    {
        bamIndexTest.close();
    }
    
    //open the optional input annotation (refined.gff)
    if (!empty(options.gtfInFile))
    {
        GffFileIn gtfInfile;
        if (!open(gtfInfile, toCString(options.gtfInFile)))
        {
            std::cerr << "Could not open " << options.gtfInFile << std::endl;
            return 1;
        }
        
        doWork(gtfInfile, gtfOutFile, countFile, matFile, statFile, options, true);
    }
    else
    {
        GffFileIn placeholder;
        doWork(placeholder, gtfOutFile, countFile, matFile, statFile, options, false);
    }
    
    
    countFile.close();
    //matFile.close();
    //statFile.close();
    
    return 0;
}

