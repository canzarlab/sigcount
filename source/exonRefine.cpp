// ==========================================================================
//                                  exonRefine
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

//#define SEQAN_ENABLE_DEBUG 1
#include <fstream>
#include <iostream>
#include <string>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/store.h>

using namespace seqan;

// if SEPARATE_STRANDS is defined, forward and reverse strands are separated and will have no common subexons
#define SEPARATE_STRANDS
#define COMPACT_NODE_IDS

// if LOCAL_SPLICEEND is defined, 'SpliceEnd' will refer to reading direction (5'->3', not forward strand position)
// makes only sense if SEPARATE_STRANDS is defined
#define LOCAL_SPLICEEND

#ifdef COMPACT_NODE_IDS
#define LARGE_NODE_ID "LargeNodeId"
#else
#define LARGE_NODE_ID "NodeId"
#endif

template <typename TValue>
inline TValue absInt(TValue t)
{
    return (t >= 0)? t: -t;
}

template <typename TValue>
inline TValue maskInt(TValue t)
{
    return t & ((1ull << (8 * sizeof(TValue) - 4)) - 1);
}

template <typename TSequence>
void removeEqualElements(TSequence &seq)
{
    typedef typename Iterator<TSequence, Standard>::Type TIter;
    
    TIter itBegin = begin(seq, Standard());
    TIter itEnd = end(seq, Standard());
    TIter src = itBegin;
    TIter dst = itBegin;
    for (; src != itEnd; ++src)
    {
        if (*dst != *src)
        {
            ++dst;
            *dst = *src;
        }
    }
    if (itBegin != itEnd)
        resize(seq, dst - itBegin + 1);
}

template <typename TSequence>
void removeEqualBoundaries(TSequence &seq)
{
    typedef typename Iterator<TSequence, Standard>::Type TIter;
    typedef typename Value<TSequence>::Type TContigExonBound;
    TContigExonBound intronBelow = 1ull << (8 * sizeof(TContigExonBound) - 2);
    TContigExonBound intronAbove = 1ull << (8 * sizeof(TContigExonBound) - 1);
    
    TIter itBegin = begin(seq, Standard());
    TIter itEnd = end(seq, Standard());
    TIter src = itBegin;
    TIter dst = itBegin;
    for (; src != itEnd; ++src)
    {
        if (maskInt(*dst) != maskInt(*src))
        {
            ++dst;
            *dst = *src;
        } else
        {
            *dst |= *src & (intronBelow | intronAbove | (intronBelow >> 2) | (intronAbove >> 2));
        }
    }
    if (itBegin != itEnd)
        resize(seq, dst - itBegin + 1);
}


template <typename TContigExonBound>
struct BoundariesLess
{
    inline bool operator() (TContigExonBound a, TContigExonBound b)
    {
        return maskInt(a) < maskInt(b);
    }
};

// Step 3a: Get Exon Boundaries
//
// Go through all annotated exons and collect their boundaries
// group according to strand and orientation.
// Afterwards sort boundaries and remove duplicated elements

template <typename TExonBoundString, typename TFragmentStore>
void getExonBoundaries(TExonBoundString &exonBounds, TFragmentStore &store)
{
    typedef typename Value<TExonBoundString>::Type TContigExonBounds;
    typedef typename  Value<TContigExonBounds>::Type TContigExonBound;
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
    typedef typename Value<TAnnotationStore>::Type TAnnotation;
    
    clear(exonBounds);
    TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
    TContigExonBound intronBelow = 1ull << (8 * sizeof(TContigExonBound) - 2);
    TContigExonBound intronAbove = 1ull << (8 * sizeof(TContigExonBound) - 1);
    
    while (!atEnd(dfsIt))
    {
        TAnnotation &anno = getAnnotation(dfsIt);
        if (anno.typeId == TFragmentStore::ANNO_EXON)
        {
            // the reverse strands are interleaved with the forward strands
            if (anno.contigId == TAnnotation::INVALID_ID)
            {
                std::cerr << "Exon annotation has no contig (parentId=" << getParentName(dfsIt) << ")." << std::endl;
                continue;
            }
            unsigned contigNum = anno.contigId;
#ifdef SEPARATE_STRANDS
            contigNum *= 2;
            if (anno.beginPos > anno.endPos) ++contigNum;
#endif
            
            if (length(exonBounds) <= contigNum)
                resize(exonBounds, contigNum + 1);
            if (anno.beginPos < anno.endPos)
            {
                appendValue(exonBounds[contigNum], anno.beginPos | intronBelow);
                appendValue(exonBounds[contigNum], anno.endPos   | intronAbove);
            } else {
                appendValue(exonBounds[contigNum], anno.endPos   | (intronBelow >> 2));
                appendValue(exonBounds[contigNum], anno.beginPos | (intronAbove >> 2));
            }
            goNextRight(dfsIt);
        }
        else
            goNext(dfsIt);
    }
    for (unsigned i = 0; i < length(exonBounds); ++i)
    {
        std::sort(begin(exonBounds[i], Standard()), end(exonBounds[i], Standard()), BoundariesLess<TContigExonBound>());
        removeEqualBoundaries(exonBounds[i]);
    }
}


// Step 3b: Refine Exon Boundaries
//
// For every exon interval get all contained boundaries and
// divide exon into subexons using these boundaries.
// Insert subexons as children of the exon (even if there is only one subexon).
// NodeIds are set according to the boundary rank in the
// sorted list of boundaries (the total rank is determined by a partial sum).

template <typename TExonBoundString, typename TFragmentStore, typename TNodeIds>
void refineExonBoundaries(TExonBoundString &exonBounds, TFragmentStore &store, TNodeIds &nodeIds)
{
    typedef typename Value<TExonBoundString>::Type TContigExonBounds;
    typedef typename  Value<TContigExonBounds>::Type TContigExonBound;
    typedef typename Iterator<TContigExonBounds, Standard>::Type TContigExonBoundIterator;
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename  Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
    typedef typename  Value<TAnnotationStore>::Type TAnnotation;
    //typedef typename  Id<TAnnotation>::Type TId;
    TContigExonBound intronBelow  = 1ull << (8 * sizeof(TContigExonBound) - 2);
    TContigExonBound intronAbove  = 1ull << (8 * sizeof(TContigExonBound) - 1);
    TContigExonBound intronBelowR = 1ull << (8 * sizeof(TContigExonBound) - 4);
    TContigExonBound intronAboveR = 1ull << (8 * sizeof(TContigExonBound) - 3);
    
    // 5'------------------------------------------>3'
    // below                                    above
    //
    //   ,--------.      ,----------.     ,---------.
    //   *        *      *          *     *         *
    //   '--------'      '----------'     '---------'
    //        ,------.      ,---.        ----.  ,-----
    //        *      *      *   *            *  *
    //        '------'      '---'        ----'  '-----
    // ===============================================
    //   ,---.,---.,-.   ,-.,---.,--.     ,--.,-.,--.
    //   * L |* B *|R*   *L|* B *| R*     *B *|-|* B*
    //   '---''---''-'   '-''---''--'     '--''-''--'
    //
    // Unless LOCAL_SPLICEEND is defined, 'L' and 'R' refer to the fwd strand.
    //
    // 'SpliceEnd' refers to the current subexon (no matter if on fwd or rev strand)
    // 'SpliceEndOther' refers to the opposite strand even if there is no annotation
    
    
    // splice tags:
    // '-' = subexon is in the middle of each exon and cannot be spliced with other subexons
    // 'B' = subexon is the exon itself or a real overlap of two subexons and can be spliced at BOTH ends
    // 'L' = subexon can only be spliced at the 5' prime end
    // 'R' = subexon can only be spliced at the 3' prime end
    char spliceChar[4] = { '-', 'L', 'R', 'B' };
    char spliceCharMirror[4] = { '-', 'R', 'L', 'B' };
    
    String<int> partialSum;
    for (unsigned i = 0, sum = 0; i < length(exonBounds); ++i)
    {
        appendValue(partialSum, sum);
        sum += length(exonBounds[i]);
    }
    
    TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
    while (!atEnd(dfsIt))
    {
        TAnnotation &anno = getAnnotation(dfsIt);
        if (anno.typeId == TFragmentStore::ANNO_EXON)
        {
            // the reverse strands are interleaved with the forward strands
            if (anno.contigId == TAnnotation::INVALID_ID)
            {
                std::cerr << "Exon annotation has no contig (parentId=" << getParentName(dfsIt) << ")." << std::endl;
                continue;
            }
            
            unsigned contigNum = anno.contigId;
#ifdef SEPARATE_STRANDS
            contigNum *= 2;
            if (anno.beginPos > anno.endPos) ++contigNum;
#endif
            
            TContigExonBoundIterator boundsBegin = begin(exonBounds[contigNum], Standard());
            TContigExonBoundIterator boundsEnd = end(exonBounds[contigNum], Standard());
            
            // determine exon orientation and order boundaries
            bool forward = true;
            TContigExonBound beginPos = anno.beginPos;
            TContigExonBound endPos = anno.endPos;
            if (beginPos > endPos)
            {
                forward = false;
                beginPos = endPos;
                endPos = anno.beginPos;
            }
            
            TContigExonBoundIterator first = std::lower_bound(boundsBegin, boundsEnd, beginPos, BoundariesLess<TContigExonBound>());
            TContigExonBoundIterator last = std::lower_bound(boundsBegin, boundsEnd, endPos, BoundariesLess<TContigExonBound>());
            
            TContigExonBound subExonBeginPos = beginPos;
            if (first != last)
            {
                TContigExonBound prevMask = *first;
                ++first;
                for (TContigExonBoundIterator it = first; it <= last; ++it)
                {
                    std::stringstream tmp;
                    TContigExonBound subExonEndPos = maskInt(*it);
                    TAnnoTreeIter childIt = createRightChild(dfsIt);
                    clearValues(childIt);
                    int nodeId;
                    TContigExonBound intronMask = (*it & (intronAbove | intronAboveR)) | (prevMask & (intronBelow | intronBelowR));
                    if (forward)
                    {       // positive strand
                        getAnnotation(childIt).beginPos = subExonBeginPos;
                        getAnnotation(childIt).endPos = subExonEndPos;
                        nodeId = (it - boundsBegin) + partialSum[contigNum];
                    } else
                    {       // negative strand
                        getAnnotation(childIt).endPos = subExonBeginPos;
                        getAnnotation(childIt).beginPos = subExonEndPos;
#ifdef SEPARATE_STRANDS
                        // for increasing nodeIds
                        nodeId = (boundsEnd - it - 1) + partialSum[contigNum];
#else
                        // nodeIds aren't increasing if we merge the strands
                        nodeId = (it - boundsBegin) + partialSum[contigNum];
#endif
                    }
                    SEQAN_ASSERT_LT(beginPos, endPos);
                    SEQAN_ASSERT_LEQ(beginPos, subExonBeginPos);
                    SEQAN_ASSERT_LEQ(subExonEndPos, endPos);
                    tmp << nodeId;
#ifdef COMPACT_NODE_IDS
                    append(nodeIds, nodeId);
#endif
                    setType(childIt, "subexon");
                    assignValueByKey(childIt, LARGE_NODE_ID, tmp.str());
                    unsigned spliceEnd = intronMask >> (8 * sizeof(TContigExonBound) - 2);
                    unsigned spliceEndR = (intronMask >> (8 * sizeof(TContigExonBound) - 4)) & 3;
                    
                    // we want spliceEnd refer to the exon's strand and spliceEndR to the opposite strand
                    if (!forward)
                        std::swap(spliceEnd, spliceEndR);
                    
                    
#ifdef LOCAL_SPLICEEND
                    if (forward)
                    {
                        assignValueByKey(childIt, "SpliceEnd", spliceChar[spliceEnd]);
                        assignValueByKey(childIt, "SpliceEndOther", spliceCharMirror[spliceEndR]);
                    }
                    else
                    {
                        assignValueByKey(childIt, "SpliceEnd", spliceCharMirror[spliceEnd]);
                        assignValueByKey(childIt, "SpliceEndOther", spliceChar[spliceEndR]);
                    }
#else
                    assignValueByKey(childIt, "SpliceEnd", spliceChar[spliceEnd]);
                    assignValueByKey(childIt, "SpliceEndOther", spliceChar[spliceEndR]);
#endif
                    
                    subExonBeginPos = subExonEndPos;
                    prevMask = *it;
                }
            }
            goNextRight(dfsIt);
        }
        else
            goNext(dfsIt);
    }
}

// Step 4: Calculate Orderings
//
// For every exon interval get all contained boundaries and
// divide exon into subexons using these boundaries.
// Fill ordering array with subexon ids for every transcript.
// NodeIds are set according to the boundary rank in the sorted list of boundaries.

template <typename TContigOrderings, typename TFragmentStore>
void calculateOrderings(TContigOrderings &orderings, TFragmentStore &store)
{
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
    typedef typename Value<TAnnotationStore>::Type TAnnotation;
    
    unsigned subExonType = 0;
    _storeAppendType(store, subExonType, "subexon");
    String<char, CStyle> str;
    Dna5String tmp;
    
    TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
    while (!atEnd(dfsIt))
    {
        TAnnotation &anno = getAnnotation(dfsIt);
        if (anno.typeId == subExonType)
        {
            unsigned mRNAId = value(nodeUp(nodeUp(dfsIt)));
            
            if (length(orderings) <= mRNAId)
                resize(orderings, mRNAId + 1);
            
            TAnnotation &anno = getAnnotation(dfsIt);
            
            unsigned i = 0;
            if (anno.beginPos < anno.endPos)
            {
                for (; i < length(orderings[mRNAId]); ++i)
                {
                    TAnnotation &anno2 = store.annotationStore[orderings[mRNAId][i]];
                    if (anno.beginPos < anno2.beginPos) break;
                }
            }
            else
            {
                for (; i < length(orderings[mRNAId]); ++i)
                {
                    TAnnotation &anno2 = store.annotationStore[orderings[mRNAId][i]];
                    if (anno.beginPos > anno2.beginPos) break;
                }
            }
            
            insertValue(orderings[mRNAId], i, value(dfsIt));
            goNextRight(dfsIt);
        }
        else
            goNext(dfsIt);
    }
}

template <typename TStream, typename TContigOrderings, typename TFragmentStore>
void writeTranscripts(TStream &target, TContigOrderings &orderings, TFragmentStore &store)
{
    CharString id;
    CharString tmp, transcript;
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename TFragmentStore::TContigPos TContigPos;
    typedef typename Value<TAnnotationStore>::Type TAnnotation;
    
    for (unsigned i = 0; i < length(orderings); ++i)
        if (!empty(orderings[i]))
        {
            id = "Locus_";
            append(id, store.annotationNameStore[store.annotationStore[i].parentId]);
            append(id, "_Transcript_");
            append(id, store.annotationNameStore[i]);
            append(id, "_Confidence_99.99");
            
            clear(transcript);
            for (unsigned j = 0; j < length(orderings[i]); ++j)
            {
                TAnnotation &anno = store.annotationStore[orderings[i][j]];
                if (_max(anno.beginPos, anno.endPos) < (TContigPos)length(store.contigStore[anno.contigId].seq))
                {
                    if (anno.beginPos < anno.endPos)
                        append(transcript, infix(store.contigStore[anno.contigId].seq, anno.beginPos, anno.endPos));
                    else
                    {
                        tmp = infix(store.contigStore[anno.contigId].seq, anno.endPos, anno.beginPos);
                        reverseComplement(tmp);
                        toLower(tmp);
                        append(transcript, tmp);
                    }
                }
                else
                {
                    std::cerr << "Genomic sequence missing for " << store.contigNameStore[anno.contigId];
                    std::cerr << " from " << anno.beginPos << " to " << anno.endPos << std::endl;
                }
                
            }
            writeRecord(target, id, transcript);
        }
}

template <typename TStream, typename TOrderings, typename TFragmentStore, typename TNodeIds>
void writeOrderings(TStream &/*target*/, TOrderings &orderings, TFragmentStore &store, TNodeIds &nodeIds)
{
    typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type TAnnotation;
    
#ifdef COMPACT_NODE_IDS
    std::sort(begin(nodeIds, Standard()), end(nodeIds, Standard()));
    removeEqualElements(nodeIds);
#endif
    
    CharString nodeId;
    String<char, CStyle> idStr;
    for (unsigned i = 0; i < length(orderings); ++i)
        if (!empty(orderings[i]))
        {
            //target << ">Locus_" << store.annotationNameStore[store.annotationStore[i].parentId];
            //target << "_Transcript_" << store.annotationNameStore[i];
            //target << "_Confidence_99.99\n";
            
            int endPos = 0;
            int lastId = -1;
            __int64 lastPos = -1;
            bool silent = false;
            for (unsigned j = 0; j < length(orderings[i]); ++j)
            {
                TAnnotation &anno = store.annotationStore[orderings[i][j]];
                //if (j > 0)
                    //target << "-(0)->";
                //if (anno.beginPos > anno.endPos)
                    //target << '-';
                if (!annotationGetValueByKey(idStr, store, anno, LARGE_NODE_ID))
                    std::cerr << "nodeId is missing for exon " << store.annotationNameStore[orderings[i][j]] << std::endl;
                int nodeId = atoi(idStr);
#ifdef COMPACT_NODE_IDS
                typename Iterator<TNodeIds, Standard>::Type it = std::lower_bound(begin(nodeIds, Standard()), end(nodeIds, Standard()), nodeId);
                nodeId = it - begin(nodeIds, Standard()) + 1;
                std::stringstream tmp;
                tmp << nodeId;
                annotationAssignValueByKey(store, anno, "NodeId", tmp.str());
#endif
                
#ifdef SEPARATE_STRANDS
                if (absInt(nodeId) <= lastId && !silent)
                {
                    std::cerr << "nodeId is not increasing in "
                    << "Locus_" << store.annotationNameStore[store.annotationStore[i].parentId]
                    << "_Transcript_" << store.annotationNameStore[i]
                    << "_Confidence_99.99" << std::endl;
                    for (unsigned k = 0; k < length(orderings[i]); ++k)
                        std::cerr << '\t' << orderings[i][k];
                    std::cerr << std::endl;
                    silent = true;
                }
#endif
                
                if (((anno.beginPos < anno.endPos && (__int64)anno.beginPos <= lastPos) ||
                     (anno.beginPos > anno.endPos && (__int64)anno.beginPos >= lastPos)) && !silent && lastId != -1)
                {
                    std::cerr << "beginPos is not monotonic in "
                    << "Locus_" << store.annotationNameStore[store.annotationStore[i].parentId]
                    << "_Transcript_" << store.annotationNameStore[i]
                    << "_Confidence_99.99" << std::endl;
                    for (unsigned k = 0; k < length(orderings[i]); ++k)
                        std::cerr << '\t' << store.annotationStore[orderings[i][k]].beginPos << '-' << store.annotationStore[orderings[i][k]].endPos;
                    std::cerr << std::endl;
                    silent = true;
                }
                //target << nodeId;
                //target << ':';
                endPos += absInt((__int64)anno.endPos - (__int64)anno.beginPos);
                //target << endPos;
                lastId = absInt(nodeId);
                lastPos = anno.beginPos;
            }
            //target << '\n';
        }
}

//temporary fix of the read records which generates a cyclic annotation
//tree if additional gene entries are given in gtf
//will be replaced when fixed in seqan
template <typename TSpec, typename TConfig>
inline void
readRecords_tmp_fix(FragmentStore<TSpec, TConfig> & fragStore,
            GffFileIn & gffFile)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (atEnd(gffFile))
        return;

    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.annotationNameStoreCache);
    refresh(fragStore.annotationTypeStoreCache);

    GffRecord record;
    IOContextGff_<TFragmentStore> ctx;

    while (!atEnd(gffFile))
    {
        readRecord(record, gffFile);
        if(record.type != "exon")
            continue;
        _readOneAnnotation(ctx, record, format(gffFile));
        _storeOneAnnotation(fragStore, ctx);
    }
    _storeClearAnnoBackLinks(fragStore.annotationStore);
    _storeCreateAnnoBackLinks(fragStore.annotationStore);
    _storeRemoveTempAnnoNames(fragStore);
}



int main(int argc, char const * argv[])
{
    typedef FragmentStore<> TFragmentStore;
    typedef String<__uint64> TContigExonBounds;
    typedef String<TContigExonBounds> TExonBounds;
    //typedef Id<TFragmentStore>::Type TId;
    
    ArgumentParser parser;
    setAppName(parser, "exonRefine");
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
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "ANNOTATIONS", true));
    std::vector<std::string> exts = GffFileIn::getFileExtensions();
    append(exts, UcscFileIn::getFileExtensions());
    setValidValues(parser, 0, exts);
    
    addOption(parser, ArgParseOption("p", "prefix", "Name for output file.", ArgParseArgument::OUTPUT_PREFIX));
    setDefaultValue(parser, "prefix", "./refined");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP]... <\\fIannotation.gff\\fP>");
    addUsageLine(parser, "[\\fIOPTIONS\\fP]... <\\fIannotation.gtf\\fP>");
    
    addDescription(parser, "exonRefine refines the set of exons into non-overlapping segments");
    
    addDescription(parser, "(c) Copyright 2010-2014 by David Weese.");
    
    //requiredArguments(parser, 2);
    
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK) return 0;
    
    //double t0 = sysTime();
    TFragmentStore store;
    
    CharString prefix;
    getOptionValue(prefix, parser, "prefix");
    //CharString orderingFileName(prefix);
    //CharString transcriptsFileName(prefix);
    CharString refinedFileName(prefix);
    //append(orderingFileName, "/ordering.fa");
    //append(transcriptsFileName, "/transcripts.fa");
    //append(refinedFileName, "./refined.gtf");
    append(refinedFileName, ".gtf");

    
    //////////////////////////////////////////////////////////////////////////////
    // Step 1: Reading the genome
    double t1 = sysTime();
    /*
     if (isSet(parser, "genome"))
     {
     CharString tmp;
     getOptionValue(tmp, parser, "genome");
     bool success = loadContigs(store, tmp, true);
     if (!success)
     {
     std::cerr << "Could not read genome " << tmp << std::endl;
     return 1;
     }
     t1 = sysTime();                std::cout << "Reading the genome took "<< t1-t0 << " seconds." << std::endl;
     }
     */
    
    //////////////////////////////////////////////////////////////////////////////
    // Step 2: Reading the annotation
    if (getArgumentValueCount(parser, 0) == 1) // Gtf/Gff?
    {
        GffFileIn annotationFile;
        if (!open(annotationFile, getArgumentValues(parser, 0)[0].c_str()))
        {
            std::cerr << "Could not open annotation " << getArgumentValues(parser, 0)[0] << std::endl;
            return 1;
        }
        readRecords_tmp_fix(store, annotationFile);
    }
    else
    {
        UcscFileIn knownGenes;
        UcscFileIn knownIsoforms;
        if (!open(knownGenes, getArgumentValues(parser, 0)[0].c_str()))
        {
            std::cerr << "Could not open knownGene " << getArgumentValues(parser, 0)[0] << std::endl;
            return 1;
        }
        if (!open(knownIsoforms, getArgumentValues(parser, 0)[1].c_str()))
        {
            std::cerr << "Could not open knownIsoforms " << getArgumentValues(parser, 0)[1] << std::endl;
            return 1;
        }
        readRecords(store, knownGenes);
        readRecords(store, knownIsoforms);
    }
    double t2 = sysTime();          std::cout << "Reading the annotation took "<< t2-t1 << " seconds." << std::endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // Step 3: Get exon boundaries and refine into subexons
    TExonBounds exonBounds;
    String<int> nodeIds;
    getExonBoundaries(exonBounds, store);
    std::cout << "Collected all exon boundaries." << std::endl;
    refineExonBoundaries(exonBounds, store, nodeIds);
    double t3 = sysTime();          std::cout << "Refining exons took "<< t3-t2 << " seconds." << std::endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // Step 4: Calculate orderings
    StringSet<String<int> > orderings;
    calculateOrderings(orderings, store);
    //double t4 = sysTime();          std::cout << "Calculating orderings took "<< t4-t3 << " seconds." << std::endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // Step 5: Write ordering file
    
    //std::ofstream orderingsFile(toCString(orderingFileName), ::std::ios_base::out | ::std::ios_base::binary);
    std::ofstream orderingsFile;
    writeOrderings(orderingsFile, orderings, store, nodeIds);
    //orderingsFile.close();
    
    double t5 = sysTime();          //std::cout << "Writing contig orderings took "<< t5-t4 << " seconds." << std::endl;
    //////////////////////////////////////////////////////////////////////////////
    // Step 6: Write transcripts
    double t6 = t5;
    /*
     if (isSet(parser, "genome"))
     {
     SeqFileOut transcriptsFile(toCString(transcriptsFileName));
     writeTranscripts(transcriptsFile, orderings, store);
     t6 = sysTime();                std::cout << "Writing spliced transcipts took "<< t6-t5 << " seconds." << std::endl;
     }
     */
    //////////////////////////////////////////////////////////////////////////////
    // Step 7: Write the annotation enhanced by subexons
    //GffFileOut file2(toCString(refinedFileName), ::std::ios_base::out | ::std::ios_base::binary);
    GffFileOut file2(toCString(refinedFileName));
    
    writeRecords(file2, store);
    double t7 = sysTime();          std::cout << "Writing the enhanced annotation took "<< t7-t6 << " seconds." << std::endl;
    
    return 0;
}
