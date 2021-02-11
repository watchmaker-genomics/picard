/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.arrays.illumina;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang.StringUtils;
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A class to represent a record (line) from an Extended Illumina Manifest [Assay] entry
 */
public class Build37ExtendedIlluminaManifestRecord extends IlluminaManifestRecord {
    protected enum Flag {
        /** Flagged by Illumina as a bad assay */
        ILLUMINA_FLAGGED,
        LIFTOVER_FAILED,
        UNSUPPORTED_GENOME_BUILD,

        /** Probe sequence not found in reference. */
        PROBE_SEQUENCE_MISMATCH,

        /** Probe sequence is on unexpected strand. */
        PROBE_SEQUENCE_STRAND_INVALID,

        /** Source sequenc not found in reference. */
        SOURCE_SEQUENCE_MISMATCH,

        /** Source sequence is invalid (contains invalid character). */
        SOURCE_SEQUENCE_INVALID,

        /** Source sequence is on unexpected strand. */
        SOURCE_SEQUENCE_STRAND_INVALID,

        /** Neither insertion nor deletion sequence found in reference. */
        INDEL_NOT_FOUND,

        /** Both insertion and deletion sequence found in reference. */
        INDEL_CONFLICT,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        SEQUENCE_MISMATCH,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        INDEL_SEQ_MISMATCH,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        INDEL_EXTENSION_ERROR,
        DUPE,
        PASS,
    }

    private final IlluminaBPMLocusEntry locusEntry;
    // TODO - switch over to using illuminaManifestRecord as a private member - do NOT extend it.
    private final IlluminaManifestRecord illuminaManifestRecord;

    private String b37Chr;
    private Integer b37Pos;
    private String snpRefAllele;
    private String snpAlleleA;
    private String snpAlleleB;
    private String rsId;
    private Flag flag = Flag.PASS;

    private Allele A;
    private Allele B;
    private Allele ref;

    // The refStrand if provided in the Illumina manifest, otherwise calculated
    // TODO - probably rename this back to refStrand once this class no longer extends IlluminaManifestRecord
    private Strand referenceStrand = null;

    private final Log log = Log.getInstance(Build37ExtendedIlluminaManifestRecord.class);

    // These are the IUPAC nucleotide codes as described here: https://www.bioinformatics.org/sms/iupac.html
    private static final String IUPAC_NUCLEOTIDE_CODES = "ACGTRYSWKMBDHVNacgtryswkmbdhvn";
    private static final String ACGT_CODES = "ACGTacgt";

    private static final String SRC_SEQ_REGEX = "([" + IUPAC_NUCLEOTIDE_CODES + "]*)\\[([" + ACGT_CODES + "-])\\/([" + ACGT_CODES + "]*)\\]([" + IUPAC_NUCLEOTIDE_CODES + "]*)";
    private static final Pattern pattern = Pattern.compile(SRC_SEQ_REGEX);

    private static final String ACGT_REGEX = "^[" + ACGT_CODES + "]+$";
    private static final Pattern ACGT_PATTERN = Pattern.compile(ACGT_REGEX);

    // Symbolics for the regex groups...
    private static final int FIVE_PRIME_SEQUENCE = 1;
    private static final int PRE_INDEL_SEQUENCE = 2;
    private static final int INDEL_SEQUENCE = 3;
    private static final int THREE_PRIME_SEQUENCE = 4;

    static final String BUILD_36 = "36";
    static final String BUILD_37 = "37";

    /**
     * This constructor is used to read records from an already created Build37ExtendedIlluminaManifestRecord file.
     * It does not work to set the Extended-specific fields
     */
    Build37ExtendedIlluminaManifestRecord(final Map<String, Integer> columnNameToIndex, final String[] line, final int index) {
        // TODO - I think this could be simplified to not necessarily rely on it including any of the Illumina Manifest records
        // but that should be done in the calling method I think
        // Not absolutely sure why this constructor needs to be here?
        super(columnNameToIndex, line, index);
        this.locusEntry = null;     // TODO - don't like this.
        this.illuminaManifestRecord = null;

        final int end = line.length;
        flag = Flag.valueOf(line[end - 1]);

        if (!isBad()) {
            b37Chr = line[end - 7];
            b37Pos = parseIntOrNull(line[end - 6]);
            snpRefAllele = line[end - 5];
            snpAlleleA = line[end - 4];
            snpAlleleB = line[end - 3];
            rsId = line[end - 2];

            A = Allele.create(snpAlleleA, snpAlleleA.equals(snpRefAllele));
            B = Allele.create(snpAlleleB, snpAlleleB.equals(snpRefAllele));
            ref = Allele.create(snpRefAllele, true);
        } else {
            b37Chr = "0";
            b37Pos = 0;
            snpRefAllele = "";
            snpAlleleA = "";
            snpAlleleB = "";
            rsId = "";

            A = Allele.NO_CALL;
            B = Allele.NO_CALL;
            ref = Allele.NO_CALL;
        }
    }

    // TODO - where is this constructor used?  Needed?
    // Probably testing.
    Build37ExtendedIlluminaManifestRecord(final IlluminaManifestRecord record,
                                   final Flag flag,
                                   final String b37Chr,
                                   final int b37Pos,
                                   final String snpRefAllele,
                                   final String snpAlleleA,
                                   final String snpAlleleB,
                                   final String rsId) {
        super(record);
        this.locusEntry = null;
        this.illuminaManifestRecord = null;
        this.flag = flag;
        this.b37Chr = b37Chr;
        this.b37Pos = b37Pos;
        this.snpRefAllele = snpRefAllele;
        this.snpAlleleA = snpAlleleA;
        this.snpAlleleB = snpAlleleB;
        this.rsId = rsId;
    }

    /**
     * This constructor is used to take a record from an Illumina Manifest and sets the Extended-specific fields
     * in preparation for writing out the Build37ExtendedIlluminaManifestRecord to file (or otherwise using it)
     */
    Build37ExtendedIlluminaManifestRecord(final IlluminaBPMLocusEntry locusEntry,
                                          final IlluminaManifestRecord record,
                                   final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                   final Map<String, File> chainFilesMap) {
        super(record);
        this.locusEntry = locusEntry;
        this.illuminaManifestRecord = record;

        boolean stringent_validation = true;

        // Check that the fields in the bpm agree with those in the (CSV) manifest.
        validateBpmLocusEntryAgainstIlluminaManifestRecord();

        // Look for entries which Illumina has marked as invalid
        if (locusEntry.chrom.equals(IlluminaManifestRecord.ILLUMINA_FLAGGED_BAD_CHR)) {
            flag = Build37ExtendedIlluminaManifestRecord.Flag.ILLUMINA_FLAGGED;
        }
//
//        if (!r.getMajorGenomeBuild().trim().equals(BUILD_36) && !r.getMajorGenomeBuild().trim().equals(BUILD_37)) {
//            flag = Build37ExtendedIlluminaManifestRecord.Flag.UNSUPPORTED_GENOME_BUILD;
//        }

        // TODO - figure out how to do this with and without liftover files
        if (!isBad()) {
            if (illuminaManifestRecord.getMajorGenomeBuild().trim().equals(BUILD_37)) {
                // no liftover needed
                b37Chr = locusEntry.chrom;
                b37Pos = locusEntry.mapInfo;
            } else {
                liftOverToBuild37(referenceFilesMap, chainFilesMap);
            }
        }

        final ReferenceSequenceFile refFile = referenceFilesMap.get(BUILD_37);

        if (!isBad()) {
            if (isSnp()) {      // TODO - will need to do this for indels too.
                setReferenceStrand(refFile, stringent_validation);
            }
        }

        if (!isBad()) {
            if (isSnp()) {
                processSnp(refFile);
            } else {
                processIndel(refFile);
            }
        }

        if (!isBad()) {     // Note that populateSnp/IndelAlleles may flag a record as bad
            A = Allele.create(snpAlleleA, snpAlleleA.equals(snpRefAllele));
            B = Allele.create(snpAlleleB, snpAlleleB.equals(snpRefAllele));
            ref = Allele.create(snpRefAllele, true);
        } else {
            A = Allele.NO_CALL;
            B = Allele.NO_CALL;
            ref = Allele.NO_CALL;
        }
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }

    public Allele getRefAllele() {
        return ref;
    }

    public Strand getReferenceStrand() { return referenceStrand; }

    public String getB37Chr() {
        return b37Chr;
    }

    public Integer getB37Pos() {
        return b37Pos;
    }

    public String getRsId() { return rsId; }

    public Boolean isBad() {
        return flag != Flag.DUPE && flag != Flag.PASS;
    }

    public Boolean isDupe() {
        return flag == Flag.DUPE;
    }

    public Flag getFlag() {
        return flag;
    }

    /**
     * Determines the chromosome and position of the record on Build 37.
     */
    private void liftOverToBuild37(final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                   final Map<String, File> chainFilesMap) {

        final File chainFileToBuild37 = chainFilesMap.get(illuminaManifestRecord.getMajorGenomeBuild());
        final LiftOver liftOver = new LiftOver(chainFileToBuild37);
        final Interval interval = new Interval(locusEntry.getChrom(), locusEntry.getMapInfo(), locusEntry.mapInfo);
        final Interval b37Interval = liftOver.liftOver(interval);

        if (b37Interval != null) {
            b37Chr = b37Interval.getContig();
            b37Pos = b37Interval.getStart();

            // Validate that the reference allele at the lifted over coordinates matches that of the original.
            String originalRefAllele = getSequenceAt(referenceFilesMap.get(BUILD_36), locusEntry.getChrom(), locusEntry.getMapInfo(), locusEntry.mapInfo);
            String newRefAllele = getSequenceAt(referenceFilesMap.get(BUILD_37), b37Chr, b37Pos, b37Pos);
            if (originalRefAllele.equals(newRefAllele)) {
                log.debug("Lifted over record " + this);
                log.debug(" From build " + illuminaManifestRecord.getMajorGenomeBuild() + " chr=" + locusEntry.getChrom() + ", position=" + locusEntry.getMapInfo()  + " To build " + BUILD_37 + " chr=" + b37Chr + ", position=" + b37Pos);
            } else {
                flag = Flag.LIFTOVER_FAILED;
                log.error("Liftover failed for record: " + this);
                log.error( " Sequence at lifted over position does not match that at original position");
            }
        } else {
            flag = Flag.LIFTOVER_FAILED;
            log.error("Liftover failed for record: " + this);
        }
    }

    /**
     * Uses source sequence to determine snpAlleleA, snpAlleleB, and Allele.
     * <p>
     */
    private void processSnp(final ReferenceSequenceFile refFile) {
        snpAlleleA = locusEntry.getSnp().substring(1, 2);
        snpAlleleB = locusEntry.getSnp().substring(3, 4);

        if (referenceStrand == Strand.NEGATIVE) {
            snpAlleleA = SequenceUtil.reverseComplement(snpAlleleA);
            snpAlleleB = SequenceUtil.reverseComplement(snpAlleleB);
        }

        //extra validation for ambiguous snps
        // TODO - is this really needed??  Or
        if (isAmbiguous()) {
            if (illuminaManifestRecord.getAlleleBProbeSeq() != null) {
                String probeAAllele = illuminaManifestRecord.getAlleleAProbeSeq().substring(illuminaManifestRecord.getAlleleAProbeSeq().length() - 1);
                String probeBAllele = illuminaManifestRecord.getAlleleBProbeSeq().substring(illuminaManifestRecord.getAlleleBProbeSeq().length() - 1);
                if (!probeAAllele.equals(snpAlleleA) && !probeBAllele.equals(snpAlleleB) && (referenceStrand == Strand.POSITIVE)) {
                    snpAlleleA = probeAAllele;
                    snpAlleleB = probeBAllele;
                }
            } else {
                // This manifest contains no Allele B Probe Sequence.  We (currently) need this for validating/trusting
                // these ambiguous SNPs, so we are flagging it.

                // This should be an error?
                log.warn("Ambiguous probe without alleleBProbeSeq!!!  Record: " + this);
            }
        }

        // TODO - should I validate this further?  How?
        snpRefAllele = getSequenceAt(refFile, b37Chr, b37Pos, b37Pos);
    }

    private void processIndel(final ReferenceSequenceFile refFile) {
        // First let's validate the probe sequence
        String alleleAProbeSeq = illuminaManifestRecord.getAlleleAProbeSeq().toUpperCase();
        if (!ACGT_PATTERN.matcher(alleleAProbeSeq).find()) {
            throw new PicardException("AlleleAProbeSeq for record: " + this + " contains non-ACGT character(s)");
        }

        // TODO - Indel processing is doing something different for finding the refStrand than SNP processing is/was.
        Strand refStrand = illuminaManifestRecord.getRefStrand();
        if (refStrand == Strand.NONE) {
            // Some Illumina manifests do not have ref_strand defined.  We will use the illumina strand instead.
            if (getIlmnStrand() == IlluminaStrand.PLUS)  {
                refStrand = Strand.POSITIVE;
            }
            else if (getIlmnStrand() == IlluminaStrand.MINUS) {
                refStrand = Strand.NEGATIVE;
            }
        }

        if (refStrand == Strand.NEGATIVE) {
            alleleAProbeSeq = SequenceUtil.reverseComplement(alleleAProbeSeq);
        } else if (refStrand != Strand.POSITIVE) {
            flag = Flag.PROBE_SEQUENCE_STRAND_INVALID;
            log.warn("Error in processIndel.  Record:" + this);
            log.warn("  AlleleAProbeSeq on unexpected strand: " + refStrand);
            return;
        }

        // Then let's validate the source sequence
        final Matcher matcher = parseSourceSeq(getSourceSeq());
        if (isSnp()) {
            throw new PicardException("This shouldn't happen");
        }
        if (!matcher.group(PRE_INDEL_SEQUENCE).equals("-")) {       // In indels it's always of the form: [-/GCA]
            throw new PicardException("Unexpected allele '-' Record: " + this);
        }

        String fivePrimeSeq = matcher.group(FIVE_PRIME_SEQUENCE).toUpperCase();
        String indelSeq = matcher.group(INDEL_SEQUENCE).toUpperCase();
        String threePrimeSeq = matcher.group(THREE_PRIME_SEQUENCE).toUpperCase();

        if (!ACGT_PATTERN.matcher(indelSeq).find()) {
            throw new PicardException("Indel sequence for record: " + this + " contains invalid (non-ACGT) character(s)");
        }
        if (getSourceStrand() == IlluminaStrand.MINUS) {
            final String temp = threePrimeSeq;
            threePrimeSeq = SequenceUtil.reverseComplement(fivePrimeSeq);
            indelSeq = SequenceUtil.reverseComplement(indelSeq);
            fivePrimeSeq = SequenceUtil.reverseComplement(temp);
        } else if (getSourceStrand() != IlluminaStrand.PLUS) {
            flag = Flag.SOURCE_SEQUENCE_STRAND_INVALID;
            log.warn("Error in processIndel.  Record: " + this);
            log.warn("  Source Sequence on unexpected strand: " + getSourceStrand());
            return;
        }

        // Get a bounding region of sequence to search for the source and probe sequences in.
        int maxLength = Math.max(fivePrimeSeq.length(), threePrimeSeq.length());
        String regionSequence = getSequenceAt(refFile, getB37Chr(),
                getB37Pos() - maxLength - alleleAProbeSeq.length(),
                getB37Pos() + maxLength + 2 * alleleAProbeSeq.length() + 2);
        SequenceAndIndex fivePrimeSequenceAndIndex = findSubsequence(fivePrimeSeq, regionSequence);
        SequenceAndIndex threePrimeSequenceAndIndex = findSubsequence(threePrimeSeq, regionSequence);

        int indexOfProbeASeq = regionSequence.indexOf(alleleAProbeSeq);

        if (indexOfProbeASeq == -1) {
            flag = Flag.PROBE_SEQUENCE_MISMATCH;
            log.warn("Error in processIndel.  Record: " + this);
            log.warn("  Couldn't find alleleAProbeSeq in reference");
            log.debug("  AlleleAProbeSeq: " + illuminaManifestRecord.getAlleleAProbeSeq());
            return;
        }
        if ((!fivePrimeSequenceAndIndex.isFound()) || (!threePrimeSequenceAndIndex.isFound())) {
            flag = Flag.SOURCE_SEQUENCE_MISMATCH;
            log.warn("Error in processIndel.  Record: " + this);
            log.warn("  Couldn't find either 5' or 3' end of source sequence in reference");
            log.debug("  5': " + fivePrimeSeq);
            log.debug("  3': " + threePrimeSeq);
            return;
        }

        final String srcSeqWithoutVariant = fivePrimeSequenceAndIndex.sequence + threePrimeSequenceAndIndex.sequence;
        final String srcSeqWithVariant = fivePrimeSequenceAndIndex.sequence + indelSeq + threePrimeSequenceAndIndex.sequence;
        int indexWithInsertion = regionSequence.indexOf(srcSeqWithVariant);
        int indexWithoutInsertion = regionSequence.indexOf(srcSeqWithoutVariant);
        boolean isDeletion = (indexWithInsertion != -1);
        boolean isInsertion = (indexWithoutInsertion != -1);
        if ((indexWithoutInsertion == -1) && (indexWithInsertion == -1)) {
            flag = Flag.INDEL_NOT_FOUND;
            log.warn("Error in processIndel.  Record: " + this);
            log.warn("  Couldn't find source sequence with or without variant in reference");
            log.debug("  w   Variant: " + srcSeqWithoutVariant);
            log.debug("  w/o Variant: " + srcSeqWithVariant);
            return;
        }

        if (isDeletion && isInsertion) {
            flag = Flag.INDEL_CONFLICT;
            log.warn("Error in processIndel.  Record: " + this);
            log.warn("  Conflict.  Both source sequence with and without variation found in reference");
            log.debug("  w   Variant:    " + srcSeqWithVariant);
            log.debug("  w/o Variant: " + srcSeqWithoutVariant);
            return;
        }

        if (isDeletion) {
            // A deletion.  Position in VCF is before the deletion
            b37Pos--;
        }
        final String refAllele = getSequenceAt(refFile, getB37Chr(), getB37Pos(), getB37Pos());
        if (isDeletion) {
            snpRefAllele = refAllele + indelSeq;
        } else {
            snpRefAllele = refAllele;
        }
        if (locusEntry.getSnp().equals("[I/D]")) {
            snpAlleleA = refAllele + indelSeq;
            snpAlleleB = refAllele;
        } else {
            snpAlleleA = refAllele;
            snpAlleleB = refAllele + indelSeq;
        }
    }

    static SequenceAndIndex findSubsequence(final String seq, final String region) {
        SequenceAndIndex sequenceAndIndex = new SequenceAndIndex(seq, -1);

        boolean isSequenceValid = ACGT_PATTERN.matcher(seq).find();
        if (isSequenceValid) {
            sequenceAndIndex = new SequenceAndIndex(seq, region.indexOf(seq));
        } else {
            int indexOfSeq = findSubsequenceWithIupacCharacters(seq, region);
            if (indexOfSeq != -1) {
                String updatedSeq = region.substring(indexOfSeq, indexOfSeq + seq.length());
                sequenceAndIndex = new SequenceAndIndex(updatedSeq, indexOfSeq);
            }
        }
        return sequenceAndIndex;
    }

    /**
     * Given a sequence with one or more IUPAC characters in it, search for that sequence
     * (using all appropriate nucleotide substitutions) in the region sequence
     *
     * @param seq sequence containing one or more IUPAC codes
     * @param region region to search
     * @return index of first subsequence, or -1
     */
    static int findSubsequenceWithIupacCharacters(final String seq, final String region) {
        if (seq.length() > region.length()) {
            return -1;
        }

        int subSequenceIndex = -1;
        byte[] seqBytes = seq.getBytes();
        byte[] regionBytes = region.getBytes();
        for (int offset = 0; offset <= (region.length() - seq.length()); offset++) {
            boolean seqMatch = true;
            for (int index = 0; index < seq.length(); index++) {
                byte seqByte = seqBytes[index];
                byte regionByte = regionBytes[offset + index];
                if (!SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(regionByte, seqByte)) {
                    seqMatch = false;
                    break;
                }
            }
            if (seqMatch) {
                subSequenceIndex = offset;
                break;
            }
        }
        return subSequenceIndex;
    }


    /**
     * Find the sequence in a reference sequence file.
     * for the contig in the range [start,stop]
     * @param refFile ReferenceSequenceFile to use
     * @param chr Contig whose subsequence to retrieve.
     * @param startPos inclusive, 1-based start of region.
     * @param endPos inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.
     */
    private static String getSequenceAt(final ReferenceSequenceFile refFile, final String chr, final int startPos, final int endPos) {
        // TODO - should I just cache the contig lengths?
        final int contigLength = refFile.getSequenceDictionary().getSequence(chr).getSequenceLength();
        int usedEndPos = Math.min(endPos, contigLength);
        return new String(refFile.getSubsequenceAt(chr, startPos, usedEndPos).getBases()).toUpperCase();
    }

    /**
     * Use regex to capture the insertion sequence and the sequence after the indel.
     * <p>
     * The source sequence is of the from V[W/X]Y, where V,W,X,Y are sequences.
     * <p>
     * - A SNP example looks like:    AGGGAGTC[A/G]GGTTGCGA
     * V     W X    Y
     * <p>
     * - A InDel example looks like:  AGCCTCGA[-/CGAA]TCACC
     * V     W   X   Y
     */
    private static Matcher parseSourceSeq(final String sourceSeq) {
        final Matcher matcher = pattern.matcher(sourceSeq);
        if (matcher.find()) {
            return matcher;
        } else {
            throw new PicardException("Could not find the pattern V[W/X]Y in the SourceSeq: " + sourceSeq);
        }
    }

    /**
     * This method sets the reference strand
     *
     * If the refStrand is provided in the Illumina manifest(s) we will use that.
     * Unless 'stringent_validation' is in use OR the refStrand is NOT provided in the Illumina manifest we will
     * attempt to calculate the reference strand by finding the probe sequence in the reference.
     *
     * @param refFile reference to use for finding the probe sequence
     * @param stringent_validation when set, the code will do more rigorous checking of the Illumina manifest.
     */
    private void setReferenceStrand(final ReferenceSequenceFile refFile, final boolean stringent_validation) {
        if (referenceStrand != null) {
            return;
        }

        referenceStrand = locusEntry.refStrand;
        if ((referenceStrand.equals(Strand.NONE)) || (stringent_validation)) {
            // Should do a warning (at the top level?)

            final String probeSeq;
            if (isAmbiguous()) {
                //ambiguous snps contain the probed base so we need to truncate the string
                probeSeq = illuminaManifestRecord.getAlleleAProbeSeq().substring(0, illuminaManifestRecord.getAlleleAProbeSeq().length() - 1);
            } else {
                probeSeq = illuminaManifestRecord.getAlleleAProbeSeq();
            }

            // Should use lifted over b37 coordinates!!!
            final String reference = getSequenceAt(refFile, b37Chr, b37Pos - probeSeq.length(), b37Pos - 1);
            final String reverseReference = SequenceUtil.reverseComplement(getSequenceAt(refFile, b37Chr, b37Pos + 1, b37Pos + probeSeq.length()));

            if (reference.equals(probeSeq)) {
                referenceStrand = Strand.POSITIVE;
            } else if (reverseReference.equals(probeSeq)) {
                referenceStrand = Strand.NEGATIVE;
            } else {
                flag = Flag.PROBE_SEQUENCE_MISMATCH;
                log.warn("Error in getStrand.  Record:" + this);
                log.warn("  Couldn't find alleleAProbeSeq in reference");
                log.debug("  AlleleAProbeSeq: " + illuminaManifestRecord.getAlleleAProbeSeq());
                log.debug("  Reference:       " + reference);
                log.debug("  Reverse Ref:     " + reverseReference);
                return;
            }
            if ((!locusEntry.refStrand.equals(Strand.NONE)) && (!referenceStrand.equals(locusEntry.refStrand))) {
                // This should be an error.
                throw new PicardException("Calculated Reference Strand differs from Reference Strand provided in Illumina Manifest");
            }
        }
    }

    private String getChrom() {
        // TODO - previously I did a '.trim()
        return locusEntry.chrom;
    }

    public void setRsId(String rsId) {
        this.rsId = rsId;
    }

    // TODO - you should have two different flags for error and dupe and only merge them later.
    public void setDupe(boolean isDupe) {
        if (!isBad()) {
            if (isDupe) {
                flag = Flag.DUPE;
            }
        }
    }


    // TODO - need an object / class to sanitize these and then pull the correct record???
    private void validateBpmLocusEntryAgainstIlluminaManifestRecord() {
        validateEntryField(locusEntry.ilmnId, illuminaManifestRecord.getIlmnId(), "ilmnId");
        validateEntryField(locusEntry.name, illuminaManifestRecord.getName(), "name");
        validateEntryField(locusEntry.ilmnStrand, illuminaManifestRecord.getIlmnStrand(), "ilmnStrand");
        validateEntryField(locusEntry.snp, illuminaManifestRecord.getSnp(), "snp");
        validateEntryField(locusEntry.chrom, illuminaManifestRecord.getChr(), "chrom");
        validateEntryField(locusEntry.ploidy, illuminaManifestRecord.getPloidy(), "ploidy");
        validateEntryField(locusEntry.species, illuminaManifestRecord.getSpecies(), "species");
        validateEntryField(locusEntry.mapInfo, illuminaManifestRecord.getPosition(), "mapInfo");
        validateEntryField(locusEntry.addressA, Integer.parseInt(illuminaManifestRecord.getAddressAId()), "addressAId");
        if (locusEntry.getVersion() == 4) {
            validateEntryField(locusEntry.alleleAProbeSeq, illuminaManifestRecord.getAlleleAProbeSeq(), "alleleAProbeSeq");
        }
        if ((locusEntry.addressB != -1) && (illuminaManifestRecord.getAddressBId() != null)) {
            validateEntryField(locusEntry.addressB, Integer.parseInt(illuminaManifestRecord.getAddressBId()), "addressBId");
        }
        if (locusEntry.getVersion() == 4) {
            validateEntryField(locusEntry.alleleBProbeSeq, illuminaManifestRecord.getAlleleBProbeSeq(), "alleleBProbeSeq");
        }
        validateEntryField(locusEntry.genomeBuild, illuminaManifestRecord.getGenomeBuild(), "genomeBuild");
        validateEntryField(locusEntry.source, illuminaManifestRecord.getSource(), "source");
        validateEntryField(locusEntry.sourceVersion, illuminaManifestRecord.getSourceVersion(), "sourceVersion");
        validateEntryField(locusEntry.sourceStrand, illuminaManifestRecord.getSourceStrand(), "sourceStrand");
        if (locusEntry.getVersion() == 4) {
            validateEntryField(locusEntry.sourceSeq, illuminaManifestRecord.getSourceSeq(), "sourceSeq");
            validateEntryField(locusEntry.topGenomicSeq, illuminaManifestRecord.getTopGenomicSeq(), "topGenomicSeq");
        }
//        if (record.getExpClusters() != null) {
//            validateEntryField(locusEntry.expClusters, Integer.parseInt(record.getExpClusters()), "expClusters");
//        }
//        validateEntryField(locusEntry.intensityOnly, record.getIntensityOnly(), "intensityOnly");
        if (locusEntry.getVersion() == 8) {
            validateEntryField(locusEntry.refStrand, illuminaManifestRecord.getRefStrand(), "refStrand");
        }
    }

    private void validateEntryField(final Object locusEntryField, final Object recordField, final String fieldName) {
        // TODO - should this just be a non-passing entry?  Probably so.
        if (!locusEntryField.equals(recordField)) {
            throw new PicardException("Field '" + fieldName + "' disagrees between BPM file (found '" + locusEntryField + "') and CSV (found: '" + recordField + "')");
        }
    }

    @Override
    public String getLine() {
        final String originalLine = super.getLine();

        final List<String> extensions = new ArrayList<>();
        extensions.add(b37Chr);
        extensions.add(b37Pos != null ? b37Pos.toString() : null);
        extensions.add(snpRefAllele);
        extensions.add(snpAlleleA);
        extensions.add(snpAlleleB);
        extensions.add(rsId);
        extensions.add(flag.name());

        return originalLine + "," + StringUtils.join(extensions, ",");
    }

    static class SequenceAndIndex {
        final String sequence;
        final int index;

        SequenceAndIndex(String sequence, int index) {
            this.sequence = sequence;
            this.index = index;
        }

        final boolean isFound() {
            return this.index != -1;
        }
    }
}
