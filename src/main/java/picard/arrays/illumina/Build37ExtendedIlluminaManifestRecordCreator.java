package picard.arrays.illumina;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import picard.PicardException;

import java.io.File;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Build37ExtendedIlluminaManifestRecordCreator {

    private final Map<String, ReferenceSequenceFile> referenceFilesMap;
    private final Map<String, File> chainFilesMap;
    // TODO - need to handle this for indels too!
    private boolean refStrandDefinedInManifest = true;

    private final Log log = Log.getInstance(Build37ExtendedIlluminaManifestRecordCreator.class);

    public static final String BUILD_36 = "36";
    public static final String BUILD_37 = "37";

    // These are the IUPAC nucleotide codes as described here: https://www.bioinformatics.org/sms/iupac.html
    public static final String IUPAC_NUCLEOTIDE_CODES = "ACGTRYSWKMBDHVNacgtryswkmbdhvn";
    public static final String ACGT_CODES = "ACGTacgt";

    private static final String SRC_SEQ_REGEX = "([" + IUPAC_NUCLEOTIDE_CODES + "]*)\\[([" + ACGT_CODES + "-])\\/([" + ACGT_CODES + "]*)\\]([" + IUPAC_NUCLEOTIDE_CODES + "]*)";
    private static final Pattern pattern = Pattern.compile(SRC_SEQ_REGEX);

    private static final String ACGT_REGEX = "^[" + ACGT_CODES + "]+$";
    private static final Pattern ACGT_PATTERN = Pattern.compile(ACGT_REGEX);

    // Symbolics for the regex groups...
    public static final int FIVE_PRIME_SEQUENCE = 1;
    public static final int PRE_INDEL_SEQUENCE = 2;
    public static final int INDEL_SEQUENCE = 3;
    public static final int THREE_PRIME_SEQUENCE = 4;

    Build37ExtendedIlluminaManifestRecordCreator(final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                          final Map<String, File> chainFilesMap) {
        this.referenceFilesMap = referenceFilesMap;
        this.chainFilesMap = chainFilesMap;
    }

    public Build37ExtendedIlluminaManifestRecord createRecord(final IlluminaBPMLocusEntry locusEntry,
                                                              final IlluminaManifestRecord illuminaManifestRecord) {

        // Check that the fields in the bpm agree with those in the (CSV) manifest.
        validateBpmLocusEntryAgainstIlluminaManifestRecord(locusEntry, illuminaManifestRecord);

        // Should I create the Build37ExtendedIlluminaManifestRecord here and set all the things in it?  Probably.
        Build37ExtendedIlluminaManifestRecord newRecord = new Build37ExtendedIlluminaManifestRecord(illuminaManifestRecord,
                Build37ExtendedIlluminaManifestRecord.Flag.PASS,
                "",
                null,
                "",
                "",
                "",
                "");

        // Look for entries which Illumina has marked as invalid
        if (locusEntry.chrom.equals(IlluminaManifestRecord.ILLUMINA_FLAGGED_BAD_CHR)) {
            newRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.ILLUMINA_FLAGGED;
            return newRecord;
        }

        // TODO - TEMP!  For testing sake I am failing all indels - only want to focus on SNPs now.
//        if (newRecord.isIndel()) {
//            newRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.ILLUMINA_FLAGGED;
//            return newRecord;
//        }
        // End TODO - TEMP!

        if (!illuminaManifestRecord.getMajorGenomeBuild().trim().equals(BUILD_36) &&
                !illuminaManifestRecord.getMajorGenomeBuild().trim().equals(BUILD_37)) {
            newRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.UNSUPPORTED_GENOME_BUILD;
            return newRecord;
        }

        // TODO - figure out how to do this with and without liftover files
        if (illuminaManifestRecord.getMajorGenomeBuild().trim().equals(BUILD_37)) {
            // no liftover needed
            newRecord.b37Chr = locusEntry.chrom;
            newRecord.b37Pos = locusEntry.mapInfo;
        } else {
            liftOverToBuild37(newRecord, locusEntry, illuminaManifestRecord);
            if (newRecord.isFail()) {
                return newRecord;
            }
        }

        final ReferenceSequenceFile refFile = referenceFilesMap.get(BUILD_37);

        if (newRecord.isSnp()) {      // TODO - will need to do this for indels too.
            setReferenceStrandForSnp(newRecord, locusEntry, illuminaManifestRecord, refFile);
        }

        if (!newRecord.isFail()) {
            if (newRecord.isSnp()) {
                processSnp(newRecord, locusEntry, illuminaManifestRecord, refFile);
            } else {
                processIndel(newRecord, locusEntry, illuminaManifestRecord, refFile);
            }
        }

        return newRecord;
    }

    /**
     * Uses source sequence to determine snpAlleleA, snpAlleleB, and Allele.
     * <p>
     */
    private void processSnp(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                            final IlluminaBPMLocusEntry locusEntry,
                            final IlluminaManifestRecord illuminaManifestRecord,
                            final ReferenceSequenceFile refFile) {
        build37ExtendedIlluminaManifestRecord.snpAlleleA = locusEntry.getSnp().substring(1, 2);
        build37ExtendedIlluminaManifestRecord.snpAlleleB = locusEntry.getSnp().substring(3, 4);

        if (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.NEGATIVE) {
            build37ExtendedIlluminaManifestRecord.snpAlleleA = SequenceUtil.reverseComplement(build37ExtendedIlluminaManifestRecord.snpAlleleA);
            build37ExtendedIlluminaManifestRecord.snpAlleleB = SequenceUtil.reverseComplement(build37ExtendedIlluminaManifestRecord.snpAlleleB);
        }

        if (build37ExtendedIlluminaManifestRecord.isAmbiguous()) {
            if (illuminaManifestRecord.getAlleleBProbeSeq() != null) {
                String probeAAllele = illuminaManifestRecord.getAlleleAProbeSeq().substring(illuminaManifestRecord.getAlleleAProbeSeq().length() - 1);
                String probeBAllele = illuminaManifestRecord.getAlleleBProbeSeq().substring(illuminaManifestRecord.getAlleleBProbeSeq().length() - 1);
                if (!probeAAllele.equals(build37ExtendedIlluminaManifestRecord.snpAlleleA) &&
                    !probeBAllele.equals(build37ExtendedIlluminaManifestRecord.snpAlleleB) &&
                        (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.POSITIVE)) {
                    build37ExtendedIlluminaManifestRecord.snpAlleleA = probeAAllele;
                    build37ExtendedIlluminaManifestRecord.snpAlleleB = probeBAllele;
                }
            } else {
                // This manifest contains no Allele B Probe Sequence.  We (currently) need this for validating/trusting
                // these ambiguous SNPs, so we are flagging it.
                build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.MISSING_ALLELE_B_PROBESEQ;
                log.warn("Error in processSnp.  Record:" + build37ExtendedIlluminaManifestRecord);
                log.warn("  Ambiguous probe without alleleBProbeSeq");
            }
        }

        // TODO - should I validate this further?  How?
        build37ExtendedIlluminaManifestRecord.snpRefAllele = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos, build37ExtendedIlluminaManifestRecord.b37Pos);
    }

    /**
     * This method sets the reference strand
     *
     * If the refStrand is provided in the Illumina manifest(s) we will use that.
     * Unless 'stringent_validation' is in use OR the refStrand is NOT provided in the Illumina manifest we will
     * attempt to calculate the reference strand by finding the probe sequence in the reference.
     *
     * @param refFile reference to use for finding the probe sequence
     */
    private void setReferenceStrandForSnp(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                                    final IlluminaBPMLocusEntry locusEntry,
                                    final IlluminaManifestRecord illuminaManifestRecord,
                                    final ReferenceSequenceFile refFile) {
        if (build37ExtendedIlluminaManifestRecord.referenceStrand != null) {
            return;
        }

        build37ExtendedIlluminaManifestRecord.referenceStrand = locusEntry.refStrand;
        if (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.NONE) {
            refStrandDefinedInManifest = false;
        }
        // TODO - NOTE: Bug2/Fix2 - I am using the Illumina-supplied refStrand if available.
        if (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.NONE) {
            // Begin Calculate the Freeseq way:
//            String alleleA = build37ExtendedIlluminaManifestRecord.snpAlleleA;
//            String alleleB = build37ExtendedIlluminaManifestRecord.snpAlleleB;
//            if (build37ExtendedIlluminaManifestRecord.getIlmnStrand() == IlluminaManifestRecord.IlluminaStrand.BOT) {
//                alleleA = SequenceUtil.reverseComplement(alleleA);
//                alleleB = SequenceUtil.reverseComplement(alleleB);
//            } else if (build37ExtendedIlluminaManifestRecord.getIlmnStrand() != IlluminaManifestRecord.IlluminaStrand.TOP) {
//                // TODO - make this a generic class of error
//                throw new PicardException("Unexpected Illumina strand value: " + build37ExtendedIlluminaManifestRecord.getIlmnStrand());
//            }
//            strand = get_strand_from_top_alleles(alleleA, alleleB, ref, win, len);

            // End calculate the freeseq way

            final String probeSeq;
            if (build37ExtendedIlluminaManifestRecord.isAmbiguous()) {
                //ambiguous snps contain the probed base so we need to truncate the string
                probeSeq = illuminaManifestRecord.getAlleleAProbeSeq().substring(0, illuminaManifestRecord.getAlleleAProbeSeq().length() - 1);
            } else {
                probeSeq = illuminaManifestRecord.getAlleleAProbeSeq();
            }

            final String reference = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos - probeSeq.length(), build37ExtendedIlluminaManifestRecord.b37Pos - 1);
            final String reverseReference = SequenceUtil.reverseComplement(getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos + 1, build37ExtendedIlluminaManifestRecord.b37Pos + probeSeq.length()));

            if (reference.equals(probeSeq)) {
                build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.POSITIVE;
            } else if (reverseReference.equals(probeSeq)) {
                build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.NEGATIVE;
            } else {
                build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.PROBE_SEQUENCE_MISMATCH;
                log.warn("Error in getStrand.  Record:" + build37ExtendedIlluminaManifestRecord);
                log.warn("  Couldn't find alleleAProbeSeq in reference");
                log.debug("  AlleleAProbeSeq: " + illuminaManifestRecord.getAlleleAProbeSeq());
                log.debug("  Reference:       " + reference);
                log.debug("  Reverse Ref:     " + reverseReference);
            }
//            if ((!locusEntry.refStrand.equals(Strand.NONE)) && (!build37ExtendedIlluminaManifestRecord.referenceStrand.equals(locusEntry.refStrand))) {
//                build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.CALC_REF_STRAND_MISMATCH;
//                log.warn("Error in getStrand.  Record:" + build37ExtendedIlluminaManifestRecord);
//                log.warn("  The calculated refStrand (" + build37ExtendedIlluminaManifestRecord.referenceStrand + ") differs from that specified in the manifest (" + locusEntry.refStrand + ")");
//            }
        }
    }

    /**
     * This method sets the reference strand
     * TODO - NOTE - This is a version that maintains compatibility with 1.6 - it does NOT use the Illumina-supplied
     *               Reference Strand.  It calculates it like we did in 1.5 and before (picard-private)
     *
     * If the refStrand is provided in the Illumina manifest(s) we will use that.
     * Unless 'stringent_validation' is in use OR the refStrand is NOT provided in the Illumina manifest we will
     * attempt to calculate the reference strand by finding the probe sequence in the reference.
     *
     * @param refFile reference to use for finding the probe sequence
     */
    private void setReferenceStrandForSnp1_6(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                                          final IlluminaBPMLocusEntry locusEntry,
                                          final IlluminaManifestRecord illuminaManifestRecord,
                                          final ReferenceSequenceFile refFile) {
        if (build37ExtendedIlluminaManifestRecord.referenceStrand != null) {
            return;
        }

        build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.NONE;
        if (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.NONE) {
            refStrandDefinedInManifest = false;
        }

        final String probeSeq;
        if (build37ExtendedIlluminaManifestRecord.isAmbiguous()) {
            //ambiguous snps contain the probed base so we need to truncate the string
            probeSeq = illuminaManifestRecord.getAlleleAProbeSeq().substring(0, illuminaManifestRecord.getAlleleAProbeSeq().length() - 1);
        } else {
            probeSeq = illuminaManifestRecord.getAlleleAProbeSeq();
        }

        final String reference = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos - probeSeq.length(), build37ExtendedIlluminaManifestRecord.b37Pos - 1);
        final String reverseReference = SequenceUtil.reverseComplement(getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos + 1, build37ExtendedIlluminaManifestRecord.b37Pos + probeSeq.length()));

        if (reference.equals(probeSeq)) {
            build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.POSITIVE;
        } else if (reverseReference.equals(probeSeq)) {
            build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.NEGATIVE;
        } else {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.PROBE_SEQUENCE_MISMATCH;
            log.warn("Error in getStrand.  Record:" + build37ExtendedIlluminaManifestRecord);
            log.warn("  Couldn't find alleleAProbeSeq in reference");
            log.debug("  AlleleAProbeSeq: " + illuminaManifestRecord.getAlleleAProbeSeq());
            log.debug("  Reference:       " + reference);
            log.debug("  Reverse Ref:     " + reverseReference);
        }
    }

    private void processIndel(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                              final IlluminaBPMLocusEntry locusEntry,
                              final IlluminaManifestRecord illuminaManifestRecord,
                              final ReferenceSequenceFile refFile) {
        // First let's validate the probe sequence
        String alleleAProbeSeq = illuminaManifestRecord.getAlleleAProbeSeq().toUpperCase();
        if (!ACGT_PATTERN.matcher(alleleAProbeSeq).find()) {
            throw new PicardException("AlleleAProbeSeq for record: " + this + " contains non-ACGT character(s)");
        }

        // TODO - Indel processing is doing something different for finding the refStrand than SNP processing is/was.
        Strand refStrand = illuminaManifestRecord.getRefStrand();
        if (refStrand == Strand.NONE) {
            // Some Illumina manifests do not have ref_strand defined.  We will use the illumina strand instead.
            if (build37ExtendedIlluminaManifestRecord.getIlmnStrand() == IlluminaManifestRecord.IlluminaStrand.PLUS)  {
                refStrand = Strand.POSITIVE;
            }
            else if (build37ExtendedIlluminaManifestRecord.getIlmnStrand() == IlluminaManifestRecord.IlluminaStrand.MINUS) {
                refStrand = Strand.NEGATIVE;
            }
        }

        if (refStrand == Strand.NEGATIVE) {
            alleleAProbeSeq = SequenceUtil.reverseComplement(alleleAProbeSeq);
        } else if (refStrand != Strand.POSITIVE) {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.PROBE_SEQUENCE_STRAND_INVALID;
            log.warn("Error in processIndel.  Record:" + build37ExtendedIlluminaManifestRecord);
            log.warn("  AlleleAProbeSeq on unexpected strand: " + refStrand);
            return;
        }

        // Then let's validate the source sequence
        final Matcher matcher = parseSourceSeq(build37ExtendedIlluminaManifestRecord.getSourceSeq());
        if (build37ExtendedIlluminaManifestRecord.isSnp()) {
            throw new PicardException("This shouldn't happen");
        }
        if (!matcher.group(PRE_INDEL_SEQUENCE).equals("-")) {       // In indels it's always of the form: [-/GCA]
            throw new PicardException("Unexpected allele '-' Record: " + build37ExtendedIlluminaManifestRecord);
        }

        String fivePrimeSeq = matcher.group(FIVE_PRIME_SEQUENCE).toUpperCase();
        String indelSeq = matcher.group(INDEL_SEQUENCE).toUpperCase();
        String threePrimeSeq = matcher.group(THREE_PRIME_SEQUENCE).toUpperCase();

        if (!ACGT_PATTERN.matcher(indelSeq).find()) {
            throw new PicardException("Indel sequence for record: " + build37ExtendedIlluminaManifestRecord + " contains invalid (non-ACGT) character(s)");
        }
        if (build37ExtendedIlluminaManifestRecord.getSourceStrand() == IlluminaManifestRecord.IlluminaStrand.MINUS) {
            final String temp = threePrimeSeq;
            threePrimeSeq = SequenceUtil.reverseComplement(fivePrimeSeq);
            indelSeq = SequenceUtil.reverseComplement(indelSeq);
            fivePrimeSeq = SequenceUtil.reverseComplement(temp);
        } else if (build37ExtendedIlluminaManifestRecord.getSourceStrand() != IlluminaManifestRecord.IlluminaStrand.PLUS) {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.SOURCE_SEQUENCE_STRAND_INVALID;
            log.warn("Error in processIndel.  Record: " + build37ExtendedIlluminaManifestRecord);
            log.warn("  Source Sequence on unexpected strand: " + build37ExtendedIlluminaManifestRecord.getSourceStrand());
            return;
        }

        // Get a bounding region of sequence to search for the source and probe sequences in.
        int maxLength = Math.max(fivePrimeSeq.length(), threePrimeSeq.length());
        String regionSequence = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.getB37Chr(),
                build37ExtendedIlluminaManifestRecord.getB37Pos() - maxLength - alleleAProbeSeq.length(),
                build37ExtendedIlluminaManifestRecord.getB37Pos() + maxLength + 2 * alleleAProbeSeq.length() + 2);
        SequenceAndIndex fivePrimeSequenceAndIndex = findSubsequence(fivePrimeSeq, regionSequence);
        SequenceAndIndex threePrimeSequenceAndIndex = findSubsequence(threePrimeSeq, regionSequence);

        int indexOfProbeASeq = regionSequence.indexOf(alleleAProbeSeq);

        if (indexOfProbeASeq == -1) {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.PROBE_SEQUENCE_MISMATCH;
            log.warn("Error in processIndel.  Record: " + build37ExtendedIlluminaManifestRecord);
            log.warn("  Couldn't find alleleAProbeSeq in reference");
            log.debug("  AlleleAProbeSeq: " + illuminaManifestRecord.getAlleleAProbeSeq());
            return;
        }
        if ((!fivePrimeSequenceAndIndex.isFound()) || (!threePrimeSequenceAndIndex.isFound())) {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.SOURCE_SEQUENCE_MISMATCH;
            log.warn("Error in processIndel.  Record: " + build37ExtendedIlluminaManifestRecord);
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
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.INDEL_NOT_FOUND;
            log.warn("Error in processIndel.  Record: " + build37ExtendedIlluminaManifestRecord);
            log.warn("  Couldn't find source sequence with or without variant in reference");
            log.debug("  w   Variant: " + srcSeqWithoutVariant);
            log.debug("  w/o Variant: " + srcSeqWithVariant);
            return;
        }

        if (isDeletion && isInsertion) {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.INDEL_CONFLICT;
            log.warn("Error in processIndel.  Record: " + build37ExtendedIlluminaManifestRecord);
            log.warn("  Conflict.  Both source sequence with and without variation found in reference");
            log.debug("  w   Variant:    " + srcSeqWithVariant);
            log.debug("  w/o Variant: " + srcSeqWithoutVariant);
            return;
        }

        if (isDeletion) {
            // A deletion.  Position in VCF is before the deletion
            build37ExtendedIlluminaManifestRecord.b37Pos--;
        }
        final String refAllele = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.getB37Chr(), build37ExtendedIlluminaManifestRecord.getB37Pos(), build37ExtendedIlluminaManifestRecord.getB37Pos());
        if (isDeletion) {
            build37ExtendedIlluminaManifestRecord.snpRefAllele = refAllele + indelSeq;
        } else {
            build37ExtendedIlluminaManifestRecord.snpRefAllele = refAllele;
        }
        if (locusEntry.getSnp().equals("[I/D]")) {
            build37ExtendedIlluminaManifestRecord.snpAlleleA = refAllele + indelSeq;
            build37ExtendedIlluminaManifestRecord.snpAlleleB = refAllele;
        } else {
            build37ExtendedIlluminaManifestRecord.snpAlleleA = refAllele;
            build37ExtendedIlluminaManifestRecord.snpAlleleB = refAllele + indelSeq;
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

    // TODO - need an object / class to sanitize these and then pull the correct record???
    private void validateBpmLocusEntryAgainstIlluminaManifestRecord(final IlluminaBPMLocusEntry locusEntry,
                                                                    final IlluminaManifestRecord illuminaManifestRecord) {
        validateEntryField(locusEntry.ilmnId, illuminaManifestRecord.getIlmnId(), "ilmnId");
        validateEntryField(locusEntry.name, illuminaManifestRecord.getName(), "name");
        validateEntryField(locusEntry.ilmnStrand, illuminaManifestRecord.getIlmnStrand(), "ilmnStrand");
        validateEntryField(locusEntry.snp, illuminaManifestRecord.getSnp(), "snp");
        validateEntryField(locusEntry.chrom, illuminaManifestRecord.getChr(), "chrom");
        validateEntryField(locusEntry.ploidy, illuminaManifestRecord.getPloidy(), "ploidy");
        // Commented out because a PsychChip had an entry where in the bpm species was 'other' and in the csv it was 'Homo Sapiens')
//        validateEntryField(locusEntry.species, illuminaManifestRecord.getSpecies(), "species");
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

    /**
     * Determines the chromosome and position of the record on Build 37.
     */
    private void liftOverToBuild37(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                                   final IlluminaBPMLocusEntry locusEntry,
                                   final IlluminaManifestRecord illuminaManifestRecord) {

        final File chainFileToBuild37 = chainFilesMap.get(illuminaManifestRecord.getMajorGenomeBuild());
        final LiftOver liftOver = new LiftOver(chainFileToBuild37);
        final Interval interval = new Interval(locusEntry.getChrom(), locusEntry.getMapInfo(), locusEntry.mapInfo);
        final Interval b37Interval = liftOver.liftOver(interval);

        if (b37Interval != null) {
            build37ExtendedIlluminaManifestRecord.b37Chr = b37Interval.getContig();
            build37ExtendedIlluminaManifestRecord.b37Pos = b37Interval.getStart();

            // Validate that the reference allele at the lifted over coordinates matches that of the original.
            String originalRefAllele = getSequenceAt(referenceFilesMap.get(BUILD_36), locusEntry.getChrom(), locusEntry.getMapInfo(), locusEntry.mapInfo);
            String newRefAllele = getSequenceAt(referenceFilesMap.get(BUILD_37), build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos, build37ExtendedIlluminaManifestRecord.b37Pos);
            if (originalRefAllele.equals(newRefAllele)) {
                log.debug("Lifted over record " + build37ExtendedIlluminaManifestRecord);
                log.debug(" From build " + illuminaManifestRecord.getMajorGenomeBuild() +
                        " chr=" + locusEntry.getChrom() +
                        ", position=" + locusEntry.getMapInfo()  +
                        " To build " + BUILD_37 +
                        " chr=" + build37ExtendedIlluminaManifestRecord.b37Chr + ", position=" + build37ExtendedIlluminaManifestRecord.b37Pos);
            } else {
                build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.LIFTOVER_FAILED;
                log.error("Liftover failed for record: " + build37ExtendedIlluminaManifestRecord);
                log.error( " Sequence at lifted over position does not match that at original position");
            }
        } else {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.LIFTOVER_FAILED;
            log.error("Liftover failed for record: " + build37ExtendedIlluminaManifestRecord);
        }
    }

    public boolean isRefStrandDefinedInManifest() {
        return refStrandDefinedInManifest;
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
