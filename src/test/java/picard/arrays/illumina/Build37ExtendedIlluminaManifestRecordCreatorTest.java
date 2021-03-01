package picard.arrays.illumina;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class Build37ExtendedIlluminaManifestRecordCreatorTest {

    @DataProvider(name = "findSubsequenceDataProvider")
    public Object[][] findSubsequenceDataProvider() {
        return new Object[][]{
                {"ANCT", "CAGCTT", 1, "AGCT"},
                {"ANCT", "CGGCTT", -1, "ANCT"},
                {"GCTT", "CGGCTT", 2, "GCTT"},
                {"GCTT", "CGGCTTCGGCTT", 2, "GCTT"},
                {"GCNT", "CGGCTT", 2, "GCTT"},
                {"GCTN", "CGGCTT", 2, "GCTT"},
                {"NAGGATTC", "CAGGATTC", 0, "CAGGATTC"},
                {"CNGGATTC", "CAGGATTC", 0, "CAGGATTC"},
                {"CANGATTC", "CAGGATTC", 0, "CAGGATTC"},
                {"CAGNATTC", "CAGGATTC", 0, "CAGGATTC"},
                {"CAGGNTTC", "CAGGATTC", 0, "CAGGATTC"},
                {"CAGGANTC", "CAGGATTC", 0, "CAGGATTC"},
                {"CAGGATNC", "CAGGATTC", 0, "CAGGATTC"},
                {"CAGGATTN", "CAGGATTC", 0, "CAGGATTC"},
                {"NNATT", "CAGGATTC", 2, "GGATT"},
                {"GNTN", "CAGGATTC", 3, "GATT"}
        };
    }

    @Test(dataProvider = "findSubsequenceDataProvider")
    public void testFindSubsequenceWithIupacCharacters(final String seq, final String region, final int expectedIndex, final String expectedSequence) {
        int index = Build37ExtendedIlluminaManifestRecordCreator.findSubsequenceWithIupacCharacters(seq, region);
        Assert.assertEquals(index, expectedIndex);
    }

    @Test(dataProvider = "findSubsequenceDataProvider")
    public void testFindSubsequence(final String seq, final String region, final int expectedIndex, final String expectedSequence) {
        Build37ExtendedIlluminaManifestRecordCreator.SequenceAndIndex sequenceAndIndex = Build37ExtendedIlluminaManifestRecordCreator.findSubsequence(seq, region);
        Assert.assertEquals(sequenceAndIndex.index, expectedIndex);
        Assert.assertEquals(sequenceAndIndex.sequence, expectedSequence);
        Assert.assertEquals((sequenceAndIndex.index != -1), sequenceAndIndex.isFound());
    }
}

