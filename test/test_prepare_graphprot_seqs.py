# from re import search
from scripttest import TestFileEnvironment

bindir = "bin/"
script = "prepare_graphprot_seqs.py"
# test file environment
testdir = "test/testenv_prepare_graphprot_seqs"
# directories relative to test file environment
bindir_rel = "../../" + bindir
datadir_rel = "../../test/"

env = TestFileEnvironment(testdir)


def test_no_such_genome_id():
    """Try using an unknown genome id for automatic chromosome limits extraction."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords.bed",
        "unknown_genome",
        datadir_rel + "artificial_genome.fa",
        expect_error=True
    )
    # assert search("chromosome_limits", run.stdout), "Error message did not contain reference to 'chromosome_limits', was :'\n{}'".format(run.stdout)
    assert(run.returncode != 0)


def test_stranded_shuffling_plus():
    """Test if negatives are shuffled on the correct strand. Shuffling plus strand peaks on a minus strand region should result in zero negatives."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords_onlyplus.bed",
        "unknown_genome",
        datadir_rel + "artificial_genome.fa",
        "--chromosome_limits", datadir_rel + "artificial_genome.limits",
        "--seq_length", "20",
        "--core_length", "10",
        "--output_file_prefix", "test_stranded_shuffling_plus",
        "--negative_site_candidate_regions_fn", datadir_rel + "artificial_candidate_region_minus.bed",
        "-v",
        expect_stderr=True,
    )
    assert "derived negative cores: 0" in run.stderr, "Error, expecting zero negative instances."


def test_stranded_shuffling_minus():
    """Test if negatives are shuffled on the correct strand. Shuffling minus strand peaks on a plus strand region should result in zero negatives."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords_onlyminus.bed",
        "unknown_genome",
        datadir_rel + "artificial_genome.fa",
        "--chromosome_limits", datadir_rel + "artificial_genome.limits",
        "--seq_length", "20",
        "--core_length", "10",
        "--output_file_prefix", "test_stranded_shuffling_minus",
        "--negative_site_candidate_regions_fn", datadir_rel + "artificial_candidate_region_plus.bed",
        "-v",
        expect_stderr=True,
    )
    assert "derived negative cores: 0" in run.stderr, "Error, expecting zero negative instances." + run.stderr + run.stdout


def test_stranded_shuffling_on_plus():
    """Test if negatives are shuffled on the correct strand. Shuffling minus two plus strand peaks and one minus strand peak on a plus strand region should result in two negatives."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords.bed",
        "unknown_genome",
        datadir_rel + "artificial_genome.fa",
        "--chromosome_limits", datadir_rel + "artificial_genome.limits",
        "--seq_length", "20",
        "--core_length", "10",
        "--output_file_prefix", "test_stranded_shuffling_on_plus",
        "--negative_site_candidate_regions_fn", datadir_rel + "artificial_candidate_region_plus.bed",
        "-v",
        expect_stderr=True,
    )
    assert "derived negative cores: 2" in run.stderr, "Error, expecting two negative instances." + run.stdout + run.stderr


def test_stranded_shuffling_minus_genomedef():
    """Test if negatives are shuffled on the correct strand. Shuffling minus strand peaks on a plus strand region should result in zero negatives."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords_onlyminus.bed",
        "hg19",
        datadir_rel + "artificial_genome.fa",
        "--seq_length", "20",
        "--core_length", "10",
        "--output_file_prefix", "test_stranded_shuffling_minus",
        "--negative_site_candidate_regions_fn", datadir_rel + "artificial_candidate_region_plus.bed",
        "-v",
        expect_stderr=True,
    )
    assert "derived negative cores: 0" in run.stderr, "Error, expecting zero negative instances." + run.stderr + run.stdout


def test_stranded_shuffling_on_plus_genomedef():
    """Test if negatives are shuffled on the correct strand. Shuffling minus two plus strand peaks and one minus strand peak on a plus strand region should result in two negatives."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords.bed",
        "hg19",
        datadir_rel + "artificial_genome.fa",
        "--seq_length", "20",
        "--core_length", "10",
        "--output_file_prefix", "test_stranded_shuffling_on_plus",
        "--negative_site_candidate_regions_fn", datadir_rel + "artificial_candidate_region_plus.bed",
        "-v",
        expect_stderr=True,
    )
    assert "derived negative cores: 2" in run.stderr, "Error, expecting two negative instances." + run.stdout + run.stderr


def test_stranded_shuffling_on_minus():
    """Test if negatives are shuffled on the correct strand. Shuffling minus two plus strand peaks and one minus strand peak on a minus strand region should result in one negative."""
    run = env.run(
        bindir_rel + script,
        datadir_rel + "artificial_coords.bed",
        "unknown_genome",
        datadir_rel + "artificial_genome.fa",
        "--chromosome_limits", datadir_rel + "artificial_genome.limits",
        "--seq_length", "20",
        "--core_length", "10",
        "--output_file_prefix", "test_stranded_shuffling_on_minus",
        "--negative_site_candidate_regions_fn", datadir_rel + "artificial_candidate_region_minus.bed",
        "-v",
        expect_stderr=True,
    )
    assert "derived negative cores: 1" in run.stderr, "Error, expecting one negative instance." + run.stdout + run.stderr
