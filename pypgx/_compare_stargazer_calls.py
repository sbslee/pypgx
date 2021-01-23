import pandas as pd
from pypgx.sdk import Results

def compare_stargazer_calls(ref_file, test_file, output_file=None):
    """Compute the concordance between two 'genotype-calls.tsv' files
    created by Stargazer.

    Parameters
    ----------
    ref_file : str
        Path to the reference or truth 'genotype-calls.tsv' file created by
        Stargazer.
    test_file : str
        Path to the test 'genotype-calls.tsv' file created by Stargazer.
    output_file : str, optional
        Path to the output file.

    Returns
    -------
    Results
        Results instance which has the following attributes: ``ref_file``,
        ``ref_samples_total``, ``test_file``, ``test_samples_total``,
        ``overlap_samples_total``, ``concordant_samples_total``,
        ``concordant_percentage``.

    """
    rf = pd.read_table(ref_file, index_col=0)
    tf = pd.read_table(test_file, index_col=0)

    rf['ref-dip'] = rf['hap1_main'] + '/' + rf['hap2_main']
    tf['test-dip'] = tf['hap1_main'] + '/' + tf['hap2_main']

    df = pd.concat([rf, tf], axis=1, join='inner')
    concordant = sum(df['ref-dip'] == df['test-dip'])
    perc = concordant / df.shape[0] * 100

    results = Results(
        ref_file=ref_file,
        ref_samples_total=rf.shape[0],
        test_file=test_file,
        test_samples_total=tf.shape[0],
        overlap_samples_total=df.shape[0],
        concordant_samples_total=concordant,
        concordant_percentage=perc,
    )

    if output_file:
        with open(output_file, 'w') as f:
            f.write(f"Reference file: {ref_file}\n")
            f.write(f"Reference samples total: {rf.shape[0]}\n")
            f.write(f"Test file: {test_file}\n")
            f.write(f"Test samples total: {tf.shape[0]}\n")
            f.write(f"Overlap samples total: {df.shape[0]}\n")
            f.write(f"Concordant samples total: {concordant}\n")
            f.write(f"Concordant percentage: {perc:.2f}\n")

    return results
