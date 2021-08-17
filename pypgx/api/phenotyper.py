from pypgx import get_activity_score

def phenotyper(target_gene, first_hap, second_hap):
    """Convert haplotype calls to a predicted phenotype.

    Different genes have different phenotypes. Many use a unit of enzyme
    activity known as an activity score (AS). For example, in the CYP2D6
    gene, the reference CYP2D6*1 allele has AS of 1 while a non-functional
    allele such as CYP2D6*4 has AS of 0. Note that in pypgx, star alleles
    with unknown function will have AS of -100; therefore, if one of the two
    haplotype calls contains an unknown function allele, the combined score
    will still be a negative number (e.g. AS=1 plus AS=-100 is still AS=-99)
    which translates to an unknown phenotype (e.g. unknown_metabolizer).

    Parameters
    ----------
    target_gene : str
        Name of the target gene.
    first_hap : str
        First haplotype call.
    second_hap : str
        Second haplotype call.

    Returns
    -------
    str
        Predicted phenotype.

    Examples
    --------

    .. code:: python

        from pypgx.phenotyper import phenotyper
        print(phenotyper("cyp2d6", "*1", "*1"))
        print(phenotyper("cyp2d6", "*1", "*4"))
        print(phenotyper("cyp2d6", "*1", "*2x2"))  # *2x2 is gene duplication.
        print(phenotyper("cyp2d6", "*4", "*5"))    # *5 is gene deletion.

    .. parsed-literal::

        normal_metabolizer
        intermediate_metabolizer
        ultrarapid_metabolizer
        poor_metabolizer

    Notes
    -----

    This method currently uses 6 different phenotyping algorithms:

    1. Default

    .. list-table::
        :widths: 50 50
        :header-rows: 1

        * - Phenotype
          - Summary
        * - unknown_function
          - AS < 0
        * - no_function
          - AS == 0
        * - decreased_function
          - 0 < AS < 2
        * - normal_function
          - AS == 2
        * - increased_function
          - AS > 2
        * - undetermined_function
          - All other cases

    2. Drug metabolizers (default)

    .. list-table::
        :widths: 50 50
        :header-rows: 1

        * - Phenotype
          - Summary
        * - unknown_metabolizer
          - AS < 0
        * - poor_metabolizer
          - AS == 0
        * - intermediate_metabolizer
          - 0 < AS <= 1.25
        * - normal_metabolizer
          - 1.25 < AS <= 2
        * - rapid_metabolizer
          - 2 < AS < 2.5
        * - ultrarapid_metabolizer
          - AS >= 2.5
        * - undetermined_metabolizer
          - All other cases

    3. Drug metabolizers - CYP2C9

    .. list-table::
        :widths: 50 50
        :header-rows: 1

        * - Phenotype
          - Summary
        * - unknown_metabolizer
          - AS < 0
        * - poor_metabolizer
          - AS <= 0.5
        * - intermediate_metabolizer
          - 0.5 < AS < 2
        * - normal_metabolizer
          - AS == 2
        * - ultrarapid_metabolizer
          - AS > 2
        * - undetermined_metabolizer
          - All other cases

    4. Drug metabolizers - CYP2D6

    .. list-table::
        :widths: 50 50
        :header-rows: 1

        * - Phenotype
          - Summary
        * - unknown_metabolizer
          - AS < 0
        * - poor_metabolizer
          - AS == 0
        * - intermediate_metabolizer
          - 0 < AS <= 1
        * - normal_metabolizer
          - 1 < AS <= 2.25
        * - ultrarapid_metabolizer
          - AS > 2.25
        * - undetermined_metabolizer
          - All other cases

    5. Drug metabolizers - DPYD

    .. list-table::
        :widths: 50 50
        :header-rows: 1

        * - Phenotype
          - Summary
        * - unknown_metabolizer
          - AS < 0
        * - poor_metabolizer
          - 0 <= AS <= 0.5
        * - intermediate_metabolizer
          - 0.5 < AS < 2
        * - normal_metabolizer
          - AS == 2
        * - ultrarapid_metabolizer
          - AS > 2
        * - undetermined_metabolizer
          - All other cases

    6. Drug transporters (default)

    .. list-table::
        :widths: 50 50
        :header-rows: 1

        * - Phenotype
          - Summary
        * - unknown_function
          - AS < 0
        * - poor_function
          - 0 <= AS <= 1
        * - decreased_function
          - 1 < AS <= 1.5
        * - normal_function
          - 1.5 < AS <= 2
        * - increased_function
          - AS > 2
        * - undetermined_function
          - All other cases

    """
    first_score = get_activity_score(target_gene, first_hap)
    second_score = get_activity_score(target_gene, second_hap)
    total_score = first_score + second_score
    if target_gene in ptcallers:
        return ptcallers[target_gene](total_score)
    else:
        return "no_phenotype"

def _phenotype_default(total):
    if total < 0:
        result = "unknown_function"
    elif total == 0:
        result = "no_function"
    elif 0 < total < 2:
        result = "decreased_function"
    elif total == 2:
        result = "normal_function"
    elif total > 2:
        result = "increased_function"
    else:
        result = "undetermined_function"
    return result

def _metabolizer_default(total):
    if total < 0:
        result = "unknown_metabolizer"
    elif total == 0:
        result = "poor_metabolizer"
    elif 0 < total <= 1.25:
        result = "intermediate_metabolizer"
    elif 1.25 < total <= 2:
        result = "normal_metabolizer"
    elif 2 < total < 2.5:
        result = "rapid_metabolizer"
    elif total >= 2.5:
        result = "ultrarapid_metabolizer"
    else:
        result = "undetermined_metabolizer"
    return result

def _metabolizer_cyp2d6(total):
    if total < 0:
        result = "unknown_metabolizer"
    elif total == 0:
        result = "poor_metabolizer"
    elif 0 < total <= 1:
        result = "intermediate_metabolizer"
    elif 1 < total <= 2.25:
        result = "normal_metabolizer"
    elif total > 2.25:
        result = "ultrarapid_metabolizer"
    else:
        result = "undetermined_metabolizer"
    return result

def _metabolizer_cyp2c9(total):
    if total < 0:
        result = "unknown_metabolizer"
    elif total <= 0.5:
        result = "poor_metabolizer"
    elif 0.5 < total < 2:
        result = "intermediate_metabolizer"
    elif total == 2:
        result = "normal_metabolizer"
    elif total > 2:
        result = "ultrarapid_metabolizer"
    else:
        result = "undetermined_metabolizer"
    return result

def _metabolizer_dpyd(total):
    if total < 0:
        result = "unknown_metabolizer"
    elif 0 <= total <= 0.5:
        result = "poor_metabolizer"
    elif 0.5 < total < 2:
        result = "intermediate_metabolizer"
    elif total == 2:
        result = "normal_metabolizer"
    elif total > 2:
        result = "ultrarapid_metabolizer"
    else:
        result = "undetermined_metabolizer"
    return result

def _transporter_default(total):
    if total < 0:
        result = "unknown_function"
    elif 0 <= total <= 1:
        result = "poor_function"
    elif 1 < total <= 1.5:
        result = "decreased_function"
    elif 1.5 < total <= 2:
        result = "normal_function"
    elif 2 < total:
        result = "increased_function"
    else:
        result = "undetermined_function"
    return result

ptcallers = {
    "abcb1": _phenotype_default,
    "cacna1s": _phenotype_default,
    "cftr": _phenotype_default,
    "cyp1a1": _metabolizer_default,
    "cyp1a2": _metabolizer_default,
    "cyp1b1": _metabolizer_default,
    "cyp2a6": _metabolizer_default,
    "cyp2a13": _metabolizer_default,
    "cyp2b6": _metabolizer_default,
    "cyp2c8": _metabolizer_default,
    "cyp2c9": _metabolizer_cyp2c9,
    "cyp2c19": _metabolizer_default,
    "cyp2d6": _metabolizer_cyp2d6,
    "cyp2e1": _metabolizer_default,
    "cyp2f1": _metabolizer_default,
    "cyp2j2": _metabolizer_default,
    "cyp2r1": _metabolizer_default,
    "cyp2s1": _metabolizer_default,
    "cyp2w1": _metabolizer_default,
    "cyp3a4": _metabolizer_default,
    "cyp3a5": _metabolizer_default,
    "cyp3a7": _metabolizer_default,
    "cyp3a43": _metabolizer_default,
    "cyp4a11": _metabolizer_default,
    "cyp4a22": _metabolizer_default,
    "cyp4b1": _metabolizer_default,
    "cyp4f2": _metabolizer_default,
    "cyp17a1": _metabolizer_default,
    "cyp19a1": _metabolizer_default,
    "cyp26a1": _metabolizer_default,
    "dpyd": _metabolizer_dpyd,
    "g6pd": _phenotype_default,
    "gstm1": _phenotype_default,
    "gstp1": _phenotype_default,
    "gstt1": _phenotype_default,
    "ifnl3": _phenotype_default,
    "nat1": _phenotype_default,
    "nat2": _phenotype_default,
    "nudt15": _metabolizer_default,
    "por": _phenotype_default,
    "ptgis": _phenotype_default,
    "ryr1": _phenotype_default,
    "slc15a2": _metabolizer_default,
    "slc22a2": _metabolizer_default,
    "slco1b1": _transporter_default,
    "slco1b3": _transporter_default,
    "slco2b1": _transporter_default,
    "sult1a1": _phenotype_default,
    "tbxas1": _phenotype_default,
    "tpmt": _metabolizer_default,
    "ugt1a1": _metabolizer_default,
    "ugt1a4": _phenotype_default,
    "ugt2b7": _phenotype_default,
    "ugt2b15": _phenotype_default,
    "ugt2b17": _phenotype_default,
    "vkorc1": _phenotype_default,
    "xpc": _phenotype_default,
}
