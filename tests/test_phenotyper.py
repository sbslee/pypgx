from pypgx.phenotyper import phenotyper

def test_phenotyper():
    assert phenotyper("cyp2d6", "*1", "*1") == "normal_metabolizer"
    assert phenotyper("cyp2d6", "*1", "*4") == "intermediate_metabolizer"
    assert phenotyper("cyp2d6", "*1", "*2x2") == "ultrarapid_metabolizer"
    assert phenotyper("cyp2d6", "*4", "*5") == "poor_metabolizer"
