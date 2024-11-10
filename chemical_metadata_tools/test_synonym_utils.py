from chemical_metadata_tools.synonym_utils import *
import pytest



def test_basic_valid_numerals():
    assert is_roman_numeral('I') is True
    assert is_roman_numeral('V') is True
    assert is_roman_numeral('X') is True
    assert is_roman_numeral('L') is True
    assert is_roman_numeral('C') is True
    assert is_roman_numeral('D') is True
    assert is_roman_numeral('M') is True

def test_complex_valid_numerals():
    assert is_roman_numeral('MCMLIV') is True  # 1954
    assert is_roman_numeral('MCMXC') is True   # 1990
    assert is_roman_numeral('MMXIV') is True   # 2014
    assert is_roman_numeral('MMMCMXCIX') is True  # 3999

def test_invalid_characters():
    assert is_roman_numeral('ABC') is False
    assert is_roman_numeral('XII!') is False
    assert is_roman_numeral('X1V') is False
    assert is_roman_numeral('') is False
    assert is_roman_numeral('  ') is False

def test_invalid_repetitions():
    assert is_roman_numeral('IIII', rule_of_four=True) is False
    assert is_roman_numeral('VV') is False
    assert is_roman_numeral('LL') is False
    assert is_roman_numeral('DD') is False

def test_invalid_subtractive_patterns():
    assert is_roman_numeral('IC') is False    # Invalid jump
    assert is_roman_numeral('XM') is False    # Invalid jump
    assert is_roman_numeral('VX') is False    # V cannot subtract

def test_max_value_constraint():
    assert is_roman_numeral('MMMM') is False  # 4000 > default max of 3999
    assert is_roman_numeral('MMMM', max_val=4000) is True  # Equal to max
    assert is_roman_numeral('MMMM', max_val=3999) is False  # Exceeds max

def test_edge_cases():
    assert is_roman_numeral(None) is False
    assert is_roman_numeral(123) is False
    assert is_roman_numeral('  IV  ') is False
    assert is_roman_numeral('iv') is False    # Must be uppercase

def test_subtractive_rules():
    assert is_roman_numeral('IV') is True     # 4
    assert is_roman_numeral('IX') is True     # 9
    assert is_roman_numeral('XL') is True     # 40
    assert is_roman_numeral('XC') is True     # 90
    assert is_roman_numeral('CD') is True     # 400
    assert is_roman_numeral('CM') is True     # 900

def test_numbers_1_to_15():
    # Single digit numbers
    assert is_roman_numeral('I') is True      # 1
    assert is_roman_numeral('II') is True     # 2
    assert is_roman_numeral('III') is True    # 3
    assert is_roman_numeral('IV') is True     # 4
    assert is_roman_numeral('V') is True      # 5
    assert is_roman_numeral('VI') is True     # 6
    assert is_roman_numeral('VII') is True    # 7
    assert is_roman_numeral('VIII') is True   # 8
    assert is_roman_numeral('IX') is True     # 9
    assert is_roman_numeral('X') is True      # 10
    assert is_roman_numeral('XI') is True     # 11
    assert is_roman_numeral('XII') is True    # 12
    assert is_roman_numeral('XIII') is True   # 13
    assert is_roman_numeral('XIV') is True    # 14
    assert is_roman_numeral('XV') is True     # 15

    # Invalid alternatives that some might try
    assert is_roman_numeral('IIII') is True  # 4 should be IV
    assert is_roman_numeral('VIIII') is True # 9 should be IX
    assert is_roman_numeral('XIIII') is True # 14 should be XIV






def get_test_cases():
    """
    Returns list of (input, expected_output) tuples for testing.
    expected_output is the string in its desired case format.
    """
    return [
        # Registry numbers/identifiers (preserve case)
        ("CAS-12345", "CAS-12345"),
        ("UNII-5QTD157FPB", "UNII-5QTD157FPB"),
        ("DTXSID12345", "DTXSID12345"),
        ("CHEBI:12345", "CHEBI:12345"),
        ("NSC-12345", "NSC-12345"),
        
        # Chemical formulas (preserve case)
        ("CH3Br", "CH3Br"),
        ("C2H5Br", "C2H5Br"),
        ("H2O", "H2O"),
        ("CH2BrCl", "CH2BrCl"),
        ("C6H8O7", "C6H8O7"),
        
        # Stereochemistry notations (preserve case)
        ("(R)-limonene", "(R)-limonene"),
        ("(S)-alanine", "(S)-alanine"),
        ("(RS)-ibuprofen", "(RS)-ibuprofen"),
        ("cis-platinum", "cis-platinum"),
        ("trans-stilbene", "trans-stilbene"),
        
        # Greek letters (preserve case)
        ("α-tocopherol", "α-tocopherol"),
        ("β-carotene", "β-carotene"),
        ("γ-cyclodextrin", "γ-cyclodextrin"),
        
        # Standard chemical names (lowercase)
        ("Ethanol", "ethanol"),
        ("METHANOL", "methanol"),
        ("Acetic Acid", "acetic acid"),
        ("Sodium Chloride", "sodium chloride"),
        
        # Tags and suffixes (preserve case)
        ("aspirin [USP]", "aspirin [USP]"),
        ("paracetamol [INN]", "paracetamol [INN]"),
        ("acetaminophen [USAN]", "acetaminophen [USAN]"),
        
        # Real examples from database (preserve case)
        ("CAS-630-08-0", "CAS-630-08-0"),
        ("NSC-12345", "NSC-12345"),
        ("pH 7.4", "pH 7.4"),
        ("25 °C", "25 °C"),
        ("E236", "E236"),
        ("UNII-7U1EE4V452", "UNII-7U1EE4V452"),
        
        # Real examples from database (should lowercase)
        ("Carbon Monoxide", "carbon monoxide"),
        ("ETHANOIC ACID", "ethanoic acid"),
        ("Formic Acid", "formic acid"),
        ("METHANOIC ACID", "methanoic acid"),
        ("Hydrochloric Acid", "hydrochloric acid"),
        ("PALMITIC ACID", "palmitic acid"),
        
        # Mixed cases requiring preservation
        ("N-methylacetamide", "N-methylacetamide"),
        ("O-ethyloxime", "O-ethyloxime"),
        ("S-methylcysteine", "S-methylcysteine"),
        ("D-glucose", "D-glucose"),
        ("L-alanine", "L-alanine"),
        
        # Temperature/concentration units (preserve case)
        ("37 °C", "37 °C"),
        ("98.6 °F", "98.6 °F"),
        ("5 mg/mL", "5 mg/mL"),
        ("10 µg/mL", "10 µg/mL"),
        
        # Real chemical synonyms from database
        # ("CH3-[CH2]14-COOH", "CH3-[CH2]14-COOH"),
        ("Carbon(II) Oxide", "carbon(II) oxide"),
        ("METHYLCHLOROFORM", "methylchloroform"),
        
        # Edge cases
        ("P-Xylene", "p-xylene"),
        ("M-Cresol", "m-cresol"),
        
        ("O-Benzenediol", "O-benzenediol"),
        ("", ""),
        (" ", " "),
        ("123", "123"),

        # Basic charges
        ("Fe+", "Fe+"),
        ("Fe-", "Fe-"),
        ("Na+", "Na+"),
        ("Cl-", "Cl-"),
        
        # Multiple charges with numbers
        ("Fe2+", "Fe2+"),
        ("Fe3+", "Fe3+"),
        ("Cu2+", "Cu2+"),
        ("Po3-", "Po3-"),
        
        # Multiple + or - signs
        ("Fe++", "Fe++"),
        ("Fe+++", "Fe+++"),
        ("Po--", "Po--"),
        ("Po---", "Po---"),
        ("Fe++++", "Fe++++"),
        ("Po----", "Po----"),
        ("Fe+++++", "Fe+++++"),
        ("Po-----", "Po-----"),
        
        # Charges in parentheses
        ("Fe(+)", "Fe(+)"),
        ("Fe(-)", "Fe(-)"),
        ("Fe(2+)", "Fe(2+)"),
        ("Fe(3+)", "Fe(3+)"),
        ("Po(1-)", "Po(1-)"),
        ("Po(2-)", "Po(2-)"),
        ("Fe(+++)", "Fe(+++)"),
        ("Po(---)", "Po(---)"),

        ("CrH4O4-1", "CrH4O4-1"),
        ("H6O6Sn-3", "H6O6Sn-3"),
        ("Hexahydroxystannate(IV)", "hexahydroxystannate(IV)"),
        ("Hexahydroxystannate(I)", "hexahydroxystannate(I)"),
        ("Hexahydroxystannate(II)", "hexahydroxystannate(II)"),

        ("UN 1828", "UN 1828"),

        ("Radium chloride (RaCl2)", "radium chloride (RaCl2)"),

        ("Lithium silicate (Li2SiO3) (6CI)", "lithium silicate (Li2SiO3) (6CI)"),

        # Basic units with numbers
        ("100 mL", "100 mL"),
        ("50 µL", "50 µL"),
        ("25 L", "25 L"),
        ("10 mg", "10 mg"),
        ("5 µg", "5 µg"),
        ("1 kg", "1 kg"),
        ("100 g", "100 g"),
        
        # Concentrations
        ("10 mg/mL", "10 mg/mL"),
        ("5 µg/mL", "5 µg/mL"),
        ("2.5 mg/L", "2.5 mg/L"),
        ("0.5 µg/L", "0.5 µg/L"),
        ("100 microg/mL", "100 microg/mL"),
        
        # Molar concentrations
        ("10 mM", "10 mM"),
        ("5 µM", "5 µM"),
        ("2.5 nM", "2.5 nM"),
        
        # Molecular weights
        ("50 kDa", "50 kDa"),
        ("10 Da", "10 Da"),
        
        # Temperature
        ("37 °C", "37 °C"),
        ("98.6 °F", "98.6 °F"),
        ("-20 °C", "-20 °C"),
        
        # pH values
        ("pH 7", "pH 7"),
        ("pH 7.4", "pH 7.4"),
        ("pH 5.5", "pH 5.5"),
        
        # Decimal numbers
        ("0.5 mL", "0.5 mL"),
        ("1.5 µg", "1.5 µg"),
        ("2.75 mg/mL", "2.75 mg/mL"),
        
        # Negative numbers
        ("-10 °C", "-10 °C"),
        ("-5 °F", "-5 °F"),
        
        # Scientific notation
        ("1e-3 mM", "1e-3 mM"),
        ("1.5e3 µL", "1.5e3 µL"),
        
        # Mixed case input
        # ("100 Ml", "100 mL"),
        ("50 Mg", "50 mg"),
        ("25 KG", "25 kg"),
        
        # Words containing unit letters (should not be modified)
        ("Light", "light"),
        ("Literature", "literature"),
        ("Long", "long"),
        ("Milligrams", "milligrams"),
        
        # Units in sentences
        ("dissolved in 100 mL", "dissolved in 100 mL"),
        ("heated to 37 °C", "heated to 37 °C"),
        ("adjusted to pH 7.4", "adjusted to pH 7.4"),
        ("concentration of 5 mg/mL", "concentration of 5 mg/mL"),
        
        # Multiple units in one string
        ("100 mL of 5 mg/mL", "100 mL of 5 mg/mL"),
        ("50 µL at 37 °C", "50 µL at 37 °C"),
        ("pH 7.4 and 25 °C", "pH 7.4 and 25 °C"),
        
        # Edge cases
        # ("100mL", "100mL"),  # No space (might want to discuss handling)
        ("mL", "mL"),  # Just the unit
        # ("L solution", "l solution"),  # L at start without number
        # ("in L", "in L"),  # L without number
        # ("THE PH IS 7", "the pH is 7"),  # pH in caps
        # ("at PH 7", "at pH 7"),  # pH in caps with value
        
        # Tricky cases with 'L'
        ("Large", "large"),
        ("Light blue", "light blue"),
        ("Low temperature", "low temperature"),
        ("5 L container", "5 L container"),
        ("in 5 L", "in 5 L"),
        
        # Range values
        ("5-10 mL", "5-10 mL"),
        ("10 - 20 mg", "10 - 20 mg"),
        ("between 5 and 10 mL", "between 5 and 10 mL"),
        
        # With parentheses
        ("(100 mL)", "(100 mL)"),
        ("(pH 7.4)", "(pH 7.4)"),
        ("(at 37 °C)", "(at 37 °C)"),
        
        # Complex combinations
        ("100 mL of pH 7.4 buffer at 37 °C", "100 mL of pH 7.4 buffer at 37 °C"),
        ("5 mg/mL in 10 mM at pH 7", "5 mg/mL in 10 mM at pH 7"),
        ("stored at -20 °C in 100 µL aliquots", "stored at -20 °C in 100 µL aliquots"),        
    ]

@pytest.mark.parametrize("input_text,expected_result", get_test_cases())
def test_fix_synonym_case(input_text: str, expected_result: str):
    """Test if the function correctly adjusts the case of chemical synonyms"""
    result = fix_synonym_case(input_text)
    assert result == expected_result, \
        f"Failed for '{input_text}': got '{result}', expected '{expected_result}'"


def test_to_debug():
    # fix_synonym_case('PALMITIC ACID')
    # fix_synonym_case("Fe(+++)")

    # fix_synonym_case("E236")
    # fix_synonym_case("CH3-[CH2]14-COOH")
    # fix_synonym_case("Hexahydroxystannate(IV)")
    # fix_synonym_case("Radium chloride (RaCl2)")
    # fix_synonym_case("Lithium silicate (Li2SiO3) (6CI)")


    # fix_chemical_case("(pH 7.4)")
    # fix_chemical_case("2.5 mg/L")
    fix_chemical_case("zinc hydroxide (Zn(OH)2)")


