import chemicals
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple



element_symbols = [i.symbol for i in chemicals.periodic_table]
_charge_patterns = None
def get_charge_patterns():
    """Get all possible charge patterns for elements"""
    global _charge_patterns
    if _charge_patterns is not None:
        return _charge_patterns
    # Basic charge signs and their repeating patterns
    signs = ['+', '-']
    numbers = list(range(1, 6))  # 1 through 5
    
    patterns = []
    
    # Add repeating signs (e.g., +, ++, +++, etc.)
    for sign in signs:
        for i in range(1, 6):
            patterns.append(sign * i)
    
    # Add number-sign combinations (e.g., 1+, 2-, +1, -2, etc.)
    for num in numbers:
        for sign in signs:
            # Number before sign (e.g., 1+)
            patterns.append(f"{num}{sign}")
            # Sign before number (e.g., +1)
            patterns.append(f"{sign}{num}")
    patterns = frozenset(patterns)
    _charge_patterns = patterns
    return patterns

# Function to check if text matches element charge pattern
def has_element_charge(t):
    # Remove any whitespace
    t = t.strip()
    charge_patterns = get_charge_patterns()
    
    # Check for element symbols with various charge notations
    for symbol in element_symbols:
        if not t.startswith(symbol):
            continue
            
        rest = t[len(symbol):]
        if not rest:  # Just the element symbol
            continue
            
        # Check bare patterns (Fe+ etc)
        if rest in charge_patterns:
            return True
            
        # Check parenthesized patterns (Fe(+) etc)
        if rest.startswith('(') and rest.endswith(')'):
            if rest[1:-1] in charge_patterns:
                return True
                
    return False

def is_roman_numeral(s: str, max_val: int = 3999, rule_of_four=False) -> bool:
    """
    Validates if a string is a valid Roman numeral and within the specified maximum value.
    
    Args:
        s: String to validate
        max_val: Maximum decimal value allowed (default 3999)
    
    Returns:
        bool: True if valid Roman numeral within range, False otherwise
    """
    if not isinstance(s, str) or not s:
        return False

    # Valid Roman numeral symbols and their values
    roman_values = {
        'I': 1,
        'V': 5,
        'X': 10,
        'L': 50,
        'C': 100,
        'D': 500,
        'M': 1000
    }

    # Check for invalid characters
    if not all(c in roman_values for c in s):
        return False

    # Check for invalid repetitions
    for char in ['V', 'L', 'D']:
        if s.count(char) > 1:
            return False
    
    if rule_of_four:
        for char in ['I', 'X', 'C', 'M']:
            if char * 4 in s:
                return False

    # Calculate value while checking subtractive rules
    total = 0
    i = 0
    while i < len(s):
        # Look ahead for subtractive pairs
        if i + 1 < len(s) and roman_values[s[i]] < roman_values[s[i + 1]]:
            # Validate subtractive pairs
            curr, next_char = s[i], s[i + 1]
            
            # Only I, X, and C can be used as subtractives
            if curr not in ['I', 'X', 'C']:
                return False
                
            # Check valid subtractive pairs
            valid_pairs = {
                'I': ['V', 'X'],
                'X': ['L', 'C'],
                'C': ['D', 'M']
            }
            
            if curr not in valid_pairs or next_char not in valid_pairs[curr]:
                return False
                
            # Calculate subtractive value
            total += roman_values[next_char] - roman_values[curr]
            i += 2
        else:
            total += roman_values[s[i]]
            i += 1

    return total <= max_val

@dataclass
class Span:
    """Represents a section of text with its position and possible nested spans"""
    start: int
    end: int
    text: str
    nested: List['Span'] = field(default_factory=list)
    
    def __post_init__(self):
        assert len(self.text) == self.end - self.start

@dataclass
class CaseEdit:
    """Represents a suggested edit to the case of a span of text"""
    span: Span
    case: str  # 'preserve', 'lower', 'upper'
    confidence: float  # 0.0 to 1.0
    
    def apply(self, text: str) -> str:
        """Apply this edit to the text"""
        if self.case == 'preserve':
            return text
        elif self.case == 'lower':
            return text.lower()
        elif self.case == 'upper':
            return text.upper()
        raise ValueError(f"Unknown case: {self.case}")

class CaseChecker:
    """Base class for implementing case checking rules"""
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        raise NotImplementedError

def find_spans(text: str) -> List[Span]:
    """
    Recursively split text into spans, handling nested brackets.
    Returns list of Span objects with nested structure.
    """
    def find_spans_recursive(text: str, offset: int = 0) -> List[Span]:
        spans = []
        current_start = 0
        bracket_stack = []
        bracket_pairs = {
            ')': '(',
            ']': '[',
            '}': '{',
            '>': '<'
        }
        
        i = 0
        while i < len(text):
            char = text[i]
            
            if char in '([{<':
                if not bracket_stack:  # Start of new bracket section
                    if i > current_start:
                        # Add text before bracket as a span
                        spans.append(Span(
                            offset + current_start,
                            offset + i,
                            text[current_start:i]
                        ))
                    current_start = i
                bracket_stack.append((char, i))
                
            elif char in ')]}>' and bracket_stack:
                matching_open = bracket_pairs[char]
                if bracket_stack[-1][0] == matching_open:
                    start_char, start_idx = bracket_stack.pop()
                    if not bracket_stack:  # End of outermost bracket group
                        # Recursively process the content inside brackets
                        inner_text = text[start_idx + 1:i]
                        inner_spans = find_spans_recursive(
                            inner_text,
                            offset + start_idx + 1
                        )
                        
                        # Create span for the entire bracket expression
                        bracket_span = Span(
                            offset + current_start,
                            offset + i + 1,
                            text[current_start:i + 1],
                            nested=inner_spans
                        )
                        spans.append(bracket_span)
                        spans.extend(inner_spans)
                        current_start = i + 1
            
            i += 1
        
        # Add final span if any
        if current_start < len(text):
            spans.append(Span(
                offset + current_start,
                offset + len(text),
                text[current_start:]
            ))
        
        return spans
    
    return find_spans_recursive(text)
class FormulaChecker(CaseChecker):
    """Checker for chemical formulas like CH3Br, H2O, etc."""
    def __init__(self):
        self.element_symbols = set(element_symbols)
        # Sort longer symbols first to avoid partial matches (e.g., 'Br' before 'B')
        self.sorted_symbols = sorted(element_symbols, key=len, reverse=True)
        self.strip_chars = '-[]'
        self.charge_patterns = get_charge_patterns()
        
    def is_chemical_formula(self, text: str) -> bool:
        """
        Check if text matches pattern of elements and numbers.
        Example valid patterns: CH3Br, H2O, C6H8O7
        """
        if not text:
            return False
            
        # Strip certain characters from ends before checking
        text = text.strip(self.strip_chars)
        if text.startswith('(') and text.endswith(')') and self.is_chemical_formula(text[1:-1]):
            return True

        pos = 0
        found_any = False
        
        while pos < len(text):
            # Try to match an element symbol at current position
            matched = False
            for symbol in self.sorted_symbols:
                if text[pos:].startswith(symbol):
                    found_any = True
                    matched = True
                    pos += len(symbol)
                    # Skip any following digits
                    while pos < len(text) and text[pos].isdigit():
                        pos += 1
                    break
            
            if not matched:
                # If we've found elements and reached the end, check if remainder is a charge
                if found_any and pos < len(text):
                    remaining = text[pos:]
                    return remaining in self.charge_patterns
                return False                
        return found_any

    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        # Check if the whole span is a chemical formula
        if self.is_chemical_formula(span.text):
            return [CaseEdit(span, 'preserve', 0.9)]
        return []

def apply_case_edits(text: str, edits: List[CaseEdit]) -> str:
    """
    Apply all case edits to the text, handling overlapping edits by confidence
    """
    # Sort edits by confidence (highest first)
    sorted_edits = sorted(edits, key=lambda e: e.confidence, reverse=True)
    
    # Convert string to list for easier manipulation
    result = list(text)
    # Track confidence of changes for each character (-1 means unmodified)
    confidences = [-1.0] * len(text)
    
    # Apply edits in order of confidence
    for edit in sorted_edits:
        span = edit.span
        
        # Check if any character in this span has been modified with higher confidence
        span_confidences = confidences[span.start:span.end]
        if any(conf >= edit.confidence for conf in span_confidences if conf >= 0):
            continue
            
        # Apply the edit and update confidences
        section = ''.join(result[span.start:span.end])
        new_section = edit.apply(section)
        result[span.start:span.end] = new_section
        confidences[span.start:span.end] = [edit.confidence] * (span.end - span.start)
    
    return ''.join(result)
# Example checkers
MAX_CHARGE = 15
class RomanNumeralChecker(CaseChecker):
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        if is_roman_numeral(span.text, MAX_CHARGE) or (span.text.startswith('(') and span.text.endswith(')') and is_roman_numeral(span.text[1:-1], MAX_CHARGE)):
            return [CaseEdit(span, 'preserve', 0.9)]
        return []

class ChemicalChargeChecker(CaseChecker):
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        if '(' in span.text and ')' in span.text:
            # Check if it's a charge notation like (II) or (2+)
            inner = span.text.strip('()')
            if (is_roman_numeral(inner, MAX_CHARGE) or 
                any(c in inner for c in '+-') or 
                (inner.isdigit() and any(c in span.text for c in '+-'))):
                return [CaseEdit(span, 'preserve', 0.9)]
        return []

class ChemicalNameChecker(CaseChecker):
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        # If it looks like a regular word (no numbers or special chars)
        if span.text.isalpha():
            return [CaseEdit(span, 'lower', 0.5)]
        return []

class DefaultLowercaseChecker(CaseChecker):
    """Default checker that suggests lowercase for any text not handled by other checkers"""
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        # Always suggest lowercase with very low confidence
        # This ensures other more specific checkers can override it
        return [CaseEdit(span, 'lower', 0.1)]

class ElementChargeChecker(CaseChecker):
    """Checker for element charges like Fe(+++) or Cu(2+)"""
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        if has_element_charge(span.text):
            # High confidence to override other checkers
            return [CaseEdit(span, 'preserve', 0.94)]
        return []

# Common scientific units and measurements that should preserve case
SCIENTIFIC_UNITS = {
    'pH',
    '°C', '°F',              # temperature
    'mg/mL', 'µg/mL',        # concentrations
    'mg/L', 'µg/L',          # concentrations
    'mL', 'µL', 'L',         # volumes
    'kg', 'mg', 'µg', 'g',   # masses 
    'mM', 'µM', 'nM',        # molar concentrations
    'kDa', 'Da',             # molecular weights
    'microg/mL',             # alternative notation
}

# Pre-compile the regex patterns for each term
UNIT_PATTERNS = {
    term: re.compile(
        r'(?:^|[\s([{])(' + re.escape(term) + r')(?=$|[\s)\]}]|[.,;:])'
    ) for term in SCIENTIFIC_UNITS
}

class ScientificUnitChecker(CaseChecker):
    """Checker for scientific units and measurements like °C, mg/mL etc"""
    
    def __init__(self):
        # Common prefixes that should preserve case
        self.special_terms = {
            'pH',
            '°C', '°F',  # temperature
            'mg/mL', 'µg/mL','mg/L', 'µg/L',  # concentrations
            'mL', 'µL', 'L',   # volumes
            'kg', 'mg', 'µg', 'g',  # masses
            'mM', 'µM', 'nM',  # molar concentrations
            'kDa', 'Da',       # molecular weights
            'microg/mL',
        }
        
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        text = span.text.strip()
        
        # Look for units at the end of the text or after a space
        for pattern in UNIT_PATTERNS.values():
            # Match if unit:
            # - is at start of string (^)
            # - follows whitespace (\s)
            # - follows punctuation ([(\[{])
            # And ensure it's followed by:
            # - end of string
            # - whitespace
            # - punctuation
            for match in re.finditer(pattern, text):
                unit = match.group(1)
                unit_start = match.start(1) + span.start
                unit_end = match.end(1) + span.start
                
                return [CaseEdit(
                    Span(unit_start, unit_end, unit),
                    'preserve',
                    0.85
                )]
                    
        return []

CI_RE = re.compile(r'^\s*\(\s*\d+\s*CI(?:\s*,\s*\d+\s*CI)*\s*\)\s*$')
class CINomenclatureChecker(CaseChecker):
    """Checker for Chemical Index nomenclature like (8CI,9CI) or (8CI, 9CI)"""
    def __init__(self):
        # Pattern explanation:
        # ^\s*        - Start of string, optional whitespace
        # \(          - Opening parenthesis
        # \s*         - Optional whitespace after opening parenthesis
        # \d+\s*CI    - One or more digits, optional whitespace, CI
        # (?:         - Start non-capturing group for additional CI terms
        #   \s*,\s*   - Comma with optional whitespace around it
        #   \d+\s*CI  - One or more digits, optional whitespace, CI
        # )*          - Zero or more additional CI terms
        # \s*         - Optional whitespace before closing parenthesis
        # \)          - Closing parenthesis
        # \s*$        - Optional whitespace, end of string
        self.ci_pattern = CI_RE
        
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        # Check if this span matches our CI pattern
        if self.ci_pattern.match(span.text.strip()):
            return [CaseEdit(span, 'preserve', 0.8)]
        return []

class NotationChecker(CaseChecker):
    """Checker for chemical notation patterns like N-, O-, D-, L-, etc."""
    def __init__(self):
        self.notation_patterns = [
            '(RS)-', '(±)-', '(R)-', '(S)-',
            'α-', 'β-', 'γ-', 'δ-', 'D-', 'L-', 'N-', 'O-', 'S-'
        ]
        # Pre-build spans for each pattern to know their lengths
        self.pattern_lengths = {pat: len(pat) for pat in self.notation_patterns}
        
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        # Check if the span starts with any of our patterns
        for pattern in self.notation_patterns:
            if span.text.startswith(pattern):
                # Create two edits:
                # 1. Preserve case for the notation part
                # 2. Lowercase for the rest with slightly lower confidence
                pattern_len = self.pattern_lengths[pattern]
                edits = [
                    CaseEdit(
                        Span(span.start, span.start + pattern_len, span.text[:pattern_len]),
                        'preserve',
                        0.9
                    )
                ]
                
                # Only add second edit if there's text after the pattern
                if pattern_len < len(span.text):
                    edits.append(
                        CaseEdit(
                            Span(span.start + pattern_len, span.end, span.text[pattern_len:]),
                            'lower',
                            0.8
                        )
                    )
                return edits
                
            # Also check for the single character versions without hyphen
            # This catches cases where the hyphen might be in a different span
            if len(pattern) == 2 and pattern[1] == '-':  # Like 'N-', 'O-', etc
                if span.text == pattern[0]:  # Just the single letter
                    return [CaseEdit(span, 'preserve', 0.9)]
                    
        return []

refrig_suffixes = {'10', '11', '110', '111', '1112a', '1113', '1114', '112', '1120', '112a',
    '113', '1130(E)', '1132a', '113a', '114', '1140', '1141', '114B2', '114a', '115', '1150', '116', 
    '12', '120', '121', '1216', '1218', '121a', '122', '1224yd(Z)', '122a', '122b', '123', '1233zd(E)', 
    '1234yf', '1234ze(E)', '123a', '123b', '124', '124a', '125', '1270', '12B1', '12B2', '12a', '13', 
    '130', '130a', '131', '131a', '131b', '132', '132a', '132b', '132bB2', '132c', '133', '1336mzz(E)', 
    '1336mzz(Z)', '133a', '133b', '134', '134a', '13B1', '13I1', '14', '140', '140a', '141', '141B2', 
    '141a', '141b', '142a', '142b', '143', '143a', '14A', '150', '150a', '151', '151a', '152', '152a',
    '160', '161', '170', '20', '21', '211', '212', '213', '214', '215', '216', '216ca', '217', '217ba', 
    '218', '22', '221', '222', '222c', '223', '223ca', '223cb', '224', '224ca', '224cb', '224cc', '225', 
    '225aa', '225ba', '225bb', '225ca', '225cb', '225cc', '225da', '225ea', '225eb', '226', '226ba', 
    '226ca', '226cb', '226da', '226ea', '227ca', '227ca2', '227ea', '227me', '22B1', '22a', '23', '231', 
    '232', '232ca', '232cb', '233', '233ca', '233cb', '233cc', '234', '234aa', '234ab', '234ba', '234bb', 
    '234bc', '234ca', '234cb', '234cc', '234cd', '234da', '234fa', '234fb', '235', '235ca', '235cb', '235cc', 
    '235da', '235fa', '236cb', '236ea', '236fa', '236me', '241', '242', '243', '243ca', '243cb', '243cc', 
    '243da', '243ea', '243ec', '244', '244ba', '244bb', '244ca', '244cb', '244cc', '244da', '244db', '244ea', 
    '244eb', '244ec', '244fa', '244fb', '245ca', '245cb', '245ea', '245eb', '245fa', '245mc', '245mf', '245qc',
    '251', '252', '252ca', '252cb', '252dc', '252ec', '253', '253ba', '253bb', '253ca', '253cb', '253ea', '253eb', 
    '253ec', '253fa', '253fb', '253fc', '254cb', '254pc', '261', '261ba', '262', '262ca', '262fa', '262fb', '263', 
    '271', '271b', '271d', '271fb', '272', '280', '281', '290', '30', '31', '32', '329ccb', '338eea', '347ccd', 
    '347mcc', '347mmy', '34E', '365mfc', '40', '400', '401A', '401B', '401C', '402A', '402B', '403A', '403B',
    '404A', '405A', '406A', '406B', '407A', '407B', '407C', '407D', '407E', '407F', '408A', '409A', '409B',
    '41', '410A', '410B', '411A', '411B', '411C', '412A', '413A', '414A', '414B', '415A', '415B', '416A',
    '417A', '417B', '417C', '418A', '419A', '419B', '420A', '421A', '421B', '422A', '422B', '422B', '422C',
    '422D', '422E', '423A', '424A', '425A', '426A', '427A', '428A', '429A', '430A', '431A', '432A', '433A', 
    '433B', '433C', '434A', '435A', '436A', '436B', '437A', '438A', '439A', '440A', '441A', '442A', 
    '443A', '444A', '445A', '448A', '449A', '452A', '452B', '453A', '454A', '454B', '454C', '455A', '456A',
    '457A', '458A', '459A', '466A', '50', '500', '501', '502', '502a', '503', '504', '505', '506', '507',
    '507A', '508', '508A', '508B', '509', '509A', '510', '510A', '511', '511A', '512A', '513A', '514A', '515B',
    '600', '600a', '601', '601a', '610', '611', '630', '631', '702', '704', '717', '718', '720', '728', 
    '729', '732', '740', '744', '744A', '744R', '764', '784', 'C316', 'C317', 'C318', 'E125', 'E134', 'E143a', 'E170'}
for r in list(refrig_suffixes):
    refrig_suffixes.add(r.lower())
class RegistryIDChecker(CaseChecker):
    """Checker for database and registry identifiers with specific format requirements"""
    
    def __init__(self):
        # Each pattern is (prefix, [separators], validation_function)
        self.registry_patterns = [
            # E numbers must be all digits
            ('E', [''], lambda rest: rest.isdigit()),
            ('Q', [''], lambda rest: rest.isdigit()),
            ('QZY', [''], lambda rest: rest.isdigit()),
            ('AKOS', [''], lambda rest: rest.isdigit()),
            ('D', [''], lambda rest: rest.isdigit()),
            ('CS', ['-'], lambda rest: rest.isdigit()),

            ('R', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('SP', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('DR', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('CFC', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('HFC', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('HCFC', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('HFO', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            ('HC', ['', '-', ' '], lambda rest: rest in refrig_suffixes),
            
            # Most registry IDs should be alphanumeric with optional dashes
            ('PDSP', ['-', ' '], lambda rest: all(c.isalnum() or c == '-'  or c == '_' for c in rest)),
            ('HY', ['-'], lambda rest: all(c.isalnum() or c == '-'  for c in rest)),
            ('UN', ['-', ' '], lambda rest: all(c.isalnum() or c == '-' for c in rest)),
            ('UNII', ['-'], lambda rest: all(c.isalnum() or c == '-' for c in rest)),
            ('DTXSID', ['', '-'], lambda rest: all(c.isalnum() or c == '-' for c in rest)),
            ('DTXCID', ['', '-'], lambda rest: all(c.isalnum() or c == '-' for c in rest)),
            ('CHEBI', [':', ' '], lambda rest: rest.replace('-', '').isalnum()),
            ('NSC', ['-', ' '], lambda rest: rest.replace('-', '').isalnum()),
            ('CAS', ['-', ' '], lambda rest: all(c.isdigit() or c == '-' for c in rest)),
            ('MFCD', ['-', ' '], lambda rest: all(c.isdigit() for c in rest)),
            ('EC', [' ', '-'], lambda rest: all(c.isalnum() or c == '-' for c in rest)),
            ('EINECS', [' ', '-'], lambda rest: all(c.isalnum() or c == '-' for c in rest)),
            ('BRN', [' ', '-'], lambda rest: rest.replace('-', '').isalnum()),
            ('HSDB', [' ', '-'], lambda rest: rest.isdigit()),
            ('CHEMBL', ['', ' ', '-'], lambda rest: rest.replace('-', '').isalnum()),
            ('MFCD', ['', ' ', '-'], lambda rest: rest.replace('-', '').isalnum()),
            ('SCHEMBL', ['', ' ', '-'], lambda rest: rest.replace('-', '').isalnum()),
            ('KBio', ['', ' ', '-'], lambda rest: rest.replace('-', '').isalnum()),
            ('SPBio', ['', ' ', '-'], lambda rest: rest.replace('-', '').isalnum()),
            ('BSPBio', ['', ' ', '-'], lambda rest: rest.replace('-', '').isalnum()),

            ('CCRIS', ['', ' ', '-'], lambda rest: rest.replace('-', '').isdigit()),
            ('STL', ['', ' ', '-'], lambda rest: rest.replace('-', '').isdigit()),
            ('NCGC', ['', ' ', '-'], lambda rest: rest.replace('-', '').isdigit()),
            ('AS', ['', ' ', '-'], lambda rest: rest.replace('-', '').isdigit()),
            ('DB', ['', ' ', '-'], lambda rest: rest.replace('-', '').isdigit()),
            ('MPCM', ['', ' ', '-'], lambda rest: rest.replace('-', '').isdigit()),
        ]
    
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        text = span.text.strip()
        
        # Check each registry pattern
        for prefix, separators, validator in self.registry_patterns:
            # If text starts with this prefix
            if text.upper().startswith(prefix):
                prefix_end = len(prefix)
                
                # Try each possible separator
                for separator in separators:
                    # Check if the separator matches
                    if separator:
                        if not text[prefix_end:].startswith(separator):
                            continue
                        id_start = prefix_end + len(separator)
                    else:
                        id_start = prefix_end
                    
                    # Get the rest of the text after prefix and separator
                    rest = text[id_start:]
                    
                    # Only preserve case if the rest matches the validation function
                    if rest and validator(rest):
                        return [CaseEdit(span, 'preserve', 0.95)]
        
        return []
class NomenclatureTagChecker(CaseChecker):
    """Checker for standard nomenclature tags like [USP], [INN], [USAN]"""
    
    def __init__(self):
        # Standard nomenclature tags that should preserve their case
        self.tags = {
            '[USP]', '[INN]', '[USAN]', '[BAN]', '[JAN]', '[WHO]',
            '[MI]', '[MART.]', '[VANDF]', '[HSDB]', '[EP]', '(MART.)',
            '[WHO-DD]'
        }
    
    def check(self, span: Span, full_text: str) -> List[CaseEdit]:
        # Clean the span text for comparison
        cleaned_text = span.text.strip()
        
        # Check if this span exactly matches one of our tags
        if cleaned_text in self.tags:
            return [CaseEdit(span, 'preserve', 0.9)]  # High confidence to preserve standard tags
            
        return []
def fix_chemical_case(text: str) -> str:
    """
    Main function to handle chemical name case fixing
    """
    # Initialize checkers
    checkers = [
        ElementChargeChecker(),
        RegistryIDChecker(),
        ScientificUnitChecker(),
        FormulaChecker(),
        NotationChecker(),
        NomenclatureTagChecker(),
        RomanNumeralChecker(),
        CINomenclatureChecker(),
        ChemicalChargeChecker(),
        ChemicalNameChecker(),
        DefaultLowercaseChecker(),
    ]
    
    # Find spans in text
    # Find spans in text and add full text span
    spans = find_spans(text)
    full_span = Span(0, len(text), text)
    spans = [full_span] + spans  # Add full span at the beginning

    # Collect all edits from all checkers
    all_edits = []
    for span in spans:
        for checker in checkers:
            edits = checker.check(span, text)
            all_edits.extend(edits)
    
    # Apply edits to get final result
    return apply_case_edits(text, all_edits)



# Example usage:
def fix_synonym_case(text):
    """
    Process a chemical synonym, converting to lowercase if appropriate.
    
    Args:
        text (str): Chemical synonym to process
        
    Returns:
        str: Processed chemical synonym
    """
    return fix_chemical_case(text)

