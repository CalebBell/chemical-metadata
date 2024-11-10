from codecs import open
import glob
import json
import os
from html.parser import HTMLParser
from joblib import Parallel, delayed
import subprocess
from pathlib import Path
import html
from chemicals import serialize_formula, nested_formula_parser, atoms_to_Hill, molecular_weight
import appdirs
import re
from sortedcontainers import SortedSet
cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)

common_chemistry_cache_dir = os.path.join(cache_dir, 'common_chemistry')
if not os.path.exists(common_chemistry_cache_dir):
    os.mkdir(common_chemistry_cache_dir)

def parse_chemical_formula(data):
    """Parse chemical formula from scifinder HTML data.
    
    Parameters
    ----------
    data : str
        HTML content from scifinder document
        
    Returns
    -------
    dict
        Dictionary containing the parsed formula data:
        - 'formula_raw': Raw formula string from document
        - 'formula_normalized': Normalized formula using Hill notation
        - 'atoms': Dictionary of atom counts
    Returns None if no valid formula is found.
    """
    try:
        # Get the section after CAS Registry Number
        section = data.split('CAS Registry Number')[1].split('<br/>')[0:3]
        
        # Look at each bold section until we find a valid formula
        for line in section:
            if '<b>' in line:
                potential_formula = line.split('<b>')[1].split('</b>')[0].strip()
                # Clean up the string
                potential_formula = potential_formula.replace('&#160;', '').strip()
                
                try:
                    # Try to parse it - if it succeeds, it's a formula
                    # atoms = nested_formula_parser(potential_formula)
                    normalized = serialize_formula(potential_formula)
                    return normalized
                except (ValueError, KeyError):
                    # If parsing fails, this bold section wasn't a formula
                    continue
                    
        # If we get here, no valid formula was found
        return None
            
    except (IndexError, KeyError):
        return None
def checkCAS(CASRN):
    try:
        check = CASRN[-1]
        CASRN = CASRN[::-1][1:]
        productsum = 0
        i = 1
        for num in CASRN:
            if num == '-':
                pass
            else:
                productsum += i*int(num)
                i += 1
        return (productsum % 10 == int(check))
    except:
        return False

def convert_pdf_to_html(pdf_path, html_path):
    """Convert a single PDF file to HTML using pdftohtml"""
    try:
        process = subprocess.Popen([
            'pdftohtml',
            '-i',
            '-enc', 'UTF-8',
            '-noframes',
            pdf_path,
            html_path
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Error converting {pdf_path}: {stderr.decode()}")
            return False
        return True
    except Exception as e:
        print(f"Exception while converting {pdf_path}: {str(e)}")
        return False

def process_pdfs_in_parallel(pdf_folder, html_folder, n_jobs=-1):
    """Process all PDFs in parallel using joblib"""
    pdf_files = list(Path(pdf_folder).glob('*.pdf'))
    html_files = list(Path(html_folder).glob('*.html'))
    
    if len(pdf_files) == len(html_files):
        print("All PDFs already converted to HTML")
        return
    
    # Create HTML folder if it doesn't exist
    Path(html_folder).mkdir(parents=True, exist_ok=True)
    
    # Prepare conversion tasks
    conversion_tasks = []
    for pdf_path in pdf_files:
        html_path = Path(html_folder) / f"{pdf_path.stem}.html"
        if not html_path.exists():
            conversion_tasks.append((str(pdf_path), str(html_path)))
    
    if not conversion_tasks:
        print("No new PDFs to convert")
        return
    
    print(f"Converting {len(conversion_tasks)} PDFs to HTML in parallel...")
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(convert_pdf_to_html)(pdf_path, html_path)
        for pdf_path, html_path in conversion_tasks
    )
    
    success_count = sum(1 for r in results if r)
    print(f"Successfully converted {success_count} of {len(conversion_tasks)} PDFs")

def parse_scifinder(data, CAS):
    d = {}
    possibilities = data.split('CAS Registry Number')[1].split('<b>')[1].split('</b>')[1].split('<br/>')
    if len(possibilities) > 2 and 'Manual Registration' in possibilities[2]:
        name = possibilities[1]
    else:
        name = possibilities[0]
    
    d['Name'] = lower_case_first_letter_name(name)

    formula = parse_chemical_formula(data)
    if formula:
        # d['Formula'] = formula # formula doesn't work with charge right now
        d['MW'] = molecular_weight(nested_formula_parser(formula))
    else:
        print(f'On CAS {CAS} and with name {name}, formula was missing or could not be parsed')
        
    if 'Alternate CAS Registry Numbers' in data:
        alt = data.split('Alternate CAS Registry Numbers')[1].split('<b>')[0].split('</b>')[1]
        alt = alt.split('<')[0]  # remove the end tag
        alt = alt.split('; ')    # If there are multiple hits - or split on ',' and include br
        d['Alternate CAS'] = alt
    if 'Deleted CAS Registry Numbers' in data:
        deleted = data.split('Deleted CAS Registry Numbers')[1].split('\n')[0].split('</b>')[1].replace('<br/>', '').replace(' ', '').split(',')
        d['Deleted CAS'] = deleted
    if 'Other Names' in data:
        hits = data.split('Other Names')[1]
        if 'Deleted CAS Registry Numbers' in hits:
            hits = hits.split('Deleted CAS Registry Numbers')[0]
        if 'Alternate CAS Registry Numbers' in hits:
            hits = hits.split('Alternate CAS Registry Numbers')[0]
        else:
            hits = hits.split('References')[0]
        hits = hits.split('\n')[0].split('</b>')[1].replace('<br/>', '')
        
        for bad_text in ['Alternate CAS Registry Numbers', 'Incompletely Defined Substance', 
                        'Coordination Compound', 'Manual Registration', 'Tabular Inorganic Substance']:
            hits = hits.split(bad_text)[0]
        hits = hits.replace('<b>', '')
        names = [i.strip() for i in hits.split(';')]
        d['Other Names'] = names
    return d

def remove_html_tags(text):
    # Replace both <sub> and <sup> tags with empty string but keep the content
    # Define patterns to clean
    patterns = [
        # Remove span tags with class
        (r'<span class="[^"]*">', ''),
        (r'</span>', ''),
        # Remove smallcaps/smallcap tags (both variants)
        (r'</?smallcaps?>', ''),
        # Remove sub/sup tags
        (r'</?sub>', ''),
        (r'</?sup>', ''),
        # Remove any other HTML tags that might appear
        (r'<[^>]+>', ''),
        # Clean up any double spaces that might result
        (r'\s+', ' '),
    ]

    cleaned = text
    for pattern, replacement in patterns:
        cleaned = re.sub(pattern, replacement, cleaned)

    # Trim any leading/trailing whitespace
    return cleaned.strip()
    
def clean_common_chemistry_synonyms(values):
    # Clean each value in the list
    return [remove_html_tags(v) for v in values]

def lower_case_first_letter_name(name):
    if len(name) > 1 and name[1].lower() == name[1]:
        return name[0].lower() + name[1:]
    return name




def parse_pdfs_in_subfolder(folder):
    h = HTMLParser()
    pdf_folder = folder / 'pdf'
    html_folder = folder / 'html'
    mol_folder = folder / 'mol'
    
    # Convert PDFs to HTML in parallel
    process_pdfs_in_parallel(pdf_folder, html_folder)
    
    # Process HTML files
    parsed_data = {}
    for html_path in html_folder.glob('*.html'):
        CAS = html_path.stem
        with open(html_path, encoding='utf-8') as f:
            text = html.unescape(f.read())
        dat = parse_scifinder(text, CAS)
        parsed_data[CAS] = dat

    # Also look through the common chemistry data online
    process_keys = list(parsed_data.keys())
    for mol_path in mol_folder.glob('*.mol'):
        CAS = mol_path.stem
        process_keys.append(CAS)
    process_keys = list(set(process_keys))

    for CAS in process_keys:
        potential_common_chemistry_path = os.path.join(common_chemistry_cache_dir, CAS)
        if os.path.exists(potential_common_chemistry_path):
            cc_data = json.loads(open(potential_common_chemistry_path).read())
            if cc_data:
                # May not have a PDF for a file
                if CAS not in parsed_data:
                    parsed_data[CAS] = {}
                if 'Name' in parsed_data[CAS]:
                    parsed_data[CAS]['Name'] = lower_case_first_letter_name(parsed_data[CAS]['Name'])
                
                # If we take a new name, make sure the old name is still a synonym
                if 'Other Names' in parsed_data[CAS]:
                    parsed_data[CAS]['Other Names'].append(parsed_data[CAS]['Name'])
                elif 'Name' in parsed_data[CAS]:
                    parsed_data[CAS]['Other Names'] = [parsed_data[CAS]['Name']]
                parsed_data[CAS]['Name'] = remove_html_tags(cc_data['name'])
                if 'synonyms' in cc_data:
                    parsed_data[CAS]['Other Names'] = list(SortedSet(clean_common_chemistry_synonyms(cc_data['synonyms']) + parsed_data[CAS].get('Other Names', [])))

                if 'inchi' in cc_data and cc_data['inchi']:
                    parsed_data[CAS]['inchi'] = cc_data['inchi']
                if 'inchiKey' in cc_data and cc_data['inchiKey']:
                    parsed_data[CAS]['inchiKey'] = cc_data['inchiKey']
                if 'canonicalSmile' in cc_data and cc_data['canonicalSmile']:
                    parsed_data[CAS]['smiles'] = cc_data['canonicalSmile']
                if 'molecularFormula' in cc_data and cc_data['molecularFormula']:
                    parsed_data[CAS]['formula'] = remove_html_tags(cc_data['molecularFormula'])
                    try:
                        parsed_data[CAS]['MW'] = molecular_weight(nested_formula_parser(parsed_data[CAS]['formula']))
                    except:
                        print(f"Failed to mass molecular formula {parsed_data[CAS]['formula']}")

        if CAS in parsed_data and 'Name' in parsed_data[CAS]:
            # force lower case names - can handle any exceptions otherwise. Make change if length > 1 and second character not a capital
            parsed_data[CAS]['Name'] = lower_case_first_letter_name(parsed_data[CAS]['Name'])

    # Write JSON output
    output_path = folder / 'Parsed CAS metadata.json'
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(parsed_data, f, indent=2, separators=(',', ': '), sort_keys=True)

if __name__ == '__main__':
    parse_pdfs_in_subfolder(Path(os.getcwd()))
