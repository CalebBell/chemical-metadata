# Justfile for chemical-metadata project
# Get the current date in YYYY_MM_DD format
date := datetime("%Y_%m_%d")
chemicals-dir := `python3 -c "import chemicals; print(chemicals.identifiers.folder)"`

default:
    @just --list --verbose

# Anions
# Remove cache for anions
clean-cache-anions:
    rm -rf ~/.cache/chemical_metadata/anions
# Parse PDFs for anions
parse-pdf-anions:
    cd anions && python3 -m chemical_metadata_tools.parse_CAS_data
# Generate new TSV for anions
generate-tsv-anions:
    cd anions && python3 generate_chemical_metadata.py mol/*.mol anion_db_{{date}}.tsv
# Compare new TSV with old TSV for anions
compare-anions:
    cd anions && python3 -m chemical_metadata_tools.compare_db "Anion db.tsv" anion_db_{{date}}.tsv
# Full update for anions (generate new db and compare)
update-anions:
    cd anions && \
    python3 generate_chemical_metadata.py mol/*.mol anion_db_{{date}}.tsv && \
    python3 -m chemical_metadata_tools.compare_db "Anion db.tsv" anion_db_{{date}}.tsv changes.txt
    cat anions/changes.txt
# Copy new anions database over old one and remove dated file
copy-anions:
    cd anions && [ -f anion_db_{{date}}.tsv ] && cp anion_db_{{date}}.tsv "Anion db.tsv" && rm anion_db_{{date}}.tsv || echo "New anion database not found"
# Cations
# Remove cache for cations
clean-cache-cations:
    rm -rf ~/.cache/chemical_metadata/cations
# Parse PDFs for cations
parse-pdf-cations:
    cd cations && python3 -m chemical_metadata_tools.parse_CAS_data
# Generate new TSV for cations
generate-tsv-cations:
    cd cations && python3 generate_chemical_metadata.py mol/*.mol cation_db_{{date}}.tsv
# Compare new TSV with old TSV for cations
compare-cations:
    cd cations && python3 -m chemical_metadata_tools.compare_db "Cation db.tsv" cation_db_{{date}}.tsv
# Full update for cations (generate new db and compare)
update-cations:
    cd cations && \
    python3 generate_chemical_metadata.py mol/*.mol cation_db_{{date}}.tsv && \
    python3 -m chemical_metadata_tools.compare_db "Cation db.tsv" cation_db_{{date}}.tsv changes.txt
    cat cations/changes.txt
# Copy new cations database over old one and remove dated file
copy-cations:
    cd cations && [ -f cation_db_{{date}}.tsv ] && cp cation_db_{{date}}.tsv "Cation db.tsv" && rm cation_db_{{date}}.tsv || echo "New cation database not found"
# Inorganic
# Remove cache for inorganic
clean-cache-inorganic:
    rm -rf ~/.cache/chemical_metadata/inorganic
# Parse PDFs for inorganic
parse-pdf-inorganic:
    cd inorganic && python3 -m chemical_metadata_tools.parse_CAS_data
# Generate new TSV for inorganic
generate-tsv-inorganic:
    cd inorganic && python3 db_preprocessor.py && python3 generate_chemical_metadata.py mol/*.mol inorganic_db_{{date}}.tsv
# Compare new TSV with old TSV for inorganic
compare-inorganic:
    cd inorganic && python3 -m chemical_metadata_tools.compare_db "Inorganic db.tsv" inorganic_db_{{date}}.tsv
# Full update for inorganic (generate new db and compare)
update-inorganic:
    cd inorganic && \
    python3 db_preprocessor.py && \
    python3 generate_chemical_metadata.py mol/*.mol inorganic_db_{{date}}.tsv && \
    python3 -m chemical_metadata_tools.compare_db "Inorganic db.tsv" inorganic_db_{{date}}.tsv changes.txt
    cat inorganic/changes.txt

# Copy new inorganic database over old one and remove dated file
copy-inorganic:
    cd inorganic && [ -f inorganic_db_{{date}}.tsv ] && cp inorganic_db_{{date}}.tsv "Inorganic db.tsv" && rm inorganic_db_{{date}}.tsv || echo "New inorganic database not found"
# Organic
# Remove cache for organic
clean-cache-organic:
    rm -rf ~/.cache/chemical_metadata/organic
# Parse PDFs for organic
parse-pdf-organic:
    cd organic && python3 -m chemical_metadata_tools.parse_CAS_data
# Generate new TSV for organic
generate-tsv-organic:
    cd organic && python3 db_preprocessor.py && python3 generate_chemical_metadata.py mol/*.mol organic_user_db_{{date}}.tsv
# Compare new TSV with old TSV for organic
compare-organic:
    cd organic && python3 -m chemical_metadata_tools.compare_db "chemical identifiers example user db.tsv" organic_user_db_{{date}}.tsv
# Full update for organic (generate new db and compare)
update-organic:
    cd organic && \
    python3 db_preprocessor.py && \
    python3 generate_chemical_metadata.py mol/*.mol organic_user_db_{{date}}.tsv && \
    python3 -m chemical_metadata_tools.compare_db "chemical identifiers example user db.tsv" organic_user_db_{{date}}.tsv changes.txt
    cat organic/changes.txt
# Copy new organic database over old one and remove dated file
copy-organic:
    cd organic && [ -f organic_user_db_{{date}}.tsv ] && cp organic_user_db_{{date}}.tsv "chemical identifiers example user db.tsv" && rm organic_user_db_{{date}}.tsv || echo "New organic database not found"

# Update all databases, main command - run this then inspect results and run `copy-all` (then optionally `install-all`)
update-all: update-anions update-cations update-inorganic update-organic
# Clean all caches
clean-all-caches: clean-cache-anions clean-cache-cations clean-cache-inorganic clean-cache-organic
# Parse all PDFs
parse-all-pdfs: parse-pdf-anions parse-pdf-cations parse-pdf-inorganic parse-pdf-organic
# Generate all TSVs
generate-all-tsvs: generate-tsv-anions generate-tsv-cations generate-tsv-inorganic generate-tsv-organic
# Compare all databases
compare-all: compare-anions compare-cations compare-inorganic compare-organic
# Remove all changes.txt files
clean-all-changes:
    rm -f anions/changes.txt cations/changes.txt inorganic/changes.txt organic/changes.txt
# Copy all new databases over old ones and remove dated files
copy-all: copy-anions copy-cations copy-inorganic copy-organic clean-all-changes


# Install individual databases to chemicals package
install-anions:
    cp anions/"Anion db.tsv" "{{chemicals-dir}}/Anion db.tsv"

install-cations:
    cp cations/"Cation db.tsv" "{{chemicals-dir}}/Cation db.tsv"

install-inorganic:
    cp inorganic/"Inorganic db.tsv" "{{chemicals-dir}}/Inorganic db.tsv"

install-organic:
    cp organic/"chemical identifiers example user db.tsv" "{{chemicals-dir}}/chemical identifiers example user db.tsv"

# Install all databases to chemicals package
install-all: install-anions install-cations install-inorganic install-organic

remove-redundant:
    python3 -m chemical_metadata_tools.remove_redundant_from_pubchem_db
