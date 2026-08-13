# AGENTS.md

Guidance for agents (and humans) working in this repository. This project builds
and maintains the chemical identifier databases (TSV files) that ship inside the
[`chemicals`](https://github.com/CalebBell/chemicals) Python package
(`chemicals.identifiers`). It is a *data curation* repo more than a software
project: most of the value is in the TSV/JSON data files, and the Python code
exists to (re)generate and validate them from source PDFs, `.mol` files, and
PubChem lookups.

This document explains what each piece does and how the pieces fit together.
It does not change, and should not be used to justify changing, the actual
pipeline — see "Ground rules" below.

## Ground rules for agents

- **Do not change the pipeline/workflow** (Justfile targets, script logic,
  column order, file naming conventions) unless explicitly asked to. Fix data,
  not process.
- **The checked-in `*.tsv` files are the source of truth** that ships to users
  of the `chemicals` package. They are *regenerated* from `mol/*.mol` +
  `Parsed CAS metadata.json` + `Good synoynms by CAS.json` (+ live PubChem
  lookups), but regeneration requires `rdkit`, `thermo`, `chemicals`,
  `pubchempy`, network access to PubChem, and (for a full `parse-pdf-*` run)
  SciFinder/CAS Common Chemistry PDFs that are not stored in git. In most
  agent sessions none of that is available, so the practical way to fix a
  single bad entry is to hand-edit the row directly in the relevant `*.tsv`
  file (and, if the error originates from a hardcoded override, also fix it
  in that directory's `generate_chemical_metadata.py` `custom_compounds`/
  `good_syns` dict so the fix survives the next real regeneration).
- **Never break the TSV column contract** (see "Row format" below) when
  editing: same tab-separated field order, empty synonym columns are fine to
  drop (not pad), and don't introduce new columns.
- CAS Registry Numbers are the primary key within each database. They must be
  unique per file and have a valid checksum (see below).
- These are curated *identifier* databases (name/formula/structure), not
  property databases — there's no boiling point, density, etc. here.

## Repository layout

```
chemical-metadata/
├── Justfile                     # all pipeline commands (see below)
├── generate_db_from_cids.py     # ad hoc: print a DB row for given PubChem CID(s)
├── opsin_2.7.0_name_to_smiles_mapping_from_chemicals_metadata.tsv
│                                 # OPSIN-derived name<->SMILES map; feeds
│                                 # iupac_names.py (currently short-circuited off)
├── chemical_metadata_tools/     # shared Python package used by every subfolder
│   ├── parse_CAS_data.py        # PDF -> HTML -> "Parsed CAS metadata.json"
│   ├── generate_db_tools.py     # mol file + JSON + PubChem -> one TSV row
│   ├── synonym_utils.py         # name/case cleanup helpers used when building synonyms
│   ├── iupac_names.py           # OPSIN name-standardization lookup (disabled, see below)
│   ├── compare_db.py            # diff two TSV generations, human-readable report
│   ├── duplicate_searcher.py    # cross-database duplicate detector (installed DBs)
│   └── remove_redundant_from_pubchem_db.py
│                                 # strip CAS numbers already covered by the
│                                 # specialized DBs out of the big PubChem DBs
├── anions/                      # anion database (charge < 0)
├── cations/                     # cation database (charge > 0)
├── inorganic/                   # neutral inorganic compounds
├── organic/                     # neutral organic compounds ("example user db")
├── duplicates/                  # output of duplicate_searcher.py (JSON reports)
└── old/                         # legacy/retired code (ChemSep import, elements,
                                  # a very old standalone `chemical_metadata` package).
                                  # Not part of the current pipeline.
```

Each of `anions/`, `cations/`, `inorganic/`, `organic/` is a self-contained
"database project" with the same shape:

```
<db>/
├── <Name> db.tsv                     # the shipped database (source of truth)
├── <db>_preferences.json             # preferred/unpreferred CAS for duplicate formula/structure groups
├── Parsed CAS metadata.json          # output of parse_CAS_data.py (from pdf/html)
├── Good synoynms by CAS.json         # output of generate_chemical_metadata.py (hand overrides)
├── generate_chemical_metadata.py     # per-db hardcoded overrides -> "Good synoynms by CAS.json"
├── mol/                              # one <CAS>.mol file per compound (from PubChem/SciFinder)
├── pdf/, html/                       # (inorganic/anions/cations) SciFinder PDF sources + converted HTML
└── old/, old_needs_review/           # superseded data kept for reference
```

Note the *db* names are inconsistent on purpose (they match what's installed
into `chemicals`):
- `anions/Anion db.tsv`
- `cations/Cation db.tsv`
- `inorganic/Inorganic db.tsv`
- `organic/chemical identifiers example user db.tsv` (this is `chemicals`'
  `PUBCHEM_EXAMPLE_DB_NAME`)

`organic/db_preprocessor.py` and `inorganic/db_preprocessor.py` are actually
each directory's `generate_chemical_metadata.py`-equivalent input-prep step —
despite the name, `db_preprocessor.py` builds `Good synoynms by CAS.json`
from a big dict of hand-curated per-CAS overrides (`custom_compounds`), while
`generate_chemical_metadata.py` in the *same* directory is the pipeline
script that actually reads `mol/*.mol` and writes the TSV. This split is
confusing — `anions/` and `cations/` fold that override step into
`generate_chemical_metadata.py` itself (no separate `db_preprocessor.py`),
`inorganic/` and `organic/` split it out. Check `Justfile` if in doubt about
which script does what for a given directory.

## Row format (the TSV contract)

Every database TSV is tab-separated, **no header row**, one compound per
line, sorted with not-yet-preferred/duplicate compounds first and
PubChem-preferred compounds last (see `sort_key` in
`generate_db_tools.ChemicalMetadataProcessor.process_files`). Columns:

| # | Field | Notes |
|---|-------|-------|
| 1 | `cid` | PubChem CID, or `-1` if none/not applicable |
| 2 | `CAS` | CAS Registry Number, hyphenated (e.g. `108-88-3`) |
| 3 | `formula` | Hill-order formula from `serialize_formula`/RDKit; ions carry a trailing charge suffix, e.g. `BiO4-3`, `MnO+3` (this is *not* a real element — don't parse charge digits as atom counts) |
| 4 | `molecular_weight` | float, rounded to 6 decimals |
| 5 | `smiles` | canonical SMILES (may be blank for some hardcoded/isotopic entries) |
| 6 | `inchi` | InChI string with the `InChI=1S/`/`InChI=1/` prefix stripped |
| 7 | `inchikey` | standard InChIKey |
| 8 | `iupac_name` | PubChem/CAS-derived name, treated as "the" IUPAC-ish name |
| 9 | `common_name` | the preferred display/common name |
| 10+ | `synonym...` | variable-length list of additional names/synonyms (may be a single empty string if there are none) |

`chemical_metadata_tools/compare_db.py` and the `chemicals` package itself
(`ChemicalMetadataDB`/`load_chemical_file` in
`remove_redundant_from_pubchem_db.py`) both parse rows this way — keep any
edit consistent with `values[0:9]` + `values[7:]`-style parsing (columns 8
and 9, i.e. `iupac_name`/`common_name`, are duplicated into the front of the
synonym list by some consumers, which is expected).

## The generation pipeline, end to end

For a given category (`anions`, `cations`, `inorganic`, `organic`), the full
refresh flow (driven by `Justfile`, all runnable via `just <recipe>`) is:

1. **`parse-pdf-<db>`** — `cd <db> && python3 -m chemical_metadata_tools.parse_CAS_data`.
   Converts every `pdf/*.pdf` (SciFinder/CAS substance detail exports) to
   `html/*.html` (via `pdftohtml`), scrapes each HTML for name, formula, MW,
   alternate/deleted CAS numbers and "Other Names", and merges in any locally
   cached CAS Common Chemistry JSON (`~/.cache/chemical_metadata/common_chemistry/<CAS>`).
   Writes `<db>/Parsed CAS metadata.json`.
2. **Hardcoded overrides** — each directory's override script
   (`generate_chemical_metadata.py` for `anions`/`cations`, or
   `db_preprocessor.py` for `inorganic`/`organic`) builds a big
   `custom_compounds`/`good_syns` dict keyed by CAS number and writes
   `<db>/Good synoynms by CAS.json`. This is where you hardcode a specific
   PubChem CID, name, SMILES, InChI, formula, MW, extra synonyms, or a
   `preferred`/not-preferred flag for a CAS number that the automatic
   pipeline gets wrong (duplicate structures resolve to one PubChem CID,
   isotopes/spin-isomers that RDKit/InChI can't round-trip, polymers to
   exclude, etc). `hardcoded_synonyms: True` means *replace* the synonym list
   instead of merging into it.
3. **`generate-tsv-<db>`** — `cd <db> && python3 generate_chemical_metadata.py mol/*.mol <db>_db_YYYY_MM_DD.tsv`.
   For every `mol/<CAS>.mol` file: read structure (preferring PubChem SDF
   metadata embedded in the file when present), round-trip through InChI for
   canonical structure/charge handling, verify the resulting formal charge
   matches what's expected for that category (`require_charge=0` for neutral
   inorganic/organic, `require_not_charge=0` for cations/anions — see
   `ChemicalMetadataProcessor.__init__`), look up additional PubChem
   synonyms/IUPAC name by CID or InChIKey (cached under
   `~/.cache/chemical_metadata/<db>/`), then combine everything with the
   SciFinder data and hardcoded overrides
   (`_combine_with_scifinder`) to produce one TSV row per compound
   (`ChemicalMetadataProcessor._format_output`). Compounds whose name
   contains "polymer"/"poly" are dropped (`remove_unwanted_compounds_after_processing`,
   with a couple of CAS-specific exceptions). CAS numbers in each script's
   `ignore_CASs` set (structures RDKit/InChI can't handle correctly) are
   skipped entirely. Output is sorted (`process_files`) so
   `preferences.json`-preferred rows sort last.
4. **`compare-<db>`** — diffs the new dated TSV against the currently-shipped
   one (`chemical_metadata_tools/compare_db.py`) and writes/prints
   `<db>/changes.txt`: compounds added/removed, per-property changes, and
   synonym added/removed lists. **Always inspect this before copying** — it's
   the main defense against silently regressing an entry.
5. **`update-<db>`** — runs steps 3+4 together (generate then compare) and
   cats `changes.txt`.
6. **`copy-<db>`** — after you've reviewed `changes.txt` and are satisfied,
   moves the dated `<db>_db_YYYY_MM_DD.tsv` over the shipped
   `<db> db.tsv`/`chemical identifiers example user db.tsv`, and moves the
   matching `_preferences.json` into place (`anion_preferences.json` etc;
   `organic` has no preferences file). This is a destructive overwrite of the
   previous database — only run it once `changes.txt` looks right.
7. **`install-<db>`** — copies the shipped TSV (+ preferences JSON) into the
   *installed* `chemicals` package's data folder
   (`python3 -c "import chemicals; print(chemicals.identifiers.folder)"`),
   so a locally installed `chemicals` immediately picks up the new data
   without a release. This does **not** touch this git repo.

`just update-all` / `copy-all` / `install-all` run all four categories
together. `just clean-all-caches` wipes the PubChem/mol-processing disk
caches (`~/.cache/chemical_metadata/<db>/`) if you need to force fresh
PubChem lookups. `just clean-all-changes` removes the `changes.txt` reports.
`just parse-all-pdfs` / `generate-all-tsvs` / `compare-all` run just that one
step across all four categories.

None of this pipeline runs automatically in CI — it's a manual, human-
(or agent-)supervised refresh process. There is no step that regenerates
`mol/*.mol` files themselves; new compounds get a new `.mol` file added by
hand (typically exported from PubChem or SciFinder) before any of the above
will pick them up.

### Other top-level tools

- **`chemical_metadata_tools.duplicate_searcher`** (`just duplicate-check`) —
  loads the databases as *installed* in the `chemicals` package (not the
  repo copies) and cross-checks PubChem CID / CAS / SMILES / InChI /
  InChIKey / formula / common name for collisions across the cation, anion,
  inorganic, "example", "small", and "large" PubChem databases, writing
  `duplicates/*.json`. A collision that's already resolved by exactly one
  `preferred_cas` entry (from the three `*_preferences.json` files) and the
  rest `unpreferred_cas` is filtered out as "already handled". Note this
  operates on the databases as installed via `just install-all`, so you must
  install before running it for it to see your latest edits.
- **`chemical_metadata_tools.remove_redundant_from_pubchem_db`**
  (`just remove-redundant`) — removes CAS numbers already present in the
  cation/anion/inorganic/example databases from the large generic PubChem
  "small"/"large" databases installed in `chemicals`, to avoid duplicate
  entries across databases. Also operates on the installed package, not this repo.
- **`generate_db_from_cids.py`** — standalone CLI, `python3 generate_db_from_cids.py <cid> [<cid> ...]`,
  prints one draft TSV-ish row (CAS column left blank) per PubChem CID by
  querying `pubchempy` directly. Useful for hand-drafting a new entry before
  it has a `.mol` file, but its output format is a slightly different/older
  shape than the real pipeline's — treat it as a starting point, not
  something to paste straight into a database file.
- **`chemical_metadata_tools/iupac_names.py`** — intended to build a
  SMILES→"best" name lookup from
  `opsin_2.7.0_name_to_smiles_mapping_from_chemicals_metadata.tsv` (picks the
  shortest, lowest-punctuation, lowercase name per SMILES). Currently
  **disabled**: `if is_debugger_active() or 1:` always takes the `or 1`
  branch and sets `iupac_standard_names = {}`, so this lookup is a no-op in
  practice. Worth knowing if `iupac_name` values look like they ignore that
  TSV — they currently do, unconditionally.
- **`chemical_metadata_tools/synonym_utils.py`** — text-processing helpers
  used while building the synonym list: `fix_synonym_case` (title-case
  cleanup that's careful around element symbols/charges/Roman numerals/
  chemical-prefix locants so it doesn't mangle formulas like `Fe2O3` or
  `(R)-`), plus supporting helpers (`has_element_charge`, `is_roman_numeral`,
  bracket/span parsing). Has a real unit test suite —
  `chemical_metadata_tools/test_synonym_utils.py` (`pytest chemical_metadata_tools/`).
- **`chemical_metadata_tools/generate_db_tools.py`** — the shared engine
  described in step 3 above (`ChemicalMetadataProcessor`), plus standalone
  helpers: `is_inchikey`/`remove_inchikeys` (strip stray InChIKeys that leak
  into synonym lists), `remove_bad_synonyms` (drops synonyms matching junk
  patterns — vendor grade/purity specs, packaging, batch/lot numbers,
  database-ID-shaped strings, etc.), `deduplicate_names` (case-insensitive
  dedup preferring the more-lowercase variant).

## Setup / running the pipeline locally

This repo's scripts depend on packages that are **not** installed in a
default agent sandbox: `rdkit`, `thermo`, `chemicals`, `pubchempy`, `joblib`,
`diskcache`, `appdirs`, `sortedcontainers`, `pandas`, `numpy`. Also
`pdftohtml` (from `poppler-utils`) is required for `parse-pdf-*`. Most of
these can be `pip install`-ed; PubChem-dependent steps additionally need
outbound network access to `pubchem.ncbi.nlm.nih.gov`, which may be slow or
rate-limited — the pipeline caches PubChem responses via `diskcache` under
`~/.cache/chemical_metadata/` specifically because repeated full-database
regeneration is expensive.

Given that, most day-to-day agent work in this repo is **direct data
curation** — reading/fixing the checked-in TSVs and their JSON inputs — not
running the full regeneration pipeline. If you do have the dependencies
available and want to verify a fix reproduces correctly, the fastest sanity
check for a single compound is `python3 generate_db_from_cids.py <cid>`
(only needs `pubchempy`) or reformulating with `rdkit` directly (formula via
`rdkit.Chem.rdMolDescriptors.CalcMolFormula(mol, True, True)`, MW via
`rdkit.Chem.Descriptors.MolWt(mol)`) rather than running the full per-category
`generate-tsv-*` recipe.

## Data-quality conventions worth knowing before editing

- **CAS checksum**: last digit = (sum of each preceding digit × its distance
  from the check digit, reading right-to-left) mod 10. A CAS number that
  fails this is either a typo or (rarely, see `inorganic/generate_chemical_metadata.py`)
  a synthetic non-CAS placeholder like `2099990000-00-0` used for a
  compound with no real CAS number — check `custom_compounds` before
  "fixing" one of those.
- **`common_name` casing**: names are generally stored lowercase-first
  (`lower_case_first_letter_name` in `parse_CAS_data.py`) unless the second
  character is itself uppercase (e.g. acronym-like names), matching CAS
  index-name style (`sulfuric acid`, not `Sulfuric Acid`).
- **Ionic formulas** encode charge as a trailing `+`/`-` run or `+`/`-`
  followed by a digit (`Fe+3`, `SO4-2`) — this is `serialize_formula`
  notation, not a literal atom count.
- **Synonym noise is deliberately filtered**, not exhaustive: vendor grade/
  purity strings, packaging, batch numbers, bare InChIKeys, and other
  low-value synonyms are stripped by `remove_bad_synonyms`/`remove_inchikeys`
  during generation. If you're hand-adding a synonym to a TSV row, keep it to
  genuine alternate names (trade names, other-language names, IUPAC
  variants) rather than something the pipeline would have filtered out.
- **Preferred/duplicate CAS pairs**: several CAS numbers can point at the
  same PubChem CID/structure (different mineral forms, hydrates, or historic
  vs. current CAS registrations of "the same" substance). These are resolved
  via `preferred: True/False` in a `custom_compounds` entry, not by deleting
  either row — both stay in the database, and `*_preferences.json` /
  `duplicate_searcher.py` use the flag to recognize the duplicate as
  intentional. Don't unilaterally delete a "duplicate-looking" row without
  checking whether it's actually one of these designed pairs. **When both
  members of such a pair share a structure, they should have byte-identical
  `cid`/`formula`/`molecular_weight`/`smiles`/`inchi`/`inchikey`/`iupac_name`
  fields** (only CAS, `common_name`, and the synonym list differ) — if they
  don't, that's a strong signal one of them regressed (this is exactly how
  the CAS 1309-38-2/1345-25-1 bug below was found: a duplicate-InChIKey scan
  across `preferred`-pair CAS numbers turned up one pair that was actually a
  copy-paste of an *unrelated* compound, not the intended structure-sharing
  pair).
- **Mixed-metal-oxide `formula` column is sometimes a reduced empirical
  ratio, not the true stoichiometry.** For dozens of `inorganic/Inorganic
  db.tsv` rows (e.g. CoTiO3, Li2MnO3, CuCr2O4, CdWO4, Bi2MoO6, B8K2O13...),
  the `formula` column stores a 1:1:1-style reduced ratio (a side effect of
  how the underlying mol file/InChI represents these compounds) while
  `iupac_name`/synonyms correctly carry the true stoichiometric formula in
  parentheses, e.g. `CoOTi` / `cobalt titanium oxide (CoTiO3)`. This is
  deliberate and consistent — do **not** "fix" the name to match the
  formula column for these. The tell for a *genuine* stoichiometry bug
  (fixed in this repo for CAS 12137-12-1, 12065-65-5, 1315-03-3) is a row
  whose *own* synonym list is internally self-contradictory — e.g. it uses
  names for two different, unrelated compounds/CAS numbers (a millerite/NiS
  row that also said "nickel sulfide (Ni3S4)", a real but different
  compound) — not just a mismatch against the reduced formula column.
- **PubChem records for extended/lattice solids (oxides, minerals) are
  sometimes simply wrong**, not just imprecisely reduced. CAS 1317-61-9 /
  1309-38-2 (Fe3O4, iron oxide/magnetite) both resolved to PubChem CID
  9816051, whose own record ("iron;tetrahydrate", formula `Fe3H8O4`, MW
  239.6, SMILES `O.O.O.O.[Fe].[Fe].[Fe]`) represents Fe3O4 as three bare
  iron atoms plus four water molecules — not a real structure. Fixed by
  hand with a proper charge-balanced ionic SMILES
  (`[Fe+2].[Fe+3].[Fe+3].[O-2].[O-2].[O-2].[O-2]`, MW 231.531), matching the
  convention already used for Fe2O3 elsewhere in the same file. If you hit
  another oxide/mineral CAS whose MW/formula looks physically wrong, check
  the raw `PUBCHEM_MOLECULAR_FORMULA`/`PUBCHEM_IUPAC_NAME` tags inside its
  `mol/<CAS>.mol` file before trusting it — low-quality PubChem CIDs like
  this one are a real, recurring failure mode for compounds without a
  well-defined molecular structure, not something the pipeline can catch
  automatically.
- **A one-off text-mangling bug hit `organic/chemical identifiers example
  user db.tsv`'s synonym lists** (not the primary `iupac_name`/`common_name`
  fields, and not the other three databases beyond a couple of stray
  instances): a handful of punctuation characters were replaced by their
  literal (mis-identified) Unicode character names or lost entirely somewhere
  in a PDF/HTML→text conversion step, e.g. a prime/apostrophe (`′`/`'`)
  became the literal text `" inverted exclamation mark"`, a comma became
  `" pound not"`, a lost synonym-list separator became `"pound>>"`, and
  alpha (`α`) became `"I+/-"`. Fixed with a handful of global substitutions
  (verified first that none collided with the literal word "compound",
  which contains "pound" as a substring) plus manual cleanup of a few
  unreconstructable vendor-spec fragments. If a future `parse-pdf-organic`
  run reintroduces this, it's almost certainly the same root cause in
  `pdftohtml`/`parse_CAS_data.parse_scifinder`'s handling of non-ASCII
  punctuation — worth fixing at the source rather than re-patching output.
- **Verified clean as of this review** (2026-08-13): no duplicate CAS
  numbers within or across the four databases, no duplicate InChIKeys
  across databases, no CAS number given as a synonym in one row that is
  actually a different row's primary CAS, no invalid CAS checksums, no
  empty `iupac_name`+`common_name` pairs, and (aside from the Fe3O4 case
  above) no formula/MW disagreement with the row's own SMILES when
  cross-checked with RDKit.
