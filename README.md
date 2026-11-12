# Docking Interaction Cleaner

A robust Python utility to parse docking reference spreadsheets, match corresponding per-project interaction results, clean amino acid residue strings, and generate structured summary and detailed Excel reports.

## Features

- Reads docking reference files containing peptide sequences, docking scores, and selected models.
- Recursively locates local interaction results across project folders.
- Cleans and formats amino acid residue identifiers using regular expressions.
- Exports a dual-sheet Excel report with ranked summaries and detailed interactions.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

Run in the current directory:
```bash
python clean_docking_interactions.py
```

Or specify custom root directories, docking spreadsheets, and output filenames:
```bash
python clean_docking_interactions.py /path/to/screening_root -d "PEPTIDE MUTATION docking.xlsx" -o Custom_Cleaned.xlsx
```
