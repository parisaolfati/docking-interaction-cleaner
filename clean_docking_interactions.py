# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import glob
import re
import pandas as pd


def clean_amino_acid_string(val_str):
    """Clean string to extract only valid amino acid formats (e.g., ALA1, GLY3)."""
    if pd.isna(val_str) or not str(val_str).strip():
        return ""
    
    # Extract valid amino acid codes with numbers
    valid_aa = re.findall(r'[A-Za-z]{3}\d+', str(val_str))
    if not valid_aa:
        tokens = [t.strip() for t in str(val_str).replace(',', ' ').split()]
        valid_aa = [t for t in tokens if re.search(r'[A-Za-z]', t) and not t.isdigit()]
        
    return ", ".join(sorted(set(valid_aa)))


def process_docking_cleaning(root_dir: Path, docking_excel_name: str, output_excel_name: str):
    """Process docking reference, map local interaction results, clean residues, and export summary reports."""
    root_dir = root_dir.resolve()
    docking_excel = root_dir / docking_excel_name
    final_output = root_dir / output_excel_name

    if not docking_excel.exists():
        raise FileNotFoundError(f"Docking reference file not found: {docking_excel}")

    print(f"Reading docking reference file from: {docking_excel}")
    df_docking = pd.read_excel(docking_excel)
    df_docking.columns = [str(c).strip() for c in df_docking.columns]

    pep_col = df_docking.columns[0]
    score_col = df_docking.columns[1]
    model_col = df_docking.columns[2]

    combined_results = []
    missing_count = 0

    print("Processing and cleaning interaction data...")

    for idx, row in df_docking.iterrows():
        peptide = str(row[pep_col]).strip()
        dock_score = row[score_col]
        selected_model_raw = str(row[model_col]).strip().lower()
        
        model_num = ''.join(filter(str.isdigit, selected_model_raw))
        if not model_num:
            model_num = '1'
        target_model_label = f"MODEL_{model_num}"

        folder_pattern = f"*{peptide}*"
        matching_folders = list(root_dir.glob(folder_pattern))
        
        found_local = False
        
        for pep_folder in matching_folders:
            if pep_folder.is_dir():
                excel_files = list(pep_folder.rglob("Detailed_Interaction_Results.xlsx"))
                
                if excel_files:
                    local_excel_path = excel_files[0]
                    try:
                        df_local = pd.read_excel(local_excel_path, sheet_name="Protein_Peptide_Pairs")
                        
                        if 'Model' in df_local.columns:
                            df_local['Model_Upper'] = df_local['Model'].astype(str).str.upper().str.strip()
                            df_filtered = df_local[df_local['Model_Upper'] == target_model_label]
                        else:
                            df_filtered = df_local
                            
                        if not df_filtered.empty:
                            found_local = True
                            for _, r in df_filtered.iterrows():
                                rec_val = str(r.get("Target_Protein", r.get("Protein_Target", ""))).strip()
                                pep_val = str(r.get("Peptide_Residue", "")).strip()
                                
                                combined_results.append({
                                    "Peptide": peptide,
                                    "Docking_Score": dock_score,
                                    "Selected_Model": target_model_label.lower(),
                                    "Receptor_Residue_ChainA": rec_val,
                                    "Peptide_Residue_ChainB": pep_val,
                                    "Interaction_Type": r.get("Interaction_Type", r.get("Type", "")),
                                    "Distance_A": r.get("Distance_A", r.get("Distance", ""))
                                })
                    except Exception as e:
                        print(f"Error reading file {local_excel_path}: {e}")
                        
                if found_local:
                    break

        if not found_local:
            missing_count += 1
            combined_results.append({
                "Peptide": peptide,
                "Docking_Score": dock_score,
                "Selected_Model": f"model_{model_num}",
                "Receptor_Residue_ChainA": "No Interactions Found",
                "Peptide_Residue_ChainB": "No Interactions Found",
                "Interaction_Type": "-",
                "Distance_A": "-"
            })

    df_final = pd.DataFrame(combined_results)

    # Build summary sheet
    summary_list = []
    for (pep, score, mod), group in df_final.groupby(["Peptide", "Docking_Score", "Selected_Model"]):
        valid = group[group["Receptor_Residue_ChainA"] != "No Interactions Found"]
        
        raw_prot = " ".join(valid["Receptor_Residue_ChainA"].dropna().astype(str))
        raw_pep = " ".join(valid["Peptide_Residue_ChainB"].dropna().astype(str))
        
        clean_prot = clean_amino_acid_string(raw_prot)
        clean_pep = clean_amino_acid_string(raw_pep)
        
        summary_list.append({
            "Peptide": pep,
            "Docking_Score": score,
            "Selected_Model": mod,
            "Total_Interactions": len(valid),
            "Receptor_Residues_ChainA": clean_prot if clean_prot else "None",
            "Peptide_Residues_ChainB": clean_pep if clean_pep else "None"
        })

    df_summary = pd.DataFrame(summary_list).sort_values(by="Docking_Score", ascending=True)

    # Save output
    with pd.ExcelWriter(final_output, engine="openpyxl") as writer:
        df_summary.to_excel(writer, sheet_name="Peptide_Summary", index=False)
        df_final.to_excel(writer, sheet_name="Detailed_Interactions", index=False)

    print("
Processing completed successfully.")
    print(f"Cleaned output saved to: {final_output.resolve()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean and summarize peptide docking interactions based on reference spreadsheets.")
    parser.add_argument(
        "root_dir",
        nargs="?",
        default=".",
        type=Path,
        help="Root directory containing screening project folders (default: current directory)."
    )
    parser.add_argument(
        "-d", "--docking",
        default="PEPTIDE MUTATION docking.xlsx",
        help="Docking reference Excel filename (default: PEPTIDE MUTATION docking.xlsx)."
    )
    parser.add_argument(
        "-o", "--output",
        default="Selected_Models_Interactions_Cleaned.xlsx",
        help="Output cleaned Excel filename (default: Selected_Models_Interactions_Cleaned.xlsx)."
    )
    args = parser.parse_args()
    
    process_docking_cleaning(args.root_dir, args.docking, args.output)
