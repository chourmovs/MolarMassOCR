# main.py

import argparse
import os
import cv2
import numpy as np
from DECIMER import predict_SMILES
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw

def preprocess_image(path, scale=1.0, blur=3, method='adaptive', blocksize=21, C=4, suffix=""):
    img = cv2.imread(path)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    if scale != 1.0:
        img = cv2.resize(img, None, fx=scale, fy=scale, interpolation=cv2.INTER_CUBIC)
    if method == 'adaptive':
        proc = cv2.adaptiveThreshold(
            img, 255, cv2.ADAPTIVE_THRESH_MEAN_C,
            cv2.THRESH_BINARY, blocksize, C
        )
    elif method == 'otsu':
        _, proc = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    else:
        proc = img
    if blur > 0:
        proc = cv2.medianBlur(proc, blur)
    out = f"{os.path.splitext(path)[0]}_prep{suffix}.png"
    cv2.imwrite(out, proc)
    return out

def display_and_save_structure(smiles, out_name=None):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("  [ERREUR] Molécule non valide.")
        return
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mass = rdMolDescriptors.CalcExactMolWt(mol)
    print(f"  Formule brute: {formula}")
    print(f"  Masse molaire: {mass:.4f} Da")
    print(f"  SMILES: {smiles}")
    # Sauvegarde l'image de la molécule
    if out_name:
        img = Draw.MolToImage(mol, size=(300, 300))
        os.makedirs(os.path.dirname(out_name), exist_ok=True)
        img.save(out_name)
        print(f"  Aperçu structure sauvé dans {out_name}")

def process_image(img_path):
    params = [
        (1.0, 21, 4, 3, 'adaptive', "_A"),
        (1.5, 21, 4, 3, 'adaptive', "_B"),
        (2.0, 21, 4, 3, 'adaptive', "_C"),
        (1.0, 31, 6, 5, 'adaptive', "_D"),
        (1.5, 11, 2, 1, 'adaptive', "_E"),
        (1.0, 0, 0, 0, 'otsu', "_F"),
    ]
    results = []
    print(f"\n== Analyse de {img_path} ==")
    for scale, blocksize, C, blur, method, suffix in params:
        prep = preprocess_image(img_path, scale=scale, blocksize=blocksize, C=C, blur=blur, method=method, suffix=suffix)
        smiles = predict_SMILES(prep)
        mol = Chem.MolFromSmiles(smiles)
        formula = rdMolDescriptors.CalcMolFormula(mol) if mol else "??"
        mass = rdMolDescriptors.CalcExactMolWt(mol) if mol else None
        results.append({"smiles": smiles, "prep_img": prep, "formula": formula, "mass": mass})
    # Affiche seulement les variantes uniques (masse et SMILES)
    seen_masses = set()
    seen_smiles = set()
    eps = 0.01
    for idx, r in enumerate(results):
        if r["mass"] is None or r["smiles"] in seen_smiles:
            continue
        already_seen = any(abs(r["mass"] - m) < eps for m in seen_masses)
        if not already_seen:
            print(f"\n--- Variante {idx+1} ---")
            display_and_save_structure(r["smiles"], out_name=f"output/{os.path.basename(img_path)}_variant{idx+1}.png")
            seen_masses.add(r["mass"])
            seen_smiles.add(r["smiles"])
    # Suppression fichiers temporaires
    for r in results:
        try:
            os.remove(r["prep_img"])
        except OSError:
            pass

def main():
    parser = argparse.ArgumentParser(
        description="Reconnaissance optique de structures chimiques à partir d'images (multi-prétraitement)."
    )
    parser.add_argument("images", nargs="+", help="Un ou plusieurs fichiers image à analyser (.png .jpg ...)")
    args = parser.parse_args()
    for img_path in args.images:
        if not os.path.isfile(img_path):
            print(f"Fichier non trouvé: {img_path}")
            continue
        process_image(img_path)
    # Nettoyage global (optionnel, mais ne touche pas output/)
    import glob
    for ext in ("*.png", "*.jpg", "*.jpeg", "*.gif"):
        for file in glob.glob(ext):
            if not file.startswith("output/"):
                try: os.remove(file)
                except OSError: pass

if __name__ == "__main__":
    main()
