import os
import sys
import tkinter as tk
from tkinter import font as tkfont
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk, ImageGrab
import tempfile
import re
import numpy as np
import requests
import gradio_client as grc
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw, Descriptors

# --- Console interne ---
class ConsoleRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget
    def write(self, msg):
        self.text_widget.configure(state='normal')
        self.text_widget.insert('end', msg)
        self.text_widget.see('end')
        self.text_widget.configure(state='disabled')  # Fix: use 'disabled'
    def flush(self): pass


# --- Gradio client (à instancier UNE seule fois) ---
client = grc.Client("https://chouchouvs-distance-smiles.hf.space/")

def enhance_image(pil_img):
    from PIL import ImageEnhance, ImageFilter
    img = ImageEnhance.Contrast(pil_img).enhance(1.3)
    img = ImageEnhance.Sharpness(img).enhance(2.0)
    img = img.filter(ImageFilter.MedianFilter(size=3))
    return img

import tempfile
from gradio_client import Client, handle_file

def predict_smiles_with_gradio(pil_img):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        pil_img.save(tmp, format="PNG")
        tmp_path = tmp.name
    result = client.predict(
        img=handle_file(tmp_path),
        api_name="/predict"
    )
    return result




atomic_weights = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
    "S": 32.06, "P": 30.974, "F": 18.998, "Cl": 35.45,
    "Br": 79.904, "I": 126.90,
}
def masse_molaire_from_formula(formula):
    pattern = r"([A-Z][a-z]?)(\d*)"
    total = 0.0
    for (elt, count) in re.findall(pattern, formula):
        count = int(count) if count else 1
        w = atomic_weights.get(elt)
        if w is None: continue
        total += w * count
    return total

def is_chemically_plausible(mol):
    if not mol: return False
    for atom in mol.GetAtoms():
        s = atom.GetSymbol()
        v = atom.GetExplicitValence()
        if s == "C" and v > 4: return False
        if s == "O" and v > 2: return False
        if s == "N" and v > 3: return False
        if s == "H" and v > 1: return False
    return True

def smiles_has_isotopes(smiles):
    return bool(re.search(r'\[\d+[A-Z][a-z]?\]', smiles))

def smiles_to_cas(smiles):
    try:
        url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
        r = requests.get(url_cid, timeout=10)
        if not r.ok or not r.text.strip().isdigit():
            return None
        cid = r.text.strip()
        url_syn = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
        r2 = requests.get(url_syn, timeout=10)
        if not r2.ok: return None
        data = r2.json()
        ids = data.get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
        for regid in ids:
            if re.match(r'^\d{2,7}-\d{2}-\d$', regid):
                return regid
        return None
    except Exception: return None

def cas_to_name(cas):
    try:
        url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{cas}/cids/TXT"
        r = requests.get(url_cid, timeout=10)
        if not r.ok or not r.text.strip().isdigit():
            return None
        cid = r.text.strip()
        url_name = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,Title/JSON"
        r2 = requests.get(url_name, timeout=10)
        if not r2.ok: return None
        data = r2.json()
        props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
        return props.get("Title") or props.get("IUPACName")
    except Exception: return None

def clean_smiles_isotopes(smiles):
    return re.sub(r'\[(\d+)([A-Z][a-z]?)\]', r'\2', smiles)

import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import font as tkfont

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Reconnaissance de masse molaire OCSR")
        self.root.geometry("950x650")
        self.root.minsize(900, 600)
        self.bg_main = "#f5f8fc"
        self.color_accent = "#4F8EF7"
        self.bg_card = "#ffffff"
        self.bg_console = "#1a2233"
        self.fg_console = "#e0e6ef"
        self.font_title = tkfont.Font(family="Segoe UI", size=17, weight="bold")
        self.font_btn = tkfont.Font(family="Segoe UI", size=12, weight="bold")
        self.font_result = tkfont.Font(family="Consolas", size=11)
        self.font_console = tkfont.Font(family="Consolas", size=10)

        # Header
        header = tk.Frame(root, bg=self.color_accent, height=54)
        header.pack(fill="x", side="top")
        lbl_title = tk.Label(header, text="Reconnaissance de masse molaire OCSR", bg=self.color_accent, fg="#fff",
                             font=self.font_title, padx=22, pady=11)
        lbl_title.pack(side="left", anchor="w")

        # MAIN BODY
        body = tk.Frame(root, bg=self.bg_main)
        body.pack(fill="both", expand=True)

        # -- ZONE IMAGE GAUCHE --
        left_frame = tk.Frame(body, bg=self.bg_main, width=310)
        left_frame.pack(side="left", fill="y", padx=(18,10), pady=18)
        left_frame.pack_propagate(False)
        self.panel_img = tk.Label(left_frame, bg=self.bg_card, bd=2, relief="groove")
        self.panel_img.pack(pady=14, padx=5, fill="both", expand=True)

        # -- ZONE DROITE (résultats + console) --
        right_frame = tk.Frame(body, bg=self.bg_main)
        right_frame.pack(side="left", fill="both", expand=True, padx=(0,22), pady=18)

        # --- HAUT : BLOC RESULTATS ("carte") ---
        card = tk.Frame(right_frame, bg=self.bg_card, bd=0, highlightbackground="#c7d4e8", highlightthickness=1)
        card.pack(fill="x", expand=False, pady=(0, 9), anchor="n")

        btn_open = tk.Button(
            card, text="📂 Ouvrir une image...", font=self.font_btn,
            bg="#e3edfa", fg="#2856b6", relief="flat", activebackground="#d2e3fa",
            activeforeground="#24447c", cursor="hand2",
            command=self.load_image
        )
        btn_open.pack(pady=8, padx=16, fill='x')

        btn_paste = tk.Button(
            card, text="📋 Coller depuis le presse-papiers (Ctrl+V)", font=self.font_btn,
            bg="#e3edfa", fg="#2856b6", relief="flat", activebackground="#d2e3fa",
            activeforeground="#24447c", cursor="hand2",
            command=self.paste_image
        )
        btn_paste.pack(pady=5, padx=16, fill='x')

        # Zone résultats
        self.text_result = tk.Text(card, width=60, height=11, font=self.font_result,
                                   bg=self.bg_card, fg="#21262c", bd=0, highlightthickness=0)
        self.text_result.pack(pady=8, padx=10)
        self.text_result.tag_configure("rouge", foreground="red")

        # Navigation variants
        nav_frame = tk.Frame(card, bg=self.bg_card)
        nav_frame.pack(pady=(0,5))
        self.btn_prev = tk.Button(
            nav_frame, text="◀", font=self.font_btn, width=3,
            bg="#dde9fc", fg="#3572b0", relief="flat", cursor="hand2",
            command=self.prev_variant
        )
        self.btn_prev.pack(side="left", padx=2)
        self.lbl_num = tk.Label(nav_frame, text="Variante 1/1", font=self.font_btn, bg=self.bg_card, fg="#3b4150")
        self.lbl_num.pack(side="left", padx=3)
        self.btn_next = tk.Button(
            nav_frame, text="▶", font=self.font_btn, width=3,
            bg="#dde9fc", fg="#3572b0", relief="flat", cursor="hand2",
            command=self.next_variant
        )
        self.btn_next.pack(side="left", padx=2)

        btn_copy_smiles = tk.Button(
            card, text="Copier SMILES dans le presse-papiers", font=self.font_btn,
            bg="#d2f2d2", fg="#25772a", activebackground="#bef0c6", relief="flat",
            cursor="hand2", command=self.copy_smiles
        )
        btn_copy_smiles.pack(pady=(7,10), padx=16, fill='x')

        # --- BAS : CONSOLE OCCUPE TOUT L'ESPACE RESTANT ---
        console_frame = tk.LabelFrame(right_frame, text="Console interne", padx=5, pady=5,
                                     bg=self.bg_main, fg="#313749", font=self.font_btn, bd=1, relief="groove",
                                     labelanchor="nw")
        console_frame.pack(fill="both", expand=True, pady=(0,0))
        self.text_console = tk.Text(console_frame, font=self.font_console,
                                    state='normal', bg=self.bg_console, fg=self.fg_console,
                                    insertbackground="white", bd=0, highlightthickness=0)
        self.text_console.pack(fill="both", expand=True)
        # Redirection de la sortie standard et erreurs vers la console interne
        sys.stdout = ConsoleRedirector(self.text_console)
        sys.stderr = ConsoleRedirector(self.text_console)


        # Afficher un log dès la création
        self.log_console("Console initialisée : prêt à analyser vos images.")


        # AJOUTE CECI
    def log_console(self, message):
            self.text_console.config(state='normal')
            self.text_console.insert('end', message + "\n")
            self.text_console.see('end')
            self.text_console.config(state='disabled')





    def process_pil_image(self, pil_img):
        results = []
        variants = [
            ("original", pil_img),
            ("enhanced", enhance_image(pil_img)),
        ]
        for suffix, img in variants:
            try:
                smiles = predict_smiles_with_gradio(img)
                smiles = clean_smiles_isotopes(smiles)
            except Exception as e:
                print(f"Erreur DECIMER distant ({suffix}) : {e}")
                smiles = ""
            mol = Chem.MolFromSmiles(smiles)
            plausible = (
                mol and is_chemically_plausible(mol) and not smiles_has_isotopes(smiles)
            )
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                mass_exact = rdMolDescriptors.CalcExactMolWt(mol)
                mass_moyenne = Descriptors.MolWt(mol)
                img_rdkit = Draw.MolToImage(mol, size=(220, 220))
            else:
                formula = "??"
                mass_exact = None
                mass_moyenne = None
                img_rdkit = None
            results.append({
                "variant": suffix,
                "smiles": smiles,
                "formula": formula,
                "mass_exact": mass_exact,
                "mass_moyenne": mass_moyenne,
                "mol_img": img_rdkit,
                "plausible": plausible,
            })
        print(f"{len(results)} variantes analysées (original + enhanced).")
        return results

    def show_results(self, pil_img):
        print("Analyse de l'image en cours...")
        results = self.process_pil_image(pil_img)
        seen_masses = set()
        eps = 0.01
        variants = []
        for idx, r in enumerate(results):
            if r["mass_exact"] is None or r["smiles"] == "":
                continue
            already = any(abs(r["mass_exact"] - m) < eps for m in seen_masses)
            if not already:
                variants.append(r)
                seen_masses.add(r["mass_exact"])
        self.variants = variants
        self.current_variant = 0
        self.show_current_variant()

    def show_current_variant(self):
        if not self.variants:
            self.panel_img.configure(image='')
            self.text_result.delete(1.0, tk.END)
            self.text_result.insert(tk.END, "Aucune structure reconnue.\n")
            self.update_nav_buttons()
            return
        var = self.variants[self.current_variant]
        if var["mol_img"]:
            img_tk = ImageTk.PhotoImage(var["mol_img"])
            self.panel_img.configure(image=img_tk)
            self.panel_img.image = img_tk
        else:
            self.panel_img.configure(image='')
        self.text_result.delete(1.0, tk.END)
        msg, warning = self.format_variant(var, self.current_variant+1)
        if warning:
            self.text_result.insert(tk.END, msg, "rouge")
        else:
            self.text_result.insert(tk.END, msg)
        self.smiles_last = var["smiles"]
        self.update_nav_buttons()

    def update_nav_buttons(self):
        total = max(1, len(self.variants))
        self.lbl_num.config(text=f"Variante {self.current_variant+1}/{total}")
        if len(self.variants) > 1:
            self.btn_prev.config(state="normal" if self.current_variant > 0 else "disabled")
            self.btn_next.config(state="normal" if self.current_variant < len(self.variants)-1 else "disabled")
        else:
            self.btn_prev.config(state="disabled")
            self.btn_next.config(state="disabled")

    def prev_variant(self):
        if self.variants and self.current_variant > 0:
            self.current_variant -= 1
            self.show_current_variant()

    def next_variant(self):
        if self.variants and self.current_variant < len(self.variants) - 1:
            self.current_variant += 1
            self.show_current_variant()

    def format_variant(self, var, idx):
        msg = f"--- Variante {idx} ---\n"
        msg += f"SMILES: {var['smiles']}\nFormule: {var['formula']}\n"
        msg += f"Masse exacte (monoisotopique via structure): {var['mass_exact']:.5f} g/mol\n"
        msg += f"Masse molaire moyenne (via structure): {var['mass_moyenne']:.2f} g/mol\n"
        mass_from_formula = masse_molaire_from_formula(var["formula"])
        msg += f"Masse molaire moyenne (calculée via formule): {mass_from_formula:.2f} g/mol\n"
        if "_cas" not in var:
            var["_cas"] = smiles_to_cas(var["smiles"])
            if var["_cas"]:
                print(f"CAS trouvé pour {var['smiles']}: {var['_cas']}")
            else:
                print(f"CAS non trouvé pour {var['smiles']}")
        cas_number = var["_cas"]
        msg += f"CAS #: {cas_number if cas_number else '(non trouvé)'}\n"
        if cas_number:
            if "_nom" not in var:
                var["_nom"] = cas_to_name(cas_number)
                if var["_nom"]:
                    print(f"Nom trouvé pour CAS {cas_number}: {var['_nom']}")
                else:
                    print(f"Nom non trouvé pour CAS {cas_number}")
            nom_courant = var["_nom"]
        else:
            nom_courant = ""
        if nom_courant:
            msg += f"Nom PubChem: {nom_courant}\n"
        warning = False
        if not var.get("plausible", True):
            msg += "⚠️ Structure possiblement invalide (valence/isotope anormal)\n"
            warning = True
        msg += "\n"
        return msg, warning

    def load_image(self):
        file_path = filedialog.askopenfilename(filetypes=[
            ("Images", "*.png;*.jpg;*.jpeg;*.bmp;*.gif"),
            ("Tous les fichiers", "*.*"),
        ])
        if not file_path: return
        try:
            pil_img = Image.open(file_path).convert("RGB")
            print(f"Image chargée : {file_path}")
            self.show_results(pil_img)
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de lire l'image:\n{e}")
            print(f"Erreur ouverture image : {e}")

    def paste_image(self):
        try:
            pil_img = ImageGrab.grabclipboard()
            if isinstance(pil_img, Image.Image):
                print("Image collée depuis le presse-papiers.")
                self.show_results(pil_img)
            else:
                messagebox.showwarning("Presse-papiers", "Aucune image trouvée dans le presse-papiers.")
                print("Presse-papiers: aucune image détectée.")
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de récupérer l'image du presse-papiers:\n{e}")
            print(f"Erreur collage image: {e}")

    def copy_smiles(self):
        if self.variants and 0 <= self.current_variant < len(self.variants):
            smiles = self.variants[self.current_variant]["smiles"]
            self.root.clipboard_clear()
            self.root.clipboard_append(smiles)
            self.root.update()
            messagebox.showinfo("SMILES copié", "La chaîne SMILES a été copiée dans le presse-papiers.")
            print(f"SMILES copié: {smiles}")
            self.log_console("Message pour test")
        else:
            messagebox.showwarning("SMILES", "Aucun SMILES à copier.")
            print("Aucun SMILES à copier.")

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
