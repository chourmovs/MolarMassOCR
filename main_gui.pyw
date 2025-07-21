import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk, ImageGrab
import glob
import re

# === Redirection de la console ===
class ConsoleRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, msg):
        self.text_widget.configure(state='normal')
        self.text_widget.insert('end', msg)
        self.text_widget.see('end')
        self.text_widget.configure(state='disabled')

    def flush(self):
        pass  # pour compatibilité

# === Fixe le chemin des modèles AVANT import des modules IA ===
MODEL_DIR = os.path.join(os.path.dirname(__file__), "models")
os.environ["DECIMER_CACHE_DIR"] = os.path.join(MODEL_DIR, "decimer")
#scale_str = f"x{scale}plus"
#model_path = os.path.join(MODEL_DIR, f"RealESRGAN_{scale_str}.pth")


# === Splashscreen immédiat, centré et large ===
splash = tk.Tk()
splash.overrideredirect(True)
splash_width, splash_height = 640, 180
screen_width = splash.winfo_screenwidth()
screen_height = splash.winfo_screenheight()
x = (screen_width // 2) - (splash_width // 2)
y = (screen_height // 2) - (splash_height // 2)
splash.geometry(f"{splash_width}x{splash_height}+{x}+{y}")
splash.config(bg="#334466")
label = tk.Label(
    splash,
    text="Reconnaissance de structure moléculaire\nChargement du moteur IA...",
    font=("Arial", 20, "bold"),
    padx=30, pady=45,
    bg="#334466", fg="white"
)
label.pack(expand=True, fill="both")
splash.update()

# Imports lents ici (après DECIMER_CACHE_DIR)
from DECIMER import predict_SMILES
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw, Descriptors
import matplotlib
matplotlib.use('Agg')
import cv2
import numpy as np
import requests

# Super-résolution (RealESRGAN doit être dans models/)
try:
    from realesrgan import RealESRGAN
    import torch
    has_sr = True
except ImportError:
    has_sr = False

splash.destroy()

# ---- Nettoyage SMILES des isotopes ----
def clean_smiles_isotopes(smiles):
    return re.sub(r'\[(\d+)([A-Z][a-z]?)\]', r'\2', smiles)

# ---- Super-résolution sans téléchargement, modèle local obligatoire ----
def super_resolve_pil(pil_img, scale=2):
    if not has_sr:
        print("RealESRGAN non disponible.")
        return pil_img
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    try:
        model = RealESRGAN(device, scale=scale)
        scale_str = f"x{scale}plus"
        model_path = os.path.join(MODEL_DIR, f"RealESRGAN_{scale_str}.pth")
        if not os.path.exists(model_path):
            print(f"Fichier modèle RealESRGAN manquant : {model_path}")
            return pil_img
        model.load_weights(model_path)
        sr_img = model.predict(pil_img)
        print(f"Super-résolution {scale_str} appliquée.")
        return sr_img
    except Exception as e:
        print("Erreur Real-ESRGAN :", e)
        return pil_img


# ---- Utilitaires chimiques ----
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
        if w is None:
            continue
        total += w * count
    return total

def is_chemically_plausible(mol):
    if not mol:
        return False
    for atom in mol.GetAtoms():
        s = atom.GetSymbol()
        v = atom.GetExplicitValence()
        if s == "C" and v > 4:
            return False
        if s == "O" and v > 2:
            return False
        if s == "N" and v > 3:
            return False
        if s == "H" and v > 1:
            return False
    return True

def smiles_has_isotopes(smiles):
    return bool(re.search(r'\[\d+[A-Z][a-z]?\]', smiles))

# ---- Recherche CAS/nom via PubChem ----
def smiles_to_cas(smiles):
    try:
        url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
        r = requests.get(url_cid, timeout=10)
        if not r.ok or not r.text.strip().isdigit():
            return None
        cid = r.text.strip()

        url_syn = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
        r2 = requests.get(url_syn, timeout=10)
        if not r2.ok:
            return None

        data = r2.json()
        ids = data.get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
        for regid in ids:
            if re.match(r'^\d{2,7}-\d{2}-\d$', regid):
                return regid
        return None
    except Exception:
        return None

def cas_to_name(cas):
    try:
        url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{cas}/cids/TXT"
        r = requests.get(url_cid, timeout=10)
        if not r.ok or not r.text.strip().isdigit():
            return None
        cid = r.text.strip()

        url_name = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,Title/JSON"
        r2 = requests.get(url_name, timeout=10)
        if not r2.ok:
            return None
        data = r2.json()
        props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
        return props.get("Title") or props.get("IUPACName")
    except Exception:
        return None

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Reconnaissance de masse molaire OCSR")
        self.root.geometry("790x540")
        self.root.resizable(False, False)

        left_frame = tk.Frame(root)
        left_frame.pack(side="left", padx=10, pady=10)

        self.panel_img = tk.Label(left_frame)
        self.panel_img.pack()

        right_frame = tk.Frame(root)
        right_frame.pack(side="right", fill="both", expand=True, padx=10)

        btn_open = tk.Button(right_frame, text="Ouvrir une image...", font=("Arial", 12), command=self.load_image)
        btn_open.pack(pady=6, fill='x')

        btn_paste = tk.Button(right_frame, text="Coller depuis le presse-papiers (Ctrl+V)", font=("Arial", 12), command=self.paste_image)
        btn_paste.pack(pady=6, fill='x')

        self.text_result = tk.Text(right_frame, width=56, height=12, font=("Consolas", 10))
        self.text_result.pack(pady=8)
        self.text_result.tag_configure("rouge", foreground="red")

        # === Console live (stdout/stderr) ===
        console_frame = tk.LabelFrame(right_frame, text="Console interne", padx=2, pady=2)
        console_frame.pack(fill="both", expand=False, pady=(0,6))
        self.text_console = tk.Text(console_frame, height=7, font=("Consolas", 9), state='disabled', bg="#21262c", fg="#e6e6e6")
        self.text_console.pack(fill="both", expand=True)
        sys.stdout = ConsoleRedirector(self.text_console)
        sys.stderr = ConsoleRedirector(self.text_console)
        print("Console interne activée. Tous les print() apparaîtront ici.")

        nav_frame = tk.Frame(right_frame)
        nav_frame.pack(pady=2)
        self.btn_prev = tk.Button(nav_frame, text="◀ Précédent", command=self.prev_variant)
        self.btn_prev.pack(side="left", padx=3)
        self.lbl_num = tk.Label(nav_frame, text="Variante 1/1", width=15)
        self.lbl_num.pack(side="left", padx=3)
        self.btn_next = tk.Button(nav_frame, text="Suivant ▶", command=self.next_variant)
        self.btn_next.pack(side="left", padx=3)

        btn_copy_smiles = tk.Button(right_frame, text="Copier SMILES dans le presse-papiers", font=("Arial", 11), command=self.copy_smiles)
        btn_copy_smiles.pack(pady=3, fill='x')

        self.variants = []
        self.current_variant = 0
        self.smiles_last = None

        self.update_nav_buttons()

        self.root.bind('<Control-v>', lambda event: self.paste_image())

    def process_pil_image(self, pil_img):
        print("Prétraitement de l'image (super-résolution si dispo)...")
        pil_img = super_resolve_pil(pil_img, scale=2)
        params = [
           # (1.0, 21, 4, 3, 'adaptive', "_A"),
            #(1.5, 21, 4, 3, 'adaptive', "_B"),
            #(2.0, 21, 4, 3, 'adaptive', "_C"),
            #(1.0, 31, 6, 5, 'adaptive', "_D"),
            #(1.5, 11, 2, 1, 'adaptive', "_E"),
            (1.0, 0, 0, 0, 'otsu', "_F"),
        ]

        temp_base = "temp_img"
        pil_img.save(f"{temp_base}.png")
        results = []
        try:
            for scale, blocksize, C, blur, method, suffix in params:
                img_cv = cv2.cvtColor(np.array(pil_img), cv2.COLOR_RGB2BGR)
                img = cv2.cvtColor(img_cv, cv2.COLOR_BGR2GRAY)
                if scale != 1.0:
                    img = cv2.resize(img, None, fx=scale, fy=scale, interpolation=cv2.INTER_CUBIC)
                if method == 'adaptive':
                    proc = cv2.adaptiveThreshold(
                        img, 255, cv2.ADAPTIVE_THRESH_MEAN_C,
                        cv2.THRESH_BINARY, blocksize, C)
                elif method == 'otsu':
                    _, proc = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
                else:
                    proc = img
                if blur > 0:
                    proc = cv2.medianBlur(proc, blur)
                preproc_path = f"{temp_base}_prep{suffix}.png"
                cv2.imwrite(preproc_path, proc)

                try:
                    smiles = predict_SMILES(preproc_path)
                    smiles = clean_smiles_isotopes(smiles)
                except Exception as e:
                    print(f"Erreur DECIMER (variante {suffix}): {e}")
                    smiles = ""
                mol = Chem.MolFromSmiles(smiles)
                plausible = (
                    mol and
                    is_chemically_plausible(mol) and
                    not smiles_has_isotopes(smiles)
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
                    "smiles": smiles,
                    "formula": formula,
                    "mass_exact": mass_exact,
                    "mass_moyenne": mass_moyenne,
                    "mol_img": img_rdkit,
                    "plausible": plausible,
                })
        finally:
            for f in glob.glob(f"{temp_base}*.png"):
                try:
                    os.remove(f)
                except Exception as e:
                    print(f"Erreur suppression temporaire {f}: {e}")
        print(f"{len(results)} variantes analysées.")
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

        # Recherche automatique du CAS et nom courant
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
        if not file_path:
            return
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
        else:
            messagebox.showwarning("SMILES", "Aucun SMILES à copier.")
            print("Aucun SMILES à copier.")

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()

