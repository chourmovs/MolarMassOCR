import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk, ImageGrab
import os
import io

from DECIMER import predict_SMILES
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors

def process_pil_image(pil_img):
    # Sauvegarde temporaire pour décimer
    temp_path = "temp_img.png"
    pil_img.save(temp_path)
    try:
        smiles = predict_SMILES(temp_path)
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None, "Structure non reconnue", None, None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mass = rdMolDescriptors.CalcExactMolWt(mol)
        img = Draw.MolToImage(mol, size=(250, 250))
        return smiles, formula, mass, img
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Reconnaissance de masse molaire OCSR")
        self.root.geometry("650x400")

        self.panel_img = tk.Label(root)
        self.panel_img.pack(side="left", padx=10, pady=10)

        frame = tk.Frame(root)
        frame.pack(side="right", fill="both", expand=True, padx=10)

        btn_open = tk.Button(frame, text="Ouvrir une image...", command=self.load_image)
        btn_open.pack(pady=10)

        btn_paste = tk.Button(frame, text="Coller depuis le presse-papiers (Ctrl+V)", command=self.paste_image)
        btn_paste.pack(pady=10)

        self.text_result = tk.Text(frame, width=45, height=15)
        self.text_result.pack(pady=10)

        self.root.bind('<Control-v>', lambda event: self.paste_image())

    def show_results(self, pil_img):
        result = process_pil_image(pil_img)
        smiles, formula, mass, mol_img = result

        # Affiche la molécule
        if mol_img:
            self.mol_photo = ImageTk.PhotoImage(mol_img)
            self.panel_img.configure(image=self.mol_photo)
            self.panel_img.image = self.mol_photo
        else:
            self.panel_img.configure(image='')

        # Résumé texte
        self.text_result.delete(1.0, tk.END)
        self.text_result.insert(tk.END, f"SMILES: {smiles}\n\nFormule: {formula}\n\nMasse molaire: {mass}\n")

    def load_image(self):
        file_path = filedialog.askopenfilename(filetypes=[
            ("Images", "*.png;*.jpg;*.jpeg;*.bmp;*.gif"),
            ("Tous les fichiers", "*.*"),
        ])
        if not file_path:
            return
        try:
            pil_img = Image.open(file_path).convert("RGB")
            self.show_results(pil_img)
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de lire l'image:\n{e}")

    def paste_image(self):
        try:
            pil_img = ImageGrab.grabclipboard()
            if isinstance(pil_img, Image.Image):
                self.show_results(pil_img)
            else:
                messagebox.showwarning("Presse-papiers", "Aucune image trouvée dans le presse-papiers.")
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de récupérer l'image du presse-papiers:\n{e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
