import tkinter as tk
from tkinter import ttk
import time
import sys
import re

# =========================
# ====== CONFIG APP =======
# =========================
APP_NAME = "Reconnaissance de masse molaire OCSR"
CURRENT_VERSION = "1.08"  # <-- Mets à jour ton numéro de version ici
GITHUB_REPO = "chourmovs/MolarMassOCR"  # owner/repo

# =============== SPLASHSCREEN ===============
class SplashScreen(tk.Toplevel):
    def __init__(self, parent, steps=7):
        super().__init__(parent)
        self.overrideredirect(True)
        w, h = 370, 140
        self.update_idletasks()
        screen_w = self.winfo_screenwidth()
        screen_h = self.winfo_screenheight()
        x = (screen_w // 2) - (w // 2)
        y = (screen_h // 2) - (h // 2)
        self.geometry(f"{w}x{h}+{x}+{y}")
        self.configure(bg="#31436a")
        self.logo = tk.Label(self, text="💡", font=("Segoe UI", 36), bg="#31436a", fg="#fff")
        self.logo.pack(pady=(16, 2))
        self.label = tk.Label(self, text="Chargement...", font=("Segoe UI", 12), bg="#31436a", fg="#fff")
        self.label.pack(pady=(0, 7))
        self.progress = ttk.Progressbar(self, orient='horizontal', length=260, mode='determinate', maximum=steps)
        self.progress.pack(pady=(0, 6))
        self.update_idletasks()
    def set_message(self, text):
        self.label.config(text=text)
        self.update_idletasks()
    def set_progress(self, val):
        self.progress['value'] = val
        self.update_idletasks()

def load_heavy_libs(splash):
    splash.set_message("Importation des modules...")
    time.sleep(0.2)
    import os; splash.set_progress(1)
    time.sleep(0.1)
    import numpy as np; splash.set_progress(2)
    time.sleep(0.1)
    import requests; splash.set_progress(3)
    time.sleep(0.1)
    from rdkit import Chem; splash.set_progress(4)
    time.sleep(0.1)
    from rdkit.Chem import rdMolDescriptors, Draw, Descriptors; splash.set_progress(5)
    time.sleep(0.1)
    import gradio_client as grc; splash.set_progress(6)
    time.sleep(0.1)
    from PIL import Image, ImageTk, ImageGrab, ImageEnhance, ImageFilter; splash.set_progress(7)
    splash.set_message("Prêt !")

# ============== REDIRECTION CONSOLE ==============
class ConsoleRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget
    def write(self, msg):
        self.text_widget.configure(state='normal')
        self.text_widget.insert('end', msg)
        self.text_widget.see('end')
        self.text_widget.configure(state='disabled')
    def flush(self): pass

# ============== POPUP D'ATTENTE (SABLIER) ==============
class BusyPopup(tk.Toplevel):
    """
    Popup modal léger avec progressbar indéterminée + messages verbeux.
    Utiliser show(msg) / update_message(msg) / close().
    """
    def __init__(self, parent, title="Veuillez patienter…"):
        super().__init__(parent)
        self.parent = parent
        self.title(title)
        self.configure(bg="#31436a")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        # Position centré parent
        self.update_idletasks()
        w, h = 420, 160
        px = parent.winfo_rootx() + parent.winfo_width()//2 - w//2
        py = parent.winfo_rooty() + parent.winfo_height()//2 - h//2
        self.geometry(f"{w}x{h}+{px}+{py}")

        # UI
        frm = tk.Frame(self, bg="#31436a")
        frm.pack(expand=True, fill="both", padx=14, pady=14)

        self.lbl_title = tk.Label(frm, text=title, font=("Segoe UI", 12, "bold"), bg="#31436a", fg="#fff")
        self.lbl_title.pack(anchor="w", pady=(2, 6))

        self.txt = tk.Text(frm, height=4, bg="#20304f", fg="#e0e6ef",
                           insertbackground="#fff", bd=0, highlightthickness=0, wrap="word")
        self.txt.pack(fill="both", expand=True, pady=(0, 10))
        self.txt.configure(state="disabled")

        self.pb = ttk.Progressbar(frm, mode="indeterminate", length=360)
        self.pb.pack(fill="x")
        self.pb.start(70)

        # empêcher fermeture (optionnel)
        self.protocol("WM_DELETE_WINDOW", lambda: None)

    def show(self, message: str):
        self.update_message(message)
        self.deiconify()
        self.update()

    def update_message(self, message: str):
        self.txt.configure(state="normal")
        self.txt.insert("end", message.strip() + "\n")
        self.txt.see("end")
        self.txt.configure(state="disabled")
        self.update_idletasks()

    def close(self):
        try:
            self.pb.stop()
        except Exception:
            pass
        self.grab_release()
        self.destroy()

# ============== APPLICATION PRINCIPALE ==============
class App:
    def __init__(self, root):
        # Imports dynamiques
        import os, requests, tempfile, threading, webbrowser, queue
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors, Draw, Descriptors
        import gradio_client as grc
        from gradio_client import handle_file
        from PIL import Image, ImageTk, ImageGrab, ImageEnhance, ImageFilter
        from tkinter import filedialog, messagebox, font as tkfont

        self.os = os
        self.requests = requests
        self.threading = threading
        self.webbrowser = webbrowser
        self.queue = queue
        self.np = np
        self.Chem = Chem
        self.rdMolDescriptors = rdMolDescriptors
        self.Draw = Draw
        self.Descriptors = Descriptors
        self.grc = grc
        self.handle_file = handle_file
        self.Image = Image
        self.ImageTk = ImageTk
        self.ImageGrab = ImageGrab
        self.ImageEnhance = ImageEnhance
        self.ImageFilter = ImageFilter
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.tkfont = tkfont
        self.tempfile = tempfile

        # Client Gradio distant
        self.hf_space_url = "https://chouchouvs-distance-smiles.hf.space/"
        self.client = grc.Client(self.hf_space_url)

        # poids atomiques
        self.atomic_weights = {
            "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
            "S": 32.06, "P": 30.974, "F": 18.998, "Cl": 35.45,
            "Br": 79.904, "I": 126.90,
        }

        # --- UI Tkinter ---
        self.root = root
        self.root.title(f"{APP_NAME} — v{CURRENT_VERSION}")
        self.root.geometry("1000x650")
        self.root.minsize(1000, 600)
        self.bg_main = "#f5f8fc"
        self.color_accent = "#4F8EF7"
        self.bg_card = "#ffffff"
        self.bg_console = "#1a2233"
        self.fg_console = "#e0e6ef"
        self.font_title = tkfont.Font(family="Segoe UI", size=17, weight="bold")
        self.font_btn = tkfont.Font(family="Segoe UI", size=12, weight="bold")
        self.font_result = tkfont.Font(family="Consolas", size=11)
        self.font_console = tkfont.Font(family="Consolas", size=10)

        # Menu
        self._build_menu()

        self.zoom_factor = 1.0
        self.min_zoom = 0.2
        self.max_zoom = 4.0
        self.img_pil_displayed = None
        self.img_tk = None
        self.pan_x = 0
        self.pan_y = 0
        self.drag_start_x = None
        self.drag_start_y = None

        # Header
        header = tk.Frame(root, bg=self.color_accent, height=54)
        header.pack(fill="x", side="top")
        lbl_title = tk.Label(header, text=APP_NAME, bg=self.color_accent, fg="#fff",
                             font=self.font_title, padx=22, pady=11)
        lbl_title.pack(side="left", anchor="w")

        # MAIN BODY
        body = tk.Frame(root, bg=self.bg_main)
        body.pack(fill="both", expand=True)
        left_frame = tk.Frame(body, bg=self.bg_main, width=410)
        left_frame.pack(side="left", fill="y", padx=(18,10), pady=18)
        left_frame.pack_propagate(False)

        self.image_container = tk.Frame(left_frame, bg=self.bg_main, width=400, height=300)
        self.image_container.pack(fill="both", expand=True)
        self.image_container.pack_propagate(False)

        self.btn_overlay = tk.Frame(self.image_container, bg='', highlightthickness=0)
        self.btn_overlay.place(relx=1.0, rely=0.0, anchor='ne', x=-2, y=2)

        btn_style = dict(
            font=("Segoe UI", 10, "bold"),
            relief="flat",
            bg="#e3edfa", fg="#3572b0",
            activebackground="#d2e3fa",
            activeforeground="#24447c",
            bd=0,
            cursor="hand2",
            padx=2, pady=0,
        )
        self.btn_hd = tk.Button(
            self.btn_overlay, text="🔎", command=self.show_full_res, **btn_style
        )
        self.btn_hd.pack(side="right", padx=(2,0))
        self.btn_reset_zoom = tk.Button(
            self.btn_overlay, text="⟳", command=self.reset_zoom, **btn_style
        )
        self.btn_reset_zoom.pack(side="right", padx=0)

        self.canvas_img = tk.Canvas(self.image_container, bg=self.bg_card, bd=2, relief="groove", highlightthickness=0)
        self.canvas_img.pack(pady=(25,10), padx=5, fill="both", expand=True)

        self.canvas_img.bind("<MouseWheel>", self.on_mousewheel_zoom)
        self.canvas_img.bind("<Button-4>", self.on_mousewheel_zoom_linux)
        self.canvas_img.bind("<Button-5>", self.on_mousewheel_zoom_linux)
        self.canvas_img.bind("<ButtonPress-1>", self.start_pan)
        self.canvas_img.bind("<B1-Motion>", self.do_pan)
        self.canvas_img.bind("<ButtonRelease-1>", self.end_pan)

        right_frame = tk.Frame(body, bg=self.bg_main)
        right_frame.pack(side="left", fill="both", expand=True, padx=(0,22), pady=18)

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

        self.text_result = tk.Text(card, width=60, height=11, font=self.font_result,
                                   bg=self.bg_card, fg="#21262c", bd=0, highlightthickness=0)
        self.text_result.pack(pady=8, padx=10)
        self.text_result.tag_configure("rouge", foreground="red")

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

        console_frame = tk.LabelFrame(right_frame, text="Console interne", padx=5, pady=5,
                                     bg=self.bg_main, fg="#313749", font=self.font_btn, bd=1, relief="groove",
                                     labelanchor="nw")
        console_frame.pack(fill="both", expand=True, pady=(0,0))
        self.text_console = tk.Text(console_frame, font=self.font_console,
                                    state='normal', bg=self.bg_console, fg=self.fg_console,
                                    insertbackground="white", bd=0, highlightthickness=0)
        self.text_console.pack(fill="both", expand=True)
        sys.stdout = ConsoleRedirector(self.text_console)
        sys.stderr = ConsoleRedirector(self.text_console)
        self.log_console(f"Console initialisée. Version: {CURRENT_VERSION}")

        # Vérification de mise à jour au démarrage (non bloquant)
        self.root.after(800, self._auto_check_update)

    # ---------- MENU ----------
    def _build_menu(self):
        m = tk.Menu(self.root)
        self.root.config(menu=m)

        menu_file = tk.Menu(m, tearoff=0)
        menu_file.add_command(label="Quitter", command=self.root.quit)
        m.add_cascade(label="Fichier", menu=menu_file)

        menu_help = tk.Menu(m, tearoff=0)
        menu_help.add_command(label="Rechercher des mises à jour…", command=self.check_updates_interactive)
        menu_help.add_command(label="Ouvrir la page GitHub", command=lambda: self._open_url(f"https://github.com/{GITHUB_REPO}"))
        menu_help.add_separator()
        menu_help.add_command(label="À propos", command=self._about)
        m.add_cascade(label="Aide", menu=menu_help)

    def _about(self):
        from tkinter import messagebox
        messagebox.showinfo("À propos", f"{APP_NAME}\nVersion {CURRENT_VERSION}\n\nRepo: {GITHUB_REPO}")

    def _open_url(self, url: str):
        try:
            self.webbrowser.open(url)
        except Exception as e:
            self.messagebox.showerror("Erreur", f"Impossible d'ouvrir le navigateur:\n{e}")

    def log_console(self, message):
        self.text_console.config(state='normal')
        self.text_console.insert('end', message + "\n")
        self.text_console.see('end')
        self.text_console.config(state='disabled')

    # ----------- ZOOM + PAN + HD --------------------
    def update_panel_img(self):
        self.canvas_img.delete("all")
        if self.img_pil_displayed is None:
            return
        pil = self.img_pil_displayed
        w, h = pil.size
        z = self.zoom_factor
        w_z, h_z = int(w * z), int(h * z)
        img_zoomed = pil.resize((w_z, h_z), self.Image.LANCZOS)
        self.img_tk = self.ImageTk.PhotoImage(img_zoomed)
        c_width = self.canvas_img.winfo_width()
        c_height = self.canvas_img.winfo_height()
        x = c_width // 2 - w_z // 2 + self.pan_x
        y = c_height // 2 - h_z // 2 + self.pan_y
        self.canvas_img.create_image(x, y, anchor="nw", image=self.img_tk)

    def on_mousewheel_zoom(self, event):
        if event.delta > 0:
            self.zoom_factor = min(self.zoom_factor * 1.13, self.max_zoom)
        else:
            self.zoom_factor = max(self.zoom_factor / 1.13, self.min_zoom)
        self.update_panel_img()

    def on_mousewheel_zoom_linux(self, event):
        if event.num == 4:
            self.zoom_factor = min(self.zoom_factor * 1.13, self.max_zoom)
        elif event.num == 5:
            self.zoom_factor = max(self.zoom_factor / 1.13, self.min_zoom)
        self.update_panel_img()

    def reset_zoom(self):
        self.zoom_factor = 1.0
        self.pan_x = 0
        self.pan_y = 0
        self.update_panel_img()

    def start_pan(self, event):
        self.drag_start_x = event.x
        self.drag_start_y = event.y

    def do_pan(self, event):
        dx = event.x - self.drag_start_x
        dy = event.y - self.drag_start_y
        self.pan_x += dx
        self.pan_y += dy
        self.drag_start_x = event.x
        self.drag_start_y = event.y
        self.update_panel_img()

    def end_pan(self, event):
        self.drag_start_x = None
        self.drag_start_y = None

    def show_full_res(self):
        if self.img_pil_displayed is None:
            return
        win = tk.Toplevel(self.root)
        win.title("Structure - Haute résolution")
        pil = self.img_pil_displayed
        img_tk = self.ImageTk.PhotoImage(pil)
        lbl = tk.Label(win, image=img_tk)
        lbl.image = img_tk
        lbl.pack()
        win.geometry(f"{pil.width}x{pil.height}")

    # ----------- LOGIQUE CHIMIE/OCSR ------------
    def enhance_image(self, pil_img):
        img = self.ImageEnhance.Contrast(pil_img).enhance(1.3)
        img = self.ImageEnhance.Sharpness(img).enhance(2.0)
        img = img.filter(self.ImageFilter.MedianFilter(size=3))
        return img

    def masse_molaire_from_formula(self, formula):
        pattern = r"([A-Z][a-z]?)(\d*)"
        total = 0.0
        for (elt, count) in re.findall(pattern, formula):
            count = int(count) if count else 1
            w = self.atomic_weights.get(elt)
            if w is None: continue
            total += w * count
        return total

    def is_chemically_plausible(self, mol):
        if not mol: return False
        for atom in mol.GetAtoms():
            s = atom.GetSymbol()
            v = atom.GetExplicitValence()
            if s == "C" and v > 4: return False
            if s == "O" and v > 2: return False
            if s == "N" and v > 3: return False
            if s == "H" and v > 1: return False
        return True

    def smiles_has_isotopes(self, smiles):
        return bool(re.search(r'\[\d+[A-Z][a-z]?\]', smiles))

    def smiles_to_cas(self, smiles):
        try:
            url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
            r = self.requests.get(url_cid, timeout=10)
            if not r.ok or not r.text.strip().isdigit():
                return None
            cid = r.text.strip()
            url_syn = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
            r2 = self.requests.get(url_syn, timeout=10)
            if not r2.ok: return None
            data = r2.json()
            ids = data.get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
            for regid in ids:
                if re.match(r'^\d{2,7}-\d{2}-\d$', regid):
                    return regid
            return None
        except Exception:
            return None

    def cas_to_name(self, cas):
        try:
            url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{cas}/cids/TXT"
            r = self.requests.get(url_cid, timeout=10)
            if not r.ok or not r.text.strip().isdigit():
                return None
            cid = r.text.strip()
            url_name = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,Title/JSON"
            r2 = self.requests.get(url_name, timeout=10)
            if not r2.ok: return None
            data = r2.json()
            props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
            return props.get("Title") or props.get("IUPACName")
        except Exception:
            return None

    def clean_smiles_isotopes(self, smiles):
        return re.sub(r'\[(\d+)([A-Z][a-z]?)\]', r'\2', smiles)

    # --------- Appels HF AVEC SABLIER & RETRIES (thread) ---------
    def predict_smiles_with_gradio(self, pil_img, verbose_cb=None, max_retries=3):
        """
        Lance l'appel Gradio en thread et bloque jusqu'au résultat mais
        permet d'envoyer des messages via verbose_cb (GUI non gelée grâce au popup).
        """
        # Sauvegarde image dans un fichier temp
        with self.tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            pil_img.save(tmp, format="PNG")
            tmp_path = tmp.name

        def v(msg):
            if verbose_cb: verbose_cb(msg)
            print(msg)

        # Petit ping pour réveiller l'espace (si endormi)
        try:
            v("Ping de l'espace HuggingFace…")
            self.client.view_api(all_endpoints=False)  # appel léger
        except Exception as e:
            v(f"Ping non déterminant (continue) : {e}")

        # Tentatives avec backoff progressif
        last_err = None
        for attempt in range(1, max_retries + 1):
            try:
                v(f"Tentative {attempt}/{max_retries} : envoi de l'image…")
                res = self.client.predict(
                    img=self.handle_file(tmp_path),
                    api_name="/predict"
                )
                v("Réponse reçue ✅")
                return res
            except Exception as e:
                last_err = e
                v(f"Erreur HF : {e}")
                if attempt < max_retries:
                    delay = 2 * attempt
                    v(f"Nouvelle tentative dans {delay}s (l'espace peut être en réveil)…")
                    time.sleep(delay)
        # Échec après retries
        raise RuntimeError(f"Échec d'appel HuggingFace après {max_retries} tentatives : {last_err}")

    # ----- Pipeline de traitement avec BusyPopup -----
    def process_pil_image(self, pil_img):
        results = []
        variants = [
            ("original", pil_img),
            ("enhanced", self.enhance_image(pil_img)),
        ]
        for suffix, img in variants:
            try:
                # Affiche popup sablier et exécute en thread
                popup = BusyPopup(self.root, title="Analyse en cours…")
                popup.show(f"Préparation de la variante: {suffix}")
                # pour ne pas bloquer : on utilise une Queue pour rapatrier résultat/erreur
                q = self.queue.Queue()

                def worker():
                    try:
                        smiles = self.predict_smiles_with_gradio(
                            img,
                            verbose_cb=lambda m: self.root.after(0, popup.update_message, m),
                            max_retries=3
                        )
                        q.put(("ok", smiles))
                    except Exception as e:
                        q.put(("err", e))

                t = self.threading.Thread(target=worker, daemon=True)
                t.start()

                # Boucle d'attente non bloquante (garde l'UI vivante)
                while t.is_alive():
                    self.root.update()
                    time.sleep(0.05)

                status, payload = q.get_nowait()
                popup.close()

                if status == "ok":
                    smiles = self.clean_smiles_isotopes(payload)
                else:
                    print(f"Erreur DECIMER distant ({suffix}) : {payload}")
                    smiles = ""
            except Exception as e:
                print(f"Erreur DECIMER ({suffix}) : {e}")
                smiles = ""

            mol = self.Chem.MolFromSmiles(smiles)
            plausible = (mol and self.is_chemically_plausible(mol) and not self.smiles_has_isotopes(smiles))
            if mol:
                formula = self.rdMolDescriptors.CalcMolFormula(mol)
                mass_exact = self.rdMolDescriptors.CalcExactMolWt(mol)
                mass_moyenne = self.Descriptors.MolWt(mol)
                img_rdkit = self.Draw.MolToImage(mol, size=(440, 440))
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
        for r in results:
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
        if not hasattr(self, "variants") or not self.variants:
            self.img_pil_displayed = None
            self.update_panel_img()
            self.text_result.delete(1.0, tk.END)
            self.text_result.insert(tk.END, "Aucune structure reconnue.\n")
            self.update_nav_buttons()
            return
        var = self.variants[self.current_variant]
        if var["mol_img"]:
            self.img_pil_displayed = var["mol_img"]
            self.pan_x = 0
            self.pan_y = 0
            self.zoom_factor = 1.0
            self.update_panel_img()
        else:
            self.img_pil_displayed = None
            self.update_panel_img()
        self.text_result.delete(1.0, tk.END)
        msg, warning = self.format_variant(var, self.current_variant+1)
        if warning:
            self.text_result.insert(tk.END, msg, "rouge")
        else:
            self.text_result.insert(tk.END, msg)
        self.smiles_last = var["smiles"]
        self.update_nav_buttons()

    def update_nav_buttons(self):
        total = max(1, len(getattr(self, "variants", [])))
        self.lbl_num.config(text=f"Variante {self.current_variant+1}/{total}")
        if total > 1:
            self.btn_prev.config(state="normal" if self.current_variant > 0 else "disabled")
            self.btn_next.config(state="normal" if self.current_variant < total-1 else "disabled")
        else:
            self.btn_prev.config(state="disabled")
            self.btn_next.config(state="disabled")

    def prev_variant(self):
        if getattr(self, "variants", None) and self.current_variant > 0:
            self.current_variant -= 1
            self.show_current_variant()

    def next_variant(self):
        if getattr(self, "variants", None) and self.current_variant < len(self.variants) - 1:
            self.current_variant += 1
            self.show_current_variant()

    def format_variant(self, var, idx):
        msg = f"--- Variante {idx} ---\n"
        msg += f"SMILES: {var['smiles']}\nFormule: {var['formula']}\n"
        msg += f"Masse exacte (monoisotopique via structure): {var['mass_exact']:.5f} g/mol\n"
        msg += f"Masse molaire moyenne (via structure): {var['mass_moyenne']:.2f} g/mol\n"
        mass_from_formula = self.masse_molaire_from_formula(var["formula"])
        msg += f"Masse molaire moyenne (calculée via formule): {mass_from_formula:.2f} g/mol\n"
        if "_cas" not in var:
            var["_cas"] = self.smiles_to_cas(var["smiles"])
            if var["_cas"]:
                print(f"CAS trouvé pour {var['smiles']}: {var['_cas']}")
            else:
                print(f"CAS non trouvé pour {var['smiles']}")
        cas_number = var["_cas"]
        msg += f"CAS #: {cas_number if cas_number else '(non trouvé)'}\n"
        if cas_number:
            if "_nom" not in var:
                var["_nom"] = self.cas_to_name(cas_number)
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

    # ----------- Actions UI chargement/coller/copie ----------
    def load_image(self):
        file_path = self.filedialog.askopenfilename(filetypes=[
            ("Images", "*.png;*.jpg;*.jpeg;*.bmp;*.gif"),
            ("Tous les fichiers", "*.*"),
        ])
        if not file_path: return
        try:
            pil_img = self.Image.open(file_path).convert("RGB")
            print(f"Image chargée : {file_path}")
            self.show_results(pil_img)
        except Exception as e:
            self.messagebox.showerror("Erreur", f"Impossible de lire l'image:\n{e}")
            print(f"Erreur ouverture image : {e}")

    def paste_image(self):
        try:
            pil_img = self.ImageGrab.grabclipboard()
            if isinstance(pil_img, self.Image.Image):
                print("Image collée depuis le presse-papiers.")
                self.show_results(pil_img)
            else:
                self.messagebox.showwarning("Presse-papiers", "Aucune image trouvée dans le presse-papiers.")
                print("Presse-papiers: aucune image détectée.")
        except Exception as e:
            self.messagebox.showerror("Erreur", f"Impossible de récupérer l'image du presse-papiers:\n{e}")
            print(f"Erreur collage image: {e}")

    def copy_smiles(self):
        if hasattr(self, 'variants') and self.variants and 0 <= self.current_variant < len(self.variants):
            smiles = self.variants[self.current_variant]["smiles"]
            self.root.clipboard_clear()
            self.root.clipboard_append(smiles)
            self.root.update()
            self.messagebox.showinfo("SMILES copié", "La chaîne SMILES a été copiée dans le presse-papiers.")
            print(f"SMILES copié: {smiles}")
            self.log_console("Message pour test")
        else:
            self.messagebox.showwarning("SMILES", "Aucun SMILES à copier.")
            print("Aucun SMILES à copier.")

    # =======================
    #   MISES À JOUR GITHUB
    # =======================
    def _parse_version(self, s: str):
        """
        Parse simple SemVer-like strings: 'v1.2.3' or '1.2.3'
        Retourne un tuple d'entiers pour comparaison.
        """
        s = s.strip().lower()
        if s.startswith("v"):
            s = s[1:]
        parts = re.split(r"[^0-9]+", s)
        nums = [int(p) for p in parts if p.isdigit()]
        while len(nums) < 3:
            nums.append(0)
        return tuple(nums[:3])

    def _fetch_latest_release(self):
        url = f"https://api.github.com/repos/{GITHUB_REPO}/releases/latest"
        headers = {"User-Agent": f"{APP_NAME}-updater"}
        r = self.requests.get(url, headers=headers, timeout=8)
        r.raise_for_status()
        data = r.json()
        tag = data.get("tag_name") or data.get("name") or ""
        html_url = data.get("html_url", f"https://github.com/{GITHUB_REPO}/releases")
        return tag, html_url

    def check_updates_interactive(self):
        try:
            tag, page = self._fetch_latest_release()
            v_remote = self._parse_version(tag)
            v_local = self._parse_version(CURRENT_VERSION)
            if v_remote > v_local:
                if self.messagebox.askyesno(
                    "Mise à jour disponible",
                    f"Nouvelle version détectée : {tag}\n\nVersion actuelle : {CURRENT_VERSION}\n\n"
                    f"Voulez-vous ouvrir la page des releases ?"
                ):
                    self._open_url(page)
            else:
                self.messagebox.showinfo("À jour", f"Aucune mise à jour disponible.\nVersion actuelle : {CURRENT_VERSION}")
        except Exception as e:
            self.messagebox.showwarning("Vérification échouée", f"Impossible de vérifier les mises à jour:\n{e}")

    def _auto_check_update(self):
        def worker():
            try:
                tag, page = self._fetch_latest_release()
                v_remote = self._parse_version(tag)
                v_local = self._parse_version(CURRENT_VERSION)
                if v_remote > v_local:
                    self.root.after(0, lambda: self._prompt_update(tag, page))
            except Exception as e:
                # silencieux au démarrage ; log console
                self.log_console(f"[MAJ] Vérification silencieuse échouée : {e}")

        t = self.threading.Thread(target=worker, daemon=True)
        t.start()

    def _prompt_update(self, tag: str, page: str):
        try:
            if self.messagebox.askyesno(
                "Mise à jour disponible",
                f"Une nouvelle version est disponible : {tag}\n"
                f"Version actuelle : {CURRENT_VERSION}\n\n"
                "Ouvrir la page GitHub pour télécharger ?"
            ):
                self._open_url(page)
        except Exception as e:
            self.log_console(f"[MAJ] Prompt update error: {e}")

# ================= MAIN =================
if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    splash = SplashScreen(root, steps=7)
    splash.update()
    load_heavy_libs(splash)
    app = App(root)  # Création de l'UI, fenêtre encore cachée

    def show_main_and_kill_splash():
        root.deiconify()
        splash.destroy()

    root.after_idle(show_main_and_kill_splash)
    root.mainloop()

