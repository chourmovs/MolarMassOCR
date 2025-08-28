# 🧪 MolarMassOCR – OCR de structures chimiques → masses molaires
---
## 🎯 Objet
MolarMassOCR est une application de bureau (Tkinter) qui :
lit une image (fichier ou presse-papiers) contenant une structure chimique,
l’envoie à un modèle OCSR distant pour obtenir un SMILES,
calcule la formule, la masse exacte (monoisotopique) et la masse molaire moyenne via RDKit,
effectue une vérification simple de plausibilité (valences, isotopes),
tente de retrouver un n° CAS et un nom via PubChem,
affiche le rendu moléculaire pour une vérification simple, propose zoom/pan, navigation entre variantes (originale & améliorée), et copie SMILES.
L’outil est interactif, sans persistance locale par défaut, et expose une console interne pour la traçabilité.


## 👥 Public visé
Utilisateurs : chimistes, analystes, ingénieurs procédés.


## 🖥️ Utilisation (utilisateurs)
```
Lancer l’application (python app.py ou binaire packagé si fourni).
Importer une structure :
📂 Ouvrir une image (PNG/JPG/BMP/GIF), ou
📋 Coller depuis le presse-papiers.
L’outil génère 2 variantes (originale + enhanced : contraste/netteté/filtre médian) et tente une reconnaissance pour chaque.
Naviguer entre variantes (◀ ▶), zoom/pan, vue HD.
Copier le SMILES en un clic.
Consulter les résultats :
SMILES, formule brute, masse exacte, masse molaire moyenne (structure & formule),
CAS (si trouvé), nom PubChem (si trouvé),
alerte si plausibilité douteuse (valences/isotopes).
Important : un message “Aucune structure reconnue” s’affiche si rien n’est exploitable.
````
---

## 🔌 Dépendances et réseau (réalité du code)
```
OCSR distant (Gradio) : https://chouchouvs-distance-smiles.hf.space/
➜ Envoi de l’image (PNG temporaire) pour obtenir le SMILES.<br>
Résolution CAS / Nom : PubChem REST (pubchem.ncbi.nlm.nih.gov)
➜ Envoi du SMILES ou d’un CAS pour récupérer CID/nom.<br>
Mode offline : non supporté en l’état (l’OCSR et PubChem nécessitent Internet).<br>
Pour un mode offline, remplacer l’endpoint OCSR par un modèle local et désactiver/respecter les appels PubChem.
```

## 🧩 Pile technique
```
UI : Tkinter (splash, console interne, canvas avec zoom/pan).
Vision & chimie : RDKit (formule, masses, rendu MolToImage).
Amélioration image : PIL (contraste, netteté, filtre médian).
Réseau : gradio_client (endpoint Hugging Face Space), requests (PubChem).
Sécurité UX : aucune écriture disque par défaut (hors fichiers temporaires image → envoi OCSR).
```

## 🔒 Sécurité digitale (focus IT & utilisateurs)
```
Ce qui sort de la machine
Vers l’endpoint OCSR (*.hf.space) :
Image de la structure (PNG) — variante originale ou « enhanced ».
Pas d’autres métadonnées explicites envoyées par l’app.
Vers PubChem :
SMILES (pour récupérer un CID), ou
CAS (pour retrouver un CID/nom).

Ce qui n’est pas fait
Pas d’upload automatique de jeux de données hors action utilisateur.
Pas de stockage local persistant par défaut (logs uniquement dans la console interne).
Pas de clés/API sensibles utilisées par défaut dans ce projet.
````

## Recommandations IT
```
Réseau : autoriser HTTPS sortant vers :
*.hf.space (Gradio/Hugging Face Space)
pubchem.ncbi.nlm.nih.gov
Proxy : si nécessaire, configurer les variables d’environnement (HTTP_PROXY, HTTPS_PROXY, NO_PROXY).
### Isolement : pour documents sensibles, envisager un self-host du modèle OCSR en interne (ou bloquer l’upload externe) et substituer PubChem par une base locale.
Traçabilité : si besoin, rediriger la console vers un logger applicatif chiffré et sans contenu des images.
```

## Risques & mitigation
```
Confidentialité IP (images de schémas, structures propriétaires) :
➜ masquer/recadrer précisément les zones nécessaires avant import,
➜ utiliser un endpoint interne si c’est critique.

Intégrité de dépendances :
➜ verrouiller requirements.txt, vérifier les hashes (voir “Intégrité” ci-dessous).

Disponibilité service externe :
➜ prévoir un fallback local ou une procédure utilisateur si l’endpoint est indisponible.

🛡️ Intégrité des données & du binaire
À la livraison (build/zip/exe)
version de MolarMassOCR,
date de build,
liste des versions clés (Python, RDKit, Pillow, gradio_client),
SHA256,
résultats de scans OK.

Côté application
Aucun écrasement de fichiers utilisateur.
Sorties affichées, pas persistées (sauf action explicite de l’utilisateur).
Calculs de masses/formules reproductibles (RDKit) à partir du SMILES retourné.
```
---
## ⚙️ Installation & lancement
```
Prérequis
Windows recommandé (tests GUI Tkinter).
Python 3.11+ (ou exécutable packagé).
Accès réseau sortant HTTPS vers *.hf.space et PubChem.

Démarrage
# Environnement Python classique
pip install rdkit-pypi pillow gradio_client requests
python app.py
Un écran Splash s’affiche pendant le chargement des bibliothèques (RDKit, PIL, gradio_client…).

🧭 Flux de traitement (simplifié)
Chargement image → création variante originale et variante améliorée.
Pour chaque variante :
a. Envoi à l’endpoint OCSR → SMILES
b. RDKit : MolFromSmiles, rendu, formule, masse exacte, masse moyenne
c. Plausibilité (valences/isotopes simples)
d. PubChem : tentative CAS, puis nom
Dé-doublonnage léger (masses proches) → navigation 1/…/N
Copie SMILES possible, logs lisibles dans la console interne.
```
🚧 Limitations connues

Dépendance réseau (OCSR + PubChem) : l’outil ne fonctionne pas hors-ligne tel quel.
Plausibilité basique (valences/isotopes) — ce n’est pas un validateur structural complet.
Qualité image : une image bruitée ou mal contrastée dégrade fortement l’OCSR.
CAS/nom : l’absence ou l’ambiguïté dans PubChem peut empêcher la résolution.


## ✅ Check-list IT / Sécurité
 Hash SHA256 du binaire/zip vérifié.
 Antivirus local (Defender / ClamAV) : OK.
 Pare-feu/Proxy : HTTPS → *.hf.space & pubchem.ncbi.nlm.nih.gov autorisés (ou endpoint interne configuré).
 Politique données sensibles : usage approuvé (ou bascule vers modèle interne).
 Requirements figés et scannés (SCA/Software Composition Analysis recommandé).
