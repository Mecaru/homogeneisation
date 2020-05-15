# Journal
----
## ToDo List
- Coder l'interface (notebook avec explications, interface user-friendly, "homogeneisation_main.ipynb")
    - Réfléchir à l'implémentation de la bdd
    - Réfléchir à l'implémentation des modèles inverses
    - Réfléchir à l'intégration de modèles auto-cohérents à l'interfaces
    - Cas matériaux poreux
    - Réfléchir à la représentation graphique de plus d'une inclusion
- Ajouter des modèles
    - Compléter le modèle de Mori-Tanaka
    - Modèles des bornes sup et inf
    - Réfléchir à l'intégration des modèles auto-cohérents avec phases
    - Réfléchir à l'importation des scripts de calcul des tenseurs d'Eshelby (contacter l'auteur pour l'implémentation sur JupyterLab)
---
## Suivi
### 27/04/2020
- Structure générale du code
- Initialisation des fichiers "journal.md", "classes.py"
- Fonctionnement des classes

### 29/04/2020
- Début de l'interface : génération d'une inclusion

### 30/04/2020
- Génération d'une microstructure
- Calcul de comportement homogénéisé avec modèles compatibles

### 04/05/2020
- Calcul automatisé depuis un fichier texte - 1ere version

### 05/05/2020
- Description des modèles avec équations

### 07/05/2020
- Ajout de la description des modèles en Markdown
- Méthode 'draw' permettant de dessiner les microstructures
- Début de l'implémentation de la définition de E et nu au lieu de K et G

### 12/05/2020 
- Comparaison de modèles entre eux
- Comparaison de modèles avec des données de fichiers .txt
- Ajout du calcul des bornes de Hashin lors du calcul du comportement homogénéisé
- Implémentation des variables E et nu en plus de K et G dans le cas isotrope
- Changement automatique du nom de l'inclusion par défaut
- Lecture des noms des modèles des fichiers input d'automatisation des simulations rendue insensible à la casse 