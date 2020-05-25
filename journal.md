# Journal
----
## ToDo List
- Coder l'interface (notebook avec explications, interface user-friendly, "homogeneisation_main.ipynb")
    - Réfléchir à l'implémentation de la bdd
    - Réfléchir à l'intégration de modèles auto-cohérents à l'interface
    - Cas matériaux poreux
    - Réfléchir à la représentation graphique de plus d'une inclusion
- Ajouter des modèles
    - Réfléchir à l'intégration des modèles auto-cohérents avec phases
- Ajouter tous les modèles simples
- Distribution des ellipsoïdes
- Calcul des tenseurs d'Eshelby
- Validation des modèles
- Comportements non isotropes avec fichier texte
- Corriger les modèles
    - Revoir le calcul de K de Mori-Tanaka
    - Eshelby en dehors des bornes
### Programme de la semaine
- Terminer l'interface des modèles inverses
    - Optimisation de la norme du vecteur de descente (pas)
    - Changement du critère d'arrêt - critère sur la variation de l'erreur plutôt que sur sa valeur
    - Interface
- Exporter les graphes en .txt
- Enlever la grid sur les plots
- fixer l'axe f entre 0 et 1
- Comportements anisotropes
- Comportements visco-élastiques "T" (ou "fréquence"), "K'", "K''", "G'", "G''" et K constant
- Patterns isotropes

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

### 13/05/2020
- Ajout d'inclusions ellipsoïdales et choix du rapport de forme

### 15/05/2020
- Dessin des ellipsoïdes lors de la génération d'une microstructure avec une seule inclusion
- Ajout des modèles d'Eshelby et différentiel au master + Test 
- Tracé de plusieurs paramètres différents lors de la comparaison de modèles

### 18/05/2020
- Modification de la structure de la comparaison de modèles
- Début de la sauvegarde des figures

### 19/05/2020
- Sauvegarde des figures dans le dossier outputs
- Ajout de données depuis un fichier texte avec plusieurs paramètres
- Idée de la structure de données pour la section modèles inverses

### 25/05/2020
- Version fonctionnelle de l'algorithme de résolution inverse à plusieurs paramètres
- Optimisation du temps de calcul - suppression des sauvegardes de fichiers