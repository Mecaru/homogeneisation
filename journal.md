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
- Organiser les dossier d'entrées et sorties (un dossier input, puis un dossier par type d'input, idem pour les outputs)
- Ajouter les graphes de validation des modèles aux descriptions de modèles
- Comportements anisotropes
- Patterns isotropes
    
### Programme de la semaine
- Intégrer les modèles autocohérents
- Suite du visco

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
- Structure théorique de l'interface de la section modèle inverse et structure des données
- Mise à jour du critère d'arrêt pour la prise en compte de la variation de l'erreur
- Exportation d'un graphe contenant uniquement des tracés de modèles en fichier texte

### 26/05/2020
- Exportation d'un graphe avec suppression automatique des données provenant de fichiers d'entrée
- Suppression des grid, axe des fractions volumiques limité à 0 - 1 (section comparaison de modèles)
- Interface modèles inverses : choix du comportement cible et début de la création d'inclusions

### 27/05/2020
- Interface modèles inverses: Création d'inclusions, de microstructure et début de la génération du dictionnaire des degrés de liberté de l'algorithme d'optimisation
- Interface modèles inverses: Génération de microstructures avec paramètres de comportement de matrice, d'inclusions ou fractions volumiques d'inclusions inconnus
- Interface modèle inverse fonctionnelle, non optimisée pour le cas de variables avec des range très différents

### 28/05/2020
- Génération d'inclusions visco-élastiques

### 29/05/2020
- Génération de microstructures visco-élastiques

### 02/06/2020
- Visco-élasticité: Interpolation des listes de fréquences
- Visco-élasticité: Début de la modification des modèles

### 03/06/2020
- Fusion de cellules
- Passage des bornes max des paramètres de comportement à 10^6
- Adaptation de Mori-Tanaka aux comportement visco-élastiques
- Adaptation d'Eshelby au visco-élastique et début de la modification du différentiel

### 08/06/2020
- Ajout du modèle différentiel visco-élastique sans tracé de la partie imaginaire de K
- Ajout du tracé de K' et K''
- Manipulation des graphes matplotlib sur le main visco

### 09/06/2020
- calculs depuis un fichier .txt: Modification du format des fichiers de sortie, section valable uniquement pour des calculs en non-visco isotrope 
- Possibilité d'entrer un nu=0.5
- Visco-élasticité: implémentation du choix de l'abscisse (fréquence ou temperature)