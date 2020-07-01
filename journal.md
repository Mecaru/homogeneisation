# Journal
----
## ToDo List
- Coder l'interface (notebook avec explications, interface user-friendly, "homogeneisation_main.ipynb")
    - Réfléchir à l'implémentation de la bdd
    - Réfléchir à l'intégration de modèles auto-cohérents à l'interface
- Ajouter des modèles
    - Réfléchir à l'intégration des modèles auto-cohérents avec phases
- Distribution des ellipsoïdes
- Calcul des tenseurs d'Eshelby
- Validation des modèles
- Corriger les modèles
    - Revoir le calcul de K de Mori-Tanaka
    - Eshelby en dehors des bornes
- Ajouter les graphes de validation des modèles aux descriptions de modèles
- Patterns isotropes
- 4 phases avec différents patterns
- Améliorer la résolution inverse (interface et soft)
- Résolution inverse en visco et avec une interphase
    
### Programme de la semaine
- Reprendre la résolution inverse avec classes_v2
- Valider différentiel
- Organiser les dossier d'entrées et sorties (un dossier input, puis un dossier par type d'input, idem pour les outputs)
- Reprendre calculs automatisés
- Notice modèle et changement behavior acceptable
- Présentation

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
- Intégration des modèles auto-cohérents au main

### 12/06/2020
- Fusion de classes.py et classes_visco.py en classes_v2.py, ajout des modèles auto-cohérents en visco, ajout de nouveaux modèles facilité

### 16/06/2020
- Fusion des cellules "calcul de comportement homogénéisé" sur les deux mains
- Visco: Ajout de modèles sur le même graphe avec tracé des différents paramètres sur différents subplots

### 17/06/2020
- Visco: changement des limites de l'axe x du graphe
- Visco: sauvegarde des graphes et ajout d'un dossier visco dans outputs
- Visco: sauvegarde des graphes en csv
- Visco: ajout des comportements en K' et K'', ou E' et E''
- Ajout de comportements anisotropes via fichiers texte sous la forme de matrices 6X6
- Passage des méthodes str et repr des classes modèles à la classe mère Model

### 18/06/2020
- Suppression des prolate et oblate, ajout d'inclusions ellipsoïdales avec deux rapports d'aspect et leur affichage
- Modification des arguments d'entrée des fonctions compute avec prise en compte des rapports d'aspect

### 19/06/2020
- Correction de bugs et version fonctionnelle des comportements anisotropes

### 22/06/2020
- Création du dossier anisotropic_behaviors
- Ajout d'une instance InclusionAndInterphase et modifications des classes

### 23/06/2020
- Ajout de la génération d'inclusions avec interphases au main
- Ajout de la génération d'inclusions avec interphases au main visco
- Dessin de plusieurs inclusions
- Dessin d'inclusions avec interphases

### 24/06/2020
- Comparaison de modèles contenant des inclusions avec interphases
- Reset de la microstructure à son état initial après la comparaison de modèles
- Correction d'un bug dû à la fonction draw lors de la génération de microstructures sans inclusions + correction d'une erreur lors du tracé des bornes de Hashin lorsque leur calcul n'est pas possible
- Ajout du couple de comportement isotrope E visco-élastique et K élastique
- Correction d'une erreur sur le nom des matrices C et S (stifness et compliance)
- Réduction de la taille des graphes et de la taille de police

### 25/06/2020
- Affichage des infos des inclusions et suppression d'inclusions, affichage des noms des inclusions plutôt que leur description dans les Dropdown (main)

### 26/06/2020
- Affichage des noms des inclusions plutôt que leur description dans la section comparaison de modèles du main
- Affichage des noms des inclusions plutôt que leur description et ajout des infos inclusions dans le main visco

### 29/06/2020
- Résolution inverse pour des inclusions sans interphase en élasticité isotrope
- Comparaison des modèles autocohérents old et new

---

### 30/06/2020
- Remplacement des tuples de fraction volumique des inclusions+interphase par des listes mutables, simplification du code associé dans comparaison de modèles
- Résolution inverse avec inclusions et interpahses

### 01/07/2020
- Fusion des cellules résolution inverse et correction d'un bug mineur