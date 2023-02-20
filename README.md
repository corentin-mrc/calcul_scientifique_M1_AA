# Calcul_scientifique_M1_AA

Projet de calcul scientifique en M1 algèbre appliquée.
Réalisé par Walid ABED et Corentin MARCOU.

## Organisation du code:

Nous avons utilisé la notion de classe. Il y a trois classes:
1. Noeud: contient un point (x, y) tel que x est dans [-a, a] et y dans [-b, b].
2. Triangle: contient un tableau de trois noeud, contient les méthodes associées à FT et à la matrice BT.
3. Maillage: contient un tableau de triangles, contient les méthodes assocées aux indices i, j, aux numéro global et aux numéro intérieur.
En plus de ces classes, le fichier boite_a_outils contient toutes les méthodes relatives à l'approximation de la solution.

Pour compiler le projet, on pourra utilisé le fichier CMakeLists.txt avec CMake ou utilisé directement g++.
L'affichage des données obtenues se fait avec Python.
