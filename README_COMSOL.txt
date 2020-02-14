https://indico.cern.ch/event/532518/contributions/2200213/attachments/1288251/1917249/COMSOL_Garfield_integration-3.pdf


1) New component
Creer un composant
Options possibles : booleans / subtraction pour soustraire un matériau
ou transformations -> réseau pour reproduire plusieurs fois un motifs (des trous par exemple)

2) New materials

3) Ajouter un maillage
Attention ! Choisir la resolution de la mesh avec attention, ca doit être le plus fin possible

4) Ajouter une physique (AC/DC électrostatique)
- une pour les vraies configurations de champs électriques
- une pour chaque electrode de lecture

5) Ajouter une etude (stationnaire)

6) Calculer

7) Exporter les donéees
- données avec V (potentiel de l'etude), V2, V3, ... (avec Vi weighting potential de l'electrode i, V=1 a l'electrode et 0 sur toutes les autres)
Attention, dans options avancées, il faut aller dans resolution --> nombre de noeuds = 2
- maillage mesh.mphtxt (!! attention a l'extension)

