# TP3_IFT3295

Emanuel Rollin - 20106951
Github: https://github.com/LipoBruh/TP3_IFT3295

## Dépendances

Penser à installer les dépendances suivantes au besoin: 
```
import numpy as np
import os
```


## Utilisation générale : 

In CLI , enter string of format :

```python plast.py -i CGTAGTCGGCTAACCAGCATAACGCTTGTAAACGTAAGAGCCC -db tRNAs.fasta -E 5 -ss 1e3 -seed '11111111111'```

Va retourner un string ressemblant à ceci : 

_____RESULTS_____
  Seed used      : 11111111111
  Pattern match  : AGCGGGGTAGAGGAATTGGTTTACTCATCAGGCTCATGACCTGAAGACTGCAGGTTCGAATCCTGTCCCCGCC
  Sequence match : AGCGGGGTAGAGGAATTGGTCGACTCATCAGGCTCATGACCTGAAGACTGCAGGTTCGAATCCTGTCCCCGCC
  ID             : >M|cat|Carica_papaya
  Fusion score   : 347
  Bit score      : 503
  E value        : 8.87221492294373e-146

Copie de ce string est aussi envoyée dans le dossier "./output/"

## PARTIE 1 

### 1.1

Voir parser.py 

Sert à extraire les reads d'un fichier de format .fasta et les formatter dans une classe Read.

### 1.2

Voir tp3.py

find_all_HSP() sert à diviser le pattern en kmers et à créer des alignements par pattern matching exact. Un alignement simple est stocké dans la classe Alignement. Plusieurs alignements pour un même Kmer peuvent être stockés dans la classe alignements, tel que pour un même read et kmer on a plusieurs indexs possibles.

Dans plast.HSP, soit un dicitonnaire, on stocke les alignements ainsi : {kmer: Alignements}

NOTE : LE BONUS EST RESPECTÉ

Lorsque vous spécifiez le input string, toute combinaisons de "0" et "1" de taille >1 devrait fonctionner.

### 1.3

Voir boyer_moore.py

La recherche exacte et l'algorithme utilisé sont découplés dans boyer_moore.py pour clarté et réutilisabilité. Boyer Moore est légèrement modifié pour supporter un seed.

### 1.4

Voir tp3.py

L'extension gloutonne et sa logique est dans la méthode plast.extend(...). S'étend aussi loin que les bornes du read et de la séquence le permettent, tant que le score demeure au dessus du seuil (>4). Les extensions sont stockés dans une classe Extended_HSP qui contient les attributs pertinents.

On n'étend que les kmers filtrés respectant :
- meilleur score pour leur séquence respective
- existe au moins un autre kmer ayant un hsp dans cette séquence

### 1.5

Voir tp3.py

La fusion utilise une classe Fusion_HSP. C'est cette dernière qui gère la fusion. On initialise une fusion avec un premier Extended_HSP et on utilise fuse_HSP(Extended_HSP) pour tenter d'ajouter une deuxième séquence à la fusion. Comme les deux extensions utilisent les indexs de leur pattern et de leur read, on peut comparer ces indexs pour vérifier s'ils se croisent. S'ils sont disjoints, on ne fait pas la fusion. Sinon, on utilise les maximums et minimums de chaque paire d'indexs pour définir la fusion.

### 1.6

Voir tp3.py

Dans la classe Fusion_HSP, les méthodes Fusion_HSP.score(), Fusion_HSP.bitscore() Fusion_HSP.e_value(), Fusion_HSP.cut_off(ss) permettent calculer et retourner les attributs de la fusion. Fusion_HSP.cut_off(ss) retourne un booléen qui facilite le filtrage si le e_value ne respecte pas le cutoff.

### 1.7

Voir tp3.py

Dans la classe PLAST, run() utilise les attributs pour exécuter le processus complet. C'est aussi la méthode qui trie et filtre basé sur le e_value et le cutoff ss. Le output correspond aux attributs de la Fusion_HSP avec le e_value le plus petit.


## PARTIE 2

### 2.1

Dans plast.py, lorsque vous utilisez le CLI, vous utilisez PLAST.run(), ce qui vous permet de rechercher un seul pattern que vous fournissez et retourne le meilleur résultat.

Dans tp3.py, si vous exécutez le script, la question 2.1 sera exécutée et le résultat stocké dans le dossier "./output".

q2_seq1 et q2_seq4 semblent retourner un résultat valide

q2_seq2 et q2_seq3 ne peuvent pas retourner de résultat avec ce seed. Facile à confirmer avec ctrl+F dans le fichier tRNAs.fasta avec un court segment de 11 pb.

### 2.2

Ajouté au fichier d'output dans "./output" pour simplifier. On appelle simplement Fusion_HSP.to_reverse() et Fusion_HSP.to_AA().

La séquence d'ARN et d'AA correspondantes sont les deux dernières lignes dans le rapport.

### 2.3

Voir pdf Q2.3.pdf

### 2.4 (BONUS)
Dépendemment du format de la graine, on peut causer des manquements:
- Un seed du format "11111111111" matche exactement le kmer de départ, implique qu'on insiste qu'au moins 11 char matchent consécutivement
- Augmenter à un seed "111111111111111111111111" pourrait devenir trop rigide et retourner aucun match, car n'accomode pas pour les gaps et force la présence d'une longue séquence dans le read. Il est commun de voir au moins une insertion ou substitution dans un pattern recherché, et augmenter la rigidité de notre pattern ne permet pas d'accomoder ces phénomènes (peu sensible, / trop spécifique). La vitesse de recherche avec Boyer-Moore est généralement accélérée avec un long pattern.
- Diminuer à un seed "11111" est peu précis et va retourner des HSP dans tous les reads. Il est très facile de trouver 5 matchs consécutifs, car pour un alphabet de 4 lettres, cette séquence sera présente 1 fois sur 1024 (plus sensible, peu précis / spécifique). Aussi, diminuer la taille du seed / pattern augmente la fréquence de matchs de Boyer-Moore et diminue la vitesse de recherche, car défavorise l'utilisation efficace des heuristiques.

### 2.5 (BONUS)
In CLI , enter string of format :

```python plast.py -i CGTAGTCGGCTAACCAGCATAACGCTTGTAAACGTAAGAGCCC -db tRNAs.fasta -E 5 -ss 1e2 -seed '111010010100110111'```

L'impact est mesuré dans l'algorithme Boyer Moore. La complexité est affectée de la longueur du pattern. Ici, les zéros correspondent à des skips, n'augment pas vraiment le temps de traitement de façon significative. 111010010100110111 est de complexité similaire à 11111111111.
- Préprocessing plus long si pattern plus long
- Meilleures heuristiques de recherche (plus rapide pour un pattern plus long)
- Complexité souvent notée ainsi :
  - Best case : O(n / m)
  - AVG case : O(n)
  - Worst case : O(n + m)

- Nuisible seulement dans le worst case. Bénéfique dans le best case. Comme les zéros n'ont pas vraiment d'impact sur Boyer-Moore avec l'implémentation actuelle, devrait augmenter la vitesse de recherche.