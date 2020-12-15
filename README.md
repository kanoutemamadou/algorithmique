# Projet 6 algorithmique: Alignement de séquences d'ADN
Etudiants:  
Mamadou KANOUTE  
Cédric Djiomou  
Minh-Quan Tran  

## Quick start

Le package `alignSeq` est un package permettant d'obtenir le meilleur alignement possible entre deux séquences d'ADN en s'inspirant de l'algorithme de **Needleman-Wunsch**.   
**Ce meilleur alignement est obtenu en prenant l'alignement qui donne le meilleur score**. 

L'algorithme utilisé **ne permet d'affirmer l'unicité de la solution** dû au fait qu'on peut avoir plusieurs alignements différents avec le même score.

L'algorithme de **Needleman-Wunsch** permet de passer d'une complexité exponentielle à une complexité **O(np)** avec **n** et **p** les longueurs des deux séquences, **contrairement à la version naïve de complexité exponentielle où on enumère tous les alignements possibles**.


## Installation du package

Pour une installation facile sous R, vous avez besoin du package devtools. Ensuite excécuter les 2 lignes de codes dans la console de R.

**devtools::install_github("kanoutemamadou/algorithmique")**   
**library(alignSeq)**


### Exemple inspiré de l'article sur Wikipédia
               
`A = GCATGCU`     
     `B =  GATTACA`    
     `d = -1`     
     `match = 1`      
     `dismatch = -1`       
     
Notons **p** la taille de la séquence A
**n** la taille de la séquence B.
Différentes étapes:
On crée une matrice de similarité S entre nos deux séquences.

`S_ij = 1 si A[i] = B[j]`           

`S_ij = -1 si A[i] != B[j]`         

`S_ij = -1 si on a un decalage, par exemple S[i] = - et S[j] = G, vice versa `                

Après on crée une matrice de taille (n+1)x(p+1) dont ses colonnes sont les éléménts de A, ses lignes les éléments de B, 

On l'appelle souvent `matrice F`.

#### Initialisation:
Les premières colonne et ligne de F sont remplies en commençant par 0, ensuite on enlève d.

F[1,1] = 0        
F[1,2] = F[1,1] + d ........ F[1,p] = F[1,p-1] + d      
F[2,1] = F[1,1] + d ....... F[n,1] = F[n-1,1] + d   


![alt text](Initialisation.png)

A partir de là, on remplit les colonnes au fur et à mesure en appliquant les formules suivantes:
```javascript
  top = F[i-1, j] + d
  
  left = F[i, j-1] + d
  
  topLeft = F[i-1, j-1] + S[A[i], B[j]]
  
  F[i,j] = max(left, top, topLeft)
```

![alt text](F_construite.png)

  **Une fois, F entièrement construite, on cherche quel alignement donne le score maximal**.     
  Pour cela, on cherche le score maximal donné.On voit bien ici que c'est `F[n,p]`. 
  
  On part de cette position et on remonte vers la position (1,1), en regardant à chaque étape à partir de quel voisin on est parti.
  
  S'il s'agit de l'élément:
  
  - diagonal, `A[i] et B[j] sont alignés`.
  
  - à la position `(i-1,j)` alors `A[i]` est aligné avec un trou, c'est à dire le symbole   -.
  
  - à la position `(i, j-1)` alors `B[j]` est aligné avec un trou, c'est à dire le symbole  -.
  
  
![alt text](Needleman-Wunsch_pairwise_sequence_alignment.png)


