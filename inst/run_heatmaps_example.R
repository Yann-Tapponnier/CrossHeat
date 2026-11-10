###############################################################################
# CrossHeat : exemple complet, TOUS les paramètres sont listés et expliqués.
# Les valeurs écrites ci-dessous sont les valeurs par défaut.
# Ne garder dans `params` que ce que tu veux changer (le reste = défaut).
# Aide : ?run_heatmaps  |  ?default_params
###############################################################################
library(CrossHeat)

# =============================================================================
# 1. DATASETS
#    Une entrée par jeu de données. Le NOM de l'entrée sert de titre au bloc.
#      obj     : objet Seurat (v5), données log-normalisées dans la couche "data"
#      species : "human" ou "mouse"
#      grp     : colonne de meta.data qui définit les groupes (colonnes/blocs)
#    Le 1er dataset est la référence pour les étiquettes de lignes.
# =============================================================================
datasets <- list(
  Dataset1_human = list(obj = seurat_Obj1, species = "human", grp = "orig.ident"),
  Dataset2_mouse = list(obj = seurat_Obj2, species = "mouse", grp = "orig.ident"),
  Dataset3_mouse = list(obj = seurat_Obj3, species = "mouse", grp = "orig.ident")
)

# =============================================================================
# 2. GENESETS
#    Symboles HUMAINS (si genesets en souris : human_to_mouse = FALSE).
#    Un vecteur par module ; l'ordre de la liste = l'ordre des blocs de lignes.
#    (Listes génériques : à remplacer par tes gènes.)
# =============================================================================
genesets <- list(
  Module1_housekeeping   = c("GAPDH", "ACTB", "B2M", "PPIA", "RPLP0", "TBP"),
  Module2_T_cell        = c("CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD4", "IL7R"),
  Module3_myeloid       = c("CD14", "LYZ", "CD68", "CSF1R", "ITGAM", "LST1", "FCGR3A")
)

# =============================================================================
# 3. PARAMETRES (tous)
# =============================================================================
params <- list(

  # ---- Vues -----------------------------------------------------------------
  # "A" = gènes x cellules | "B" = gènes x groupes | "C" = modules x groupes
  # Toute combinaison : "A", "B", c("A","C"), c("A","B","C")
  versions       = c("A", "B", "C"),   # + "D" : dot plot (couleur = moyenne z-scorée, taille = % de cellules)

  # ---- Orthologues (utilisés seulement si `map` n'est pas fourni) ------------
  # (si TOUS les datasets ont la même espèce : aucune conversion, genesets = symboles de cette espèce)
  ortho_method   = "gprofiler",         # "gprofiler", "homologene" ou "babelgene"
  ortho_non121   = "keep_both_species", # orthologues non 1:1 (argument non121_strategy d'orthogene)
  human_to_mouse = TRUE,                # TRUE : genesets en symboles humains ; FALSE : en symboles souris
  map_file       = NULL,                # ex. "map_orthologues.rds" pour sauvegarder la table ; NULL = non

  # ---- Gènes manquants --------------------------------------------------------
  keep_missing   = FALSE,               # FALSE : seuls les gènes présents dans TOUS les datasets sont tracés
                                        # TRUE  : tous les gènes sont gardés ; si un gène manque dans un
                                        #   dataset (ou n'a pas d'ortholog), sa ligne est toute grise (na_color)
                                        #   et son barplot vaut 0. Les gènes concernés sont listés dans la console.
                                        #   Si tu changes keep_missing, ne réutilise PAS out$map.

  # ---- Transformation ---------------------------------------------------------
  clip           = 2,                   # écrêtage du z-score à [-clip, clip] ; NULL = aucun
                                        #   (fixe aussi les bornes de l'échelle de couleur)
  min_n          = 10,                  # vues B et C : nb min de cellules exprimant le gène dans un
                                        #   groupe, sinon case grise (NA)

  # ---- Vue A (cellules) -------------------------------------------------------
  n_cells        = 200,                 # nb de cellules tirées au hasard par groupe
  seed           = 1,                   # graine du tirage

  # ---- Mise en page et texte ---------------------------------------------------
  flexible       = FALSE,               # TRUE : les blocs s'étirent avec la fenêtre (largeurs relatives ;
                                        #   col_w et les 60 mm de la vue A sont alors ignorés)
  show_n         = TRUE,                # barplot horizontal "cells expr." à droite (vues A et B)
  title          = "Module expression", # titre de la figure ; NULL = pas de titre
  title_version  = TRUE,                # ajoute le nom de la vue au titre
  title_size     = 14,                  # taille de police du titre
  col_w          = 8,                   # largeur d'une colonne en mm (vues B et C, si flexible = FALSE)
  gene_size      = 8,                   # taille de police des noms de gènes (lignes)
  col_size       = 9,                   # taille de police des noms de colonnes (B, C)
                                        #   et des labels de groupes (A)

  # ---- Couleurs et tiles ------------------------------------------------------
  colors         = c("navy", "white", "firebrick"),  # échelle de -clip à +clip : 2 couleurs ou plus, réparties
                                        #   régulièrement (nb impair : couleur du milieu = 0). Ou un seul mot :
                                        #   "pheatmap" (défaut de pheatmap = RdYlBu inversée, 7 couleurs),
                                        #   ou un nom RColorBrewer ("RdBu", "rev_Spectral" = inversée)
  color_range    = NULL,                # NULL : échelle de -clip à +clip ; ou c(min, max), ex. c(-1, 3)
                                        #   (échelle asymétrique ; garder clip >= max)
  na_color       = "grey85",            # couleur des cases NA (masquées par min_n, ou gène sans variance)
  bar_color      = "grey40",            # couleur du barplot "cells expr."
  module_tiles   = FALSE,               # TRUE : tile coloré par module, sur le côté gauche du 1er bloc
  group_tiles    = FALSE,               # TRUE : tiles colorés par groupe, au-dessus de chaque bloc
  module_colors  = NULL,                # vecteur NOMMÉ, ex. c(Module1_housekeeping = "grey30"); NULL = auto
  group_colors   = NULL,                # vecteur NOMMÉ, ex. c(ctrl = "grey60", treated = "firebrick"); NULL = auto
                                        #   (même nom de groupe = même couleur dans tous les blocs ;
                                        #    noms manquants = couleur automatique)
  tile_size      = 3,                   # épaisseur des tiles (mm)

  # ---- Vue D (dot plot) ---------------------------------------------------------
  # couleur = moyenne par groupe (z-score par gène) ; taille du point = proportion de cellules
  # du groupe qui expriment le gène. min_n = 0 pour dessiner tous les points.
  row_h          = 8,                   # hauteur d'une ligne en mm (si flexible = FALSE)
  dot_scale      = 0.9,                 # diamètre max du point (100 % des cellules), fraction de la case

  # ---- Export PNG -----------------------------------------------------------
  save_png       = FALSE,               # TRUE : un PNG par vue
  out_dir        = ".",                 # dossier de sortie (créé si besoin)
  file_prefix    = "heatmap",           # fichiers : <out_dir>/<file_prefix>_<vue>.png
  png_width      = 12,                  # largeur en pouces
  png_height     = 7,                   # hauteur en pouces (augmenter si beaucoup de gènes)
  png_res        = 300                  # résolution en dpi
)

# =============================================================================
# 4. LANCEMENT
#    out$map : table d'orthologues (gene_human, gene_mouse, module, label)
#    out$res : matrices tracées, une entrée par vue
# =============================================================================
out <- run_heatmaps(datasets, genesets, params)

# Relancer avec d'autres réglages SANS refaire g:Profiler (argument `map`) :
out2 <- run_heatmaps(datasets, genesets,
                     modifyList(params, list(versions = "B", clip = 3, gene_size = 6,
                                             save_png = TRUE, out_dir = "figures",
                                             file_prefix = "CrossHeat")),
                     map = out$map)
# Si tu changes `genesets` ou `datasets`, ne pas réutiliser out$map.

# Dot plot (vue D) : palette "dotplot" sur l'échelle -1 -> 3
out_dot <- run_heatmaps(datasets, genesets,
                        modifyList(params, list(versions = "D", colors = "dotplot",
                                                color_range = c(-1, 3), clip = 3, min_n = 0)),
                        map = out$map)

# Garder TOUS les gènes (lignes grises là où un gène manque) ; min_n = 0 pour ne masquer aucune case.
# Pas de `map = out$map` ici : la table de gènes doit être recalculée.
out_all <- run_heatmaps(datasets, genesets,
                        modifyList(params, list(keep_missing = TRUE, min_n = 0)))

# Blocs qui s'étirent + couleurs perso + tiles :
out3 <- run_heatmaps(datasets, genesets,
                     modifyList(params, list(flexible = TRUE,
                                             colors = "pheatmap",
                                             module_tiles = TRUE, group_tiles = TRUE)),
                     map = out$map)

# Une seule espèce (ici un seul dataset souris) : pas de conversion d'orthologues,
# les genesets sont alors des symboles SOURIS.
out4 <- run_heatmaps(list(Dataset2_mouse = datasets$Dataset2_mouse),
                     list(Module2_T_cell = c("Cd3d", "Cd3e", "Cd3g", "Cd2")),
                     list(versions = c("A", "B")))
