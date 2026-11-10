#' Default parameters
#'
#' Returns the full list of parameters used by [run_heatmaps()]. Pass only the
#' values you want to change in the `params` argument of [run_heatmaps()]; every
#' other parameter keeps the default listed below.
#'
#' @details
#' **Views**
#' \describe{
#'   \item{`versions`}{Default `c("A", "B", "C")`. Views to draw, any combination:
#'     `"A"` genes x cells, `"B"` genes x groups, `"C"` modules x groups.}
#' }
#'
#' **Orthologs** (used only when the datasets are of different species and `map`
#' is not supplied to [run_heatmaps()]; if all datasets have the same species no
#' conversion is done)
#' \describe{
#'   \item{`ortho_method`}{Default `"gprofiler"`. Ortholog source: `"gprofiler"`,
#'     `"homologene"` or `"babelgene"`.}
#'   \item{`ortho_non121`}{Default `"keep_both_species"`. How to treat non one-to-one
#'     orthologs (argument `non121_strategy` of [orthogene::convert_orthologs()]).}
#'   \item{`human_to_mouse`}{Default `TRUE`. `TRUE` if `genesets` are human symbols,
#'     `FALSE` if they are mouse symbols.}
#'   \item{`map_file`}{Default `NULL`. Path of an `.rds` file where the gene table
#'     is saved (e.g. `"map_orthologues.rds"`); `NULL` = no saving.}
#' }
#'
#' **Missing genes**
#' \describe{
#'   \item{`keep_missing`}{Default `FALSE`: only the genes present in **all** datasets are
#'     drawn (the dropped genes are listed in the console). `TRUE`: keep every gene of
#'     `genesets`; in a dataset where a gene is absent (or has no ortholog) its row is
#'     entirely `na_color` (grey) and its bar is 0. The genes concerned are listed in the
#'     console. If you change `keep_missing`, do not reuse `out$map` from a previous call:
#'     let the table be recomputed.}
#' }
#'
#' **Transformation**
#' \describe{
#'   \item{`clip`}{Default `2`. Clipping of the per-gene z-score to `[-clip, clip]`;
#'     `NULL` = no clipping. It also sets the limits of the colour scale (2 if `NULL`).}
#'   \item{`min_n`}{Default `10`. Views B and C: minimum number of cells expressing a gene
#'     in a group; below it the tile is `na_color` (`NA`).}
#' }
#'
#' **View A (cells)**
#' \describe{
#'   \item{`n_cells`}{Default `200`. Number of cells sampled at random per group
#'     (all cells if the group is smaller).}
#'   \item{`seed`}{Default `1`. Random seed of the sampling.}
#' }
#'
#' **Layout and text**
#' \describe{
#'   \item{`flexible`}{Default `FALSE`. `TRUE`: the blocks stretch with the plot window
#'     (widths in relative units, `col_w` and the 60 mm of view A are ignored).
#'     `FALSE`: fixed widths in mm.}
#'   \item{`show_n`}{Default `TRUE`. Horizontal bar plot of the number of expressing
#'     cells (whole dataset) on the right of each block. Views A and B only.}
#'   \item{`title`}{Default `"Module expression"`. Title of the figure; `NULL` = no title.}
#'   \item{`title_version`}{Default `TRUE`. Append the view name to the title.}
#'   \item{`title_size`}{Default `14`. Font size of the title.}
#'   \item{`col_w`}{Default `8`. Width of one column in mm. Views B and C, when
#'     `flexible = FALSE`.}
#'   \item{`gene_size`}{Default `8`. Font size of the gene names (rows).}
#'   \item{`col_size`}{Default `9`. Font size of the column names (B, C) and of the
#'     group labels under the blocks (A).}
#' }
#'
#' **Colours and tiles**
#' \describe{
#'   \item{`colors`}{Default `c("navy", "white", "firebrick")`. Heatmap colours for
#'     low (`-clip`), middle (0) and high (`+clip`) values. Exactly 3 colours.}
#'   \item{`na_color`}{Default `"grey85"`. Colour of the `NA` tiles (masked by `min_n`,
#'     or genes with no variance).}
#'   \item{`bar_color`}{Default `"grey40"`. Colour of the bar plot of expressing cells.}
#'   \item{`module_tiles`}{Default `FALSE`. Add a coloured tile on the left side of the
#'     first block, one colour per module.}
#'   \item{`group_tiles`}{Default `FALSE`. Add coloured tiles on top of every block,
#'     one colour per group (the same group name has the same colour in all blocks).}
#'   \item{`module_colors`}{Default `NULL` (automatic). Named vector of colours for the
#'     modules, e.g. `c(Module1 = "grey30", Module2 = "#E69F00")`. Names missing from
#'     the vector get an automatic colour.}
#'   \item{`group_colors`}{Default `NULL` (automatic). Named vector of colours for the
#'     groups, same logic as `module_colors`.}
#'   \item{`tile_size`}{Default `3`. Thickness of the tiles in mm.}
#' }
#'
#' **PNG export**
#' \describe{
#'   \item{`save_png`}{Default `FALSE`. Save one PNG per view.}
#'   \item{`out_dir`}{Default `"."`. Output folder, created if needed.}
#'   \item{`file_prefix`}{Default `"heatmap"`. Files are named
#'     `<out_dir>/<file_prefix>_<view>.png`, e.g. `heatmap_A.png`.}
#'   \item{`png_width`, `png_height`}{Default `12` and `7`. Size in inches.}
#'   \item{`png_res`}{Default `300`. Resolution in dpi.}
#' }
#' @return A named list of parameters.
#' @export
default_params <- function() list(
  versions       = c("A", "B", "C"),  # A: gènes x cellules | B: gènes x groupes | C: modules x groupes
  # orthologues
  ortho_method   = "gprofiler",       # "gprofiler", "homologene" ou "babelgene"
  ortho_non121   = "keep_both_species",  # stratégie orthogene pour les orthologues non 1:1
  human_to_mouse = TRUE,
  map_file       = NULL,              # ex. "map_orthologues.rds" pour sauvegarder la table
  # gènes manquants
  keep_missing   = FALSE,             # TRUE : garder les gènes absents d'un dataset (ligne grise)
  # transformation
  clip           = 2,                 # écrêtage du z-score ; NULL = pas d'écrêtage
  min_n          = 10,                # B/C : nb min de cellules exprimant le gène dans un groupe
  # version A
  n_cells        = 200,               # cellules tirées par groupe
  seed           = 1,
  # mise en page et texte
  flexible       = FALSE,             # TRUE : les blocs s'étirent avec la fenêtre
  show_n         = TRUE,              # barplot "cells expr." à droite (A et B)
  title          = "Module expression",  # NULL = pas de titre
  title_version  = TRUE,
  title_size     = 14,
  col_w          = 8,                 # largeur (mm) d'une colonne en B/C (si flexible = FALSE)
  gene_size      = 8,                 # taille de police des noms de gènes (lignes)
  col_size       = 9,                 # taille de police des noms de colonnes / groupes
  # couleurs et tiles
  colors         = c("navy", "white", "firebrick"),  # bas, milieu (0), haut
  na_color       = "grey85",
  bar_color      = "grey40",
  module_tiles   = FALSE,             # tile coloré par module, sur le côté
  group_tiles    = FALSE,             # tile coloré par groupe, au-dessus
  module_colors  = NULL,              # vecteur nommé, ex. c(Module1 = "red"); NULL = auto
  group_colors   = NULL,              # vecteur nommé, ex. c(ctrl = "grey"); NULL = auto
  tile_size      = 3,                 # épaisseur des tiles (mm)
  # PNG
  save_png       = FALSE,
  out_dir        = ".",
  file_prefix    = "heatmap",
  png_width      = 12,
  png_height     = 7,
  png_res        = 300
)

###############################################################################
# Ce qui se calcule UNE fois : table de gènes (+ orthologues) et gènes communs
###############################################################################
#' @noRd
sym_of <- function(species) if (species == "human") "gene_human" else "gene_mouse"

#' Build the gene table shared by all datasets
#'
#' If the datasets are of different species, calls [Orthologue_Table()] to map
#' the gene symbols. If **all datasets have the same species**, no conversion is
#' done: the gene sets are used as given and must be symbols of that species.
#' By default only the genes present in **all** datasets are kept; with
#' `keep_missing = TRUE` (see [default_params()]) every gene is kept and the
#' missing ones are drawn as grey rows. The resulting table defines the row
#' order of every heatmap.
#'
#' @param genesets Named list of gene symbols, one vector per module. Human
#'   symbols when the datasets mix species (see `human_to_mouse` in
#'   [default_params()]); symbols of the common species otherwise.
#' @param datasets Named list of datasets, see [run_heatmaps()].
#' @param P Parameter list, see [default_params()].
#' @return A data frame with the column `gene_human` and/or `gene_mouse`, plus
#'   `module` and `label`.
#' @export
prepare_map <- function(genesets, datasets, P = default_params()) {
  species <- unique(vapply(datasets, function(d) d$species, character(1)))
  bad <- setdiff(species, c("human", "mouse"))
  if (length(bad)) stop("`species` must be \"human\" or \"mouse\", not: ", paste(bad, collapse = ", "))

  if (length(species) == 1) {
    # un seul type d'espèce : pas de conversion d'orthologues
    message("All datasets are ", species, ": no ortholog conversion, genesets used as given.")
    st  <- unique(stack(genesets))
    map <- data.frame(gene = as.character(st$values),
                      module = factor(as.character(st$ind), levels = names(genesets)),
                      stringsAsFactors = FALSE)
    names(map)[1] <- paste0("gene_", species)
    map <- map[order(map$module), , drop = FALSE]
    rownames(map) <- NULL
  } else {
    map <- Orthologue_Table(genesets, Method = P$ortho_method,
                            human_to_mouse = P$human_to_mouse,
                            non121_strategy = P$ortho_non121)

    if (isTRUE(P$keep_missing)) {
      # gènes sans ortholog : Orthologue_Table() les perd, on les remet (NA dans l'autre espèce)
      from <- if (P$human_to_mouse) "human" else "mouse"
      to   <- if (P$human_to_mouse) "mouse" else "human"
      st   <- unique(stack(genesets))
      lost <- st[!as.character(st$values) %in% map[[paste0("gene_", from)]], , drop = FALSE]
      if (nrow(lost)) {
        add <- data.frame(a = as.character(lost$values),
                          b = factor(as.character(lost$ind), levels = names(genesets)),
                          c = NA_character_, stringsAsFactors = FALSE)
        names(add) <- c(paste0("gene_", from), "module", paste0("gene_", to))
        map <- rbind(map, add[, names(map), drop = FALSE])
      }
      # ordre = ordre des genesets
      map <- map[order(map$module, match(map[[paste0("gene_", from)]], unlist(genesets))), ]
      rownames(map) <- NULL
    }
  }
  if (!is.null(P$map_file)) saveRDS(map, P$map_file)
  print(map)

  # présence de chaque gène dans chaque dataset (NA = pas d'ortholog = absent)
  present <- lapply(datasets, function(d) {
    g <- map[[sym_of(d$species)]]
    !is.na(g) & g %in% rownames(d$obj)
  })
  everywhere <- Reduce(`&`, present)

  if (isTRUE(P$keep_missing)) {
    # on garde tout ; les lignes absentes d'un dataset seront grises
    message("Genes kept but missing in at least one dataset (grey rows there):")
    print(map[!everywhere, ])
  } else {
    # on garde seulement les gènes présents dans TOUS les datasets
    message("Genes dropped (missing in at least one dataset):"); print(map[!everywhere, ])
    map <- map[everywhere, , drop = FALSE]
  }

  map$module <- factor(map$module, levels = names(genesets))
  lab_col <- if ("gene_human" %in% names(map)) "gene_human" else "gene_mouse"
  lab <- map[[lab_col]]
  other <- setdiff(intersect(c("gene_human", "gene_mouse"), names(map)), lab_col)
  if (anyNA(lab) && length(other)) lab[is.na(lab)] <- map[[other]][is.na(lab)]
  map$label <- make.unique(lab)               # étiquettes de lignes
  map                                         # définit l'ordre des lignes pour tous les datasets
}

###############################################################################
# Extraction
###############################################################################
# A : une colonne par cellule (tirage au hasard par groupe)
#' @noRd
extract_cells <- function(d, map, n_cells = 200, seed = 1) {
  g    <- map[[sym_of(d$species)]]
  ok   <- !is.na(g) & g %in% rownames(d$obj)      # gène présent dans ce dataset
  grp  <- d$obj@meta.data[[d$grp]]
  expr <- GetAssayData(d$obj, layer = "data")[unique(g[ok]), , drop = FALSE]

  # nb de cellules exprimant le gène sur TOUT le dataset (pas seulement le tirage) ; 0 si absent
  n_total <- numeric(length(g))
  n_total[ok] <- Matrix::rowSums(expr > 0)[g[ok]]

  set.seed(seed)
  cells <- unlist(lapply(split(colnames(d$obj), grp),
                         \(x) sample(x, min(n_cells, length(x)))))
  # gènes absents : ligne entière à NA (case grise)
  mat <- matrix(NA_real_, nrow = length(g), ncol = length(cells),
                dimnames = list(map$label, cells))
  if (any(ok)) mat[ok, ] <- as.matrix(expr[, cells, drop = FALSE])[g[ok], , drop = FALSE]

  g_cells <- as.character(grp[match(cells, colnames(d$obj))])
  list(mat = mat,
       group   = factor(g_cells, levels = unique(g_cells)),   # ordre d'apparition = ordre des blocs
       n_total = unname(n_total))
}

# B / C : une colonne par groupe (moyenne) + masque par groupe
#' @noRd
extract_avg <- function(d, map, min_n = 10) {
  g    <- map[[sym_of(d$species)]]
  ok   <- !is.na(g) & g %in% rownames(d$obj)      # gène présent dans ce dataset
  grp  <- d$obj@meta.data[[d$grp]]
  expr <- GetAssayData(d$obj, layer = "data")[unique(g[ok]), , drop = FALSE]
  idx  <- split(seq_along(grp), grp)

  # gènes absents : ligne entière à NA (case grise), 0 cellule
  avg    <- matrix(NA_real_, nrow = length(g), ncol = length(idx),
                   dimnames = list(map$label, names(idx)))
  n_expr <- matrix(0, nrow = length(g), ncol = length(idx),
                   dimnames = list(map$label, names(idx)))
  n_total <- numeric(length(g))

  if (any(ok)) {
    # moyenne en espace non-log (comme AverageExpression), puis log1p
    a <- sapply(idx, \(i) Matrix::rowMeans(expm1(expr[, i, drop = FALSE])))
    n <- sapply(idx, \(i) Matrix::rowSums(expr[, i, drop = FALSE] > 0))
    a <- log1p(a)
    avg[ok, ]    <- a[g[ok], , drop = FALSE]
    n_expr[ok, ] <- n[g[ok], , drop = FALSE]
    n_total[ok]  <- Matrix::rowSums(expr > 0)[g[ok]]
  }

  avg[n_expr < min_n] <- NA          # trop peu de cellules positives -> case grise

  list(mat = avg, group = NULL, n_expr = n_expr, n_total = n_total)
}

###############################################################################
# Transformations
###############################################################################
#' @noRd
zrow <- function(mat, clip = NULL) {
  z <- t(scale(t(mat)))              # scale ignore les NA
  z[is.nan(z)] <- NA                 # NA et non 0 : le masque reste visible
  if (!is.null(clip)) z <- pmax(pmin(z, clip), -clip)
  z
}

# C : score par module = moyenne des z-scores (remplaçable par UCell)
#' @noRd
module_score <- function(zmat, d, map) {
  g  <- map[[sym_of(d$species)]]
  ok <- !duplicated(paste(map$module, g))
  z  <- zmat[ok, , drop = FALSE]
  s  <- t(sapply(levels(map$module), \(m)
                 colMeans(z[map$module[ok] == m, , drop = FALSE], na.rm = TRUE)))
  s[is.nan(s)] <- NA
  s
}

# échelle de couleur : bas / milieu (0) / haut
#' @noRd
make_col_fun <- function(clip, colors = c("navy", "white", "firebrick")) {
  if (length(colors) != 3)
    stop("`colors` needs exactly 3 colours: low, middle (0), high. Got ", length(colors))
  lim <- if (is.null(clip)) 2 else clip
  colorRamp2(c(-lim, 0, lim), colors)
}

# couleurs nommées pour les tiles : automatiques, surchargées par `user` (vecteur nommé)
#' @noRd
resolve_colors <- function(values, user = NULL, palette = "Set 3") {
  values <- unique(as.character(values))
  out <- stats::setNames(grDevices::hcl.colors(max(length(values), 1), palette)[seq_along(values)],
                         values)
  if (!is.null(user)) {
    if (is.null(names(user)))
      stop("Custom tile colours must be a NAMED vector, e.g. c(group1 = \"red\", group2 = \"blue\")")
    hit <- intersect(names(user), values)
    out[hit] <- unname(user[hit])
    miss <- setdiff(values, names(user))
    if (length(miss))
      message("No colour given for: ", paste(miss, collapse = ", "), " (automatic colour used)")
  }
  out
}

###############################################################################
# Plot d'un dataset
###############################################################################
#' @noRd
plot_block <- function(mat, title, map, col_fun, P, first = TRUE, by_gene = TRUE,
                       col_group = NULL, n_total = NULL,
                       group_cols = NULL, module_cols = NULL) {
  # --- barplot du nb de cellules exprimant le gène (à droite)
  ra <- NULL
  if (by_gene && P$show_n && !is.null(n_total)) {
    args <- list(anno_barplot(n_total, which = "row", width = unit(15, "mm"),
                              gp = gpar(fill = P$bar_color, col = NA)),
                 annotation_label = "cells\nexpr.")
    names(args)[1] <- paste0("n_", title)          # nom unique par dataset
    ra <- do.call(rowAnnotation, args)
  }

  # --- tiles de modules (côté gauche, 1er bloc seulement)
  la <- NULL
  if (isTRUE(P$module_tiles) && first) {
    rows_mod <- if (by_gene) as.character(map$module) else rownames(mat)
    args <- list(anno_simple(rows_mod, col = module_cols[unique(rows_mod)],
                             width = unit(P$tile_size, "mm")),
                 show_annotation_name = FALSE, show_legend = FALSE)
    names(args)[1] <- "module_tile"
    la <- do.call(rowAnnotation, args)
  }

  # --- tiles de groupes (au-dessus) : un bloc coloré par groupe (A) ou une case par colonne (B, C)
  ta <- NULL
  if (isTRUE(P$group_tiles)) {
    if (!is.null(col_group)) {
      labs <- levels(droplevels(col_group))
      anno <- anno_block(gp = gpar(fill = group_cols[labs], col = NA),
                         height = unit(P$tile_size, "mm"))
    } else {
      anno <- anno_simple(colnames(mat), col = group_cols[colnames(mat)],
                          height = unit(P$tile_size, "mm"))
    }
    args <- list(anno, show_annotation_name = FALSE, show_legend = FALSE)
    names(args)[1] <- paste0("gtile_", title)
    ta <- do.call(HeatmapAnnotation, args)
  }

  # --- mode A : un label par groupe sous chaque bloc de cellules
  ba <- NULL
  if (!is.null(col_group)) {
    lab_gp <- gpar(fontsize = P$col_size)
    labs   <- levels(droplevels(col_group))
    args <- list(anno_block(gp = gpar(fill = NA, col = NA),      # sans cadre ni fond
                            labels = labs,
                            labels_gp = lab_gp,
                            labels_rot = 90,                      # écriture verticale
                            height = max_text_width(labs, gp = lab_gp) + unit(2, "mm")),
                 show_annotation_name = FALSE)
    names(args)[1] <- paste0("grp_", title)
    ba <- do.call(HeatmapAnnotation, args)
  }

  # --- largeur : fixe (mm) ou relative (flexible)
  w <- if (isTRUE(P$flexible)) {
    if (is.null(col_group)) unit(ncol(mat), "null") else unit(1, "null")
  } else {
    if (is.null(col_group)) unit(ncol(mat) * P$col_w, "mm") else unit(60, "mm")
  }

  Heatmap(mat, name = "z-score", col = col_fun, na_col = P$na_color,
          cluster_rows = FALSE, cluster_columns = FALSE,
          row_split    = if (by_gene) map$module else NULL,
          row_gap      = unit(4, "mm"), row_title_rot = 0,
          column_split = col_group, column_gap = unit(1, "mm"),
          show_row_names    = first,
          show_column_names = is.null(col_group),
          row_names_side = "left", row_names_gp = gpar(fontsize = P$gene_size),
          column_names_gp = gpar(fontsize = P$col_size),
          column_title = title, border = TRUE,
          left_annotation   = la,
          right_annotation  = ra,
          top_annotation    = ta,
          bottom_annotation = ba,
          width = w)
}

###############################################################################
# Assemblage d'une version (+ PNG optionnel)
###############################################################################
#' @noRd
build <- function(extract, transform, datasets, map, P, version_label,
                  by_gene = TRUE, use_group = FALSE, file = NULL, ...) {
  col_fun <- make_col_fun(P$clip, P$colors)
  res <- lapply(datasets, \(d) { e <- extract(d, map, ...); list(e = e, z = transform(e$mat, d)) })

  # couleurs de tiles communes à tous les blocs (même nom = même couleur)
  all_groups <- unique(unlist(lapply(res, \(r)
                                     if (use_group) levels(r$e$group) else colnames(r$z))))
  group_cols  <- resolve_colors(all_groups, P$group_colors, "Set 3")
  module_cols <- resolve_colors(levels(map$module), P$module_colors, "Dark 3")

  hts <- lapply(seq_along(res), \(i)
                plot_block(res[[i]]$z, names(res)[i], map = map, col_fun = col_fun, P = P,
                           first       = (i == 1),
                           by_gene     = by_gene,
                           col_group   = if (use_group) res[[i]]$e$group else NULL,
                           n_total     = res[[i]]$e$n_total,
                           group_cols  = group_cols,
                           module_cols = module_cols))
  ht <- Reduce(`+`, hts)

  main <- P$title
  if (!is.null(main) && P$title_version) main <- paste0(main, " - ", version_label)

  draw_it <- function()
    draw(ht, ht_gap = unit(5, "mm"),
         column_title = main,
         column_title_gp = gpar(fontsize = P$title_size, fontface = "bold"))

  draw_it()
  if (!is.null(file)) {
    png(file, width = P$png_width, height = P$png_height, units = "in", res = P$png_res)
    draw_it()
    dev.off()
    message("PNG saved: ", file)
  }
  invisible(res)
}

###############################################################################
# Point d'entrée
#   datasets : liste(nom = list(obj, species, grp))
#   genesets : liste(module = c(symboles))
#   params   : liste partielle, complétée par default_params()
#   map      : table déjà calculée (sinon calculée ici) -> évite de refaire gprofiler
# Renvoie list(map, res) en invisible
###############################################################################
#' Cross-species heatmaps of gene modules
#'
#' Draws side-by-side heatmaps of the same gene modules across several Seurat
#' objects (human and/or mouse). When the datasets mix species, mouse symbols
#' are obtained through an ortholog table computed once with
#' [Orthologue_Table()]. When all datasets have the same species, no conversion
#' is done and the gene sets are used as given.
#'
#' Three views are available:
#' \itemize{
#'   \item **A** genes x cells (random sample of cells per group),
#'   \item **B** genes x groups (average expression, with a mask for groups where
#'     fewer than `min_n` cells express the gene),
#'   \item **C** modules x groups (mean z-score per module).
#' }
#'
#' @param datasets Named list. Each element is a list with `obj` (Seurat object,
#'   log-normalised data in the `data` layer), `species` (`"human"` or `"mouse"`)
#'   and `grp` (name of the metadata column that defines the groups).
#'   The first dataset is the reference for row labels. A single dataset is allowed.
#' @param genesets Named list of gene symbols, one vector per module. **Human**
#'   symbols when the datasets mix species; symbols of the common species when all
#'   datasets have the same species (no conversion in that case).
#' @param params Named list of parameters to override, see [default_params()].
#' @param map Gene table returned by a previous call (`out$map`). Avoids
#'   recomputing the orthologs.
#' @return Invisibly, a list with `map` (gene table) and `res`
#'   (one element per view, with the matrices that were plotted).
#' @examples
#' \dontrun{
#' # One list element per dataset: Seurat object, species, grouping column.
#' # The list names ("Dataset1_human", ...) are used as titles of the blocks.
#' datasets <- list(
#'   Dataset1_human = list(obj = seurat_Obj1, species = "human", grp = "orig.ident"),
#'   Dataset2_mouse = list(obj = seurat_Obj2, species = "mouse", grp = "orig.ident"),
#'   Dataset3_mouse = list(obj = seurat_Obj3, species = "mouse", grp = "orig.ident")
#' )
#'
#' # Gene sets: HUMAN symbols, one vector per module (replace with your own genes)
#' genesets <- list(
#'   Module1_housekeeping   = c("GAPDH", "ACTB", "B2M", "PPIA", "RPLP0", "TBP"),
#'   Module2_T_cell        = c("CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD4"),
#'   Module3_myeloid       = c("CD14", "LYZ", "CD68", "CSF1R", "ITGAM", "LST1")
#' )
#'
#' # Minimal call: all parameters at their default values
#' out <- run_heatmaps(datasets, genesets)
#'
#' # Change a few parameters
#' out <- run_heatmaps(datasets, genesets,
#'                     params = list(versions = c("A", "C"), clip = 3, gene_size = 6,
#'                                   flexible = TRUE, colors = c("blue", "white", "red"),
#'                                   module_tiles = TRUE, group_tiles = TRUE))
#'
#' # Re-run with other settings without recomputing the orthologs
#' out2 <- run_heatmaps(datasets, genesets,
#'                      params = list(versions = "B", save_png = TRUE),
#'                      map = out$map)
#'
#' # Keep ALL genes, even those missing in a dataset (grey rows there).
#' # No `map = out$map` here: the gene table must be recomputed.
#' out3 <- run_heatmaps(datasets, genesets,
#'                      params = list(keep_missing = TRUE, min_n = 0))
#'
#' # One species only (here a single mouse dataset): no ortholog conversion,
#' # genesets are MOUSE symbols
#' run_heatmaps(list(Dataset2_mouse = datasets$Dataset2_mouse),
#'              list(Module_T_cell = c("Cd3d", "Cd3e", "Cd2", "Cd8a")))
#'
#' # Every parameter, explained: see ?default_params, the README, or
#' # system.file("examples", "run_heatmaps_example.R", package = "CrossHeat")
#' }
#' @export
run_heatmaps <- function(datasets, genesets, params = list(), map = NULL) {
  P <- modifyList(default_params(), params)
  unknown <- setdiff(names(params), names(default_params()))
  if (length(unknown)) warning("Unknown parameter(s) ignored: ", paste(unknown, collapse = ", "))
  if (is.null(map)) map <- prepare_map(genesets, datasets, P)
  if (isTRUE(P$save_png)) dir.create(P$out_dir, showWarnings = FALSE, recursive = TRUE)

  png_file <- function(v) if (isTRUE(P$save_png))
    file.path(P$out_dir, paste0(P$file_prefix, "_", v, ".png")) else NULL

  one <- function(v) switch(v,
    A = build(extract_cells, \(m, d) zrow(m, P$clip), datasets, map, P,
              "A: genes x cells", use_group = TRUE, file = png_file("A"),
              n_cells = P$n_cells, seed = P$seed),
    B = build(extract_avg, \(m, d) zrow(m, P$clip), datasets, map, P,
              "B: genes x groups", file = png_file("B"),
              min_n = P$min_n),
    C = build(extract_avg, \(m, d) module_score(zrow(m, P$clip), d, map), datasets, map, P,
              "C: modules x groups", by_gene = FALSE, file = png_file("C"),
              min_n = P$min_n),
    stop("Unknown version: ", v, " (choose among A, B, C)"))

  res <- setNames(lapply(P$versions, one), P$versions)
  invisible(list(map = map, res = res))
}
