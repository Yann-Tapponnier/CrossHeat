##############################################################################################################################
#################################variante de CLAUDE pour conservé les gene en module #########################################
##############################################################################################################################
# Persmet et garder les Orthologue des 2 espèces = cool pour le heatmap:
# Script indépendant : source("Orthologue_Table.R") au début de la session.
# Nécessite le package orthogene (BiocManager::install("orthogene")).
Orthologue_Table <- function(gene_lists,
                             Method = c("gprofiler", "homologene", "babelgene"),
                             human_to_mouse = TRUE,
                             non121_strategy = "keep_both_species") {
  Method <- match.arg(Method)
  from <- if (human_to_mouse) "human" else "mouse"
  to   <- if (human_to_mouse) "mouse" else "human"
  
  genes_df <- unique(stack(gene_lists))
  names(genes_df) <- c("gene_in", "module")
  
  res <- orthogene::convert_orthologs(
    gene_df         = unique(genes_df$gene_in),
    input_species   = from,
    output_species  = to,
    method          = Method,
    non121_strategy = non121_strategy,
    gene_output     = "columns"
  )
  res <- res[, c("input_gene", "ortholog_gene")]
  
  out <- merge(genes_df, res, by.x = "gene_in", by.y = "input_gene")
  out$module <- factor(out$module, levels = names(gene_lists))
  out <- out[order(out$module), ]
  names(out)[names(out) == "gene_in"]       <- paste0("gene_", from)
  names(out)[names(out) == "ortholog_gene"] <- paste0("gene_", to)
  rownames(out) <- NULL
  out
}
