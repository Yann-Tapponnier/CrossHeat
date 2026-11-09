#' Ortholog table for gene modules
#'
#' Converts a list of gene sets between human and mouse with
#' [orthogene::convert_orthologs()] and returns one row per gene.
#'
#' @param gene_lists Named list of gene symbols, one vector per module.
#' @param Method Ortholog source: `"gprofiler"`, `"homologene"` or `"babelgene"`.
#' @param human_to_mouse `TRUE` if `gene_lists` are human symbols (default),
#'   `FALSE` if they are mouse symbols.
#' @param non121_strategy How to handle non one-to-one orthologs, passed to
#'   [orthogene::convert_orthologs()].
#' @return A data frame with columns `gene_human`, `gene_mouse` and `module`
#'   (a factor ordered like `gene_lists`).
#' @export
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
