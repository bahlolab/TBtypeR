
make_tbprofiler_bc <- function(
    output = 'TBtypeR.TBprofiler_barcode.bed',
    panel = TBtypeR::tbt_panel
)
{
  panel %>%
    rename_panel_tbprofiler() %>%
    transmute(
      chrom = 'Chromosome',
      start = pos - 1L,
      end = pos,
      lineage = phylotype,
      allele = if_else(genotype == 0L, ref, alt),
      X1 = str_extract(phylotype, '[^\\.]+'),
      X2 = 'None',
      X3 = 'None'
    ) %>%
    write_tsv(output, col_names = F)
}

rename_panel_tbprofiler <- function(panel) {
  # fastlin expect parent to be strict prefix
  orig <-
    panel %>%
    select(phylo = phylotype, parent = parent_phylotype) %>%
    distinct()

  new <- with(orig, str_remove(phylo, str_c('^', parent, '\\.')) %>%
                str_remove_all('\\.') %>%
                str_remove_all('_'))
  new <- str_c(new, '[', seq_along(new), ']')
  names(new) <- orig$phylo
  mod <-
    orig %>%
    mutate(phylo = unname(new[phylo]),
           parent = if_else(!is.na(new[parent]),
                            unname(new[parent]),
                            parent))

  parents <- filter(mod, parent == 'root') %>% pull(phylo)
  while (TRUE) {
    childs <- filter(mod, parent %in% parents) %>% pull(phylo)
    if (length(childs) == 0) {
      break
    }
    new <-
      mod %>%
      with(if_else(
        phylo %in% childs,
        str_c(parent, '.', phylo),
        phylo
      ))
    names(new) <- mod$phylo

    mod <-
      mod %>%
      mutate(phylo = new[phylo],
             parent = if_else(!is.na(new[parent]), new[parent], parent))

    parents <- unname(new[childs])
  }

  mod <-
    mod %>%
    mutate(across(everything(), \(x) str_remove_all(x, '\\[[0-9]+\\]')))

  old <- with(orig, unique(c(phylo, parent)))
  new <- with(mod, unique(c(phylo, parent)))
  names(new) <- old

  panel %>%
    mutate(phylotype = new[phylotype],
           parent_phylotype = new[parent_phylotype])

}
