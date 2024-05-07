# Refine SNP panel (barcode) based on observed data

tbrefine <- function(
    gds = NULL,
    vcf = NULL,
    sample_phylotypes = NULL,
    add_new = TRUE,
    remove_old = TRUE,
    panel = TBtypeR::tbt_panel,
    min_ac = 2L,
    min_snp_node = 5L,
    max_var_missing = 0.05,
    min_ext_fold = 5,
    max_global_discordand_rate = 0.01,
    max_local_discordant_rate = 0.05,
    min_discordant_sample = 2L,
    ...) {

  # check args
  assert_that(
    is_gds(gds) || is_scalar_character(vcf),
    is.data.frame(sample_phylotypes) || is.null(sample_phylotypes),
    is_scalar_integerish(min_ac) && min_ac > 0,
    is_scalar_integerish(min_snp_node) && min_snp_node > 0,
    is_scalar_double(max_var_missing),
    is_scalar_double(min_ext_fold),
    is_scalar_double(max_global_discordand_rate),
    is_scalar_double(max_local_discordant_rate),
    is_scalar_integerish(min_discordant_sample) && min_discordant_sample > 0
  )

  if (is_scalar_character(vcf)) {
    gds <- vcf_to_gds(vcf)
  }

  if (is.null(sample_phylotypes )) {
    sample_phylotypes <-
      tbtype(gds = gds, ...) %>%
      filter_tbtype(max_phylotypes = 2) %>%
      filter(n_phy == 1) %>%
      unnest_mixtures() %>%
      select(sample_id, phylotype = mix_phylotype)
  }

  sample_phylotypes <-
    sample_phylotypes %>%
    semi_join(panel, by = 'phylotype')

  new_panel_snps <- tibble()
  rm_snps <- character()

  if (add_new) {
    # add new SNPs consistent with existing panel

    # exclude panel sites
    seqResetFilter(gds)
    seqSetFilter(gds, !seqGetData(gds, 'position') %in% tbt_panel$pos)
    seqSetFilter(gds, sample.id = sample_phylotypes$sample_id)

    # get sample genotypes
    gt12 <- seqGetData(gds, 'genotype')
    gt <- gt12[1,,]
    gt[gt != gt12[2,,]] <- NA_integer_
    rownames(gt) <- seqGetData(gds, 'sample.id')
    colnames(gt) <- seqGetData(gds, 'variant.id')
    rm(gt12)

    # filter missing sites
    keep <- colSums(!is.na(gt)) / nrow(gt) > max_var_missing
    seqSetFilter(gds, variant.id = seqGetData(gds, 'variant.id')[keep])
    gt <- gt[, keep]

    cand_snps <-
      setNames(0:3, 0:3) %>%
      map_dfc(\(x) colSums(gt == x, na.rm = T)) %>%
      mutate(vid = as.integer(colnames(gt)), .before = `0`)  %>%
      pivot_longer(-vid, names_to = 'aid', values_to = 'count') %>%
      mutate(aid = as.integer(aid)) %>%
      arrange(vid, desc(count)) %>%
      group_by(vid) %>%
      mutate(freq = count / sum(count)) %>%
      dplyr::slice(1:2) %>%
      mutate(allele = c('a', 'b')) %>%
      ungroup() %>%
      pivot_wider(names_from = allele, values_from = c(aid, count, freq)) %>%
      mutate(freq_other = 1 - freq_a - freq_b) %>%
      filter(freq_a < 1,
             count_b >= min_ac,
             freq_b > min_ext_fold * (1-freq_a-freq_b)) %>%
      select(-starts_with(c('count', 'freq'))) %>%
      inner_join(
        variantInfo(gds) %>%
          as_tibble() %>%
          rename(vid = variant.id,
                 chrom = chr),
        by = join_by(vid)
      ) %>%
      # only allow A allele == reference allele
      filter(aid_a == 0L) %>%
      mutate(alt = map2_chr(alt, aid_b, \(alt, i) { str_split(alt, ',', n=3, T)[i] }))

    # new SNPs present in at most one phylotype
    cand_snp_samp_phy <-
      cand_snps %>%
      mutate(sample_id = map2(vid, aid_b, function(vid, aid_b) {
        names(which(gt[, as.character(vid)] == aid_b))
      })) %>%
      select(vid, sample_id) %>%
      unnest(sample_id) %>%
      inner_join(sample_phylotypes, by = 'sample_id') %>%
      chop(sample_id) %>%
      group_by(vid) %>%
      filter(n() == 1) %>%
      ungroup() %>%
      unnest(sample_id)

    new_snps <-
      cand_snp_samp_phy %>%
      split.data.frame(.$phylotype) %>%
      map_dfr(function(data) {

        root_phylotype <- data$phylotype[1]

        snp_tbl <-
          data %>%
          select(vid, sample_id) %>%
          mutate(gt = 1L) %>%
          complete(vid, sample_id, fill=list(gt = 0L))

        phylo <-
          snp_tbl %>%
          pivot_wider(
            names_from = sample_id,
            values_from = gt) %>%
          mutate(ref = 0L) %>%
          as.data.frame() %>%
          column_to_rownames('vid') %>%
          t() %>%
          dist() %>%
          hclust() %>%
          ape::as.phylo() %>%
          ape::root('ref') %>%
          as_tibble() %>%
          mutate(label = if_else(is.na(label), str_c('n', node), label)) %>%
          ape::as.phylo()

        node_sample <-
          phylo %>%
          as_tibble() %>%
          select(node) %>%
          mutate(
            sample_id =
              phangorn::Descendants(phylo, node) %>%
              map(~ na.omit(TBtypeR:::node_to_label(phylo, .))),
            depth = lengths(phangorn::Ancestors(phylo, node))
          ) %>%
          unnest(sample_id)

        node_vars <-
          node_sample %>%
          inner_join(snp_tbl, relationship = 'many-to-many', by = 'sample_id') %>%
          count(node, depth, vid, gt) %>%
          pivot_wider(names_from = gt,
                      names_prefix = 'n_',
                      values_from = n,
                      values_fill = list(n = 0L)) %>%
          mutate(n_0 = `if`(exists('n_0'), n_0, 0L)) %>%
          group_by(vid) %>%
          # limit to clades with most nodes containing
          filter(n_1 == max(n_1)) %>%
          ungroup() %>%
          filter(n_1 > 0, n_0 == 0) %>%
          select(-n_1, -n_0) %>%
          arrange_all() %>%
          chop(vid) %>%
          group_by(vid) %>%
          filter(depth == max(depth)) %>%
          ungroup() %>%
          mutate(label = node_to_label(phylo, node)) %>%
          filter(lengths(vid) >= min_snp_node)

        if (nrow(node_vars) == 0) {
          return(tibble())
        }

        node_labels <-
          condense_phylo(
            phylo,
            setdiff(seq_len(treeio::Nnode2(phylo)),
                    c(node_vars$node, treeio::rootnode(phylo))
            )
          ) %>%
          (function(x) {
            as_tibble(x) %>%
              mutate(
                phylotype =
                  label_phylo(x, root_label = '') %>%
                  as_tibble() %>%
                  pull(label) %>%
                  str_remove('^\\.'),
                parent_phylotype = phylotype[match(parent, node)]) %>%
              select(label, phylotype, parent_phylotype)
          }) %>%
          mutate(across(ends_with('phylotype'), ~ str_c(root_phylotype, '[', ., ']'))) %>%
          mutate(parent_phylotype = str_remove(parent_phylotype, '\\[\\]'))

        inner_join(
          node_labels,
          node_vars,
          by = 'label'
        ) %>%
          select(phylotype, parent_phylotype, vid)
      })

    new_panel_snps <-
      new_snps %>%
      unnest(vid) %>%
      inner_join(cand_snps, by = 'vid') %>%
      mutate(genotype = 1L,
             reference = 'TBrefineR') %>%
      select(chrom, pos, ref, alt, genotype, phylotype, parent_phylotype, reference)

  }

  if (remove_old) {
    # remove SNPs inconsistent with observed data

    # exclude panel sites
    seqResetFilter(gds)
    seqSetFilter(gds, seqGetData(gds, 'position') %in% tbt_panel$pos)
    seqSetFilter(gds, sample.id = sample_phylotypes$sample_id)

    # get sample genotypes
    gt12 <- seqGetData(gds, 'genotype')
    gt <- gt12[1,,]
    gt[gt != gt12[2,,]] <- -1L
    rownames(gt) <- seqGetData(gds, 'sample.id')
    colnames(gt) <-
      tibble(pos = seqGetData(gds, 'position')) %>%
        left_join(panel_with_vid(panel) %>% select(vid, pos),
                  by = 'pos') %>%
        pull(vid)
    rm(gt12)

    # encode as 0/1
    as_tibble(variantInfo(gds)) %>%
      select(pos, ref, alt) %>%
      separate_rows(alt) %>%
      group_by(pos, ref) %>%
      mutate(GT = seq_along(alt)) %>%
      ungroup() %>%
      inner_join(panel_with_vid(panel), by = c('pos', 'ref', 'alt')) %>%
      select(GT, vid) %>%
      chop(vid) %>%
      pwalk(function(GT, vid) {
        gt[, vid][gt[, vid] == GT] <<- 1L
      })
    gt[gt > 1] <- -1L

    # global_discordance rates

    lin_geno <- panel_to_geno(panel)[, colnames(gt)]

    gt_exp <- gt
    gt_exp[] <- 0L

    pwalk(sample_phylotypes, function(sample_id, phylotype) {
      gt_exp[sample_id, ] <<- lin_geno[phylotype, ]
    })

    rm_global <-
      which(colSums(gt != gt_exp, na.rm = T) / colSums(!is.na(gt)) > max_global_discordand_rate) %>%
      names()

    # local_discordance
    phylo <- panel_to_phylo(panel)

    rm_local <-
      phylo %>%
      as_tibble() %>%
      filter(label != 'root') %>%
      mutate(descendent =
               Descendants(phylo, node, 'all') %>%
               map( ~ node_to_label(phylo, .))
             ) %>%
      select(phylotype = label, descendent) %>%
      inner_join(
        select(panel_with_vid(panel), vid, pos, genotype, phylotype),
        by = 'phylotype') %>%
      unnest(descendent) %>%
      inner_join(
        select(sample_phylotypes, descendent=phylotype, sample_id),
        by = 'descendent',
        relationship = 'many-to-many') %>%
      select(phylotype, vid, pos, sample_id, genotype) %>%
      distinct() %>%
      group_by(sample_id) %>%
      mutate(correct = genotype == gt[sample_id[1], vid]) %>%
      ungroup() %>%
      count(phylotype, vid, correct) %>%
      filter(!is.na(correct)) %>%
      group_by(phylotype, vid) %>%
      summarise(
        n_incorrect = sum(n[!correct]),
        n = sum(n),
        discordance = n_incorrect / n) %>%
      filter(
        n_incorrect > min_discordant_sample,
        discordance > max_local_discordant_rate) %>%
      pull(vid)

    rm_snps <- union(rm_local, rm_global)
  }


  if (nrow(new_panel_snps)) {
    n_new <- n_distinct(new_panel_snps$phylotype)
    message('Added ', nrow(new_panel_snps), ' SNPs and ',  n_new, ' strains to panel')
  }
  if (length(rm_snps)) {
    message('Removed ', length(rm_snps), ' SNPs from panel')
  }

  refined_panel <-
    panel_with_vid(panel) %>%
    filter(!vid %in% rm_snps) %>%
    select(-vid) %>%
    bind_rows(new_panel_snps) %>%
    arrange_all()

}

