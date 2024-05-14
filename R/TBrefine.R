# Refine SNP panel (barcode) based on observed data

#' @importFrom SeqArray seqGetData seqSetFilter seqResetFilter
#' @importFrom purrr map_dfc pwalk
#' @importFrom dplyr count
tbrefine <- function(
    gds = NULL,
    vcf = NULL,
    sample_phylotypes = NULL,
    add_new = TRUE,
    remove_old = TRUE,
    panel = TBtypeR::tbt_panel,
    MIN_AC = 2L,
    MIN_SNP_NODE = 5L,
    MIN_VAF_NODE = 0.95,
    MAX_VAR_MISSING = 0.05,
    MIN_C_FOLD = 5,
    GLOBAL_MAX_ERROR_RATE = 0.01,
    LOCAL_MAX_ERROR_RATE = 0.05,
    LOCAL_MIN_SAMPLE_ERROR = 2L,
    ...) {

  # check args
  assert_that(
    is_gds(gds) || is_scalar_character(vcf),
    is.data.frame(sample_phylotypes) || is.null(sample_phylotypes),
    is_scalar_integerish(MIN_AC) && MIN_AC > 0,
    is_scalar_integerish(MIN_SNP_NODE) && MIN_SNP_NODE > 0,
    is_scalar_double(MAX_VAR_MISSING),
    is_scalar_double(MIN_C_FOLD),
    is_scalar_double(GLOBAL_MAX_ERROR_RATE),
    is_scalar_double(LOCAL_MAX_ERROR_RATE),
    is_scalar_integerish(LOCAL_MIN_SAMPLE_ERROR) && LOCAL_MIN_SAMPLE_ERROR > 0
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
  new_parent_phylo <- tibble(phylotype = character(), new_parent_phylotype = character())
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
    keep <- colSums(!is.na(gt)) / nrow(gt) > MAX_VAR_MISSING
    seqSetFilter(gds, variant.id = seqGetData(gds, 'variant.id')[keep])
    gt <- gt[, keep]

    # set of biallelic SNPs
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
             count_b >= MIN_AC,
             freq_b > MIN_C_FOLD * (1-freq_a-freq_b)) %>%
      select(-starts_with(c('count', 'freq'))) %>%
      inner_join(
        variantInfo(gds) %>%
          as_tibble() %>%
          rename(vid = variant.id,
                 chrom = chr),
        by = 'vid'
      ) %>%
      # only allow A allele == reference allele
      filter(aid_a == 0L) %>%
      mutate(alt = map2_chr(alt, aid_b, \(alt, i) { str_split(alt, ',', n=3, T)[i] }))

    # binarise sample gts at candidate sites
    gt_cand <- gt[, as.character(cand_snps$vid)]

    cand_snps %>%
      select(vid, aid_b) %>%
      mutate(vid = as.character(vid)) %>%
      chop(vid) %>%
      pwalk(function(vid, aid_b) {
        b <- which(gt_cand[, vid] == aid_b)
        ot <- setdiff(which(gt_cand[, vid] > 0 ), b)
        if (length(ot)) {
          gt_cand[, vid][ot] <<- 0L
        }
        gt_cand[, vid][b] <<- 1L
      })

    # create sample phylo to get ancestors etc
    sample_phylo <-
      panel %>%
      bind_rows(
        sample_phylotypes %>%
          transmute(chrom = 'chrom',
                    pos = 1L,
                    ref = 'A',
                    alt = 'A',
                    genotype = 1L,
                    parent_phylotype = phylotype,
                    phylotype = sample_id)
      ) %>%
      panel_to_phylo()

    # samples for each branch
    branch_sample <-
      sample_phylotypes %>%
      select(sample_id) %>%
      mutate(node = label_to_node(sample_phylo, sample_id)) %>%
      mutate(branch = phangorn::Ancestors(sample_phylo, node, 'all')) %>%
      select(-node) %>%
      unnest(branch) %>%
      chop(sample_id) %>%
      mutate(branch = node_to_label(sample_phylo, branch))

    # AC/AF by branch
    branch_ac <- matrix(
      0L, nrow = nrow(branch_sample), ncol = nrow(cand_snps),
      dimnames = list(branch = branch_sample$branch, variant = cand_snps$vid))
    branch_af <- branch_ac

    branch_sample %>%
      pwalk(function(branch, sample_id){
        branch_ac[branch, ] <<- colSums(gt_cand[sample_id, , drop=F], na.rm = TRUE)
        branch_af[branch, ] <<-  branch_ac[branch, ] / colSums(!is.na(gt_cand[sample_id, ,drop=F]))
      })

    branch_snps <-
      tibble(branch = rownames(branch_ac)) %>%
      mutate(
        node = label_to_node(sample_phylo, branch),
        depth = lengths(phangorn::Ancestors(sample_phylo, node, 'all'))
      ) %>%
      mutate(snp_data = map(branch, function(b) {
        tibble(vid = names(which(branch_ac[b, ] > 0)),
               ac = branch_ac[b, vid],
               af = branch_af[b, vid]) %>%
          mutate(vid = as.integer(vid))
      })) %>%
      unnest(snp_data) %>%
      arrange(vid, branch) %>%
      group_by(vid) %>%
      mutate(is_root = ac == max(ac),
             is_root = is_root & depth == max(depth[is_root])) %>%
      filter(any(is_root)) %>%
      filter(depth >= depth[is_root]) %>%
      rename(sub_branch = branch,
             sub_node = node) %>%
      mutate(branch = sub_branch[is_root],
             node = sub_node[is_root]) %>%
      ungroup() %>%
      select(-is_root) %>%
      chop(-c(branch, node, vid)) %>%
        chop(-node) %>%
        mutate(
          children =
            Children(sample_phylo, node) %>%
            map(~ node_to_label(sample_phylo, .)) %>%
            map(setdiff, sample_phylotypes$sample_id)
        ) %>%
        unchop(-c(node, children)) %>%
      mutate(children_obs = map2(children, sub_branch, ~ sort(intersect(.x, unlist(.y))))) %>%
      mutate(case = case_when(
        lengths(children_obs) == 0                 ~ 'no_children',
        lengths(children_obs) == lengths(children) ~ 'all_children',
        lengths(children_obs)  < lengths(children) ~ 'some_children'
      ))

    ##### CASE TIP_NODE/NO_CHILDREN #####
    # filter new nodes with MIN_SNP_NODE
    # perform heirachical clustering of SNPS in each branch
    # remove SNPs inconsistent with clusters
    new_tip_snps <-
      branch_snps %>%
      filter(case == 'no_children') %>%
      select(branch, vid) %>%
      arrange(branch, vid) %>%
      mutate(vid = as.character(vid)) %>%
      chop(vid) %>%
      filter(lengths(vid) > MIN_SNP_NODE) %>%
      mutate(data = map2(vid, seq_along(vid),  function(vid, i) {

        default <- tibble(phylotype = character(),
                          parent_phylotype = character(),
                          vid = list())

        snp_mat <- gt_cand[which(rowSums(gt_cand[, vid], na.rm=T) > 0), vid, drop = FALSE]
        snp_mat <- rbind(
          matrix(0L, ncol = ncol(snp_mat), dimnames = list('ref', colnames(snp_mat))),
          snp_mat
        )
        if (any(is.na(snp_mat))) {
          # filter too many missing
          missingness <- colSums(is.na(snp_mat)) / nrow(snp_mat)
          snp_mat <- snp_mat[, missingness < MAX_VAR_MISSING, drop = F]
          if (ncol(snp_mat) == 0) {
            return(default)
          }
        }
        if (any(is.na(snp_mat))) {
          # impute remaining missing
          imp <- missMDA::imputePCA(snp_mat)$completeObs
          missing <- which(is.na(snp_mat))
          snp_mat[missing] <- imp[missing]  %>% (\(x) { if_else(x < 0.5, 0L, 1L) })
        }

        snp_tbl <-
          as.data.frame(snp_mat) %>%
          rownames_to_column('sample_id') %>%
          pivot_longer(-sample_id, names_to = 'vid', values_to = 'gt')

        phylo <-
          snp_mat %>%
          dist() %>%
          hclust() %>%
          ape::as.phylo() %>%
          ape::root('ref') %>%
          as_tibble() %>%
          mutate(label = if_else(is.na(label), str_c('n', node), label)) %>%
          tidytree::as.phylo()

        # filter variants for monophyletic clusters
        node_vars <-
          phylo %>%
          as_tibble_shh() %>%
          select(node) %>%
          mutate(
            sample_id =
              phangorn::Descendants(phylo, node) %>%
              map(~ na.omit(TBtypeR:::node_to_label(phylo, .))),
            depth = lengths(phangorn::Ancestors(phylo, node))
          ) %>%
          unnest(sample_id) %>%
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
          ungroup()  %>%
          filter(n_1 > 0, n_0 == 0) %>%
          select(-n_1, -n_0) %>%
          arrange_all() %>%
          chop(vid) %>%
          group_by(vid) %>%
          filter(depth == max(depth)) %>%
          ungroup() %>%
          mutate(label = node_to_label(phylo, node)) %>%
          filter(lengths(vid) >= MIN_SNP_NODE)

        if (nrow(node_vars) == 0) {
          return(default)
        }

        node_labels <-
          condense_phylo(
            phylo,
            setdiff(
              seq_len(treeio::Nnode2(phylo)),
              c(node_vars$node, treeio::rootnode(phylo), seq_len(treeio::Ntip(phylo))))
          ) %>%
          (function(x) {
            as_tibble_shh(x) %>%
              mutate(
                phylotype =
                  label_phylo(
                    x,
                    root_label = if_else(treeio::rootnode(phylo) %in% node_vars$node, '1', ''),
                    symbols = letters,
                    sep = '',
                    node_labels = seq_len(treeio::Ntip(phylo)) %>% setNames(., node_to_label(phylo, .))
                  ) %>%
                  as_tibble_shh() %>%
                  pull(label) %>%
                  str_remove('^\\.'),
                parent_phylotype = phylotype[match(parent, node)]) %>%
              select(label, phylotype, parent_phylotype)
          }) %>%
          mutate(across(ends_with('phylotype'), ~ str_c('[', ., ']')))

        inner_join(
          node_labels,
          node_vars,
          by = 'label') %>%
          select(phylotype, parent_phylotype, vid) %>%
          mutate(vid = map(vid, as.integer))

      })) %>%
      select(-vid) %>%
      unnest(data)

    if (nrow(new_tip_snps)) {
      new_panel_snps <-
        new_tip_snps %>%
        mutate(across(ends_with('phylotype'), ~ str_c(branch, .))) %>%
        mutate(parent_phylotype = str_remove(parent_phylotype,'\\[\\]')) %>%
        select(phylotype, parent_phylotype, vid)
    }


    ##### CASE ALL CHILDREN #####
    # filter MIN_VAF_NODE in all sub_branches
    new_node_snps <-
      branch_snps %>%
      filter(case == 'all_children') %>%
      filter(map_lgl(af, ~ all(. >= MIN_VAF_NODE))) %>%
      select(branch, vid) %>%
      chop(vid) %>%
      rename(phylotype = branch) %>%
      inner_join(
        panel %>%
          select(phylotype, parent_phylotype) %>%
          distinct(),
        by = 'phylotype'
      )

    if (nrow(new_node_snps)) {
      new_panel_snps <- bind_rows(new_panel_snps, new_node_snps)
    }

    ##### CASE SOME CHILDREN #####
    # filter MIN_VAF_NODE in all sub_branches
    # filter for MIN_SNP_NODE
    # filter for non-overlapping supersets
    # add ancestral node(s)

    new_ancestor_snps <-
      branch_snps %>%
      filter(case == 'some_children') %>%
      filter(map_lgl(af, ~ all(. >= MIN_VAF_NODE))) %>%
      select(vid, branch, sub_branch, depth, ac) %>%
      unnest(c(sub_branch, depth, ac)) %>%
      group_by(vid) %>%
      filter(depth == 1 + suppressWarnings(min(depth))) %>%
      ungroup() %>%
      select(-depth) %>%
      chop(c(sub_branch, ac)) %>%
      mutate(ac = map_int(ac, sum)) %>%
      chop(c(ac, vid)) %>%
      mutate(ac = map_int(ac, sum)) %>%
      arrange(branch, vid) %>%
      filter(lengths(vid) > MIN_SNP_NODE) %>%
      group_by(branch) %>%
      filter(valid_nodesets(sub_branch, ac)) %>%
      mutate(labels = label_nodesets(sub_branch)) %>%
      ungroup() %>%
      unnest(labels) %>%
      mutate(across(ends_with('phylotype'), ~ str_c(branch, '[', ., ']'))) %>%
      mutate(parent_phylotype = str_remove(parent_phylotype, '\\[\\]$')) %>%
      select(phylotype, parent_phylotype, vid, sub_branch)

    if (nrow(new_ancestor_snps)) {
      new_parent_phylo <-
        new_ancestor_snps %>%
        select(new_parent_phylotype = phylotype,
               phylotype = sub_branch) %>%
        mutate(n = lengths(phylotype)) %>%
        unnest(phylotype) %>%
        group_by(phylotype) %>%
        slice(which.min(n)) %>%
        select(-n)

      new_panel_snps <-
        bind_rows(
          new_panel_snps,
          select(new_ancestor_snps, -sub_branch)
        )
    }

    new_panel_snps <-
      new_panel_snps %>%
      unnest(vid) %>%
      inner_join(cand_snps, by = 'vid') %>%
      mutate(genotype = 1L,
             reference = 'TBrefineR') %>%
      select(chrom, pos, ref, alt, genotype, phylotype, parent_phylotype, reference)
  }

  if (remove_old) {
    # remove SNPs inconsistent with observed data

    # select panel sites
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
      which(colSums(gt != gt_exp, na.rm = T) / colSums(!is.na(gt)) > GLOBAL_MAX_ERROR_RATE) %>%
      names()

    panel_phylo <- panel_to_phylo(panel)

    # local_discordance
    rm_local <-
      panel_phylo %>%
      as_tibble() %>%
      filter(label != 'root') %>%
      mutate(descendant =
               Descendants(panel_phylo, node, 'all') %>%
               map( ~ node_to_label(panel_phylo, .))
             ) %>%
      select(phylotype = label, descendant) %>%
      inner_join(
        select(panel_with_vid(panel), vid, pos, genotype, phylotype),
        by = 'phylotype') %>%
      unnest(descendant) %>%
      inner_join(
        select(sample_phylotypes, descendant=phylotype, sample_id),
        by = 'descendant',
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
        n_incorrect > LOCAL_MIN_SAMPLE_ERROR,
        discordance > LOCAL_MAX_ERROR_RATE) %>%
      pull(vid)

    rm_snps <- union(rm_local, rm_global)
  }


  if (length(rm_snps)) {
    message('Removed ', length(rm_snps), ' SNPs from panel')
  }
  x <- new_panel_snps %>% filter(phylotype %in% panel$phylotype)
  if (nrow(x)) {
    message('Added ', nrow(x), ' SNPs to ', dplyr::n_distinct(x$phylotype), ' existing strains in panel')
  }
  x <- new_panel_snps %>% filter(!phylotype %in% panel$phylotype)
  if (nrow(new_panel_snps)) {
    message('Added ', nrow(x), ' SNPs and ',  dplyr::n_distinct(x$phylotype), ' new strains to panel')
  }

  refined_panel <-
    panel_with_vid(panel) %>%
    filter(!vid %in% rm_snps) %>%
    select(-vid) %>%
    bind_rows(new_panel_snps) %>%
    left_join(new_parent_phylo, by = 'phylotype') %>%
    mutate(
      parent_phylotype = if_else(
        is.na(new_parent_phylotype),
        parent_phylotype,
        new_parent_phylotype)
    ) %>%
    select(-new_parent_phylotype) %>%
    arrange_all()

  # check validity
  panel_to_phylo(refined_panel, T)

  return(refined_panel)

}

# function to determine maximal compatible set of subbranch nodesets
valid_nodesets <- function(nodesets, weights = rep(1L, length(nodesets))) {

  if (length(nodesets) < 2) { return(rep(TRUE, length(nodesets))) }

  superset_id <- function(subsets) {
    if (length(subsets) == 1) { return(1L)}
    supersets <-
      reduce(subsets, function(x, y) {
        X <- `if`(is.list(x), x, list(x))
        if (any(y %in% unlist(X))) {
          map(X, ~ `if`(any(y %in% .), union(y, .), .))
        } else {
          c(X, list(y))
        }
      })
    map_int(subsets, \(x) which(map_lgl(supersets, \(y) any(x %in% y))))
  }

  # return T/F for length of .nodesets

  recursion <- function(.nodesets) {
    if (length(.nodesets) == 1L) { return(TRUE) }
    ssid <- superset_id(.nodesets)
    if (dplyr::n_distinct(ssid) == 1L) {
      # quick check for validity
      largest <- which(lengths(.nodesets) == max(lengths(.nodesets)))[1]
      if (all(map_lgl(.nodesets, ~ all(. %in% .nodesets[[largest]])))) {
        val <- all(recursion(.nodesets[-largest]))
        return(rep(val, length(.nodesets)))
      } else {
        return(rep(FALSE, length(.nodesets)))
      }
    } else {
      val <- map_lgl(seq_len(max(ssid)), ~all(recursion(.nodesets[ssid == .])))
      val[ssid]
    }
  }

  # preseve the most weights
  valid <- recursion(nodesets)
  if (any (!valid)) {
    combs <-
      tibble(keep = combinations(which(!valid), zero_len_as_na = F)) %>%
      filter(lengths(keep) > 0, lengths(keep) < max(lengths(keep))) %>%
      mutate(weight = map_dbl(keep, ~ sum(weights[.]))) %>%
      arrange(desc(weight))
    nodesets_keep <- integer(0)
    for (k in combs$keep) {
      if (all(recursion(nodesets[k]))) {
        nodesets_keep <- k
        break
      }
    }
    valid[nodesets_keep] <- TRUE
  }
  return(valid)
}

label_nodesets <- function(nodeset) {

  if(length(nodeset) == 1) {
    return(tibble(phylotype = 'A', parent_phylotype = ''))
  } else if(length(nodeset) == 0) {
    return(tibble(phylotype = character(), parent_phylotype = character()))
  }

  X <-
    tibble(nodeset) %>%
    mutate(rn = row_number()) %>%
    arrange(lengths(nodeset)) %>%
    mutate(i = row_number()) %>%
    mutate(parent = purrr::map2_int(i, nodeset, function(i, ns) {
      if (i == length(nodeset)) { return( NA_integer_ )}
      map_lgl(
        nodeset[seq.int(i+1, length(nodeset))],
        function(NS) all(ns %in% NS)) %>%
        which() %>%
        first() %>%
        (\(x) i + x)
    })) %>%
    mutate(phylotype = str_c('node', i),
           parent_phylotype = str_c('node', parent) %>% replace_na('root'))

  X %>%
    panel_to_phylo(minimal_checks = T) %>%
    (function(x) {
      inner_join(
        as_tibble_shh(x) %>%
          select(node, label),
        label_phylo(x, symbols = LETTERS, root_label = '', sep = '') %>%
          as_tibble_shh() %>%
          select(node, new_label = label),
        by = 'node'
      ) %>%
        with(setNames(new_label, label))
    }) %>%
    (function(x) mutate(X, across(ends_with('phylotype'), ~x[.]))) %>%
    arrange(rn) %>%
    select(phylotype, parent_phylotype)
}

