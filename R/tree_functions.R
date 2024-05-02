# TODO
#   tidy functions, add imports, remove junk

tips <- function(phylo) {
  seq_len(treeio::Nnode2(phylo)) %>%
    {
      .[treeio::isTip(phylo, .)]
    }
}

phylo_labels <- function(phylo) {
  node_to_label(phylo, seq_len(treeio::Nnode2(phylo)))
}

inner_nodes <- function(phylo) {
  seq_len(treeio::Nnode2(phylo)) %>%
    {
      .[!treeio::isTip(phylo, .)]
    }
}

# add tbl_tree to class of dataframe
#' @importFrom rlang is_integerish is_character
as_tbl_tree <- function(df) {
  assert_that(
    is.data.frame(df),
    all(c("parent", "node", "label") %in% colnames(df)),
    is_integerish(df$parent),
    is_integerish(df$node),
    is_character(df$label)
  )

  class(df) <- union("tbl_tree", class(df))
  df
}

# convert a dataframe to phylo using tidy_tree
df2phylo <- function(df) {
  tidytree::as.phylo(as_tbl_tree(df))
}

# drop all tips from a phylo object
drop_all_tips <- function(phylo) {
  assert_that(is_phylo(phylo))

  parent <- phylo$edge[, 1]
  node <- phylo$edge[, 2]
  tip_set <- node[!(node %in% parent)]
  condense_phylo(phylo, tip_set)
}

# condense phylo object removing selected nodes
#' @importFrom rlang is_integerish
#' @importFrom dplyr filter select mutate bind_rows slice
#' @importFrom magrittr "%>%"
condense_phylo <- function(phylo, nodes_rm) {
  assert_that(
    is_phylo(phylo),
    is_integerish(nodes_rm)
  )

  tbl_tree <-
    with(as_tibble(phylo), {
      for (i in nodes_rm) {
        nd <- node[i]
        par <- parent[i]
        ch_i <- which(parent == nd)
        parent[ch_i] <- par
        branch.length[ch_i] <- branch.length[ch_i] + branch.length[i]
      }
      tibble(
        parent = parent,
        node = node,
        branch.length = branch.length,
        label = label
      ) %>%
        slice(-nodes_rm)
    })

  new_tip <-
    filter(tbl_tree, !node %in% parent) %>%
    pull(node)

  new_inner <- setdiff(tbl_tree$node, new_tip)

  new_order <-
    filter(tbl_tree, node %in% new_inner) %>%
    select(parent, node) %>%
    as.matrix() %>%
    t() %>%
    c() %>%
    unique()

  bind_rows(
    filter(tbl_tree, node %in% new_tip) %>%
      mutate(
        node = seq_along(new_tip),
        parent = length(new_tip) + match(parent, new_order)
      ),
    filter(tbl_tree, node %in% new_inner) %>%
      mutate(
        parent = length(new_tip) + match(parent, new_order),
        node = length(new_tip) + match(node, new_order)
      )
  ) %>%
    df2phylo()
}

#' @importFrom dplyr mutate rowwise ungroup n slice pull group_by select left_join distinct arrange
#' @importFrom tidyr gather
#' @importFrom magrittr "%>%"
collapse_phylotypes <- function(phylo, geno, min_dist = 1L, include_dist = TRUE) {
  assert_that(
    is_integer(geno) && is.matrix(geno) && max(geno, na.rm = T) == 1L && min(geno, na.rm = T) == 0L,
    is_phylo(phylo),
    setequal(rownames(geno), c(phylo$tip.label, phylo$node.label)),
    is_scalar_integerish(min_dist) && min_dist > 0L
  )

  phylo_dist <- phylo_geno_dist(phylo, geno, min_dist = min_dist, as_phylo = TRUE)

  node_dist <- ape::dist.nodes(phylo_dist)

  nodes_ident <-
    which(node_dist < min_dist & lower.tri(node_dist), arr.ind = T) %>%
    as_tibble() %>%
    mutate(id = seq_len(n())) %>%
    gather(row, col, key = "key", value = "node") %>%
    select(-key) %>%
    {
      left_join(., select(., node) %>% distinct() %>% mutate(depth = n_ancestor(phylo, node)), "node")
    } %>%
    group_by(id) %>%
    arrange(desc(depth), node) %>%
    slice(1) %>%
    ungroup() %>%
    pull(node) %>%
    unique()

  phylo_r1 <- `if`(include_dist, phylo_dist, phylo)
  phylo_r2 <- `if`(
    length(nodes_ident) > 0L,
    condense_phylo(phylo_r1, nodes_ident),
    phylo_r1
  )

  return(phylo_r2)
}

#' @importFrom dplyr mutate rowwise ungroup n slice pull group_by select left_join distinct arrange
#' @importFrom magrittr "%>%"
phylo_geno_dist <- function(phylo, geno, min_dist = 1L, as_phylo = FALSE) {
  # check args
  assert_that(
    is_integer(geno) && is.matrix(geno) && max(geno, na.rm = T) == 1L && min(geno, na.rm = T) == 0L,
    is_phylo(phylo),
    setequal(rownames(geno), c(phylo$tip.label, phylo$node.label)),
    is_bool(as_phylo)
  )

  phylo_dist <-
    tidytree::as_tibble(phylo) %>%
    mutate(parent_label = node_to_label(phylo, parent)) %>%
    rowwise() %>%
    mutate(branch.length = `if`(
      parent_label == "root", 1L,
      sum(geno[label, ] != geno[parent_label, ], na.rm = TRUE)
    )) %>%
    ungroup() %>%
    mutate(branch.length = replace(branch.length, node == parent, NA)) %>%
    as_tbl_tree() %>%
    tidytree::as.phylo()

  if (as_phylo) {
    return(phylo_dist)
  } else {
    return(ape::dist.nodes(phylo_dist))
  }
}

# return phylo with new labels
#' @importFrom rlang is_integerish is_dictionaryish is_string
#' @importFrom dplyr filter select mutate bind_rows slice
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
label_phylo <- function(phylo,
                        node_labels = integer(0L),
                        sep = "|",
                        anc_sep = "-A",
                        root_label = "root") {
  assert_that(
    is_phylo(phylo),
    is_integerish(node_labels) & is_dictionaryish(node_labels),
    all(node_labels %in% phylo$edge),
    is_string(sep)
  )

  rn <- treeio::rootnode(phylo)

  if (!rn %in% node_labels) {
    node_labels <- c(magrittr::set_names(rn, root_label), node_labels)
  }

  in_tree <-
    as_tibble(phylo) %>%
    mutate(
      label_new = names(node_labels)[match(node, node_labels)],
      children = map(node, function(nd) treeio::child(phylo, nd)),
      depth = n_ancestor(phylo, node),
      nchild = lengths(children),
      is_anc = (nchild == 1) & (!node %in% node_labels)
    )

  anc <-
    filter(in_tree, is_anc) %>%
    mutate(child = unlist(children)) %>%
    {
      `if`(
        nrow(.) > 0,
        with(., {
          level <- rep(1L, length(node))
          descendent <- child
          for (dp in unique(sort(setdiff(depth, max(depth)), T))) {
            ii <- which(depth == dp) %>% purrr::keep(~ child[.] %in% node)
            level[ii] <- level[match(child[ii], node)] + 1L
            descendent[ii] <- descendent[match(child[ii], node)]
          }
          tibble(label_old = label, level = level, descendent = descendent)
        }),
        tibble(label_old = character(), level = integer(), descendent = integer())
      )
    }

  new_labs <-
    condense_phylo(phylo, filter(in_tree, is_anc) %>% pull(node)) %>%
    {
      mutate(as_tibble(.),
        depth = n_ancestor(., node),
        children = map(node, function(nd) Children(., nd))
      )
    } %>%
    left_join(select(in_tree, label, label_new, nchild), "label") %>%
    with({
      label_new[parent == node] <- root_label
      for (dp in seq_len(max(depth) + 1L) - 1L) {
        slen <- length(label_new)
        for (i in which(depth == dp)) {
          ch <- children[[i]] %>% discard(~ !is.na(label_new[.]))
          label_new[ch] <- str_c(label_new[i], sep = sep, seq_along(ch))
        }
      }
      tibble(label_old = label, label_new = label_new)
    })

  in_tree %>%
    select(parent, node, branch.length, label_old = label) %>%
    left_join(new_labs, "label_old") %>%
    (function(df) {
      select(df, label_desc = label_new, descendent = node) %>%
        right_join(anc, "descendent") %>%
        mutate(label_anc = str_c(label_desc, anc_sep, level)) %>%
        select(label_old, label_anc) %>%
        {
          left_join(df, ., "label_old")
        }
    }) %>%
    mutate(label_new = if_else(is.na(label_new), label_anc, label_new)) %>%
    select(parent, node, branch.length, label = label_new) %>%
    df2phylo()
}

get_tips_by_nodelab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    filter(label %in% node_lab) %>%
    mutate(tips = get_tip_list(phylo, node)) %>%
    with(set_names(tips, label))
}

child_lab_by_lab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    (function(tree_tbl) {
      filter(tree_tbl, label %in% node_lab) %>%
        mutate(children = map(node, function(nd) tree_tbl$label[match(child(phylo, nd), tree_tbl$node)])) %>%
        with(set_names(children, label))
    })
}

child_dist_by_lab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    (function(tree_tbl) {
      filter(tree_tbl, label %in% node_lab) %>%
        mutate(distance = map(node, function(nd) tree_tbl$branch.length[match(child(phylo, nd), tree_tbl$node)])) %>%
        with(set_names(distance, label))
    })
}

parent_lab_by_lab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    (function(tree_tbl) {
      filter(tree_tbl, label %in% node_lab) %>%
        mutate(parents = map_chr(node, function(nd) tree_tbl$label[match(treeio::parent(phylo, nd), tree_tbl$node)])) %>%
        with(set_names(parents, label))
    })
}

# recursive function to count ancestors
n_ancestor <- function(phylo, node) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]

  n_ancestor_rec <- function(nd) {
    i <- match(nd, child)
    if (is.na(i)) {
      return(0L)
    }
    return(n_ancestor_rec(parent[i]) + 1L)
  }
  vapply(node, n_ancestor_rec, integer(1))
}


# recursive function to number of tips bellow a node
n_offspring <- function(phylo, node) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- !{
    child %in% parent
  }

  n_offspring_rec <- function(nd) {
    i <- which(nd == parent)
    vapply(i, function(j) {
      if (is_tip[j]) {
        return(1L)
      } else {
        return(n_offspring_rec(child[j]))
      }
    }, integer(1)) %>% sum()
  }
  vapply(node, n_offspring_rec, integer(1))
}

# count number of children
n_child <- function(phylo, node) {
  parent <- phylo$edge[, 1]
  parent <- parent[parent %in% node]
  vapply(node, function(n) sum(parent == n), integer(1L))
}

# recursive function to number of tips bellow a node
node_n_tip <- function(phylo, node) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- !{
    child %in% parent
  }

  node_n_tip_rec <- function(nd) {
    i <- which(nd == parent)
    vapply(i, function(j) {
      if (is_tip[j]) {
        return(1L)
      } else {
        return(m_node_n_tip_rec(child[j]))
      }
    }, integer(1)) %>% sum()
  }
  m_node_n_tip_rec <- memoise(node_n_tip_rec)
  vapply(node, m_node_n_tip_rec, integer(1))
}

# recursive function to get list of tips for each node
get_tip_list <- function(phylo, node, as_label = TRUE) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- !{
    child %in% parent
  }
  label <- phylo$tip.label

  get_tip_list_rec <- function(nd) {
    i <- which(nd == parent)

    map(i, function(j) {
      if (is_tip[j]) {
        return(j)
      } else {
        return(get_tip_list_rec(child[j]))
      }
    }) %>% unlist()
  }
  m_get_tip_list_rec <- memoise(get_tip_list_rec)
  l <- lapply(node, m_get_tip_list_rec)

  if (as_label) {
    lapply(l, function(i) label[i])
  } else {
    lapply(l, function(i) child[i])
  }
}

label_to_node <- function(phylo, labels) {
  tidytree::as_tibble(phylo) %>%
    with(node[match(labels, label)])
}

node_to_label <- function(phylo, nodes) {
  tidytree::as_tibble(phylo) %>%
    with(label[match(nodes, node)])
}

get_child_list <- function(phylo, node) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- !{
    child %in% parent
  }
  label <- phylo$tip.label

  get_child_list_rec <- function(nd) {
    i <- which(nd == parent)

    map(i, function(j) {
      if (is_tip[j]) {
        return(child[j])
      } else {
        return(c(child[j], m_get_child_list_rec(child[j])))
      }
    }) %>% unlist()
  }
  m_get_child_list_rec <- memoise(get_child_list_rec)
  map(node, m_get_child_list_rec)
}

get_children <- function(phylo, node, depth = 1L) {
  assert_that(
    is_phylo(phylo),
    is_scalar_integerish(node),
    is_scalar_integerish(depth) & depth >= 0L,
    node %in% phylo$edge
  )

  parent <- phylo$edge[, 1]
  children <- phylo$edge[, 2]

  get_children_rec <- function(nd, depth) {
    ch <- children[which(nd == parent)]
    if ((depth > 1L) & (length(ch) > 0L)) {
      map(ch, function(ch_) {
        get_children_rec(ch_, depth - 1L)
      }) %>%
        unlist() %>%
        {
          c(ch, .)
        }
    } else {
      ch
    }
  }
  if (depth > 0L) {
    get_children_rec(node, depth)
  } else {
    integer(0)
  }
}

get_child_list_by_lab <- function(phylo, labels) {
  label_to_node(phylo, labels) %>%
    get_child_list(phylo, .) %>%
    map(~ node_to_label(phylo, .)) %>%
    set_names(labels)
}

dist_to_tip <- function(phylo, node) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- !{
    child %in% parent
  }

  dist_to_tip_rec <- function(nd) {
    i <- which(nd == parent)

    map_int(i, function(j) {
      if (is_tip[j]) {
        return(1L)
      } else {
        return(m_dist_to_tip_rec(child[j]) + 1L)
      }
    }) %>% min()
  }
  m_dist_to_tip_rec <- memoise(dist_to_tip_rec)
  map_int(node, m_dist_to_tip_rec)
}


get_nodes <- function(phylo) {
  phylo$edge[, 1] %>%
    sort() %>%
    unique()
}

get_tips <- function(phylo) {
  phylo$edge[, 2] %>%
    sort() %>%
    setdiff(get_nodes(phylo))
}

get_root <- function(phylo) {
  setdiff(unique(phylo$edge[, 1]), phylo$edge[, 2])
}

# recursive function to get ancestors
get_anc_list <- function(phylo, node, as_label = FALSE) {
  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]

  if (as_label) {
    labels <- c(phylo$tip.label, phylo$node.label)
    get_anc_rec <- function(nd) {
      p <- parent[match(nd, child)]
      if (is.na(p)) {
        return(character(0))
      }
      return(c(get_anc_rec(p), labels[p]))
    }
  } else {
    get_anc_rec <- function(nd) {
      p <- parent[match(nd, child)]
      if (is.na(p)) {
        return(integer(0))
      }
      return(c(get_anc_rec(p), p))
    }
  }
  m_get_anc_rec <- memoise::memoise(get_anc_rec)
  map(node, m_get_anc_rec)
}

