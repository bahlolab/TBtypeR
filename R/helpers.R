
is_proportion <- function(x, na.rm = T) {
  all((0 <= x) & (1 >= x), na.rm = na.rm)
}

#' @importFrom rlang is_scalar_double is_scalar_integer
is_scalar_proportion <- function(x, na.rm = T) {
  (is_scalar_double(x) || is_scalar_integer(x)) && is_proportion(x, na.rm = na.rm)
}


is_phylo <- function(x) {
  inherits(x, "phylo")
}

#' @importFrom rlang is_double
force_to_interval <- function(x, min_val = 0, max_val = 1) {
  assert_that(is_double(x))
  x[x < min_val] <- min_val
  x[x > max_val] <- max_val
  return(x)
}

#' @importFrom rlang is_double is_scalar_double
c_spike_in <- function(p0, p1) {
  assert_that(
    is_double(p0),
    is_scalar_double(p1)
  )
  return(c(p0 * (1 - p1), p1))
}



permute_genotypes <- function(phy_gts, n_perm) {
  rate <- rowSums(phy_gts) / ncol(phy_gts)
  # "permuted" genotypes background
  vapply(
    seq_len(nrow(phy_gts)),
    function(i) {
      rbinom(n_perm, 1, rate[i])
    },
    integer(n_perm)
  ) %>% t()
}

#' @importFrom dplyr mutate select summarise starts_with row_number
#' @importFrom tidyr chop pivot_longer
#' @importFrom magrittr set_colnames
permute_similar <- function(states, new_gts, n_perm) {
  assert_that(
    is.factor(states),
    is_integerish(new_gts),
    length(states) == length(new_gts),
    is_scalar_integerish(n_perm)
  )

  state_smry <-
    tibble(state = states, new_gt = new_gts) %>%
    mutate(site = row_number()) %>%
    group_by(state) %>%
    summarise(
      N = sum(new_gt),
      n_site = length(site),
      site = list(site)
    )

  perm <- matrix(0L, nrow = length(states), ncol = n_perm)
  ones <- state_smry %>%
    filter(N == n_site) %>%
    pull(site) %>%
    unlist()
  if (!is.null(ones)) {
    perm[ones, ] <- 1L
  }
  for (i in with(state_smry, which(N > 1 & N < n_site))) {
    for (j in seq_len(n_perm)) {
      perm[sample(state_smry$site[[i]], state_smry$N[i]), j] <- 1L
    }
  }
  return(perm)
}


#' @importFrom rlang is_scalar_double is_scalar_integer
binom_likelihood <- function(x, size, p, err, by_site = FALSE) {
  # check args
  assert_that(
    is_integerish(x),
    is_integerish(size),
    length(x) == length(size),
    length(x) == length(p) || length(p) == 1L,
    is_proportion(p),
    is_proportion(err),
    length(err) == 1 || length(err) == length(x)
  )
  pe <- p * (1 - err) + (1 - p) * (err)
  res <- dbinom(x, size, pe, log = TRUE)
  if (by_site) {
    return(res)
  } else {
    return(sum(res, na.rm = TRUE))
  }
}

mix_likelihood <- function(data, model, error_rate, by_site = FALSE) {
  assert_that(
    all(c("bac", "dp") %in% colnames(data)),
    nrow(data) == length(model)
  )

  with(data, binom_likelihood(bac, dp, model, error_rate, by_site = by_site))
}

mix_model <- function(gts, mix_prop) {
  assert_that(
    is.matrix(gts),
    is_proportion(mix_prop),
    ncol(gts) == length(mix_prop)
  )

  force_to_interval(colSums(t(gts) * mix_prop), min_val = 0, max_val = 1)
}

#' @importFrom phangorn Children
parent_and_children <- function(phylo, node) {
  unlist(c(treeio::parent(phylo, node), unlist(Children(phylo, node))))
}

node_dist_range <- function(node, phylo, node_dist, max_dist,
                            inclusive = TRUE,
                            nodes_exclude = integer()) {
  nodes <- which(node_dist[node, ] <= max_dist)

  if (inclusive) {
    nodes <- union(nodes, parent_and_children(phylo, node))
  }

  nodes <- setdiff(nodes, union(node, nodes_exclude))

  dist_rem <- max_dist - node_dist[node, nodes]

  if (inclusive && any(dist_rem > 0)) {
    nodes <-
      map(
        which(dist_rem > 0),
        function(i) {
          node_dist_range(
            node = nodes[i],
            phylo = phylo,
            node_dist = node_dist,
            max_dist = dist_rem[i],
            inclusive = inclusive,
            nodes_exclude = c(nodes, node, nodes_exclude)
          )
        }
      ) %>%
      unlist() %>%
      union(., nodes)
  }
  return(sort(nodes))
}


#' @importFrom dplyr if_else
flip_at <- function(x, at, states = c(0L, 1L)) {
  x %>% replace(at, if_else(x[at] == states[1], states[2], states[1]))
}

# all pairwise interger combinations
comb2_int <- function(x, names = c("a", "b")) {
  if (length(x) == 1) {
    n <- x
  } else {
    n <- length(x)
  }
  n <- as.integer(n)
  assert_that(n > 1)
  dplyr::tibble(a = rep(seq_len(n - 1), seq.int(n - 1, 1))) %>%
    dplyr::mutate(
      b = {
        i <- seq_along(a) + 1
        o <- c(0, cumsum((n - 2):1))
        (i - o[a]) %>% as.integer()
      }
    ) %>%
    {
      if (length(x) == 1) {
        .
      } else {
        dplyr::mutate(., a = x[a], b = x[b])
      }
    } %>%
    magrittr::set_colnames(names)
}

log_add <- function(x, y) {
  .max <- pmax(x, y)
  .min <- pmin(x, y)
  .max + log1p(exp(.min - .max))
}

#' @importFrom purrr reduce
log_sum <- function(x) {
  reduce(sort(x, T), function(x, y) x + log1p(exp(y - x)))
}

insert_cols <- function(m1, m2, at) {
  assert_that(
    is.matrix(m1),
    (is.matrix(m2) && nrow(m1) == nrow(m2)) || (is.vector(m2) && nrow(m1) == length(m2)),
    nrow(m1) == nrow(m2),
    is_scalar_integerish(at) && at > 0 && at <= ncol(m1) + 1
  )

  cbind(
    m1[, seq_len(at - 1)],
    `if`(is.matrix(m2), m2, matrix(m2)),
    m1[, at + seq_len(ncol(m1) - at + 1) - 1]
  )
}


#' @importFrom phangorn Ancestors Descendants Children
#' @importFrom purrr map
exclusions <- function(nodes,
                       phylo,
                       node_dist,
                       exclude_parent = FALSE,
                       exclude_child = FALSE,
                       exclude_ancestor = FALSE,
                       exclude_descendant = FALSE,
                       exclude_distance = NULL,
                       exclude_inner = FALSE,
                       exclude_self = TRUE,
                       exclude_root = TRUE) {
  exclude <-
    c(
      `if`(exclude_parent, treeio::parent(phylo, nodes), integer()),
      `if`(exclude_child, unlist(Children(phylo, nodes)), integer()),
      `if`(exclude_ancestor, unlist(Ancestors(phylo, nodes)), integer()),
      `if`(exclude_descendant, unlist(Descendants(phylo, nodes)), integer()),
      `if`(exclude_root, treeio::rootnode(phylo), integer()),
      `if`(exclude_self, nodes, integer()),
      `if`(exclude_inner, inner_nodes(phylo), integer()),
      `if`(
        !is.null(exclude_distance),
        unlist(map(nodes, ~ which(node_dist[., ] <= exclude_distance))),
        integer()
      )
    ) %>%
    sort() %>%
    unique()

  return(exclude)
}

hash_set <- function(.keys = character()) {
  # hashmap::hashmap(as.character(.keys), rep(TRUE, length(.keys)))
  env <- new.env(hash = TRUE)
  for (k in .keys) {
    env[[k]] <- TRUE
  }
  env
}

hs_add <- function(env, keys) {
  # hs[[keys]] <- rep(TRUE, length(keys))
  # hs
  for (k in keys) {
    env[[k]] <- TRUE
  }
  env
}

#' @importFrom tidyr replace_na
#' @importFrom purrr map_lgl
hs_contains <- function(env, keys) {
  # assert_that(inherits(hs, 'Rcpp_Hashmap'))
  # replace_na(hs[[keys]], FALSE)
  keys %in% names(env)
}

combinations <- function(x, zero_len_as_na = TRUE) {
  n <- length(x)
  map(seq_len(2^n) - 1, function(y) {
    v <- x[head(as.logical(intToBits(y)), n)]
    if (zero_len_as_na && length(v) == 0) {
      as(NA, class(x[1]))
    } else {
      v
    }
  })
}



pdist <- function(mat, threads = 1) {

  dist_func_ptr <- RcppXPtrUtils::cppXPtr(
    "double customDist(const arma::mat &A, const arma::mat &B) {
        uint64_t d = 0;
        for (arma::uword idx = 0; idx < A.size(); ++idx) {
          if (A.at(idx) != B.at(idx) && std::isfinite(A.at(idx)) && std::isfinite(B.at(idx))) {
                ++d;
          }
        }
        return d;
     }",
    depends = c("RcppArmadillo"))

  parallelDist::parDist(mat, method="custom", func = dist_func_ptr, threads = threads)
}


