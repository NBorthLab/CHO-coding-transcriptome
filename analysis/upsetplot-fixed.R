UpSettingBug <- function(
  m,
  comb_col = "black",
  pt_size = unit(3, "mm"),
  lwd = 2,
  bg_col = "#F0F0F0",
  bg_pt_col = "#CCCCCC",
  set_order = order(set_size(m), decreasing = TRUE),
  comb_order = if (attr(m, "param")$set_on_rows) {
  order.comb_mat(m[set_order, ], decreasing = TRUE)
  } else {
  order.comb_mat(m[, set_order], decreasing = TRUE)
  },
  top_annotation = upset_top_annotation(m),
  right_annotation = upset_right_annotation(m),
  left_annotation = NULL,
  row_names_side = "left",
  ...
) {
  param <- attr(m, "param")
  set_on_rows <- param$set_on_rows
  mode <- param$mode

  m2 <- m

  class(m2) <- "matrix"

  pt_size <- pt_size
  lwd <- lwd

  if (!is.null(left_annotation)) {
    if (missing(right_annotation)) {
      right_annotation <- NULL
    }
    if (missing(row_names_side)) {
      row_names_side <- "right"
    }
  }

  if(length(bg_col) == 1) bg_col <- c(bg_col, "white")

  n_set <- length(set_size(m))
  bg_col <- rep(bg_col, times = n_set)
  bg_col <- bg_col[seq_len(n_set)]

  if(set_on_rows) {
    n_comb = ncol(m)
    if(length(comb_col == 1)) comb_col = rep(comb_col, n_comb)

    layer_fun = function(j, i, x, y, w, h, fill) {
      nr = round(1/as.numeric(h[1]))
      nc = round(1/as.numeric(w[1]))
      subm = matrix(pindex(m2, i, j), nrow = nr, byrow = FALSE)
      for(k in seq_len(nr)) {
        grid.rect(y = k/nr, height = 1/nr, just = "top",
          gp = gpar(fill = bg_col[k], col = NA))
      }
      grid.points(x, y, size = pt_size, pch = 16,
        gp = gpar(col = ifelse(pindex(m2, i, j), comb_col[j], bg_pt_col)))
      jj = unique(j)
      for(k in seq_len(nc)) {
        if(sum(subm[, k]) >= 2) {
          i_min = min(which(subm[, k] > 0))
          i_max = max(which(subm[, k] > 0))
          grid.lines(c(k - 0.5, k - 0.5)/nc, (nr - c(i_min, i_max) + 0.5)/nr, gp = gpar(col = comb_col[jj[k]], lwd = lwd))
        }
      }
    }

    # check top annotation
    # if it is specified by upset_top_annotation and gp(col) is not set
    ra = top_annotation
    if(length(ra) == 1) {
      ta_call = substitute(top_annotation)
      ta_call = as.list(ta_call)
      if("upset_top_annotation" %in% as.character(ta_call[[1]])) {
        if(!"gp" %in% names(as.list(ta_call))) {
          ra@anno_list[[1]]@fun@var_env$gp$fill = comb_col
          ra@anno_list[[1]]@fun@var_env$gp$col = comb_col
        }
      }
    }

    ht = Heatmap(m2, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(type = "none"),
      layer_fun = layer_fun, show_heatmap_legend = FALSE,
      top_annotation = ra,
      right_annotation = right_annotation, left_annotation = left_annotation,
      row_names_side = row_names_side, col = c("0" = bg_pt_col, "1" = comb_col[1]),
      row_order = set_order, column_order = comb_order, ...)
  } else {
    n_comb = nrow(m)
    if(length(comb_col == 1)) comb_col = rep(comb_col, n_comb)

    layer_fun = function(j, i, x, y, w, h, fill) {
      nr = round(1/as.numeric(h[1]))
      nc = round(1/as.numeric(w[1]))
      subm = matrix(pindex(m2, i, j), nrow = nr, byrow = FALSE)
      for(k in seq_len(nc)) {
        grid.rect(x = k/nc, width = 1/nc, just = "right", gp = gpar(fill = bg_col[k], col = NA))
      }
      grid.points(x, y, size = pt_size, pch = 16, gp = gpar(col = ifelse(pindex(m2, i, j), comb_col[i], "#CCCCCC")))
      ii = unique(i)
      for(k in seq_len(nr)) {
        if(sum(subm[k, ]) >= 2) {
          i_min = min(which(subm[k, ] > 0))
          i_max = max(which(subm[k, ] > 0))
          grid.lines((c(i_min, i_max) - 0.5)/nc, (nr - c(k ,k) + 0.5)/nr, gp = gpar(col = comb_col[ii[k]], lwd = lwd))
        }
      }
    }

    ra = right_annotation
    if(length(ra) == 1) {
      ta_call = substitute(right_annotation)
      ta_call = as.list(ta_call)
      if(as.character(ta_call[[1]]) == "upset_right_annotation") {
        if(!"gp" %in% names(as.list(ta_call))) {
          ra@anno_list[[1]]@fun@var_env$gp$fill = comb_col
          ra@anno_list[[1]]@fun@var_env$gp$col = comb_col
        }
      }
    }
    la = left_annotation
    if(length(la) == 1) {
      ta_call = substitute(left_annotation)
      ta_call = as.list(ta_call)
      if(as.character(ta_call[[1]]) == "upset_left_annotation") {
        if(!"gp" %in% names(as.list(ta_call))) {
          la@anno_list[[1]]@fun@var_env$gp$fill = comb_col
          la@anno_list[[1]]@fun@var_env$gp$col = comb_col
        }
      }
    }
    ht = Heatmap(m2, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(type = "none"),
      layer_fun = layer_fun, show_heatmap_legend = FALSE,
      top_annotation = top_annotation,
      right_annotation = ra, left_annotation = la, col = c("0" = bg_pt_col, "1" = comb_col[1]),
      row_order = comb_order, column_order = set_order, ...)
  }
  ht@heatmap_param$type = "UpSet"
  attr(ht, "UpSet") = TRUE
  ht
}
