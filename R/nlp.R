#' Function to determine the depth to bedrock based on litholog intervals that
#' are labelled as either 'Bedrock' or 'Surficial'
#'
#' This function finds the first occurrence of 'Bedrock' that is below any other
#' intervals labelled as 'Surficial' in each log. As such, this represents a
#' maximum estimate of bedrock depth; it avoids the problems with glaciotectonic
#' rafts of bedrock, but can potentially overestimate bedrock depth if some
#' intervals are misclassified as surficial with the bedrock strata.
#'
#' @param lithologs tibble of litholog data containing the column '.pred_class'
#'   which as two factor levels, 'Bedrock' and 'Surficial'.
#' @param option, either c("last", "first")
#'
#' @return tibble containing the bedrock depths per well, with 'gicwellid'
#'   and '.bedrock_dep' columns.
#' @export
pick_bedrock <- function(lithologs, option = c("last", "first")) {
  option <- match.arg(option)

  if (option == "last") {
    # get the maximum depth of any units that are classified as surficial
    max_surf <- lithologs |>
      filter(.pred_class == "Surficial") |>
      group_by(gicwellid) |>
      summarise(minv = max(int_top_dep, na.rm = TRUE))

    # take the top of the next interval beneath any surficial as the bedrock top
    ypred <- left_join(lithologs, max_surf, by = join_by("gicwellid"))

    ypred <- ypred |>
      group_by(gicwellid) |>
      filter(
        (int_top_dep > minv | is.na(minv)),
        .pred_class == "Bedrock"
      ) |>
      slice_head(n = 1) |>
      ungroup()
  } else {
    ypred <- lithologs |>
      filter(.pred_class == "Bedrock") |>
      group_by(gicwellid) |>
      slice_head(n = 1) |>
      ungroup()
  }

  ypred <- ypred |>
    rename(.bedrock_dep = "int_top_dep") |>
    select(gicwellid, .bedrock_dep)

  return(ypred)
}


add_projected_coords <- function(obj) {
  obj <- obj |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326", keepgeom = TRUE) |>
    terra::project("epsg:3402")

  obj$x <- terra::crds(obj)[, 1]
  obj$y <- terra::crds(obj)[, 2]
  obj <- obj |>
    as.data.frame() |>
    dplyr::as_tibble()
  return(obj)
}

#' Snap bedrock picks to the log intervals and assign a unit column to either
#' 'Bedrock' or 'Surficial'
#'
#' Bedrock top picks, often made in 'Viewlog' are not aligned precisely with the
#' litholog intervals; they often float slightly above or below the interval
#' boundaries. This function finds the closest interval boundary (top or bottom)
#' depth to the pick and moves the pick to align precisely with the boundary.
#'
#' @param lithologs a data.frame/tibble of litholog intervals with 'int_top_dep'
#'   representing the interval top depths, and 'bedrock_dep' defining the picked
#'   bedrock depth in the log.
#'
#' @return a tibble with the modified bedrock depths
#' @export
lithologs_snap <- function(lithologs) {
  snapped_bedrock <- lithologs |>
    dplyr::group_by(.data$gicwellid) |>
    dplyr::mutate(diff = abs(.data$int_top_dep - .data$bedrock_dep)) |>
    dplyr::filter(diff == min(.data$diff)) |>
    dplyr::slice_head() |>
    dplyr::ungroup() |>
    dplyr::mutate(bedrock_dep_snapped = dplyr::if_else(
      abs(.data$int_top_dep - .data$bedrock_dep) <
        abs(.data$int_bot_dep - .data$bedrock_dep),
      .data$int_top_dep,
      .data$int_bot_dep
    )) |>
    dplyr::select("gicwellid", "bedrock_dep_snapped")

  lithologs <- lithologs |>
    dplyr::left_join(snapped_bedrock) |>
    dplyr::mutate(
      bedrock_dep = as.numeric(.data$bedrock_dep),
      bedrock_dep_snapped = as.numeric(.data$bedrock_dep_snapped),
      bedrock_dep = dplyr::if_else(
        !is.na(.data$bedrock_dep_snapped),
        .data$bedrock_dep_snapped,
        .data$bedrock_dep
      )
    ) |>
    dplyr::select(-"bedrock_dep_snapped")

  return(lithologs)
}

#' Assign litholog intervals to 'Bedrock' or 'Surficial' units
#'
#' @param lithologs tibble of litholog data with 'int_top_dep' and 'bedrock_dep'
#'   columns
#'
#' @return tibble with a new 'unit' column containing a factor with the levels =
#'   c("Bedrock", "Surficial")
#' @export
lithologs_assign <- function(lithologs) {
  lithologs |>
    dplyr::mutate(
      unit = dplyr::if_else(
        .data$int_top_dep >= .data$bedrock_dep,
        "Bedrock",
        "Surficial"
      ),
      unit = factor(.data$unit, levels = c("Bedrock", "Surficial"))
    )
}

#' Join the bedrock picks with the lithologs
#'
#' The ‘gicwellid’ field in the lithologs will represent ‘LOC_NAME_ALT2’ for
#' logs coming from ABGEOL, and will represent the GIC_WELL_ID for those coming
#' from AWWID, and also contains UWIs and other local names etc. This was
#' because original GFM pick compilation had a single gicwellid column that used
#' LOC_NAME_ALT2 as the identifier for picks made from ABGEOL, so it is the only
#' column that can be used to relate the picks back to the logs.
#'
#' Unfortunately, this identifier is not entirely unique and there are some
#' duplicated ids.
#'
#' This function joins the picks with the lithologs based on 'gicwellid', but also
#' removes joined locations where the pick coordinates are not within 500 m of
#' the litholog coordinates. This should not remove more than 3000 picks out of
#' ~ 130,000 otherwise an error is raised.
#'
#' @param lithologs a tibble of compiled lithologs
#' @param picks a tibble of compiled picks
#'
#' @return the joined data
#' @export
lithologs_join <- function(lithologs, picks) {
  # convert coordinates
  lithologs <- add_projected_coords(lithologs)
  picks <- add_projected_coords(picks)

  # join logs with picks
  picks_df <- picks |>
    dplyr::select(c("gicwellid", "x", "y", "bedrock_dep")) |>
    dplyr::rename(picks_x = "x", picks_y = "y")

  lithologs_labelled <- dplyr::left_join(
    lithologs,
    picks_df,
    by = dplyr::join_by("gicwellid")
  )

  # set joined columns to NA if the coordinates do not match within tolerance
  lithologs_labelled <- lithologs_labelled |>
    dplyr::mutate(
      dist = sqrt((.data$x - .data$picks_x)^2 + (.data$y - .data$picks_y)^2),
      bedrock_dep = dplyr::if_else(
        .data$dist < 1000,
        .data$bedrock_dep,
        NA_real_
      )
    ) |>
    dplyr::select(-c("picks_x", "picks_y", "x", "y", "dist"))

  return(lithologs_labelled)
}

#' Preprocessing recipe for NLP model
#'
#' @param data tibble of training data. Needs to be grouped by the 'gicwellid'
#'   column
#'
#' @return a recipe object
#' @export
recipe_nlp <- function(data) {
  rec <- data |>
    recipes::recipe(unit ~ .) |>
    recipes::update_role(
      dplyr::matches("gicwellid"),
      old_role = "predictor",
      new_role = "identifier"
    ) |>
    recipes::step_mutate_at(
      recipes::all_logical_predictors(),
      fn = as.integer
    )

  # add indicator variables for missing data
  rec <- rec |>
    recipes::step_indicate_na(material, material_desc, colour)

  # set infinite values to nan
  rec <- rec |>
    recipes::step_mutate_at(
      recipes::all_integer_predictors(),
      fn = ~ ifelse(is.infinite(.), NA_integer_, .)
    ) |>
    recipes::step_mutate_at(
      recipes::all_numeric_predictors(),
      fn = ~ ifelse(is.infinite(.), NA_real_, .)
    )

  # tokenization
  rec <- rec |>
    textrecipes::step_tokenize(
      material, material_desc,
      options = list(strip_numeric = TRUE, lowercase = TRUE, strip_punct = TRUE)
    ) |>
    textrecipes::step_tokenize(
      colour,
      token = "ngrams",
      options = list(n = 3L, n_min = 1L, lowercase = TRUE)
    ) |>
    textrecipes::step_stopwords(
      material, material_desc,
      stopword_source = "snowball"
    ) |>
    textrecipes::step_tokenmerge(
      material, material_desc,
      prefix = "material"
    ) |>
    textrecipes::step_tf(material) |>
    textrecipes::step_tf(colour) |>
    textrecipes::step_clean_names(starts_with("tf_")) |>
    recipes::step_rm(recipes::has_role("identifier"))

  return(rec)
}

#' Helper function to add additional predictors to each water well
#'
#' @param data tibble of water well data
#'
#' @return tibble of water well data with added terms
#' @export
add_features <- function(data) {
  gravel_regex <- paste(gravel_terms, collapse = "|")
  missing_regex <- paste(nan_terms, collapse = "|")

  # add per interval features
  # binary variables for terms that indicate missing descriptions
  interval_data <- data |>
    lazy_dt() |>
    mutate(
      material = tolower(material),
      material_desc = tolower(material_desc),
      colour = tolower(colour)
    )

  # text-pattern matching terms
  interval_data <- interval_data |>
    mutate(
      matches_gravel = pmax(
        as.integer(str_detect(material, !!gravel_regex)),
        as.integer(str_detect(material_desc, !!gravel_regex)),
        na.rm = TRUE
      ),
      matches_till = pmax(
        as.integer(str_detect(material, "till|drift|glacial")),
        as.integer(str_detect(material_desc, "till|drift|glacial")),
        na.rm = TRUE
      )
    )

  # add grouped by well features features
  # borehole descriptive features
  grouped_data <- interval_data |>
    group_by(gicwellid) |>
    mutate(
      int_thk = int_bot_dep - int_top_dep,
      n_intervals = n(),
      max_depth = max(int_bot_dep, na.rm = TRUE),
      rel_depth = int_top_dep / max(int_bot_dep, na.rm = TRUE),
      int_thk_rel = int_thk / max_depth,
      waterbearing = as.integer(waterbearing)
    )

  # proportions of text-pattern matching terms over sliding windows
  grouped_data <- grouped_data |>
    mutate(
      # prop of gravel-like terms beneath current interval
      gr_w1 = slide_dbl(matches_gravel, .before = 1, .after = 1, .f = mean),
      gr_w2 = slide_dbl(matches_gravel, .before = 2, .after = 2, .f = mean),
      gr_w3 = slide_dbl(matches_gravel, .before = 3, .after = 3, .f = mean),
      gr_w5 = slide_dbl(matches_gravel, .before = 5, .after = 5, .f = mean),

      # prop of till-like terms beneath current interval
      till_w1 = slide_dbl(matches_till, .before = 1, .after = 1, .f = mean),
      till_w2 = slide_dbl(matches_till, .before = 2, .after = 2, .f = mean),
      till_w3 = slide_dbl(matches_till, .before = 3, .after = 3, .f = mean),
      till_w5 = slide_dbl(matches_till, .before = 5, .after = 5, .f = mean)
    )

  # occurrences of text-pattern matching terms within each borehole log
  grouped_data <- grouped_data |>
    mutate(
      max_grvl = max(int_top_dep[matches_gravel == 1]),
      max_till = max(int_top_dep[matches_till == 1])
    )

  # return ungrouped data
  grouped_data |>
    ungroup() |>
    as_tibble() |>
    select(-contains("matches"))
}

#' Terms that indicate missing descriptions
#'
#' @export
nan_terms <- c(
  "unknown", "fill", "undetermined", "no recovery", "other",
  "no data", "predrilled", "old", "well", "backfill",
  "see comments", "unspecified", "unreadable", "lost",
  "synthetic", "comments", "notes", "^formation$"
)

#' Common gravel unit terms
#'
#' @export
gravel_terms <- c(
  "gravel",
  "pebbl",
  "rocks",
  "\\bstones\\b",
  "boulder",
  "cobbl",
  "conglomerate",
  "chert"
)
