

library(tidyverse)
library(sp)
library(sf)
library(data.table)

`%!in%` <- compose(`!`, `%in%`)

load_data_hook <- function(raw_path = NULL,
                           mesh_path = "./data/",
                           show_plot = F, done = T){
  
  if(done){
    cell_data= read.csv(paste0(mesh_path, "2D_hook_data.csv"))%>%
      select(-X)
    cell_data$id_cell[cell_data$sampl == "S1" & cell_data$pool == "P3" & 
                        cell_data$Gen == "Col0" & cell_data$id_cell == 68] = 9000
    cell_data$id_cell[cell_data$sampl == "S1" & cell_data$pool == "P3" & 
                        cell_data$Gen == "Col0" & cell_data$id_cell == 194] = 9001.5
  }else{
  # Generate 2D cell data from MorphoGrapX mesh
  path = raw_path # "./data/raw/"
  fls = list.files(raw_path)
  fls = fls[grepl(".csv", fls)]
  
  hook = list()
  for (fl in fls){
    print(fl)
    df = load_cell_data(paste0(path,fl))
    
    df = df %>% 
      filter(y < 1000, id_cell > 0) %>%
      mutate(id_point = paste0(round(x), ";", round(y))) %>%
      filter(!duplicated(id_point))

    if(nrow(df)>3){
      df_hook = get_cell_polygons(df, file_name = fl)

      hook[[length(hook) + 1]] <- df_hook 
    }
  }
  hook_tibble <- as_tibble(bind_rows(hook))
  
  mba = hook_tibble%>%
    filter(id_cell != -1)%>%
    group_by(id_cell, Gen, sampl, pool)%>%
    distinct(basal_apical, euc, atan)%>%
    summarise(n = n(),
              mba = round(mean(basal_apical)),
              m_euc = mean(euc),
              m_atan = mean(atan), .groups = "drop")

  hook_tibble <- hook_tibble %>% 
    left_join(mba %>%select(-n),by = c("id_cell", "Gen", "sampl", "pool") )
  
  # 1. compute per-cell statistics with center and angle
  cell_rate <- hook_tibble %>%
    dplyr::group_by(Gen, pool, sampl, time, id_cell) %>%
    dplyr::summarise(
      mx = mean(x), my = mean(y),
      var_area = mean(area),
      var_major_l = mean(major_axis_length),
      var_minor_l = mean(minor_axis_length),
      var_circ = mean(circularity),
      .groups = "drop"
    )
  
  cell_rate <- cell_rate %>%
    arrange(Gen, pool, sampl, id_cell, time) %>%
    group_by(Gen, pool, sampl, id_cell) %>%
    mutate(
      across(starts_with("var"),
             ~ .x - lag(.x),
             .names = "diff_{.col}")
    ) %>%
    mutate(
      base_area    = ifelse(any(time == 0), var_area[time == 0][1], NA_real_),
      base_major_l = ifelse(any(time == 0), var_major_l[time == 0][1], NA_real_),
      base_minor_l = ifelse(any(time == 0), var_minor_l[time == 0][1], NA_real_),
      base_circ    = ifelse(any(time == 0), var_circ[time == 0][1], NA_real_)
    ) %>%
    mutate(
      norm_diff_var_area    = diff_var_area    / base_area,
      norm_diff_var_major_l = diff_var_major_l / base_major_l,
      norm_diff_var_minor_l = diff_var_minor_l / base_minor_l,
      norm_diff_var_circ    = diff_var_circ    / base_circ
    ) %>%
    ungroup()
  
  # 3. attach back to raw tibble
  cell_data <- hook_tibble %>%
    left_join(cell_rate%>%filter(time == 8)%>%select(-time),
              by = c("Gen","pool","sampl","id_cell"))
  }
  
  if(show_plot){
    cell_data %>%
      filter(id_cell != -1)%>%
      ggplot(aes(x,y))+
      geom_polygon(aes(fill = diff_var_area, group = factor(id_cell)))+
      coord_fixed()+
      facet_grid(Gen*time~sampl*pool)+
      theme_dark()+
      viridis::scale_fill_viridis(option = "H")+
      labs(fill = "Cell area\ngrowth [µm²]")
    }# Example
  
  return(cell_data)
}

make_corresp_matrix <- function(cell_data){
  cell_data_ordered <- ordering_cells(cell_data%>%filter(id_cell != -1, time == 0))
  
  outside = cell_data%>%filter(id_cell == -1)%>%mutate(cut = NA, inner_outer = NA, basal_apical2 = NA)
  cell_data_ordered = rbind(cell_data_ordered, 
                            outside%>%mutate(id_plant = paste0(Gen, pool, sampl)))
  
  cell_data_binned <- cell_data_ordered %>%
    filter(time == 0) %>%
    mutate(
      finner_outer = floor(inner_outer),
      finner_outer = ifelse(finner_outer<1,1, ifelse(finner_outer>8,8, finner_outer)),
      fbasal_apical2 = floor(basal_apical2),
      fbasal_apical2 = ifelse(fbasal_apical2<1,1, ifelse(fbasal_apical2>8,8, fbasal_apical2)),
      inner_bin = ifelse(finner_outer<4,4, finner_outer ),
      inner_bin = ifelse(inner_bin>=5,"outer", "inner" ),
      basal_bin = ifelse(fbasal_apical2<4,4, fbasal_apical2 ),
      basal_bin = ifelse(basal_bin>=5,"apical", "basal" ),
      bin = paste0(inner_bin, "-",basal_bin),
      g = paste0(Gen,bin))
  
  corr_matrix = cell_data_binned %>%
    group_by(id_plant,time, id_cell)%>%
    mutate(
      x1 = x,
      y1 = y,
      x2 = lead(x, default = first(x)),  # connect last to first
      y2 = lead(y, default = first(y))
    ) %>%
    ungroup()
  
  corr_matrix_dcr = corr_matrix%>%filter(Gen == "dcr")
  dcr_prop = corr_matrix_dcr %>%
    dplyr::group_by(id_cell, finner_outer, fbasal_apical2)%>%
    dplyr::distinct(diff_var_major_l, diff_var_minor_l,
                    major_axis_length, minor_axis_length)%>%
    dplyr::group_by(finner_outer, fbasal_apical2)%>%
    dplyr::summarise(L_dcr = mean(diff_var_major_l/major_axis_length),
                     l_dcr = mean(diff_var_minor_l/minor_axis_length), .groups = "drop")
  
  
  corr_matrix_cyp = corr_matrix%>%filter(Gen == "cyp")
  cyp_prop = corr_matrix_cyp %>%
    dplyr::group_by(id_cell, finner_outer, fbasal_apical2)%>%
    dplyr::distinct(diff_var_major_l, diff_var_minor_l,
                    major_axis_length, minor_axis_length)%>%
    dplyr::group_by(finner_outer, fbasal_apical2)%>%
    dplyr::summarise(L_cyp = mean(diff_var_major_l/major_axis_length),
                     l_cyp = mean(diff_var_minor_l/minor_axis_length), .groups = "drop")
  
  corr_matrix_col = corr_matrix%>%filter(Gen == "Col0")
  col_prop = corr_matrix_col %>%
    dplyr::group_by(id_cell, finner_outer, fbasal_apical2)%>%
    dplyr::distinct(diff_var_major_l, diff_var_minor_l,
                    major_axis_length, minor_axis_length)%>%
    dplyr::group_by(finner_outer, fbasal_apical2)%>%
    dplyr::summarise(L_col = mean(diff_var_major_l/major_axis_length),
                     l_col = mean(diff_var_minor_l/minor_axis_length), .groups = "drop")
  
  corr_matrix <- corr_matrix %>% 
    left_join(col_prop, by = c("finner_outer", "fbasal_apical2")) %>% 
    left_join(dcr_prop, by = c("finner_outer", "fbasal_apical2")) %>% 
    left_join(cyp_prop, by = c("finner_outer", "fbasal_apical2")) 
  
  return(corr_matrix)
}


quadrant_stat_analysis <- function(corr_matrix){
  # Stat t-test on strain 
  corr_matrix = corr_matrix%>%filter(bin != "NA-NA")%>%
    distinct(norm_diff_var_area,norm_diff_var_major_l, Gen, bin, id_cell, basal_bin, inner_bin)
  
  lm_strain <- lm(norm_diff_var_area ~ Gen*bin, 
                  data = corr_matrix)
  lm_elong <- lm(norm_diff_var_major_l ~ Gen*bin, 
                 data = corr_matrix)
  
  lm_elong <- lm(norm_diff_var_major_l ~ bin, 
                 data = corr_matrix%>%filter(Gen == "Col0"))
  
  emm_strain = emmeans::emmeans(lm_strain, ~ Gen * bin)
  emm_elonga = emmeans::emmeans(lm_elong, ~ Gen*bin)
  
  emm_elong_df = as.data.frame(emm_elonga)
  
  emm_df <- as.data.frame(emm_strain)
  
  ggplot() +
    geom_boxplot(aes(x = bin, norm_diff_var_area, group = bin, fill = Gen), alpha = 0.5,
                 data =corr_matrix)+
    theme_bw()+
    ylim(-1,4)
  
  cld_elong <- multcomp::cld(emm_elonga, Letters_elong = letters, adjust = "sidak") %>%
    as.data.frame() %>%
    rename(group_label = .group)
  corrmat_elong <- corr_matrix%>%
    left_join(cld_elong, by = c("bin", "Gen"))%>%
    mutate(group = bin)
  
  pl = ggplot(corrmat_elong, aes(x = inner_bin, y = norm_diff_var_major_l, fill = emmean)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_text(    aes(x = inner_bin, y = 2, label = group_label),
      vjust = 0,
      fontface = "bold"
    ) +
    viridis::scale_fill_viridis(option = "turbo") +
    labs(
      x = "Genotype × bin",
      y = "Strain",
      title = "genotype × bin combinations"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    )+
    facet_wrap(basal_bin~Gen)
  
  pl = ggplot(corrmat_elong%>%filter(Gen == "Col0"), aes(x = inner_bin, y = norm_diff_var_major_l, fill = emmean)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_text(    aes(x = inner_bin, y = 2, label = group_label),
                  vjust = 0,
                  fontface = "bold"
    ) +
    viridis::scale_fill_viridis(option = "turbo") +
    labs(
      x = "Genotype × bin",
      y = "Strain",
      title = "genotype × bin combinations"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    )+
    facet_wrap(basal_bin~Gen)
  pl
  ggsave("./data/out/fig/Gen_bin_elong_boxplot.svg")
  
  corr_matrix = corr_matrix %>% left_join(emm_df, by = c("Gen", "bin"))
  
  print(corr_matrix %>% 
    ggplot() +
    geom_boxplot(aes(x = inner_bin, norm_diff_var_area, fill = emmean))+
    theme_bw()+
    facet_wrap(basal_bin~Gen)+
    labs(y = 'strain [-]', fill = "mean strain [-]")+
    viridis::scale_fill_viridis(option = "turbo"))
  
  res_df <- corr_matrix %>%
    mutate(
      resid_strain = residuals(lm_strain),
      resid_elong = residuals(lm_elong)
    )
  
  car::leveneTest(resid_strain ~ interaction(Gen,bin), data = res_df, center = median)
  
  # simple linear model on absolute residuals
  lm_resid_strain <- lm(abs(resid_strain) ~ Gen * bin, data = res_df)
  lm_resid_elong <- lm(abs(resid_elong) ~ Gen * bin, data = res_df)
  
  emm_res_strain <- emmeans::emmeans(lm_resid_strain, ~ Gen * bin)
  emm_res_elong <- emmeans::emmeans(lm_resid_elong, ~ Gen * bin)
  
  # Compact letter display (grouping letters)
  cld_res_strain <- multcomp::cld(emm_res_strain, Letters_strain = letters, adjust = "sidak") %>%
    as.data.frame() %>%
    rename(group_label = .group)
  
  cld_res_elong <- multcomp::cld(emm_res_elong, Letters_elong = letters, adjust = "sidak") %>%
    as.data.frame() %>%
    rename(group_label = .group)
  
  res_df_strain <- res_df %>%
    left_join(cld_res_strain, by = c("Gen", "bin"))%>%
    mutate(group = interaction(Gen, bin))
  
  res_df_elong <- res_df %>%
    left_join(cld_res_elong, by = c("Gen", "bin"))%>%
    mutate(group = interaction(Gen, bin))
  
  ggplot(res_df_strain, aes(x = inner_bin, y = abs(resid_strain), fill = Gen)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
    geom_text(
      aes(x = inner_bin, y = 2, label = group_label),
      vjust = 0,
      fontface = "bold"
    ) +
    theme_minimal(base_size = 14) +
    scale_fill_viridis_d() +
    labs(
      x = "Genotype × bin",
      y = "absolute residuals",
      title = "Variance amplitude across genotype × bin combinations"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    )+
    facet_wrap(basal_bin~Gen)+
    ylim(0,2.5)
  
  ggplot(res_df_elong, aes(x = inner_bin, y = abs(resid_strain), fill = Gen)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
    geom_text(
      aes(x = inner_bin, y = 2, label = group_label),
      vjust = 0,
      fontface = "bold"
    ) +
    theme_minimal(base_size = 14) +
    scale_fill_viridis_d() +
    labs(
      x = "Genotype × bin",
      y = "absolute residuals",
      title = "Variance amplitude across genotype × bin combinations"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    )+
    facet_wrap(basal_bin~Gen)+
    ylim(0,2.5)
  
}

prep_geo <- function(corr_matrix, filtering = "Col0P3S1"){
  data = corr_matrix%>%
    filter(id_plant == filtering)%>%
    mutate(id_cell = ifelse(id_cell == -1, max(id_cell)+1, id_cell))
    
  return(data)
  
}
  


load_cell_data <- function(path){
  
  if(!is.na(str_match(path, ".txt")[1])){
    df = read.delim(path, header = F, sep = " ", fill = T)
    df = df[-1,2:5]
    colnames(df) <- c("x", "y", "z", "id_cell")
    df$idx = 1:nrow(df)
    
    bef = which(is.na(df$id_cell))[1]
    df = df %>%
      filter(idx < bef)
    write.csv(df, file = str_replace_all(path, ".txt", ".csv"))
    file.remove(path)
  }else if(!is.na(str_match(path, ".csv")[1])){
    df = read.csv(path)
  }else{
    df = as_tibble(matrix(c(0), ncol = 5))
    colnames(df) <- c("x", "y", "z", "id_cell", "idx")
    df = df%>%filter(idx != 0 )
  }

  return(df)
}
  

get_cell_polygons <- function(df_hook, verbatim = FALSE, alpha = 1.5, file_name = "P1_Col0_S1_T0.csv") {

  
  name = unlist(str_split(file_name, "_"))
  df_hook$pool = name[1]
  df_hook$Gen = name[2]
  df_hook$sampl = name[3]
  df_hook$time = parse_number(name[4])
  
  df_hook = df_hook%>%
    mutate(x=ifelse((pool == "P2" & Gen == "Col0" & sampl == "S2") |
                      (pool == "P3" & Gen == "Col0" & sampl == "S1") |
                      (pool == "P1" & Gen == "cyp" & sampl == "S1") |
                      (pool == "P2" & Gen == "cyp" & sampl == "S1") |
                      (pool == "P2" & Gen == "cyp" & sampl == "S2") |
                      (pool == "P1" & Gen == "dcr" & sampl == "S1") | 
                      (pool == "P3" & Gen == "Col0" & sampl == "S2"), -x,x)
    )
  
  
  id_cell_vector <- sort(unique(df_hook$id_cell))
  id_cell_vector <- id_cell_vector[id_cell_vector > 0]
  
  cell_polygons <- list()
  multi_polygons <- list()
  multi_ids = NULL
  for (i in id_cell_vector) {
    
    if(i == id_cell_vector[round(length(id_cell_vector)/2)]){
      print("50% done for this mesh...")
    }
    tmp <- df_hook %>%
      filter(id_cell == i) %>%
      mutate(id_point = paste0(round(x, 3), ";", round(y, 3))) %>%
      filter(!duplicated(id_point))
    
    if (nrow(tmp) > 3) {
      ## try alphashape
      my.ashape <- try(alphahull::ashape(x = tmp$x, y = tmp$y, alpha = alpha), silent = TRUE)
      while (!inherits(my.ashape, "ashape")) {
        tmp <- tmp %>% mutate(x = jitter(x, 0.004), y = jitter(y, 0.004))
        my.ashape <- try(alphahull::ashape(x = tmp$x, y = tmp$y, alpha = alpha), silent = TRUE)
      }
      
      # extract lines
      a <- data.frame(my.ashape$edges)[, c("x1", "y1", "x2", "y2")]
      l <- lapply(1:nrow(a), function(j) {
        sf::st_linestring(matrix(as.numeric(a[j, ]), ncol = 2, byrow = TRUE))
      })
      my_multilinestring <- sf::st_sfc(l, crs = 2056) |> sf::st_sf()
      
      if (sf::st_is_empty(sf::st_union(my_multilinestring) |> sf::st_polygonize())) {
        sf_points <- sf::st_as_sf(tmp, coords = c("x", "y"), crs = 2056)
        concave_hull <- concaveman::concaveman(sf_points, concavity = 2.5)
        r_poly <- smoothr::smooth(concave_hull, method = "ksmooth", smoothness = 5)
      } else {
        alphapoly <- sf::st_union(my_multilinestring) |> sf::st_polygonize() |> sf::st_collection_extract()
        r_poly <- smoothr::smooth(alphapoly, method = "ksmooth", smoothness = 5)
      }
      
      # keep largest polygon if multiple
      areas <- as.numeric(sf::st_area(r_poly))
      if (length(areas) > 1) {
        r_poly <- r_poly[which.max(areas)]
      }
      
      # measure properties
      coords <- sf::st_coordinates(r_poly)[, 1:2]
      ## simplify for output
      simpl <- sf::st_simplify(r_poly, dTolerance = 0.1, preserveTopology = TRUE)

      pol = pol_to_data(simpl)%>%
        mutate(id_cell = i)
      cell_polygons[[length(cell_polygons) + 1]] <- pol
      multi_ids <- c(multi_ids, i)
      multi_polygons[[length(multi_polygons) + 1]] <- simpl
      
      
    }
  }
  
  # Git the outline of the general topology

  cells_sfc <- do.call(c, multi_polygons)   # combine into an sfc object
  cells_sf  <- cells_sf <- st_sf(
    id_cell = multi_ids,
    geometry = cells_sfc
  )
  
  # Compute centroids
  df_sf_geom <- cells_sf %>%
    mutate(centroid = st_centroid(geometry))
  # Extract coordinates
  centroids <- st_coordinates(df_sf_geom$centroid)
  # Add centroid x/y and mean x/y
  df_sf_geom <- df_sf_geom %>%
    mutate(mx = centroids[,1],
           my = centroids[,2]
    )
  cat("Click on the plot to set a new center, or press ESC to skip.\n")
  
  x11()
  plot(df_sf_geom$mx, df_sf_geom$my)
  new_center = locator(n = 1)
  dev.off()
  
  # --- OPTION 2: TYPE new center manually ---
  if (is.null(new_center)) {
    cx <- as.numeric(readline("Enter new center x: "))
    cy <- as.numeric(readline("Enter new center y: "))
  } else {
    cx <- new_center$x
    cy <- new_center$y
  }
  
  # Shift coordinates
  df_sf_geom <- df_sf_geom %>%
    mutate(
      mx_shifted = mx - cx,
      my_shifted = my - cy
    )
  
  index_base = as_tibble(basal_to_apical(df_sf_geom))
  
  outline <- st_union(cells_sf)
  cells_buf <- st_buffer(cells_sf, 2)  # 
  
  cells_union <- st_union(cells_buf)
  walls = pol_to_data(cells_union)%>%
    mutate(id_cell = -1)
  
  cell_polygons[[length(cell_polygons) + 1]] <- walls
  cells_df <- bind_rows(cell_polygons)
  cells_df = left_join(cells_df, index_base%>%select(id_cell, dist_along,basal_apical,cx,cy, atan, euc), by = c("id_cell"))
  
  cells_df$pool = name[1]
  cells_df$Gen = name[2]
  cells_df$sampl = name[3]
  cells_df$time = parse_number(name[4])
  
  return(cells_df)
}


basal_to_apical <- function(df_sf_geom){


  df <- df_sf_geom %>%
    mutate(euc   = sqrt((mx_shifted)^2 + (my_shifted)^2),
           atan  = atan2(mx_shifted, my_shifted),   # shift 
           atan0 = ifelse(atan > 0, pi+atan, atan)
    )
  
  df %>%
    ggplot(aes(mx_shifted, my_shifted))+
    geom_point(aes(colour = atan))+coord_fixed()
  df = df %>%
    arrange(-atan)
  
  fit <- smooth.spline(df$atan, df$euc, spar = 1.1)
  
  
  # sample finely along the spiral
  theta_seq <- seq(min(df$atan)*0.95, max(df$atan)*0.95, length.out = 2000)
  r_seq <- predict(fit, theta_seq)$y
  
  # back to Cartesian coordinates
  curve_pts <- tibble(
    t = 1:length(theta_seq),
    x = r_seq * sin(theta_seq),
    y = r_seq * cos(theta_seq)
  )
  curve_pts = curve_pts %>%
    mutate(dx = c(0, diff(x)),
           dy = c(0, diff(y)),
           seg_len = sqrt(dx^2 + dy^2),
           arc_len = cumsum(seg_len))
  
  print(ggplot() +
    geom_point(aes(mx_shifted, my_shifted), data = df, color = "black", alpha = 0.6) +
    geom_path(aes(x, y), data = curve_pts, color = "orange", size = 1.2) +
    coord_fixed())
  
  total_len <- max(curve_pts$arc_len)
  
  df <- df %>%
    mutate(
      cx = st_coordinates(centroid)[,1],
      cy = st_coordinates(centroid)[,2],
    ) %>%
    rowwise() %>%
    mutate(
      nearest = {
        dists <- (curve_pts$x - cx)^2 + (curve_pts$y - cy)^2
        which.min(dists)
      },
      dist_along = curve_pts$arc_len[nearest],
      basal_apical = (dist_along / total_len) * 10
    ) %>%
    ungroup()
  return(df)
}

  
pol_to_data <- function(pol, verbatim = F){
  
  coords <- sf::st_coordinates(pol)[, 1:2]
  
  # Perform PCA
  pca <- prcomp(coords, center = TRUE, scale. = FALSE)
  # Extract PCA results
  pc1 <- pca$rotation[, 1]  # Major axis direction
  pc2 <- pca$rotation[, 2]  # Minor axis direction
  # Compute axis lengths (project points onto principal components)
  proj_coords <- as.matrix(coords) %*% pca$rotation
  major_axis_length <- diff(range(proj_coords[, 1]))  # Along PC1
  minor_axis_length <- diff(range(proj_coords[, 2]))  # Along PC2
  if(verbatim){
    # Plot cell polygon with PCA axes
    pl = ggplot() +
      geom_polygon(aes(x = coords[, 1], y = coords[, 2]), fill = "lightblue", color = "blue") +
      geom_segment(aes(x = pca$center[1], y = pca$center[2],
                       xend = pca$center[1] + pc1[1] * major_axis_length / 2,
                       yend = pca$center[2] + pc1[2] * major_axis_length / 2),
                   color = "red", size = 1.2) +
      geom_segment(aes(x = pca$center[1], y = pca$center[2],
                       xend = pca$center[1] + pc2[1] * minor_axis_length / 2,
                       yend = pca$center[2] + pc2[2] * minor_axis_length / 2),
                   color = "blue", size = 1.2) +
      coord_equal()
    print(pl)
  }
  
  perimeter = as.numeric(sf::st_length(sf::st_cast(pol, "MULTILINESTRING")))
  area = as.numeric( sf::st_area(pol))
  circularity = 4*pi*(area/perimeter^2)
  
  
  pol_tibble <- tibble(
    x = coords[, 1], y = coords[, 2],
    major_axis_length,
    minor_axis_length,
    perimeter = perimeter,
    area = area,
    circularity = circularity
  )
  return(pol_tibble)
}


write_geo <- function(data, dim = 2, path_geo, Celldomain =F, make_base = F){
  prefix = data$id_plant[1]
  
  scenario = building_strain_cases(data%>%filter(id_cell != -1))
  write_case(scenario, prefix)
  
  `%!in%` <- compose(`!`, `%in%`)
  date = Sys.time()
  x1 = paste0('// Gmsh project created on ', date,'\nSetFactory("OpenCASCADE");\n//+\n')
  
  if(make_base){
    data_p = prep_base(data, dir_path = "./data/in/growth_data/bc/")
  }
  pl = data_p %>%
    filter(id_cell == -1)%>%
    ggplot(aes(x,y))+
    geom_polygon()+coord_fixed()
  
  print(pl)
  
  data = data_p %>%
    mutate(id_cell = ifelse(id_cell == -1, max(id_cell)+1, id_cell))
  
  data = data %>%filter(!is.na(x))%>%
    group_by(time, id_cell)%>%
    mutate(
      x1 = x,
      y1 = y,
      x2 = lead(x, default = first(x)),  # connect last to first
      y2 = lead(y, default = first(y))
    ) %>%
    ungroup()
  
  if(dim == 2){
    data$z = 0
    data$z1 = 0
    data$z2 = 0
  }
  
  center = c(mean(data$x), mean(data$y))
  central = data %>% filter(id_cell != max(id_cell)) %>% 
    dplyr::group_by(id_cell) %>% 
    dplyr::summarise(mx = mean(x), my = mean(y), .groups = "drop")%>%
    mutate(euc = sqrt((mx-center[1])^2+(my-center[2])^2))
  center_cell_id = central$id_cell[central$euc == min(central$euc)][1] 
  
  k = h = j = 1
  txt = x1
  for(i in sort(unique(data$id_cell))){
    tmp = data%>%filter(id_cell == i)
    i_point = unique_point_id(tmp)
    
    i_point$idx = seq(k,-1+k+nrow(i_point),1)
    i_point$idx2 = c(i_point$idx[-1],i_point$idx[1])
    
    if(i < max(data$id_cell)){
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',1.0};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nSurface(',
                    h,") = {", j,'};\n//+\n', collapse = "")
      if(Celldomain){
        Physical = paste0("Physical Surface(",h,") = {",h,"};\n")
      }else{
        Physical = paste0("//Physical Surface(",h,") = {",h,"};\n")
      }
      
      txt = paste0(txt, dots, Lines, Surf, Physical)
      if(i == center_cell_id){
        fix_curve = paste0('Physical Curve("fix", 2) = {',paste0(i_point$idx, collapse = ", "),'};\n')
        fix_id= i_point$idx
      }
    }else{
      
      all_inner_id = 1:(k-1)
      all_inner_id = all_inner_id[all_inner_id %!in% fix_id]
      Physical_curve = paste0('Physical Curve("inner", 1) = {',paste0( all_inner_id  , collapse = ", "),'};\n')
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',0.5};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nPlane Surface(',
                    h,") = {", paste0(sort(seq(1,j,2), decreasing = T), collapse = ", "),'};\n//+\n', collapse = "")
      Physical = paste0("Physical Surface(0) = {",h,"};\n")
      outer_curve = paste0('Physical Curve("outer", 3) = {',paste0(i_point$idx, collapse = ", "),'};\n')
      txt = paste0(txt, dots, Lines, Surf, Physical, Physical_curve, fix_curve, outer_curve)
    }
    
    h = h+1
    j = j+2
    k = max(i_point$idx)+1
  }
  write(txt, path_geo)
}


#' @title remove duplicated points
#'
#'
#' @param tmp data to filter
#' @keywords GMSH Geo
#' @export
#'
#'

unique_point_id <- function(tmp){
  tmp$id_point = paste0('x:',tmp$x1,'_y:',tmp$y1,'_z:',tmp$z1)
  tmp$id_point2 = paste0('x:',tmp$x2,'_y:',tmp$y2,'_z:',tmp$z2)
  i_point = tibble(id_point = unique(c(tmp$id_point, tmp$id_point2)))
  
  x_coor = unlist(str_split(unlist(str_split(i_point$id_point, pattern = 'x:')),pattern = "_y"))
  y_coor = unlist(str_split(unlist(str_split(i_point$id_point, pattern = 'y:')),pattern = "_z"))
  z_coor = unlist(str_split(i_point$id_point, pattern = 'z:'))
  i_point$x = parse_number(x_coor[seq(2,length(x_coor),3)])
  i_point$y = parse_number(y_coor[seq(2,length(y_coor),3)])
  i_point$z = parse_number(z_coor[seq(2,length(z_coor),2)])
  
  i_point
  return(i_point)
}


read_geo <- function(geo_path, out_geo = NULL){
  geo = read_lines(geo_path)
  point = geo[grepl("Point", geo)]
  point = unlist(str_split(point, pattern = "="))
  point = point[seq(2, length(point), by = 2)]
  pointX = as_data_frame(t(matrix(parse_number(unlist(str_split(point, pattern = ","))), nrow = 4)))
  Lines = geo[grepl("Line", geo)]
  Lines = unlist(str_split(Lines, pattern = "="))
  Lines = Lines[seq(2, length(Lines), by = 2)]
  LinesX = as_data_frame(t(matrix(parse_number(unlist(str_split(Lines, pattern = ","))), nrow = 2)))
  Surf = geo[grepl("Curve Loop", geo)]
  Surf = unlist(str_split(Surf, pattern = "="))
  Surf = Surf[seq(2, length(Surf), by = 2)]
  
  data = NULL
  for(i in 1:length(Surf)){
    i_l = parse_number(unlist(str_split(Surf[i], pattern = ",")))
    l_p = LinesX[abs(i_l), ]
    x_point = tibble(x1 = ifelse(i_l>=0,pointX$V1[l_p$V1],pointX$V1[l_p$V2]),
                     y1 = ifelse(i_l>=0,pointX$V2[l_p$V1],pointX$V2[l_p$V2]),
                     x2 = ifelse(i_l>=0,pointX$V1[l_p$V2],pointX$V1[l_p$V1]),
                     y2 = ifelse(i_l>=0,pointX$V2[l_p$V2], pointX$V2[l_p$V1]))%>%
      mutate(id_cell = i)
    ggplot(x_point)+geom_polygon(aes(x1, y1))+coord_fixed()
    data = rbind(data, x_point)
  }
  Phys = paste0("Physical Surface(",1:length(Surf),") = {",1:length(Surf),"};")
  geo = c(geo, Phys)
  if(!is.null(out_geo)){
    write(geo, out_geo)
  }

  
  return(data)
}


read_skeleton <- function(path= "./data/Branch_P1_col0_S1_T0.csv"){
  data = read.csv(path)
  data%>%
    ggplot()+
    geom_segment(aes(x = V1.x, y = -V1.y, xend = V2.x, yend = -V2.y))+coord_fixed()
  points = paste0("Point(",1:nrow(data),") = {",data$V1.x/1.2,", ",-data$V1.y/1.2,", ",data$V1.z,', 1};')
  init = c('// Gmsh project created on Thu Dec 12 10:58:32 2024', 'SetFactory("OpenCASCADE");')
  geo = c(init, points)
  write(geo, "./data/Branch_P1_col0_S1_T0.geo")
  
}

calculate_area <- function(df) {
  df_sf <- df %>%
    group_by(id_cell) %>%
    summarise(area = as.numeric(st_area(st_polygon(list(
      rbind(cbind(x, y), cbind(x[1], y[1]))  # Ensuring the polygon is closed
    )))), .groups = "drop") 
}

plot_hsd_boxplot <- function(data, x = "X", y = "Y", qqplot = TRUE) {
  
  # Dynamically evaluate x and y
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)
  
  # Summarize: compute max(y) by group x
  data_summarized <- data %>%
    group_by(!!x_sym) %>%
    summarize(Max_y = max(!!y_sym, na.rm = TRUE)) %>%
    ungroup()
  
  # Perform ANOVA and HSD test
  formula_str <- as.formula(paste(y, "~", x))
  hsd <- agricolae::HSD.test(aov(formula_str, data = data), x, group = TRUE)
  
  # Convert HSD result into dataframe
  hsd_df <- data.frame(hsd$groups) %>%
    mutate(!!x_sym := row.names(hsd$groups)) %>%
    rename(groups = groups) %>%
    select(-matches(y))
  
  # Merge summarized data with HSD groups
  data_summarized <- left_join(data_summarized, hsd_df, by = x)
  
  # MAIN BOXPLOT
  p_box <- ggplot(data, aes(x = reorder(!!x_sym, !!y_sym, FUN = median), y = !!y_sym)) +
    geom_boxplot(aes(fill = !!x_sym), alpha = 0.6) +
    geom_text(data = data_summarized,
              aes(x = !!x_sym, y = Max_y * 1.05, label = groups),
              vjust = 0) +
    theme_classic() +
    viridis::scale_fill_viridis(discrete = TRUE) +
    ylab(y) +
    xlab(x) +
    ggtitle(paste("HSD group comparison for", y, "by", x))+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  # QQPLOT for residuals
  model <- aov(formula_str, data = data)
  residuals_model <- residuals(model)
  p_qq <- ggplot(data.frame(resid = residuals_model), aes(sample = resid)) +
    stat_qq() +
    stat_qq_line() +
    theme_classic() +
    ggtitle("QQ plot of ANOVA residuals")
  
  # Combine plots
  if(qqplot){
    gridExtra::grid.arrange(p_box, p_qq, ncol = 2)
  }
  else{
    p_box
  }
  
}

growth_tensor <- function(data, h=8*60){
  data = data%>%mutate(sig_xx = 0, sig_xy = 0, sig_yy = 0)
  for(i in unique(data$id_cell)){
    coords = data%>%filter(id_cell == i)%>%select(x,y)
    # Compute the covariance matrix
    cov_matrix <- cov(coords)
    eigen_decomp <- eigen(cov_matrix)
    eigenvectors <- eigen_decomp$vectors
    eigenvalues <- c(data$norm_diff_var_major_l[data$id_cell == i][1]/h, data$norm_diff_var_minor_l[data$id_cell == i][1]/h)
    eigenvalues[eigenvalues < 0] = 0
    # Construct the tensor matrix for each cell
    tensor_matrix <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
    data$sig_xx[data$id_cell == i] = tensor_matrix[1,1]
    data$sig_xy[data$id_cell == i] = tensor_matrix[1,2]
    data$sig_yy[data$id_cell == i] = tensor_matrix[2,2]
  }
  return(data)
}

building_strain_cases <- function(data) {
  
  # --- Helper function to apply a given growth rule and compute tensors ---
  make_case <- function(data, name, major_expr, minor_expr) {
    case <- data %>%
      mutate(
        norm_diff_var_major_l = !!rlang::parse_expr(major_expr),
        norm_diff_var_minor_l = !!rlang::parse_expr(minor_expr)
      )
    case <- growth_tensor(case)
    case %>% mutate(growth = name)
  }
  
  # --- Define growth scenarios as a list ---
  scenarios <- list(
    heter_ani  = list("norm_diff_var_major_l", "norm_diff_var_minor_l"),
    slow_ani   = list("norm_diff_var_major_l*0.8", "norm_diff_var_minor_l*0.8"),
    high_ani   = list("norm_diff_var_major_l*1.2", "norm_diff_var_minor_l*1.2"),
    dcr = list("L_dcr", "l_dcr"),
    cyp = list("L_cyp", "l_cyp"),
    col = list("L_col", "l_col")
  )
  
  # --- Build the standard cases automatically ---
  case_list <- purrr::imap(scenarios, ~make_case(data, .y, .x[[1]], .x[[2]]))
  
  # Extract hetergrowth_ani for downstream computations
  case_hetergrowth_ani <- dplyr::filter(case_list$heter_ani, TRUE)
  
  # --- Compute mean strains (anisotropic reference) ---
  mean_strain_bin <- case_hetergrowth_ani %>%
    group_by(id_cell,bin) %>%
    distinct(norm_diff_var_major_l,norm_diff_var_minor_l)
  
  make_bin_growth <- function(data, mean_strain_bin, target_bins, ref_bins, name) {
    data %>%
      mutate(
        norm_diff_var_major_l = ifelse(
          bin %in% target_bins,
          mean(mean_strain_bin$norm_diff_var_major_l[mean_strain_bin$bin %in% ref_bins], na.rm = TRUE),
          norm_diff_var_major_l
        ),
        norm_diff_var_minor_l = ifelse(
          bin %in% target_bins,
          mean(mean_strain_bin$norm_diff_var_minor_l[mean_strain_bin$bin %in% ref_bins], na.rm = TRUE),
          norm_diff_var_minor_l
        )
      ) %>%
      growth_tensor() %>%
      mutate(growth = name)
  }
  
  case_homoggrowth_bin_apiave <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-apical", "outer-apical"),
    ref_bins = c("inner-apical", "outer-apical"),
    name = "homog_apiave"
  )
  
  case_homoggrowth_bin_apiinn <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-apical","outer-apical"),
    ref_bins = c("inner-apical"),
    name = "homog_apiinn"
  )
  
  case_homoggrowth_bin_apiout <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-apical","outer-apical"),
    ref_bins = c("outer-apical"),
    name = "homog_apiout"
  )
  
  # --- and similarly for basal side ---
  case_homoggrowth_bin_basave <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-basal", "outer-basal"),
    ref_bins = c("inner-basal", "outer-basal"),
    name = "homog_basave"
  )
  
  case_homoggrowth_bin_basinn <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-basal", "outer-basal"),
    ref_bins = c("inner-basal"),
    name = "homog_basinn"
  )
  
  case_homoggrowth_bin_basout <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-basal", "outer-basal"),
    ref_bins = c("outer-basal"),
    name = "homog_basout"
  )
  
  # --- Derived homogeneous/isotropic scenarios ---
  case_homoggrowth_ani <- case_hetergrowth_ani %>%
    mutate(
      norm_diff_var_major_l = mean(mean_strain_bin$norm_diff_var_major_l, na.rm = TRUE),
      norm_diff_var_minor_l = mean(mean_strain_bin$norm_diff_var_minor_l, na.rm = TRUE)
    ) %>%
    growth_tensor() %>%
    mutate(growth = "homog_ani")
  
  
  case_homoggrowth_iso <- case_hetergrowth_ani %>%
    mutate(
      diff_var_major_l = mean(
        c(mean_strain_bin$norm_diff_var_major_l ,
          mean_strain_bin$norm_diff_var_minor_l ), na.rm = TRUE),
      diff_var_minor_l = diff_var_major_l
    ) %>%
    growth_tensor() %>%
    mutate(growth = "homog_iso")
  
  case_hetergrowth_iso <- case_hetergrowth_ani %>%
    mutate(
      norm_diff_var_major_l = (norm_diff_var_major_l + norm_diff_var_minor_l) / 2,
      norm_diff_var_minor_l = (norm_diff_var_major_l + norm_diff_var_minor_l) / 2
    ) %>%
    growth_tensor() %>%
    select(1:ncol(case_homoggrowth_ani)) %>%
    mutate(growth = "heter_iso")
  
  # --- Combine everything ---
  all_cases <- bind_rows(
    case_list)
  all_cases = rbind(all_cases,
                    case_homoggrowth_iso,
                    case_homoggrowth_ani,
                    case_hetergrowth_iso,
                    case_homoggrowth_bin_apiave,
                    case_homoggrowth_bin_apiinn,
                    case_homoggrowth_bin_apiout,
                    case_homoggrowth_bin_basave,
                    case_homoggrowth_bin_basinn,
                    case_homoggrowth_bin_basout
  )
  
  return(all_cases)
}

building_growth_cases <- function(data) {
  
  # --- Helper function to apply a given growth rule and compute tensors ---
  make_case <- function(data, name, major_expr, minor_expr) {
    case <- data %>%
      mutate(
        diff_var_major_l = !!rlang::parse_expr(major_expr),
        diff_var_minor_l = !!rlang::parse_expr(minor_expr)
      )
    case <- growth_tensor(case)
    case %>% mutate(growth = name)
  }
  
  # --- Define growth scenarios as a list ---
  scenarios <- list(
    heter_ani  = list("diff_var_major_l/1000", "diff_var_minor_l/1000"),
    slow_ani   = list("diff_var_major_l/2000", "diff_var_minor_l/2000"),
    high_ani   = list("diff_var_major_l/500", "diff_var_minor_l/500"),
    dcr = list("(L_dcr*major_axis_length)/1000", "(l_dcr*minor_axis_length)/1000"),
    cyp = list("(L_cyp*major_axis_length)/1000", "(l_cyp*minor_axis_length)/1000"),
    col = list("(L_col*major_axis_length)/1000", "(l_col*minor_axis_length)/1000")
  )
  
  # --- Build the standard cases automatically ---
  case_list <- purrr::imap(scenarios, ~make_case(data, .y, .x[[1]], .x[[2]]))
  
  # Extract hetergrowth_ani for downstream computations
  case_hetergrowth_ani <- dplyr::filter(case_list$heter_ani, TRUE)
  
  # --- Compute mean strains (anisotropic reference) ---
  mean_strain <- case_hetergrowth_ani %>%
    group_by(id_cell) %>%
    distinct(diff_var_major_l,diff_var_minor_l,base_major_l, base_minor_l) %>%
    mutate(maxda = diff_var_major_l / base_major_l, minda = diff_var_minor_l / base_minor_l)
  
  mean_strain_bin <- case_hetergrowth_ani %>%
    group_by(id_cell, bin) %>%
    distinct(diff_var_major_l,diff_var_minor_l,base_major_l, base_minor_l) %>%
    mutate(maxda = diff_var_major_l / base_major_l, minda = diff_var_minor_l / base_minor_l)
  
  make_bin_growth <- function(data, mean_strain_bin, target_bins, ref_bins, name) {
    data %>%
      mutate(
        diff_var_major_l = ifelse(
          bin %in% target_bins,
          mean(mean_strain_bin$maxda[mean_strain_bin$bin %in% ref_bins], na.rm = TRUE) * base_major_l,
          diff_var_major_l
        ),
        diff_var_minor_l = ifelse(
          bin %in% target_bins,
          mean(mean_strain_bin$minda[mean_strain_bin$bin %in% ref_bins], na.rm = TRUE) * base_minor_l,
          diff_var_minor_l
        )
      ) %>%
      growth_tensor() %>%
      mutate(growth = name)
  }
  
  case_homoggrowth_bin_apiave <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-apical", "outer-apical"),
    ref_bins = c("inner-apical", "outer-apical"),
    name = "homog_apiave"
  )
  
  case_homoggrowth_bin_apiinn <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-apical","outer-apical"),
    ref_bins = c("inner-apical"),
    name = "homog_apiinn"
  )
  
  case_homoggrowth_bin_apiout <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-apical","outer-apical"),
    ref_bins = c("outer-apical"),
    name = "homog_apiout"
  )
  
  # --- and similarly for basal side ---
  case_homoggrowth_bin_basave <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-basal", "outer-basal"),
    ref_bins = c("inner-basal", "outer-basal"),
    name = "homog_basave"
  )
  
  case_homoggrowth_bin_basinn <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-basal", "outer-basal"),
    ref_bins = c("inner-basal"),
    name = "homog_basinn"
  )
  
  case_homoggrowth_bin_basout <- make_bin_growth(
    data = case_hetergrowth_ani,
    mean_strain_bin = mean_strain_bin,
    target_bins = c("inner-basal", "outer-basal"),
    ref_bins = c("outer-basal"),
    name = "homog_basout"
  )
  
  # --- Derived homogeneous/isotropic scenarios ---
  case_homoggrowth_ani <- case_hetergrowth_ani %>%
    mutate(
      diff_var_major_l = mean(mean_strain$maxda, na.rm = TRUE)*base_major_l,
      diff_var_minor_l = mean(mean_strain$minda, na.rm = TRUE)*base_minor_l
    ) %>%
    growth_tensor() %>%
    mutate(growth = "homog_ani")
  

  case_homoggrowth_iso <- case_hetergrowth_ani %>%
    mutate(
      diff_var_major_l = mean(
        c(mean_strain$maxda*mean_strain$base_major_l ,
          mean_strain$minda*mean_strain$base_minor_l ), na.rm = TRUE),
      diff_var_minor_l = diff_var_major_l
    ) %>%
    growth_tensor() %>%
    mutate(growth = "homog_iso")
  
  case_hetergrowth_iso <- case_hetergrowth_ani %>%
    mutate(
      diff_var_major_l = (diff_var_major_l + diff_var_minor_l) / 2,
      diff_var_minor_l = (diff_var_major_l + diff_var_minor_l) / 2
    ) %>%
    growth_tensor() %>%
    select(1:ncol(case_homoggrowth_ani)) %>%
    mutate(growth = "heter_iso")
  
  # --- Combine everything ---
  all_cases <- bind_rows(
    case_list)
  all_cases = rbind(all_cases,
    case_homoggrowth_iso,
    case_homoggrowth_ani,
    case_hetergrowth_iso,
    case_homoggrowth_bin_apiave,
    case_homoggrowth_bin_apiinn,
    case_homoggrowth_bin_apiout,
    case_homoggrowth_bin_basave,
    case_homoggrowth_bin_basinn,
    case_homoggrowth_bin_basout
  )
  
  return(all_cases)
}

# Define the standardization function
standardize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Function to apply standardization to all numeric columns in a dataframe
standardize_dataframe <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      standardize(col)
    } else {
      col  # Keep non-numeric columns unchanged
    }
  })
  return(df)
}

# Function to add rank within each group
add_group_rank <- function(df, group_col, value_col, desc = T) {
  if(desc){
    df <- df %>%
      group_by(.data[[group_col]]) %>%
      mutate(rank = rank(desc(.data[[value_col]]), ties.method = "first")) %>%
      ungroup()
  }else{
    df <- df %>%
      group_by(.data[[group_col]]) %>%
      mutate(rank = rank(.data[[value_col]], ties.method = "first")) %>%
      ungroup()
  }

  
  return(df)
}

ordering_cells_archive <- function(data){
  data_T0 = data %>%
    filter(time == 0)%>%
    mutate(id_plant = paste0(Gen,pool,sampl))
  
  out_order = NULL
  for(i in unique(data_T0$id_plant)){
    tmp_data = data_T0%>%filter(id_plant == i)
    tmp = tibble(id_cell = unique(tmp_data$id_cell))
    
    dist_mat= dist(tmp)
    hclust_avg <- hclust(dist_mat, method = 'average')
    k = parse_number(str_sub(paste0(max(tmp_data$id_cell)),1,2))
    if (k > 20){k =11}
    cut_avg <- cutree(hclust_avg, k)
    tmp$cut = cut_avg
    tmp = left_join(tmp, tmp_data%>%group_by(id_cell)%>%distinct(atan0), by = "id_cell")
    tmp_i = add_group_rank(tmp, "cut", "atan0")%>%
      select(-atan0)
    
    out_order= rbind(out_order, tmp_i %>% mutate(id_plant = i))
  }
  
  data_T0 = left_join(data_T0,
                   out_order, by = c("id_plant", "id_cell"))
  data_T8 = left_join(data %>%
                        filter(time == 8)%>%
                        mutate(id_plant = paste0(Gen,pool,sampl)),
                      out_order, by = c("id_plant", "id_cell"))
  return(rbind(data_T0, data_T8))
}


write_case <- function(data, prefix = "P2",dir_path = "~/GitHub/SRobertGroup/2D-HypoHookFEM/data/in/growth_data/"){
  idlabel_young = tibble(ID_label = seq(0,length(unique(data$id_cell))), Young_value = 10)
  # Overwrite the Young modulus for a first case
  idlabel_young$Young_value[1] = 20
  
  for(i_case in unique(data$growth)){
    tmp = data %>%filter(growth == i_case)
    idlabel_growth = left_join(idlabel_young,
                               tmp%>%
                                 dplyr::distinct(id_cell, sig_xx, sig_xy, sig_yy)%>%
                                 arrange(id_cell)%>%
                                 mutate(ID_label = 1:length(unique(data$id_cell))), by = 'ID_label')%>%
      mutate(across(everything(), ~ tidyr::replace_na(.x, 0)))
    
    write.csv(idlabel_growth,paste0(dir_path,prefix,"_",i_case,".csv"))
    
  }
}


ordering_cells <- function(data) {
  
  # helper for scaling
  scale_to_10 <- function(x) {
    if (length(unique(x)) == 1) return(rep(10, length(x)))
    (x - min(x)) / (max(x) - min(x)) * 9
  }
  
  order_one_plant <- function(df) {
    tmp_data <- df %>% filter(time == 0)
    
    # clustering on cells
    tmp <- tibble(id_cell = unique(tmp_data$id_cell))%>%
      arrange(id_cell)
    dist_mat <- dist(tmp)
    hclust_avg <- hclust(dist_mat, method = "average")
    k <- readr::parse_number(str_sub(paste0(max(tmp_data$id_cell)), 1, 2))
    if (k > 20) k <- 11
    cut_avg <- cutree(hclust_avg, k)
    tmp$cut <- cut_avg
    
    # merge back atan
    tmp <- left_join(tmp,
                     tmp_data %>% group_by(id_cell) %>% distinct(m_atan, x, y, diff_var_area, mx, my),
                     by = "id_cell")

    
    # compute whether each cell is above/below the line
    tmp_i <- tmp %>%
      mutate(inner_outer = scale_to_10(cut))
    
    # compute whether each cell is above/below the line
    tmp_i$s_atan = scale_to_10(tmp_i$m_atan)
    tmp_i <- tmp_i %>%
      mutate(basal_apical2 = s_atan
      ) %>%
      select(id_cell, cut, inner_outer, basal_apical2)%>%
      distinct(id_cell, .keep_all = TRUE)   # ensure one row per cell
    
    # join back to both T0 and T8
    df <- left_join(df, tmp_i, by = "id_cell")
    return(df)
  }
  
  data %>%
    mutate(id_plant = paste0(Gen, pool, sampl)) %>%
    group_by(id_plant) %>%
    group_modify(~ order_one_plant(.x)) %>%
    ungroup()
}



find_bottom_bc <- function(data = tibble(x = numeric(), y = numeric()),
         r = 1) {
  
  set_of_points = data%>%
    filter(id_cell != -1) %>% 
    group_by(id_cell, cut, fbasal_apical2)%>%
    summarise(y = min(y), .groups = "drop")%>%
    filter(fbasal_apical2 == min(fbasal_apical2))
  set_of_points = set_of_points %>%
    filter( cut %in% c(7,3,12,9))%>%
    left_join(data, by = c("id_cell","cut","fbasal_apical2" ,"y"))
  

  
  # Convert to sf objects — no CRS assumption (planar space)
  points_sf <- st_as_sf(set_of_points, coords = c("x", "y"), crs = NA_crs_)
  data_sf   <- st_as_sf(data%>%filter(id_cell == -1), coords = c("x", "y"), crs = NA_crs_)
  
  # Helper: make a closed circle polygon
  make_circle <- function(x, y, r, n = 100) {
    angles <- seq(0, 2 * pi, length.out = n)
    coords <- cbind(x + r * cos(angles), y + r * sin(angles))
    coords <- rbind(coords, coords[1, ])  # ensure closure
    st_polygon(list(coords))
  }
  
  # Build circle polygons around each point
  coords <- st_coordinates(points_sf)
  circle_list <- map(seq_len(nrow(coords)),
                     ~ make_circle(coords[.x, 1], coords[.x, 2], r))
  
  circles <- st_sfc(circle_list, crs = NA_crs_)
  buffers <- st_sf(id = seq_along(circles), geometry = circles)
  
  # Spatial join: keep only points that fall within any circle
  result <- st_coordinates(st_join(data_sf, buffers, join = st_intersects, left = FALSE))
  result = as_tibble(result)
  colnames(result) = c("x","y")
  
  print(data%>%
    ggplot(aes(x,y))+
    geom_polygon(aes(group = id_cell), fill= "white", colour = "black")+coord_fixed()+
    geom_point(aes(x,y), data = result, colour = "red"))
  
  return(set_of_points %>%select(x,y))
}

project_point_to_line <- function(x, y, m, b) {
  # point on line
  x0 <- 0
  y0 <- b
  
  # direction vector
  dx <- 1
  dy <- m
  
  # vector from P0 to P
  vx <- x - x0
  vy <- y - y0
  
  # dot products
  t <- (vx * dx + vy * dy) / (dx * dx + dy * dy)
  
  # projected point
  x_proj <- x0 + t * dx
  y_proj <- y0 + t * dy
  
  list(x_proj = x_proj, y_proj = y_proj)
}


prep_base <- function(data, dir_path = NULL){
  

  tmp = data%>%filter(fbasal_apical2 == 1)
  
  center_estimates <- tmp %>%
    distinct(cx,cy,m_atan, m_euc)%>%
    mutate(
      x0 = cx - m_euc * sin(m_atan),
      y0 = cy - m_euc * cos(m_atan)
    )%>%
    summarise(x0 = mean(x0), y0 = mean(y0))
  
  tmp = tmp %>%
    mutate(atan = atan2(y-center_estimates$y0, x-center_estimates$x0))
  
  fit = aov(y~x, tmp%>%group_by(inner_outer)%>%summarise(y = min(y), x = mean(x)))
  m = fit$coefficients[2]
  
  k=10
  while(tmp$atan[tmp$y == min(tmp$y)]<0){
    center_estimates <- tmp %>%
      distinct(cx,cy,m_atan, m_euc)%>%
      mutate(
        x0 = cx - m_euc * sin(m_atan),
        y0 = cy - m_euc * cos(m_atan)
      )%>%
      summarise(x0 = mean(x0), y0 = mean(y0)-k)
    
    tmp = tmp %>%
      mutate(atan = atan2(y-center_estimates$y0, x-center_estimates$x0))
    k = k +1
  }
  
  tmp = tmp%>%mutate(b = y-m*x)
  b = min(tmp$b)-5
  
  angle_limits = tmp %>%dplyr::group_by(id_cell)%>%
    summarise(
      max_a = max(atan, na.rm = TRUE),
      min_a = min(atan, na.rm = TRUE),
      .groups = "drop"
    )
  selected_pts <- tmp %>%
    inner_join(angle_limits, by = "id_cell") %>%
    filter(atan == max_a | atan == min_a)
  
  hull_indices <- chull(selected_pts$x, selected_pts$y)
  convex_hull <- selected_pts[hull_indices, ]
  
  x11()
  plot(convex_hull$x, convex_hull$y, asp = 1)
  abline(coef = c(b,m))
  pts <- locator(n = 2)
  dev.off()
  
  
  # click exactly two points on the polygon
  selected_points <- sapply(1:2, function(i) {
    which.min(sqrt((convex_hull$x - pts$x[i])^2 + (convex_hull$y - pts$y[i])^2))
  })
  
  projections <- lapply(1:2, function(i) {
    project_point_to_line(
      x = convex_hull$x[selected_points[i]],
      y = convex_hull$y[selected_points[i]],
      m = m,
      b = b-3
    )
  })
  # convert to a dataframe
  pol_proj <- data.frame(
    x  = convex_hull$x[selected_points],
    y  = convex_hull$y[selected_points],
    x_proj = sapply(projections, `[[`, "x_proj"),
    y_proj = sapply(projections, `[[`, "y_proj")
  ) %>%
    arrange(x_proj)
  
  x_end = c(pol_proj$x[1], pol_proj$x_proj[1], pol_proj$x_proj[2], pol_proj$x[2], pol_proj$x[1])
  y_end = c(pol_proj$y[1], pol_proj$y_proj[1], pol_proj$y_proj[2], pol_proj$y[2], pol_proj$y[1])
  
  end = st_polygon(list(as.matrix(cbind(x = x_end, y = y_end))))
  
  # sf obj of hull
  hull_poly <- st_polygon(list(as.matrix(rbind(convex_hull[, c("x", "y")], convex_hull[1, c("x", "y")]) )))
  
  ggplot()+geom_sf(data = end)+
    geom_abline(slope = m, intercept = b-3)
  
  hull_union <- st_union(st_sfc(hull_poly), st_sfc(end))%>% st_sf()
  print(hull_union)
  hull_union_df <- st_coordinates(hull_union) %>% as.data.frame()
  names(hull_union_df) <- c("x", "y", "id")
  
  
  outside = data %>%
    filter(id_cell == -1)
  wall = st_polygon(list(as.matrix(rbind(outside[, c("x", "y")], outside[1, c("x", "y")]) )))
  
  wall = st_union(wall,st_union(st_sfc(hull_poly), st_sfc(end)))
  wall = nngeo::st_remove_holes(wall)
  
  ggplot()+geom_sf(data = wall)
  
  simpl <- sf::st_simplify(wall, dTolerance = 0.1, preserveTopology = TRUE)
  hull_union_df <- st_coordinates(simpl) %>% as.data.frame()
  names(hull_union_df) <- c("x", "y", "id")
  
  if(nrow(hull_union_df) < nrow(outside)){
    data[data$id_cell == -1, c("x", "y")][1:nrow(hull_union_df),] = hull_union_df[, c("x", "y")]
    data[data$id_cell == -1, c("x", "y")][(nrow(hull_union_df)+1):nrow(outside),] = NA
    data = data %>%filter(!is.na(x))%>%
      group_by(id_plant,time, id_cell)%>%
      mutate(
        x1 = x,
        y1 = y,
        x2 = lead(x, default = first(x)),  # connect last to first
        y2 = lead(y, default = first(y))
      ) %>%
      ungroup()
  }
  bc = tibble(m = m, b = b-3)
  prefix = data$id_plant[1]
  write.csv(bc, paste0(dir_path,prefix,"_bc.csv"))
  return(data)
}




