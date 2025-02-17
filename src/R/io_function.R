

library(tidyverse)
library(sp)
library(sf)
library(data.table)

`%!in%` <- compose(`!`, `%in%`)

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
  
  

get_cell_polygons <- function(df_hook, verbatim = F){
  
  id_cell_vector = unique(df_hook$id_cell)
  id_cell_vector = sort(id_cell_vector[id_cell_vector > 0])
  
  Hook_cell = Wall_agregate = NULL
  for(i in id_cell_vector){
    
    tmp = df_hook %>% filter(id_cell == i) %>% 
      mutate(id_point = paste0(round(x, 3), ";", round(y,3))) %>% 
      filter(!duplicated(id_point))

    if (nrow(tmp)>3){
      aa =1.5
      # Get the alphahull shape and convert it to dataframe
      my.ashape = try(alphahull::ashape(x= tmp$x, y = tmp$y, alpha = aa), silent = T)
      while(!class(my.ashape) == "ashape"){
        print("err in ashape function")
        pl = ggplot(data = tmp)+geom_point(aes(x,y))+coord_fixed()
        print(i)
        print(pl)
        tmp = tmp %>% mutate(x = jitter(x, 0.004), y = jitter(y, 0.004))
        my.ashape = try(alphahull::ashape(x= tmp$x, y = tmp$y, alpha = aa), silent = T)
      }
      
      # converted to sf polygons
      a <- data.frame(my.ashape$edges)[,c( 'x1', 'y1', 'x2', 'y2')]
      l <- sf::st_linestring(matrix(as.numeric(a[1,]), ncol=2, byrow = T))
      for(j in 2:nrow(a)){
        l <- c(l, sf::st_linestring(matrix(as.numeric(a[j,]), ncol=2, byrow = T)))
      }
      my_multilinestring = sf::st_sf(geom = sf::st_sfc(l), crs = 2056)
      
      if(sf::st_is_empty(sf::st_union(my_multilinestring)%>% sf::st_polygonize())){
        sf_points <- sf::st_as_sf(tmp, coords = c("x", "y"), crs =  2056)
        concave_hull <- concaveman::concaveman(sf_points, concavity = 2.5, length_threshold = 0)
        r_poly_smooth <- smoothr::smooth(concave_hull, method = "ksmooth", smoothness = 5)
        plot(r_poly_smooth)
      }else{
        alphapoly <-  sf::st_union(my_multilinestring)%>% sf::st_polygonize() %>% sf::st_collection_extract()
        r_poly_smooth <- smoothr::smooth(alphapoly, method = "ksmooth", smoothness = 10)
      }   
      
      are = as.numeric( sf::st_area(r_poly_smooth))
      if(length(are) != 1){
        plot(r_poly_smooth)
        r_poly_smooth = r_poly_smooth[[which(are == max(are))]]
      }
      # Get the coordinates
      coords <- sf::st_coordinates(r_poly_smooth)[, 1:2]

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
      
      perimeter = as.numeric(sf::st_length(sf::st_cast(r_poly_smooth, "MULTILINESTRING")))
      area = as.numeric( sf::st_area(r_poly_smooth))
      circularity = 4*pi*(area/perimeter^2)
      
      simpl = sf::st_simplify(r_poly_smooth, dTolerance = 0.5)
      coords <- sf::st_coordinates(simpl)[, 1:2]
      
      pol = tibble(x = coords[,1], y = coords[,2], 
                   id_cell = i)
      pol = pol %>% 
            mutate(x1 = x,
                   y1 = y,
                   x2 = c(pol$x[-1], pol$x[1]),
                   y2 = c(pol$y[-1], pol$y[1]),
                   major_axis_length,
                   minor_axis_length,
                   perimeter, area, circularity)

      Hook_cell = rbind(Hook_cell, pol)
      
      # Initialize wall
      upWall = scale_polygon(simpl, 1.05)
      wallcoords <- sf::st_coordinates(upWall)[, 1:2]
      wall =  tibble(x = wallcoords[,1], y = wallcoords[,2])
      Wall_agregate = rbind(Wall_agregate, wall)
    }
  }
  
  return(list(Hook_cell, Wall_agregate))
}


# Function to scale down a polygon towards its center of mass
scale_polygon <- function(polygon, scale_factor = 0.8) {
  # Calculate the centroid (center of mass) of the polygon
  centroid <- sf::st_centroid(polygon)
  
  # Extract the coordinates of the polygon and the centroid
  polygon_coords <- sf::st_coordinates(polygon)[, 1:2]
  centroid_coords <- sf::st_coordinates(centroid)[1, 1:2]
  
  x_c = (polygon_coords[,1]-centroid_coords[1])*scale_factor +centroid_coords[1]
  y_c = (polygon_coords[,2]-centroid_coords[2])*scale_factor +centroid_coords[2]
  
  # Create a new scaled polygon
  scaled_polygon <- sf::st_polygon(list(matrix( c(x_c, y = y_c), ncol = 2)))
  
  # Return the scaled polygon as an sf object with the same CRS as the original
  return(sf::st_sfc(scaled_polygon, crs = st_crs(polygon)))
}



write_geo <- function(data, dim = 2, path_geo, averaging = FALSE, extrude=0){
  
  date = Sys.time()
  x1 = paste0('// Gmsh project created on ', date,'\nSetFactory("OpenCASCADE");\n//+\n')
  
  if(dim == 2){
    data$z = 0
    data$z1 = 0
    data$z2 = 0
  }
  
  k = h = j = 1
  txt = x1
  for(i in sort(unique(data$id_cell))){
    tmp = data%>%filter(id_cell == i)
    i_point = unique_point_id(tmp)
    
    i_point$idx = seq(k,-1+k+nrow(i_point),1)
    i_point$idx2 = c(i_point$idx[-1],i_point$idx[1])
    
    dots = paste0("Point(",i_point$idx,") = {",i_point$x, ' ,',i_point$y,' ,',i_point$z, ',1.0};\n//+\n', collapse = "")
    Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
    Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nSurface(',
                  h,") = {", j,'};\n//+\n', collapse = "")
    Physical = paste0("Physical Surface(",h,") = {",h,"};")
    
    txt = paste0(txt, dots, Lines, Surf, Physical)
    
    h = h+1
    j = j+2
    k = max(i_point$idx)+1
  }
  write(txt, path_geo)
}

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


growth_tensor <- function(data){
  data = data%>%mutate(sig_xx = 0, sig_xy = 0, sig_yy = 0)
  for(i in unique(data$id_cell)){
    coords = data%>%filter(id_cell == i)%>%select(x,y)
    # Compute the covariance matrix
    cov_matrix <- cov(coords)
    eigen_decomp <- eigen(cov_matrix)
    eigenvectors <- eigen_decomp$vectors
    eigenvalues <- c(data$diff_var_major_l[data$id_cell == i][1], data$diff_var_minor_l[data$id_cell == i][1])
    eigenvalues[eigenvalues < 0] = 0
    # Construct the tensor matrix for each cell
    tensor_matrix <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
    data$sig_xx[data$id_cell == i] = tensor_matrix[1,1]
    data$sig_xy[data$id_cell == i] = tensor_matrix[1,2]
    data$sig_yy[data$id_cell == i] = tensor_matrix[2,2]
  }
  return(data)
}


building_growth_cases <- function(data){
  # Localized growth with high strain anisotropy
  # Wild type ground truth
  case_hetergrowth_ani = data%>%mutate(diff_var_major_l =  diff_var_major_l/1000, # ifelse(diff_var_major_l<0,0,
                                       diff_var_minor_l =  diff_var_minor_l/1000) # ifelse(diff_var_minor_l<0,0,
  case_hetergrowth_ani = growth_tensor(case_hetergrowth_ani)
  
  case_dcr = data%>%mutate(diff_var_major_l =  (L_dcr*major_axis_length)/1000, # ifelse(diff_var_major_l<0,0,
                           diff_var_minor_l =  (l_dcr*minor_axis_length)/1000) # ifelse(diff_var_minor_l<0,0,
  case_dcr = growth_tensor(case_dcr)
  
  case_cyp = data%>%mutate(diff_var_major_l =  (L_cyp*major_axis_length)/1000, # ifelse(diff_var_major_l<0,0,
                           diff_var_minor_l =  (l_cyp*minor_axis_length)/1000) # ifelse(diff_var_minor_l<0,0,
  case_cyp = growth_tensor(case_cyp)
  
  case_col = data%>%mutate(diff_var_major_l =  (L_col*major_axis_length)/1000, # ifelse(diff_var_major_l<0,0,
                           diff_var_minor_l =  (l_col*minor_axis_length)/1000) # ifelse(diff_var_minor_l<0,0,
  case_col = growth_tensor(case_col)
  
  mean_strain = case_hetergrowth_ani%>%dplyr::group_by(id_cell)%>%
    dplyr::summarise(maxd = mean(diff_var_major_l), 
                     mind = mean(diff_var_minor_l),
                     ma = mean(area), .groups = "drop")%>%
    mutate(maxda = maxd/ma,
           minda = mind/ma)
  
  # distributed growth with anisotropic strain rate
  case_homoggrowth_ani = case_hetergrowth_ani%>%
    mutate(diff_var_major_l = mean(mean_strain$maxda*mean(mean_strain$ma, na.rm = T), na.rm = T),
           diff_var_minor_l = mean(mean_strain$minda*mean(mean_strain$ma, na.rm = T), na.rm = T))
  case_homoggrowth_ani = growth_tensor(case_homoggrowth_ani)
  
  # distributed growth with isotropic strain rate
  case_homoggrowth_iso = case_hetergrowth_ani%>%
    mutate(diff_var_major_l =  mean(c(mean_strain$maxda*mean(mean_strain$ma, na.rm = T), 
                                      mean_strain$minda*mean(mean_strain$ma, na.rm = T)), na.rm = T),
           diff_var_minor_l =  mean(c(mean_strain$maxda*mean(mean_strain$ma, na.rm = T), 
                                      mean_strain$minda*mean(mean_strain$ma, na.rm = T)), na.rm = T))
  case_homoggrowth_iso = growth_tensor(case_homoggrowth_iso)
  
  # Localized growth with isotropic strain rate
  case_hetergrowth_iso = left_join(case_hetergrowth_ani, mean_strain, by = "id_cell")%>%
    mutate(diff_var_major_l = (maxd+ mind)/2,
           diff_var_minor_l = (maxd+mind)/2)
  case_hetergrowth_iso = growth_tensor(case_hetergrowth_iso)[,1:ncol(case_homoggrowth_ani)]
  
  case = rbind(case_homoggrowth_iso%>%mutate(growth = "homog_iso"),
               case_homoggrowth_ani%>%mutate(growth = "homog_ani"),
               case_hetergrowth_iso%>%mutate(growth = "heter_iso"),
               case_hetergrowth_ani%>%mutate(growth = "heter_ani"),
               case_dcr %>%mutate(growth = "dcr"),
               case_cyp %>%mutate(growth = "cyp"),
               case_col %>%mutate(growth = "Col0"))
  return(case)
  
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

ordering_cells <- function(data){
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
  idlabel_young = tibble(ID_label = seq(0,length(unique(data$id_cell))), Young_value = 0.001)
  # Overwrite the Young modulus for a first case
  idlabel_young$Young_value[1] = 0.002
  
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



