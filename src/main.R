

########################################
# Hook analysis
# generation of mesh
# Author of the file: Adrien Heymans
# Date: oct 2025
#
# Function in io_function
########################################

setwd("~/GitHub/SRobertGroup/2D-HypoHookFEM/")
source("./src/R/io_function.R")

# Load all cell data
cell_data = load_data_hook()

# Example
cell_data %>%
  filter(time == 0, id_cell != -1)%>%
  ggplot(aes(x,y))+
  geom_polygon(aes(fill = diff_var_major_l/base_major_l, group = factor(id_cell)))+
  coord_fixed()+
  facet_grid(Gen*time~sampl*pool)+
  theme_dark()+
  viridis::scale_fill_viridis(option = "H")+
  labs(fill = "Cell length\ngrowth [µm]")


# Generate corresponding matrix
# transferring growth prop of corresponding matrix 
corresponding_matrix = make_corresp_matrix(cell_data)

plot_data = corresponding_matrix %>%
  group_by(finner_outer, fbasal_apical2, Gen)%>%
  summarise(strain_rate = mean(norm_diff_var_area, na.rm = T),
            elonga_rate = mean(norm_diff_var_major_l, na.rm = T), 
            L_col = mean(L_col, na.rm = T),
            L_cyp = mean(L_cyp, na.rm = T),
            L_dcr = mean(L_dcr, na.rm = T),.groups = "drop")%>%
  mutate(L = ifelse(Gen == "Col0", L_col, ifelse(Gen == "cyp", L_cyp, L_dcr)))


corresponding_matrix %>%
  filter(Gen == "Col0", bin != "NA-NA")%>%
  distinct(inner_bin, basal_bin, norm_diff_var_major_l)%>%
  ggplot()+
  geom_boxplot(aes(inner_bin, norm_diff_var_major_l))+
  facet_wrap(~basal_bin)+theme_classic()



# Plot on pseudo hook the growth gradient
ggplot(plot_data%>%filter(Gen== "Col0"), aes(x = fbasal_apical2, y = as.numeric(finner_outer), fill = elonga_rate)) +
  geom_tile() +
  facet_wrap(~Gen) +
  coord_radial(start = -0.532 * pi, end = 0.2 * pi, inner.radius = 0.4)+
  viridis::scale_fill_viridis(option = "turbo", limits = c(0.18,1.7)) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

  ggplot(plot_data%>%filter(Gen== "Col0"), aes(x = fbasal_apical2, y = as.numeric(finner_outer), fill = elonga_rate)) +
  geom_tile() +
  facet_wrap(~Gen) +
  coord_radial(start = -0.532 * pi, end = 0.2 * pi, inner.radius = 0.4)+
  viridis::scale_fill_viridis(option = "turbo", limits = c(0.18,1.7)) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )





ggsave("./data/out/fig/fig01_elong_corrmat.svg")

plot_data = corresponding_matrix %>%
  group_by(basal_bin, inner_bin, Gen)%>%
  summarise(strain_rate = mean(norm_diff_var_area, na.rm = T),
            elonga_rate = mean(norm_diff_var_major_l, na.rm = T), .groups = "drop")%>%
  mutate(inner_bin = ifelse(inner_bin == "inner", 1, 2),
         basal_bin = ifelse(basal_bin == "basal", 1, 2))

plot_hsd_boxplot(data = corresponding_matrix%>%filter(Gen == "Col0")%>%
                   filter(!is.na(basal_bin)), x = "bin", 
                 y = "norm_diff_var_major_l",qqplot = F)
plot_hsd_boxplot(data = corresponding_matrix%>%filter(Gen == "dcr")%>%
                   filter(!is.na(basal_bin)), x = "bin", 
                 y = "norm_diff_var_major_l",qqplot = F)
plot_hsd_boxplot(data = corresponding_matrix%>%filter(Gen == "cyp")%>%
                   filter(!is.na(basal_bin)), x = "bin", 
                 y = "norm_diff_var_major_l",qqplot = F)


corresponding_matrix %>% 
  filter(!is.na(inner_bin))%>%
        ggplot() +
        geom_boxplot(aes(x = basal_bin, norm_diff_var_major_l, fill = factor(inner_bin)))+
        theme_bw()+
        facet_wrap(~Gen)+
        labs(y = 'strain [-]', fill = "mean strain [-]")

ggsave("~/GitHub/SRobertGroup/2D-HypoHookFEM/box_plot_normelon_col.svg")

# Plot on pseudo hook the growth gradient
ggplot(plot_data, aes(x = as.numeric(basal_bin), y = as.numeric(inner_bin), fill = elonga_rate)) +
  geom_tile() +
  facet_wrap(~Gen) +
  coord_radial(start = -0.532 * pi, end = 0.2 * pi, inner.radius = 0.4)+
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

plot_data = corresponding_matrix %>%
  group_by(finner_outer, fbasal_apical2, Gen)%>%
  summarise(strain_rate = mean(norm_diff_var_area, na.rm = T),
            elonga_rate = mean(norm_diff_var_major_l, na.rm = T), .groups = "drop")


# Plot on pseudo hook the growth gradient
ggplot(plot_data, aes(x = as.numeric(fbasal_apical2), y = as.numeric(finner_outer), fill = elonga_rate)) +
  geom_tile() +
  facet_wrap(~Gen) +
  coord_radial(start = -0.532 * pi, end = 0.2 * pi, inner.radius = 0.4)+
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

ggsave("./elong_corr_mat.png")

# Stat test on quadrant
# Strain rate, normalized elongation rate, 
quadrant_stat_analysis(corresponding_matrix) # t-test + Levene test


Col0_P3S1 = corresponding_matrix%>%
  filter(id_plant == "Col0P3S1", id_cell != -1)

Col0_P3S2 = corresponding_matrix%>%
  filter(id_plant == "Col0P3S2", id_cell != -1)

# Construct the growth scenarios
WT_scenario = building_strain_cases(data = corresponding_matrix%>%
                                      filter(Gen == "Col0", id_cell != -1))


# Example
WT_scenario %>%  filter(id_cell != max(WT_scenario$id_cell))%>%
  ggplot(aes(x,y))+
  geom_polygon(aes(fill = norm_diff_var_major_l, group = factor(id_cell)))+
  coord_fixed()+
  facet_wrap(~growth)+  scale_fill_viridis_c(option = "turbo")

plot_data = WT_scenario %>%
  group_by(finner_outer, fbasal_apical2, growth)%>%
  summarise(strain_rate = mean(norm_diff_var_area, na.rm = T),
            elonga_rate = mean(norm_diff_var_major_l, na.rm = T), .groups = "drop")

ggplot(plot_data%>%
         filter(growth %in% c("homog_apiinn", "homog_apiout", "homog_basinn", "homog_basout")),
       aes(x = as.numeric(fbasal_apical2), y = as.numeric(finner_outer), fill = elonga_rate)) +
  geom_tile() +
  facet_wrap(~growth) +
  coord_radial(start = -0.532 * pi, end = 0.2 * pi, inner.radius = 0.4)+
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

ggsave("~/GitHub/SRobertGroup/2D-HypoHookFEM/heatmap_elong_api_base.svg")


for(i in unique(corresponding_matrix$id_plant)){
  message(i)
  print(corresponding_matrix %>%
    filter(time == 0, id_cell != -1, id_plant == i)%>%
    ggplot(aes(x,y))+
    geom_polygon(aes(fill = norm_diff_var_area, group = factor(id_cell)))+
    coord_fixed()+
    facet_grid(Gen*time~sampl*pool)+
    theme_dark()+
    viridis::scale_fill_viridis(option = "H")+
    labs(fill = "Cell length\ngrowth [µm]"))
  
  
  data_i = corresponding_matrix%>%
    filter(id_plant == i)
  prefix = data_i$id_plant[1]
  
  scenario = building_strain_cases(data_i%>%filter(id_cell != -1))
  
  write_case(scenario, prefix)
  
  # write_geo(data_i, 
  #           path_geo = paste0("./data/in/mesh/",i,".geo"),
  #           Celldomain = T,
  #           make_base = T)
}



csv_path = "./data/out/raw_csv/"

files <- list.files(csv_path, recursive = TRUE, full.names = TRUE)

data <- files |>
  map_dfr(function(f) {
    
    parts <- str_split(f, .Platform$file.sep)[[1]]
    fo <- parts[length(parts) - 2]
    gr_prop <- parts[length(parts) - 1]
    fl <- basename(f)
    
    read.csv(f) |>
      transmute(
        x = Points.0,
        y = Points.1,
        dx = Displacement.Vector.0,
        dy = Displacement.Vector.1,
        id = vtkOriginalPointIds,
        base = row_number(),
        time = parse_number(fl),
        growth = gr_prop,
        name = fo
      )
  })

data = data %>% 
  mutate(xx = x +dx, yy = y +dy)

data %>% filter(name=="Col0P3S1", growth=="col") %>% 
  ggplot(aes(xx,yy, colour = time))+geom_point()

get_lm_segment <- function(df) {
  fit <- lm(y ~ x, data = df)
  
  xr <- range(df$x, na.rm = TRUE)
  tibble(x1= xr[1],x2=xr[2],
           y1=predict(fit, newdata = data.frame(x = xr[1])),
           y2=predict(fit, newdata = data.frame(x = xr[2])),
           slope = coef(fit)[2],
           intercept = coef(fit)[1])
}


cluster_data = data
cluster_data$vector = NA

for (name_i in unique(cluster_data$name)){
  df = cluster_data%>%filter(name == name_i)
  dist_mat= dist(df[, c("x", "y")])
  hclust_avg <- hclust(dist_mat, method = 'average')
  cut <- cutree(hclust_avg, 2)
  cluster_data$vector[cluster_data$name == name_i] = cut
}

# Which vector is apical and which is basal?
cluster_data = cluster_data %>% 
  group_by(name, growth) %>% 
  mutate(vector_type = ifelse(x < mean(x), "Basal", "Apical")) %>% 
  ungroup()

cluster_data%>%
  ggplot(aes(xx,yy, colour = vector_type))+geom_point()+facet_wrap(growth~name)


vector_data <- cluster_data %>%
  mutate(x = xx, y = yy) %>%
  group_by(time, name, growth, vector_type) %>%
  summarize(
    get_lm_segment(cur_data()),
    .groups = "drop"
  )

data = vector_data%>%transmute(sample_name = name,
                               growth_properties = growth,
                               side = vector_type,
                               time,
                               x1, x2,
                               y1, y2,
                               slope, intercept)


# read csv from folder
path_kinematics = "./data/out/kinematics/"
fls <- list.files(path_kinematics)
fls = fls[grepl(".csv", fls)]

# ex file name Col0P3S1_api_in_angles.csv
# name = "Col0P3S1"
# growth_label = "api_in"
data_kinematics <- map_dfr(
  fls,
  function(f) {
    read.csv(file.path(path_kinematics, f)) %>%
      mutate(
        file_base = sub("_angles\\.csv$", "", f),
        name = sub("_.*$", "", file_base),
        growth_label = sub("^[^_]*_", "", file_base)
      )
  }
)
head(data_kinematics)

# Figure 1 
# name == "Col0P3S2"
# growth_label %in% c("api_in", "api_out", "base_in", "base_out", "Col")
# first data point 0, 180°
ref_angle = data_kinematics %>%
group_by(name, growth_label) %>% filter(time == min(time))%>%
mutate(ref_angle = -final_angle)%>%ungroup()%>%select(name, growth_label, ref_angle)

data_kinematics = data_kinematics%>%left_join(ref_angle, by = c("name", "growth_label"))
data_kinematics%>%filter(name == "Col0P3S2", growth_label %in% c("api_in", "api_out", "base_in", "base_out", "Col"))%>%
mutate(angle = 180 - final_angle - ref_angle)%>%
ggplot(aes(time/60, angle, colour = growth_label)) +
  geom_smooth(size = 1.2, alpha = 0.8, se = F) +
  facet_wrap(~name) +
  theme_classic() +
  labs(
    x = "Time [h]",
    y = "Hook Angle [°]",
    colour = "Growth:",
    title = "Apical Hook Angle Kinetics"
  ) +
  scale_color_viridis_d(option = "turbo")


# Figure 2: Combined Facets
# Define the groups for each panel
panel_a_files <- c("Col0P3S2_Col", "Col0P3S2_dcr", "Col0P3S2_cyp")
panel_b_files <- c("Col0P3S2_Col", "dcrP2S1_dcr", "cypP3S1_cyp", "dcrP2S1_Col", "cypP3S1_Col")


# Create subset for Panel A
df_a <- data_kinematics %>%
  filter(file_base %in% panel_a_files) %>%
  mutate(panel = "Panel A")

# Create subset for Panel B
df_b <- data_kinematics %>%
  filter(file_base %in% panel_b_files) %>%
  mutate(panel = "Panel B")

# Combine and Plot
bind_rows(df_a, df_b) %>%
  mutate(angle = 180 - final_angle - ref_angle) %>%
  ggplot(aes(time/60, angle, colour = file_base)) +
  geom_smooth(size = 1.2, alpha = 0.8, se = FALSE, method = "loess", span = 2) +
  facet_wrap(~panel, scales = "free_x") + 
  theme_classic() +
  labs(
    x = "Time [h]",
    y = "Hook Angle [°]",
    colour = "Growth:",
    title = "Apical Hook Angle Kinetics"
  ) +
  scale_color_viridis_d(option = "turbo") +
  ylim(150, 260)

ggsave("./kinematics_figure2_combined.svg", height = 6, width = 12)

