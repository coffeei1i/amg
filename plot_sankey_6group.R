

rm(list=ls())
library(networkD3)

read_name="A02782A3"
library(ggplot2)
library(ggalluvial)  
source('sankey_func.R')
sankey_color=c(
  ExN_Avp_Ebf3 = '#ff87ad',
  ExN_Bace2_Met = '#ffb9cd',
  ExN_Car12_Krt20 = '#db6292',
  ExN_Car12_Met = '#ff7091',
  ExN_Ccn3_Dcn = '#ffa389',
  ExN_Ebf2_Aox3 = '#ff8015',
  ExN_Grp_Hgf = '#ff9f28',
  ExN_Lhx5_Reln = '#ffb30a',
  ExN_Lhx9_Abi3bp = '#ec4c50',
  ExN_Otp_Sim1 = '#ff2e7e',
  ExN_Otp_Stc1 = '#ff4345',
  ExN_Rorb_Plch1 = '#e52b5d',
  ExN_Rspo1_Col18a1 = '#d5231f',
  ExN_Rspo2_Npsr1 = '#e77a7b',
  ExN_Sim1_Strip2 = '#ff1e83',
  ExN_Tfap2c_Pou3f1 = '#ff02b5',
  InN_Abcb5_Pax6 = '#d8f075',
  InN_Apob_Prkcq = '#a6e8ab',
  InN_Calcr_Esr2 = '#7AE6AB',
  InN_Drd2_Adora2a = '#B8FFCA',
  InN_Esr2_Cyp19a1 = '#ADE6A6',
  InN_Gal_Dock8 = '#82AD7D',
  InN_Htr3a_Frem1 = '#00979D',
  InN_Htr3a_Glyctk = '#00A809',
  InN_Isl1_Kcnh8 = '#26BF64',
  InN_Lamp5_Id2 = '#73ca96',
  InN_Lamp5_Ndnf = '#00A79D',
  InN_Lhx8_Th = '#2F8C4D',
  InN_Ndnf_Ppp1r1c = '#669D6A',
  InN_Nxph2_Foxp2 = '#588169',
  InN_Prkcd_Chodl = '#5f9365',
  InN_Pvalb_Adamts5 = '#005C07',
  InN_Pvalb_Vipr2 = '#008F1F',
  InN_Sst_Edaradd = '#53879D',
  InN_Sst_Myzap = '#53A39D',
  InN_Sst_Ntn1 = '#174F61',
  InN_Sst_Vipr2 = '#46877b',
  InN_Tshr_Qrfprl = '#7F9922',
  InN_Tshz1_Drd1 = '#C2E32C',
  InN_Tshz1_Lamp5 = '#96E32C',
  InN_Vip_Crh = '#74996b',
  Non_Astro = '#898584',
  Non_COP = '#d6d2d5',
  Non_Endo = '#81705a',
  Non_Fibro = '#a9a9a9',
  Non_Micro = '#7f605c',
  Non_Mural = '#808080',
  Non_OPC = '#7e6b6d',
  Non_Oligo = '#907b77',
  tha = 'blue',
  BA = "#eb6100",
  BMAp = "#e83828",
  COAp = "#fff100",
  PA = "#e60012",
  IAL = "#2ea7e0",
  CEA = "#036eb8",
  BMAa = "#f39800",
  COAa = "#ffe200",
  MEAav = "#ffe200",  
  MEAad = "#00913a",
  MEApv = "#8fc31f",
  MEApd = "#006934",
  LA="#f4b739"
)
setwd('/data')
read_name="A02782A4"

st_obj_final = dior::read_h5(file=paste0(read_name,"_for_python.h5"))
st_obj_final@assays[["counts"]]=NULL
st_obj_final@meta.data[["orig.ident"]]=read_name


read_name="A02782A4"
st_obj_cellbin = dior::read_h5(file=paste0(read_name,"_for_python.h5"))
length(names(table(st_obj_cellbin@meta.data[["predicted_classes4"]])))
unique_cell_types <- names(table(st_obj_cellbin@meta.data[["predicted_classes4"]]))



for (read_name in c("A02782A3","A02781E1")){
  st_obj_cellbin = dior::read_h5(file=paste0(read_name,"_for_python.h5"))
  
  
  
  st_obj_cellbin@assays[["counts"]]=NULL
  st_obj_cellbin@meta.data[["orig.ident"]]=read_name
  st_obj_final=merge(st_obj_final,st_obj_cellbin)
  
}
st_obj_cellbin=st_obj_final
non_indices <- which(!grepl("^Non", st_obj_cellbin@meta.data$predicted_classes4))
st_obj_cellbin=st_obj_cellbin[,non_indices]

table(st_obj_cellbin@meta.data[["orig.ident"]])


st_obj_cellbin=subset(st_obj_cellbin,group2!="MEAPV")


cell_types <- st_obj_cellbin@meta.data[["predicted_classes4"]]
brain_regions <- st_obj_cellbin@meta.data[["group2"]]

data_frame <- data.frame(cell_type = cell_types, brain_region = brain_regions)
library(dplyr)

freq_table <- table(data_frame$cell_type, data_frame$brain_region)

freq_df <- as.data.frame.table(freq_table)
names(freq_df) <- c("cell_type", "brain_region", "freq")


filtered_freq_df <- freq_df %>%
  filter(!grepl("^Group", brain_region),  
         !grepl("Unknown", brain_region))  



filtered_freq_df$brain_region[filtered_freq_df$brain_region == "IAL" | filtered_freq_df$brain_region == "IAM"] <- "IAL"


table(filtered_freq_df$brain_region)


filtered_freq_df <- filtered_freq_df[filtered_freq_df$brain_region != "MEAPV",]


table(filtered_freq_df$brain_region)


complete_gene_order <- c( "ExN_Grp_Hgf","ExN_Rspo2_Npsr1", "ExN_Ccn3_Dcn", "ExN_Bace2_Met", 
                          "ExN_Car12_Krt20","ExN_Car12_Met", "ExN_Sim1_Strip2", "InN_Tshz1_Lamp5", 
                          "InN_Tshz1_Drd1", "InN_Sst_Vipr2", "InN_Drd2_Adora2a", "InN_Prkcd_Chodl", 
                          "InN_Isl1_Kcnh8", "InN_Gal_Dock8", "ExN_Avp_Ebf3", "ExN_Otp_Stc1", 
                          "ExN_Otp_Sim1", "ExN_Rspo1_Col18a1", "ExN_Rorb_Plch1", "ExN_Ebf2_Aox3", 
                          "ExN_Lhx9_Abi3bp", "ExN_Lhx5_Reln", "ExN_Tfap2c_Pou3f1", "InN_Tshr_Qrfprl", 
                          "InN_Apob_Prkcq", "InN_Lhx8_Th", "InN_Nxph2_Foxp2", "InN_Abcb5_Pax6", 
                          "InN_Calcr_Esr2", "InN_Esr2_Cyp19a1", "InN_Sst_Myzap", "InN_Sst_Ntn1", 
                          "InN_Sst_Edaradd", "InN_Pvalb_Adamts5", "InN_Pvalb_Vipr2", "InN_Lamp5_Ndnf", 
                          "InN_Lamp5_Id2", "InN_Vip_Crh", "InN_Htr3a_Glyctk", "InN_Htr3a_Frem1", 
                          "InN_Ndnf_Ppp1r1c"
)
brain_region_order <- c("LA", "BA", "BMAp", "COAp", "PA", "IAL", "CEAL", "CEA","BMAa", "COAa", 
                        "MEAav", "MEAad", "MEApv", "MEApd")


filtered_freq_df$cell_type <- factor(filtered_freq_df$cell_type, levels = complete_gene_order)
filtered_freq_df$brain_region <- factor(filtered_freq_df$brain_region, levels = brain_region_order)



total_freqs <- filtered_freq_df %>%
  group_by(cell_type) %>%
  summarize(total_freq = sum(freq))


filtered_freq_df <- filtered_freq_df %>%
  left_join(total_freqs, by = "cell_type") %>%
  mutate(frequency_percent = freq / total_freq)
i=5

  frequency_threshold <- i * 5
  

  filtered_freq_df2 <- filtered_freq_df %>%
    filter(frequency_percent > 0.05 +0.01*i)
  
  print(length(rownames(filtered_freq_df)))
  print(length(rownames(filtered_freq_df2)))
  
  p <- ggplot(data = filtered_freq_df2,
              aes(axis1 = cell_type, axis2 = brain_region, y = freq, fill = cell_type)) +
    geom_alluvium(aes(width = 0.5)) +  
    geom_stratum(width = 0.1) +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_fill_manual(values = sankey_color) +  
    theme_minimal() +  
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(),  
      legend.position = "none"  
    ) +
    ggtitle(paste("Sankey Diagram of Cell Types and Brain Regions: Threshold", frequency_threshold))  # 添加图表标题
  



filtered_freq_df2 <- filtered_freq_df %>%
  filter(frequency_percent > 0.05 +0.01*0)







p<-ggplot(data = filtered_freq_df2,
          aes(axis1 = cell_type, axis2 = brain_region, y = freq, fill = cell_type)) +
  geom_alluvium(aes(width = 0.5)) +  
  geom_stratum(width = 0.1) +  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values = sankey_color) + 
  theme_minimal() +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    legend.position = "none"  
  ) +
  ggtitle("Sankey Diagram of Cell Types and Brain Regions")  
p


highlight_cell_types <- c("ExN_Grp_Hgf", "ExN_Rspo2_Npsr1", "ExN_Ccn3_Dcn", "ExN_Bace2_Met", "ExN_Car12_Krt20", "ExN_Car12_Met", "ExN_Sim1_Strip2")
sankey_plot <- create_sankey_diagram(
  data = filtered_freq_df2,
  highlight_cell_types = highlight_cell_types,  
  sankey_colors = sankey_color,  
  plot_title = "Sankey Diagram of Cell Types and Brain Regions"
)




highlight_cell_types <- c("InN_Tshz1_Lamp5", "InN_Tshz1_Drd1", "InN_Sst_Vipr2", "InN_Drd2_Adora2a", "InN_Prkcd_Chodl", "InN_Isl1_Kcnh8")
sankey_plot <- create_sankey_diagram(
  data = filtered_freq_df2,
  highlight_cell_types = highlight_cell_types,  # 示例细胞类型
  sankey_colors = sankey_color,  # 假设你已经定义了颜色向量
  plot_title = "Sankey Diagram of Cell Types and Brain Regions"
)


ggsave(filename = paste0("/data2/xyh_desktop/xyh_desktop/yb/202504desktop/20250407/sankey/group2_0.5.pdf"), plot = sankey_plot, width = 10, height = 20)


#group3

highlight_cell_types <- c("InN_Gal_Dock8", "ExN_Avp_Ebf3", "ExN_Otp_Stc1", "ExN_Otp_Sim1", "ExN_Rspo1_Col18a1", "ExN_Rorb_Plch1", "ExN_Ebf2_Aox3", "ExN_Lhx9_Abi3bp", "ExN_Lhx5_Reln", "ExN_Tfap2c_Pou3f1")
sankey_plot <- create_sankey_diagram(
  data = filtered_freq_df2,
  highlight_cell_types = highlight_cell_types,  
  sankey_colors = sankey_color,  
  plot_title = "Sankey Diagram of Cell Types and Brain Regions"
)


#group4

highlight_cell_types <- c("InN_Tshr_Qrfprl", 
                          "InN_Apob_Prkcq", "InN_Lhx8_Th", "InN_Nxph2_Foxp2", "InN_Abcb5_Pax6", 
                          "InN_Calcr_Esr2", "InN_Esr2_Cyp19a1")
sankey_plot <- create_sankey_diagram(
  data = filtered_freq_df2,
  highlight_cell_types = highlight_cell_types,  
  sankey_colors = sankey_color,  
  plot_title = "Sankey Diagram of Cell Types and Brain Regions"
)

highlight_cell_types <- c("InN_Sst_Myzap", "InN_Sst_Ntn1", 
                          "InN_Sst_Edaradd", "InN_Pvalb_Adamts5", "InN_Pvalb_Vipr2")
sankey_plot <- create_sankey_diagram(
  data = filtered_freq_df2,
  highlight_cell_types = highlight_cell_types,  
  sankey_colors = sankey_color, 
  plot_title = "Sankey Diagram of Cell Types and Brain Regions"
)


#group6

highlight_cell_types <- c("InN_Lamp5_Ndnf", 
                          "InN_Lamp5_Id2", "InN_Vip_Crh", "InN_Htr3a_Glyctk", "InN_Htr3a_Frem1", 
                          "InN_Ndnf_Ppp1r1c")
sankey_plot <- create_sankey_diagram(
  data = filtered_freq_df2,
  highlight_cell_types = highlight_cell_types,  
  sankey_colors = sankey_color,  
  plot_title = "Sankey Diagram of Cell Types and Brain Regions"
)





