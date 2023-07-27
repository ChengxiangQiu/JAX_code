
##########################################
### parameters used for making 3D UMAP ###
##########################################

t1 = list(family = 'Helvetica',
          size = 25,
          color = "black")
t2 = list(family = 'Helvetica',
          size = 15,
          color = "grey")

####################################
### major trajectory color plate ###
####################################

major_trajectory_color_plate = c("Neuroectoderm_and_glia"             = "#f96100",
                                 "Intermediate_neuronal_progenitors"  = "#2e0ab7",
                                 "Eye_and_other"                      = "#00d450",
                                 "Ependymal_cells"                    = "#b75bff",
                                 "CNS_neurons"                        = "#e5c000",
                                 "Mesoderm"                           = "#bb46c5",
                                 "Definitive_erythroid"               = "#dc453e",
                                 "Epithelial_cells"                   = "#af9fb6",
                                 "Endothelium"                        = "#00a34e",
                                 "Muscle_cells"                       = "#ffa1f5",
                                 "Hepatocytes"                        = "#185700",
                                 "White_blood_cells"                  = "#7ca0ff",
                                 "Neural_crest_PNS_glia"              = "#fff167",
                                 "Adipocytes"                         = "#7f3e39",
                                 "Primitive_erythroid"                = "#ffa9a1",
                                 "Neural_crest_PNS_neurons"           = "#b5ce92",
                                 "T_cells"                            = "#ff9d47",
                                 "Lung_and_airway"                    = "#02b0d1",
                                 "Intestine"                          = "#ff007a",
                                 "B_cells"                            = "#01b7a6",
                                 "Olfactory_sensory_neurons"          = "#e6230b",
                                 "Cardiomyocytes"                     = "#643e8c",
                                 "Oligodendrocytes"                   = "#916e00",
                                 "Mast_cells"                         = "#005361",
                                 "Megakaryocytes"                     = "#3f283d",
                                 "Testis_and_adrenal"                 = "#585d3b")

######################################
### Timepoints between E8.5 and P0 ###
######################################

day_color_plate  =  c("E8.0"  = "#5E4FA2", "E8.25"  = "#545BA8", "E8.5"  =  "#4A68AE", "E8.75"  = "#4075B4",
                      "E9.0"  = "#3682BA", "E9.25"  = "#398FB9", "E9.5"  =  "#449DB4", "E9.75"  = "#50AAAE",
                      "E10.0" = "#5CB7A9", "E10.25" = "#69C3A4", "E10.5" =  "#78C9A4", "E10.75" = "#88CFA4",
                      "E11.0" = "#98D5A4", "E11.25" = "#A7DBA4", "E11.5" =  "#B5E1A1", "E11.75" = "#C3E69F",
                      "E12.0" = "#D0EC9C", "E12.25" = "#DDF199", "E12.5" =  "#E8F59B", "E12.75" = "#EDF8A4",
                      "E13.0" = "#F3FAAD", "E13.25" = "#F9FCB6", "E13.5" =  "#FFFFBF", "E13.75" = "#FEF7B3",
                      "E14.0" = "#FEF0A7", "E14.25" = "#FEE99B", "E14.333" ="#FEE28F", "E14.75" = "#FDD985",
                      "E15.0" = "#FDCD7B", "E15.25" = "#FDC272", "E15.5" =  "#FDB768", "E15.75" = "#FCAB5F",
                      "E16.0" = "#FA9C58", "E16.25" = "#F88D52", "E16.5" =  "#F67E4B", "E16.75" = "#F46F44",
                      "E17.0" = "#EE6445", "E17.25" = "#E75947", "E17.5" =  "#E04F4A", "E17.75" = "#D9444D",
                      "E18.0" = "#CF384D", "E18.25" = "#C32A4A", "E18.5" =  "#B71C47", "E18.75" = "#AA0E44", 
                      "P0"    = "#9E0142")

#################################
### Somite counts color plate ###
#################################

somite_color_plate = rev(c("#fde725", "#eae51a", "#d5e21a", "#c0df25", "#a8db34", 
                           "#93d741", "#7fd34e", "#6ccd5a", "#58c765", "#48c16e", 
                           "#3aba76", "#2eb37c", "#25ab82", "#20a386", "#1e9c89", 
                           "#1f948c", "#228c8d", "#25848e", "#287d8e", "#2b758e", 
                           "#2e6d8e", "#32658e", "#365d8d", "#3a548c", "#3e4a89", 
                           "#424186", "#453882", "#472e7c", "#482374", "#48186a", 
                           "#470d60", "#440154"))

names(somite_color_plate) = paste0(c(0, 2:12, 14:18, 20:34), " somites")

####################################
### Posterior embryo color plate ###
####################################

posterior_embryo_color_plate = c("Notochord"                           = "#1F77B4",
                                 "Nodal cilia"                         = "#FF7F0E", 
                                 "NMPs and spinal cord progenitors"    = "#2CA02C",
                                 "Gut"                                 = "#D62728",
                                 "Mesodermal progenitors (Tbx6+)"      = "#9467BD")

#########################
### renal color plate ###
#########################

renal_color_plate = c("Posterior intermediate mesoderm"    = "#A6CEE3",
                      "Anterior intermediate mesoderm"     = "#1F78B4",
                      "Ureteric bud"                       = "#B2DF8A",
                      "Collecting duct intercalated cells" = "#33A02C",
                      "Collecting duct principal cells"    = "#FB9A99",
                      "Metanephric mesenchyme"             = "#E31A1C",
                      "Nephron progenitors"                = "#FDBF6F",
                      "Podocytes"                          = "#FF7F00",
                      "Proximal tubule cells"              = "#CAB2D6",
                      "Ascending loop of Henle"            = "#6A3D9A",
                      "Distal convoluted tubule"           = "#FFFF99",
                      "Connecting tubule"                  = "#B15928")


##########################################
### lateral plate mesoderm color plate ###
##########################################

LPM_color_plate = c("Splanchnic mesoderm"                   = "#ca47cc",
                    "Lung mesenchyme"                       = "#6dd251",
                    "Hepatic mesenchyme"                    = "#663fc6",
                    "Renal stromal cells"                   = "#f8f434",
                    "Cardiopharyngeal mesoderm"             = "#3f2b71",
                    "Vascular smooth muscle cells"          = "#c1bf5a",
                    "Vascular smooth muscle cells (Pparg+)" = "#9f7ada",
                    "Meninges"                              = "#6ad39e",
                    "Airway smooth muscle cells"            = "#dc4989",
                    "Amniotic mesoderm"                     = "#587f3a",
                    "Allantois"                             = "#8f3e7f",
                    "Extraembryonic mesoderm"               = "#d58d40",
                    "Somatic mesoderm"                      = "#7d9cd9",
                    "Renal pericytes and mesangial cells"   = "#d94e38",
                    "Gut mesenchyme"                        = "#79bec6",
                    "Foregut mesenchyme"                    = "#983b40",
                    "Gastrointestinal smooth muscle cells"  = "#7f5d35",
                    "Gonad progenitor cells"                = "#42222e",
                    "Sertoli cells"                         = "#d399af",
                    "Granulosa cells"                       = "#516aca",
                    "Proepicardium"                         = "#5b6180",
                    "Mesothelial cells"                     = "#3a5440")

LPM_E85_color_plate = c("Anterior intermediate mesoderm"               = "#1F78B4",
                        "First heart field"                            = "#d85091",
                        "Lateral plate and intermediate mesoderm"      = "#7d5a3b",
                        "Mesodermal progenitors (Tbx6+)"               = "#38503a",
                        "Posterior intermediate mesoderm"              = "#A6CEE3",
                        "Second heart field"                           = "#3d2130",
                        "LPM:Allantois"                                = "#8f3e7f",
                        "LPM:Amniotic mesoderm"                        = "#587f3a",
                        "LPM:Cardiopharyngeal mesoderm"                = "#3f2b71",
                        "LPM:Extraembryonic mesoderm"                  = "#d58d40",
                        "LPM:Foregut mesenchyme"                       = "#983b40",
                        "LPM:Gut mesenchyme"                           = "#79bec6",
                        "LPM:Hepatic mesenchyme"                       = "#663fc6",
                        "LPM:Proepicardium"                            = "#5b6180",
                        "LPM:Renal pericytes and mesangial cells"      = "#d94e38",
                        "LPM:Somatic mesoderm"                         = "#7d9cd9",
                        "LPM:Splanchnic mesoderm"                      = "#ca47cc",
                        "NMPs and spinal cord progenitors"             = "#cb566d",
                        "LPM:Gonad progenitor cells"                   = "#42222e",
                        "LPM:Lung mesenchyme"                          = "#6dd251",
                        "LPM:Meninges"                                 = "#6ad39e",
                        "LPM:Vascular smooth muscle cells"             = "#c1bf5a",
                        "LPM:Mesothelial cells"                        = "#3a5440")


#######################
### eye color plate ###
#######################

eye_color_plate = c("Eye field"                            = '#6d8840', 
                    "Retinal pigment cells"                = '#c95493',
                    "Naive retinal progenitor cells"       = "#A6CEE3",
                    "Retinal progenitor cells"             = "#1F78B4",
                    "Bipolar precursor cells"              = "#B2DF8A",
                    "Ciliary margin cells"                 = "#33A02C",
                    "Photoreceptor precursor cells"        = "#FB9A99",
                    "Cone precursor cells"                 = "#E31A1C",
                    "Rod precursor cells"                  = "#FDBF6F",
                    "Retinal ganglion cells"               = "#FF7F00",
                    "PV-containing retinal ganglion cells" = "#CAB2D6",
                    "Amacrine/Horizontal precursor cells"  = "#6A3D9A",
                    "Amacrine cells"                       = "#FFFF99",
                    "Horizontal cells"                     = "#B15928",
                    "Cholinergic amacrine cells"           = "#808080")

#################################
### neuroectoderm color plate ###
#################################

neuroectoderm_color_plate = c("Spinal dI1 interneurons"                = '#636EFA',
                              "Spinal dI2 interneurons"                = '#EF553B',
                              "Spinal dI3 interneurons"                = '#00CC96',
                              "Spinal dI4 interneurons"                = '#AB63FA',
                              "Spinal dI5 interneurons"                = '#FFA15A',
                              "Spinal dI6 interneurons"                = '#19D3F3',
                              "Spinal V0 interneurons"                 = '#FF6692',
                              "Spinal V1 interneurons"                 = '#B6E880',
                              "Spinal V2a interneurons"                = '#FECB52',
                              "Spinal V2b interneurons"                = '#1F77B4',
                              "Spinal V3 interneurons"                 = '#FF7F0E',
                              "Di/mesencephalon glutamatergic neurons" = '#2CA02C',
                              "Striatal projection neurons"            = '#9467BD',
                              "Precerebellar neurons"                  = '#8C564B',
                              "Di/mesencephalon GABAergic neurons"     = '#E377C2',
                              "Hypothalamic Sim1 neurons"              = '#7F7F7F',
                              "Midbrain dopaminergic neurons"          = '#BCBD22',
                              "Spinal cord motor neuron progenitors"   = '#17BECF',
                              "Spinal cord motor neurons"              = "#17BECF",
                              
                              "Anterior floor plate"                   = "#9a632b",
                              "Anterior roof plate"                    = "#b1c4d9",
                              "Diencephalon"                           = "#9e3136",
                              "Dorsal telencephalon"                   = "#93d6b6",
                              "Hindbrain"                              = "#8f366c",
                              "Hypothalamus"                           = "#848b38",
                              "Hypothalamus (Sim1+)"                   = "#52288d",
                              "Midbrain"                               = "#d3c398",
                              "Midbrain-hindbrain boundary"            = "#3d2945",
                              "Floorplate and p3 domain"               = "#b898c9",
                              "Posterior roof plate"                   = "#363d27",
                              "Spinal cord/r7/r8"                      = "#467d4e",
                              "Telencephalon"                          = "#d98f90",
                                  
                              "Astrocytes"                             = "#504d84",
                              "Cajal-Retzius cells"                    = "#cfd98e",
                              "Choroid plexus"                         = "#692135",
                              "Cranial motor neurons"                  = "#65e3be",
                              "Eye field"                              = "#c15a34",
                              "GABAergic cortical interneurons"        = "#e8afd6",
                              "Intermediate progenitors"               = "#45251f",
                              "Neural progenitor cells (Neurod1+)"     = "#d7c1b7",
                              "Neurons (Slc17a8+)"                     = "#437a7d",
                              "Thalamic neuronal precursors"           = "#5f5b28")

astrocytes_color_plate = c("VA1 astrocytes"      = "#9ac25e",
                           "VA2 astrocytes"      = "#b35948",
                           "VA3 astrocytes"      = "#b98b3a",
                           "Anterior astrocytes" = "#8d50a9")

neuron_day_color_plate = rev(c("#fcffa4", "#f2ea69", "#f9cb35", "#fcac11", "#f98e09", 
                               "#f1731d", "#e45a31", "#d24644", "#bc3754", "#a32c61", 
                               "#8a226a", "#71196e", "#57106e", "#3d0965", "#210c4a", "#0b0724", "#000004"))

names(neuron_day_color_plate) = c("E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
                                  "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
                                  "E12.25", "E12.5", "E12.75")



################################
### Birth-series color plate ###
################################

birth_color_plate = c("E18.75"       = "#7aa457",
                      "P0"           = "#cb6a49",
                      "Csection_0m"  = "#0d0887",
                      "Csection_20m" = "#7e03a8",
                      "Csection_40m" = "#cc4778",
                      "Csection_60m" = "#f89540",
                      "Csection_80m" = "#f0f921",
                      "NatBirth"     = "#02b0d1")


######################################
### Re-annotation P-S + E8.5b data ###
######################################

gastrulation_color_plate = c("Allantois"                        = "#8f3e7f",
                             "Amniotic mesoderm"                = "#587f3a",
                             "Extraembryonic mesoderm"          = "#d58d40",
                             "Gut mesenchyme"                   = "#79bec6",
                             "Hindbrain"                        = "#8f366c",
                             "Primitive erythroid cells"        = "#c96d44",
                             "Blood progenitors"                = "#5d69d7",
                             "Cardiopharyngeal mesoderm"        = "#e68028",
                             "Definitive ectoderm"              = "#4e7ac1",
                             "Embryonic visceral endoderm"      = "#d2482a",
                             "Epiblast"                         = "#58abe3",
                             "Floor plate"                      = "#d13e4b",
                             "Gut"                              = "#54c465",
                             "Hematoendothelial progenitors"    = "#c843ae",
                             "Intermediate mesoderm"            = "#77b635",
                             "Midbrain"                         = "#d77fdc",
                             "Neural crest"                     = "#d3a52d",
                             "Paraxial mesoderm (Tbx6-)"        = "#9652c6",
                             "Parietal endoderm"                = "#a8ad38",
                             "Primitive streak"                 = "#7d5fa4",
                             "Spinal cord"                      = "#388a35",
                             "Amniotic ectoderm"                = "#dc4074",
                             "Anterior primitive streak"        = "#5cbf8c",
                             "Cardiogenic mesoderm"             = "#d64694",
                             "CLE and NMPs"                     = "#41bcb3",
                             "Definitive endoderm"              = "#994f2d",
                             "Endothelium"                      = "#af99e2",
                             "Extraembryonic ectoderm"          = "#956f2c",
                             "Extraembryonic visceral endoderm" = "#ae5c82",
                             "Forebrain"                        = "#9cb36a",
                             "Lateral plate mesoderm"           = "#e389b6",
                             "Nascent mesoderm"                 = "#318766",
                             "Notochord"                        = "#e48380",
                             "Paraxial mesoderm (Tbx6+)"        = "#737129",
                             "Primordial germ cells"            = "#ac4b55",
                             "Surface ectoderm"                 = "#cea365")


