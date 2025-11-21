# mada-intestinal-parasites-prev

This repo contains code for calculating and visualizing the prevalence of multiple intestinal parasites across diverse environments in Madagascar. The data are obtained from the following three studies.

| Study | Citation |
|---------------|---------------------------------------------------------|
| MAHERY NE Cohort | Golden CD, Anjaranirina EJG, Fernald LCH, Hartl DL, Kremen C, Milner DA Jr, Ralalason DH, Ramihantaniarivo H, Randriamady H, Rice BL, Vaitla B, Volkman SK, Vonona MA, Myers SS. Cohort Profile: The Madagascar Health and Environmental Research (MAHERY) study in north-eastern Madagascar. Int J Epidemiol. 2017 Dec 1;46(6):1747-1748d. doi: 10.1093/ije/dyx071. PMID: 29040632; PMCID: PMC5837654. |
| MAHERY-Antongil (a.k.a. Darwin) NE Cohort | Golden CD, Borgerson C, Rice BL, Allen LH, Anjaranirina EJG, Barrett CB, Boateng G, Gephart JA, Hampel D, Hartl DL, Knippenberg E, Myers SS, Ralalason DH, Ramihantaniarivo H, Randriamady H, Shahab-Ferdows S, Vaitla B, Volkman SK, Vonona MA. Cohort Description of the Madagascar Health and Environmental Research-Antongil (MAHERY-Antongil) Study in Madagascar. Front Nutr. 2019 Jul 19;6:109. doi: 10.3389/fnut.2019.00109. PMID: 31428615; PMCID: PMC6690017. |
| CRS Cross-Sectional Study | Golden CD, Rice BL, Randriamady HJ, Vonona AM, Randrianasolo JF, Tafangy AN, Andrianantenaina MY, Arisco NJ, Emile GN, Lainandrasana F, Mahonjolaza RFF, Raelson HP, Rakotoarilalao VR, Rakotomalala AANA, Rasamison AD, Mahery R, Tantely ML, Girod R, Annapragada A, Wesolowski A, Winter A, Hartl DL, Hazen J, Metcalf CJE. Study Protocol: A Cross-Sectional Examination of Socio-Demographic and Ecological Determinants of Nutrition and Disease Across Madagascar. Front Public Health. 2020 Sep 17;8:500. doi: 10.3389/fpubh.2020.00500. PMID: 33042943; PMCID: PMC7527467. |

Materials are split into the following folders:

-   `Cleaned_Data`: Cleaned data ready for analysis

-   `Code`: Code to calculate overall prevalence across all datasets (`1_get_overall_prevs.R`), calculate age- and sex-specific prevalences (`2_age_sex_prevs.R`), and calculate region-specific prevalences (`3_region_prevs.R`). Code for summarizing the model results is also available in this folder.

-   `Model_Outputs`: Figures and tables produced using the outputs from the analysis models.
