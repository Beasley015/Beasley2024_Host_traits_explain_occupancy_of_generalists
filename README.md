# Ectoparasite life history traits influence occupancy patterns at varying organizational scales

## **Abstract**

Ectoparasites are exposed to a ‘dual’ environment: the conditions on an individual host and the external environment. However, variation in ectoparasite life history traits, such as the portion of the life cycle spent on-host, leads to differences in selective pressure exerted by each environment. Parasites that spend most of the life cycle on-host are pressured to undergo increased host specialization, leading to differences in host specificity and occupancy patterns compared to ephemeral parasites which only contact the host to feed. Using data from small mammals and ectoparasites in Vermont, I used a multi-scale MSOM to 1) estimate ectoparasite occupancy on individual hosts nested within geographic sites, 2) calculate the Bayesian R2 at the site and host levels of the model to determine the variation in occupancy explained by each level, and 3) compared number of host species and the R2 values across different life history categories. Life history category was significantly associated with host specificity and host-level Bayesian R2: parasites which spend their full life cycle on the host parasitized fewer host species and had significantly more variation in occupancy explained by host-level covariates than ephemeral parasites. However, there were no significant differences in site-level R2 between life history categories, and no significant associations between site-level R2 and host specificity, suggesting that additional factors may play a role in structuring small mammal/ectoparasite communities.

## **Data**

All data used in the study was collected by the author from April-September 2020 in Chittenden County, Vermont. Data from related studies collected in 2019 and data collected from UVM Mammalogy field trips may also be included in some data files.

### MammRawData.csv

Contains data from small mammals captured in Sherman live traps. Key to column abbreviations:

Habitat: Coarse descriptor of habitat where trapping transect was located. Farm = 
actively farmed land, Field = old field, Forest = second-growth forest.

County: County where transect was located, no spaces.

Site: Location of trapping transect. StM = St. Michael's College (Colchester, VT), Jericho = University of Vermont Jericho Research Forest (Jericho, VT), Audubon = Green Mountain Audubon Center (Huntington, VT), Elmer = Elmer Farm (near Middlebury, VT), RiverBerry = River Berry Farm (Fairfax, VT), Barr = Barr Hill Natural Area (Greensboro, VT), MBR = Marsh-Billings-Rockefeller National Park, Buck = Helen Buckner Memorial Preserve (West Haven, VT), Dummer = Fort Dummer State Park, Hort = University of Vermont Horticultural Research Center (South Burlington, VT), Butternut = Butternut Hill Natural Area (North Hero, VT), Intervale = Intervale Center (Burlington, VT)

Date: date (day/month/year)

Day: Survey number, values range from 1-3

Techs: Initials of technicians collecting data. Author = EB

Station: Trap station along the linear transect

Bait: Bait type. PB = peanut butter + sunflower seeds, seed = sunflower seeds only

Tag: For individals marked with ear tags, tag number. For shrews, portion of the body
where fur was clipped for marking. LH = left hind, RH = right hind, CH = center hind

Genus: Genus of captured mammal

Species: species epithet of captured mammal

Abbrev: Abbreviation of species name. First two letters of genus + first two letters
of species epithet

Sex: F = female, M = male

Repro: reproductive status of captured individual. ns = non-scrotal male, s = scrotal 
male, sm = female with small nipples, lg = female with large nipples

Mass: mass in grams

TL: tail length in millimeters

RHFL: right hind foot length in millimeters

EL: ear length in millimeters. Measurement was taken from the left ear unless ear was
damaged

Ecto: Sample number of ectoparasites, if any were collected

DNA: Sample number for ear tissue sample, if any was collected

### Nocap_sites.csv

A list of sites and trap days on which no mammals were captured. This file makes it easier to set up the data for analysis. See above for column abbreviations. 

### EctosAll.csv

Data for ectoparasites collected on both live and deceased small mammals collected in summer 2020. Key to column abbreviations:

SampleNo.: ID number that matches a sample collected from a live mammal that was released after data collection (format E###) or a mammal that was found dead and prepared and deposited as a voucher specimen in UVM's Zadock Thompson Zoological Collections (format AAA###). 

Order: Order of an individual ectoparasite.

Family: Family of an individual ectoparasite.

Genus: Genus of an individual ectoparasite.

Species: Species of an individual ectoparasite.

### VegRawData.csv

Vegetation data sampled on every 3rd trapping station along the small mammal trapping transect. Contains data from this project and a related project from 2019. Key to column abbreviations:

Habitat: Coarse descriptor of habitat where trapping transect was located. Farm = 
actively farmed land, Field = old field, Forest = second-growth forest.

County: County where transect was located, no spaces.

Site: Location of trapping transect. StM = St. Michael's College (Colchester, VT), Jericho = University of Vermont Jericho Research Forest (Jericho, VT), Audubon = Green Mountain Audubon Center (Huntington, VT), Elmer = Elmer Farm (near Middlebury, VT), RiverBerry = River Berry Farm (Fairfax, VT), Barr = Barr Hill Natural Area (Greensboro, VT), MBR = Marsh-Billings-Rockefeller National Park, Buck = Helen Buckner Memorial Preserve (West Haven, VT), Dummer = Fort Dummer State Park, Hort = University of Vermont Horticultural Research Center (South Burlington, VT), Butternut = Butternut Hill Natural Area (North Hero, VT), Intervale = Intervale Center (Burlington, VT)

Date: date (day/month/year)

Techs: Initials of technicians collecting data. Author = EB

Station: Trap station along the linear transect

Canopy: Index of canopy cover as measured by a spherical convex densiometer.

Weins10 to Weins60: Vertical complexity as measured using a Weins pole. Vertical complexity is measured by counting the number of times living vegitation touches each 10-centimeter interval makred on the pole. The numbers in the column names indicate the interval: e.g. Weins10 is the 0-10cm interval.

%Tree to %DeadVeg: Percent ground cover composition within a grid measuring 0.5m by 0.5m. 

Other: Data that does not fit into any other column. Primarily includes percentages of relatively rare ground cover types (e.g. moss, water)

### FleaClassifications.csv

Flea classifications were gleaned from the literature (see Appendix S1, Table 1 of the manuscript for literature cited). 

## **Code**

EctoMultiScale.R - Full code for all analysis associated with the manuscript, including data cleaning, multi-scale MSOM, model selection, Bayesian R^2^, and figure creation. Code for saving intermediate outputs is included but is commented out.
