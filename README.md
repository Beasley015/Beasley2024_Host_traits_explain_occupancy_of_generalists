# HOST TRAITS EXPLAIN MORE VARIATION IN OCCUPANCY OF “GENERALISTS” THAN “SPECIALISTS” DUE TO STRONG HOST PREFERENCES AMONG GENERALISTS

## **Abstract**

The range of hosts a parasite can successfully occupy is partially determined by their niche breadth, the set of environmental conditions necessary to maintain a stable population. Niche breadth is often quantified using host specificity, which encompasses the number of host species a parasite can exploit and the parasite’s distribution among its hosts. Parasites with a wider niche breadth can potentially occupy more host species and are often more evenly distributed among their hosts than parasites with a narrower niche breadth. However, parasites interact with potential hosts within the context of a geographic locality and the set of environmental characteristics associated with it.  The extent to which environmental filters associated with host individuals and their geographic context explains variation in occupancy of parasites, and the extent to which variation in occupancy is associated with host range and specificity, is poorly understood. Using data from small mammals and ectoparasites in Vermont, I used a multi-scale, multi-species occupancy model (MSOM) to 1) estimate ectoparasite occupancy at 10 geographic sites and on individual hosts within each site, 2) quantify the variation in occupancy explained by the site and host levels of the model using Bayesian R2, and 3) evaluate associations between explained variation and host range of ectoparasites. For ectoparasites collected from at least 4 different host species, I calculated structural specificity to determine the distribution of these parasites across their hosts, and β-specificity to evaluate changes in host use across habitats. Host range was significantly associated with host-level Bayesian R2: generalist parasites had more variation in occupancy explained by host-level covariates than specialist parasites. This result may be explained by differences in structural specificity: many generalists disproportionally occurred on a single host species, suggesting that host characteristics act as habitat filters for these parasites. There were no significant associations between site-level Bayesian R2 and host specificity. However, many generalists demonstrated high β-specificity, suggesting these parasites may “switch” hosts depending on host availability. These results highlight that the terms “specialist” and “generalist” are context-dependent and may not accurately describe the niche breadth of parasite taxa. Understanding variation in host specificity as it pertains to potential habitat filters may be important for predicting which parasites are able to bypass host filters and “jump” to a novel host, which has implications for the surveillance and management of vector-borne diseases.

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

### UpdatedTracks.kml
Polyline file of trapping transects. Used to generate maps and test for spatial autocorrelation between sites.

## **Code**

EctoMultiScale.R - Full code for all analysis associated with the manuscript, including data cleaning, multi-scale MSOM, model selection, Bayesian R^2^, and figure creation. Code for saving intermediate outputs is included but is commented out.
